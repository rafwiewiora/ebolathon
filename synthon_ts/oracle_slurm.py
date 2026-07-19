"""Docking oracle backend #4 -- your OWN SLURM cluster (CPU->GPU funnel).

Runs docking on an HPC with SLURM instead of the cloud backends (Rowan / muni),
which are storage/queue-limited. One `score()` call becomes ONE SLURM **array
job** (one array task per `array_chunk` ligands), so a whole round of the
synthon-TS screen fans out across every CPU/GPU node you have.

Engines (all free; a CPU->GPU funnel; selectable via `engine=`):

* ``qvina2`` -- CPU / cheap bulk Vina (matches the CDK2/1HCK calibration).
* ``smina``  -- CPU Vina/Vinardo (same PDBQT prep as qvina2; ``--scoring vinardo``).
* ``unidock`` -- GPU throughput (Uni-Dock: GPU Vina/Vinardo, thousands/GPU).
* ``gnina``  -- GPU accuracy (GNINA 1.3: Vina + CNN rescoring; SDF ligands).

Oracle contract (identical to the muni/Rowan backends)::

    oracle.score(smiles_list, template_smiles=None) -> {smiles: score}

lower score = better (Vina/Vinardo kcal/mol); failed ligands are omitted.

What runs WHERE
---------------
* Login node (this process): receptor prep, RDKit ligand prep, sbatch, poll,
  score parsing. Assumes a **shared filesystem** (the norm on HPC) so compute
  nodes read the ligand files and write the pose/score files this process reads
  back.
* Compute nodes (the array tasks): the chosen engine, over its chunk, using the
  rendered ``slurm/*.sbatch`` template.

Cluster specifics (partitions, GPU gres, module/conda activation) are all
configurable and marked ``TODO: set for your cluster``.
"""
from __future__ import annotations

import glob
import math
import os
import re
import shutil
import subprocess
import time
import uuid

from .core import Target
from .pocket import _fetch_structure  # reuse the RCSB download / local-file helper

try:  # RDKit is required for ligand prep; guarded so import of the module never fails
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _HAVE_RDKIT = True
except Exception:  # pragma: no cover
    _HAVE_RDKIT = False

try:  # Meeko is the preferred SMILES->PDBQT path; obabel is the fallback
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    _HAVE_MEEKO = True
except Exception:  # pragma: no cover
    _HAVE_MEEKO = False


# --------------------------------------------------------------------------- #
# Engine metadata
# --------------------------------------------------------------------------- #
_GPU_ENGINES = {"unidock", "gnina"}
_LIGAND_FMT = {"qvina2": "pdbqt", "smina": "pdbqt", "unidock": "pdbqt", "gnina": "sdf"}
_TEMPLATE = {"qvina2": "cpu_qvina.sbatch", "smina": "cpu_qvina.sbatch",
             "unidock": "gpu_unidock.sbatch", "gnina": "gpu_gnina.sbatch"}
_ENGINE_BIN = {"qvina2": "qvina2", "smina": "smina"}
_SLURM_DIR = os.path.join(os.path.dirname(__file__), "slurm")


# --------------------------------------------------------------------------- #
# Ligand prep (RDKit embed-only -> PDBQT/SDF). Module-level = locally testable.
# --------------------------------------------------------------------------- #
def embed_3d(smiles: str, seed: int = 0xF00D):
    """SMILES -> a single 3D conformer with RDKit ETKDGv3.

    HARD RULE (user global): **embed only -- never run MMFF/UFF force-field
    optimization**. We add Hs, embed one conformer, and stop. A random-coords
    retry covers hard embeds. Returns an RDKit Mol (with Hs) or None."""
    if not _HAVE_RDKIT:
        raise RuntimeError("RDKit is required for ligand prep (pip install rdkit)")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(m, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(m, params) != 0:
            return None
    # NO force-field optimization here, by design (embed-only rule).
    return m


def ligand_to_sdf(smiles: str, path: str, seed: int = 0xF00D) -> bool:
    """Embed `smiles` and write a 3D SDF (for gnina). Returns True on success."""
    m = embed_3d(smiles, seed)
    if m is None:
        return False
    w = Chem.SDWriter(path)
    w.write(m)
    w.close()
    return os.path.exists(path) and os.path.getsize(path) > 0


def ligand_to_pdbqt(smiles: str, path: str, seed: int = 0xF00D,
                    obabel_bin: str = "obabel") -> bool:
    """Embed `smiles` and write a PDBQT (for qvina2/smina/unidock).

    Prefers Meeko (`MoleculePreparation` + `PDBQTWriterLegacy`); falls back to
    Open Babel if Meeko isn't installed. Returns True on success."""
    m = embed_3d(smiles, seed)
    if m is None:
        return False
    if _HAVE_MEEKO:
        try:
            prep = MoleculePreparation()
            setups = prep.prepare(m)
            setup = setups[0] if isinstance(setups, (list, tuple)) else setups
            pdbqt, ok, _err = PDBQTWriterLegacy.write_string(setup)
            if ok and pdbqt:
                with open(path, "w") as fh:
                    fh.write(pdbqt)
                return os.path.getsize(path) > 0
        except Exception:  # noqa: BLE001 -- fall through to obabel
            pass
    return _obabel_from_mol(m, path, obabel_bin)


def _obabel_from_mol(mol, out_path: str, obabel_bin: str = "obabel") -> bool:
    """Write `mol` to a temp SDF (RDKit) then convert to `out_path` with obabel.
    obabel builds the AutoDock torsion tree + Gasteiger charges for .pdbqt."""
    if not shutil.which(obabel_bin):
        return False
    import tempfile
    tmp = tempfile.NamedTemporaryFile("w", suffix=".sdf", delete=False)
    tmp.close()
    try:
        w = Chem.SDWriter(tmp.name)
        w.write(mol)
        w.close()
        subprocess.run([obabel_bin, tmp.name, "-O", out_path],
                       capture_output=True, text=True)
        return os.path.exists(out_path) and os.path.getsize(out_path) > 0
    finally:
        os.unlink(tmp.name)


# --------------------------------------------------------------------------- #
# Score parsing (module-level = locally testable against saved engine output)
# --------------------------------------------------------------------------- #
_VINA_RE = re.compile(r"REMARK\s+VINA\s+RESULT:\s*(-?\d+\.?\d*)")


def parse_vina_pdbqt(text: str):
    """Best (lowest) score from a Vina-family output PDBQT.

    qvina2/smina/unidock write one ``REMARK VINA RESULT: <score> <rmsd_lb>
    <rmsd_ub>`` per pose (mode 1 first = best). We take the min to be safe.
    Returns a float or None if no result line is present."""
    scores = [float(s) for s in _VINA_RE.findall(text or "")]
    return min(scores) if scores else None


def parse_gnina_sdf(text: str):
    """Best (lowest) ``minimizedAffinity`` from a gnina output SDF.

    gnina writes each pose as an SDF record with a ``> <minimizedAffinity>``
    property block (Vina kcal/mol, lower = better) plus CNN scores. We scan all
    poses and return the minimum minimizedAffinity, or None."""
    if not text:
        return None
    lines = text.splitlines()
    vals = []
    for i, ln in enumerate(lines):
        if "minimizedAffinity" in ln and ln.lstrip().startswith(">"):
            for j in range(i + 1, min(i + 3, len(lines))):
                s = lines[j].strip()
                if s:
                    try:
                        vals.append(float(s.split()[0]))
                    except ValueError:
                        pass
                    break
    return min(vals) if vals else None


# --------------------------------------------------------------------------- #
# Receptor prep helpers
# --------------------------------------------------------------------------- #
def _strip_receptor_pdb(text: str) -> str:
    """Keep only protein ATOM/TER records -- drop HETATM (waters, ligands, ions,
    buffer). A pragmatic default; if your pocket needs a metal/cofactor kept,
    prep the receptor PDBQT yourself and pass its path as the target protein.
    TODO: revisit per target if a cofactor is catalytically required."""
    keep = []
    for ln in text.splitlines():
        if ln.startswith("ATOM") or ln.startswith("TER"):
            keep.append(ln)
    keep.append("END")
    return "\n".join(keep) + "\n"


# --------------------------------------------------------------------------- #
# The oracle
# --------------------------------------------------------------------------- #
class SlurmDockingOracle:
    """`score(smiles_list) -> {smiles: score}` via one SLURM array job per call.

    Parameters
    ----------
    target : `core.Target`. ``target.protein`` is a PDB id (downloaded from
        RCSB), a local ``.pdb``/``.cif`` path, or a pre-prepared ``.pdbqt`` path
        (used as-is). ``target.pocket`` is ``[[cx,cy,cz],[sx,sy,sz]]``.
        ``target.exhaustiveness`` and ``target.scoring_function`` drive docking.
    engine : ``qvina2`` | ``smina`` | ``unidock`` | ``gnina``.
    partition : SLURM partition for CPU engines (qvina2/smina).
        TODO: set for your cluster.
    gpu_partition : SLURM partition for GPU engines (unidock/gnina).
        TODO: set for your cluster.
    time_limit : ``#SBATCH --time`` per array task.
    array_chunk : ligands per array task (one task docks this many).
    poll_interval : seconds between ``squeue`` polls.
    workdir : shared-filesystem scratch dir (ligands, sbatch, poses, scores).
    env_activate : shell snippet prepended to every job (module load / micromamba
        activate). TODO: set for your cluster.
    log : logging callable.
    """

    def __init__(self, target: Target, engine: str = "qvina2",
                 partition: str | None = None, gpu_partition: str | None = None,
                 time_limit: str = "02:00:00", array_chunk: int = 200,
                 poll_interval: int = 30, workdir: str = "slurm_work",
                 env_activate: str | None = None, log=print):
        if engine not in _LIGAND_FMT:
            raise ValueError(f"unknown engine {engine!r}; choose from {sorted(_LIGAND_FMT)}")
        self.t = target
        self.engine = engine
        self.partition = partition
        self.gpu_partition = gpu_partition
        self.time_limit = time_limit
        self.array_chunk = max(1, int(array_chunk))
        self.poll_interval = poll_interval
        self.workdir = os.path.abspath(workdir)
        self.env_activate = env_activate
        self.log = log
        self._receptor_pdbqt: str | None = None
        # {smiles: (best_pose_file, best_score)} -- populated as rounds score,
        # consumed by write_top_hits_sdf_gz for lean pose output.
        self._best_pose: dict = {}

    # -- oracle contract ---------------------------------------------------
    def score(self, smiles_list, template_smiles=None) -> dict:
        """Dock `smiles_list` on the cluster; return {smiles: best_score}.
        `template_smiles` is ignored (the box defines the site), matching the
        muni/Rowan batch backends."""
        smiles_list = list(smiles_list)
        if not smiles_list:
            return {}
        receptor = self._prepare_receptor()
        run_dir = os.path.join(self.workdir, f"run_{uuid.uuid4().hex[:8]}")
        os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(run_dir, "out"), exist_ok=True)

        idx2smi, idx2file = self._prep_ligands(run_dir, smiles_list)
        if not idx2smi:
            self.log("[slurm] no ligands prepared this batch; nothing to submit")
            return {}
        n_tasks = self._write_chunks(run_dir, idx2smi, idx2file)
        self._render_sbatch(run_dir, receptor, n_tasks)
        self._submit_and_wait(run_dir, n_tasks)
        scores = self._collect(run_dir, idx2smi)
        self.log(f"[slurm] scored {len(scores)}/{len(idx2smi)} ligands "
                 f"(engine={self.engine}, {n_tasks} array tasks)")
        return scores

    # -- receptor ----------------------------------------------------------
    def _prepare_receptor(self) -> str:
        """PDB(id/path) -> receptor.pdbqt, cached once in `workdir`.

        If `target.protein` is already a `.pdbqt`, use it verbatim. Otherwise
        fetch (RCSB id or local file via `pocket._fetch_structure`), strip
        HETATM for PDB input, and convert with Meeko `mk_prepare_receptor.py`
        (if on PATH) else Open Babel."""
        if self._receptor_pdbqt and os.path.exists(self._receptor_pdbqt):
            return self._receptor_pdbqt
        prot = self.t.protein
        if isinstance(prot, str) and prot.lower().endswith(".pdbqt") and os.path.exists(prot):
            self._receptor_pdbqt = os.path.abspath(prot)
            return self._receptor_pdbqt

        os.makedirs(self.workdir, exist_ok=True)
        rec = os.path.join(self.workdir, "receptor.pdbqt")
        if os.path.exists(rec):
            self._receptor_pdbqt = rec
            return rec

        text, fmt = _fetch_structure(prot)
        src = os.path.join(self.workdir, f"receptor_src.{fmt}")
        if fmt == "pdb":
            text = _strip_receptor_pdb(text)
        with open(src, "w") as fh:
            fh.write(text)

        if not self._receptor_via_meeko(src, rec) and not self._receptor_via_obabel(src, rec):
            raise RuntimeError(
                "receptor prep failed: need Meeko's mk_prepare_receptor.py or "
                "obabel on PATH (install via synthon_ts/slurm/install.sh). "
                f"source structure: {src}")
        self.log(f"[slurm] prepared receptor -> {rec}")
        self._receptor_pdbqt = rec
        return rec

    def _receptor_via_meeko(self, src: str, rec: str) -> bool:
        exe = shutil.which("mk_prepare_receptor.py")
        if not exe:
            return False
        stem = rec[:-6] if rec.endswith(".pdbqt") else rec
        # Meeko CLI flags vary across versions; try the common invocations.
        for args in ([exe, "--read_pdb", src, "-o", stem, "-p"],
                     [exe, "-i", src, "-o", stem, "-p"]):
            try:
                subprocess.run(args, capture_output=True, text=True)
                if os.path.exists(rec) and os.path.getsize(rec) > 0:
                    return True
            except Exception:  # noqa: BLE001
                continue
        return False

    def _receptor_via_obabel(self, src: str, rec: str) -> bool:
        obabel = shutil.which("obabel")
        if not obabel:
            return False
        # -xr = rigid receptor PDBQT (Gasteiger charges added for .pdbqt output)
        subprocess.run([obabel, src, "-O", rec, "-xr"], capture_output=True, text=True)
        return os.path.exists(rec) and os.path.getsize(rec) > 0

    # -- ligands -----------------------------------------------------------
    def _prep_ligands(self, run_dir: str, smiles_list):
        """Embed + write one ligand file per SMILES. Returns
        ({index: smiles}, {index: relpath}) for the ones that prepped OK."""
        ligdir = os.path.join(run_dir, "ligands")
        os.makedirs(ligdir, exist_ok=True)
        fmt = _LIGAND_FMT[self.engine]
        idx2smi, idx2file = {}, {}
        for i, smi in enumerate(smiles_list):
            path = os.path.join(ligdir, f"lig_{i}.{fmt}")
            ok = (ligand_to_sdf(smi, path) if fmt == "sdf"
                  else ligand_to_pdbqt(smi, path))
            if ok:
                idx2smi[i] = smi
                idx2file[i] = os.path.relpath(path, run_dir)
            else:
                self.log(f"[slurm] ligand prep failed (idx {i}): {smi[:70]}")
        self.log(f"[slurm] prepped {len(idx2smi)}/{len(smiles_list)} ligands ({fmt})")
        return idx2smi, idx2file

    def _write_chunks(self, run_dir: str, idx2smi, idx2file) -> int:
        """Write one chunk manifest per array task. Returns the task count.

        Uni-Dock wants a bare ligand-index file (one path/line); the per-ligand
        engines (qvina2/smina/gnina) want ``<index> <path>`` lines."""
        indices = list(idx2smi)
        n_tasks = math.ceil(len(indices) / self.array_chunk)
        for t in range(n_tasks):
            sub = indices[t * self.array_chunk:(t + 1) * self.array_chunk]
            with open(os.path.join(run_dir, f"chunk_{t + 1}.txt"), "w") as fh:
                for idx in sub:
                    if self.engine == "unidock":
                        fh.write(f"{idx2file[idx]}\n")
                    else:
                        fh.write(f"{idx} {idx2file[idx]}\n")
        return n_tasks

    # -- sbatch render + submit -------------------------------------------
    def _render_sbatch(self, run_dir: str, receptor: str, n_tasks: int) -> str:
        """Fill the engine's template and write `run_dir/job.sbatch`."""
        gpu = self.engine in _GPU_ENGINES
        partition = (self.gpu_partition if gpu else self.partition)
        if not partition:
            partition = "TODO_SET_PARTITION"
            self.log("[slurm] WARNING: no partition set -- sbatch will fail until "
                     "you pass --partition / --gpu-partition")
        gres = ("#SBATCH --gres=gpu:1        # TODO: adjust GPU count/type for your cluster"
                if gpu else "# (CPU engine -- no --gres requested)")
        env = self.env_activate or (
            "# TODO: add env activation for your cluster, e.g.\n"
            "# module load cuda/12.2\n"
            "# eval \"$(/path/to/micromamba shell hook -s bash)\" && micromamba activate synthon-dock")

        (cx, cy, cz), (sx, sy, sz) = self.t.pocket
        repl = {
            "PARTITION": partition,
            "TIME": self.time_limit,
            "ARRAY": f"1-{n_tasks}",
            "GRES": gres,
            "ENV_ACTIVATE": env,
            "RUN_DIR": run_dir,
            "RECEPTOR": receptor,
            "CX": cx, "CY": cy, "CZ": cz,
            "SX": sx, "SY": sy, "SZ": sz,
            "EXH": int(self.t.exhaustiveness),
            "SCORING": self.t.scoring_function,
            "ENGINE_BIN": _ENGINE_BIN.get(self.engine, self.engine),
            "EXTRA_FLAGS": self._extra_flags(),
        }
        with open(os.path.join(_SLURM_DIR, _TEMPLATE[self.engine])) as fh:
            text = fh.read()
        for k, v in repl.items():
            text = text.replace("{{" + k + "}}", str(v))
        job = os.path.join(run_dir, "job.sbatch")
        with open(job, "w") as fh:
            fh.write(text)
        return job

    def _extra_flags(self) -> str:
        if self.engine == "smina" and self.t.scoring_function == "vinardo":
            return "--scoring vinardo"
        if self.engine == "gnina":
            return "--cnn_scoring rescore"
        return ""

    def _submit_and_wait(self, run_dir: str, n_tasks: int) -> str:
        sbatch = shutil.which("sbatch")
        if not sbatch:
            raise RuntimeError(
                "sbatch not found -- SlurmDockingOracle must run on a SLURM login "
                "node. (The sbatch + ligand files were written to "
                f"{run_dir}; you can submit them by hand there.)")
        proc = subprocess.run([sbatch, "job.sbatch"], cwd=run_dir,
                              capture_output=True, text=True)
        m = re.search(r"Submitted batch job (\d+)", proc.stdout)
        if not m:
            raise RuntimeError(f"sbatch failed: {proc.stdout.strip()} {proc.stderr.strip()}")
        job_id = m.group(1)
        self.log(f"[slurm] submitted array job {job_id} ({n_tasks} tasks)")
        self._wait(job_id)
        return job_id

    def _wait(self, job_id: str) -> None:
        squeue = shutil.which("squeue") or "squeue"
        while True:
            q = subprocess.run([squeue, "-j", str(job_id), "-h", "-o", "%T"],
                               capture_output=True, text=True)
            states = [s for s in q.stdout.split() if s]
            if not states:
                break
            self.log(f"[slurm] job {job_id}: {len(states)} task(s) active "
                     f"({','.join(sorted(set(states)))})")
            time.sleep(self.poll_interval)
        self.log(f"[slurm] job {job_id} complete")

    # -- collect -----------------------------------------------------------
    def _collect(self, run_dir: str, idx2smi) -> dict:
        """Parse the best score per ligand from the engine output; also cache the
        best pose file per SMILES for lean pose export."""
        outdir = os.path.join(run_dir, "out")
        out = {}
        for idx, smi in idx2smi.items():
            if self.engine == "gnina":
                f = os.path.join(outdir, f"out_{idx}.sdf")
                parser = parse_gnina_sdf
            elif self.engine == "unidock":
                f = os.path.join(outdir, f"lig_{idx}_out.pdbqt")
                parser = parse_vina_pdbqt
            else:
                f = os.path.join(outdir, f"out_{idx}.pdbqt")
                parser = parse_vina_pdbqt
            if not os.path.exists(f):
                continue
            try:
                with open(f) as fh:
                    sc = parser(fh.read())
            except Exception:  # noqa: BLE001
                continue
            if sc is None:
                continue
            out[smi] = sc
            prev = self._best_pose.get(smi)
            if prev is None or sc < prev[1]:
                self._best_pose[smi] = (f, sc)
        return out


# --------------------------------------------------------------------------- #
# Lean pose output: ONE gzipped multi-molecule SDF (not dozens of loose PDBs)
# --------------------------------------------------------------------------- #
def _pose_file_to_mol(path: str):
    """Read the best pose from an engine output file as an RDKit Mol.
    SDF (gnina) is read directly; PDBQT (vina-family) is converted with obabel
    (first/best pose only). Returns a Mol or None."""
    if not _HAVE_RDKIT:
        return None
    if path.endswith(".sdf"):
        supp = Chem.SDMolSupplier(path, sanitize=False)
        for m in supp:
            if m is not None:
                return m
        return None
    obabel = shutil.which("obabel")
    if not obabel:
        return None
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
    tmp.close()
    try:
        # -f 1 -l 1 = first (best) pose only
        subprocess.run([obabel, path, "-O", tmp.name, "-f", "1", "-l", "1"],
                       capture_output=True, text=True)
        supp = Chem.SDMolSupplier(tmp.name, sanitize=False)
        for m in supp:
            if m is not None:
                return m
    finally:
        os.unlink(tmp.name)
    return None


def write_top_hits_sdf_gz(hits, oracle: SlurmDockingOracle, out_path: str,
                          log=print):
    """Write the top hits' docked poses as ONE gzipped multi-mol SDF.

    `hits` is a best-first list of ``(smiles, score)``. Each written record
    carries ``Score``, ``SMILES`` and ``Rank`` as SDF properties. Uses the pose
    files the oracle cached during scoring (no re-dock). Deliberately ONE
    ``top_hits.sdf.gz`` -- not dozens of loose PDBs -- to keep the repo lean."""
    if not _HAVE_RDKIT:
        log("[slurm] RDKit unavailable; skipping pose export")
        return None
    best = getattr(oracle, "_best_pose", {})
    mols = []
    for rank, (smi, score) in enumerate(hits, 1):
        entry = best.get(smi)
        if not entry:
            continue
        m = _pose_file_to_mol(entry[0])
        if m is None:
            continue
        m.SetProp("_Name", f"rank{rank}")
        m.SetProp("Score", f"{score:.3f}")
        m.SetProp("SMILES", smi)
        m.SetProp("Rank", str(rank))
        mols.append(m)
    if not mols:
        log("[slurm] no cached poses to export")
        return None
    import gzip
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
    tmp.close()
    try:
        w = Chem.SDWriter(tmp.name)
        for m in mols:
            w.write(m)
        w.close()
        with open(tmp.name, "rb") as fin, gzip.open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    finally:
        os.unlink(tmp.name)
    log(f"[slurm] wrote {len(mols)} poses -> {out_path}")
    return out_path
