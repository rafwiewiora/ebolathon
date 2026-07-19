"""Export docked poses of the top hits to disk for live inspection in PyMOL.

The screen's oracle only returns *scores* (and the muni backend never produces
poses at all). To actually LOOK at how the winners bind, we re-dock the top hits
through the direct Rowan single-``docking`` workflow, pull the best pose geometry
back, and write:

    <out_dir>/receptor.pdb        prepared receptor (once, grey cartoon)
    <out_dir>/rank{N}_{score}.pdb one docked ligand pose per hit (sticks)
    <out_dir>/top_hits.sdf        best-effort combined SDF (RDKit bond perception)
    <out_dir>/hits.csv            rank, smiles, score, mmgbsa, posebusters_valid
    <out_dir>/view.pml            PyMOL script — load + colour + label + zoom

Pose geometry comes back from ``retrieve_calculation_molecules(pose_uuid)`` as a
list of one molecule dict whose ``atoms`` are ``{atomic_number, position:[x,y,z],
mass}`` in the RECEPTOR coordinate frame (verified live). We convert those atoms
straight to PDB ``HETATM`` records; PyMOL infers bonds by distance, so no bond
list is needed.

Requires ``ROWAN_API_KEY`` regardless of the screen's scoring backend (muni batch
docking has no poses), and a prepared-protein UUID in ``target.protein``.
"""
from __future__ import annotations

import csv
import os
import shutil
import tempfile

import rowan

try:
    from stjames.workflows.docking import VinaSettings
except Exception:  # pragma: no cover
    from rowan import VinaSettings  # type: ignore

try:
    from stjames import ELEMENT_SYMBOL  # {atomic_number: "C", ...}
except Exception:  # pragma: no cover
    ELEMENT_SYMBOL = {}

from .core import Target


# --------------------------------------------------------------------------- #
# atoms -> PDB
# --------------------------------------------------------------------------- #
def _symbol(z: int) -> str:
    return ELEMENT_SYMBOL.get(z, "C") if ELEMENT_SYMBOL else "C"


def atoms_to_pdb(atoms: list, resname: str = "LIG", chain: str = "X") -> str:
    """Render a list of ``{atomic_number, position}`` atoms as PDB HETATM records.

    Element symbols go in cols 77-78 (the field PyMOL trusts); bonds are left to
    PyMOL's distance-based perception. A per-element serial makes atom names
    unique (``C1``, ``C2``, ``N1`` ...)."""
    lines = []
    per_elem: dict = {}
    for i, a in enumerate(atoms, 1):
        z = a.get("atomic_number", 6)
        x, y, z_ = a["position"]
        sym = _symbol(z)
        per_elem[sym] = per_elem.get(sym, 0) + 1
        name = f"{sym}{per_elem[sym]}"[:4]
        # strict PDB columns: name 13-16, altLoc 17, resName 18-20, chain 22,
        # resSeq 23-26, xyz 31-54, occ/temp 55-66, element 77-78.
        lines.append(
            "HETATM%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
            % (i, name, "", resname[:3], chain[:1], 1, "", x, y, z_, 1.00, 0.00, sym)
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def _atoms_to_rdkit_block(atoms: list, name: str):
    """Best-effort RDKit mol (with perceived bonds) from pose atoms, for SDF.
    Returns an RDKit Mol or None if RDKit / bond perception is unavailable."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Geometry import Point3D
    except Exception:
        return None
    try:
        rw = Chem.RWMol()
        conf = Chem.Conformer(len(atoms))
        for i, a in enumerate(atoms):
            idx = rw.AddAtom(Chem.Atom(int(a.get("atomic_number", 6))))
            x, y, z = a["position"]
            conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
        m = rw.GetMol()
        m.AddConformer(conf, assignId=True)
        rdDetermineBonds.DetermineBonds(m, charge=0)  # perception only, no FF opt
        m.SetProp("_Name", name)
        return m
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# Rowan helpers
# --------------------------------------------------------------------------- #
def _download_receptor(protein_uuid: str, dest_pdb: str, log=print) -> None:
    """Download the prepared receptor PDB to `dest_pdb`. ``download_pdb_file``
    writes ``<protein-name>.pdb`` INTO a directory, so we stage in a temp dir and
    move the single file to the requested path."""
    prot = rowan.retrieve_protein(protein_uuid)
    with tempfile.TemporaryDirectory() as td:
        prot.download_pdb_file(path=td)
        pdbs = [f for f in os.listdir(td) if f.lower().endswith(".pdb")]
        if not pdbs:
            raise RuntimeError(f"no PDB written for protein {protein_uuid}")
        shutil.move(os.path.join(td, pdbs[0]), dest_pdb)
    log(f"[poses] wrote receptor {dest_pdb}")


def _dock_one_pose(smiles: str, target: Target, max_poses: int, log=print):
    """Dock `smiles` once; return the best pose as
    ``(atoms, score, mmgbsa, posebusters_valid)`` or None on failure."""
    settings = VinaSettings(executable=target.executable,
                            scoring_function=target.scoring_function,
                            exhaustiveness=target.exhaustiveness, max_poses=max_poses)
    wf = rowan.submit_docking_workflow(
        protein=target.protein, pocket=target.pocket,
        initial_molecule=rowan.Molecule.from_smiles(smiles),
        docking_settings=settings, name="synthon-TS pose export", max_credits=200)
    wf.wait_for_result(poll_interval=5)
    wf.fetch_latest(in_place=True)
    scores = (wf.model_dump().get("data") or {}).get("scores") or []
    best = None
    for rec in scores:
        s = rec.get("score")
        if isinstance(s, (int, float)) and (best is None or s < best["score"]):
            best = rec
    if best is None:
        return None
    mols = rowan.retrieve_calculation_molecules(best["pose"])
    if not mols:
        return None
    return (mols[0]["atoms"], float(best["score"]),
            best.get("mmgbsa_score"), best.get("posebusters_valid"))


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
def export_top_poses(hits, target: Target, out_dir: str, api_key: str,
                     top_k: int = 10, log=print, cache: dict | None = None) -> dict:
    """Re-dock the top hits and write receptor + poses + view.pml into `out_dir`.

    Parameters
    ----------
    hits : ranked list of ``(smiles, score)`` (best first). Only the score from
        the pose re-dock is written; the incoming score is used only for ordering.
    target : a ``Target`` whose ``protein`` is a *prepared Rowan UUID* and whose
        ``pocket`` is the docking box. (For the muni backend, prepare a Rowan UUID
        separately — muni batch docking produces no poses.)
    out_dir : output directory (created if needed).
    api_key : Rowan API key (required; poses always go through direct Rowan).
    cache : optional ``{smiles: pose_tuple}`` dict reused across calls so live
        per-round exports never re-dock a hit already posed.

    Returns a summary dict; the key output for the user is ``out_dir/view.pml``.
    """
    os.makedirs(out_dir, exist_ok=True)
    rowan.api_key = api_key or os.environ["ROWAN_API_KEY"]
    cache = cache if cache is not None else {}

    receptor_path = os.path.join(out_dir, "receptor.pdb")
    if not os.path.exists(receptor_path):
        _download_receptor(target.protein, receptor_path, log=log)

    hits = list(hits)[:top_k]
    rows = []          # (rank, smiles, dock_score, mmgbsa, pb_valid, pdb_filename)
    rdkit_mols = []
    for rank, (smi, _incoming) in enumerate(hits, 1):
        pose = cache.get(smi)
        if pose is None:
            log(f"[poses] docking rank {rank}: {smi[:60]}")
            try:
                pose = _dock_one_pose(smi, target, max_poses=4, log=log)
            except Exception as e:  # noqa: BLE001
                log(f"[poses] rank {rank} dock failed: {str(e)[:140]}")
                pose = None
            if pose is not None:
                cache[smi] = pose
        if pose is None:
            continue
        atoms, dock_score, mmgbsa, pb_valid = pose
        fname = f"rank{rank}_{dock_score:+.2f}.pdb"
        with open(os.path.join(out_dir, fname), "w") as fh:
            fh.write(atoms_to_pdb(atoms, resname="LIG"))
        rows.append((rank, smi, dock_score, mmgbsa, pb_valid, fname))
        m = _atoms_to_rdkit_block(atoms, name=f"rank{rank}")
        if m is not None:
            m.SetProp("rank", str(rank)); m.SetProp("score", f"{dock_score:.2f}")
            m.SetProp("SMILES", smi)
            rdkit_mols.append(m)
        log(f"[poses]   rank {rank}: dock={dock_score:+.2f} pb={pb_valid} -> {fname}")

    _write_csv(os.path.join(out_dir, "hits.csv"), rows)
    _write_sdf(os.path.join(out_dir, "top_hits.sdf"), rdkit_mols, log=log)
    _write_pml(out_dir, rows, target.pocket, log=log)
    log(f"[poses] exported {len(rows)} poses -> {out_dir}  (open {out_dir}/view.pml)")
    return {"out_dir": out_dir, "n_poses": len(rows),
            "view_pml": os.path.join(out_dir, "view.pml"),
            "receptor": receptor_path}


def _write_csv(path: str, rows) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "smiles", "score", "mmgbsa", "posebusters_valid", "pose_file"])
        for rank, smi, score, mmgbsa, pb, fname in rows:
            w.writerow([rank, smi, f"{score:.3f}",
                        "" if mmgbsa is None else mmgbsa,
                        "" if pb is None else pb, fname])


def _write_sdf(path: str, mols, log=print) -> None:
    if not mols:
        # leave a stub so downstream tooling / the user isn't surprised by absence
        open(path, "w").close()
        return
    try:
        from rdkit import Chem
        w = Chem.SDWriter(path)
        for m in mols:
            w.write(m)
        w.close()
    except Exception as e:  # noqa: BLE001
        log(f"[poses] SDF write skipped: {str(e)[:120]}")


def _write_pml(out_dir: str, rows, pocket, log=print) -> None:
    """Write a robust, re-loadable PyMOL script. Uses ``cd`` to the out dir so the
    relative loads work no matter where PyMOL was launched, and ``reinitialize``
    so re-running it (live refresh) never stacks duplicate objects."""
    abs_dir = os.path.abspath(out_dir)
    (cx, cy, cz), _ = pocket
    lines = [
        "# ================================================================",
        "# synthon-TS docked poses - live view",
        "#",
        "# LIVE REFRESH: while a screen is running it rewrites these files",
        "# after every round. To pull in the newest poses, in the PyMOL",
        "# command line just re-run this script:",
        "#     PyMOL>  @view.pml",
        "# (or use the menu: File > Reload All). 'reinitialize' below keeps",
        "# re-loads clean - no duplicated/stacked objects.",
        "# ================================================================",
        "reinitialize",
        "bg_color white",
        f"cd {abs_dir}",
        "",
        "# --- receptor ---",
        "load receptor.pdb, receptor",
        "hide everything, receptor",
        "show cartoon, receptor",
        "color grey70, receptor",
        "set cartoon_transparency, 0.15, receptor",
        "",
        "# --- docked ligand poses (best score first) ---",
    ]
    pose_objs = []
    for rank, smi, score, mmgbsa, pb, fname in rows:
        obj = f"pose_rank{rank}"
        pose_objs.append(obj)
        lines += [
            f"load {fname}, {obj}",
            f"show sticks, {obj}",
            f"util.cbag('{obj}')",   # colour by element, green carbons
            f'label first {obj}, "{rank}: {score:+.2f}"',
        ]
    if pose_objs:
        grp = " ".join(pose_objs)
        lines += [
            "",
            "# --- pocket context + framing ---",
            f"group poses, {grp}",
            f"select pocket_res, byres (receptor within 5 of (poses))",
            "show lines, pocket_res",
            "set label_size, 18",
            "set label_color, black",
            f"zoom (poses), 6",
        ]
    else:
        # no poses yet (e.g. very first round) - at least frame the box centre
        lines += [
            "pseudoatom box_center, pos=[%.3f, %.3f, %.3f]" % (cx, cy, cz),
            "zoom box_center, 15",
        ]
    lines += ["", "set ray_shadows, 0", "# end view.pml", ""]
    with open(os.path.join(out_dir, "view.pml"), "w") as fh:
        fh.write("\n".join(lines))
    log(f"[poses] wrote {out_dir}/view.pml")
