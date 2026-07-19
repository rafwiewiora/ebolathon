"""Docking oracle backend #2 — direct Rowan SDK (no muni).

Uses the user's own Rowan account + credits via the single ``docking`` workflow
(``rowan.submit_docking_workflow``), one Rowan workflow per ligand, fired off in
concurrent waves so Rowan runs them in parallel.

IMPORTANT ACCESS NOTE (verified live 2026-07-19 on this account)
----------------------------------------------------------------
* ``batch_docking`` is feature-gated: ``submit_batch_docking_workflow`` returns
  ``400 You do not have access to this feature``. Do NOT use it here.
* ``analogue_docking`` needs a bound template pose and FAILS on arbitrary
  anchors — not usable for the synthon loop.
* Single ``docking`` IS accessible and is the backend used below. It needs a
  target the account can resolve. A **prepared protein UUID** always works; a
  bare PDB id ("1HCK") may or may not resolve at runtime depending on the
  account/server state, so we prepare a protein once up front and dock against
  its UUID (see ``prepare_protein_uuid``). Preparation is done with
  ``rowan.create_protein_from_pdb_id(...)`` + ``Protein.prepare()`` (the
  ``pocket_detection`` workflow was observed to fail, so we do NOT rely on it).

Score schema: a finished single dock exposes ``wf.model_dump()["data"]["scores"]``,
a list of per-pose records (``DockingScore``: ``score``, ``pose``, ``strain``,
``rmsd``, ``mmgbsa_score``, ...). ``score`` is the Vina affinity in kcal/mol
(lower = better). We take the best (lowest) ``score`` across poses per ligand. A
failed dock has an empty ``scores`` list and is omitted from the result.

Auth: ``rowan.api_key`` from env ROWAN_API_KEY (never hardcode).
"""
from __future__ import annotations

import os
import re
import time

import rowan

# VinaSettings lives in stjames; Rowan re-exports it as rowan.VinaSettings too.
try:
    from stjames.workflows.docking import VinaSettings
except Exception:  # pragma: no cover
    from rowan import VinaSettings  # type: ignore

from .core import Target

_UUID_RE = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.I
)


def _looks_like_uuid(s: str) -> bool:
    return bool(_UUID_RE.match(s.strip()))


# --------------------------------------------------------------------------- #
# Protein preparation
# --------------------------------------------------------------------------- #
def prepare_protein_uuid(pdb_or_uuid: str, api_key: str | None = None,
                         max_credits: int = 200, log=print) -> str:
    """Resolve a docking target to a prepared-protein UUID that single ``docking``
    can use on THIS account.

    * If ``pdb_or_uuid`` already looks like a UUID, it is returned unchanged.
    * Otherwise it is treated as a PDB id: a protein is created from the PDB id
      and prepared (sanitize + add H etc.), and the prepared UUID is returned.

    ``max_credits`` is currently informational (protein prep is cheap, ~2 cr).
    """
    rowan.api_key = api_key or os.environ["ROWAN_API_KEY"]
    if _looks_like_uuid(pdb_or_uuid):
        log(f"[rowan-oracle] using protein UUID as-is: {pdb_or_uuid}")
        return pdb_or_uuid
    log(f"[rowan-oracle] preparing protein from PDB id {pdb_or_uuid!r} ...")
    protein = rowan.create_protein_from_pdb_id(
        pdb_or_uuid, name=f"{pdb_or_uuid}-synthon-ts")
    protein.prepare()  # blocks with internal polling; sets sanitized state
    log(f"[rowan-oracle] prepared protein {pdb_or_uuid} -> uuid={protein.uuid} "
        f"(sanitized={getattr(protein, 'sanitized', '?')})")
    return protein.uuid


# --------------------------------------------------------------------------- #
# Oracle
# --------------------------------------------------------------------------- #
class RowanDockingOracle:
    """Single-``docking`` oracle: one Rowan workflow per ligand, submitted in
    concurrent waves so Rowan parallelises them.

    ``score(smiles_list, template_smiles=None) -> {smiles: best_score}`` where
    lower is better (Vina). Ligands whose dock fails are omitted.

    Cost/quality (measured live, 1HCK/CDK2, qvina2/vina/exh=8/max_poses=4):
      * ``do_csearch=False`` (default): ~2 credits, ~2 min/ligand. Roscovitine
        best score -7.5. Vina samples poses itself from the embedded conformer.
      * ``do_csearch=True``: ~6.4 credits, ~6.4 min/ligand. Roscovitine -7.7.
        Marginally better; rarely worth ~3x cost for a 100s-of-docks loop.
    """

    def __init__(self, target: Target, api_key: str | None = None,
                 protein_uuid: str | None = None, max_poses: int = 4,
                 max_credits: int = 200, poll_interval: int = 5,
                 max_inflight: int = 8, do_csearch: bool = False,
                 name: str = "synthon-TS dock", log=print):
        rowan.api_key = api_key or os.environ["ROWAN_API_KEY"]
        self.t = target
        self.max_poses = max_poses
        self.max_credits = max_credits
        self.poll = poll_interval
        self.max_inflight = max_inflight
        self.do_csearch = do_csearch
        self.name = name
        self.log = log
        # Resolve the target to a prepared-protein UUID exactly once.
        self.protein_uuid = protein_uuid or prepare_protein_uuid(
            target.protein, api_key=rowan.api_key, max_credits=max_credits, log=log)

    # -- settings ----------------------------------------------------------
    def _vina_settings(self) -> "VinaSettings":
        return VinaSettings(
            executable=self.t.executable,            # 'qvina2'
            scoring_function=self.t.scoring_function,  # 'vina' | 'vinardo'
            exhaustiveness=self.t.exhaustiveness,
            max_poses=self.max_poses,
        )

    # -- public API --------------------------------------------------------
    def score(self, smiles_list, template_smiles=None) -> dict:
        # template_smiles is unused for single docking (each ligand is docked
        # independently into the pocket box); accepted for contract parity.
        out: dict = {}
        smiles_list = list(smiles_list)
        for i in range(0, len(smiles_list), self.max_inflight):
            wave = smiles_list[i:i + self.max_inflight]
            out.update(self._dock_wave(wave, i))
        return out

    # -- one concurrent wave ----------------------------------------------
    def _dock_wave(self, wave, offset) -> dict:
        settings = self._vina_settings()
        submitted = []  # (smiles, workflow) or (smiles, None) on submit failure
        for smi in wave:
            try:
                mol = rowan.Molecule.from_smiles(smi)
                wf = rowan.submit_docking_workflow(
                    protein=self.protein_uuid,
                    pocket=self.t.pocket,
                    initial_molecule=mol,
                    docking_settings=settings,
                    do_csearch=self.do_csearch,
                    name=f"{self.name} #{offset + len(submitted)}",
                    max_credits=self.max_credits,
                )
                submitted.append((smi, wf))
            except Exception as e:
                self.log(f"[rowan-oracle] submit failed for {smi[:40]}: "
                         f"{str(e)[:140]}")
                submitted.append((smi, None))
        self.log(f"[rowan-oracle] wave submitted: "
                 f"{sum(1 for _, w in submitted if w)}/{len(wave)} ligands")

        out = {}
        for smi, wf in submitted:
            if wf is None:
                continue
            try:
                wf.wait_for_result(poll_interval=self.poll)
                wf.fetch_latest(in_place=True)
                sc = _best_score(wf.model_dump())
            except Exception as e:
                self.log(f"[rowan-oracle] dock failed for {smi[:40]}: "
                         f"{str(e)[:140]}")
                sc = None
            if sc is not None:
                out[smi] = sc
                self.log(f"[rowan-oracle]   {sc:+.2f}  {smi[:60]}")
            else:
                self.log(f"[rowan-oracle]   (no score) {smi[:60]}")
        return out


# --------------------------------------------------------------------------- #
# Score parsing
# --------------------------------------------------------------------------- #
def _best_score(dump: dict):
    """Best (lowest) Vina score from a single-``docking`` model_dump.

    Expected shape: ``dump["data"]["scores"]`` is a list of per-pose records,
    each a dict with a numeric ``score`` field (Vina affinity, kcal/mol; lower
    better). Returns the min score, or ``None`` if the dock failed / produced no
    scores. Parsing is defensive across minor schema shifts."""
    data = dump.get("data")
    if not isinstance(data, dict):
        return None
    scores = data.get("scores")
    if not isinstance(scores, list) or not scores:
        return None
    vals = []
    for rec in scores:
        v = None
        if isinstance(rec, dict):
            for k in ("score", "docking_score", "affinity", "best_score"):
                if isinstance(rec.get(k), (int, float)):
                    v = float(rec[k])
                    break
        elif isinstance(rec, (int, float)):
            v = float(rec)
        if v is not None:
            vals.append(v)
    return min(vals) if vals else None
