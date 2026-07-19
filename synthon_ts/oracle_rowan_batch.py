"""Docking oracle backend #3 — direct Rowan ``batch_docking`` (QVina2).

High-throughput bulk scoring on your OWN Rowan account: one workflow docks many
ligands, ~0.1 cr/ligand (vs ~1.5-2 cr for single docking), dozens per job in
seconds. Verified working with a PREPARED protein UUID (2026-07-19). This is the
path for screening THOUSANDS of molecules; single docking (oracle_rowan.py) is
~15-20x slower/pricier and only needed when you want poses per ligand.

IMPORTANT: like all batch docking, this is **scores-only — no poses**. Pose
export re-docks the top-K through single docking (poses.export_top_poses handles
that automatically: the batch oracle exposes no pose_cache, so export falls back
to a fresh dock for just the winners).

Auth: ``rowan.api_key`` from env ROWAN_API_KEY. Needs a prepared protein UUID
(bare PDB id fails ``get_protein 400``); ``prepare_protein_uuid`` handles it.
"""
from __future__ import annotations

import os

import rowan

from .core import Target
from .oracle_rowan import prepare_protein_uuid


class RowanBatchDockingOracle:
    """``score(smiles_list) -> {smiles: best_score}`` via Rowan ``batch_docking``.

    Ligands are split into chunks of ``chunk`` and submitted in concurrent waves
    of ``max_inflight`` jobs (submit the whole wave, then collect) so Rowan runs
    them in parallel. Lower score = better; failed ligands are omitted."""

    def __init__(self, target: Target, api_key: str | None = None,
                 protein_uuid: str | None = None, chunk: int = 60,
                 max_inflight: int = 6, max_credits: int = 2000,
                 poll_interval: int = 5, name: str = "synthon-TS batch dock",
                 log=print):
        rowan.api_key = api_key or os.environ["ROWAN_API_KEY"]
        self.t = target
        self.chunk = chunk
        self.max_inflight = max_inflight
        self.max_credits = max_credits
        self.poll = poll_interval
        self.name = name
        self.log = log
        # No pose_cache: batch docking returns no poses (export re-docks top-K).
        self.protein_uuid = protein_uuid or prepare_protein_uuid(
            target.protein, api_key=rowan.api_key, log=log)

    def score(self, smiles_list, template_smiles=None) -> dict:
        smiles_list = list(smiles_list)
        chunks = [smiles_list[i:i + self.chunk]
                  for i in range(0, len(smiles_list), self.chunk)]
        out: dict = {}
        for w in range(0, len(chunks), self.max_inflight):
            wave = chunks[w:w + self.max_inflight]
            submitted = []
            for ch in wave:
                try:
                    wf = rowan.submit_batch_docking_workflow(
                        smiles_list=ch, protein=self.protein_uuid, pocket=self.t.pocket,
                        executable=self.t.executable,
                        scoring_function=self.t.scoring_function,
                        exhaustiveness=self.t.exhaustiveness,
                        name=self.name, max_credits=self.max_credits)
                    submitted.append((ch, wf))
                except Exception as e:  # noqa: BLE001
                    self.log(f"[rowan-batch] submit failed ({len(ch)} ligs): {str(e)[:140]}")
            self.log(f"[rowan-batch] wave: {len(submitted)}/{len(wave)} jobs "
                     f"({sum(len(c) for c, _ in submitted)} ligands)")
            for ch, wf in submitted:
                try:
                    wf.wait_for_result(poll_interval=self.poll)
                    wf.fetch_latest(in_place=True)
                    data = wf.model_dump().get("data") or {}
                    scores = data.get("best_scores") or []
                    sent = data.get("initial_smiles_list") or ch
                    n = 0
                    for smi, sc in zip(sent, scores):
                        if isinstance(sc, (int, float)):
                            out[smi] = float(sc)
                            n += 1
                    self.log(f"[rowan-batch]   +{n} scored (total {len(out)})")
                except Exception as e:  # noqa: BLE001
                    self.log(f"[rowan-batch] job failed ({len(ch)} ligs): {str(e)[:140]}")
        return out
