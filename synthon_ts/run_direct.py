#!/usr/bin/env python
"""Direct-API synthon-TS screen — onepot + Rowan SDK only, NO muni.

Fully portable: drop this repo on any machine with Python 3.11+, `pip install
onepot rowan-python rdkit`, set two env vars, and run.

    export ONEPOT_API_KEY=...        # onepot direct key
    export ROWAN_API_KEY=...         # your Rowan account key
    python -m synthon_ts.run_direct --query "<SMILES>" --max-docks 120

Docking oracle = Rowan single `docking` (one workflow per ligand, fired in
concurrent waves). NB: Rowan `batch_docking` is feature-gated on this account
(400 "you do not have access"), and `analogue_docking` needs a bound template
pose, so single docking is the backend. The target PDB is prepared ONCE into a
Rowan protein UUID up front, then reused for every dock. See oracle_rowan.py.
"""
from __future__ import annotations

import argparse
import json
import os

from .core import LoopConfig, Target, diverse_top, run_loop
from .oracle_rowan import RowanDockingOracle, prepare_protein_uuid

# CDK2 / PDB 1HCK ATP-site box (Rowan Pocketeer, rank-1 druggability pocket).
DEFAULT_POCKET = [[98.49406480789185, 94.42433500289917, 98.22893345355988],
                  [22.622045516967773, 24.517674446105957, 19.841248273849487]]


def main():
    ap = argparse.ArgumentParser(description="Direct-API synthon-TS screen (onepot + Rowan)")
    ap.add_argument("--query", required=True, help="seed query SMILES")
    ap.add_argument("--protein", default="1HCK", help="PDB id (auto-prepared) or Rowan UUID")
    ap.add_argument("--pocket", default=None, help="JSON [[cx,cy,cz],[sx,sy,sz]]; default CDK2 box")
    ap.add_argument("--scoring", default="vina", choices=["vina", "vinardo"])
    ap.add_argument("--executable", default="qvina2", choices=["qvina2", "qvina-w", "vina"])
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--max-poses", type=int, default=4, help="poses kept per ligand")
    ap.add_argument("--csearch", action="store_true",
                    help="run Rowan conformer search per ligand (~3x slower/costlier, "
                         "marginally better scores; default off — vina samples poses itself)")
    ap.add_argument("--max-inflight", type=int, default=8,
                    help="ligands submitted concurrently per wave")
    ap.add_argument("--max-docks", type=int, default=120)
    ap.add_argument("--max-rounds", type=int, default=6)
    ap.add_argument("--seed-hits", type=int, default=60)
    ap.add_argument("--round-hits", type=int, default=30)
    ap.add_argument("--anchors", type=int, default=3)
    ap.add_argument("--max-credits", type=int, default=200, help="Rowan credit cap per single dock")
    ap.add_argument("--out", default="synthon_ts_direct_hits.json")
    args = ap.parse_args()

    api_key = os.environ["ROWAN_API_KEY"]
    pocket = json.loads(args.pocket) if args.pocket else DEFAULT_POCKET

    # Prepare the protein exactly once (PDB id -> prepared UUID); reuse everywhere.
    protein_uuid = prepare_protein_uuid(args.protein, api_key=api_key)

    target = Target(protein=protein_uuid, pocket=pocket, executable=args.executable,
                    scoring_function=args.scoring, exhaustiveness=args.exhaustiveness)
    cfg = LoopConfig(query_smiles=args.query, seed_max_results=args.seed_hits,
                     round_max_results=args.round_hits, n_anchors_per_round=args.anchors,
                     max_rounds=args.max_rounds, max_docks=args.max_docks)
    oracle = RowanDockingOracle(target, api_key=api_key, protein_uuid=protein_uuid,
                                max_poses=args.max_poses, max_credits=args.max_credits,
                                max_inflight=args.max_inflight,
                                do_csearch=args.csearch)

    res = run_loop(oracle, target, cfg, onepot_key=os.environ["ONEPOT_API_KEY"])
    _report(res, cfg, args.out)


def _report(res, cfg, out):
    ranked = res["ranked"]
    print(f"\n=== synthon-TS (direct) done: {res['n_docks']} docks, "
          f"{len(res['pool'])} scored ===")
    print("Top scaffold-diverse hits (lower score = better):")
    for p in diverse_top(ranked, cfg.top_k):
        print(f"  {p.score:+6.2f}  {p.smiles}")
    json.dump(res["top"], open(out, "w"), indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
