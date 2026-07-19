#!/usr/bin/env python
"""Muni-backed synthon-TS screen — retrieval via onepot direct API, docking via
the muni CLI (`rowan_batch_docking`, QVina2). Runs on the shared muni workspace
and bills muni credits (~0.05 cr/ligand). No Rowan key needed.

    export ONEPOT_API_KEY=...        # onepot direct key (muni onepot is broken)
    muni page use "<your trial page>"      # so jobs land on the right page
    python -m synthon_ts.run_muni --query "<SMILES>" --max-docks 200

Uses the CDK2/1HCK ATP-site box + prepared protein UUID by default (the
calibration target). Override with --protein / --pocket for another target.
"""
from __future__ import annotations

import argparse
import json
import os

from .core import LoopConfig, Target, diverse_top, run_loop
from .oracle_muni import MuniBatchDockingOracle

# CDK2 / PDB 1HCK ATP-site box (Rowan Pocketeer, rank-1 pocket) + the prepared
# protein UUID from muni's pocket-detection job (thread it so docking uses the
# same prepared structure the box was computed against).
DEFAULT_POCKET = [[98.49406480789185, 94.42433500289917, 98.22893345355988],
                  [22.622045516967773, 24.517674446105957, 19.841248273849487]]
DEFAULT_PROTEIN_UUID = "f4901ad5-9715-4a6b-b9b9-f48165cb5b6d"


def main():
    ap = argparse.ArgumentParser(description="Muni-backed synthon-TS screen")
    ap.add_argument("--query", required=True, help="seed query SMILES")
    ap.add_argument("--protein", default=DEFAULT_PROTEIN_UUID,
                    help="muni/Rowan prepared-protein UUID, or a PDB id")
    ap.add_argument("--pocket", default=None, help="JSON [[cx,cy,cz],[sx,sy,sz]]; default CDK2 box")
    ap.add_argument("--executable", default="qvina2", choices=["qvina2", "vina"])
    ap.add_argument("--scoring", default="vina", choices=["vina", "vinardo", "ad4"])
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--max-docks", type=int, default=200)
    ap.add_argument("--max-rounds", type=int, default=6)
    ap.add_argument("--seed-hits", type=int, default=60)
    ap.add_argument("--round-hits", type=int, default=40)
    ap.add_argument("--anchors", type=int, default=3)
    ap.add_argument("--coarse", action="store_true",
                    help="coarse anchor-based attribution (no per-analog decompose; "
                         "cheaper onepot). Default is precise per-product attribution.")
    ap.add_argument("--page-id", default=None, help="muni page id for trial jobs")
    ap.add_argument("--out", default="synthon_ts_muni_hits.json")
    args = ap.parse_args()

    pocket = json.loads(args.pocket) if args.pocket else DEFAULT_POCKET
    target = Target(protein=args.protein, pocket=pocket, executable=args.executable,
                    scoring_function=args.scoring, exhaustiveness=args.exhaustiveness)
    cfg = LoopConfig(query_smiles=args.query, seed_max_results=args.seed_hits,
                     round_max_results=args.round_hits, n_anchors_per_round=args.anchors,
                     max_rounds=args.max_rounds, max_docks=args.max_docks,
                     precise_attribution=not args.coarse)
    oracle = MuniBatchDockingOracle(target, page_id=args.page_id,
                                    name="synthon-TS dock")

    res = run_loop(oracle, target, cfg, onepot_key=os.environ["ONEPOT_API_KEY"])
    _report(res, cfg, args.out)


def _report(res, cfg, out):
    ranked = res["ranked"]
    print(f"\n=== synthon-TS (muni) done: {res['n_docks']} docks, "
          f"{len(res['pool'])} scored ===")
    print("Top scaffold-diverse hits (lower score = better):")
    for p in diverse_top(ranked, cfg.top_k):
        print(f"  {p.score:+6.2f}  {p.smiles}")
    json.dump(res["top"], open(out, "w"), indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
