#!/usr/bin/env python
"""One-command synthon-TS screen: target PDB + a pocket ligand -> ranked hits +
docked poses you can open live in PyMOL.

    python -m synthon_ts.run \
        --pdb 1HCK --pocket-ligand ATP \
        --query "CCC(CO)Nc1nc(NCc2ccccc2)c2ncn(C(C)C)c2n1" \
        --backend direct --out-dir runs/cdk2 --pymol

What it does, end to end:
  1. Derives the docking box from the named pocket ligand in the PDB
     (`pocket.box_from_ligand`) - no manual box needed.
  2. Prepares the protein (direct: Rowan `create_protein_from_pdb_id` + prepare
     -> UUID; muni: passes the PDB id / UUID + box).
  3. Seeds the loop - either from `--query` (a known binder) OR, when `--query`
     is omitted, by sampling drug-like, CORE-native molecules from onepot
     `sample_space` (great when you have NO known binder).
  4. Runs the synthon-aware Thompson-Sampling screen (`core.run_loop`), dropping
     off-profile analogs with a local RDKit drug-like filter before docking.
  5. Exports docked poses of the current top hits after every round AND at the
     end (`poses.export_top_poses`) so `view.pml` refreshes live in PyMOL.

Keys via env only: ROWAN_API_KEY (always needed for pose export; also for the
direct docking backend) and ONEPOT_API_KEY (retrieval + sample_space).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from collections import defaultdict

from onepot import Client

from .core import LoopConfig, Target, diverse_top, run_loop
from .filters import DrugLikeLimits, make_filter
from .pocket import box_from_ligand
from .poses import export_top_poses

PYMOL_BIN = "/Users/bb/.local/share/mamba/envs/cadd-pymol/bin/pymol"


# --------------------------------------------------------------------------- #
# Protein prep
# --------------------------------------------------------------------------- #
def _resolve_rowan_uuid(pdb: str, api_key: str, log=print) -> str:
    """Resolve `pdb` (PDB id, local .pdb file, or Rowan UUID) to a prepared Rowan
    protein UUID (needed for pose export and the direct docking backend)."""
    import rowan
    from .oracle_rowan import _looks_like_uuid, prepare_protein_uuid
    rowan.api_key = api_key
    if _looks_like_uuid(pdb):
        return pdb
    if os.path.exists(pdb):
        log(f"[run] uploading local protein file {pdb} to Rowan ...")
        prot = rowan.upload_protein(name=os.path.basename(pdb), file_path=pdb)
        prot.prepare()
        log(f"[run] prepared uploaded protein -> {prot.uuid}")
        return prot.uuid
    return prepare_protein_uuid(pdb, api_key=api_key, log=log)


# --------------------------------------------------------------------------- #
# sample_space seeding
# --------------------------------------------------------------------------- #
def _sample_seeds(client: Client, count: int, strategy: str, limits: DrugLikeLimits,
                  mw_min: float, qed_min: float, seed: int | None,
                  max_price: int | None, log=print):
    """Draw `count` drug-like, CORE-native seed molecules from onepot sample_space.

    The property window is native to sample_space (NOTE: onepot's `clogp` is
    Crippen logP, used here as a proxy for the cLogD 1-5 target). sample_space has
    NO price/supplier param, but each molecule carries `price_usd`, so we
    over-sample and apply the price cut LOCALLY (supplier risk isn't returned at
    draw time; downstream analog searches enforce it natively). Returns
    (list[smiles], used_seed) - log `used_seed` to replay the exact draw."""
    props = {
        "molecular_weight": {"min": mw_min, "max": limits.mw_max},
        "clogp": {"min": limits.logp_min, "max": limits.logp_max},
        "tpsa": {"max": limits.tpsa_max},
        "hbd": {"max": limits.hbd_max},
        "hba": {"max": limits.hba_max},
        "rotatable_bonds": {"max": limits.rotatable_max},
        "qed": {"min": qed_min},
    }
    # over-sample so enough survive the local price cut ($125/$295 two-tier)
    draw = count * 4 if max_price is not None else count
    resp = client.sample_space(count=draw, strategy=strategy, properties=props,
                               seed=seed, include_properties=True)
    used_seed = resp["seed"]
    mols = resp["molecules"]
    kept = mols
    if max_price is not None:
        kept = [m for m in mols if (m.get("price_usd") is None or m["price_usd"] <= max_price)]
        dropped = len(mols) - len(kept)
        if dropped:
            log(f"[run] price filter dropped {dropped}/{len(mols)} seeds (> ${max_price})")
    kept = kept[:count]
    smis = [m["smiles"] for m in kept]
    log(f"[run] sample_space kept {len(smis)} seed(s) (strategy={strategy}, "
        f"seed={used_seed}; replay with --sample-seed {used_seed}):")
    for m in kept:
        log(f"[run]   seed (${m.get('price_usd')}): {m['smiles']}")
    return smis, used_seed


# --------------------------------------------------------------------------- #
# Convergence metrics (data-only; a dashboard consumes convergence.json)
# --------------------------------------------------------------------------- #
def _convergence_round(ranked, rnd: int, top_k: int) -> dict:
    """Serialize the synthon-space state implied by the scored pool `ranked`
    (best-first Products). Reward = -docking_score (higher = better). Rebuilds
    per-(reaction_class, bb_index, bb) reward stats from each product's attributed
    synthons - i.e. the bandit state, recomputed from the pool it was built on."""
    scores = [p.score for p in ranked if p.score is not None]
    best_score = min(scores) if scores else None
    topk = sorted(scores)[:top_k]
    topk_mean = (sum(topk) / len(topk)) if topk else None

    # (rxn, idx, bb) -> list of (reward, weight)
    arms: dict = defaultdict(list)
    for p in ranked:
        if p.reward is None or not (p.reaction_class and p.bbs):
            continue
        weights = p.bb_weights or {}
        for idx, bb in p.bbs:
            w = weights.get(idx, 1.0)
            if w > 0:
                arms[(p.reaction_class, idx, bb)].append((p.reward, float(w)))

    positions: dict = defaultdict(lambda: defaultdict(list))
    n_confident = 0
    for (rxn, idx, bb), rows in arms.items():
        sw = sum(w for _, w in rows)
        mean = sum(r * w for r, w in rows) / sw if sw else 0.0
        cnt = len(rows)
        if cnt >= 2:
            n_confident += 1
            rs = [r for r, _ in rows]
            m = sum(rs) / len(rs)
            std = math.sqrt(sum((r - m) ** 2 for r in rs) / (len(rs) - 1))
        else:
            std = 0.0
        positions[rxn][str(idx)].append(
            {"bb": bb, "mean_reward": round(mean, 4), "count": cnt, "std": round(std, 4)})
    # rank synthons per position by mean_reward desc, keep top ~12
    for rxn in positions:
        for idx in positions[rxn]:
            positions[rxn][idx].sort(key=lambda d: d["mean_reward"], reverse=True)
            positions[rxn][idx] = positions[rxn][idx][:12]

    return {
        "round": rnd, "n_docks": len(scores),
        "best_score": None if best_score is None else round(best_score, 4),
        "topk_mean": None if topk_mean is None else round(topk_mean, 4),
        "top_k": top_k,
        "scores": [round(s, 4) for s in scores],
        "positions": {r: dict(v) for r, v in positions.items()},
        "n_synthons_seen": len(arms),
        "n_synthons_confident": n_confident,
    }


def _write_convergence(path: str, meta: dict, history: list) -> None:
    with open(path, "w") as fh:
        json.dump({**meta, "rounds": history}, fh, indent=2)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="One-command synthon-TS screen with pose export for PyMOL")
    # target / pocket
    ap.add_argument("--pdb", required=True, help="PDB id (auto-downloaded) or local .pdb/.cif path")
    ap.add_argument("--pocket-ligand", required=True, help="resname of the bound pocket ligand (e.g. ATP)")
    ap.add_argument("--pocket-chain", default=None, help="optional chain id to disambiguate the ligand")
    ap.add_argument("--padding", type=float, default=8.0, help="Angstroms padding around the ligand box")
    # seeding
    ap.add_argument("--query", default=None,
                    help="seed SMILES (a known binder). Omit to seed from onepot sample_space.")
    ap.add_argument("--sample-seeds", type=int, default=3, help="# seeds to draw when --query is omitted")
    ap.add_argument("--sample-strategy", choices=["diverse", "random"], default="diverse")
    ap.add_argument("--sample-seed", type=int, default=None, help="uint32 to replay a sample_space draw")
    ap.add_argument("--mw-min", type=float, default=300.0, help="sample_space MW floor (seeds only)")
    ap.add_argument("--qed-min", type=float, default=0.5, help="sample_space QED floor (seeds only)")
    # drug-like window (drives BOTH sample_space seeds AND the local analog filter)
    ap.add_argument("--mw-max", type=float, default=550.0)
    ap.add_argument("--logp-min", type=float, default=1.0)
    ap.add_argument("--logp-max", type=float, default=5.0)
    ap.add_argument("--hbd-max", type=int, default=5)
    ap.add_argument("--hba-max", type=int, default=10)
    ap.add_argument("--tpsa-max", type=float, default=140.0)
    ap.add_argument("--rot-max", type=int, default=8)
    ap.add_argument("--no-druglike-filter", action="store_true",
                    help="disable the local RDKit drug-like filter on retrieved analogs")
    # cost / supplier (native onepot search filters; local price cut on seeds)
    ap.add_argument("--max-price", type=int, default=200,
                    help="max analog price USD (onepot tiers are $125/$295; 200 keeps the $125 tier)")
    ap.add_argument("--max-supplier-risk", default="low",
                    choices=["low", "medium", "high"],
                    help="max supplier risk for analogs (native onepot search filter)")
    # backend / docking
    ap.add_argument("--backend", choices=["direct", "muni", "slurm"], default="direct")
    ap.add_argument("--scoring", default="vina", choices=["vina", "vinardo"])
    ap.add_argument("--executable", default="qvina2", choices=["qvina2", "qvina-w", "vina"])
    ap.add_argument("--exhaustiveness", type=int, default=8)
    # SLURM backend (your own HPC — massive CPU/GPU parallelism, no cloud queue)
    ap.add_argument("--engine", default="qvina2",
                    choices=["qvina2", "smina", "unidock", "gnina"],
                    help="SLURM engine: CPU->GPU funnel (qvina2/smina CPU; unidock GPU throughput; gnina GPU+CNN)")
    ap.add_argument("--partition", default=None,
                    help="SLURM partition for CPU engines (qvina2/smina). TODO: set for your cluster")
    ap.add_argument("--gpu-partition", default=None,
                    help="SLURM partition for GPU engines (unidock/gnina). TODO: set for your cluster")
    ap.add_argument("--slurm-time", default="02:00:00", help="SLURM #SBATCH --time per array task")
    ap.add_argument("--array-chunk", type=int, default=200, help="ligands per SLURM array task")
    ap.add_argument("--slurm-poll", type=int, default=30, help="seconds between squeue polls")
    ap.add_argument("--slurm-workdir", default="slurm_work",
                    help="shared-filesystem scratch dir for SLURM ligands/poses/scores")
    ap.add_argument("--env-activate", default=None,
                    help="shell snippet prepended to each SLURM job (module load / micromamba activate). "
                         "TODO: set for your cluster")
    ap.add_argument("--dock-mode", choices=["single", "batch"], default="single",
                    help="direct backend only: 'single' returns poses (~1.5-2 cr/lig, "
                         "8 concurrent, slow); 'batch' is scores-only but ~15-20x "
                         "cheaper/faster — use it to screen THOUSANDS (poses re-dock top-K)")
    ap.add_argument("--batch-chunk", type=int, default=60, help="ligands per batch_docking job")
    ap.add_argument("--max-inflight", type=int, default=8, help="concurrent docking jobs/waves")
    ap.add_argument("--page-id", default=None, help="muni page id (muni backend only)")
    # budget
    ap.add_argument("--max-docks", type=int, default=120)
    ap.add_argument("--max-rounds", type=int, default=6)
    ap.add_argument("--seed-hits", type=int, default=60)
    ap.add_argument("--round-hits", type=int, default=30)
    ap.add_argument("--anchors", type=int, default=3)
    ap.add_argument("--coarse", action="store_true", help="coarse attribution (cheaper onepot)")
    ap.add_argument("--belief-store", default=None,
                    help="shared synthon-belief JSON: loaded as priors at start, updated "
                         "at end, so campaigns COMPOUND across runs instead of restarting.")
    ap.add_argument("--seed-from-best", type=int, default=0, metavar="N",
                    help="also seed from the top-N global hits (runs/top10.csv) so new "
                         "runs elaborate proven winners, not just random draws.")
    # output
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--top-k", type=int, default=10, help="# top hits to export as poses")
    ap.add_argument("--no-poses", action="store_true",
                    help="skip pose export entirely (scores + convergence + hits.csv only). "
                         "Use with --backend muni, or when Rowan pose storage is unavailable.")
    ap.add_argument("--pymol", action="store_true", help="auto-launch PyMOL on the output")
    return ap


def _load_dotenv():
    """Optional convenience: load a local .env into os.environ if python-dotenv is
    installed (`pip install python-dotenv`). No-op otherwise — then export the keys
    yourself (`set -a; source .env; set +a`). See .env.example."""
    try:
        from dotenv import load_dotenv
        load_dotenv()
    except Exception:
        pass


def _seeds_from_leaderboard(n):
    """Top-n accumulated global-best SMILES (runs/top10.csv) to seed a new run from
    proven winners — the cross-run 'use everything we've learned' seeding."""
    import csv
    path = "runs/top10.csv"
    if not os.path.exists(path):
        return []
    out = []
    try:
        for r in csv.DictReader(open(path)):
            smi = r.get("smiles")
            if smi:
                out.append(smi)
            if len(out) >= n:
                break
    except Exception:
        pass
    return out


def main(argv=None):
    _load_dotenv()
    args = _build_parser().parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    onepot_key = os.environ["ONEPOT_API_KEY"]
    # Rowan is needed for the `direct` docking backend and for Rowan-based pose
    # export (direct + muni). The `slurm` backend docks AND poses locally on your
    # cluster, so it needs NO Rowan key.
    rowan_key = os.environ.get("ROWAN_API_KEY")
    need_rowan = (args.backend == "direct"
                  or (args.backend == "muni" and not args.no_poses))
    if need_rowan and not rowan_key:
        sys.exit("ROWAN_API_KEY is required for the direct backend and Rowan pose export.")

    limits = DrugLikeLimits(mw_max=args.mw_max, logp_min=args.logp_min,
                            logp_max=args.logp_max, hbd_max=args.hbd_max,
                            hba_max=args.hba_max, tpsa_max=args.tpsa_max,
                            rotatable_max=args.rot_max)

    # 1) box from the pocket ligand (raw PDB, ligand still present) --------
    print(f"[run] deriving box from ligand {args.pocket_ligand!r} in {args.pdb} ...")
    box = box_from_ligand(args.pdb, args.pocket_ligand, chain=args.pocket_chain,
                          padding=args.padding)
    print(f"[run] box center={[round(v,2) for v in box[0]]} size={[round(v,2) for v in box[1]]}")

    # 2) seeds: --query OR sample_space -----------------------------------
    onepot_client = Client(api_key=onepot_key, base_url="https://api.onepot.ai")
    if args.query:
        seed_smiles = [args.query]
        print(f"[run] seeding from --query: {args.query}")
    else:
        seed_smiles, _used = _sample_seeds(
            onepot_client, args.sample_seeds, args.sample_strategy, limits,
            args.mw_min, args.qed_min, args.sample_seed, args.max_price)
        if not seed_smiles:
            sys.exit("sample_space returned no seeds; relax the property window.")
    # cross-run knowledge: also seed from the accumulated global-best hits
    if args.seed_from_best:
        best = _seeds_from_leaderboard(args.seed_from_best)
        if best:
            seed_smiles = best + [s for s in seed_smiles if s not in best]
            print(f"[run] +{len(best)} seeds from accumulated best hits (top10.csv)")

    # 3) protein prep + docking oracle ------------------------------------
    if args.backend == "direct":
        rowan_uuid = _resolve_rowan_uuid(args.pdb, rowan_key)
        dock_target = Target(protein=rowan_uuid, pocket=box, executable=args.executable,
                             scoring_function=args.scoring, exhaustiveness=args.exhaustiveness)
        if args.dock_mode == "batch":
            from .oracle_rowan_batch import RowanBatchDockingOracle
            oracle = RowanBatchDockingOracle(dock_target, api_key=rowan_key,
                                             protein_uuid=rowan_uuid, chunk=args.batch_chunk,
                                             max_inflight=args.max_inflight, name="synthon-TS dock")
            print(f"[run] direct BATCH docking (scores-only; poses re-dock top-K); "
                  f"chunk={args.batch_chunk}, {args.max_inflight} jobs concurrent")
        else:
            from .oracle_rowan import RowanDockingOracle
            oracle = RowanDockingOracle(dock_target, api_key=rowan_key, protein_uuid=rowan_uuid,
                                        max_inflight=args.max_inflight, name="synthon-TS dock")
        pose_target = dock_target  # same prepared receptor
    elif args.backend == "slurm":
        from .oracle_slurm import SlurmDockingOracle
        dock_target = Target(protein=args.pdb, pocket=box, executable=args.engine,
                             scoring_function=args.scoring, exhaustiveness=args.exhaustiveness)
        oracle = SlurmDockingOracle(
            dock_target, engine=args.engine, partition=args.partition,
            gpu_partition=args.gpu_partition, time_limit=args.slurm_time,
            array_chunk=args.array_chunk, poll_interval=args.slurm_poll,
            workdir=args.slurm_workdir, env_activate=args.env_activate)
        print(f"[run] SLURM backend: engine={args.engine}, array-chunk={args.array_chunk}, "
              f"partition={args.gpu_partition if args.engine in ('unidock','gnina') else args.partition}")
        # SLURM produces poses locally (gzipped SDF) — no Rowan receptor needed.
        pose_target = None
    else:
        from .oracle_muni import MuniBatchDockingOracle
        dock_target = Target(protein=args.pdb, pocket=box, executable=args.executable,
                             scoring_function=args.scoring, exhaustiveness=args.exhaustiveness)
        oracle = MuniBatchDockingOracle(dock_target, page_id=args.page_id,
                                        name="synthon-TS dock")
        # muni produces no poses -> prepare a separate Rowan receptor for export,
        # UNLESS --no-poses (e.g. the Rowan account can't store poses right now).
        if args.no_poses:
            pose_target = None
        else:
            print("[run] muni backend: preparing a Rowan receptor for pose export ...")
            rowan_uuid = _resolve_rowan_uuid(args.pdb, rowan_key)
            pose_target = Target(protein=rowan_uuid, pocket=box, executable=args.executable,
                                 scoring_function=args.scoring, exhaustiveness=args.exhaustiveness)

    # 4) loop config ------------------------------------------------------
    cfg = LoopConfig(query_smiles=seed_smiles[0], seed_smiles=seed_smiles,
                     seed_max_results=args.seed_hits, round_max_results=args.round_hits,
                     n_anchors_per_round=args.anchors, max_rounds=args.max_rounds,
                     max_docks=args.max_docks, top_k=args.top_k,
                     precise_attribution=not args.coarse,
                     max_price=args.max_price, max_supplier_risk=args.max_supplier_risk)
    mol_filter = None if args.no_druglike_filter else make_filter(limits)

    # 5) live pose export + convergence metrics via round callback --------
    pose_cache: dict = {}
    conv_history: list = []
    conv_meta = {"target": args.pdb, "query_or_seeds": seed_smiles,
                 "backend": args.backend}
    conv_path = os.path.join(args.out_dir, "convergence.json")

    def round_callback(ranked, rnd):
        # convergence metrics (data-only, overwrite each round)
        conv_history.append(_convergence_round(ranked, rnd, args.top_k))
        _write_convergence(conv_path, conv_meta, conv_history)
        # live poses
        hits = [(p.smiles, p.score) for p in diverse_top(ranked, args.top_k)]
        if args.no_poses or pose_target is None:
            _write_hits_csv(os.path.join(args.out_dir, "hits.csv"), hits)   # scores only
        elif hits:
            print(f"[run] round {rnd}: exporting {len(hits)} live poses ...")
            export_top_poses(hits, pose_target, args.out_dir, rowan_key,
                             top_k=args.top_k, cache=pose_cache,
                             pose_refs=getattr(oracle, "pose_cache", None))

    res = run_loop(oracle, dock_target, cfg, onepot_key=onepot_key,
                   round_callback=round_callback, mol_filter=mol_filter,
                   belief_store=args.belief_store)

    # final export --------------------------------------------------------
    if not conv_history:  # seed-only run (no elaboration rounds fired)
        conv_history.append(_convergence_round(res["ranked"], 0, args.top_k))
        _write_convergence(conv_path, conv_meta, conv_history)
    final_hits = [(p.smiles, p.score) for p in diverse_top(res["ranked"], args.top_k)]
    if args.no_poses or pose_target is None:
        _write_hits_csv(os.path.join(args.out_dir, "hits.csv"), final_hits)
    elif final_hits:
        print(f"[run] final: exporting {len(final_hits)} poses ...")
        export_top_poses(final_hits, pose_target, args.out_dir, rowan_key,
                         top_k=args.top_k, cache=pose_cache,
                         pose_refs=getattr(oracle, "pose_cache", None))

    # SLURM backend: lean pose output — ONE gzipped multi-mol SDF from the
    # engine's own poses (no Rowan). Skipped with --no-poses.
    if args.backend == "slurm" and not args.no_poses and final_hits:
        from .oracle_slurm import write_top_hits_sdf_gz
        print(f"[run] final: exporting {len(final_hits)} SLURM poses -> top_hits.sdf.gz ...")
        write_top_hits_sdf_gz(final_hits, oracle,
                              os.path.join(args.out_dir, "top_hits.sdf.gz"))

    _report(res, cfg, args.out_dir)

    if args.pymol:
        _launch_pymol(args.out_dir)


def _write_hits_csv(path, hits):
    """Minimal scores-only hits.csv (no poses) — same columns the dashboard/gen_top10
    expect, with empty pose fields — so muni / --no-poses runs still populate the
    leaderboard."""
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "smiles", "score", "mmgbsa", "posebusters_valid", "pose_file"])
        for i, (smi, score) in enumerate(hits, 1):
            w.writerow([i, smi, f"{score:.3f}", "", "", ""])


def _report(res, cfg, out_dir):
    ranked = res["ranked"]
    print(f"\n=== synthon-TS done: {res['n_docks']} docks, {len(res['pool'])} scored ===")
    print("Top scaffold-diverse hits (lower score = better):")
    for p in diverse_top(ranked, cfg.top_k):
        print(f"  {p.score:+6.2f}  {p.smiles}")
    json.dump(res["top"], open(os.path.join(out_dir, "hits.json"), "w"), indent=2)
    view = os.path.join(out_dir, "view.pml")
    print(f"\nWrote hits + poses to {out_dir}")
    print(f"Open the docked poses in PyMOL:\n    {PYMOL_BIN} {view}")
    print("Live refresh while a run is going: in PyMOL run  @view.pml  (or File > Reload All)")


def _launch_pymol(out_dir):
    view = os.path.join(out_dir, "view.pml")
    if not os.path.exists(PYMOL_BIN):
        print(f"[run] --pymol: PyMOL not found at {PYMOL_BIN}; skipping launch.")
        return
    if not os.path.exists(view):
        print(f"[run] --pymol: no {view} to open; skipping launch.")
        return
    print(f"[run] launching PyMOL on {view} ...")
    subprocess.Popen([PYMOL_BIN, view], stdout=subprocess.DEVNULL,
                     stderr=subprocess.DEVNULL)


if __name__ == "__main__":
    main()
