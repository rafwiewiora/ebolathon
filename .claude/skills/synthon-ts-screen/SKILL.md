---
name: synthon-ts-screen
description: Run the one-command synthon-aware Thompson-Sampling screen over onepot CORE with a docking oracle. Use when asked to virtually screen / find hits / grow chemistry against a protein target using this repo — given a target PDB + a pocket-defining ligand (and optionally a seed molecule). Produces docked poses for PyMOL + convergence metrics. This is the CURRENT pipeline (supersedes the older synthon-thompson-screen / shape-screen skills).
---

# Synthon-TS screen (one command)

The working pipeline is the `synthon_ts` package. Run it with one command; it
derives the docking box from a pocket ligand, preps the receptor, seeds from a
known binder or onepot `sample_space`, then loops **retrieve → dock → decompose
→ steer** and writes top-hit poses + convergence data.

## Prereqs
- Python env with `onepot`, `rowan-python`, `rdkit` (here: `/Users/bb/.local/share/mamba/envs/researchia/bin/python`).
- Env keys (never hardcode): `ONEPOT_API_KEY`, and `ROWAN_API_KEY` for the direct backend (`cat ~/.config/ebolathon_rowan_key`).
- PyMOL to view poses: `/Users/bb/.local/share/mamba/envs/cadd-pymol/bin/pymol`.

## Run
```bash
cd <repo>
python -m synthon_ts.run \
  --pdb 1HCK --pocket-ligand ATP \
  --query "<seed SMILES>" \            # omit + use --sample-seeds N to seed from sample_space
  --backend direct --out-dir runs/trial1 \
  --max-docks 8 --seed-hits 6 --round-hits 4 --anchors 1 --max-rounds 1 --top-k 3
```
- `--pocket-ligand` = **where** the pocket is (its coordinates → box). `--query` = **what** chemistry to grow from. Different roles.
- `--backend direct` = Rowan, returns poses, uses the Rowan credits (default, simplest). `--backend muni` = cheap batch docking, **scores only** (no poses).

## Outputs (`--out-dir`, refreshed each round → open in PyMOL live)
`receptor.pdb`, `rank{N}_{score}.pdb` poses, `top_hits.sdf`, `hits.csv`,
`convergence.json`, and `view.pml` (`pymol runs/trial1/view.pml`).

## Gotchas
- Receptor is auto-prepared from the PDB id to a Rowan protein UUID — a **bare PDB id fails at dock time** (`get_protein 400`), which is why prep is mandatory and handled for you.
- Direct docks are ~1–2 min each; keep `--max-docks` small for a first trial.
- Drug-like + price/supplier filters are on by default (MW≤550, cLogP 1–5, TPSA≤140, HBD≤5, HBA≤10, rotB≤8, QED≥0.5; `--max-price 200`, `--max-supplier-risk low`). onepot exposes cLogP (proxy for cLogD). `--no-druglike-filter` to disable.
- Credit-free dry run of the loop: `python -m synthon_ts.selftest`.

## How it differs from true Thompson Sampling
We have onepot (a retriever), not CORE's building-block catalog + reaction SMIRKS.
So instead of *choosing reagents and building* the exact product, we *steer a
Tanimoto window and decompose whatever onepot returns*. Because we decompose
every docked analog, credit assignment is exact per-synthon; the difference is
the acting half (retrieve-in-a-neighborhood vs construct-anything). Full design
+ backends in `synthon_ts/README.md`.
