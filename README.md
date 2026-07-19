# Ebolathon — synthon-aware screen over onepot CORE

A one-command virtual screen over onepot's **CORE** makeable chemical space, with
a **docking oracle** as the fitness function. Give it a target PDB + the ligand
that marks the pocket; it derives the box, preps the receptor, seeds the search
(from a known binder *or* a random draw of makeable CORE), then loops
**retrieve → dock → decompose → steer** — concentrating docking on the building
blocks that score well. It writes **top-hit docked poses for PyMOL** and
**per-round convergence metrics**.

The deliverable lives in [`synthon_ts/`](synthon_ts/) — see
[`synthon_ts/README.md`](synthon_ts/README.md) for the full design, the two
docking backends, and how it relates to (and differs from) true Thompson Sampling.

## Quick start

```bash
# 1. environment
pip install onepot rowan-python rdkit          # + PyMOL to view poses

# 2. keys (env only — never commit these)
export ONEPOT_API_KEY=<your-onepot-key>
export ROWAN_API_KEY=<your-rowan-key>

# 3. run — target PDB + the ligand that defines the pocket + a seed
python -m synthon_ts.run \
    --pdb 1HCK --pocket-ligand ATP \
    --query "CCC(CO)Nc1nc(NCc2ccccc2)c2ncn(C(C)C)c2n1" \
    --backend direct --out-dir runs/trial1 \
    --max-docks 8 --seed-hits 6 --round-hits 4 --anchors 1 --max-rounds 1 --top-k 3
```

**No known binder?** Drop `--query` and add `--sample-seeds 3` — it seeds from
onepot `sample_space` (property-filtered makeable CORE molecules).

## What you give it

| Flag | Meaning |
|---|---|
| `--pdb` | PDB id (auto-downloaded) or a local `.pdb`/`.cif` |
| `--pocket-ligand` | resname of the bound ligand that **marks where the pocket is** (only its location is used) |
| `--query` | seed SMILES — the molecule the **chemical search** grows from. Omit → seed from `sample_space` |
| `--backend` | `direct` (Rowan; returns poses; uses your Rowan credits) or `muni` (cheap batch; scores only) |

`--pocket-ligand` = *where to dock*; `--query` = *what chemistry to start from*. They're different roles (see synthon_ts/README.md).

## Filters (defaults)

- **Drug-like window:** MW ≤ 550, cLogP 1–5, TPSA ≤ 140, HBD ≤ 5, HBA ≤ 10, rot-bonds ≤ 8, QED ≥ 0.5.
  Native to `sample_space` for seeds; applied locally (RDKit) to retrieved analogs before docking. (onepot exposes cLogP, not cLogD — cLogP is used as the proxy.)
- **Cost/supply:** `--max-price 200` ($125 tier only) + `--max-supplier-risk low`, native on every onepot search.

## Outputs (in `--out-dir`, refreshed each round)

- `receptor.pdb`, `rank{N}_{score}.pdb` docked poses, `top_hits.sdf`, `hits.csv`
- `view.pml` — open live with `pymol runs/trial1/view.pml`
- `convergence.json` — per-round score trajectory + per-position synthon leaderboards (the space pruning to the best building blocks)

## Notes for an agent running this

- **Backend `direct` is the simple default:** Rowan-only, one credit pool, poses native. Only use `muni` if you want ultra-cheap bulk scoring of thousands of molecules (then poses cost a small Rowan re-dock of the top hits).
- Direct docks are ~1–2 min each — keep `--max-docks` small for a first trial.
- The receptor is prepared to a Rowan protein UUID automatically from the PDB id; you never handle the UUID.
- Credit-free dry run of the loop logic: `python -m synthon_ts.selftest`.
