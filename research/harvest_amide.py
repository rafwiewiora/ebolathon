#!/usr/bin/env python
"""Reaction harvest for rxn_a46303e0 (amide coupling).

Take the top-scoring acid synthons (pos0) x top amine synthons (pos1) the bandit
surfaced, enumerate the full amide-coupling grid, drug-like-filter, flag which
products are NOVEL (never docked) vs already-scored, and check makeability via
onepot exact_lookup. Writes runs/harvest/candidates.{csv,smi} for docking.

This is the "exploitation harvest": did the best pieces combine into something
better than anything the sampling pass actually retrieved?
"""
import csv
import glob
import json
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

import sys
sys.path.insert(0, ".")
from synthon_ts.filters import DrugLikeLimits, passes_druglike

RXN = "rxn_a46303e0"
N_ACID, N_AMINE = 14, 14          # top-N per position
AMIDE = AllChem.ReactionFromSmarts("[C:1](=[O:2])[OX2H1].[N;!$(NC=O):3]>>[C:1](=[O:2])[N:3]")
ACID_PAT = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
AMINE_PAT = Chem.MolFromSmarts("[NX3;H1,H2;!$(NC=O)]")


def top_synthons():
    best = {}
    for f in glob.glob("runs/*/convergence.json"):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        if not d.get("rounds"):
            continue
        for idx, arr in d["rounds"][-1].get("positions", {}).get(RXN, {}).items():
            for s in arr:
                k = (int(idx), s["bb"])
                best[k] = max(best.get(k, -9e9), s["mean_reward"])
    acids = [(r, bb) for (p, bb), r in best.items() if p == 0
             and Chem.MolFromSmiles(bb) and Chem.MolFromSmiles(bb).HasSubstructMatch(ACID_PAT)]
    amines = [(r, bb) for (p, bb), r in best.items() if p == 1
              and Chem.MolFromSmiles(bb) and Chem.MolFromSmiles(bb).HasSubstructMatch(AMINE_PAT)]
    acids.sort(reverse=True); amines.sort(reverse=True)
    return acids[:N_ACID], amines[:N_AMINE]


def enumerate_grid(acids, amines):
    limits = DrugLikeLimits()
    rows = {}   # canonical smiles -> record
    for ar, a in acids:
        am = Chem.MolFromSmiles(a)
        for nr, n in amines:
            nm = Chem.MolFromSmiles(n)
            for prod in AMIDE.RunReactants((am, nm)):
                try:
                    m = prod[0]; Chem.SanitizeMol(m)
                    smi = Chem.MolToSmiles(m)
                except Exception:
                    continue
                if not passes_druglike(smi, limits):
                    continue
                # keep the best (acid_reward + amine_reward) provenance if dup
                score = ar + nr
                if smi not in rows or score > rows[smi]["combo_reward"]:
                    rows[smi] = {"smiles": smi, "acid": a, "amine": n,
                                 "acid_reward": round(ar, 2), "amine_reward": round(nr, 2),
                                 "combo_reward": round(score, 2)}
    return list(rows.values())


def already_scored():
    seen = set()
    for f in glob.glob("runs/*/hits.csv"):
        try:
            for r in csv.DictReader(open(f)):
                m = Chem.MolFromSmiles(r.get("smiles", ""))
                if m:
                    seen.add(Chem.MolToSmiles(m))
        except Exception:
            pass
    return seen


def check_makeable(smiles_list):
    """onepot exact_lookup: which enumerated products are actually purchasable/CORE."""
    try:
        from onepot import Client
        c = Client(api_key=os.environ["ONEPOT_API_KEY"])
        makeable = {}
        for i in range(0, len(smiles_list), 50):
            chunk = smiles_list[i:i + 50]
            resp = c.search(chunk, exact_lookup=True)
            for q in resp.get("queries", []):
                res = q.get("results", [])
                makeable[q["query_smiles"]] = (res[0].get("price_usd") if res else None)
        return makeable
    except Exception as e:
        print(f"[harvest] makeability check skipped: {str(e)[:100]}")
        return {}


def main():
    os.makedirs("runs/harvest", exist_ok=True)
    acids, amines = top_synthons()
    print(f"[harvest] top {len(acids)} acids x {len(amines)} amines")
    grid = enumerate_grid(acids, amines)
    seen = already_scored()
    for r in grid:
        r["novel"] = r["smiles"] not in seen
    grid.sort(key=lambda r: r["combo_reward"], reverse=True)
    print(f"[harvest] {len(grid)} unique drug-like products "
          f"({sum(r['novel'] for r in grid)} novel / never docked)")

    make = check_makeable([r["smiles"] for r in grid])
    for r in grid:
        r["price_usd"] = make.get(r["smiles"])
        r["makeable"] = r["price_usd"] is not None
    n_make = sum(r["makeable"] for r in grid)
    print(f"[harvest] {n_make}/{len(grid)} makeable via onepot exact_lookup")

    with open("runs/harvest/candidates.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["smiles", "combo_reward", "acid_reward",
                                           "amine_reward", "novel", "makeable",
                                           "price_usd", "acid", "amine"])
        w.writeheader()
        for r in grid:
            w.writerow(r)
    # SMILES for docking — prefer makeable, else all
    dock = [r["smiles"] for r in grid if r["makeable"]] or [r["smiles"] for r in grid]
    with open("runs/harvest/candidates.smi", "w") as fh:
        fh.write("\n".join(dock) + "\n")
    print(f"[harvest] wrote runs/harvest/candidates.csv + candidates.smi "
          f"({len(dock)} to dock)")
    print("[harvest] top combos by summed synthon reward:")
    for r in grid[:8]:
        tag = "NOVEL" if r["novel"] else "seen "
        mk = f"${r['price_usd']}" if r["makeable"] else "not-makeable"
        print(f"  [{tag}] combo={r['combo_reward']:+.1f} {mk:12} {r['smiles']}")


if __name__ == "__main__":
    main()
