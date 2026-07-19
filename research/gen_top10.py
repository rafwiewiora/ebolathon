#!/usr/bin/env python
"""Aggregate the best docked molecules across ALL runs into a live top-10
leaderboard (runs/TOP10.md + runs/top10.csv). Deduped by canonical SMILES,
best (lowest) docking score kept. No timestamp on purpose — the file only
changes when the actual top-10 changes, so the auto-push loop commits only real
updates. Run: python research/gen_top10.py"""
from __future__ import annotations

import csv
import glob
import os

try:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except Exception:  # pragma: no cover
    Chem = None


def _canon(smi: str) -> str:
    if Chem is None:
        return smi
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else smi


def main():
    rows = []
    for f in sorted(glob.glob("runs/*/hits.csv")):
        run = os.path.basename(os.path.dirname(f))
        try:
            with open(f) as fh:
                for r in csv.DictReader(fh):
                    try:
                        score = float(r["score"])
                    except (KeyError, ValueError, TypeError):
                        continue
                    rows.append({"smiles": r.get("smiles", ""), "score": score,
                                 "mmgbsa": r.get("mmgbsa", ""),
                                 "pb": r.get("posebusters_valid", ""), "run": run})
        except Exception:
            continue

    best: dict = {}
    for r in rows:
        key = _canon(r["smiles"])
        if key not in best or r["score"] < best[key]["score"]:
            best[key] = r
    top = sorted(best.values(), key=lambda x: x["score"])[:10]
    n_runs = len({r["run"] for r in rows})

    lines = [
        "# Top 10 molecules — Ebola glycoprotein (T0R site)",
        "",
        f"Best docking score wins (lower = better). {len(best)} unique molecules "
        f"across {n_runs} run(s). Auto-updated as the screen runs.",
        "",
        "| # | score | MMGBSA | PoseBusters | source run | SMILES |",
        "|---|------:|-------:|:-----------:|-----------|--------|",
    ]
    for i, r in enumerate(top, 1):
        lines.append(f"| {i} | {r['score']:.2f} | {r['mmgbsa']} | {r['pb']} | "
                     f"{r['run']} | `{r['smiles']}` |")
    with open("runs/TOP10.md", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    with open("runs/top10.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "score", "mmgbsa", "posebusters_valid", "run", "smiles"])
        for i, r in enumerate(top, 1):
            w.writerow([i, f"{r['score']:.3f}", r["mmgbsa"], r["pb"], r["run"], r["smiles"]])

    print(f"wrote runs/TOP10.md — {len(top)} molecules "
          f"(best {top[0]['score']:.2f})" if top else "wrote runs/TOP10.md — no hits yet")


if __name__ == "__main__":
    main()
