#!/usr/bin/env python
"""Aggregate live screen results into:
  - runs/TOP10.md / runs/top10.csv  (human leaderboard)
  - runs/dashboard.json             (single feed the GitHub Pages dashboard fetches)

Deduped by canonical SMILES, best (lowest) docking score kept. No wall-clock
timestamp in the files on purpose — they change only when the data changes, so
the auto-push loop commits only real updates (the dashboard shows its own
client-side refresh time). Run: python research/gen_top10.py"""
from __future__ import annotations

import csv
import glob
import json
import os

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except Exception:  # pragma: no cover
    Chem = None


def _canon(smi: str) -> str:
    if Chem is None:
        return smi
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else smi


def _svg(smi: str, w: int = 220, h: int = 150) -> str:
    """Compact transparent-background 2D depiction (RDKit) for inline embedding.
    Black atoms/bonds — the dashboard wraps these in a white card so they read in
    both light and dark themes. Returns '' if RDKit is missing or parse fails."""
    if Chem is None or not smi:
        return ""
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return ""
        d = rdMolDraw2D.MolDraw2DSVG(w, h)
        o = d.drawOptions()
        o.clearBackground = False          # transparent — card supplies the white
        o.bondLineWidth = 1
        o.padding = 0.08
        rdMolDraw2D.PrepareAndDrawMolecule(d, m)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        # strip the XML/doctype header so it drops straight into HTML
        i = svg.find("<svg")
        return svg[i:].strip() if i >= 0 else ""
    except Exception:
        return ""


def _load_hits():
    rows = []
    for f in sorted(glob.glob("runs/*/hits.csv")):
        run = os.path.basename(os.path.dirname(f))
        if "baseline" in run:      # random-baseline control, not a strategy hit
            continue
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
    return rows


def _top10(rows):
    best: dict = {}
    for r in rows:
        key = _canon(r["smiles"])
        if key not in best or r["score"] < best[key]["score"]:
            best[key] = r
    return sorted(best.values(), key=lambda x: x["score"])[:10], len(best)


def _top_bb_score(rnd, k=5):
    """Per-round synthon (BB) 'score' for the trajectory: the mean reward of the
    top-k synthon arms this round, negated into score space (lower = better) so it
    plots on the same axis as the molecule scores."""
    rews = []
    for _rxn, pos in rnd.get("positions", {}).items():
        for _idx, arr in pos.items():
            for s in arr:
                rews.append(s["mean_reward"])
    rews.sort(reverse=True)
    top = rews[:k]
    return round(-sum(top) / len(top), 2) if top else None


def _primary_run():
    """Pick the run to feature: prefer a live 'screen'/'batch' run that has data,
    else the run with the most docked molecules."""
    cands = []
    for f in sorted(glob.glob("runs/*/convergence.json")):
        run = os.path.basename(os.path.dirname(f))
        try:
            d = json.load(open(f))
        except Exception:
            continue
        rounds = d.get("rounds") or []
        if not rounds:
            continue
        best = rounds[-1].get("best_score")
        if best is None:
            continue
        docked = rounds[-1].get("n_docks", 0)
        live = any(k in run for k in ("screen", "batch"))
        cands.append((live, -best, docked, run, d, rounds))
    if not cands:
        return None
    # feature the live run with the BEST (most negative) score, then most docked
    cands.sort(key=lambda c: (c[0], c[1], c[2]), reverse=True)
    _, _, _, run, d, rounds = cands[0]
    last = rounds[-1]
    positions = last.get("positions", {})
    # positions per reaction (for the "which slot" indicator) = max bb_index + 1
    rxn_npos = {rxn: max((int(k) for k in posd), default=0) + 1
                for rxn, posd in positions.items()}
    syn = []
    for rxn, pos in positions.items():
        for idx, arr in pos.items():
            for s in arr:
                syn.append({"reward": s["mean_reward"], "count": s["count"],
                            "rxn": rxn, "pos": int(idx), "n_pos": rxn_npos.get(rxn, 2),
                            "bb": s["bb"], "svg": _svg(s["bb"])})
    syn.sort(key=lambda x: x["reward"], reverse=True)
    return {
        "run": run,
        "seeds": len(d.get("query_or_seeds", [])),
        "docked": last.get("n_docks", 0),
        "synthons_seen": last.get("n_synthons_seen", 0),
        "synthons_confident": last.get("n_synthons_confident", 0),
        "best_score": last.get("best_score"),
        "topk_mean": last.get("topk_mean"),
        "trajectory": [{"round": r.get("round"), "best": r.get("best_score"),
                        "topk_mean": r.get("topk_mean"), "n_docks": r.get("n_docks"),
                        "top_bb_score": _top_bb_score(r)}
                       for r in rounds],
        "top_synthons": syn[:8],
    }


def main():
    rows = _load_hits()
    top, n_unique = _top10(rows)
    n_runs = len({r["run"] for r in rows})

    # --- TOP10.md / top10.csv ---
    lines = ["# Top 10 molecules — Ebola glycoprotein (T0R site)", "",
             f"Best docking score wins (lower = better). {n_unique} unique molecules "
             f"across {n_runs} run(s). Auto-updated as the screen runs.", "",
             "| # | score | MMGBSA | PoseBusters | source run | SMILES |",
             "|---|------:|-------:|:-----------:|-----------|--------|"]
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

    # --- campaign totals across ALL runs (the whole campaign, not one run) ---
    total_docked = 0
    n_campaigns = 0
    for f in glob.glob("runs/*/convergence.json"):
        run = os.path.basename(os.path.dirname(f))
        if "baseline" in run:
            continue
        try:
            cd = json.load(open(f))
        except Exception:
            continue
        if not cd.get("rounds"):
            continue
        n = cd["rounds"][-1].get("n_docks", 0)
        if n > 0:
            total_docked += n
            n_campaigns += 1

    # --- cumulative synthon knowledge: the belief store (learned across ALL runs) ---
    cum_syn = []
    n_arms = 0
    try:
        bs = json.load(open("runs/belief_store.json"))
        n_arms = len(bs.get("obs", []))
        for rxn, idx, bb, rows in bs.get("obs", []):
            sw = sum(w for _, w in rows)
            if sw <= 0:
                continue
            m = sum(r * w for r, w in rows) / sw
            cum_syn.append({"reward": round(m, 2), "n": round(sw, 1), "rxn": rxn,
                            "pos": int(idx), "bb": bb})
        cum_syn.sort(key=lambda x: (x["reward"], x["n"]), reverse=True)
        cum_syn = cum_syn[:8]
        for s in cum_syn:
            s["svg"] = _svg(s["bb"])
    except Exception:
        pass

    # --- dashboard.json (feed for the live GitHub Pages dashboard) ---
    dash = {
        "primary": _primary_run(),
        "cumulative_synthons": cum_syn,
        "totals": {"docked": total_docked, "campaigns": n_campaigns,
                   "best": (top[0]["score"] if top else None), "unique": n_unique,
                   "synthons_learned": n_arms},
        "n_unique": n_unique,
        "n_runs": n_runs,
        "top10": [{"rank": i, "score": round(r["score"], 2), "smiles": r["smiles"],
                   "run": r["run"], "pb": r["pb"], "svg": _svg(r["smiles"])}
                  for i, r in enumerate(top, 1)],
    }
    with open("runs/dashboard.json", "w") as fh:
        json.dump(dash, fh, indent=2)

    best = f"{top[0]['score']:.2f}" if top else "n/a"
    print(f"wrote TOP10.md + dashboard.json — {len(top)} molecules, best {best}")


if __name__ == "__main__":
    main()
