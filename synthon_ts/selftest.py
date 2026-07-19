"""Credit-free dry run of the whole loop with mocked onepot + docking.

    python -m synthon_ts.selftest
"""
from __future__ import annotations

from .core import LoopConfig, Target, diverse_top, run_loop
from .mocks import MockOnePot, MockOracle

try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")   # mock SMILES are intentionally fake
except Exception:
    pass


def main():
    target = Target(protein="1HCK", pocket=[[98.49, 94.42, 98.23], [22.6, 24.5, 19.8]])
    cfg = LoopConfig(query_smiles="c1ccc(NC(=O)c2ccccc2)cc1", seed_max_results=8,
                     round_max_results=8, n_anchors_per_round=2, max_rounds=4,
                     patience=2, max_docks=120, top_k=10, random_seed=1)
    res = run_loop(MockOracle(), target, cfg, client=MockOnePot())
    print(f"\n=== selftest OK: {res['n_docks']} docks, {len(res['pool'])} scored ===")
    for p in diverse_top(res["ranked"], 5):
        print(f"  {p.score:+.2f}  rxn={p.reaction_class}  {p.smiles[:44]}")
    assert res["n_docks"] > 0 and res["pool"], "loop produced no scored products"
    print("assertions passed.")


if __name__ == "__main__":
    main()
