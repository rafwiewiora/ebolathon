"""Offline mocks so the loop logic can be exercised without spending onepot or
docking credits. Not used in production runs.

`MockOnePot.search` matches the REAL onepot response shape: `results` analogs
carry only {smiles, similarity, price_usd} (NO synthons), and the QUERY carries a
`decompositions` block {reaction_class, bbs}. Decompositions are deterministic
per molecule (stable across calls) and drawn from a small BB vocabulary so
different products share synthons — that gives the bandit real signal to
accumulate, and makes the precise-attribution decompose step meaningful."""
from __future__ import annotations

import random
import zlib

_RXNS = ["rxn_amide01", "rxn_suzuki02"]
_BB_VOCAB = 12  # small pool per position -> synthons recur across products


def _h(s: str) -> int:
    """Stable (process-independent) small hash."""
    return zlib.crc32(s.encode())


class MockOnePot:
    """Fakes onepot's search/decompose in the real shape. `search` accepts a
    list of query SMILES and returns one query block each; with decompose=True
    every query block gets a deterministic `decompositions` entry. Honours
    bb_filters only loosely (they don't change the mock analog set)."""

    def __init__(self, seed=0, analogs_per_query=8):
        self.rng = random.Random(seed)
        self.analogs_per_query = analogs_per_query

    def _decomp(self, smi: str) -> dict:
        h = _h(smi)
        rxn = _RXNS[h % len(_RXNS)]
        bb0 = f"BB0_{h % _BB_VOCAB}c1ccccc1"
        bb1 = f"BB1_{(h // _BB_VOCAB) % _BB_VOCAB}c1ccncc1"
        return {"reaction_class": rxn,
                "bbs": [{"bb_index": 0, "smiles": bb0},
                        {"bb_index": 1, "smiles": bb1}]}

    def _query_block(self, qs: str, max_results: int, decompose: bool) -> dict:
        n = min(max_results, self.analogs_per_query)
        results = []
        for i in range(n):
            # first result is the query itself (similarity 1.0), like real onepot
            smi = qs if i == 0 else f"{qs}#{i}"
            results.append({
                "smiles": smi,
                "inchikey": f"MOCK{_h(smi) % 10**10:010d}-N",
                "supplier_risk": "low",
                "has_non_us_bbs": False,
                "similarity": 1.0 if i == 0 else round(self.rng.uniform(0.3, 0.95), 3),
                "price_usd": self.rng.choice([125, 295, 1000]),
            })
        block = {"query_smiles": qs, "query_inchikey": f"MOCKQ{_h(qs) % 10**8:08d}-N",
                 "results": results}
        if decompose:
            block["decompositions"] = [self._decomp(qs)]
        return block

    def search(self, smiles_list, max_results=20, decompose=False,
               bb_filters=None, **kw):
        queries = [self._query_block(qs, max_results, decompose) for qs in smiles_list]
        return {"queries": queries,
                "credits_used": len(smiles_list),
                "credits_remaining": None}


class MockOracle:
    """Cheap deterministic docking surrogate: score ~ -(len + N-count) with
    reproducible per-SMILES noise. Lower is better, matching the Vina
    convention."""

    def __init__(self, seed=0):
        self.seed = seed

    def score(self, smiles_list, template_smiles=None):
        out = {}
        for s in smiles_list:
            base = -(0.03 * len(s) + 0.4 * s.lower().count("n"))
            noise = random.Random(_h(s) ^ self.seed).gauss(0, 0.3)
            out[s] = round(base + noise, 2)
        return out
