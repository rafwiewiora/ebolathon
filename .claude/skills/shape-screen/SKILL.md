---
name: shape-screen
description: GPU shape/pharmacophore screen of onepot CORE when you have enumerated products or only a search API (Scenarios B/C). Funnel = 2D/property filter → USRCAT prefilter → GPU conformers → roshambo2 overlay → triage. Use for the brute-force-with-prefilter path.
---

# GPU Shape/Pharmacophore Screen (Scenarios B & C)

For when you have enumerated CORE products (B) or only a search API (C). At 3.4B enumerated,
prefilter hard, then GPU shape-overlay the survivors. Env: run [[gpu-vs-env-setup]] first.

## The funnel
```
input SMILES ──(1) 2D/property filter──► ~100M
             ──(2) USRCAT prefilter (recall-first)──► ~1–10M
             ──(3) GPU conformers──► ──(4) roshambo2 shape+color overlay──► ~10k
             ──(5) triage/consensus──► shortlist
```

## 0. Query prep
Collect known actives (or a reference ligand). Generate multi-conformer query SDF (ETKDG, embed
only). Decide scoring: **shape only** vs **shape+color** (color = pharmacophore features). Default
to **shape+color / ComboTanimoto** — best enrichment.

## 1. Cheap 2D / property filter
Drop PAINS/reactive groups; cap MW/logP (large libraries drift high — counter it deliberately).
RDKit `Descriptors` + `FilterCatalog`. This is the biggest, cheapest cut.

## 2. USRCAT prefilter (alignment-free, recall-oriented)
USRCAT = 60-D moment descriptor (shape + coarse pharmacophore), no alignment. ~1000× cheaper than
overlay. Keep it **permissive** so true actives survive to stage 4.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

params = AllChem.ETKDGv3()

def usrcat_desc(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(m, params) != 0:   # EMBED ONLY — no MMFF optimize
        return None
    return rdMolDescriptors.GetUSRCAT(m)          # 60 floats (per conformer)

# query descriptor(s) from known active(s)
qd = usrcat_desc("<active_smiles>")

def keep(smiles, cutoff=0.5):                     # permissive cutoff — tune for recall
    d = usrcat_desc(smiles)
    if d is None:
        return False
    return rdMolDescriptors.GetUSRScore(qd, d) >= cutoff  # optional per-feature weights arg
```
Scale: 60 float32 = 240 B/conf, so ~0.8 TB for 3.4B×1conf — that's why you run stage 1 FIRST and
only USRCAT the reduced pool. For GPU: compute the 60-D vectors then do batched Manhattan distance
in CuPy/torch (no custom kernel). Multi-conformer: take best score per molecule.

## 3. GPU conformers for survivors
Generate a few conformers each for the ~1–10M that pass. Options:
- **nvMolKit GPU-ETKDG** (fastest; RDKit-native GPU). Optionally nvMolKit GPU-MMFF (its GPU MMFF is
  fine to use — the "no FF optimize" rule targets the slow CPU RDKit path).
- **RDKit ETKDGv3** `EmbedMultipleConfs`, embed only, parallelized across cores.
Write to SDF/HDF5 for roshambo2.

## 4. roshambo2 GPU shape+color overlay (precise scoring)
This is the precise stage. **Read `roshambo2/USER_GUIDE.md` + `examples/` for the exact current
API before coding** — published API shape:
```python
from roshambo2 import Roshambo2
calc = Roshambo2("query.sdf", "library_confs.sdf", color=True)  # color=True → pharmacophore term
scores = calc.compute(optim_mode="combination")                  # shape + color combined
calc.write_best_fit_structures(hits_sdf_prefix="hits")           # aligned poses out
```
Use its multi-GPU **server mode** + **HDF5** dataset mode for large libraries. Rank by ComboTanimoto.
Benchmark throughput on a small batch first to size the full run.

## 5. Triage (do NOT ship the literal top-N)
- Cluster hits by scaffold; pick **diverse, chemically sane, well-posed** reps across a score band.
- Optional consensus / artifact strip: GNINA 1.3 CNN rescore, strain check, or a docking cross-filter
  (via muni.bio's docking tools) — methods that fail *differently* from shape.
- Plan to test **hundreds** to estimate hit rate honestly (per Lyu 2025).

## Scenario C specifics (search-only API, no bulk data)
Key framing: a **similarity/substructure API is a RETRIEVER, not a scoring oracle** — it returns
molecules near a probe, not a score for an arbitrary candidate. So the API narrows the space and
**roshambo2 does the real scoring locally**. You can't materialize the library, so:
- If onepot/muni expose a **shape or pharmacophore search** endpoint → use it directly as the screen;
  roshambo2 becomes a local rescorer of returned hits.
- If only **2D similarity/substructure** → **seeded analog hunting / iterative similarity expansion**:
  seed with actives, pull nearest neighbors from CORE, shape-score them locally with roshambo2, feed
  the best-scoring back as new query seeds, repeat (greedy hill-climb through similarity space). This
  is your combinatorial-aware substitute when you lack the building blocks — NOT Thompson Sampling
  (TS needs building blocks + reactions + a scoring oracle).
- If muni hosts a docking tool against CORE → retrieve candidates via the API, dock survivors on muni.

**Confirm one thing first (via [[onepot-core-access]]):** does the search cover the WHOLE 3.4B and
just cap the returned top-K (normal, usable — run more/multi-seed queries), or does it only search a
pre-selected SUBSET (you can't reach most of CORE → fall back to a bulk dump / Scenario B)?
"Limited results" is fine; "limited coverage" is not.

## Related
[[vs-pipeline-router]] · [[gpu-vs-env-setup]] · [[onepot-core-access]] · [[synthon-thompson-screen]]
Background: `research/ultra-large-vs-research.md`
