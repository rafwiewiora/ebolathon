---
name: synthon-thompson-screen
description: Combinatorial-aware screen of onepot CORE when you have the building-block lists + reaction SMIRKS (Scenario A). Thompson Sampling over reagents×reactions with a GPU shape (roshambo2) or docking oracle, plus a synthon-docking variant. Use when you can enumerate CORE yourself.
---

# Synthon / Thompson-Sampling Screen (Scenario A)

Use when onepot gives you the **building-block lists + the 7 reaction SMIRKS/SMARTS**. Then CORE is a
searchable reactions×reagents space and you can run Thompson Sampling (Klarich/Walters) or synthon
docking without brute-forcing all 3.4B. Env: run [[gpu-vs-env-setup]] first.

## Why TS fits CORE
CORE = 7 reactions × 320,108 building blocks. TS treats **each reagent as a bandit "arm"** with a
belief over the scores of products containing it; it samples arms, assembles + scores a full product,
updates beliefs, and concentrates on promising reagent combinations — evaluating <1% of the space.
Because it scores **full products** (not isolated fragments) it avoids the fragment-context problem
that plagues fragment-first docking.

## Inputs you need from onepot
1. One **reagent SMILES file per reaction component** (grouped by the reaction's SMARTS, as CORE is
   built). e.g. `rxn3_componentA.smi`, `rxn3_componentB.smi`.
2. The **reaction SMARTS/SMIRKS** for each of the 7 reactions.
(These are described in arXiv 2601.12603; ask onepot for the machine-readable files — see
[[onepot-core-access]].)

## TS setup (github.com/PatWalters/TS)
Config is JSON. Run: `python ts_main.py config.json`. Key fields:
```json
{
  "reaction_smarts": "<SMARTS for one CORE reaction>",
  "reagent_file_list": ["rxnN_componentA.smi", "rxnN_componentB.smi"],
  "evaluator_class_name": "Roshambo2Evaluator",
  "evaluator_arg": {"query_sdf": "query.sdf", "color": true},
  "num_warmup_trials": 3,
  "num_ts_iterations": 1000,
  "ts_mode": "maximize",
  "results_filename": "ts_hits.csv"
}
```
Run one TS job **per reaction** (each reaction is its own combinatorial space), then merge/rank hits.

## The oracle — IMPORTANT
TS's built-in evaluators: `FPEvaluator` (2D fingerprint), `MWEvaluator`, and `ROCSEvaluator`
(**needs OpenEye — commercial**). For a free GPU path, **write a custom evaluator** subclassing TS's
`Evaluator` that scores a product SMILES with roshambo2:

```python
from evaluators import Evaluator          # TS base class
from rdkit import Chem
from rdkit.Chem import AllChem
# from roshambo2 import Roshambo2

class Roshambo2Evaluator(Evaluator):
    def __init__(self, input_dict):
        self.query_sdf = input_dict["query_sdf"]
        self.color = input_dict.get("color", True)
        self.num_evals = 0
    def evaluate(self, mol):               # mol = an assembled product (RDKit Mol)
        self.num_evals += 1
        m = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv3()) != 0:   # embed only, no FF optimize
            return -1.0
        # write m to a temp SDF, run roshambo2 vs self.query_sdf, return best ComboTanimoto
        # calc = Roshambo2(self.query_sdf, tmp_sdf, color=self.color)
        # return float(best_combo_tanimoto(calc.compute(optim_mode="combination")))
        ...
```
Confirm TS's exact base-class name/method and roshambo2's API against their repos before finalizing.
Note: per-product conformer+overlay is the cost; batch products per TS iteration if the API allows.

## Alternatives within Scenario A
- **Self-enumerate focused sublibraries**: once TS surfaces the best reagents per position, enumerate
  just those combinations with the SMIRKS and hand them to `shape-screen` for a dense final pass.
- **Synthon docking (V-SYNTHES style)**: dock capped synthons, keep those whose growth vector points
  into open pocket, elaborate + re-dock. Use if you have a target structure and want structure-based
  (not shape) scoring. Docking via GNINA 1.3 (GPU) or muni.bio's docking tools.
- **Docking oracle for TS**: swap the Roshambo2Evaluator for a docking evaluator (GNINA/Vina) if
  structure-based is preferred over shape.

## Triage
Same as `shape-screen` stage 5: diversity over top-score, artifact/consensus filter, test hundreds.
TS can mode-collapse onto one chemotype — enforce scaffold diversity (or use Enhanced TS
roulette-wheel selection) so you don't get 1000 near-identical hits.

## Related
[[vs-pipeline-router]] · [[gpu-vs-env-setup]] · [[shape-screen]] · [[onepot-core-access]]
Background: `research/ultra-large-vs-research.md`
