---
name: vs-pipeline-router
description: Start here for any onepot CORE virtual-screening task. Routes to the right pipeline skill based on what data access you have (building blocks + reactions vs enumerated dump vs search-only API), and lists what to request from onepot/muni.
---

# Virtual Screening Pipeline Router (onepot CORE)

Use this to pick the pipeline. The whole strategy branches on **what access you have to onepot CORE** (~3.4B enumerated compounds; 7 reaction classes × 320,108 building blocks; arXiv 2601.12603). Confirm access first with the [[onepot-core-access]] skill.

## Decision tree

```
Do you have the building-block lists + reaction SMIRKS/SMARTS?
├── YES → Scenario A → skill: synthon-thompson-screen
│         (Thompson Sampling or synthon docking; you can enumerate focused sublibraries yourself)
│
└── NO → Can you get a bulk / subset SMILES dump of enumerated products?
         ├── YES → Scenario B → skill: shape-screen
         │         (2D/property filter → USRCAT prefilter → GPU conformers → roshambo2 overlay → triage)
         │
         └── NO, only a search API (similarity/substructure/shape/docking via onepot or muni)
                  → Scenario C → skill: shape-screen (search-API section)
                    (the API RETRIEVES candidates that a LOCAL scorer (roshambo2) then ranks —
                     seeded analog hunting / iterative similarity expansion; or run muni's hosted
                     shape/docking tool against CORE. A similarity/substructure API is NOT a TS
                     scoring oracle — it returns neighbors of a probe, not a score of a candidate.)
```

## Key judgement call: do you even need Thompson Sampling?

TS earns its keep on **non-enumerated** spaces of 10¹⁰–10¹¹ (Enamine REAL Space). onepot CORE is
**already enumerated at ~3.4B** — only ~70× smaller than REAL. So:

- **If you get the BBs + reactions (Scenario A):** TS / synthon docking is elegant and cheap. Do it.
- **If you only get products or search (B/C):** don't force TS. A hard prefilter + GPU shape screen
  brute-forces 3.4B fine. Use `shape-screen`.

Either way the endpoint is a **GPU shape/pharmacophore screen with roshambo2** — that's the build goal.

## Standard funnel (all scenarios converge here)

```
3.4B ──2D/property filter──► ~100M ──USRCAT (recall-first)──► ~1–10M
     ──GPU conformers (nvMolKit/ETKDG)──► ──roshambo2 shape+color overlay──► ~10k
     ──artifact strip / consensus──► shortlist ──test HUNDREDS (not 20)──► hits
```

Guiding principles (from the literature dossier at `research/ultra-large-vs-research.md`):
score enriches not ranks; consensus across orthogonal error modes; diversity over top-score;
test enough to know your hit rate; the assay is ground truth.

## Related skills
- [[onepot-core-access]] — figure out which scenario you're in; questions to send onepot/muni
- [[gpu-vs-env-setup]] — micromamba env (rdkit + nvmolkit + roshambo2 + cupy)
- [[shape-screen]] — Scenario B/C pipeline (prefilter → roshambo2)
- [[synthon-thompson-screen]] — Scenario A pipeline (TS / synthon docking)

Full background: `research/ultra-large-vs-research.md` in this repo.
