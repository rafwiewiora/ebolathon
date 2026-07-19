---
name: onepot-core-access
description: Determine what access you have to onepot CORE / muni.bio and which screening scenario you're in. Covers what onepot CORE and muni.bio are, the access matrix, what each level unlocks, and the exact questions to send onepot/muni.
---

# onepot CORE + muni.bio Access

Run this first — the whole pipeline choice ([[vs-pipeline-router]]) depends on what access you get.

## What the pieces are
- **onepot CORE** (onepot.ai; arXiv 2601.12603): an **enumerated** reaction-based chemical space,
  ~3.4B compounds (product page shows ~2.7B — likely current released subset), built from **7
  med-chem reaction classes × 320,108 building blocks** (SMARTS grouping + SMIRKS enumeration + ML
  feasibility scoring). onepot is a synthesis company (AI chemist "Phil", robotic synthesis, ~5-day
  turnaround, ~$125 entry tier per compound). Access = request a link, agree to CORE terms; their
  terms mention delivering the dataset and **subsets via link, API, file transfer, or email**.
- **muni.bio**: an agentic compute platform ("where molecules meet agents") exposing **~55 tool-models**
  (structure prediction, protein design, small-molecule design, **docking**, etc.) via a code-first
  CLI, GPU-managed. It's where you *run* tools; CORE is the *library*.
- **Relationship**: NOT publicly documented. Unknown whether muni hosts a CORE-search tool or they're
  just partners. **Confirm this directly.**

## Access matrix → scenario
| What you get | Scenario | Pipeline |
|---|---|---|
| Building-block SMILES files (per reaction component) + reaction SMIRKS | **A** | [[synthon-thompson-screen]] (TS / synthon docking / self-enumeration) |
| Bulk or subset SMILES dump of enumerated products | **B** | [[shape-screen]] (prefilter → roshambo2) |
| Search-only API (2D similarity / substructure / shape / docking) | **C** | [[shape-screen]] §"Scenario C" (API as oracle / SpaceLight expansion / muni docking) |
| 3D conformers provided? | any | if not, generate locally (nvMolKit/ETKDG) |

**Highest-leverage single ask:** the building-block + reaction-SMIRKS files. That one artifact
unlocks Scenario A (TS, synthon docking) AND lets you self-enumerate any focused sublibrary for B.
The decomposition is already published (7 reactions, 320,108 BBs), so it's not secret.

## Questions to send onepot / muni
1. Can we get the **building-block lists (SMILES, per reaction component) and the 7 reaction
   SMIRKS/SMARTS** in machine-readable form? (unlocks synthon/TS methods + self-enumeration)
2. Failing that, can we get a **bulk or filtered subset dump** of enumerated product SMILES (and IDs)?
   What subset sizes / filters are supported?
3. What exactly does the **onepot API** do — retrieve/filter enumerated products, or run **searches**?
   If search: which modalities (2D similarity, substructure, **shape**, **pharmacophore**, docking)?
   Auth, rate limits, batch size?
4. Does **muni.bio** expose a tool that runs **against CORE** (a shape or docking model pointed at the
   3.4B), or do we bring our own compute and only pull candidates from onepot?
5. Are **3D conformers** provided, or do we generate them ourselves?
6. Any restrictions in the CORE terms on **bulk descriptor computation / local storage** of the space?

## How to probe programmatically (once you have credentials)
- Hit the API root / docs endpoint; enumerate available routes.
- Test a tiny similarity/substructure query with a known active; inspect the response schema (does it
  return reagent/reaction provenance per hit? that partially reconstructs Scenario A).
- On muni: list the CLI's available tool-models; check for shape/pharmacophore/docking tools and
  whether any accept a CORE handle as the library.

## Related
[[vs-pipeline-router]] · [[shape-screen]] · [[synthon-thompson-screen]] · [[gpu-vs-env-setup]]
Background: `research/ultra-large-vs-research.md`
