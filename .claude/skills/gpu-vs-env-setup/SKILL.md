---
name: gpu-vs-env-setup
description: Set up the GPU virtual-screening environment (RDKit + nvMolKit GPU conformers + roshambo2 GPU shape overlay + CuPy + Thompson Sampling). Use before running any shape-screen or synthon-thompson-screen pipeline.
---

# GPU Virtual Screening Environment Setup

Builds the toolchain for GPU shape/pharmacophore screening of onepot CORE.

## Hard constraints (this user)
- **Package management: ONLY `/Users/bb/.local/bin/micromamba`. NEVER conda or mamba.**
- **RDKit 3D: EMBED conformers only** (`EmbedMolecule`/`EmbedMultipleConfs`). **Do NOT run
  MMFF/UFF force-field optimization** unless strictly necessary — it's a known slowdown.
  (nvMolKit's *GPU* MMFF is a separate fast path and is fine.)

## Requirements
- NVIDIA GPU, compute capability ≥ 7.0 (V100+), CUDA toolkit ≥ 12.5, driver ≥ 560.28.
- `nvcc`, a C++ compiler, and CMake (≥3.26) that can find CUDA — needed to build roshambo2.

## 1. Base env + RDKit + nvMolKit + CuPy

```bash
MM=/Users/bb/.local/bin/micromamba
$MM create -n vs -c conda-forge python=3.11 rdkit cupy pytest -y
$MM run -n vs python -c "import rdkit; print('rdkit', rdkit.__version__)"

# nvMolKit = the "NVIDIA RDKit": GPU ETKDG conformers + GPU MMFF + GPU fingerprints/Tanimoto
$MM install -n vs -c conda-forge nvmolkit -y
$MM run -n vs python -c "import nvmolkit; print('nvmolkit ok')"
```
nvMolKit needs a CUDA PyTorch present; if the conda-forge package doesn't pull it, install a
CUDA torch wheel into the env first. Verify GPU is visible: `$MM run -n vs nvidia-smi`.

## 2. roshambo2 (GPU Gaussian shape+color overlay) — build from source

roshambo2 is the core overlay engine (MIT, github.com/molecularinformatics/roshambo2). Unlike
roshambo v1 it does NOT need a custom-compiled RDKit.

```bash
cd /Users/bb/repos/hackathons/ebolathon
git clone https://github.com/molecularinformatics/roshambo2.git
cd roshambo2
# Their environment.yaml pins deps; create with micromamba (NOT conda):
/Users/bb/.local/bin/micromamba env create -n roshambo2 -f environment.yaml
/Users/bb/.local/bin/micromamba run -n roshambo2 pip install .
# sanity check
cd test && /Users/bb/.local/bin/micromamba run -n roshambo2 pytest -q
```
NOTE: the exact Python API (class/args/optim_mode) is documented in the repo's
`USER_GUIDE.md` and `examples/README.md` — **read those before writing screen code**; the
`shape-screen` skill's snippet is based on the published API but confirm signatures there.

Decision: keep roshambo2 in its own env (its environment.yaml is opinionated) and shell out to
it, OR try installing it into the `vs` env. Separate env is safer to start.

## 3. Thompson Sampling (only for Scenario A)

```bash
cd /Users/bb/repos/hackathons/ebolathon
git clone https://github.com/PatWalters/TS.git
# Its built-in ROCSEvaluator needs OpenEye (commercial). For a free GPU path you will write a
# custom Evaluator that calls roshambo2 — see the synthon-thompson-screen skill.
```

## 4. Smoke test (end to end, tiny)
- Generate conformers for 2–3 SMILES with RDKit ETKDG (embed only).
- Run roshambo2 on a 2-molecule query vs 10-molecule library, confirm you get shape+color scores
  and aligned output SDF.
- Time it; record conformers/sec and overlays/sec to size the real run.

## Related
[[vs-pipeline-router]] · [[shape-screen]] · [[synthon-thompson-screen]]
