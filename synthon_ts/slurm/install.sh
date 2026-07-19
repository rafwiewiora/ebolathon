#!/usr/bin/env bash
# ===========================================================================
# Install the SLURM docking-engine environment for the synthon-TS screen.
#
# Uses micromamba ONLY (never conda/mamba). Run this once on your cluster (on a
# node/login shell that can see conda-forge + bioconda). The engines + prep
# tools then live in the env named below; point run.py --env-activate at it.
#
# NOTE: gnina and unidock are GPU (CUDA) engines and are Linux-only on
# bioconda. On macOS / non-CUDA nodes install just the CPU engines (qvina,
# smina) + prep tools; comment out the GPU lines. If bioconda lacks a
# prebuilt gnina/unidock for your CUDA/toolchain, fall back to the upstream
# release binaries (github.com/gnina/gnina, github.com/dptech-corp/Uni-Dock)
# and drop them on PATH in your --env-activate snippet instead.
# ===========================================================================
set -euo pipefail

# --- parameters (edit for your cluster) ------------------------------------
MICROMAMBA="${MICROMAMBA:-/Users/bb/.local/bin/micromamba}"   # TODO: micromamba path on your cluster
ENV_NAME="${ENV_NAME:-synthon-dock}"                          # TODO: env name (must match --env-activate)

CHANNELS="-c conda-forge -c bioconda"

echo ">> Creating micromamba env '${ENV_NAME}' with docking engines + prep tools ..."
"${MICROMAMBA}" create -y -n "${ENV_NAME}" ${CHANNELS} \
  python=3.10 \
  rdkit \
  numpy \
  \
  `#### ligand/receptor prep ####` \
  meeko \
  openbabel \
  \
  `#### CPU / cheap bulk -- Vina-family (matches the qvina2 calibration) ####` \
  qvina \
  smina \
  \
  `#### GPU / throughput -- Uni-Dock (GPU Vina/Vinardo, thousands of ligands/GPU) ####` \
  unidock \
  \
  `#### GPU / accuracy -- GNINA 1.3 (Vina + CNN rescoring) ####` \
  gnina

echo ">> Done."
echo ">> Activate with:  ${MICROMAMBA} activate ${ENV_NAME}"
echo ">> In run.py pass e.g.:  --env-activate 'eval \"\$(${MICROMAMBA} shell hook -s bash)\" && micromamba activate ${ENV_NAME}'"
