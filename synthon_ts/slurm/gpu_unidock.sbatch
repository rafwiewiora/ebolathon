#!/bin/bash
# ===========================================================================
# synthon-TS SLURM docking task -- GPU throughput engine (Uni-Dock)
#
# Uni-Dock (DP Technology) is GPU Vina/Vinardo: it batches the WHOLE array-task
# chunk onto one GPU, thousands of ligands/GPU. Rendered by
# synthon_ts.oracle_slurm.SlurmDockingOracle -- {{TOKEN}}s substituted at submit.
# ===========================================================================
#SBATCH --job-name=synthon-dock-unidock
#SBATCH --partition={{PARTITION}}          # TODO: set for your cluster (GPU partition)
#SBATCH --time={{TIME}}
#SBATCH --array={{ARRAY}}
#SBATCH --output=logs/task_%A_%a.out
#SBATCH --error=logs/task_%A_%a.err
{{GRES}}

set -uo pipefail

# --- environment activation (module load cuda / micromamba activate) -------
{{ENV_ACTIVATE}}

RUN_DIR="{{RUN_DIR}}"
cd "$RUN_DIR"
mkdir -p out logs

RECEPTOR="{{RECEPTOR}}"
CX={{CX}}; CY={{CY}}; CZ={{CZ}}
SX={{SX}}; SY={{SY}}; SZ={{SZ}}
EXH={{EXH}}
SCORING="{{SCORING}}"            # vina | vinardo

# The chunk manifest is a Uni-Dock ligand-index file: one ligand .pdbqt path
# per line. Uni-Dock writes each result to out/<ligand_basename>_out.pdbqt
# (e.g. lig_42_out.pdbqt) carrying REMARK VINA RESULT lines.
CHUNK="chunk_${SLURM_ARRAY_TASK_ID}.txt"

unidock \
  --receptor "$RECEPTOR" \
  --ligand_index "$CHUNK" \
  --center_x $CX --center_y $CY --center_z $CZ \
  --size_x $SX --size_y $SY --size_z $SZ \
  --scoring "$SCORING" \
  --exhaustiveness $EXH \
  --search_mode balance \
  --dir out/ \
  > "logs/unidock_${SLURM_ARRAY_TASK_ID}.log" 2>&1 \
  || echo "unidock batch failed for task ${SLURM_ARRAY_TASK_ID}" >&2
