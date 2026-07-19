#!/bin/bash
# ===========================================================================
# synthon-TS SLURM docking task -- GPU accuracy engine (GNINA 1.3)
#
# GNINA = AutoDock Vina + CNN rescoring on GPU. Docks per ligand (lower GPU
# throughput than Uni-Dock, higher accuracy). Rendered by
# synthon_ts.oracle_slurm.SlurmDockingOracle -- {{TOKEN}}s substituted at submit.
# ===========================================================================
#SBATCH --job-name=synthon-dock-gnina
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
EXTRA_FLAGS="{{EXTRA_FLAGS}}"    # e.g. "--cnn_scoring rescore"

# Each chunk manifest line is:  <index> <ligand_sdf_path>
# Output out/out_<index>.sdf carries <minimizedAffinity> (Vina kcal/mol, lower
# = better -- matches the qvina2 calibration) plus <CNNaffinity>/<CNNscore>.
CHUNK="chunk_${SLURM_ARRAY_TASK_ID}.txt"

while read -r IDX LIG; do
  [ -z "${IDX:-}" ] && continue
  gnina \
    --receptor "$RECEPTOR" \
    --ligand "$LIG" \
    --center_x $CX --center_y $CY --center_z $CZ \
    --size_x $SX --size_y $SY --size_z $SZ \
    --exhaustiveness $EXH $EXTRA_FLAGS \
    --out "out/out_${IDX}.sdf" \
    > "logs/gnina_${IDX}.log" 2>&1 \
    || echo "gnina failed for ligand ${IDX}" >&2
done < "$CHUNK"
