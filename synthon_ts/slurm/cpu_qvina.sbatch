#!/bin/bash
# ===========================================================================
# synthon-TS SLURM docking task -- CPU / Vina-family engine (qvina2 or smina)
#
# Rendered by synthon_ts.oracle_slurm.SlurmDockingOracle: every {{TOKEN}} below
# is substituted at submit time. #SBATCH lines flagged TODO are cluster-specific
# -- the oracle fills PARTITION/TIME/ARRAY/GRES from your run.py flags, but you
# should still sanity-check them for your site.
# ===========================================================================
#SBATCH --job-name=synthon-dock-cpu
#SBATCH --partition={{PARTITION}}          # TODO: set for your cluster (CPU partition)
#SBATCH --time={{TIME}}
#SBATCH --array={{ARRAY}}
#SBATCH --cpus-per-task=1                   # TODO: bump for more Vina CPU threads (see --cpu below)
#SBATCH --output=logs/task_%A_%a.out
#SBATCH --error=logs/task_%A_%a.err
{{GRES}}

set -uo pipefail

# --- environment activation (module load / micromamba activate) ------------
{{ENV_ACTIVATE}}

RUN_DIR="{{RUN_DIR}}"
cd "$RUN_DIR"
mkdir -p out logs

RECEPTOR="{{RECEPTOR}}"
ENGINE_BIN="{{ENGINE_BIN}}"       # qvina2 or smina
EXTRA_FLAGS="{{EXTRA_FLAGS}}"     # e.g. "--scoring vinardo" for smina/vinardo, "" for qvina2
CX={{CX}}; CY={{CY}}; CZ={{CZ}}
SX={{SX}}; SY={{SY}}; SZ={{SZ}}
EXH={{EXH}}

# Each chunk manifest line is:  <index> <ligand_pdbqt_path>
CHUNK="chunk_${SLURM_ARRAY_TASK_ID}.txt"

while read -r IDX LIG; do
  [ -z "${IDX:-}" ] && continue
  "$ENGINE_BIN" \
    --receptor "$RECEPTOR" \
    --ligand "$LIG" \
    --center_x $CX --center_y $CY --center_z $CZ \
    --size_x $SX --size_y $SY --size_z $SZ \
    --exhaustiveness $EXH $EXTRA_FLAGS \
    --out "out/out_${IDX}.pdbqt" \
    > "logs/dock_${IDX}.log" 2>&1 \
    || echo "dock failed for ligand ${IDX}" >&2
done < "$CHUNK"
