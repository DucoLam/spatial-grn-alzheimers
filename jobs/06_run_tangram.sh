#!/bin/bash
# Usage: sbatch --job-name=tg_<prefix> 06_run_tangram.sh "<cell_type>" <prefix>
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --job-name=run_tangram
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/tangram_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/tangram_%j.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
CELL_TYPE="${1}"
PREFIX="${2}"
if [ -z "$CELL_TYPE" ] || [ -z "$PREFIX" ]; then
    echo "ERROR: Usage: sbatch 06_run_tangram.sh '<cell_type>' <prefix>"
    exit 1
fi

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
echo "Tangram run: cell_type='$CELL_TYPE'  prefix='$PREFIX'"
echo "Job ID: $SLURM_JOB_ID  Node: $(hostname)  Start: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/06_run_tangram.py \
    --cell_type "$CELL_TYPE" --prefix "$PREFIX"

echo "Done. End: $(date)"
