#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=01:30:00
#SBATCH --job-name=sq2_purityk
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_purityk_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_purityk_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
set -e
DFLAM=/tudelft.net/staff-umbrella/ScReNI/dflam
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "start $(date)  node $(hostname)"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind "$HOME:$HOME" \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml python -u "$DFLAM/src/analysis/sq2_purity_acrossk.py"
echo "done $(date)"
