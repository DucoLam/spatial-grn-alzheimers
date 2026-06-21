#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=03:59:00
#SBATCH --job-name=sq2_pcsens
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_pcsens_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_pcsens_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
DFLAM=/tudelft.net/staff-umbrella/ScReNI/dflam
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "start $(date) $(hostname)"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind "$HOME:$HOME" \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml python -u "$DFLAM/src/analysis/sq2_pc_sensitivity.py"
echo "done $(date)"
