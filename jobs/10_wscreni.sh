#!/bin/bash
# Usage: sbatch 10_wscreni.sh <prefix> [n_jobs]
# Example: sbatch 10_wscreni.sh l23_it 8
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --job-name=wscreni
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/wscreni_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/wscreni_%j.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
PREFIX="${1}"
N_JOBS="${2:-8}"

if [ -z "$PREFIX" ]; then
    echo "ERROR: Usage: sbatch 10_wscreni.sh <prefix> [n_jobs]"
    exit 1
fi

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "Starting wscreni job: prefix='$PREFIX'  n_jobs=$N_JOBS"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec   --writable-tmpfs   --bind /tudelft.net:/tudelft.net   --bind $HOME:$HOME   "$SIF"   pixi run --manifest-path /opt/app/pixi.toml   python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/10_wscreni.py     --prefix "$PREFIX"     --n_jobs $N_JOBS

echo "Done"
echo "End time: $(date)"
