#!/bin/bash
# Usage: sbatch 08_subsample.sh <cell_type> <prefix>
# Example: sbatch --job-name=sub_l23it 08_subsample.sh 'L2/3 IT' l23_it
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=01:00:00
#SBATCH --job-name=subsample
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/subsample_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/subsample_%j.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
CELL_TYPE="${1}"
PREFIX="${2}"

if [ -z "$CELL_TYPE" ] || [ -z "$PREFIX" ]; then
    echo "ERROR: Usage: sbatch 08_subsample.sh <cell_type> <prefix>"
    exit 1
fi

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "Starting subsample job: cell_type='$CELL_TYPE'  prefix='$PREFIX'"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec   --writable-tmpfs   --bind /tudelft.net:/tudelft.net   --bind $HOME:$HOME   "$SIF"   pixi run --manifest-path /opt/app/pixi.toml   python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/08_subsample.py     --cell_type "$CELL_TYPE"     --prefix "$PREFIX"

echo "Done"
echo "End time: $(date)"
