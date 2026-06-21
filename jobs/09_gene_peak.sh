#!/bin/bash
# Usage: sbatch 09_gene_peak.sh <prefix>
# Example: sbatch 09_gene_peak.sh l23_it
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --job-name=gene_peak
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/gene_peak_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/gene_peak_%j.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
PREFIX="${1}"

if [ -z "$PREFIX" ]; then
    echo "ERROR: Usage: sbatch 09_gene_peak.sh <prefix>"
    exit 1
fi

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "Starting gene_peak job: prefix='$PREFIX'"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec   --writable-tmpfs   --bind /tudelft.net:/tudelft.net   --bind $HOME:$HOME   "$SIF"   bash -c "pip install --quiet --user rdata pyfaidx &&     pixi run --manifest-path /opt/app/pixi.toml     python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/09_gene_peak.py       --prefix $PREFIX"

echo "Done"
echo "End time: $(date)"
