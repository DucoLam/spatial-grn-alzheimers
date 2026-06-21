#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --job-name=gene_peak_astrocytes
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/09_gene_peak_astrocytes/output.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/09_gene_peak_astrocytes/output.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
mkdir -p /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/09_gene_peak_astrocytes

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container_0-1-3.sif

echo "Starting gene_peak_astrocytes job"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Deps (rdata, pyfaidx, MOODS) are native in container_0-1-3.sif

apptainer exec \
  --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --bind $HOME:$HOME \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/09_gene_peak_astrocytes.py

echo "Done"
echo "End time: $(date)"
