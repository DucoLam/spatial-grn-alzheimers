#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --job-name=wscreni_astrocytes
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/10_wscreni_astrocytes/output.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/10_wscreni_astrocytes/output.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
mkdir -p /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/10_wscreni_astrocytes

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "Starting wscreni_astrocytes job"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec \
  --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --bind $HOME:$HOME \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/10_wscreni_astrocytes.py

echo "Done"
echo "End time: $(date)"
