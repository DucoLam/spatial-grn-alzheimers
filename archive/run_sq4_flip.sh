#!/bin/bash
#SBATCH --job-name=sq4_cosdiag
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB
#SBATCH --time=00:30:00
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq4_cosdiag_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq4_cosdiag_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  /tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/archive/_sq4_flip_check.py
