#!/bin/bash
#SBATCH --job-name=merf_don
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=00:20:00
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/merf_don_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/merf_don_%j.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
apptainer exec --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --bind $HOME:$HOME \
  /tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/archive/_merfish_donors.py
