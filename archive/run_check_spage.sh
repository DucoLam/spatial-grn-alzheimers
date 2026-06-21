#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=96GB
#SBATCH --time=00:20:00
#SBATCH --job-name=check_spage
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/check_spage_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/check_spage_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/archive/_check_spage.py
