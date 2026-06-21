#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=00:20:00
#SBATCH --job-name=plot_grns
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/plot_grns_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/plot_grns_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/12_plot_grns.py
