#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=00:15:00
#SBATCH --job-name=plot_cdf
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/plot_cdf_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/plot_cdf_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/06_plot_cdf.py
