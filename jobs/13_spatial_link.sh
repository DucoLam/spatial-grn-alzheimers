#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --time=06:00:00
#SBATCH --job-name=spatial_link
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/spatial_link_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/spatial_link_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container_0-1-3.sif
echo "start $(date)  node $(hostname)"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/13_spatial_link.py
echo "done $(date)"
