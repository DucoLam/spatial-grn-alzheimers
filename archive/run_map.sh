#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=00:30:00
#SBATCH --job-name=map_l23
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/map_l23_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/map_l23_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
DFLAM=/tudelft.net/staff-umbrella/ScReNI/dflam
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
echo "start $(date)  node $(hostname)"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind "$HOME:$HOME" \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml python -u "$DFLAM/archive/_map_l23.py"
echo "done $(date)"
