#!/bin/bash
#SBATCH --job-name=sq4_morans_agg
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=00:20:00
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq4_morans_agg_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq4_morans_agg_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
set -e
DFLAM=/tudelft.net/staff-umbrella/ScReNI/dflam
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "aggregate start $(date)  node $(hostname)"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind "$HOME:$HOME" \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u "$DFLAM/src/analysis/sq4_regulator_morans.py" --aggregate
echo "aggregate done $(date)"
