#!/bin/bash
#SBATCH --job-name=hvg_selection
#SBATCH --output=slurm/01_hvg_selection/output.out
#SBATCH --error=slurm/01_hvg_selection/output.err
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
PROJECT=/tudelft.net/staff-umbrella/ScReNI/bsc-screni
SCRIPT=/tudelft.net/staff-umbrella/ScReNI/dflam/bsc-screni
SIF=$PROJECT/container_0-1-3.sif

apptainer exec --writable-tmpfs --pwd /opt/app --containall \
  --bind "$SCRIPT/src/:/opt/app/src/" \
  --bind "$PROJECT/pixi.toml:/opt/app/pixi.toml" \
  --bind "$PROJECT/data/:/opt/app/data/" \
  --bind "$PROJECT/output/:/opt/app/output/" \
  --env PYTHONPATH=/opt/app/src \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -m screni.data.feature_selection
