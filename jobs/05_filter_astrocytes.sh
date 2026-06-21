#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=128GB
#SBATCH --time=00:30:00
#SBATCH --job-name=filter_astrocytes
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm/05_filter_astrocytes/output.out
#SBATCH --error=slurm/05_filter_astrocytes/output.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

apptainer exec \
  --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --env PYTHONPATH=/opt/app/src \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/05_filter_astrocytes.py
