#!/bin/bash
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=256GB
#SBATCH --time=00:30:00
#SBATCH --job-name=subsample_rna_hvg_only
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm/02_subsample_rna_hvg_only/output.out
#SBATCH --error=slurm/02_subsample_rna_hvg_only/output.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

apptainer exec \
  --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --env PYTHONPATH=/opt/app/src \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/02_subsample_rna_hvg_only.py
