#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=00:15:00
#SBATCH --job-name=test_c013
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/test_c013_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/test_c013_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container_0-1-3.sif
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME   "$SIF" pixi run --manifest-path /opt/app/pixi.toml   python -c 'import MOODS.scan,MOODS.tools,pyfaidx,rdata,anndata,scipy,sklearn; print("ALL IMPORTS OK; MOODS", getattr(MOODS,"__version__","?"))'
