#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=00:10:00
#SBATCH --job-name=insp_sec
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/insp_sec_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/insp_sec_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net   /tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif   pixi run --manifest-path /opt/app/pixi.toml   python -u /tudelft.net/staff-umbrella/ScReNI/dflam/archive/_inspect_sections.py
