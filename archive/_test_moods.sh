#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=00:15:00
#SBATCH --job-name=test_moods
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/test_moods_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/test_moods_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
echo '### check toolchain ###'
apptainer exec --bind /tudelft.net:/tudelft.net "$SIF" bash -lc 'which gcc g++ swig cc 2>/dev/null; gcc --version 2>/dev/null | head -1'
echo '### try: swig wheel then MOODS-python ###'
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME   "$SIF" pixi run --manifest-path /opt/app/pixi.toml   bash -c 'pip install --user swig && pip install --user --no-build-isolation MOODS-python' 2>&1 | tail -20
echo '### import test ###'
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME   "$SIF" pixi run --manifest-path /opt/app/pixi.toml   python -c 'import MOODS.scan, MOODS.tools; print("MOODS IMPORT OK")'
