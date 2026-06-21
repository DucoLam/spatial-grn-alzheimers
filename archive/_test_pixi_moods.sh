#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=00:20:00
#SBATCH --job-name=test_pxm
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/test_pxm_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/test_pxm_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
WORK=/tudelft.net/staff-umbrella/ScReNI/dflam/_pixi_moods_test
rm -rf "$WORK"; mkdir -p "$WORK"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME "$SIF" bash -lc "
cd '$WORK'
export PIXI_HOME='$WORK/.pixi_home'
pixi init . 2>&1 | tail -2
pixi add python=3.11 2>&1 | tail -3
echo '--- try conda-forge moods ---'
pixi add moods 2>&1 | tail -8 || echo 'conda-forge moods FAILED'
echo '--- try bioconda moods ---'
pixi project channel add bioconda 2>&1 | tail -2
pixi add moods 2>&1 | tail -8 || echo 'bioconda moods FAILED'
echo '--- import test ---'
pixi run python -c 'import MOODS.scan,MOODS.tools; print(\"MOODS OK via pixi\")' 2>&1 | tail -3
"
