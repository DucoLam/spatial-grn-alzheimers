#!/bin/bash
# Usage: sbatch 10_wscreni_pool.sh <prefix> [n_jobs]   (e.g. pool4 16)
#SBATCH --partition=general
#SBATCH --qos=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --time=4-00:00:00
#SBATCH --job-name=wscreni_pool
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/wscreni_pool_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/wscreni_pool_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
PREFIX="${1:-pool4}"
NJOBS="${2:-16}"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container_0-1-3.sif
echo "10_wscreni start: $(date)  prefix=$PREFIX  n_jobs=$NJOBS  job $SLURM_JOB_ID"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # let n_jobs do the parallelism
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/10_wscreni.py --prefix "$PREFIX" --n_jobs "$NJOBS"
echo "Done: $(date)"
