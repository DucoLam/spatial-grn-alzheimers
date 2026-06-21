#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=03:59:00
#SBATCH --job-name=sq2_analysis
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_analysis_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_analysis_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
set -e
DFLAM=/tudelft.net/staff-umbrella/ScReNI/dflam
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# analysis/<SQ>/<test> directory tree (one dir per SQ, tests under SQ2)
mkdir -p "$DFLAM"/analysis/SQ1 \
         "$DFLAM"/analysis/SQ2/variogram \
         "$DFLAM"/analysis/SQ2/mantel \
         "$DFLAM"/analysis/SQ2/morans_i \
         "$DFLAM"/analysis/SQ2/grn_clustering \
         "$DFLAM"/analysis/SQ3 \
         "$DFLAM"/analysis/SQ4

echo "start $(date)  node $(hostname)"
run () {
  echo "===== $1 ====="
  apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind "$HOME:$HOME" \
    "$SIF" pixi run --manifest-path /opt/app/pixi.toml python -u "$DFLAM/src/analysis/$1"
}
run sq2_variogram.py
run sq2_mantel.py
run sq2_morans_i.py
run sq2_grn_clustering.py
echo "done $(date)"
