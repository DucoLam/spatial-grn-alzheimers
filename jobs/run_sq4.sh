#!/bin/bash
#SBATCH --job-name=sq4_morans
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=03:59:00
#SBATCH --array=0-49%10
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq4_morans_%A_%a.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq4_morans_%A_%a.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
set -e
DFLAM=/tudelft.net/staff-umbrella/ScReNI/dflam
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif
N_SHARDS=50          # MUST equal the array width (0-49 => 50)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p "$DFLAM"/analysis/SQ4/shards

echo "start $(date)  node $(hostname)  shard ${SLURM_ARRAY_TASK_ID}/${N_SHARDS}"
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind "$HOME:$HOME" \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u "$DFLAM/src/analysis/sq4_regulator_morans.py" \
  --shard "$SLURM_ARRAY_TASK_ID" --n-shards "$N_SHARDS"
echo "done $(date)  shard ${SLURM_ARRAY_TASK_ID}"
