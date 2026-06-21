#!/bin/bash
# Sharded wScReNI: sbatch 10_wscreni_shard.sh <prefix>   (16 shards in parallel)
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB
#SBATCH --time=12:00:00
#SBATCH --job-name=wscreni_shard
#SBATCH --array=0-15
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/wscreni_shard_%A_%a.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/wscreni_shard_%A_%a.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
PREFIX="${1:-pool4}"
NSHARDS=16
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container_0-1-3.sif
echo "wscreni shard $SLURM_ARRAY_TASK_ID/$NSHARDS  prefix=$PREFIX  $(date)  node $(hostname)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/10_wscreni.py \
    --prefix "$PREFIX" --n_jobs 8 --n_shards $NSHARDS --shard_id $SLURM_ARRAY_TASK_ID
echo "shard $SLURM_ARRAY_TASK_ID done: $(date)"
