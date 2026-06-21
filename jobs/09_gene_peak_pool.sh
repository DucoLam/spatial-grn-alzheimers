#!/bin/bash
# Usage: sbatch 09_gene_peak_pool.sh <prefix>   (e.g. pool4)
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=96GB
#SBATCH --time=24:00:00
#SBATCH --job-name=gene_peak_pool
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/gene_peak_pool_%j.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/gene_peak_pool_%j.err
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
PREFIX="${1:-pool4}"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container_0-1-3.sif   # has MOODS + rdata + pyfaidx natively
echo "09_gene_peak start: $(date)  prefix=$PREFIX  job $SLURM_JOB_ID"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net --bind $HOME:$HOME \
  "$SIF" pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/09_gene_peak.py --prefix "$PREFIX"
echo "Done: $(date)"
