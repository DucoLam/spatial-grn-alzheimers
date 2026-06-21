#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=3:00:00
#SBATCH --job-name=plot_gene_expression
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm/07_plot_gene_expression/output.out
#SBATCH --error=slurm/07_plot_gene_expression/output.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "Starting plot gene expression job"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Working dir: $(pwd)"
echo "Start time: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec \
  --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --env PYTHONPATH=/opt/app/src:$(pwd)/SpaGE \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/07_plot_gene_expression.py

echo "Plot gene expression job done"
echo "End time: $(date)"
