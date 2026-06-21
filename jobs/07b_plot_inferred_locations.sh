#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB
#SBATCH --time=00:30:00
#SBATCH --job-name=plot_inferred_140
#SBATCH --mail-user=dflam@student.tudelft.nl
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/07b_plot_inferred_140/output.out
#SBATCH --error=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/07b_plot_inferred_140/output.err

export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
mkdir -p /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/07b_plot_inferred_140
mkdir -p /tudelft.net/staff-umbrella/ScReNI/dflam/plots/07b_inferred_locations_140

SIF=/tudelft.net/staff-umbrella/ScReNI/bsc-screni/container.sif

echo "Starting plot_inferred_locations job"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start time: $(date)"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

apptainer exec \
  --writable-tmpfs \
  --bind /tudelft.net:/tudelft.net \
  --bind $HOME:$HOME \
  "$SIF" \
  pixi run --manifest-path /opt/app/pixi.toml \
  python -u /tudelft.net/staff-umbrella/ScReNI/dflam/src/pipeline/07b_plot_inferred_locations.py

echo "Done"
echo "End time: $(date)"
