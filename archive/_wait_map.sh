#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
JID=12618920
LOG=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/map_l23_${JID}.out
while squeue -j $JID -h -o %T 2>/dev/null | grep -qE 'PENDING|RUNNING|CONFIGURING|COMPLETING'; do
  sleep 10
done
echo "==== job left queue ===="
cat "$LOG" 2>/dev/null
ls -l /tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/grn_clustering/map_l23_H21019_k2468.png 2>/dev/null
