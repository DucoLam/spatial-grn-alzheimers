#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
LOGS=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs
while squeue -j 12659191,12659192,12659193,12659194 -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 20
done
echo "==== all four SQ2 jobs left the queue ===="
for j in 12659191 12659192 12659193 12659194; do
  echo "----- job $j -----"
  tail -n 6 $LOGS/*_$j.out 2>/dev/null
done
