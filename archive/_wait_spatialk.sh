#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
L=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_spatialk_12659292.out
while squeue -j 12659292 -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 20
done
echo "==== spatial-k sweep 12659292 finished ===="
cat $L 2>/dev/null
echo "---- err ----"
tail -n 6 /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_spatialk_12659292.err 2>/dev/null
