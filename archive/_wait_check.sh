#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
while squeue -j 12659322 -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 15
done
echo "==== check_spage 12659322 finished ===="
cat /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/check_spage_12659322.out 2>/dev/null
echo "---- err ----"
tail -n 8 /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/check_spage_12659322.err 2>/dev/null
