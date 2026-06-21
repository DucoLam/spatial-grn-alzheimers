#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
L=/tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_wjac_12659270.out
while squeue -j 12659270 -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 15
done
echo "==== wjaccard job 12659270 finished ===="
cat $L 2>/dev/null
echo "---- err ----"
tail -n 6 /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq2_wjac_12659270.err 2>/dev/null
echo "---- output csvs ----"
ls -l /tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/variogram/variogram_wjaccard* 2>/dev/null
