#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
while squeue -j 12659586,12659587 -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 45
done
echo "==== SQ2-pcsens (12659586) + SQ3 (12659587) finished ===="
echo "----- pc_sensitivity -----"
cat /tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ2/pc_sensitivity/pc_sensitivity.csv 2>/dev/null
echo "----- SQ3 correlations -----"
cat /tudelft.net/staff-umbrella/ScReNI/dflam/analysis/SQ3/sq3_correlations.csv 2>/dev/null
echo "----- sq3 err tail -----"
tail -n 6 /tudelft.net/staff-umbrella/ScReNI/dflam/slurm/logs/sq3_12659587.err 2>/dev/null
