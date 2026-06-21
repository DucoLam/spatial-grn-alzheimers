#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
JOBS="12659514,12659515,12659516,12659517"
while squeue -j $JOBS -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 60
done
echo "==== all 4 cross-donor jobs left the queue ===="
cd /tudelft.net/staff-umbrella/ScReNI/dflam/data/validation
for f in crossdonor_astrocyte crossdonor_l23_it crossdonor_l4_it crossdonor_oligo; do
  echo "----- ${f} -----"
  cat ${f}.csv 2>/dev/null || echo "MISSING ${f}.csv"
done
