#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
JOBS="12659502,12659503,12659504,12659505,12659506,12659507,12659508,12659509"
while squeue -j $JOBS -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 60
done
echo "==== all 8 panel-validation jobs left the queue ===="
cd /tudelft.net/staff-umbrella/ScReNI/dflam/data/validation
for f in panel100_astrocyte panel100_l23_it panel100_l4_it panel100_oligo panel500_astrocyte panel500_l23_it panel500_l4_it panel500_oligo; do
  echo "----- ${f} -----"
  grep -E 'argmax|cell_type' ${f}_persection.csv 2>/dev/null || echo "MISSING ${f}_persection.csv"
done
