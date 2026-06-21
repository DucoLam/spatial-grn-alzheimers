#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
while squeue -j 12659517,12659655,12659505,12659509 -h -o %T 2>/dev/null | grep -qE 'PEND|RUN|CONFIG|COMPLETING'; do
  sleep 60
done
echo "==== final batch (oligo panels, oligo cross-donor, composition) finished ===="
cd /tudelft.net/staff-umbrella/ScReNI/dflam
echo "--- oligo gene panels (argmax) ---"
grep -H argmax data/validation/panel100_oligo_persection.csv data/validation/panel500_oligo_persection.csv 2>/dev/null
echo "--- oligo cross-donor ---"
cat data/validation/crossdonor_oligo.csv 2>/dev/null
echo "--- composition corr (significant after FDR) ---"
awk -F, 'NR==1 || $7<0.1' analysis/SQ3/sq3_composition_corr.csv 2>/dev/null | head -40
echo "--- composition files ---"
ls -1 analysis/SQ3/ | grep composition
