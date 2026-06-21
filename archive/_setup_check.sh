#!/bin/bash
export PYTHONPATH="/tudelft.net/staff-umbrella/ScReNI/dflam/src:$PYTHONPATH"
cd /tudelft.net/staff-umbrella/ScReNI/dflam
sed -i 's/\r$//' run_sq2_analysis.sh
echo "=== syntax check ==="
for f in analysis_common.py sq2_variogram.py sq2_mantel.py sq2_morans_i.py sq2_grn_clustering.py; do
  if python3 -m py_compile "$f"; then echo "ok   $f"; else echo "FAIL $f"; fi
done
echo "=== import check (container) ==="
apptainer exec --writable-tmpfs --bind /tudelft.net:/tudelft.net \
  bsc-screni/container.sif pixi run --manifest-path /opt/app/pixi.toml \
  python /tudelft.net/staff-umbrella/ScReNI/dflam/archive/_check_imports.py
