#!/bin/bash
PY=python
S="$(dirname "$0")/linresp.py"
export OMP_NUM_THREADS=1
$PY $S --beta 2.3 --L 4 --T 8 --therm 5 --meas 40 || exit 1   # smoke
for b in 2.0 2.2 2.3 2.4 2.5 2.6; do
  $PY $S --beta $b --L 8 --T 16 --meas 4000 --out res_$b.json > log_$b.txt 2>&1 &
done
wait; cat res_*.json
