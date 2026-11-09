#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
$PY gauss_angle.py selftest
mkdir -p out; : > out/angle_L6.jsonl
for cfg in "3 6 6 0" "4 6 8 0" "3 6 3 0"; do
  $PY gauss_angle.py $cfg | tee -a out/angle_L6.jsonl
done
