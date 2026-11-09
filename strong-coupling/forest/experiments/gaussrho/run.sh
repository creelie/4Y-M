#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
$PY gauss_rho.py selftest
mkdir -p out; : > out/rho.jsonl
for cfg in "3 6 3 0,0,0,0" "4 6 4 0,0,0,0" "3 6 6 0,0,0,0" "3 6 6 1,1,1,1" "3 8 6 0,0,0,0" "4 6 8 0,0,0,0" "4 6 8 2,2,2,2" "3 8 9 0,0,0,0"; do
  $PY gauss_rho.py $cfg | tee -a out/rho.jsonl
done
