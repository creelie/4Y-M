#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
$PY gauss_cb.py selftest
mkdir -p out; : > out/cb.jsonl
for cfg in "6 2" "8 2" "6 3" "8 3"; do
  $PY gauss_cb.py $cfg | tee -a out/cb.jsonl
done
