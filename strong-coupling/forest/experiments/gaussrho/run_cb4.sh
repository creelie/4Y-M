#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
$PY gauss_cb.py selftest
mkdir -p out; : > out/cb4.jsonl
for cfg in "8 4" "9 4"; do
  $PY gauss_cb.py $cfg | tee -a out/cb4.jsonl
done
