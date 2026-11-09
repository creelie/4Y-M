#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
$PY gauss_k0.py selftest
mkdir -p out; : > out/k0.jsonl
for n in 2 3 4 6 8 10 12 16 20 24; do $PY gauss_k0.py $n 0.5 1 2 4 | tee -a out/k0.jsonl; done
