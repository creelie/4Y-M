#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
$PY gauss_massU.py selftest
mkdir -p out; : > out/massU_t.jsonl
for ml in 0.5 1 2 4; do for n in 2 3 4 6 8 10; do $PY gauss_massU.py $n $ml 9 t | tee -a out/massU_t.jsonl; done; done
