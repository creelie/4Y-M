#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
mkdir -p out; : > out/twoscale.jsonl
for n in 4 8 12 16; do for ml in 1 4; do $PY gauss_twoscale.py $n $ml 16 | tee -a out/twoscale.jsonl; done; done
