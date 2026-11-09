#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
mkdir -p out; : > out/kscan.jsonl
$PY gauss_kscan.py selftest | tee -a out/kscan.jsonl
for n in 8 12 16; do for ml in 1 4; do $PY gauss_kscan.py $n $ml 24 | tee -a out/kscan.jsonl; done; done
