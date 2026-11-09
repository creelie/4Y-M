#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
gcc -O3 -march=native -o pathcov pathcov.c -lm
./pathcov selftest
mkdir -p out1
for bw in 0.5 1.5 2.3; do for L in 4 6 8; do
  ./pathcov 4 $L $bw id 200 400 2 $RANDOM out1/b4_L${L}_bw${bw} &
done; done
wait
$PY analyze.py out1
