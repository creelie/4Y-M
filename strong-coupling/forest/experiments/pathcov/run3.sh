#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
gcc -O3 -march=native -o pathcov pathcov.c -lm
./pathcov selftest | tail -1
mkdir -p out3
for L in 6 8; do
  for bw in 3.0 4.0; do ./pathcov 4 $L $bw id 200 3000 2 $RANDOM out3/b4_L${L}_bw${bw}_id & done
  ./pathcov 4 $L 2.3 haar 200 3000 2 $RANDOM out3/b4_L${L}_bw2.3_haar &
done
wait
$PY analyze2.py out3 out3/summary.json
