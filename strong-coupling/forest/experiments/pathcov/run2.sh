#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
gcc -O3 -march=native -o pathcov pathcov.c -lm
./pathcov selftest | tail -1
mkdir -p out2
$PY analyze2.py control out2/control.json
for bw in 1.5 2.3; do
  for L in 4 6 8; do ./pathcov 4 $L $bw id 200 3000 2 $RANDOM out2/b4_L${L}_bw${bw} & done
  ./pathcov 6 4 $bw id 200 3000 2 $RANDOM out2/b6_L4_bw${bw} &
done
wait
$PY analyze2.py out2 out2/summary.json
