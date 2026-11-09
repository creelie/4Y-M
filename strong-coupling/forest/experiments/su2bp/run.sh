#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
gcc -O2 -o su2bp su2bp.c -lm
mkdir -p out; rm -f out/*.dat
{
echo "6 5 6 400 3 8 3 1 out/selftest_L6n6.dat"
echo "32 7.5 6 800 5 12 10 7 out/thin_b7.5_n6_s10.dat"
for n in 2 3 4 6; do echo "32 5 $n 1200 5 12 4 $RANDOM out/b5_n${n}.dat"; done
for n in 3 4 6 9; do echo "32 7.5 $n 1200 5 12 4 $RANDOM out/b7.5_n${n}.dat"; done
for n in 4 6 8 12; do echo "32 10 $n 1200 5 12 4 $RANDOM out/b10_n${n}.dat"; done
} | xargs -P 14 -L 1 bash -c './su2bp $0 $1 $2 $3 $4 $5 $6 $7 > $8'
$PY analyze.py out/*.dat | tee out/summary.jsonl
