#!/bin/bash
set -e
cd "$(dirname "$0")"
PY=python
gcc -O2 -o u1bp u1bp.c -lm
mkdir -p out; rm -f out/*.dat
{
echo "4 1.5 4 600 3 8 3 1 out/selftest_L4n4.dat"
echo "32 1.9 4 1500 5 12 10 7 out/thin_b1.9_n4_s10.dat"
for b in 1.6 1.9 2.2; do for n in 2 3 4 6 8; do echo "32 $b $n 2000 5 12 4 $RANDOM out/b${b}_n${n}.dat"; done; done
} | xargs -P 18 -L 1 bash -c './u1bp $0 $1 $2 $3 $4 $5 $6 $7 > $8'
$PY analyze.py out/*.dat | tee out/summary.jsonl
