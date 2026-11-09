#!/bin/bash
# Higher-statistics SU(2)_3 rerun: 5000 configs, fourth beta=15 (L=48 keeps L*g^2=12.8 as at beta=10).
set -e
cd "$(dirname "$0")"
PY=python
gcc -O2 -o su2bp su2bp.c -lm
mkdir -p outhi; rm -f outhi/*.dat
{
echo "6 5 6 1500 3 8 3 11 outhi/selftest_L6n6.dat"
for n in 2 3 4 6; do echo "32 5 $n 5000 5 12 4 $RANDOM outhi/b5_n${n}.dat"; done
for n in 3 4 6 9; do echo "32 7.5 $n 5000 5 12 4 $RANDOM outhi/b7.5_n${n}.dat"; done
for n in 4 6 8 12; do echo "32 10 $n 5000 5 12 4 $RANDOM outhi/b10_n${n}.dat"; done
for n in 6 9 12 18; do echo "48 15 $n 5000 5 12 4 $RANDOM outhi/b15_n${n}.dat"; done
} | xargs -P 17 -L 1 bash -c './su2bp $0 $1 $2 $3 $4 $5 $6 $7 > $8'
$PY analyze.py outhi/*.dat | tee outhi/summary.jsonl
