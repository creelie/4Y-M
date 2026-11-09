#!/bin/bash
# SU(2)_3 at fixed physical box L*g^2 = 12.8: beta 5,7.5 (re-done at small L), 20, 25.
set -e
cd "$(dirname "$0")"
PY=python
gcc -O2 -o su2bp su2bp.c -lm
mkdir -p outfb; rm -f outfb/*.dat
{
for n in 2 3 4 6; do echo "16 5 $n 5000 5 12 4 $RANDOM outfb/b5_n${n}.dat"; done
for n in 3 4 6 9; do echo "24 7.5 $n 5000 5 12 4 $RANDOM outfb/b7.5_n${n}.dat"; done
for n in 8 12 16 24; do echo "64 20 $n 5000 5 12 4 $RANDOM outfb/b20_n${n}.dat"; done
for n in 10 15 20 30; do echo "80 25 $n 4000 5 12 4 $RANDOM outfb/b25_n${n}.dat"; done
} | xargs -P 16 -L 1 bash -c './su2bp $0 $1 $2 $3 $4 $5 $6 $7 > $8'
$PY analyze.py outfb/*.dat | tee outfb/summary.jsonl
