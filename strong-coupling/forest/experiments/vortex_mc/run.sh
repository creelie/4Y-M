#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
python vmc.py 4 3000
