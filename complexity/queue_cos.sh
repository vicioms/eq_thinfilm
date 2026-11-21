#!/bin/bash
# Gaussian (FLRT) disorder, piecewise-linear interpolation (M nodes per period)
set -eo pipefail
cd "$(dirname "$0")"
python3 -u run_levels.py --model cos --A 0.6 --M 16 --K 8 12 16 20 24 28 32 --nsamp 4000 | tee log_cos.txt
python3 -u run_levels.py --model cos --A 0.6 --M 32 --K 8 12 16 20 24 28 --nsamp 2000 | tee -a log_cos.txt
python3 -u run_levels.py --model cos --A 0.6 --M 16 --K 36 40 --nsamp 2000 | tee -a log_cos.txt
echo QUEUE_DONE | tee -a log_cos.txt
