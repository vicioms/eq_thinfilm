#!/bin/bash
# direct runs (one enumeration per force, f_c by bisection), used by verify_tilt.py
# to check the tilt-symmetry mapping; P = 4, kappa = 1
set -eo pipefail
cd "$(dirname "$0")"
python3 -u run.py --A 0.6 --K 8 12 16 20 24 28 32 --nsamp 2000 --fmax 0.3 --nf 16 | tee log_A0.6.txt
python3 -u run.py --A 0.6 --K 36 40 --nsamp 600 --fmax 0.3 --nf 16 --seed0 0 | tee -a log_A0.6.txt
python3 -u run.py --A 1.0 --K 8 12 16 20 24 28 --nsamp 2000 --fmax 0.5 --nf 16 | tee log_A1.txt
python3 -u run.py --A 0.3 --K 8 16 24 32 40 48 --nsamp 2000 --fmax 0.15 --nf 16 | tee log_A0.3.txt
echo QUEUE_DONE
