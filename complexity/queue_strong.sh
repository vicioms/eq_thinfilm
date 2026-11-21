#!/bin/bash
# strong-pinning runs (k = A > 1), force grid up to A
set -eo pipefail
cd "$(dirname "$0")"
python3 -u run_levels.py --A 2.0 --K 6 8 10 12 14 16 18 20 --nsamp 2000 --fmaxfrac 1.0 | tee log_strong.txt
python3 -u run_levels.py --A 3.0 --K 6 8 10 12 14 16 --nsamp 2000 --fmaxfrac 1.0 | tee -a log_strong.txt
python3 -u run_levels.py --A 6.0 --K 4 5 6 7 8 9 10 --nsamp 1000 --fmaxfrac 1.0 | tee -a log_strong.txt
echo QUEUE_DONE | tee -a log_strong.txt
