#!/bin/bash
# production runs for the note (level method, all forces per enumeration)
set -eo pipefail
cd "$(dirname "$0")"
python3 -u run_levels.py --A 0.6 --K 8 12 16 20 24 28 32 36 40 44 --nsamp 4000 | tee log_levels.txt
python3 -u run_levels.py --A 1.0 --K 8 12 16 20 24 28 32 --nsamp 4000 | tee -a log_levels.txt
python3 -u run_levels.py --A 0.3 --K 8 16 24 32 40 48 56 --nsamp 4000 | tee -a log_levels.txt
python3 -u run_levels.py --A 0.6 --K 48 --nsamp 2000 | tee -a log_levels.txt
python3 -u run_levels.py --A 1.0 --K 36 --nsamp 2000 | tee -a log_levels.txt
python3 -u run_levels.py --A 0.3 --K 64 --nsamp 2000 | tee -a log_levels.txt
echo QUEUE_DONE | tee -a log_levels.txt
