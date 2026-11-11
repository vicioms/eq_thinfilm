#!/bin/bash
# Finite-size-scaling scan on a SLURM cluster, CPU version: each array task
# runs its shard with one single-threaded worker per allocated core.
#
# Every array task runs a disjoint slice (--shard) of the same deterministic
# task list and writes runs/fss_cpu/results_shard<ID>.jsonl. Resubmitting the
# same script skips (L, a, seed) already present in any results*.jsonl, so
# timed-out jobs can simply be resubmitted.
#
# Results go to runs/fss_cpu (not runs/fss) so that the CPU and GPU arrays can
# run at the same time without writing to the same shard files.
#
#   cd bubble_sim && mkdir -p logs && sbatch slurm_fss_cpu.sh
#   python fss_analyze.py --dir runs/fss_cpu     # once (some of) the shards are done
#
# (logs/ must exist before sbatch: SLURM opens the output file before the
#  script runs.)
#
# Cost (one core, Delta = 0.08): L=64 ~10 s, 128 ~35 s, 256 ~3.5 min, 512 ~25 min
# per sample, i.e. ~95 core-hours for the whole list below.

#SBATCH --job-name=bubble_fss_cpu
#SBATCH --partition=medium             # your CPU partition (default partition if omitted)
#SBATCH --array=0-15                 # NSHARDS = 16 (keep in sync below)
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/fss_cpu_%A_%a.out

set -euo pipefail
NSHARDS=16

# --- environment (edit) ------------------------------------------------------
source /home/schimmenti/miniconda3/etc/profile.d/conda.sh
conda activate torchcompile
# needs: numpy scipy matplotlib numba
# ---------------------------------------------------------------------------

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

python fss_scan.py \
    --backend cpu --workers "$SLURM_CPUS_PER_TASK" \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS" \
    --outdir runs/fss_cpu \
    --aniso 0 0.25 0.5 0.75 1 \
    --L 64 128 256 512 \
    --samples 64:256 128:128 256:64 512:32 \
    --Delta 0.08 --w 3 --burgers-angle 0 \
    --dh 0.001 --n-bisect 10
