#!/bin/bash
# Finite-size-scaling scan on a SLURM cluster: one GPU per array task.
#
# Every array task runs a disjoint slice (--shard) of the same deterministic
# task list and writes runs/fss/results_shard<ID>.jsonl. Resubmitting the same
# script skips (L, a, seed) already present in any results*.jsonl, so
# timed-out jobs can simply be resubmitted.
#
#   sbatch slurm_fss.sh
#   python fss_analyze.py            # once (some of) the shards are done
#
# Adapt the partition / account / modules to your cluster.

#SBATCH --job-name=bubble_fss
#SBATCH --array=0-15                 # NSHARDS = 16 (keep in sync below)
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fss_%A_%a.out
##SBATCH --partition=gpu
##SBATCH --account=XXX

set -euo pipefail
NSHARDS=16

# --- environment (edit) ------------------------------------------------------
# module load cuda/12.4
# source ~/miniforge3/etc/profile.d/conda.sh && conda activate bubble
# needs: numpy scipy matplotlib numba numba-cuda  (pip install numba-cuda)
# ---------------------------------------------------------------------------

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
nvidia-smi --query-gpu=name,memory.total --format=csv

python fss_scan.py \
    --backend cuda --gpus 1 \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS" \
    --aniso 0 0.25 0.5 0.75 1 \
    --L 256 512 1024 2048 \
    --samples 256:128 512:64 1024:32 2048:16 \
    --Delta 0.08 --w 3 --burgers-angle 0 \
    --dh 0.001 --n-bisect 10
