#!/bin/bash
# Finite-size-scaling scan on a SLURM cluster: one GPU per array task.
#
# Every array task runs a disjoint slice (--shard) of the same deterministic
# task list and writes $OUT/results_shard<ID>.jsonl. Resubmitting the same
# script skips (L, a, seed) already present in any results*.jsonl, so
# timed-out jobs can simply be resubmitted.
#
#   cd bubble_sim && mkdir -p logs && sbatch slurm_fss.sh
#   python fss_shapes.py --dir runs/fss_D0.04 --no-plot   # on the cluster, then push shapes.npz
#   python fss_analyze.py --dir runs/fss_D0.04             # once (some of) the shards are done
#
# Delta = 0.04 (half the first CPU run): at L = 128 the bubble is pinned with
# h_d ~ 0.03 and no overhangs; the initial bubble (R0 = L/6) collapses at h = 0
# for small L, so h_start is raised in steps of --h-step 0.005 until it is pinned.
#
# Check the GPU backend once before the first submission (on a GPU node):
#   python test_gpu_backend.py
#
# (logs/ must exist before sbatch: SLURM opens the output file before the
#  script runs.)

#SBATCH --job-name=bubble_fss
#SBATCH --partition=graphic
#SBATCH --array=0-15                 # NSHARDS = 16 (keep in sync below)
#SBATCH --gres=gpu:1
#SBATCH --constraint="[a100|h200]"
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/fss_%A_%a.out

set -euo pipefail
NSHARDS=16
DELTA=0.04
OUT=runs/fss_D${DELTA}

# --- environment (edit) ------------------------------------------------------
module load cuda/12.8
source /home/schimmenti/miniconda3/etc/profile.d/conda.sh
conda activate torchcompile
# needs: numpy scipy matplotlib numba numba-cuda  (pip install numba-cuda)
# ---------------------------------------------------------------------------

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
nvidia-smi --query-gpu=name,memory.total --format=csv

python fss_scan.py \
    --backend cuda --gpus 1 \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS" \
    --outdir "$OUT" \
    --aniso 0 0.25 0.5 0.75 1 \
    --L 128 256 512 1024 2048 \
    --samples 128:256 256:128 512:64 1024:32 2048:16 \
    --Delta "$DELTA" --w 3 --burgers-angle 0 \
    --h-start 0 --h-step 0.005 \
    --dh 0.001 --n-bisect 10
