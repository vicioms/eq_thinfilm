#!/bin/bash
# Finite-size-scaling scan on a SLURM cluster: one GPU per array task.
#
#   cd bubble_sim && mkdir -p logs
#   sbatch slurm_fss.sh        # R0 = L/6       -> runs/fss_D0.04      (extends the existing run)
#   sbatch slurm_fss.sh 32     # fixed R0 = 32  -> runs/fss_D0.04_R32
#
# then, on the cluster, for each output directory:
#   python fss_shapes.py --dir runs/fss_D0.04 --no-plot --workers 8
#   git add runs/fss_D0.04/shapes.npz runs/fss_D0.04/results_shard*.jsonl
#
# Every array task runs a disjoint slice (--shard) of the same deterministic
# task list and writes $OUT/results_shard<ID>.jsonl. (L, a, seed) already in
# any results*.jsonl are skipped: the L/6 series only runs the new seeds and
# L = 4096, and timed-out jobs can simply be resubmitted.
#
# Delta = 0.04, L = 128 ... 4096. Cost on an A100 (median per sample, R0 = L/6):
# L = 128-1024: 20-75 s, 2048: ~4 min, 4096: ~20 min (estimated) -> ~150 GPU-h
# per series, ~10 h per shard. Configurations are stored as float16:
# ~33 MB per L = 4096 sample, ~10 GB per series in $OUT/configs.
#
# With a fixed R0 the bubble is the same size at every L: h_d then contains
# the same Laplace term at every L, and the box size only changes how much
# disorder the runaway bubble samples. R0 = 32 is pinned below h_d at Delta = 0.04.
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
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/fss_%A_%a.out

set -euo pipefail
NSHARDS=16
DELTA=0.04
R0=${1:-}                              # empty: R0 = L/6
if [ -n "$R0" ]; then
    OUT=runs/fss_D${DELTA}_R${R0}
    R0_ARG=(--r0 "$R0")
else
    OUT=runs/fss_D${DELTA}
    R0_ARG=()
fi
echo "output: $OUT  R0: ${R0:-L/6}"

# --- environment (edit) ------------------------------------------------------
module load cuda/12.8
source /data/biophys/schimmenti/uv_dir/mhn/bin/activate
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
    --L 128 256 512 1024 2048 4096 \
    --samples 128:512 256:512 512:256 1024:128 2048:64 4096:32 \
    --Delta "$DELTA" --w 3 --burgers-angle 0 ${R0_ARG[@]+"${R0_ARG[@]}"} \
    --h-start 0 --h-step 0.005 \
    --dh 0.001 --n-bisect 10
