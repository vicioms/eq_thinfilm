#!/bin/bash
# L = 8192 bubbles (R0 = L/6), added to the existing series in runs/fss_D0.04.
#
#   cd bubble_sim && mkdir -p logs
#   sbatch slurm_fss_8192.sh
#
# Same protocol as slurm_fss.sh (Delta = 0.04, w = 3, b along x, h ramp from 0,
# dh = 0.001, 10 bisections), so the analysis treats L = 8192 as one more size.
# Purpose: at a = 1 the edge flanks of the critical bubble advance only 9, 13, 17
# sites from R0 at L = 1024, 2048, 4096 (the screw tips 204, 413, 759); this
# checks how the flank advance, the threshold gap to edge walls and the
# large-scale zeta continue one size up. a = 0 is the isotropic control.
#
# 3 a x 16 seeds = 48 tasks, one per array task (NSHARDS = 48). (L, a, seed)
# already in any results*.jsonl are skipped, so timed-out tasks can be
# resubmitted as they are. Results are appended to results_shard<ID>.jsonl of the
# series (the L <= 4096 jobs are finished, so no two jobs write the same file).
#
# Cost on an A100, extrapolated from 2048 -> 4096 (x5.7 at a = 1, x7.7 at a = 0):
# median per sample ~2 h (a = 1) to ~4.5 h (a = 0), worst ~10 h -> ~150 GPU-h.
# Memory: 8192^2 float64 fields ~0.5 GB each on the GPU; configurations are
# stored as float16, ~130 MB per sample, ~6 GB in total in $OUT/configs.
# --max-sweeps 8e6: a relaxation that hits the cap counts as depinned (biasing
# h_d low); none hit 2e6 up to L = 4096, the cap only matters if one runs long.
#
# then, on the cluster (re-extracts the whole series, all L):
#   python analyze_bubbles.py --extract runs/fss_D0.04 --workers 8
#   git add runs/fss_D0.04/contours.npz runs/fss_D0.04/results_shard*.jsonl
# and anywhere:
#   python analyze_bubbles.py runs/fss_D0.04
#
# (logs/ must exist before sbatch: SLURM opens the output file before the
#  script runs.)

#SBATCH --job-name=bubble_8192
#SBATCH --partition=graphic
#SBATCH --array=0-47                 # NSHARDS = 48 (keep in sync below)
#SBATCH --gres=gpu:1
#SBATCH --constraint="[a100|h200]"
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/fss8192_%A_%a.out

set -euo pipefail
NSHARDS=48
DELTA=0.04
OUT=runs/fss_D${DELTA}
echo "output: $OUT  L = 8192  R0: L/6"

# --- environment (edit) ------------------------------------------------------
module load cuda/12.8
source /data/biophys/schimmenti/uv_dir/mhn/bin/activate
# ---------------------------------------------------------------------------

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
nvidia-smi --query-gpu=name,memory.total --format=csv

python fss_scan.py \
    --backend cuda --gpus 1 \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS" \
    --outdir "$OUT" \
    --aniso 0 0.5 1 \
    --L 8192 \
    --samples 8192:16 \
    --Delta "$DELTA" --w 3 --burgers-angle 0 \
    --h-start 0 --h-step 0.005 \
    --dh 0.001 --n-bisect 10 \
    --max-sweeps 8000000
