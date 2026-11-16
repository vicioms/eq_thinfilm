#!/bin/bash
# Flat-wall depinning scan (KFR strip geometry) on a SLURM cluster, one GPU per array task.
#
#   cd bubble_sim && mkdir -p logs
#   sbatch slurm_fss_flat.sh 0      # b along x: wall line || b  (screw wall) -> runs/flat_D0.04_b0_A2
#   sbatch slurm_fss_flat.sh 90     # b along y: wall line _|_ b (edge wall)  -> runs/flat_D0.04_b90_A2
#   sbatch slurm_fss_flat.sh 90 4   # optional 2nd argument: aspect, box L x (aspect L), default 2
#
# The wall and the box never rotate (wall along x, motion along y); the angle
# rotates b in the noise spectrum only.
# Any other angle (e.g. 80, 85 for edge walls tilted by 10 or 5 deg) runs a = 1
# only, for the tilt dependence of the threshold (Tang-Kardar-Dhar: around the
# hard direction F_c(s) - F_c(0) ~ -|s|^1.57 for qKPZ).
#
# At a = 0 the noise is isotropic and the same seed gives the same (rotated)
# realization for every angle, so a = 0 is run only with angle 0.
#
# then, on the cluster, for each output directory:
#   python analyze_flat.py --extract runs/flat_D0.04_b0_A2 --workers 8
#   git add runs/flat_D0.04_b0_A2/walls.npz runs/flat_D0.04_b0_A2/results_shard*.jsonl
# and anywhere:
#   python analyze_flat.py runs/flat_D0.04_b0_A2 runs/flat_D0.04_b90_A2
#
# Cost, measured (angle 0, aspect 2, first shards): ~30 s per sample at
# L = 256, ~55 s at 512, ~85 s at 1024, ~11-13 min at 2048; about 36 GPU-h
# for angle 0 (3 values of a) and ~25 GPU-h for 90 (2 values of a). Cost
# grows about linearly with the aspect.
# Restarting skips finished (L, a, seed).

#SBATCH --job-name=flat_fss
#SBATCH --partition=graphic
#SBATCH --array=0-15                 # NSHARDS = 16 (keep in sync below)
#SBATCH --gres=gpu:1
#SBATCH --constraint="[a100|h200]"
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/flat_%A_%a.out

set -euo pipefail
NSHARDS=16
DELTA=0.04
ANGLE=${1:?usage: sbatch slurm_fss_flat.sh <burgers angle: 0 | 90 | ...> [aspect]}
ASPECT=${2:-2}
OUT=runs/flat_D${DELTA}_b${ANGLE}_A${ASPECT}
case "$ANGLE" in
    0)  ANISO=(0 0.5 1) ;;
    90) ANISO=(0.5 1) ;;
    *)  ANISO=(1) ;;
esac
echo "output: $OUT  burgers angle: $ANGLE  aspect: $ASPECT  a: ${ANISO[*]}"

# --- environment (edit) ------------------------------------------------------
module load cuda/12.8
source /data/biophys/schimmenti/uv_dir/mhn/bin/activate
# needs: numpy scipy matplotlib numba numba-cuda contourpy
# ---------------------------------------------------------------------------

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs
nvidia-smi --query-gpu=name,memory.total --format=csv

python fss_scan.py \
    --geometry flat --y0-frac 0.1 --aspect "$ASPECT" \
    --backend cuda --gpus 1 \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS" \
    --outdir "$OUT" \
    --aniso "${ANISO[@]}" \
    --L 256 512 1024 2048 \
    --samples 256:256 512:128 1024:64 2048:32 \
    --Delta "$DELTA" --w 3 --burgers-angle "$ANGLE" \
    --h-start 0 --h-step 0.005 \
    --dh 0.001 --n-bisect 10
