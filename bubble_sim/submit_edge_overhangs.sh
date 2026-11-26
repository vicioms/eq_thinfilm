#!/bin/bash
# Edge-wall overhang runs proposed in the phase-field wall note ("Overhangs at depinning in
# phase-field walls"): a = 1, b at 90 deg (edge wall), w = 3, Ly = 2 Lx, full configurations saved.
#
#   cd bubble_sim && bash submit_edge_overhangs.sh           # submits the three jobs below
#   bash submit_edge_overhangs.sh --dry-run                 # prints the sbatch commands only
#
#   runs/flat_D0.08_b90_A2   Delta = 0.08, L = 512 (64 walls), 1024 (32 walls)
#   runs/flat_D0.16_b90_A2   Delta = 0.16, L = 512 (64 walls), 1024 (32 walls)
#   runs/flat_D0.04_b90_A2   Delta = 0.04, adds L = 4096 (32 walls) to the existing edge series
#
# If sizeable overhangs appear at Delta = 0.08 or 0.16, extend to L = 2048 and 4096 with e.g.
#   SAMPLES="2048:32 4096:32" ANISO=1 sbatch slurm_fss_flat.sh 90 2 0.08
# Restarting skips finished (L, a, seed), so the same command resumes a job.
#
# Cost: ~4 h per wall at L = 4096, Delta = 0.04 (8x the 2048 time), i.e. 2 walls per shard of 16,
# within the 1-day limit of slurm_fss_flat.sh. The Delta = 0.08 and 0.16 runs have half the walls of the
# Delta = 0.04 edge series up to L = 1024 (~30 min per shard there), so they should be short.
#
# After the jobs, on the cluster, for each output directory:
#   python analyze_flat.py --extract runs/flat_D0.08_b90_A2 --workers 8
#   git add runs/flat_D0.08_b90_A2/walls.npz runs/flat_D0.08_b90_A2/results_shard*.jsonl
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
run=(sbatch)
[ "${1:-}" = "--dry-run" ] && run=(echo sbatch)

SAMPLES="512:64 1024:32" ANISO=1 "${run[@]}" slurm_fss_flat.sh 90 2 0.08
SAMPLES="512:64 1024:32" ANISO=1 "${run[@]}" slurm_fss_flat.sh 90 2 0.16
SAMPLES="4096:32"        ANISO=1 "${run[@]}" slurm_fss_flat.sh 90 2 0.04
