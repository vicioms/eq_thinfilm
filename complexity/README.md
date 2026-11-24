# Quenched complexity of a driven elastic line

Working note on the number of equilibria of an elastic chain in a random potential at force f. It compares the sample-averaged ⟨ln N⟩ with the Kac–Rice (annealed) result ln⟨N⟩ of Fyodorov, Le Doussal, Rosso and Texier, *Ann. Phys.* **397**, 1 (2018), [arXiv:1703.10066](https://arxiv.org/abs/1703.10066).

Read [quenched_complexity.pdf](quenched_complexity.pdf).

## Method

Model and counting:
- **Chain:** free ends, zero mass, κ = 1, K sites.
- **Shooting:** equilibria are the zeros of one shooting function. The tangent of the shooting map gives det H (Gelfand–Yaglom) and the Morse index (Sturm).
- **Exact enumeration:** with a piecewise-linear pinning force, the zeros are enumerated exactly by splitting the shooting interval at every kink. Two disorders are used: a triangle wave with a random phase per site, and the piecewise-linear interpolation of FLRT's Gaussian ξ₁cos u + ξ₂sin u.
- **Statistical tilt symmetry:** substitute ũ_i = u_i + f i(i+1)/2. Then N(f) is, in law, the number of crossings of the level fK by the zero-force shooting function W_K(a) = Σ V_i′(u_i), and f_c(K) is distributed as max_a W_K / K. So one f = 0 enumeration gives every force and the threshold.

Annealed quantities:
- **Triangle model:** exact at every K. This covers ⟨N_K(f)⟩ = P p_{S_K}(Kf) ⟨|det(−Δ_N + U)|⟩, the generalised Lyapunov exponent Λ(q), the minima rate Λ_st and the Cramér function I(f).
- **Gaussian model:** the exact asymptotic rate, from a tilted transfer operator.

## Main results

Disorder strengths: triangle wave with L_c = 2.2, 1.4 and 1, and Gaussian disorder with L_c = 1.4.

- **Annealed constants.** C = Λ(1) L_c ≈ 0.456, as in FLRT. The annealed minima constant is C_st ≈ 0.355; FLRT did not compute it.
- **Zero force, fixed amplitude.** For the triangle wave, quenched = annealed within about 1%, for all equilibria and for minima. Typical equilibria carry the tilted stability exponent Λ′(1), not the Lyapunov exponent γ.
- **Zero force, Gaussian disorder.** Σ_q(0) is 2% below annealed, because the random site amplitude is shared by all branches. A phase-annealed, amplitude-quenched operator product predicts the gap with no free parameter: 0.0067, against the measured 0.006–0.008.
- **Finite force.** Quenched = annealed within 1–3% up to f ≈ 0.6–0.7 f_c, with no O(f²) difference. Beyond about 0.75 f_c the quenched complexity falls steeply below the annealed one and vanishes at f_c.
- **Threshold ratio, exactly.** No-passing relaxation gives the exact sample threshold of chains up to K = 1024 at any pinning strength. f_c^ann / f_c = 1.59 ± 0.02, 1.55 and 1.51 at L_c = 2.2, 1.4 and 1. It falls monotonically to the strong-pinning curve (1.17 at k = 6, where it matches the minimax). Extrapolated to the continuum it is 1.60–1.70, consistent with FLRT's dynamical 1.64. Enumeration gives 1.65–1.72 for Gaussian disorder at L_c = 1.4.
- **Threshold fluctuations.** The sample threshold fluctuates as K^(−1/2), as the front of a branching process in a time-random environment does.
- **Mechanism.** The branching structure lives in an environment shared by all branches. The bulk of the population averages that environment over u, and the leading edge does not.
- **Theory.** Σ_q is the Legendre transform of a quenched pressure φ_q(τ), and f_c is its large-τ slope. Conditional Jensen gives exact bounds Σ_q ≤ Σ_sa ≤ Σ_ann. The Jensen gap is O(τ⁴), i.e. o(f²), while the population stays spread in u, and becomes large once it localises (onset ≈ 0.75 f_c). For Gaussian disorder, amplitude quenching lowers the threshold bound from 0.486 to 0.423.
- **Torus transfer operators.** The same population dynamics is a product of random positive operators on the torus (u, w) mod P, which can be iterated at K = ∞. It reproduces the exact quenched pressure through the onset. It shows that the gap is the one-step Jensen term, i.e. the phase variance of the population's growth. It also confirms f_c and the jump as the grid is refined.
- **Replicated Kac–Rice.** ⟨Z^n⟩ of the counting partition function is an exact n-replica transfer operator on the offsets between replicas. Locked replicas close the second-moment window at x coth x = e^{Λ(1)} (x = τA), which is f_ann/√2 in the continuum. Below that, the pair operator, Var ln Z and the measured gap satisfy the replica Gaussian closure Λ − φ_q = σ²/2 = (Λ₂ − 2Λ)/2 to within errors, including the lattice gap at f = 0. The closure gives only a third of the drop from f_ann to f_c, and one-step RSB with independent clusters gives none. So f_ann/f_c comes from correlations between clusters of near-critical strands that meet again in the shared environment.
- **Rigorous bounds below the annealed threshold.** Fractional moments with a change of measure give φ_q ≤ Λ_Q + d(P‖Q) for any tilted disorder law Q. Tilting the law of the phase second differences, for which Λ_Q is an exact transfer operator, proves f_c < f_ann. The bound is 0.1074, 0.2621 and 0.4887 at A = 0.3, 0.6 and 1 with one step of memory, closing 2.5–20% of the gap. Local tilts converge slowly at weak pinning.
- **Block tilts.** Penalising each block of ℓ phases by its own quenched partition function Z_b (computed exactly over polygons of entry states) turns the bound into φ_q ≤ K⁻¹ ln⟨Z ∏ Z_b⁻¹⟩ + ⟨ln Z_b⟩/ℓ. The first term is the growth rate of planted strands. At ℓ = 16 (7–16 L_c) this closes about half of the measured pressure gap at weak pinning and up to 80% at strong pinning, against 7–35% for single-site tilts. So the gap lives in phase correlations over about ten Larkin lengths. For f_c alone it gives f_c ≤ 0.104, 0.239 and 0.444 (10%, 33% and 45% of the way from f_ann to the exact f_c) at A = 0.3, 0.6, 1; as τ → ∞ it reduces to the subadditive block-maximum bound, so the relaxation thresholds stay the accurate route to f_c. The bound holds for every ℓ, but the planted estimator has a tail ~e^{γℓ} and visibly fails at ℓ = 32; only population-stable points (ℓ ≤ 16, τA ≤ 3) are reported.
- **Jump at threshold.** The quenched complexity does not vanish at f_c. At a force O(1/K) below each sample's threshold there are still e^{cK} equilibria, with c ≈ Λ(1)/2, and exponentially many minima.
- **Strong pinning (solvable).** As k → ∞, c → ln 2, because each site can sit on either side of its force peak. Also c_min → 0 and f_c → A − P/2, with sample fluctuations P/√(12K). Enumeration up to k = 6 shows c crossing over from 0.22/L_c toward ln 2.
- **Exact threshold for k > 2.** Minima have every site on a rising branch, so the depinning threshold is a lattice minimax over period choices. It is verified on 540 samples to 3·10⁻⁹. A linear-time dynamic programme extends it to K = 320 and predicts f_c = A − P/2 + 0.90 P/k. The annealed bound overshoots by a constant 0.316 P.
- **What c is.** f_c is the maximal ergodic average of the pinning force, and c is the entropy of the force-maximising invariant measure. In line with the Ledrappier–Young formula, c = d·λ* within 10%. Here λ* is the Lyapunov exponent of the critical orbit, and d ≈ 0.28 is the box-counting dimension of the near-critical set in the shooting parameter; both are measured independently.

## Contents

- `qcount.py`: exact enumeration (numba). `count` runs one force at a time. `count_levels` (triangle) and `count_levels_pl` (general piecewise-linear) give all forces from one f = 0 enumeration.
- `annealed.py`: exact annealed quantities. These are the transfer operators for Λ(q), Λ_st and the tilted Λ(τ), plus the Lyapunov exponent, the Irwin–Hall density and the Cramér function.
- `run_levels.py`: production runs. For each sample it records N(f), N_min(f), ln|det H|, the u-localisation of the equilibria and max W/K. Results go to `data_levels/`, which is git-ignored.
- `queue_levels.sh`, `queue_cos.sh`: the production runs used in the note.
- `run.py`, `queue.sh`: direct runs, one enumeration per force with f_c by bisection, used only to check the tilt mapping.
- `verify_tilt.py`: the exact and statistical checks of the tilt mapping, written to `tilt_check.json`.
- `semi_annealed.py`: the amplitude-quenched rate Λ_q(τ) for Gaussian disorder, written to `semi_annealed_cos.json`.
- `profile_center.py`: the centre and width of the f = 0 level profile, written to `profile_center_A0.6.json`.
- `analyze.py`: writes `statistics.json`, `table_rates.tex`, `table_threshold.tex` and `fig1`–`fig4` (PDF and PNG).
- `torus.py`, `torus_runs.py`: transfer operators on the torus and the runs used in the note, written to `torus.json` and `torus_decomposition.json`.
- `lyap_leader.py`: the stretching rate of the critical configuration, written to `lyap_leader.json`.
- `queue_strong.sh`, `minimax.py`: the strong-pinning runs (k = 2, 3, 6) and the minimax threshold, written to `minimax.json`.
- `strong_exact.py`, `strong_scan.py`: the exact threshold formula for k > 2 (lattice minimax and dynamic programme), its test, and the infinite-chain thresholds, written to `strong_exact.json` and `strong_scan.json`.
- `replica.py`: the pair moment Λ₂(τ) of the replicated Kac–Rice operator by population dynamics (`--check` tests n = 1 against the exact Λ), written to `replica.json`. `theory.py --replica-only` writes `replica_summary.json` and `fig10`.
- `fracbound.py`, `fracbound2.py`: the change-of-measure bounds with iid and Markov tilts of the phase second differences, written to `fracbound.json` and `fracbound2.json`.
- `blockbound.py`: the block-tilt bounds, with Z_b traced exactly over convex polygons and the planted-strand population at two sizes per point (`--check` tests λ = 0 against the exact Λ), written to `blockbound.json`. `theory.py --block-only` writes `block_summary.json` and `fig11`.
- `relax.py`: exact sample thresholds by no-passing (Middleton) forward relaxation plus bisection, K ≤ 1024 at k = 0.3–6, written to `relax.json`. `theory.py --relax-only` extrapolates them (`relax_summary.json`) and draws `fig9`.
- `fractal_top.py`: box counting of the near-critical set in the shooting parameter, written to `fractal_top.json`.
- `theory.py`: computes the quenched pressure from the level counts, its frozen asymptote and the near-threshold rates, and evaluates the f_sa bound. It writes `theory.json`, `table_theory.tex` and `fig5`–`fig11`.
- `main.tex`: the note.

## Reproduce

Requires Python 3, NumPy, SciPy, Numba, mpmath and Matplotlib. The production queues take about 1 h on 10 cores, and the direct runs about 1 h more.

```bash
./queue_levels.sh
./queue_cos.sh
python3 profile_center.py --A 0.6 --K 8 12 16 20 24 28 32
./queue.sh
python3 verify_tilt.py
python3 semi_annealed.py
python3 analyze.py
python3 torus_runs.py
python3 lyap_leader.py
python3 fractal_top.py
./queue_strong.sh
python3 minimax.py
python3 strong_exact.py
python3 strong_scan.py
python3 relax.py
python3 replica.py
python3 fracbound.py
python3 fracbound2.py
python3 blockbound.py
python3 theory.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
cp main.pdf quenched_complexity.pdf
```
