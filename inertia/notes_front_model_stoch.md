# Stochastic block front model (front_model_stoch.py)

## Model (function `sim2`)
Line coarse-grained into blocks of size ell = c*m = sqrt(m) (heat-attenuation length; 8 at m=64, 16 at m=256; periodic). Per block: u, v, phonon energies e_+, e_-, T = e_+ + e_-.
- f_i = f + K (u_{i-1} + u_{i+1} - 2 u_i), K = kappa/ell^2 = 1/m.
- m dv = [f_i - v - F(v,T)] dt + sqrt(2 T dt / nu) xi,  v >= 0.   (F = sc.Line(m).F, tabulated 40x41; noise from the block's own T)
- d e_pm = [-/+ (c/ell)(e_pm,i - e_pm,i-+1) - e_pm/m + v F/2] dt   (power balance T = m v F in a uniform hot state; upwind phonon transport at c).
- Stuck = v ~ 0: F(v->0) is large, so a block stays stuck unless noise kicks v over the unstable root. Release by active neighbours (invasion) comes from (i) the heat they send in (raises T, so noise) and (ii) the elastic load K du. No rate is imposed by hand; no input from (i).
- Active block: v > 0.3 v_min-ish (0.3 * initial hot v). dt = 0.02 m.
Check: uniform hot state reproduces the closure (f_min = 0.2006 vs 0.2008, v(f_min) = 0.07).
Parameters: ell = sqrt(m), K = 1/m, dt are fixed by the physics. One free parameter: nu (independent modes per block, noise variance T^2-scale/nu; I scanned 4, 10, 20, 60; nu = 4 used for the main runs). Closure tables are exact inputs.
Dropped earlier variants (not shipped as results): T-only model without u (front pinned by lattice, zero > 1.1 f_min); quenched block force offsets (hot state destroyed by heterogeneity).

## Results (m = 64 unless noted; r = f/f_min; s/c from the slope of the active-block count, 3 seeds, N = 3000 blocks)
| r | nu=4 | nu=10 | nu=20 | nu=60 | m=256, nu=4 |
|---|---|---|---|---|---|
| 0.96 | -12 | -9..-14 | collapse | collapse | -4..-8 |
| 1.00 | -6.5 | -22 | -17..-37 | -36..-48 | -7 |
| 1.03 | ~0 (+0.5,-0.2,-0.6) | -13 | -12 | -0.9 | ~0 (+0.5,-0.2,-0.5) |
| 1.05 | ~0 (0.5,0,0) | -1.3 | -0.8 | -0.84 | ~0.1..0.6 |
| 1.07 | 0.1..0.6 | -0.4 | -0.8 | -0.8 | +0.6 |
| 1.09 | +0.67 | ~0 | ~-0.5 | -0.8 | +0.7 |
| 1.12 | +0.70 | +0.70 | +0.7/-0.6 | ~0 | +0.7 |
| 1.20 | +0.78 | +0.78 | +0.77 | +0.76 | +0.78 |
(large negative values are the whole hot domain dying, not a front speed.)

Bulk, all-hot start, m=64, nu=4, N=4096 blocks, to 3000 m: rho(0)/late rho: r=0.98: 0.44 -> 0 within 50m; 1.00: 0.70 -> 0 (half-time 75m); 1.02: half-time 275m; 1.03: 725m; 1.04: 2025m (still decaying at 3000m); 1.06: stationary 0.905; 1.08: 0.924. m=256: 1.00 dies (50m); 1.03: half-time 1950m; 1.05: 0.927; 1.07: 0.943.

## Verdict
(1) Edge-speed zero: at nu=4 it sits at 1.03-1.05 f_min (m=64 and 256 alike); at nu>=10 at 1.08-1.12. Noise moves it down but never below f_min. Measured: 0.975-0.98. FAIL, but the zero moves the right way as noise grows; for nu < 4 the hot interior dies even above f_min.
(2) Patchy phase: BOTH sides are wrong in detail. Below the threshold (1.02-1.04 f_min) the activity does decay slowly, with half-times that grow steeply (275m, 725m, 2025m) as the threshold ~1.045 is approached from below. That is a long-lived phase, but the decay is a sigmoid (linear shrinking of hot domains at |s|~ 0.5-1 c) rather than exponential, and sits above f_min, not at 0.974 f_min. Above the threshold the stationary activity is ~0.9, not 0.4, and there are no patches (uniform hot state with a few cold blocks).
(3) Transition: clearly discontinuous: rho jumps 0 -> 0.9 over 0.01-0.02 f_min, front speed jumps from ~0 to ~0.7 c. No DP behaviour (no check of the t^-0.16 law is meaningful: the active state never has low density). Consistent with a bistable hot/cold system where s changes sign at a Maxwell-like point; the measured rho = 0.4 at 0.989 f_min indicates a different mixed (patchy) steady state that the model does not have.
Saturation 0.78 c at 1.2 f_min (measured 0.84-0.93).
What the model captures: closure-consistent hot state, edge speed saturating at ~c, long-lived metastable sliding with divergent lifetime as the threshold is approached from below, noise lowering the zero. What it misses: the zero below f_min, patchy steady state, DP-like continuous onset. The likely missing ingredient is quenched heterogeneity that keeps hot patches alive below f_min (one-loop F(v->0, T=0) being much larger than the true f_c may matter too); my quenched-offset attempt was crude and destroyed the hot state.
Figure: fig_front_stoch.png (edge speed, rho(t), late-time activity). Data: data/stoch_*.npz.
