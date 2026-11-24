# Hybrid closure (sc_hybrid.py)

## Definition
Fixed v, kappa = eta = A = 1. Unknown B(tau). Given B:
1. C_F(tau) = Delta_B(v tau) (Gaussian-smeared cubic B-spline), S_F(w) = 2 int C_F cos(w tau) dtau.
2. Shape: S0(w) = S_F(w) Im chi_loc(w)/(eta w)   (driven spectrum, as in sc_noneq.py).
3. F_pin(B) = -int R(tau) Delta'_{B(tau)}(v tau) dtau (same as sc.fpin, T = 1 with B in place of T b).
4. Amplitude: mu = v F_pin(B) / P2[S0], P2[S] = eta int dw/pi w^2 S(w).
5. B_new(tau) = (2/pi) int dw mu S0(w) (1 - cos w tau).
Fixed point B = B_new (damped iteration, alpha = 0.5, start from B = 0 at the largest v, then continuation downward in v; tolerance 1e-6 relative on B for v tau < 40). At the fixed point eta <du'^2> = eta int dw/pi w^2 (mu S0) = v F_pin exactly (P2/P1 = 1 to 1e-15 by construction), so B is the driven-spectrum shape with the equilibrium-closure-style power-balance amplitude. mu = 1 recovers sc_noneq. At large v, B -> 0, P2 -> linresp.power, mu -> 1 (checked: mu = 1.05 at v = 0.3, m = 256; f_hyb = f_eq = 0.3077 to 4 digits at v = 0.3, and the two coincide within 1e-3 for v >~ 0.1).
Grid: 25 velocities geomspace(0.005, 0.3), dtau = 0.1, tau_max = max(40 m, 1e4), pad 4. m = 4096 not run (not cheap).
Converged everywhere (max residual 1e-6, 21-34 iterations). One branch only (from B = 0 with continuation); other fixed points not searched.

## Flow curves
| m | f_min eq. closure (v_min) | f_hyb | measured threshold (0.97-0.98 f_min_eq) |
|---|---|---|---|
| 256 | 0.1396 (0.046) | no minimum: monotone, 0.1173 at v = 0.005 (grid edge) | ~0.136-0.137 |
| 1024 | 0.0943 (0.033) | 0.0867 at v = 0.0113 (shallow min, 0.0869 at 0.0087-0.0140) | ~0.092 |
For m = 256 extra points at lower v (vmin = 0.001): f = 0.1194, 0.1175, 0.1164, 0.1154 at v = 0.01, 0.005, 0.0025, 0.00125: still slowly decreasing ~1e-3 per halving (a log-like creep), no clean threshold, value ~0.115 and below. So the hybrid threshold is ~0.84 f_min_eq (m = 256) and ~0.92 f_min_eq (m = 1024), i.e. 14-16% and 6-8% BELOW the measured thresholds, whereas the equilibrium closure is 2-3% above.

Sliding branch vs measured (mean v over seeds, data/flow_m*.npz), f_meas / f_hyb / f_eq:
- m = 256: v = 0.0068: 0.1266 / 0.1182 / 0.229; v = 0.055: 0.1413 / 0.1367 / 0.1404; v = 0.091: 0.1560 / 0.1557 / 0.1552; v >= 0.12 all three agree to < 1e-3 (3 decimals).
- m = 1024: v = 0.0132: 0.0897 / 0.0868 / 0.1083; v = 0.061: 0.1053 / 0.1050 / 0.1046; v >= 0.09 all agree to 1e-3.
Hybrid is better than eq. below v_min (eq. rises steeply, as it must since T ~ v^2 and no self-sustained noise; the measured curve is nearly flat, which the hybrid reproduces within 3-7%), equivalent for v >~ 0.06, and slightly worse at v near v_min for m = 256 (0.1367 vs 0.1413 measured; eq. 0.1404).

## B(tau) vs measured (closure_check, v = 0.1158 for m = 256, 0.0899 for m = 1024), ratio to measured
| v tau | m=256 hyb | m=256 eq | m=1024 hyb | m=1024 eq |
|---|---|---|---|---|
| <1 | 1.57 | 1.54 | 1.63 | 1.61 |
| 1-3 | 1.48 | 1.39 | 1.59 | 1.55 |
| 3-10 | 1.21 | 0.89 | 1.35 | 1.20 |
| 10-30 | 1.01 | 0.62 | 1.14 | 0.87 |
| 30-100 | 0.84 | 0.49 | 0.93 | 0.68 |
| >100 | 0.83 | 0.48 | 0.80 | 0.57 |
Long-lag shape is much better than the eq. closure (0.8-1.0 vs 0.5-0.7), but the short-lag excess (1.5-1.6x) is NOT cured and is slightly worse than sc_noneq (1.3 / 1.16), because power balance forces mu = 1.2-1.4 > 1 at these velocities (B_drv alone would carry too little power). The measured B(0) ratio puzzle (~0.65 T b) remains. Fp at the measured runs: 0.0541 (meas 0.0542) for m = 256; 0.0312 (meas 0.0301) for m = 1024.

## Energy consistency
By construction P2/P1 = 1 at every v. But constructing the amplitude by balance does not remove the self-sustained noise: as v decreases, mu falls ~ v (0.093 at v = 0.005, m = 256; 0.027 at 0.00125) while Bmax stays 500-1800 (m = 256 / 1024) and F_pin does not go to 0 (F_pin ~ 0.11 at v -> 0). The state is hot: P1 = v F_pin -> 0 is balanced by P2 -> 0 only through a tiny mu multiplying a shape whose absolute scale is large (S_F is ~ saturated because B is huge but the shape itself carries a w^2-moment that is v-independent). So the hybrid trades the violated energy balance (sc_noneq, P2/P1 up to 100) for a spurious nonzero zero-velocity pinning force: a finite friction as v -> 0 (apparent threshold ~0.115 for m = 256, below true f0 = 0.137), not a pinned state.

## Verdict
Mixed; overall NOT better as a threshold theory, and not a clean replacement for the equilibrium closure. Better: the energy balance is exact (fixes sc_noneq), long-lag B(tau) shape (0.8-1.0 vs 0.5-0.7 of measured), and sliding flow curve below v_min (flat, within 3-7% vs eq. +80%/+20% off). Worse: threshold 6-16% too low (eq. is 2-3% too high), no sharp f_min/v_min for m = 256, short-lag B still 1.5x too large, and the low-v state is a hot, self-sustained one (F_pin(v -> 0) finite, f slowly creeping down). Equivalent to eq. closure for v >~ 0.1 (3 decimals). Caveats: single branch from B = 0 with continuation; only m = 256, 1024; grid min of m = 256 at the edge v = 0.005.
Data: data/sc_hybrid.npz (vs; Fp_, f_, P1_, P2_, mu_, it_, err_, Bmax_, f_eq_, hyb_at_meas_, eq_at_meas_, meas_v_, meas_f_, fmin_, vmin_, cmp_* per m).
