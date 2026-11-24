# Non-equilibrium Gaussian closure (sc_noneq.py)

## Method
Fixed v, kappa = eta = A = 1. Iterate (damping 0.5, start B = 0, tolerance 1e-6 on relative change of B for v tau < 40; 22-27 iterations, always converged):
1. C_F(tau) = Delta_{B(tau)}(v tau), exact Gaussian smearing of the cubic B-spline (4th difference of E[(x+sd z)_+^3]).
2. S_F(w) = 2 int C_F cos(w tau) dtau.
3. S_du = S_F Im chi_loc/(eta w).
4. B(tau) = (2/pi) int dw S_du (1 - cos w tau).
All transforms are DCT-IV on a midpoint tau grid (dtau = 0.1, tau_max = max(40 m, 1e4); m = 4096: dtau = 0.2) and a midpoint w grid with zero-padding (pad 4; 2 for m = 4096) to resolve the w^{-3/2} tail.
F_pin = -int R Delta'_{B}(v tau) (sc.fpin, same grid as R). P1 = v F_pin, P2 = eta int dw/pi w^2 S_du.
Two flow curves: f1 = eta v + F_pin (as requested), and f2 = eta v + P2/v (energy-based, the exact identity if P2 were exact).

## Validation (B = 0, no smearing)
P2 vs linresp.power: m = 256, v = 0.02/0.1158/0.3: 6.7187e-3 vs 6.7186e-3, 5.0330e-3 vs 5.0330e-3, 2.1262e-3 vs 2.1262e-3; m = 1024: agreement to 2e-6 relative (4.0337e-3 vs 4.0337e-3 etc.). Error < 0.002%.

## Result 1: the closure does NOT give a threshold from F_pin; the energy check fails at low v
P2/P1 (v = 0.005, 0.0105, 0.0322, 0.1, 0.3 approx; indices 0,10,15,20,25,29 of the grid):
- m = 256: 96.9, 7.7, 2.5, 1.02, 0.78, 0.95
- m = 1024: 48.4, 4.2, 1.44, 0.73, 0.93, 0.99
- m = 4096: 24.3, 2.25, 0.85, 0.83, 0.98, 1.00
At the measured-run velocities P2/P1 = 0.80 (m = 256, v = 0.1158), 0.71 (m = 1024, v = 0.0899), 0.79 (m = 4096, v = 0.080). So the energy balance holds to 5-30% only for v >~ 0.1 and fails badly (factor 2-100) at low v.
Cause: at low v, P2 saturates at a v-independent value (~6.1e-3 for m = 256, close to the radiation limit <F^2>/(2 sqrt(kappa m)) = 6.9e-3): the converged B is huge (B_max ~ 4000 at v = 0.005) and the force seen by a site decorrelates through the line's own fluctuations, i.e. a self-sustained noise that the equilibrium closure lacks. But the pinning-force expression (single-R, one-loop) does not see this power: P1 = v F_pin stays tiny.
Consequently f1(v) = eta v + F_pin is monotonically increasing in v on the whole grid; its "minimum" is at the grid edge v = 0.005 (f1 = 0.0176, 0.0171, 0.0160 for m = 256, 1024, 4096; F_pin(0.005) = 0.0126, 0.0121, 0.0110), no threshold in the usual sense. I regard f1 at v < ~0.05 as unreliable.

## Result 2: f2 = eta v + P2/v has a minimum
| m | eq. closure f_min | non-eq f2_min (v_min) | measured f0 | non-eq f1 at measured v vs measured F_pin |
|---|---|---|---|---|
| 256 | 0.1396 | 0.1526 (0.084) | 0.137 | 0.0561 vs 0.0542 (+3.5%) |
| 1024 | 0.0943 | 0.1022 (0.055) | 0.092 | 0.0312 vs 0.0301 (+3.6%) |
| 4096 | 0.0626 | 0.0674 (0.042) | - | 0.00869 vs 0.01004 (-13%) |
f2_min overshoots the eq. closure and measurement by ~8-11% (grid-edge values of the min; log grid spacing ~12%, minimum is shallow). So the non-eq closure raises the threshold in the direction the B-rescaling test suggested (+10%), but the overshoot is in the wrong direction relative to the simulated f0, and f2 relies on P2 which is only consistent at v >~ 0.1, near v_min (0.04-0.08) the P2/P1 ratio is 0.8-1.4 depending on m. Treat as indicative only.

## B(tau) comparison (B_noneq / B_measured, data/closure_check_m*_f*.npz)
| v tau range | m = 256 (v = 0.1158) | m = 1024 (v = 0.0899) | m = 4096 (v = 0.0800) |
|---|---|---|---|
| < 1 | 1.31 | 1.16 | 1.04 |
| 1-3 | 1.23 | 1.13 | 1.03 |
| 3-10 | 0.98 | 0.94 | 0.94 |
| 10-30 | 0.80 | 0.75 | 0.75 |
| 30-100 | 0.66 | 0.61 | 0.69 |
| > 100 | 0.65 | 0.52 | 0.56 |
For comparison, the equilibrium closure B_cl/B_meas at the first / last tau points is 1.55 / 0.48 (256), 1.62 / 0.56 (1024), 1.54 / 0.72 (4096). The non-eq closure removes most of the short-time excess (1.55 -> 1.3, 1.6 -> 1.16, 1.5 -> 1.04) and moves the shape toward the data, but it still underestimates the long-time B by a factor ~0.5-0.65 (measured B is 1.5-2 times larger), so the long-time (anomalous) growth is not captured. Arrays: cmp_* keys in data/sc_noneq.npz.

## Caveats
- Sites are treated as independent (no spatial correlation of the force on neighbours; the true measured B includes collective effects), and the closure for F_pin uses the same R, not a B-dependent vertex.
- P2 != P1 at low v is a genuine inconsistency of the closure, not numerical: P2 is stable w.r.t. dtau, padding, tau_max (checked at m = 256 with tau_max = 15000 vs 1e4: same values to 3 digits); B_max exceeding 100 at low v makes C_F decorrelate within v tau << 1.
- f_min values are grid values (30 points, geomspace 0.005-0.3), not interpolated.
- Fixed point found by iteration from B = 0 only; other fixed points not searched for.
Data: data/sc_noneq.npz (vs; f1_, Fp_, P1_, P2_, f2_, it_, err_ per m; cmp_* with the measured runs).
