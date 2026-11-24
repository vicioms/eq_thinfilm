# Coarse-grained front model (front_model.py)

## Method
Lattice u_i (dx=1, kappa=eta=1): m u_tt + u_t = u_xx + f - Phi_i. Phonon energies e_pm (upwind advection at c=1/sqrt(m)),
d_t e_pm +- c d_x e_pm = -e_pm/m + v Phi/2 (sliding sites only), T = e_+ + e_-. Sliding sites: Phi = F(u_t, T_i), closure
F tabulated from sc.Line(m).F on a 25 x 20 (log v in [1e-3,1.5], T in [0,4]) grid, bilinear interpolation (log v, T).
Stuck sites: u_t=0, Phi = u_xx+f, released when |u_xx+f| > f_s(T); sliding sites with |u_t|<v_lo and drive<f_s re-stick.
f_s choices: exp (f_c exp(-T/T0), T0=1 base, 0.3 variant), closure (F(v_lo,T)), const (f_c=0.308). Initial condition:
300 hot sites at closure velocity v(f) (v_min if f<f_min) with T=T(v), rest pinned, T=0; u=0. Front = rightmost sliding site;
s from linear fit over t in [0.2,1]*15m. dt = sqrt(m)/4 (c dt = 1/4); halving dt changes nothing except noise at r<=1.05 (below).
Tables: data/front_table_m*.npz (f_min = 0.1397, 0.0943 reproduced); scans: data/front_scan_m*.npz.

## s/c vs f/f_min (base: v_lo=1e-3, f_s=f_c exp(-T))
| f/fmin | m=256 model | m=1024 model |
|---|---|---|
| 0.95 | nan (all stuck) | -1.32 |
| 1.00 | -1.02 | -1.19 |
| 1.04 | -1.03 | -1.17 |
| 1.06 | -1.14 | -1.10 |
| 1.07 | +0.48 | -1.02 |
| 1.08 | +0.76 | +0.73 |
| 1.10 | +0.86 | +0.85 |
| 1.15 | +0.94 | +0.95 |
| 1.25 | +0.97 | +0.98 |

Measured (statistics3.json): m=256: f/fmin 0.859,0.917,0.96,1.002,1.045,1.11,1.217 -> s/c -5.3,-2.2,-0.32,+0.32,+0.66,+0.79,+0.84 (f0=0.981);
m=1024: 0.848,0.912,0.954,0.997,1.039,1.113,1.272 -> -4.2,-2.2,-0.35,+0.35,+0.76,+0.88,+0.92 (f0=0.975).

## Sensitivity
- f_s(T) (exp T0=1, closure-based, constant f_c): negligible effect (s identical within ~0.02 c). T0=0.3 changes only the retreating branch (-1.0 -> ~-0.1..-0.8 c) and the zero stays at ~1.06-1.07.
- v_lo: large effect on the sub-threshold side. v_lo=1e-2: front frozen (s=0) for f/fmin <= 1.07 (sites re-stick before they can accelerate), sign change not meaningful; v_lo=1e-4 (m=256) gives a spurious slow advance (+0.2..0.4 c) even at 0.9 fmin. m=1024 is insensitive to v_lo 1e-3 vs 1e-4. Saturation branch (f>=1.1 fmin) insensitive to everything.
- Front location s>0 branch is converged in dt; the retreat branch at 0.95-1.05 for m=256 is noisy (dt/2 gives -1.5..-1.0 at 0.95-1.0).

## Verdict
Partial. Positive: (i) saturation is reproduced quantitatively: 0.86/0.85 c at 1.1 fmin, 0.94-0.95 c at 1.15, 0.97-0.98 c at 1.25 (measured 0.84/0.93 and 0.88-0.92 c at 1.11-1.27; model slightly high, ~5%; the saturation is essentially the phonon-heat front, upwind advection of T triggering the next site, so it is set by c); (ii) s changes sign, and robustly (independent of f_s choice).
Negative: (i) the zero is at f0/fmin ~ 1.07 (between 1.06 and 1.08) for both masses, ABOVE f_min, not 0.98 below it; (ii) the transition is a near step (-1.0 c -> +0.5..0.75 c within 0.01-0.02 fmin) not a linear crossing through zero with |s|~0.35 c at f0 +- 0.02 fmin; (iii) the retreating branch is ~ -1 c (cooling front at the sound speed: the hot region simply loses its heat/ sliding at the phonon speed) rather than the measured gentle -0.3 .. -5 c linear growth; (iv) results depend on v_lo for the below-threshold side.
Diagnosis: the model has a sharp heat-front mechanism (advance happens only if the phonon heat arriving at the pinned site can lower f_s / sustain sliding, a threshold, giving a discontinuous s ~ c jump), while the measured s(f) is smooth and crosses zero slightly below f_min, indicating a front whose speed is controlled by a continuous balance (e.g. nucleation/avalanche of the pinned side by the elastic load u_xx, a finite-width, noisy front) that the deterministic uniform-disorder closure (no disorder realisation, T-independent stuck threshold) does not contain. The shift of f0 above f_min is plausibly from the neglected lattice-disorder fluctuations that release stuck sites below the mean-field f_s. Not tuned further.
