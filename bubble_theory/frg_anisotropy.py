"""Threshold anisotropy of flat walls: FRG predictions vs the flat runs (notes_frg.md).

  1. Scale-dependent anisotropy of the force felt by a straight wall segment of length ell,
         R(ell) = [Var_ell(edge) - Var_ell(screw)] / [Var_ell(edge) + Var_ell(screw)] / a,
     computed from the noise spectrum of bubble_depinning.dislocation_spectrum (R -> 1 for ell >> w).
  2. Straight-wall FRG scaling h_c(alpha) = A(alpha)^(2/3) h_c(A = 1), A = 1 + a cos 2 alpha.
  3. Exact symmetries: h_edge(a) = h_screw(-a), h_45(a) = h_45(-a); even/odd decomposition.
  4. One-parameter form h = h0(L) [1 + a C(L) cos 2 alpha]^(2/3): C from edge and screw, a = 0.5 and 1.
  5. FSS amplitudes B(alpha, a) of <h_d> = h_inf + B L^(-3/4).
  6. TKD crossover length of edge walls (a = 1) with the measured, scale-dependent C(ell).
  7. Direct test of orientation smearing: <cos 2 theta> of the critical walls' tangents.

Usage (from the repo root):
  python3 bubble_theory/frg_anisotropy.py [--runs bubble_sim/runs] [--A A2]
"""
import argparse
import collections
import glob
import json
import os
import sys

import numpy as np
from scipy.integrate import dblquad
from scipy.optimize import curve_fit

W_NOISE = 3.0          # noise core width w (fss_scan default)
W_WALL = np.sqrt(2.0)  # wall width: Phi = tanh(zeta / sqrt 2) for c = eps0 = 1
LS = [256, 512, 1024, 2048]


# ---------------------------------------------------------------- 1. R(ell)
def seg_variance(ell, edge):
    """Variance of the force on a straight segment of length ell, a = 1 spectrum, up to a constant.
    S ~ (k.b)^2/|k| exp(-k^2 w^2); edge: b || normal; screw: b || tangent.
    The wall profile filters k_n by rho_hat; the segment average filters k_t by sinc^2."""
    def rho_hat(kn):
        x = 0.5 * np.pi * kn * W_WALL
        return 1.0 if abs(x) < 1e-9 else x / np.sinh(x)

    def f(kn, kt):
        k = np.hypot(kn, kt)
        if k == 0:
            return 0.0
        kb = kn if edge else kt
        s = np.sinc(kt * ell / (2 * np.pi))  # np.sinc(x) = sin(pi x)/(pi x)
        return kb ** 2 / k * np.exp(-(k * W_NOISE) ** 2) * rho_hat(kn) ** 2 * s ** 2

    K = 4.0 / W_NOISE
    return dblquad(f, 0, K, 0, K, epsabs=1e-12)[0]


def R_of_ell(ells):
    out = []
    for ell in ells:
        e, s = seg_variance(ell, True), seg_variance(ell, False)
        out.append((e - s) / (e + s))
    return np.array(out)


# ---------------------------------------------------------------- data
def load_thresholds(runs, A):
    H = {}
    for b in ("b0", "b45", "b90"):
        d = collections.defaultdict(list)
        for f in glob.glob(os.path.join(runs, f"flat_D0.04_{b}_{A}", "results*.jsonl")):
            for line in open(f):
                r = json.loads(line)
                d[(r["aniso"], r["L"])].append(r["h_d"])
        for (a, L), v in d.items():
            H[(b, a, L)] = (np.mean(v), np.std(v, ddof=1) / np.sqrt(len(v)))
    return H


def tangent_cos2(curve, ell):
    xy = np.asarray(curve, float)
    ds = np.hypot(*np.diff(xy, axis=0).T)
    s = np.r_[0, np.cumsum(ds)]
    S = np.arange(0, s[-1], 0.5)
    x, y = np.interp(S, s, xy[:, 0]), np.interp(S, s, xy[:, 1])
    k = max(1, int(round(ell / 0.5)))
    return np.mean(np.cos(2 * np.arctan2(y[k:] - y[:-k], x[k:] - x[:-k])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="bubble_sim/runs")
    ap.add_argument("--A", default="A2", help="aspect tag of the flat runs (A2 or A4)")
    ap.add_argument("--no-walls", action="store_true", help="skip the tangent measurement (needs walls.npz)")
    args = ap.parse_args()

    print("== 1. anisotropy of the force on a straight segment, R(ell) (-> 1 for ell >> w = 3)")
    ells = [1, 3, 6, 12, 24, 48, 96]
    for ell, r in zip(ells, R_of_ell(ells)):
        print(f"   ell = {ell:3d}: R = {r:.3f}")

    H = load_thresholds(args.runs, args.A)
    g = lambda b, a, L: H[(b, a, L)][0]
    print(f"\n== 2. straight-wall FRG scaling vs data ({args.A}); h0 = isotropic (a = 0) threshold at the same L")
    print("   L     h0      a   edge/h0 [pred (1+a)^2/3]   screw/h0 [pred (1-a)^2/3]   45/h0 [pred 1]")
    for L in LS:
        h0 = g("b0", 0.0, L)
        for a in (0.5, 1.0):
            s45 = f"{g('b45', a, L) / h0:.3f}" if ("b45", a, L) in H else "  -  "
            print(f"   {L:4d}  {h0:.5f}  {a:.1f}  {g('b90', a, L) / h0:.3f} [{(1 + a) ** (2 / 3):.3f}]"
                  f"             {g('b0', a, L) / h0:.3f} [{(1 - a) ** (2 / 3):.3f}]              {s45}")

    print("\n== 3. exact symmetries: odd part (edge - screw)/2, even part (edge + screw)/2 - h0, 45° offset h45 - h0")
    for L in LS:
        h0 = g("b0", 0.0, L)
        parts = []
        for a in (0.5, 1.0):
            e, s = g("b90", a, L), g("b0", a, L)
            parts.append(f"a={a}: odd {0.5 * (e - s):+.5f} even {0.5 * (e + s) - h0:+.5f}")
        d45 = g("b45", 1.0, L) - h0
        print(f"   L={L:4d}  " + "  ".join(parts) + f"  | a=1: h45-h0 {d45:+.5f}")

    print("\n== 4. one-parameter form h = h0 [1 + a C cos 2alpha]^(2/3): C from each run")
    print("   L      C_edge(a=.5) C_screw(a=.5) C_edge(a=1) C_screw(a=1)")
    Cs = []
    for L in LS:
        h0 = g("b0", 0.0, L)
        row = []
        for a in (0.5, 1.0):
            row += [((g("b90", a, L) / h0) ** 1.5 - 1) / a, (1 - (g("b0", a, L) / h0) ** 1.5) / a]
        Cs.append(np.mean(row))
        print(f"   {L:4d}  " + "  ".join(f"{c:11.3f}" for c in row))
    Cs = np.array(Cs)
    y, lnA = np.polyfit(np.log(LS), np.log(Cs), 1)
    print(f"   mean C(L) ~ L^{y:.2f}")

    print("\n== 5. FSS amplitudes, <h_d> = h_inf + B L^(-3/4)  (B - B(a=0) should be odd and ~linear in a cos 2alpha)")
    f = lambda L, hi, B: hi + B * np.asarray(L, float) ** -0.75
    fits = {}
    for key in [("b0", 0.0), ("b0", 0.5), ("b0", 1.0), ("b90", 0.5), ("b90", 1.0), ("b45", 1.0)]:
        yv = np.array([H[(key[0], key[1], L)][0] for L in LS])
        ev = np.array([H[(key[0], key[1], L)][1] for L in LS])
        p, cov = curve_fit(f, LS, yv, sigma=ev, absolute_sigma=True, p0=[yv[-1], 0.1])
        chi2 = np.sum(((yv - f(LS, *p)) / ev) ** 2)
        fits[key] = p
        print(f"   {key[0]:>3s} a={key[1]:.1f}: h_inf = {p[0]:.5f}±{np.sqrt(cov[0, 0]):.5f}  B = {p[1]:+.4f}  chi2/dof = {chi2 / 2:.1f}")
    B0 = fits[("b0", 0.0)][1]
    for a in (0.5, 1.0):
        print(f"   a={a}: B_edge - B0 = {fits[('b90', a)][1] - B0:+.4f},  B_screw - B0 = {fits[('b0', a)][1] - B0:+.4f}")
    print(f"   a=1.0: B_45 - B0 = {fits[('b45', 1.0)][1] - B0:+.4f}  (even in a: expected ~a^2)")

    print("\n== 6. TKD crossover for edge walls (a = 1), with the measured threshold anisotropy C(ell)")
    print("   g(ell) = |h_c''(ell)| r_f (ell/L_c)^zeta / (c/mu),  h_c''/h_c = -(8/3) a C / (1 + a C),")
    print("   C(ell) = min(R(ell), fitted power law of section 4): the threshold anisotropy of a segment")
    print("   cannot exceed the anisotropy of the force it feels (power law extrapolated below ell = 256)")
    cmu, rf, Lc, zeta, hc = 0.47, 3.0, 6.0, 1.25, fits[("b0", 0.0)][0]
    grid = np.array([6, 12, 24, 48, 96, 192, 384])
    Cg = np.minimum(R_of_ell(grid), np.exp(lnA) * grid ** y)
    gl = (8 / 3) * Cg / (1 + Cg) * hc * rf * (grid / Lc) ** zeta / cmu
    for ell, c, gg in zip(grid, Cg, gl):
        print(f"   ell = {ell:4d}: C = {c:.3f}  g = {gg:.2f}")
    print("   (g scales as r_f L_c^-zeta: L_c and r_f are Larkin estimates with O(1) prefactors)")
    print(f"   ell_x (g = 1): {np.exp(np.interp(0.0, np.log(gl), np.log(grid))):.0f} sites")

    if args.no_walls:
        return
    print("\n== 7. direct orientation smearing: <cos 2 theta> of the critical walls' tangents (chord over ell = 6)")
    print("   orientation-smearing model: edge needs <cos2theta> = C, screw needs -<cos2theta> = C")
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bubble_sim"))
    from analysis.data import load_flat
    for b, sign in (("b90", +1), ("b0", -1)):
        W = load_flat(os.path.join(args.runs, f"flat_D0.04_{b}_{args.A}"))
        for a in (0.5, 1.0):
            row = []
            for L in LS:
                idx = np.nonzero((W["L"] == L) & (W["aniso"] == a))[0][:64]
                row.append(np.mean([tangent_cos2(W["curve"][i], 6.0) for i in idx]))
            name = "edge " if b == "b90" else "screw"
            print(f"   {name} a={a}: " + "  ".join(f"L={L}: {r:+.3f}" for L, r in zip(LS, row)))


if __name__ == "__main__":
    main()
