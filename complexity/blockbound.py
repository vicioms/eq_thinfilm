"""Block coarse-grained change of measure: upper bounds on phi_q(tau) and f_c
for the triangle model (P = 4, slope k = A), with tilts that see whole blocks.

Cut the chain into blocks of l consecutive second differences V (iid uniform
under P).  For any positive block function Z_b(V_block) take
    dQ/dP = prod_b Z_b^{-lam} / <Z_b^{-lam}> .
Then <ln Z>_P <= ln <Z>_Q + D(P||Q) becomes, per site,
    phi_q(tau) <= B = K^-1 ln < Z prod_b Z_b^{-lam} >_P + (lam / l) <ln Z_b>_P ,
any multiplicative constant in Z_b cancelling.  The bound holds for every
choice of (l, lam, Z_b); only the Monte Carlo error of B is not controlled.

Z_b (production, nlines = -1): the quenched partition function of the block
alone, P^-2 int d theta_0 d theta_-1 e^{taub sum_j g(theta_j)} ||M||_F, computed
exactly: the block map is affine on convex polygons of entry states, which are
traced and integrated in closed form (log_zb_poly).  Earlier variants kept for
comparison: a fixed grid of M entry states (nlines = 0) and exact tracing along
nlines free-end shooting lines (nlines > 0).  Pruning (delta) drops pieces below
e^{-delta} of the running total; it is deterministic, so Z_b stays a function of V.

First term: by stationarity <Z F(V)> = P <e^{tau W(0)} |psi(0)| F(V)>, and along
the strand a = 0 the theta_i are iid uniform, the V's being functions of them.
So it is the growth rate of a population of planted strands: theta_i drawn from
the tilted law e^{tau g}, weight = tangent growth, and at each block end the
factor Z_b(V_block)^{-lam}.  lam = 0 gives back the exact Lambda(tau).
Caveat: the block weight / Z_b is heavy-tailed (up to ~e^{gamma l}), so the
population estimate is biased low unless N >> e^{gamma l}; production runs use
two population sizes per point and theory.py keeps only stable points.

    python3 blockbound.py --check   (lam = 0 against the exact Lambda)
    python3 blockbound.py [--l 4,8,16] [--A 0.3,0.6,1] [--tauA 1,2,3,4] [--N 2000,8000]  -> blockbound.json
"""
import json
import time
import os
import sys
import numpy as np
from numba import njit, prange

P = 4.0


@njit(cache=True, inline="always")
def _g(x, A):
    y = x / (0.5 * P)
    n = np.floor(y)
    t = y - n
    if (int(n) & 1) == 0:
        return -A + 2.0 * A * t, 1.0
    return A - 2.0 * A * t, -1.0


@njit(cache=True)
def log_zb(V, A, k, taub, tan, e0, e1):
    """ln Z_b(V) - taub*A*l: log-mean over the M entry states."""
    M = e0.shape[0]
    l = V.shape[0]
    lw = np.empty(M)
    for m in range(M):
        th = e0[m]; thm = e1[m]
        # two tangent vectors (Frobenius norm of the block product)
        p1 = 1.0; q1 = 0.0; p2 = 0.0; q2 = 1.0
        acc = 0.0
        for j in range(l):
            g0, s0 = _g(th, A)
            tn = (2.0 * th - thm + g0 - V[j]) % P
            thm = th; th = tn
            gj, sj = _g(th, A)
            acc += taub * (gj - A)
            if tan:
                c = 2.0 + k * sj
                a1 = c * p1 - q1; q1 = p1; p1 = a1
                a2 = c * p2 - q2; q2 = p2; p2 = a2
                if j % 8 == 7:
                    nr = np.sqrt(p1 * p1 + q1 * q1 + p2 * p2 + q2 * q2)
                    p1 /= nr; q1 /= nr; p2 /= nr; q2 /= nr
                    acc += np.log(nr)
        if tan:
            acc += 0.5 * np.log(p1 * p1 + q1 * q1 + p2 * p2 + q2 * q2)
        lw[m] = acc
    mx = lw.max()
    s = 0.0
    for m in range(M):
        s += np.exp(lw[m] - mx)
    return mx + np.log(s / M)


@njit(cache=True)
def _split(a, th, thp, S, n, out_a, out_th, out_thp, out_S):
    """Insert vertices where th crosses a kink of g (multiples of P/2);
    every field is affine in a between consecutive vertices."""
    h = 0.5 * P
    m = 0
    cap = out_a.shape[0]
    for i in range(n - 1):
        if m >= cap - 2:
            return -1
        out_a[m] = a[i]; out_th[m] = th[i]; out_thp[m] = thp[i]; out_S[m] = S[i]; m += 1
        t0 = th[i]; t1 = th[i + 1]
        if t1 > t0:
            c = np.floor(t0 / h) + 1.0
            while c * h < t1:
                if m >= cap - 2:
                    return -1
                r = (c * h - t0) / (t1 - t0)
                out_a[m] = a[i] + r * (a[i + 1] - a[i]); out_th[m] = c * h
                out_thp[m] = thp[i] + r * (thp[i + 1] - thp[i]); out_S[m] = S[i] + r * (S[i + 1] - S[i])
                m += 1; c += 1.0
        elif t1 < t0:
            c = np.ceil(t0 / h) - 1.0
            while c * h > t1:
                if m >= cap - 2:
                    return -1
                r = (c * h - t0) / (t1 - t0)
                out_a[m] = a[i] + r * (a[i + 1] - a[i]); out_th[m] = c * h
                out_thp[m] = thp[i] + r * (thp[i + 1] - thp[i]); out_S[m] = S[i] + r * (S[i + 1] - S[i])
                m += 1; c -= 1.0
    out_a[m] = a[n - 1]; out_th[m] = th[n - 1]; out_thp[m] = thp[n - 1]; out_S[m] = S[n - 1]; m += 1
    return m


@njit(cache=True)
def _gval(x, A):
    y = x / (0.5 * P)
    n = np.floor(y)
    t = y - n
    if (int(n) & 1) == 0:
        return -A + 2.0 * A * t
    return A - 2.0 * A * t


@njit(cache=True)
def log_zb_curve(V, A, taub, nlines, cap):
    """Exact block partition function on free-end shooting lines:
    entry (theta_0, theta_-1) = (a + x, a), a in [0, P), x on nlines offsets;
    Z_b = mean_x int da e^{taub sum_{j=1..l} (g(theta_j) - A)} |d theta_l / da|.
    The map is piecewise affine, so the image curves are traced exactly."""
    l = V.shape[0]
    a = np.empty(cap); th = np.empty(cap); thp = np.empty(cap); S = np.empty(cap)
    b_a = np.empty(cap); b_th = np.empty(cap); b_thp = np.empty(cap); b_S = np.empty(cap)
    logs = np.empty(nlines)
    for il in range(nlines):
        x = (il + 0.5) * P / nlines
        a[0] = 0.0; th[0] = x; thp[0] = 0.0; S[0] = 0.0
        a[1] = P; th[1] = x + P; thp[1] = P; S[1] = 0.0
        n = 2
        bad = False
        for j in range(l):
            n = _split(a, th, thp, S, n, b_a, b_th, b_thp, b_S)
            if n < 0:
                bad = True; break
            # map: theta_{j+1} = 2 theta_j - theta_{j-1} + g(theta_j) - V_{j+1}
            for i in range(n):
                t = b_th[i]
                tn = 2.0 * t - b_thp[i] + _gval(t, A) - V[j]
                b_thp[i] = t; b_th[i] = tn
            n = _split(b_a, b_th, b_thp, b_S, n, a, th, thp, S)
            if n < 0:
                bad = True; break
            for i in range(n):
                S[i] += taub * (_gval(th[i], A) - A)
        if bad:
            logs[il] = -1e300
            continue
        mx = -1e300
        for i in range(n):
            if S[i] > mx:
                mx = S[i]
        tot = 0.0
        for i in range(n - 1):
            d = abs(th[i + 1] - th[i])
            s0 = S[i] - mx; s1 = S[i + 1] - mx
            ds = s1 - s0
            if abs(ds) < 1e-9:
                e = np.exp(0.5 * (s0 + s1))
            else:
                e = (np.exp(s1) - np.exp(s0)) / ds
            tot += d * e
        logs[il] = mx + np.log(tot / P)
    mx = logs.max()
    s = 0.0
    for il in range(nlines):
        s += np.exp(logs[il] - mx)
    return mx + np.log(s / nlines)


@njit(cache=True)
def _split_p(a, th, thp, S, al, n, oa, oth, othp, oS, oal):
    """As _split, with a live flag per segment (al[i]: segment i -> i+1).
    Dead segments are not refined, and runs of dead segments are merged."""
    h = 0.5 * P
    m = 0
    cap = oa.shape[0]
    for i in range(n - 1):
        if m >= cap - 2:
            return -1
        if al[i] == 0 and m > 0 and oal[m - 1] == 0:
            continue  # extend the previous dead segment
        oa[m] = a[i]; oth[m] = th[i]; othp[m] = thp[i]; oS[m] = S[i]; oal[m] = al[i]; m += 1
        if al[i] == 0:
            continue
        t0 = th[i]; t1 = th[i + 1]
        if t1 != t0:
            if t1 > t0:
                c = np.floor(t0 / h) + 1.0; dc = 1.0
            else:
                c = np.ceil(t0 / h) - 1.0; dc = -1.0
            while (c * h - t1) * dc < 0.0:
                if m >= cap - 2:
                    return -1
                r = (c * h - t0) / (t1 - t0)
                oa[m] = a[i] + r * (a[i + 1] - a[i]); oth[m] = c * h
                othp[m] = thp[i] + r * (thp[i + 1] - thp[i]); oS[m] = S[i] + r * (S[i + 1] - S[i])
                oal[m] = 1; m += 1; c += dc
    oa[m] = a[n - 1]; oth[m] = th[n - 1]; othp[m] = thp[n - 1]; oS[m] = S[n - 1]; oal[m] = 0; m += 1
    return m


@njit(cache=True)
def log_zb_prune(V, A, taub, nlines, cap, delta):
    """log_zb_curve with deterministic pruning: after each site a segment is
    dropped when its partial contribution |(d theta_j, d theta_j-1)| e^{max S} is below e^{-delta}
    times the running total.  Still a fixed function of V."""
    l = V.shape[0]
    a = np.empty(cap); th = np.empty(cap); thp = np.empty(cap); S = np.empty(cap); al = np.zeros(cap, np.int8)
    b_a = np.empty(cap); b_th = np.empty(cap); b_thp = np.empty(cap); b_S = np.empty(cap); b_al = np.zeros(cap, np.int8)
    pot = np.empty(cap)
    logs = np.empty(nlines)
    for il in range(nlines):
        x = (il + 0.5) * P / nlines
        a[0] = 0.0; th[0] = x; thp[0] = 0.0; S[0] = 0.0; al[0] = 1
        a[1] = P; th[1] = x + P; thp[1] = P; S[1] = 0.0; al[1] = 0
        n = 2
        bad = False
        for j in range(l):
            n = _split_p(a, th, thp, S, al, n, b_a, b_th, b_thp, b_S, b_al)
            if n < 0:
                bad = True; break
            for i in range(n):
                t = b_th[i]
                tn = 2.0 * t - b_thp[i] + _gval(t, A) - V[j]
                b_thp[i] = t; b_th[i] = tn
            n = _split_p(b_a, b_th, b_thp, b_S, b_al, n, a, th, thp, S, al)
            if n < 0:
                bad = True; break
            for i in range(n):
                S[i] += taub * (_gval(th[i], A) - A)
            if delta > 0.0:
                mx = -1e300
                for i in range(n - 1):
                    if al[i]:
                        d = np.hypot(th[i + 1] - th[i], thp[i + 1] - thp[i])
                        pot[i] = max(S[i], S[i + 1]) + (np.log(d) if d > 0 else -1e300)
                        if pot[i] > mx:
                            mx = pot[i]
                tot = 0.0
                for i in range(n - 1):
                    if al[i]:
                        tot += np.exp(pot[i] - mx)
                thr = mx + np.log(tot) - delta
                for i in range(n - 1):
                    if al[i] and pot[i] < thr:
                        al[i] = 0
        if bad:
            logs[il] = -1e300
            continue
        mx = -1e300
        for i in range(n - 1):
            if al[i]:
                if S[i] > mx:
                    mx = S[i]
                if S[i + 1] > mx:
                    mx = S[i + 1]
        tot = 0.0
        for i in range(n - 1):
            if al[i] == 0:
                continue
            d = abs(th[i + 1] - th[i])
            s0 = S[i] - mx; s1 = S[i + 1] - mx
            ds = s1 - s0
            if abs(ds) < 1e-9:
                e = np.exp(0.5 * (s0 + s1))
            else:
                e = (np.exp(s1) - np.exp(s0)) / ds
            tot += d * e
        logs[il] = mx + np.log(tot / P)
    mx = logs.max()
    s = 0.0
    for il in range(nlines):
        s += np.exp(logs[il] - mx)
    return mx + np.log(s / nlines)


@njit(parallel=True, cache=True)
def mean_log_zb_prune(nb, l, A, taub, nlines, cap, delta, seed):
    out = np.empty(nb)
    for b in prange(nb):
        np.random.seed(seed + b)
        V = np.random.random(l) * P
        out[b] = log_zb_prune(V, A, taub, nlines, cap, delta)
    return out


@njit(parallel=True, cache=True)
def _block_prune(Vb, A, taub, nlines, cap, delta, lam, w):
    N = Vb.shape[0]
    for j in prange(N):
        w[j] -= lam * log_zb_prune(Vb[j], A, taub, nlines, cap, delta)


@njit(parallel=True, cache=True)
def mean_log_zb(nb, l, A, k, taub, tan, e0, e1, seed):
    out = np.empty(nb)
    for b in prange(nb):
        np.random.seed(seed + b)
        V = np.random.random(l) * P
        out[b] = log_zb(V, A, k, taub, tan, e0, e1)
    return out


@njit(parallel=True, cache=True)
def mean_log_zb_curve(nb, l, A, taub, nlines, cap, seed):
    out = np.empty(nb)
    for b in prange(nb):
        np.random.seed(seed + b)
        V = np.random.random(l) * P
        out[b] = log_zb_curve(V, A, taub, nlines, cap)
    return out


@njit(parallel=True, cache=True)
def _block_curve(Vb, A, taub, nlines, cap, lam, w):
    N = Vb.shape[0]
    for j in prange(N):
        w[j] -= lam * log_zb_curve(Vb[j], A, taub, nlines, cap)



# ---------------------------------------------------------------------------
# Exact 2D block partition function: convex polygons in the entry square
# ---------------------------------------------------------------------------

@njit(cache=True)
def _clip(ex, ey, th, thp, S, nv, c, sgn, ox, oy, oth, othp, oS):
    """Sutherland-Hodgman: keep the part of a convex polygon with sgn*(th - c) >= 0."""
    m = 0
    MAXV = ox.shape[0]
    for i in range(nv):
        j = (i + 1) % nv
        di = sgn * (th[i] - c); dj = sgn * (th[j] - c)
        if di >= 0.0:
            if m >= MAXV:
                return -1
            ox[m] = ex[i]; oy[m] = ey[i]; oth[m] = th[i]; othp[m] = thp[i]; oS[m] = S[i]; m += 1
        if (di >= 0.0) != (dj >= 0.0):
            r = di / (di - dj)
            if m >= MAXV:
                return -1
            ox[m] = ex[i] + r * (ex[j] - ex[i]); oy[m] = ey[i] + r * (ey[j] - ey[i])
            oth[m] = c; othp[m] = thp[i] + r * (thp[j] - thp[i]); oS[m] = S[i] + r * (S[j] - S[i])
            m += 1
    return m


@njit(cache=True)
def _area(ex, ey, nv):
    s = 0.0
    for i in range(nv):
        j = (i + 1) % nv
        s += ex[i] * ey[j] - ex[j] * ey[i]
    return 0.5 * abs(s)


@njit(cache=True)
def _dd2(a, b):
    """(e^b - e^a)/(b - a), stable."""
    d = b - a
    if abs(d) < 1e-6:
        return np.exp(0.5 * (a + b)) * (1.0 + d * d / 24.0)
    return (np.exp(b) - np.exp(a)) / d


@njit(cache=True)
def _tri_exp(s1, s2, s3):
    """Integral of e^{affine} over the unit simplex = divided difference f[s1,s2,s3]."""
    a = s1; b = s2; c = s3
    if a > b:
        a, b = b, a
    if b > c:
        b, c = c, b
    if a > b:
        a, b = b, a
    if c - a < 1e-5:
        m = (a + b + c) / 3.0
        return 0.5 * np.exp(m)
    return (_dd2(b, c) - _dd2(a, b)) / (c - a)


@njit(cache=True)
def _poly_int(ex, ey, S, nv):
    """Integral of e^{S} (S affine) over a convex polygon, by a fan of triangles."""
    tot = 0.0
    for i in range(1, nv - 1):
        ar = 0.5 * abs((ex[i] - ex[0]) * (ey[i + 1] - ey[0]) - (ex[i + 1] - ex[0]) * (ey[i] - ey[0]))
        if ar > 0.0:
            tot += 2.0 * ar * _tri_exp(S[0], S[i], S[i + 1])
    return tot


@njit(cache=True)
def log_zb_poly(V, A, taub, cap, delta, etilt):
    """Exact Z_b(V) = P^-2 int_{[0,P)^2} d theta_0 d theta_-1
         e^{taub sum_j (g(theta_j) - A)} ||M||_F ,
    j = 1..l (and j = -1, 0 if etilt), M the tangent product (Frobenius norm).
    The block map is affine on convex polygons, which are traced exactly;
    a polygon is dropped when its weight is below e^{-delta} of the total
    (deterministic, so Z_b stays a fixed function of V)."""
    l = V.shape[0]
    h = 0.5 * P
    k = 4.0 * A / P
    MAXV = 20
    PX = np.empty((2, cap, MAXV)); PY = np.empty((2, cap, MAXV)); PT = np.empty((2, cap, MAXV))
    PTP = np.empty((2, cap, MAXV)); PS = np.empty((2, cap, MAXV))
    NV = np.zeros((2, cap), np.int64); MM = np.empty((2, cap, 4)); LM = np.empty((2, cap))
    t1x = np.empty(MAXV); t1y = np.empty(MAXV); t1t = np.empty(MAXV); t1p = np.empty(MAXV); t1s = np.empty(MAXV)
    pot = np.empty(cap)
    # four initial squares (slabs of theta_0 and theta_-1)
    cur = 0; n = 0
    for ix in range(2):
        for iy in range(2):
            x0 = ix * h; y0 = iy * h
            xs = (x0, x0 + h, x0 + h, x0); ys = (y0, y0, y0 + h, y0 + h)
            for v in range(4):
                PX[cur, n, v] = xs[v]; PY[cur, n, v] = ys[v]
                PT[cur, n, v] = xs[v]; PTP[cur, n, v] = ys[v]
                # g is affine on each square: evaluate with the slab's own branch
                sx = 0.0; sy = 0.0
                if etilt:
                    tx = (xs[v] - x0) / h; ty = (ys[v] - y0) / h
                    gx = (-A + 2 * A * tx) if ix == 0 else (A - 2 * A * tx)
                    gy = (-A + 2 * A * ty) if iy == 0 else (A - 2 * A * ty)
                    sx = taub * (gx - A) + taub * (gy - A)
                PS[cur, n, v] = sx
            NV[cur, n] = 4
            MM[cur, n, 0] = 1.0; MM[cur, n, 1] = 0.0; MM[cur, n, 2] = 0.0; MM[cur, n, 3] = 1.0
            LM[cur, n] = 0.0
            n += 1
    overflow = False
    for j in range(l):
        for phase in range(2):
            # split every polygon along the kinks of g at the current theta
            nxt = 1 - cur; m2 = 0
            for p in range(n):
                nv = NV[cur, p]
                tmin = 1e300; tmax = -1e300
                for v in range(nv):
                    t = PT[cur, p, v]
                    if t < tmin:
                        tmin = t
                    if t > tmax:
                        tmax = t
                mlo = int(np.floor(tmin / h)); mhi = int(np.floor(tmax / h))
                if mhi > mlo and tmax == mhi * h:
                    mhi -= 1
                for ms in range(mlo, mhi + 1):
                    if m2 >= cap:
                        overflow = True; break
                    if mlo == mhi:
                        nn = nv
                        for v in range(nv):
                            PX[nxt, m2, v] = PX[cur, p, v]; PY[nxt, m2, v] = PY[cur, p, v]
                            PT[nxt, m2, v] = PT[cur, p, v]; PTP[nxt, m2, v] = PTP[cur, p, v]; PS[nxt, m2, v] = PS[cur, p, v]
                    else:
                        nn = _clip(PX[cur, p], PY[cur, p], PT[cur, p], PTP[cur, p], PS[cur, p], nv, ms * h, 1.0,
                                   t1x, t1y, t1t, t1p, t1s)
                        if nn < 3:
                            continue
                        nn = _clip(t1x, t1y, t1t, t1p, t1s, nn, (ms + 1) * h, -1.0,
                                   PX[nxt, m2], PY[nxt, m2], PT[nxt, m2], PTP[nxt, m2], PS[nxt, m2])
                        if nn < 3:
                            continue
                        if _area(PX[nxt, m2], PY[nxt, m2], nn) <= 0.0:
                            continue
                    NV[nxt, m2] = nn
                    for q in range(4):
                        MM[nxt, m2, q] = MM[cur, p, q]
                    LM[nxt, m2] = LM[cur, p]
                    m2 += 1
            cur = nxt; n = m2
            if phase == 0:
                # map theta_j -> theta_{j+1} (affine on each polygon), tangent update
                for p in range(n):
                    nv = NV[cur, p]
                    tc = 0.0
                    for v in range(nv):
                        tc += PT[cur, p, v]
                    tc /= nv
                    _, sg = _g(tc, A)
                    ms = np.floor(tc / h)
                    for v in range(nv):
                        t = PT[cur, p, v]
                        tt = (t - ms * h) / h
                        gv = (-A + 2 * A * tt) if sg > 0 else (A - 2 * A * tt)
                        tn = 2.0 * t - PTP[cur, p, v] + gv - V[j]
                        PTP[cur, p, v] = t; PT[cur, p, v] = tn
                    c = 2.0 + k * sg
                    m11 = MM[cur, p, 0]; m12 = MM[cur, p, 1]; m21 = MM[cur, p, 2]; m22 = MM[cur, p, 3]
                    a11 = c * m11 - m21; a12 = c * m12 - m22
                    MM[cur, p, 0] = a11; MM[cur, p, 1] = a12; MM[cur, p, 2] = m11; MM[cur, p, 3] = m12
                    nr = np.sqrt(a11 * a11 + a12 * a12 + m11 * m11 + m12 * m12)
                    for q in range(4):
                        MM[cur, p, q] /= nr
                    LM[cur, p] += np.log(nr)
            else:
                # add taub * g(theta_{j+1}) (affine on each polygon now)
                for p in range(n):
                    nv = NV[cur, p]
                    tc = 0.0
                    for v in range(nv):
                        tc += PT[cur, p, v]
                    tc /= nv
                    _, sg = _g(tc, A)
                    ms = np.floor(tc / h)
                    for v in range(nv):
                        tt = (PT[cur, p, v] - ms * h) / h
                        gv = (-A + 2 * A * tt) if sg > 0 else (A - 2 * A * tt)
                        PS[cur, p, v] += taub * (gv - A)
        # prune
        if delta > 0.0 and n > 1:
            mx = -1e300
            for p in range(n):
                nv = NV[cur, p]
                smax = -1e300
                for v in range(nv):
                    if PS[cur, p, v] > smax:
                        smax = PS[cur, p, v]
                ar = _area(PX[cur, p], PY[cur, p], nv)
                pot[p] = smax + LM[cur, p] + (np.log(ar) if ar > 0 else -1e300)
                if pot[p] > mx:
                    mx = pot[p]
            tot = 0.0
            for p in range(n):
                tot += np.exp(pot[p] - mx)
            thr = mx + np.log(tot) - delta
            m2 = 0
            for p in range(n):
                if pot[p] >= thr:
                    if m2 != p:
                        nv = NV[cur, p]
                        for v in range(nv):
                            PX[cur, m2, v] = PX[cur, p, v]; PY[cur, m2, v] = PY[cur, p, v]
                            PT[cur, m2, v] = PT[cur, p, v]; PTP[cur, m2, v] = PTP[cur, p, v]; PS[cur, m2, v] = PS[cur, p, v]
                        NV[cur, m2] = nv
                        for q in range(4):
                            MM[cur, m2, q] = MM[cur, p, q]
                        LM[cur, m2] = LM[cur, p]
                    m2 += 1
            n = m2
    # integrate
    lv = np.empty(n)
    mx = -1e300
    for p in range(n):
        nv = NV[cur, p]
        smax = -1e300
        for v in range(nv):
            if PS[cur, p, v] > smax:
                smax = PS[cur, p, v]
        for v in range(nv):
            t1s[v] = PS[cur, p, v] - smax
        I = _poly_int(PX[cur, p], PY[cur, p], t1s, nv)
        fr = np.sqrt(MM[cur, p, 0] ** 2 + MM[cur, p, 1] ** 2 + MM[cur, p, 2] ** 2 + MM[cur, p, 3] ** 2)
        lv[p] = (np.log(I * fr) if I > 0 else -1e300) + smax + LM[cur, p]
        if lv[p] > mx:
            mx = lv[p]
    s = 0.0
    for p in range(n):
        s += np.exp(lv[p] - mx)
    if overflow:
        mx += 0.0  # truncated sum: still a fixed function of V
    return mx + np.log(s / (P * P))


@njit(parallel=True, cache=True)
def mean_log_zb_poly(nb, l, A, taub, cap, delta, etilt, seed):
    out = np.empty(nb)
    for b in prange(nb):
        np.random.seed(seed + b)
        V = np.random.random(l) * P
        out[b] = log_zb_poly(V, A, taub, cap, delta, etilt)
    return out


@njit(parallel=True, cache=True)
def _block_poly(Vb, A, taub, cap, delta, etilt, lam, w):
    N = Vb.shape[0]
    for j in prange(N):
        w[j] -= lam * log_zb_poly(Vb[j], A, taub, cap, delta, etilt)


@njit(parallel=True, cache=True)
def _site(th, thm, gm, p, q, Vb, pos, A, k, tau, u1, u2, w):
    """Advance every walker by one site.  th, thm: theta_{i-1}, theta_{i-2}
    on entry; gm = g(theta_{i-1}).  Draws theta_i from the tilted law."""
    N = th.shape[0]
    ea = np.exp(-tau * A); eb = np.exp(tau * A)
    for j in prange(N):
        if tau > 1e-12:
            gi = np.log(ea + u1[j] * (eb - ea)) / tau
        else:
            gi = -A + 2 * A * u1[j]
        if u2[j] < 0.5:
            ti = (gi + A) / (2 * A) * 0.5 * P; si = 1.0
        else:
            ti = 0.5 * P + (A - gi) / (2 * A) * 0.5 * P; si = -1.0
        Vb[j, pos] = (2.0 * th[j] - thm[j] + gm[j] - ti) % P
        thm[j] = th[j]; th[j] = ti; gm[j] = gi
        a_ = (2.0 + k * si) * p[j] - q[j]
        q[j] = p[j]; p[j] = a_
        nr = np.sqrt(a_ * a_ + q[j] * q[j])
        p[j] /= nr; q[j] /= nr
        w[j] += np.log(nr)


@njit(parallel=True, cache=True)
def _block(Vb, A, k, taub, tan, e0, e1, lam, w):
    N = Vb.shape[0]
    for j in prange(N):
        w[j] -= lam * log_zb(Vb[j], A, k, taub, tan, e0, e1)


def entries(M, seed=12345):
    r = np.random.default_rng(seed)
    return r.uniform(0, P, M), r.uniform(0, P, M)


def growth(A, tau, l, lam, taub, tan=True, M=1024, N=20_000, nblocks=60, burn=10, seed=0, nlines=0, cap=1 << 15, delta=0.0):
    """K^-1 ln <Z prod Z_b^-lam> by population dynamics (per site)."""
    rng = np.random.default_rng(seed)
    k = 4.0 * A / P
    e0, e1 = entries(M)
    th = rng.uniform(0, P, N); thm = rng.uniform(0, P, N)
    gm = np.array([_g(x, A)[0] for x in th])
    ang = rng.uniform(0, np.pi, N); p = np.cos(ang); q = np.sin(ang)
    Vb = np.zeros((N, l))
    x = tau * A
    lmg = np.log(np.sinh(x) / x) if x > 1e-12 else 0.0
    acc = []
    for b in range(nblocks):
        w = np.zeros(N)
        for pos in range(l):
            _site(th, thm, gm, p, q, Vb, pos, A, k, tau, rng.random(N), rng.random(N), w)
        if lam != 0.0:
            if nlines < 0:
                _block_poly(Vb, A, taub, cap, delta, nlines == -2, lam, w)
            elif nlines > 0 and delta > 0:
                _block_prune(Vb, A, taub, nlines, cap, delta, lam, w)
            elif nlines > 0:
                _block_curve(Vb, A, taub, nlines, cap, lam, w)
            else:
                _block(Vb, A, k, taub, tan, e0, e1, lam, w)
        m = w.max(); ww = np.exp(w - m); mean = ww.mean()
        if b >= burn:
            acc.append((np.log(mean) + m) / l + lmg)
        c = np.cumsum(ww); c /= c[-1]
        idx = np.minimum(np.searchsorted(c, (rng.random() + np.arange(N)) / N), N - 1)
        th, thm, gm, p, q = th[idx], thm[idx], gm[idx], p[idx], q[idx]
    acc = np.array(acc)
    nb = 10
    blk = acc[: len(acc) // nb * nb].reshape(nb, -1).mean(1)
    return float(acc.mean()), float(blk.std() / np.sqrt(nb))


def elogzb(A, l, taub, tan=True, M=1024, nb=40_000, seed=777, nlines=0, cap=1 << 15, delta=0.0):
    k = 4.0 * A / P
    e0, e1 = entries(M)
    if nlines < 0:
        z = mean_log_zb_poly(nb, l, A, taub, cap, delta, nlines == -2, seed)
    elif nlines > 0 and delta > 0:
        z = mean_log_zb_prune(nb, l, A, taub, nlines, cap, delta, seed)
    elif nlines > 0:
        z = mean_log_zb_curve(nb, l, A, taub, nlines, cap, seed)
    else:
        z = mean_log_zb(nb, l, A, k, taub, tan, e0, e1, seed)
    assert z.min() > -1e200, "curve capacity exceeded"
    return float(z.mean()), float(z.std() / np.sqrt(nb))


def bound(A, tau, l, lam, taub, tan=True, M=1024, nlines=0, cap=1 << 15, nb=40_000, delta=0.0, **kw):
    G, eG = growth(A, tau, l, lam, taub, tan, M, nlines=nlines, cap=cap, delta=delta, **kw)
    if lam == 0.0:
        return G, eG
    E, eE = elogzb(A, l, taub, tan, M, nb=nb, nlines=nlines, cap=cap, delta=delta)
    return G + lam / l * E, float(np.hypot(eG, lam / l * eE))


def lam_exact(A, tau):
    from annealed import gle
    x = tau * A
    return gle(1.0, A) + (np.log(np.sinh(x) / x) if x > 1e-8 else 0.0)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    if "--check" in sys.argv:
        for A, tA in ((0.6, 0.0), (0.6, 2.0), (0.3, 1.5)):
            tau = tA / A
            g, e = growth(A, tau, 4, 0.0, 0.0)
            print(f"A={A} tauA={tA}: pop {g:.4f}+-{e:.4f}  exact {lam_exact(A, tau):.4f}", flush=True)
        raise SystemExit
    # production: exact polygon Z_b (mode -1), mild pruning, two population sizes
    res = json.load(open("blockbound.json")) if os.path.exists("blockbound.json") else {}
    def opt(name, default):
        if name in sys.argv:
            j = sys.argv.index(name)
            return [float(x) for x in sys.argv[j + 1].split(",")]
        return default
    ells = [int(x) for x in opt("--l", [4, 8, 16])]
    As = opt("--A", [0.3, 0.6, 1.0])
    tAs = opt("--tauA", [1.0, 2.0, 3.0, 4.0])
    Ns = [int(x) for x in opt("--N", [2000, 8000])]
    MODE, DELTA, CAP = -1, 30.0, 1 << 15
    for l in ells:
        for A in As:
            rows = res.setdefault(str(A), [])
            done = {(r["l"], round(r["tauA"], 3)) for r in rows}
            for tA in tAs:
                if (l, round(tA, 3)) in done:
                    continue
                tau = tA / A
                t0 = time.time()
                E, eE = elogzb(A, l, tau, nb=8000, nlines=MODE, cap=CAP, delta=DELTA)
                L = lam_exact(A, tau)
                row = dict(l=l, tauA=tA, tau=tau, lam=1.0, Lambda=L, ElnZb=E / l + tA, ElnZb_err=eE / l, B=[], B_err=[], N=Ns)
                for N in Ns:
                    G, eG = growth(A, tau, l, 1.0, tau, N=N, nblocks=max(24, 384 // l), burn=6, nlines=MODE, cap=CAP, delta=DELTA)
                    row["B"].append(G + E / l); row["B_err"].append(float(np.hypot(eG, eE / l)))
                rows.append(row)
                print(f"A={A} l={l} tauA={tA}: B(N)={[round(b, 4) for b in row['B']]} +- {[round(e, 4) for e in row['B_err']]} "
                      f"Lambda={L:.4f} B/tau={row['B'][-1]/tau:.4f} ({time.time()-t0:.0f}s)", flush=True)
                json.dump(res, open("blockbound.json", "w"), indent=1)
