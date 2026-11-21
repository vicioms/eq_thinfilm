"""Exact enumeration of the equilibria of a driven elastic chain.

Model (kappa = 1, zero mass, free ends):

    H(u) = 1/2 sum_{i=0}^{K-2} (u_{i+1} - u_i)^2 + sum_{i=0}^{K-1} [V_i(u_i) - f u_i]

The derivative g_i(u) = V_i'(u) is a triangle wave of amplitude A and
period P with an independent uniform phase phi_i at every site, so that at
any fixed u the pair (g_i, V_i'') is uniform on [-A, A] times a fair sign
+-k, k = 4A/P, and the two are independent.  This is the hypothesis under
which the Kac-Rice average factorises (FLRT 2018, Eq. 19-21).

Equilibria are found by shooting from the free end: u_0 = a, w_{-1} = 0,

    w_i = w_{i-1} + g_i(u_i) - f,      u_{i+1} = u_i + w_i,

and the free end at site K-1 requires F(a) = w_{K-1}(a) = 0.  Since the
disorder is P-periodic in u, equilibria modulo the global shift by P are in
one-to-one correspondence with the zeros of F on [0, P).

Because g_i is piecewise linear, u_i(a) and w_{i-1}(a) are piecewise affine
in a.  We split the shooting interval at every kink of every g_i (depth-first)
so that on each leaf the whole trajectory is affine; F then has at most one
zero per leaf.  The count is exact up to floating point, and the local
parametrisation (each leaf stores its own origin and length) keeps the
round-off from being amplified by the chaotic forward map.

On every leaf the tangent (psi_i, phi_i) = (du_i/da, dw_i/da) is constant and
gives the leading principal minors of the Hessian -Delta_N + U:
psi_1, ..., psi_{K-1}, phi_{K-1}.  The number of sign changes in
(1, psi_1, ..., psi_{K-1}, phi_{K-1}) is the Morse index (Sturm), and
|phi_{K-1}| = |det Hessian| (Gelfand-Yaglom).
"""
import numpy as np
from numba import njit

STACK = 1 << 16


@njit(cache=True)
def _tri(u, ph, A, P):
    """Triangle wave value and slope sign at u (phase ph)."""
    half = 0.5 * P
    x = (u - ph) / half
    n = np.floor(x)
    t = x - n
    if (int(n) & 1) == 0:
        return -A + 2.0 * A * t, 1.0
    return A - 2.0 * A * t, -1.0


@njit(cache=True)
def count(phase, A, P, f, kmax_index, stop_at_first, nlogbins, logmax):
    """Enumerate all equilibria for one sample.

    Returns (n_tot, n_nodes, index_hist[kmax_index+1], logdet_hist[nlogbins],
    sum_logdet, overflow_flag).  index_hist counts equilibria by Morse index
    (the last bin collects index >= kmax_index); logdet_hist bins
    ln|det H| / K on [-logmax, logmax].
    """
    K = phase.shape[0]
    k = 4.0 * A / P
    half = 0.5 * P
    # node: depth, len, u0, w0 (w_{i-1}), psi, phi(prev), nneg, lastsign
    st_d = np.empty(STACK, np.int64)
    st_L = np.empty(STACK)
    st_u = np.empty(STACK)
    st_w = np.empty(STACK)
    st_ps = np.empty(STACK)
    st_ph = np.empty(STACK)
    st_n = np.empty(STACK, np.int64)
    st_s = np.empty(STACK)
    top = 0
    st_d[0] = 0
    st_L[0] = P
    st_u[0] = 0.0
    st_w[0] = 0.0
    st_ps[0] = 1.0
    st_ph[0] = 0.0
    st_n[0] = 0
    st_s[0] = 1.0
    top = 1
    ntot = 0
    nodes = 0
    ihist = np.zeros(kmax_index + 1, np.int64)
    lhist = np.zeros(nlogbins, np.int64)
    sumlog = 0.0
    overflow = False
    while top > 0:
        top -= 1
        i = st_d[top]
        L = st_L[top]
        u0 = st_u[top]
        w0 = st_w[top]
        ps = st_ps[top]
        phm = st_ph[top]
        nn = st_n[top]
        ls = st_s[top]
        nodes += 1
        ph_i = phase[i]
        # kinks of g_i crossed by u(t) = u0 + ps t, t in [0, L]
        ua = u0
        ub = u0 + ps * L
        xa = (ua - ph_i) / half
        xb = (ub - ph_i) / half
        if ps > 0:
            m_lo = np.floor(xa) + 1.0
            m_hi = np.ceil(xb) - 1.0
            step = 1.0
        else:
            m_lo = np.ceil(xa) - 1.0
            m_hi = np.floor(xb) + 1.0
            step = -1.0
        nk = int(round((m_hi - m_lo) * step)) + 1
        if nk < 0:
            nk = 0
        # pieces: [t_j, t_{j+1}] with t_0 = 0, t_{nk+1} = L
        t_prev = 0.0
        u_prev = u0
        for j in range(nk + 1):
            if j < nk:
                mk = m_lo + step * j
                uk = ph_i + half * mk
                t_next = (uk - u0) / ps
                if t_next > L:
                    t_next = L
            else:
                t_next = L
                uk = ub
            Lp = t_next - t_prev
            if Lp > 0.0:
                # slope of g on this piece: evaluate at midpoint
                umid = 0.5 * (u_prev + uk)
                gmid, sg = _tri(umid, ph_i, A, P)
                U = sg * k  # V_i'' on this piece
                g1 = gmid + U * (u_prev - umid)  # V_i'(u) at the piece start
                w1 = w0 + phm * t_prev + g1 - f
                phi = phm + U * ps  # phi_i = phi_{i-1} + U_i psi_i
                if i == K - 1:
                    # free end: zero of w_{K-1}(t) = w1 + phi (t - t_prev)
                    if phi != 0.0:
                        tz = -w1 / phi
                        if tz >= 0.0 and tz < Lp:
                            # last minor: det(-Delta_N + U) = phi_{K-1}
                            sdet = 1.0 if phi > 0 else -1.0
                            nidx = nn + (1 if sdet != ls else 0)
                            ntot += 1
                            if nidx >= kmax_index:
                                ihist[kmax_index] += 1
                            else:
                                ihist[nidx] += 1
                            ld = np.log(abs(phi)) / K
                            sumlog += ld
                            b = int((ld + logmax) / (2 * logmax) * nlogbins)
                            if b < 0:
                                b = 0
                            if b >= nlogbins:
                                b = nlogbins - 1
                            lhist[b] += 1
                            if stop_at_first:
                                return ntot, nodes, ihist, lhist, sumlog, overflow
                else:
                    psn = ps + phi
                    un = u_prev + w1
                    sn = 1.0 if psn > 0 else -1.0
                    nnn = nn + (1 if sn != ls else 0)
                    if top >= STACK:
                        overflow = True
                        return ntot, nodes, ihist, lhist, sumlog, overflow
                    st_d[top] = i + 1
                    st_L[top] = Lp
                    st_u[top] = un
                    st_w[top] = w1
                    st_ps[top] = psn
                    st_ph[top] = phi
                    st_n[top] = nnn
                    st_s[top] = sn
                    top += 1
            t_prev = t_next
            u_prev = uk
    return ntot, nodes, ihist, lhist, sumlog, overflow


def draw_phases(rng, K, P):
    return rng.uniform(0.0, P, size=K)


def run_sample(phase, A, P, f, kmax_index=8, stop_at_first=False):
    return count(phase, A, P, f, kmax_index, stop_at_first, 200, 3.0)


@njit(cache=True)
def count_levels(phase, A, P, levels, nub):
    """One f = 0 enumeration, all forces at once (statistical tilt symmetry).

    With w~_i = w_i + f (i+1) and u~_i = u_i + f i (i+1)/2 the shooting map at
    force f becomes the f = 0 map with phases shifted by the deterministic
    amount f i (i+1)/2, which leaves their (uniform, independent) law
    unchanged.  The free-end condition w_{K-1} = 0 becomes W_K(a) = f K, with
    W_K(a) = w_{K-1}(a) at f = 0.  Hence, in law and for each f separately,

        N_K(f) = #{a in [0, P): W_K(a) = f K},   f_c(K) = max_a W_K(a) / K.

    The Hessian (and so Morse index and det) is the same in both frames.

    levels: sorted array of W values (= f K).  Returns
      N[l], Nmin[l]   crossings of level l (all, Morse index 0)
      slog[l]         sum of ln|det H| / K over the crossings
      ub[l, nub]      histogram of u_{K-1} mod P at the crossings
      Wmax, Wmin      extreme values of W_K over the circle
      nodes, overflow
    """
    K = phase.shape[0]
    k = 4.0 * A / P
    half = 0.5 * P
    nl = levels.shape[0]
    N = np.zeros(nl, np.int64)
    Nmin = np.zeros(nl, np.int64)
    slog = np.zeros(nl)
    ub = np.zeros((nl, nub), np.int64)
    Wmax = -1e300
    Wmin = 1e300
    st_d = np.empty(STACK, np.int64)
    st_L = np.empty(STACK)
    st_u = np.empty(STACK)
    st_w = np.empty(STACK)
    st_ps = np.empty(STACK)
    st_ph = np.empty(STACK)
    st_n = np.empty(STACK, np.int64)
    st_s = np.empty(STACK)
    st_d[0] = 0
    st_L[0] = P
    st_u[0] = 0.0
    st_w[0] = 0.0
    st_ps[0] = 1.0
    st_ph[0] = 0.0
    st_n[0] = 0
    st_s[0] = 1.0
    top = 1
    nodes = 0
    overflow = False
    while top > 0:
        top -= 1
        i = st_d[top]
        L = st_L[top]
        u0 = st_u[top]
        w0 = st_w[top]
        ps = st_ps[top]
        phm = st_ph[top]
        nn = st_n[top]
        ls = st_s[top]
        nodes += 1
        ph_i = phase[i]
        ua = u0
        ubd = u0 + ps * L
        xa = (ua - ph_i) / half
        xb = (ubd - ph_i) / half
        if ps > 0:
            m_lo = np.floor(xa) + 1.0
            m_hi = np.ceil(xb) - 1.0
            step = 1.0
        else:
            m_lo = np.ceil(xa) - 1.0
            m_hi = np.floor(xb) + 1.0
            step = -1.0
        nk = int(round((m_hi - m_lo) * step)) + 1
        if nk < 0:
            nk = 0
        t_prev = 0.0
        u_prev = u0
        for j in range(nk + 1):
            if j < nk:
                mk = m_lo + step * j
                uk = ph_i + half * mk
                t_next = (uk - u0) / ps
                if t_next > L:
                    t_next = L
            else:
                t_next = L
                uk = ubd
            Lp = t_next - t_prev
            if Lp > 0.0:
                umid = 0.5 * (u_prev + uk)
                gmid, sg = _tri(umid, ph_i, A, P)
                U = sg * k
                g1 = gmid + U * (u_prev - umid)
                w1 = w0 + phm * t_prev + g1
                phi = phm + U * ps
                if i == K - 1:
                    w2 = w1 + phi * Lp
                    lo = min(w1, w2)
                    hi = max(w1, w2)
                    if hi > Wmax:
                        Wmax = hi
                    if lo < Wmin:
                        Wmin = lo
                    if phi != 0.0:
                        l0 = np.searchsorted(levels, lo)
                        l1 = np.searchsorted(levels, hi)
                        if l1 > l0:
                            sdet = 1.0 if phi > 0 else -1.0
                            nidx = nn + (1 if sdet != ls else 0)
                            ld = np.log(abs(phi)) / K
                            for l in range(l0, l1):
                                tz = (levels[l] - w1) / phi
                                uz = u_prev + ps * tz
                                b = int(np.floor((uz - P * np.floor(uz / P)) / P * nub))
                                if b >= nub:
                                    b = nub - 1
                                N[l] += 1
                                if nidx == 0:
                                    Nmin[l] += 1
                                slog[l] += ld
                                ub[l, b] += 1
                else:
                    psn = ps + phi
                    un = u_prev + w1
                    sn = 1.0 if psn > 0 else -1.0
                    nnn = nn + (1 if sn != ls else 0)
                    if top >= STACK:
                        overflow = True
                        return N, Nmin, slog, ub, Wmax, Wmin, nodes, overflow
                    st_d[top] = i + 1
                    st_L[top] = Lp
                    st_u[top] = un
                    st_w[top] = w1
                    st_ps[top] = psn
                    st_ph[top] = phi
                    st_n[top] = nnn
                    st_s[top] = sn
                    top += 1
            t_prev = t_next
            u_prev = uk
    return N, Nmin, slog, ub, Wmax, Wmin, nodes, overflow


@njit(cache=True)
def count_levels_pl(G, off, P, levels, nub):
    """count_levels for a general periodic piecewise-linear force.

    Site i has V_i'(u) linear between the nodes u = off[i] + m h, h = P / M,
    with values G[i, m mod M].  The triangle wave is M = 2.  Used here for the
    piecewise-linear interpolation of FLRT's Gaussian xi1 cos u + xi2 sin u
    disorder (random off[i] keeps the law shift-invariant).
    """
    K, M = G.shape
    h = P / M
    nl = levels.shape[0]
    N = np.zeros(nl, np.int64)
    Nmin = np.zeros(nl, np.int64)
    slog = np.zeros(nl)
    ub = np.zeros((nl, nub), np.int64)
    Wmax = -1e300
    Wmin = 1e300
    st_d = np.empty(STACK, np.int64)
    st_L = np.empty(STACK)
    st_u = np.empty(STACK)
    st_w = np.empty(STACK)
    st_ps = np.empty(STACK)
    st_ph = np.empty(STACK)
    st_n = np.empty(STACK, np.int64)
    st_s = np.empty(STACK)
    st_d[0] = 0
    st_L[0] = P
    st_u[0] = 0.0
    st_w[0] = 0.0
    st_ps[0] = 1.0
    st_ph[0] = 0.0
    st_n[0] = 0
    st_s[0] = 1.0
    top = 1
    nodes = 0
    overflow = False
    while top > 0:
        top -= 1
        i = st_d[top]
        L = st_L[top]
        u0 = st_u[top]
        w0 = st_w[top]
        ps = st_ps[top]
        phm = st_ph[top]
        nn = st_n[top]
        ls = st_s[top]
        nodes += 1
        o = off[i]
        ubd = u0 + ps * L
        xa = (u0 - o) / h
        xb = (ubd - o) / h
        if ps > 0:
            m_lo = np.floor(xa) + 1.0
            m_hi = np.ceil(xb) - 1.0
            step = 1.0
        else:
            m_lo = np.ceil(xa) - 1.0
            m_hi = np.floor(xb) + 1.0
            step = -1.0
        nk = int(round((m_hi - m_lo) * step)) + 1
        if nk < 0:
            nk = 0
        t_prev = 0.0
        u_prev = u0
        for j in range(nk + 1):
            if j < nk:
                mk = m_lo + step * j
                uk = o + h * mk
                t_next = (uk - u0) / ps
                if t_next > L:
                    t_next = L
            else:
                t_next = L
                uk = ubd
            Lp = t_next - t_prev
            if Lp > 0.0:
                umid = 0.5 * (u_prev + uk)
                x = (umid - o) / h
                mc = np.floor(x)
                m0 = int(mc) % M
                if m0 < 0:
                    m0 += M
                m1 = (m0 + 1) % M
                U = (G[i, m1] - G[i, m0]) / h
                g1 = G[i, m0] + U * (u_prev - (o + h * mc))
                w1 = w0 + phm * t_prev + g1
                phi = phm + U * ps
                if i == K - 1:
                    w2 = w1 + phi * Lp
                    lo = min(w1, w2)
                    hi = max(w1, w2)
                    if hi > Wmax:
                        Wmax = hi
                    if lo < Wmin:
                        Wmin = lo
                    if phi != 0.0:
                        l0 = np.searchsorted(levels, lo)
                        l1 = np.searchsorted(levels, hi)
                        if l1 > l0:
                            sdet = 1.0 if phi > 0 else -1.0
                            nidx = nn + (1 if sdet != ls else 0)
                            ld = np.log(abs(phi)) / K
                            for l in range(l0, l1):
                                tz = (levels[l] - w1) / phi
                                uz = u_prev + ps * tz
                                b = int(np.floor((uz - P * np.floor(uz / P)) / P * nub))
                                if b >= nub:
                                    b = nub - 1
                                N[l] += 1
                                if nidx == 0:
                                    Nmin[l] += 1
                                slog[l] += ld
                                ub[l, b] += 1
                else:
                    psn = ps + phi
                    un = u_prev + w1
                    sn = 1.0 if psn > 0 else -1.0
                    nnn = nn + (1 if sn != ls else 0)
                    if top >= STACK:
                        overflow = True
                        return N, Nmin, slog, ub, Wmax, Wmin, nodes, overflow
                    st_d[top] = i + 1
                    st_L[top] = Lp
                    st_u[top] = un
                    st_w[top] = w1
                    st_ps[top] = psn
                    st_ph[top] = phi
                    st_n[top] = nnn
                    st_s[top] = sn
                    top += 1
            t_prev = t_next
            u_prev = uk
    return N, Nmin, slog, ub, Wmax, Wmin, nodes, overflow


def cos_nodes(rng, K, M, sigma, P=2 * np.pi):
    """Nodes of the PL interpolation of V_i'(u) = xi1 cos u + xi2 sin u
    (xi ~ N(0, sigma^2)), on the grid off_i + m P / M with a uniform offset."""
    h = P / M
    xi = rng.normal(0.0, sigma, size=(K, 2))
    off = rng.uniform(0.0, h, size=K)
    um = off[:, None] + h * np.arange(M)[None, :]
    G = xi[:, :1] * np.cos(um) + xi[:, 1:] * np.sin(um)
    return G, off


@njit(cache=True)
def record_level(phase, A, P, level, maxrec):
    """f = 0 enumeration (triangle model) recording every crossing of the
    single level W_K = level: returns (u_{K-1} at the crossing, Morse index,
    ln|det H|, n, ln|psi_{K-1}|), psi = du_{K-1}/da the tangent along the
    shooting curve.  Configurations are rebuilt backward with configs()."""
    K = phase.shape[0]
    k = 4.0 * A / P
    half = 0.5 * P
    uz_out = np.empty(maxrec)
    idx_out = np.empty(maxrec, np.int64)
    ld_out = np.empty(maxrec)
    lp_out = np.empty(maxrec)
    nrec = 0
    st_d = np.empty(STACK, np.int64)
    st_L = np.empty(STACK)
    st_u = np.empty(STACK)
    st_w = np.empty(STACK)
    st_ps = np.empty(STACK)
    st_ph = np.empty(STACK)
    st_n = np.empty(STACK, np.int64)
    st_s = np.empty(STACK)
    st_d[0] = 0; st_L[0] = P; st_u[0] = 0.0; st_w[0] = 0.0
    st_ps[0] = 1.0; st_ph[0] = 0.0; st_n[0] = 0; st_s[0] = 1.0
    top = 1
    while top > 0:
        top -= 1
        i = st_d[top]; L = st_L[top]; u0 = st_u[top]; w0 = st_w[top]
        ps = st_ps[top]; phm = st_ph[top]; nn = st_n[top]; ls = st_s[top]
        ph_i = phase[i]
        ubd = u0 + ps * L
        xa = (u0 - ph_i) / half
        xb = (ubd - ph_i) / half
        if ps > 0:
            m_lo = np.floor(xa) + 1.0; m_hi = np.ceil(xb) - 1.0; step = 1.0
        else:
            m_lo = np.ceil(xa) - 1.0; m_hi = np.floor(xb) + 1.0; step = -1.0
        nk = int(round((m_hi - m_lo) * step)) + 1
        if nk < 0:
            nk = 0
        t_prev = 0.0
        u_prev = u0
        for j in range(nk + 1):
            if j < nk:
                uk = ph_i + half * (m_lo + step * j)
                t_next = (uk - u0) / ps
                if t_next > L:
                    t_next = L
            else:
                t_next = L
                uk = ubd
            Lp = t_next - t_prev
            if Lp > 0.0:
                umid = 0.5 * (u_prev + uk)
                gmid, sg = _tri(umid, ph_i, A, P)
                U = sg * k
                g1 = gmid + U * (u_prev - umid)
                w1 = w0 + phm * t_prev + g1
                phi = phm + U * ps
                if i == K - 1:
                    if phi != 0.0:
                        tz = (level - w1) / phi
                        if tz >= 0.0 and tz < Lp and nrec < maxrec:
                            sdet = 1.0 if phi > 0 else -1.0
                            uz_out[nrec] = u_prev + ps * tz
                            idx_out[nrec] = nn + (1 if sdet != ls else 0)
                            ld_out[nrec] = np.log(abs(phi))
                            lp_out[nrec] = np.log(abs(ps))
                            nrec += 1
                else:
                    psn = ps + phi
                    sn = 1.0 if psn > 0 else -1.0
                    st_d[top] = i + 1; st_L[top] = Lp; st_u[top] = u_prev + w1
                    st_w[top] = w1; st_ps[top] = psn; st_ph[top] = phi
                    st_n[top] = nn + (1 if sn != ls else 0); st_s[top] = sn
                    top += 1
            t_prev = t_next
            u_prev = uk
    return uz_out[:nrec], idx_out[:nrec], ld_out[:nrec], nrec, lp_out[:nrec]



@njit(cache=True)
def record_level_a(phase, A, P, level, maxrec):
    """As record_level, but also returns the absolute shooting parameter a of
    every crossing (leaf origins are tracked; local t is measured in units
    of a since psi_0 = 1).  Original docstring:
    f = 0 enumeration (triangle model) recording every crossing of the
    single level W_K = level: returns (u_{K-1} at the crossing, Morse index,
    ln|det H|, n, ln|psi_{K-1}|), psi = du_{K-1}/da the tangent along the
    shooting curve.  Configurations are rebuilt backward with configs()."""
    K = phase.shape[0]
    k = 4.0 * A / P
    half = 0.5 * P
    uz_out = np.empty(maxrec)
    idx_out = np.empty(maxrec, np.int64)
    ld_out = np.empty(maxrec)
    lp_out = np.empty(maxrec)
    a_out = np.empty(maxrec)
    nrec = 0
    st_a = np.empty(STACK)
    st_a[0] = 0.0
    st_d = np.empty(STACK, np.int64)
    st_L = np.empty(STACK)
    st_u = np.empty(STACK)
    st_w = np.empty(STACK)
    st_ps = np.empty(STACK)
    st_ph = np.empty(STACK)
    st_n = np.empty(STACK, np.int64)
    st_s = np.empty(STACK)
    st_d[0] = 0; st_L[0] = P; st_u[0] = 0.0; st_w[0] = 0.0
    st_ps[0] = 1.0; st_ph[0] = 0.0; st_n[0] = 0; st_s[0] = 1.0
    top = 1
    while top > 0:
        top -= 1
        i = st_d[top]; L = st_L[top]; u0 = st_u[top]; w0 = st_w[top]
        ps = st_ps[top]; phm = st_ph[top]; nn = st_n[top]; ls = st_s[top]
        a0 = st_a[top]
        ph_i = phase[i]
        ubd = u0 + ps * L
        xa = (u0 - ph_i) / half
        xb = (ubd - ph_i) / half
        if ps > 0:
            m_lo = np.floor(xa) + 1.0; m_hi = np.ceil(xb) - 1.0; step = 1.0
        else:
            m_lo = np.ceil(xa) - 1.0; m_hi = np.floor(xb) + 1.0; step = -1.0
        nk = int(round((m_hi - m_lo) * step)) + 1
        if nk < 0:
            nk = 0
        t_prev = 0.0
        u_prev = u0
        for j in range(nk + 1):
            if j < nk:
                uk = ph_i + half * (m_lo + step * j)
                t_next = (uk - u0) / ps
                if t_next > L:
                    t_next = L
            else:
                t_next = L
                uk = ubd
            Lp = t_next - t_prev
            if Lp > 0.0:
                umid = 0.5 * (u_prev + uk)
                gmid, sg = _tri(umid, ph_i, A, P)
                U = sg * k
                g1 = gmid + U * (u_prev - umid)
                w1 = w0 + phm * t_prev + g1
                phi = phm + U * ps
                if i == K - 1:
                    if phi != 0.0:
                        tz = (level - w1) / phi
                        if tz >= 0.0 and tz < Lp and nrec < maxrec:
                            sdet = 1.0 if phi > 0 else -1.0
                            uz_out[nrec] = u_prev + ps * tz
                            idx_out[nrec] = nn + (1 if sdet != ls else 0)
                            ld_out[nrec] = np.log(abs(phi))
                            lp_out[nrec] = np.log(abs(ps))
                            a_out[nrec] = a0 + t_prev + tz
                            nrec += 1
                else:
                    psn = ps + phi
                    sn = 1.0 if psn > 0 else -1.0
                    st_d[top] = i + 1; st_L[top] = Lp; st_u[top] = u_prev + w1
                    st_w[top] = w1; st_ps[top] = psn; st_ph[top] = phi
                    st_n[top] = nn + (1 if sn != ls else 0); st_s[top] = sn
                    st_a[top] = a0 + t_prev
                    top += 1
            t_prev = t_next
            u_prev = uk
    return uz_out[:nrec], idx_out[:nrec], ld_out[:nrec], nrec, lp_out[:nrec], a_out[:nrec]



def configs(phase, A, P, level, uz):
    """Rebuild equilibrium configurations (zero-force frame) backward from the
    final state (u_{K-1}, w_{K-1} = level).  Returns u[n, K] and the residual
    |w_{-1}| (should be ~0)."""
    K = len(phase)
    half = 0.5 * P
    n = len(uz)
    u = np.empty((n, K))
    u[:, K - 1] = uz
    w = np.full(n, float(level))          # w_{K-1}
    for i in range(K - 1, -1, -1):
        x = (u[:, i] - phase[i]) / half
        m = np.floor(x)
        t = x - m
        g = np.where(m.astype(np.int64) % 2 == 0, -A + 2 * A * t, A - 2 * A * t)
        w = w - g                          # w_{i-1}
        if i > 0:
            u[:, i - 1] = u[:, i] - w
    return u, np.abs(w)
