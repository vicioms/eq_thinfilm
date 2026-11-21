"""Exact depinning threshold for k = |V''| > 2 (triangle model).

A site on the falling side of its force peak has V'' = -k, so the Hessian
-Delta_N + U has a diagonal entry 2 - k (1 - k at a free end) < 0 there, and
the configuration is unstable as soon as k > 2.  So for k > 2 the minima
are exactly the configurations with every site on a rising branch, where
V_i'(u) = A + k (u - q_i), q_i = p_i + n_i P the peak (p_i = phi_i + P/2),
u_i - q_i in [-P/2, 0].  Force balance Delta_N^2 u = V' - f is then linear:

    x = u - q = y - (A - f)/k,   y = (k - Delta_N^2)^{-1} Delta_N^2 q,

using (k - Delta_N^2)^{-1} 1 = 1/k.  A lift n is valid iff
max y - min y <= P/2, and then for all A - f >= k max y.  The depinning
threshold is the largest force with a stable equilibrium, so

    f_c = A - min_n  k max_i y_i(n)       (exact for k > 2).

Lifts are enumerated through the slopes w_j = q_{j+1} - q_j in d_j + P Z,
with curvatures |Delta^2 q| <= L (L = 2P is generous; checked against
exact thresholds).

    python3 strong_exact.py   -> strong_exact.json
"""
import json
import numpy as np
from numba import njit
from multiprocessing import Pool
from qcount import count

P = 4.0


@njit(cache=True)
def _solve(kk, c):
    """(k - Delta_N^2)^{-1} c, tridiagonal: diag k+2 (k+1 at ends), off -1."""
    K = c.shape[0]
    a = np.empty(K); b = np.empty(K)
    diag = np.full(K, kk + 2.0)
    diag[0] = kk + 1.0
    diag[K - 1] = kk + 1.0
    # Thomas algorithm, off-diagonals -1
    a[0] = diag[0]
    b[0] = c[0]
    for i in range(1, K):
        m = -1.0 / a[i - 1]
        a[i] = diag[i] + m          # diag - (-1)(-1)/a = diag - 1/a
        b[i] = c[i] - m * b[i - 1]  # c + b/a
    y = np.empty(K)
    y[K - 1] = b[K - 1] / a[K - 1]
    for i in range(K - 2, -1, -1):
        y[i] = (b[i] + y[i + 1]) / a[i]
    return y


@njit(cache=True)
def min_beta(d, kk, L):
    """min over lifts of k max y, among valid lifts (max y - min y <= P/2)."""
    K = d.shape[0] + 1
    best = 1e300
    # iterative DFS over slopes w_0 .. w_{K-2}
    nmax = int(2 * L / P) + 3
    cand = np.empty((K - 1, nmax))
    ncand = np.zeros(K - 1, np.int64)
    pos = np.zeros(K - 1, np.int64)
    w = np.zeros(K - 1)
    c = np.zeros(K)

    j = 0
    # candidates for w_0: curvature c_0 = w_0 in [-L, L]
    lo = -L
    m0 = np.ceil((lo - d[0]) / P)
    n = 0
    while d[0] + (m0 + n) * P <= L and n < nmax:
        cand[0, n] = d[0] + (m0 + n) * P
        n += 1
    ncand[0] = n
    pos[0] = 0
    while j >= 0:
        if pos[j] >= ncand[j]:
            j -= 1
            if j >= 0:
                pos[j] += 1
            continue
        w[j] = cand[j, pos[j]]
        if j == K - 2:
            # complete: curvatures
            c[0] = w[0]
            for i in range(1, K - 1):
                c[i] = w[i] - w[i - 1]
            c[K - 1] = -w[K - 2]
            if abs(c[K - 1]) <= L:
                y = _solve(kk, c)
                ymax = y.max(); ymin = y.min()
                if ymax - ymin <= P / 2 and kk * ymax < best:
                    best = kk * ymax
            pos[j] += 1
            continue
        # candidates for w_{j+1} with curvature w_{j+1} - w_j in [-L, L]
        jn = j + 1
        m0 = np.ceil((w[j] - L - d[jn]) / P)
        n = 0
        while d[jn] + (m0 + n) * P <= w[j] + L and n < nmax:
            cand[jn, n] = d[jn] + (m0 + n) * P
            n += 1
        ncand[jn] = n
        pos[jn] = 0
        j = jn
    return best


def fc_predicted(phase, A, L=2 * P):
    p = phase + P / 2
    d = np.diff(p)
    return A - min_beta(d, A * 4 / P, L)


def fc_exact(phase, A):
    """Sample threshold by bisection on 'an equilibrium exists' (count)."""
    lo, hi = 0.0, A
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if count(phase, A, P, mid, 8, True, 10, 3.0)[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def work(args):
    seed, K, A = args
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0, P, K)
    return fc_exact(ph, A), fc_predicted(ph, A)



@njit(cache=True)
def min_beta_bb(d, kk, L, G, best0):
    """Branch and bound version of min_beta.  G = (k - Delta_N^2)^{-1}.
    Lower bound for a partial lift (curvatures c_0..c_j fixed): for i <= j,
    k y_i >= k (sum_{m<=j} G_im c_m - L sum_{m>j} G_im)."""
    K = d.shape[0] + 1
    best = best0
    nmax = int(2 * L / P) + 3
    cand = np.empty((K - 1, nmax))
    ncand = np.zeros(K - 1, np.int64)
    pos = np.zeros(K - 1, np.int64)
    w = np.zeros(K - 1)
    c = np.zeros(K)
    part = np.zeros((K, K))   # part[j, i] = sum_{m<=j} G_im c_m
    tail = np.zeros((K, K))   # tail[j, i] = sum_{m>j} G_im
    for i in range(K):
        s = 0.0
        for jj in range(K - 1, -1, -1):
            tail[jj, i] = s
            s += G[i, jj]
    # candidates for w_0 (descending: low curvature lifts tend to be large w)
    m0 = np.floor((L - d[0]) / P)
    n = 0
    while d[0] + (m0 - n) * P >= -L and n < nmax:
        cand[0, n] = d[0] + (m0 - n) * P
        n += 1
    ncand[0] = n
    pos[0] = 0
    j = 0
    while j >= 0:
        if pos[j] >= ncand[j]:
            j -= 1
            if j >= 0:
                pos[j] += 1
            continue
        w[j] = cand[j, pos[j]]
        c[j] = w[j] if j == 0 else w[j] - w[j - 1]
        # partial sums up to curvature index j
        lb = -1e300
        for i in range(K):
            v = (part[j - 1, i] if j > 0 else 0.0) + G[i, j] * c[j]
            part[j, i] = v
            if i <= j:
                t = kk * (v - L * tail[j, i])
                if t > lb:
                    lb = t
        if lb >= best:
            pos[j] += 1
            continue
        if j == K - 2:
            c[K - 1] = -w[K - 2]
            if abs(c[K - 1]) <= L:
                ymax = -1e300; ymin = 1e300
                for i in range(K):
                    yi = part[j, i] + G[i, K - 1] * c[K - 1]
                    if yi > ymax:
                        ymax = yi
                    if yi < ymin:
                        ymin = yi
                if ymax - ymin <= P / 2 and kk * ymax < best:
                    best = kk * ymax
            pos[j] += 1
            continue
        jn = j + 1
        m0 = np.floor((w[j] + L - d[jn]) / P)
        n = 0
        while d[jn] + (m0 - n) * P >= w[j] - L and n < nmax:
            cand[jn, n] = d[jn] + (m0 - n) * P
            n += 1
        ncand[jn] = n
        pos[jn] = 0
        j = jn
    return best


def green(K, kk):
    D = np.diag(np.full(K, -2.0)) + np.diag(np.ones(K - 1), 1) + np.diag(np.ones(K - 1), -1)
    D[0, 0] = -1.0
    D[-1, -1] = -1.0
    return np.linalg.inv(kk * np.eye(K) - D)


def fc_predicted_bb(phase, A, L=2 * P):
    p = phase + P / 2
    d = np.diff(p)
    K = len(phase)
    kk = A * 4 / P
    return A - min_beta_bb(d, kk, L, green(K, kk), 1e300)


@njit(cache=True)
def _feasible_dp(e, kk, beta, Gb, r, L, WI):
    """Is there a lift with k y_i <= beta for all i (banded kernel Gb[i, r+m-i])?
    Curvatures c_j = e_j + n_j P (|c_j| <= L) for j < K-1, c_{K-1} = -w_{K-2}.
    State: running slope index (w on a grid of step P/... tracked exactly as
    sum) and the last 2r curvature choice indices."""
    K = e.shape[0] + 1
    nopt = int(2 * L / P) + 2
    # options per site: values e_j + n P within [-L, L]
    opt = np.full((K - 1, nopt), np.nan)
    nop = np.zeros(K - 1, np.int64)
    for j in range(K - 1):
        n0 = np.ceil((-L - e[j]) / P)
        q = 0
        while e[j] + (n0 + q) * P <= L and q < nopt:
            opt[j, q] = e[j] + (n0 + q) * P
            q += 1
        nop[j] = q
    H = 2 * r                      # history length
    base = nopt
    nh = 1
    for _ in range(H):
        nh *= base
    # w is a sum of curvatures: w_j = sum_{m<=j} c_m.  Its fractional part
    # mod P is fixed per j, so index w by its integer offset in [-WI, WI].
    nw = 2 * WI + 1
    cur = np.zeros((nw, nh), np.bool_)
    nxt = np.zeros((nw, nh), np.bool_)
    fr = np.zeros(K)               # fractional reference of w_j
    s = 0.0
    for j in range(K - 1):
        s += e[j]
        fr[j] = s
    # initial: choose c_0
    for q in range(nop[0]):
        wi = int(round((opt[0, q] - fr[0]) / P)) + WI
        if 0 <= wi < nw:
            cur[wi, q] = True
    hist = np.empty(H, np.int64)
    for j in range(1, K - 1):
        nxt[:, :] = False
        for wi in range(nw):
            for h in range(nh):
                if not cur[wi, h]:
                    continue
                # decode history: h = sum hist[t] base^t, hist[0] = most recent (c_{j-1})
                hh = h
                for t in range(H):
                    hist[t] = hh % base
                    hh //= base
                for q in range(nop[j]):
                    cj = opt[j, q]
                    # constraint at site i = j - r (all curvatures i-r..i+r known)
                    i = j - r
                    ok = True
                    if i >= 0:
                        acc = Gb[i, r + (j - i)] * cj
                        for t in range(H):
                            m = j - 1 - t
                            if m < 0 or m < i - r:
                                break
                            acc += Gb[i, r + (m - i)] * opt[m, hist[t]]
                        if kk * acc > beta:
                            ok = False
                    if not ok:
                        continue
                    w_new = (wi - WI) * P + fr[j - 1] + cj
                    wn = int(round((w_new - fr[j]) / P)) + WI
                    if wn < 0 or wn >= nw:
                        continue
                    hn = q
                    mult = base
                    for t in range(H - 1):
                        hn += hist[t] * mult
                        mult *= base
                    nxt[wn, hn] = True
        cur, nxt = nxt, cur
    # final curvature c_{K-1} = -w_{K-2}; check the last r+1 sites
    for wi in range(nw):
        for h in range(nh):
            if not cur[wi, h]:
                continue
            hh = h
            for t in range(H):
                hist[t] = hh % base
                hh //= base
            wlast = (wi - WI) * P + fr[K - 2]
            cl = -wlast
            if abs(cl) > L:
                continue
            ok = True
            for i in range(max(0, K - 1 - r), K):
                acc = 0.0
                for m in range(max(0, i - r), min(K, i + r + 1)):
                    if m == K - 1:
                        cm = cl
                    else:
                        t = K - 2 - m
                        if t >= H:
                            continue
                        cm = opt[m, hist[t]]
                    acc += Gb[i, r + (m - i)] * cm
                if kk * acc > beta:
                    ok = False
                    break
            if ok:
                return True
    return False


def fc_dp(phase, A, r=2, L=2 * P, WI=None):
    K = len(phase)
    p = phase + P / 2
    d = np.diff(p)
    e = np.empty(K - 1)
    e[0] = d[0]
    e[1:] = d[1:] - d[:-1]
    kk = A * 4 / P
    G = green(K, kk)
    Gb = np.zeros((K, 2 * r + 1))
    for i in range(K):
        for m in range(max(0, i - r), min(K, i + r + 1)):
            Gb[i, r + m - i] = G[i, m]
    if WI is None:
        WI = int(3 * np.sqrt(K)) + 4
    lo, hi = 0.0, 3 * P
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if _feasible_dp(e, kk, mid, Gb, r, L, WI):
            hi = mid
        else:
            lo = mid
    return A - hi


if __name__ == "__main__":
    out = {}
    for A, Ks in [(3.0, [6, 8, 10]), (6.0, [6, 8, 10])]:
        for K in Ks:
            with Pool(9) as p:
                r = np.array(p.map(work, [(31_000_000 + 1000 * K + s, K, A) for s in range(90)]))
            diff = r[:, 1] - r[:, 0]
            out[f"A{A:g}_K{K}"] = dict(exact=r[:, 0].tolist(), predicted=r[:, 1].tolist())
            print(f"k={A:g} K={K}: <fc exact>={r[:,0].mean():.4f}  <fc predicted>={r[:,1].mean():.4f}  "
                  f"max|diff|={np.abs(diff).max():.2e}  median|diff|={np.median(np.abs(diff)):.2e}", flush=True)
    json.dump(out, open("strong_exact.json", "w"), indent=1)
