"""Annealed (Kac-Rice) quantities for the triangle-wave chain of qcount.py.

Exact finite-K mean number of equilibria (kappa = 1):

    <N_K(f)> = P * p_{S_K}(K f) * <|D_K|>,

S_K = sum of K iid uniform[-A, A] (Irwin-Hall), D_K = det(-Delta_N + U),
U_i = +-k iid fair signs, k = 4A/P.  D_K is the second component of
M_{K-1} ... M_0 (1, 0)^T with M_U = [[1+U, 1], [U, 1]] acting on
(psi_i, phi_{i-1}).

Rates:  Lambda(q) = lim (1/K) ln <|D_K|^q>  (generalised Lyapunov exponent),
        I(f) = Cramer function of uniform[-A, A],
        Sigma_ann(f) = Lambda(1) - I(f).
"""
import numpy as np
import mpmath as mp
from scipy.optimize import brentq


def _mats(k):
    return [np.array([[1 + U, 1.0], [U, 1.0]]) for U in (k, -k)]


def _angle_ops(k, n):
    """Grid on the projective line theta in [0, pi); for each U return the
    image angle and the stretching factor of the unit vector v_theta."""
    th = (np.arange(n) + 0.5) * np.pi / n
    v = np.stack([np.cos(th), np.sin(th)])
    out = []
    for M in _mats(k):
        x = M @ v
        r = np.hypot(x[0], x[1])
        th2 = np.mod(np.arctan2(x[1], x[0]), np.pi)
        out.append((th2, r))
    return th, out


def _interp(h, th2, n):
    # periodic linear interpolation of h (sampled at cell centres) at th2
    x = th2 / (np.pi / n) - 0.5
    i0 = np.floor(x).astype(int)
    t = x - i0
    return (1 - t) * h[i0 % n] + t * h[(i0 + 1) % n]


def mean_abs_det(K, k, q=1.0, n=20000):
    """<|D_K|^q> by backward transfer on the angle (exact up to grid)."""
    th, ops = _angle_ops(k, n)
    H = np.abs(np.sin(th)) ** q  # h_K(v) = |v_2|^q, homogeneous of degree q
    logscale = 0.0
    for _ in range(K):
        Hn = 0.5 * sum(r ** q * _interp(H, th2, n) for th2, r in ops)
        s = Hn.max()
        H = Hn / s
        logscale += np.log(s)
    # evaluate at theta = 0 (initial vector (1, 0))
    h0 = _interp(H, np.array([0.0]), n)[0]
    return np.exp(logscale) * h0, logscale + np.log(h0)


def mean_abs_det_enum(K, k, q=1.0):
    """Exact enumeration over 2^K sign sequences (small K check)."""
    v = np.array([[1.0], [0.0]])
    Ms = _mats(k)
    for _ in range(K):
        v = np.concatenate([M @ v for M in Ms], axis=1)
    return np.mean(np.abs(v[1]) ** q)


def gle(q, k, n=4000, iters=4000, tol=1e-13):
    """Lambda(q) = log of the top eigenvalue of the angle transfer operator."""
    th, ops = _angle_ops(k, n)
    H = np.ones(n)
    lam = 0.0
    for it in range(iters):
        Hn = 0.5 * sum(r ** q * _interp(H, th2, n) for th2, r in ops)
        s = Hn.max()
        Hn /= s
        if np.max(np.abs(Hn - H)) < tol:
            H = Hn
            lam = np.log(s)
            break
        H = Hn
        lam = np.log(s)
    return lam


def lyapunov(k, K=2_000_000, seed=0):
    """Furstenberg Lyapunov exponent gamma = Lambda'(0) by direct iteration."""
    rng = np.random.default_rng(seed)
    U = np.where(rng.random(K) < 0.5, k, -k)
    x, y = 1.0, 0.0
    s = 0.0
    for Ui in U:
        x, y = (1 + Ui) * x + y, Ui * x + y
        r = np.hypot(x, y)
        s += np.log(r)
        x /= r
        y /= r
    return s / K


def irwin_hall_density(K, A, x, dps=60):
    """Density of S_K = sum of K iid uniform[-A, A] at x (exact)."""
    mp.mp.dps = dps
    # rescale to sum of K uniform[0,1]: y = (x + K A) / (2A)
    y = mp.mpf(x + K * A) / (2 * A)
    if y <= 0 or y >= K:
        return 0.0
    s = mp.mpf(0)
    for j in range(int(mp.floor(y)) + 1):
        s += (-1) ** j * mp.binomial(K, j) * (y - j) ** (K - 1)
    s /= mp.factorial(K - 1)
    return float(s / (2 * A))


def cramer_uniform(f, A):
    """I(f) = sup_l [l f - ln(sinh(l A)/(l A))] for uniform[-A, A]."""
    x = abs(f) / A
    if x == 0:
        return 0.0
    if x >= 1:
        return np.inf
    # solve coth(t) - 1/t = x
    def lang(t):
        return t / 3 - t ** 3 / 45 if t < 1e-3 else 1 / np.tanh(t) - 1 / t
    t = brentq(lambda t: lang(t) - x, 1e-12, 1e9)
    # ln(sinh t / t), stable for large t
    lst = t - np.log(2 * t) + np.log1p(-np.exp(-2 * t)) if t > 1e-3 else t * t / 6
    return t * x - lst


def annealed_mean(K, A, P, f, n=20000):
    k = 4 * A / P
    D, logD = mean_abs_det(K, k, 1.0, n)
    p = irwin_hall_density(K, A, K * f)
    return P * p * D


def sigma_ann(f, A, P, lam1=None):
    k = 4 * A / P
    if lam1 is None:
        lam1 = gle(1.0, k)
    return lam1 - cramer_uniform(f, A)


def f_ann(A, P, lam1=None):
    k = 4 * A / P
    if lam1 is None:
        lam1 = gle(1.0, k)
    return brentq(lambda f: lam1 - cramer_uniform(f, A), 1e-9, A * (1 - 1e-8))


def gle_stable(k, n=6000, iters=20000, tol=1e-13):
    """Annealed growth rate of <|D_K| 1(all leading minors > 0)>, i.e. of the
    mean number of local minima (up to the force factor).  Oriented vectors
    (psi, phi) with psi > 0 are parametrised by theta in (-pi/2, pi/2); a step
    that makes psi_{j+1} <= 0 is killed (Sturm: a negative eigenvalue)."""
    th = -0.5 * np.pi + (np.arange(n) + 0.5) * np.pi / n
    v = np.stack([np.cos(th), np.sin(th)])
    ops = []
    for M in _mats(k):
        x = M @ v
        r = np.hypot(x[0], x[1])
        alive = x[0] > 0
        th2 = np.arctan2(x[1], x[0])
        ops.append((th2, r * alive))
    def interp(h, t2):
        xx = (t2 + 0.5 * np.pi) / (np.pi / n) - 0.5
        xx = np.clip(xx, 0, n - 1)
        i0 = np.minimum(np.floor(xx).astype(int), n - 2)
        t = xx - i0
        return (1 - t) * h[i0] + t * h[i0 + 1]
    H = np.ones(n)
    lam = 0.0
    for it in range(iters):
        Hn = 0.5 * sum(r * interp(H, t2) for t2, r in ops)
        s = Hn.max()
        Hn /= s
        lam = np.log(s)
        if np.max(np.abs(Hn - H)) < tol:
            break
        H = Hn
    return lam


def gle_curve(k, qs, n=4000):
    return np.array([gle(q, k, n) for q in qs])


# ---------------------------------------------------------------------------
# Piecewise-linear interpolation of FLRT's Gaussian disorder (qcount.cos_nodes)
# ---------------------------------------------------------------------------

def pl_cos_site_law(sigma, M, nR=48, nth=96, nt=24, P=2 * np.pi):
    """Quadrature for the joint law of (V', V'') at a uniform point for the
    PL interpolation (M nodes per period) of V' = R cos(u - alpha), R Rayleigh.
    Returns weights, V', U arrays (weights sum to 1)."""
    h = P / M
    # Rayleigh(sigma): R = sigma sqrt(2 E), E ~ Exp(1); Gauss-Laguerre in E
    xe, we = np.polynomial.laguerre.laggauss(nR)
    R = sigma * np.sqrt(2 * xe)
    th = (np.arange(nth) + 0.5) * 2 * np.pi / nth
    t = (np.arange(nt) + 0.5) / nt
    Rg, Tg, tg = np.meshgrid(R, th, t, indexing="ij")
    W = np.broadcast_to(we[:, None, None], Rg.shape) / (nth * nt)
    Vp = Rg * ((1 - tg) * np.cos(Tg) + tg * np.cos(Tg + h))
    U = Rg * (np.cos(Tg + h) - np.cos(Tg)) / h
    return W.ravel(), Vp.ravel(), U.ravel()


def gle_tilted(law, tau, nU=240, n=3000, iters=3000, tol=1e-12):
    """Lambda(tau) = lim (1/K) ln < exp(tau sum V') |D_K| > for iid (V', U)
    pairs with joint quadrature law (W, V', U): top eigenvalue of the angle
    transfer operator with U binned (bin weights carry exp(tau V'))."""
    W, Vp, U = law
    shift = tau * Vp.max() if tau > 0 else 0.0
    w = W * np.exp(tau * Vp - shift)   # log-shifted to avoid overflow
    edges = np.quantile(U, np.linspace(0, 1, nU + 1))
    idx = np.clip(np.searchsorted(edges, U, side="right") - 1, 0, nU - 1)
    wb = np.bincount(idx, w, nU)
    Ub = np.bincount(idx, w * U, nU) / np.maximum(wb, 1e-300)
    th = (np.arange(n) + 0.5) * np.pi / n
    v = np.stack([np.cos(th), np.sin(th)])
    ops = []
    for Ui, wi in zip(Ub, wb):
        if wi <= 0:
            continue
        M_ = np.array([[1 + Ui, 1.0], [Ui, 1.0]])
        x = M_ @ v
        r = np.hypot(x[0], x[1])
        ops.append((np.mod(np.arctan2(x[1], x[0]), np.pi), wi * r))
    H = np.ones(n)
    lam = 0.0
    for _ in range(iters):
        Hn = sum(r * _interp(H, t2, n) for t2, r in ops)
        s = Hn.max()
        Hn /= s
        lam = np.log(s)
        if np.max(np.abs(Hn - H)) < tol:
            break
        H = Hn
    return lam + shift


def sigma_ann_pl(v, law, taus=None):
    """Sigma_ann(v) = inf_tau [Lambda(tau) - tau v] (exact annealed rate)."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda t: gle_tilted(law, t) - t * v,
                          bounds=(0.0, 60.0), method="bounded",
                          options={"xatol": 1e-4})
    return res.fun, res.x


def amplitude_quenched_rate(sigma, M, tau=0.0, n=1000, nth=64, nt=16, K=400_000,
                            seed=0, P=2 * np.pi, nRgrid=256):
    """'Phase-annealed, amplitude-quenched' rate for the PL Gaussian model:

        Lambda_q(tau) = lim (1/K) < ln || T_{R_K} ... T_{R_1} || >,

    T_R the angle transfer operator averaged over the phase (and position in
    the cell) at fixed site amplitude R, with weight exp(tau V').  The
    amplitudes R_i ~ Rayleigh(sigma) are shared by every branch at site i;
    the phases are averaged by a population spread uniformly in u.  For a
    fixed amplitude (triangle model) this equals the annealed Lambda(1).
    Returns (Lambda_q, its standard error from 20 blocks, Lambda_ann check).
    """
    rng = np.random.default_rng(seed)
    h = P / M
    th_c = (np.arange(nth) + 0.5) * 2 * np.pi / nth        # cell phase
    tt = (np.arange(nt) + 0.5) / nt                         # position in cell
    c = (np.cos(th_c + h) - np.cos(th_c)) / h               # U / R
    vp = (1 - tt[None, :]) * np.cos(th_c)[:, None] + tt[None, :] * np.cos(th_c + h)[:, None]  # V'/R
    ang = (np.arange(n) + 0.5) * np.pi / n
    v = np.stack([np.cos(ang), np.sin(ang)])
    # R grid: Rayleigh quantiles (stratified sampling of R by index)
    qs = (np.arange(nRgrid) + 0.5) / nRgrid
    Rg = sigma * np.sqrt(-2 * np.log(1 - qs))
    ops = []
    for R in Rg:
        idx_list, wl_list, wr_list, r_list = [], [], [], []
        wtau = np.exp(tau * R * vp).mean(1)                 # E_t e^{tau V'} per cell phase
        for U, wt in zip(R * c, wtau):
            Mx = np.array([[1 + U, 1.0], [U, 1.0]]) @ v
            r = np.hypot(Mx[0], Mx[1]) * wt / nth
            th2 = np.mod(np.arctan2(Mx[1], Mx[0]), np.pi)
            x = th2 / (np.pi / n) - 0.5
            i0 = np.floor(x).astype(int)
            t = x - i0
            idx_list.append(i0); wl_list.append((1 - t) * r); wr_list.append(t * r)
        ops.append((np.array(idx_list), np.array(wl_list), np.array(wr_list)))

    def apply(op, H):
        i0, wl, wr = op
        return (wl * H[i0 % n] + wr * H[(i0 + 1) % n]).sum(0)

    # annealed check: average operator over the R grid
    H = np.ones(n)
    for _ in range(300):
        Hn = sum(apply(op, H) for op in ops) / nRgrid
        s = Hn.max(); H = Hn / s
    lam_ann = np.log(s)
    # quenched product (backward iteration of h, equivalent rate)
    H = np.ones(n)
    for _ in range(200):
        H = apply(ops[rng.integers(nRgrid)], H); H /= H.max()
    logs = np.empty(K)
    for i in range(K):
        H = apply(ops[rng.integers(nRgrid)], H)
        s = H.max(); H /= s; logs[i] = np.log(s)
    blocks = logs.reshape(20, -1).mean(1)
    return logs.mean(), blocks.std() / np.sqrt(20), lam_ann
