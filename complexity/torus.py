"""Quenched pressure from the population dynamics on the torus.

In the zero-force frame the shooting curve evolves under the random twist maps
T_i(u, w) = (u + w', w'), w' = w + V_i'(u).  Only (u mod P, w mod P) matter
for the future, so the tilted curve-length measure

    rho_i = sum over curve points of e^{tau W_i} |tangent|  (weights)

lives on the torus [0, P)^2 and Z_K(tau) = Z_0 prod_i m_i with
m_i = |rho_{i+1}| / |rho_i|.  The quenched pressure is the time average

    phi_q(tau) = lim_K K^-1 sum_i ln m_i     (one long run, ergodicity).

Discretisation: rho on an n x n grid, transported by the exact map with
bilinear deposit (sub-cell sampling s x s).  The curve's tangent direction is
carried as a line field (doubled angle, weighted average per cell): strands
through a point align with the unstable direction, so one direction per
cell suffices.  The grid smooths structure below h = P/n, i.e. it acts as
weak noise; results are checked for convergence in n.

    from torus import pressure
    phi, stderr = pressure(A=0.6, tau=0.0, n=256, steps=4000)
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _run(A, P, tau, n, s, steps, burn, seed, record_loc):
    np.random.seed(seed)
    k = 4.0 * A / P
    half = 0.5 * P
    h = P / n
    rho = np.ones((n, n)) / (n * n)
    c2 = np.ones((n, n))     # cos(2 theta) of the line field, theta = 0 (along u)
    s2 = np.zeros((n, n))
    nrho = np.empty((n, n))
    nc2 = np.empty((n, n))
    ns2 = np.empty((n, n))
    logs = np.empty(steps)
    loc = np.empty(steps)
    sub = (np.arange(s) + 0.5) / s
    for it in range(steps):
        ph = np.random.random() * P
        nrho[:, :] = 0.0
        nc2[:, :] = 0.0
        ns2[:, :] = 0.0
        tot = 0.0
        for i in range(n):
            for j in range(n):
                r0 = rho[i, j]
                if r0 <= 0.0:
                    continue
                # line direction of the cell
                th = 0.5 * np.arctan2(s2[i, j], c2[i, j])
                cx = np.cos(th)
                cy = np.sin(th)
                for a in range(s):
                    u = (i + sub[a]) * h
                    x = (u - ph) / half
                    mm = np.floor(x)
                    t = x - mm
                    if (int(mm) & 1) == 0:
                        g = -A + 2.0 * A * t
                        U = k
                    else:
                        g = A - 2.0 * A * t
                        U = -k
                    vx = (1.0 + U) * cx + cy
                    vy = U * cx + cy
                    rr = np.sqrt(vx * vx + vy * vy)
                    wgt0 = r0 / (s * s) * np.exp(tau * g) * rr
                    th2 = np.arctan2(vy, vx)
                    cc = np.cos(2 * th2)
                    ss = np.sin(2 * th2)
                    for b in range(s):
                        w = (j + sub[b]) * h
                        w2 = w + g
                        u2 = u + w2
                        # bilinear deposit on the periodic grid (cell centres)
                        fu = u2 / h - 0.5
                        fw = w2 / h - 0.5
                        iu = int(np.floor(fu))
                        iw = int(np.floor(fw))
                        du = fu - iu
                        dw = fw - iw
                        iu0 = iu % n
                        iw0 = iw % n
                        iu1 = (iu + 1) % n
                        iw1 = (iw + 1) % n
                        for (p, q, wb) in ((iu0, iw0, (1 - du) * (1 - dw)), (iu1, iw0, du * (1 - dw)),
                                           (iu0, iw1, (1 - du) * dw), (iu1, iw1, du * dw)):
                            ww = wgt0 * wb
                            nrho[p, q] += ww
                            nc2[p, q] += ww * cc
                            ns2[p, q] += ww * ss
                        tot += wgt0
        logs[it] = np.log(tot)
        # localisation of the population in u: participation of the u-marginal
        if record_loc:
            pu = nrho.sum(axis=1) / tot
            loc[it] = n * (pu * pu).sum()   # 1 if uniform in u, up to n if localised
        for i in range(n):
            for j in range(n):
                rho[i, j] = nrho[i, j] / tot
                c2[i, j] = nc2[i, j]
                s2[i, j] = ns2[i, j]
    return logs, loc


def pressure(A, tau, n=256, s=2, steps=3000, burn=300, seed=1, record_loc=False):
    logs, loc = _run(A, 4.0, tau, n, s, steps, burn, seed, record_loc)
    x = logs[burn:]
    nb = 20
    bl = x[: len(x) // nb * nb].reshape(nb, -1).mean(1)
    out = (x.mean(), bl.std() / np.sqrt(nb))
    if record_loc:
        return out + (loc[burn:].mean(),)
    return out


@njit(cache=True)
def _growth_over_phases(rho, c2, s2, A, P, tau, nphi):
    """m(phi) = sum_cells rho e^{tau g(u - phi)} |M_U(u-phi) v_theta| for a
    quadrature of nphi phases: the one-step growth factor as a function of the
    fresh phase, at fixed population (rho, line field)."""
    n = rho.shape[0]
    k = 4.0 * A / P
    half = 0.5 * P
    h = P / n
    out = np.zeros(nphi)
    for q in range(nphi):
        ph = (q + 0.5) * P / nphi
        tot = 0.0
        for i in range(n):
            u = (i + 0.5) * h
            x = (u - ph) / half
            mm = np.floor(x)
            t = x - mm
            if (int(mm) & 1) == 0:
                g = -A + 2.0 * A * t
                U = k
            else:
                g = A - 2.0 * A * t
                U = -k
            eg = np.exp(tau * g)
            for j in range(n):
                r0 = rho[i, j]
                if r0 <= 0.0:
                    continue
                th = 0.5 * np.arctan2(s2[i, j], c2[i, j])
                cx = np.cos(th)
                cy = np.sin(th)
                vx = (1.0 + U) * cx + cy
                vy = U * cx + cy
                tot += r0 * eg * np.sqrt(vx * vx + vy * vy)
        out[q] = tot
    return out


@njit(cache=True)
def _run_decomp(A, P, tau, n, s, steps, seed, nphi, every):
    """As _run, and every `every` steps records ln <m>_phi and <ln m>_phi
    over a phase quadrature (exact conditional averages at fixed past)."""
    np.random.seed(seed)
    k = 4.0 * A / P
    half = 0.5 * P
    h = P / n
    rho = np.ones((n, n)) / (n * n)
    c2 = np.ones((n, n))
    s2 = np.zeros((n, n))
    nrho = np.empty((n, n)); nc2 = np.empty((n, n)); ns2 = np.empty((n, n))
    logs = np.empty(steps)
    nrec = steps // every
    lnE = np.empty(nrec); Eln = np.empty(nrec); var2 = np.empty(nrec)
    r = 0
    sub = (np.arange(s) + 0.5) / s
    for it in range(steps):
        if it % every == 0 and r < nrec:
            mphi = _growth_over_phases(rho, c2, s2, A, P, tau, nphi)
            Em = mphi.mean()
            lnE[r] = np.log(Em)
            Eln[r] = np.log(mphi).mean()
            var2[r] = 0.5 * ((mphi / Em - 1.0) ** 2).mean()
            r += 1
        ph = np.random.random() * P
        nrho[:, :] = 0.0; nc2[:, :] = 0.0; ns2[:, :] = 0.0
        tot = 0.0
        for i in range(n):
            for j in range(n):
                r0 = rho[i, j]
                if r0 <= 0.0:
                    continue
                th = 0.5 * np.arctan2(s2[i, j], c2[i, j])
                cx = np.cos(th); cy = np.sin(th)
                for a in range(s):
                    u = (i + sub[a]) * h
                    x = (u - ph) / half
                    mm = np.floor(x)
                    t = x - mm
                    if (int(mm) & 1) == 0:
                        g = -A + 2.0 * A * t; U = k
                    else:
                        g = A - 2.0 * A * t; U = -k
                    vx = (1.0 + U) * cx + cy
                    vy = U * cx + cy
                    rr = np.sqrt(vx * vx + vy * vy)
                    wgt0 = r0 / (s * s) * np.exp(tau * g) * rr
                    th2 = np.arctan2(vy, vx)
                    cc = np.cos(2 * th2); ss = np.sin(2 * th2)
                    for b in range(s):
                        w = (j + sub[b]) * h
                        w2 = w + g
                        u2 = u + w2
                        fu = u2 / h - 0.5; fw = w2 / h - 0.5
                        iu = int(np.floor(fu)); iw = int(np.floor(fw))
                        du = fu - iu; dw = fw - iw
                        iu0 = iu % n; iw0 = iw % n; iu1 = (iu + 1) % n; iw1 = (iw + 1) % n
                        for (p, q, wb) in ((iu0, iw0, (1 - du) * (1 - dw)), (iu1, iw0, du * (1 - dw)),
                                           (iu0, iw1, (1 - du) * dw), (iu1, iw1, du * dw)):
                            ww = wgt0 * wb
                            nrho[p, q] += ww; nc2[p, q] += ww * cc; ns2[p, q] += ww * ss
                        tot += wgt0
        logs[it] = np.log(tot)
        for i in range(n):
            for j in range(n):
                rho[i, j] = nrho[i, j] / tot
                c2[i, j] = nc2[i, j]; s2[i, j] = ns2[i, j]
    return logs, lnE[:r], Eln[:r], var2[:r]


def decomposition(A, tau, n=192, s=2, steps=2500, burn=300, seed=3, nphi=48, every=5):
    """phi_q, the memory term <Lambda - ln<m>_phi>, the one-step Jensen term
    <ln<m>_phi - <ln m>_phi> and its second-order estimate <Var_phi m / 2<m>^2>."""
    logs, lnE, Eln, var2 = _run_decomp(A, 4.0, tau, n, s, steps, seed, nphi, every)
    b = burn // every
    return dict(phi=float(logs[burn:].mean()), lnEm=float(lnE[b:].mean()),
                jensen=float((lnE[b:] - Eln[b:]).mean()), jensen2=float(var2[b:].mean()))
