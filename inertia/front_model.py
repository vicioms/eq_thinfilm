"""Coarse-grained front model: damped lattice line + phonon-transported temperature.
m u_tt + u_t = u_xx + f - Phi_i ;  d_t e_pm +- c d_x e_pm = -(1/m) e_pm + v_i Phi_i/2 ; T=e_+ + e_-.
Phi_i = F(u_t,T_i) (closure, tabulated) for sliding sites; stuck sites hold until drive > f_s(T).
Usage: python front_model.py  (builds tables, runs scans, writes data/front_*.npz)
"""
import numpy as np, os, sys, json
from numba import njit, prange
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
VG = np.geomspace(1e-3, 1.5, 25)
TG = np.concatenate([[0.0], np.geomspace(2e-3, 4.0, 19)])

def table(m):
    fn = f'{D}/front_table_m{int(m)}.npz'
    if os.path.exists(fn):
        z = np.load(fn); return z['F'], z['fmin'], z['vmin'], z['vs'], z['Tv']
    import sc
    L = sc.Line(m)
    F = np.array([[L.F(v, T) for T in TG] for v in VG])
    vs = np.geomspace(0.01, 0.5, 30)
    fl, Tv, _ = L.flow(vs)
    i = np.argmin(fl)
    np.savez(fn, F=F, fmin=fl[i], vmin=vs[i], vs=vs, fl=fl, Tv=Tv)
    return F, fl[i], vs[i], vs, Tv

@njit(cache=True)
def Fint(F, lv, T, VG, TG):
    # bilinear in log v, T (clamped)
    x = (lv - np.log(VG[0])) / (np.log(VG[-1]) - np.log(VG[0])) * (len(VG) - 1)
    x = min(max(x, 0.0), len(VG) - 1.000001); i = int(x); a = x - i
    j = np.searchsorted(TG, T) - 1
    j = min(max(j, 0), len(TG) - 2)
    b = min(max((T - TG[j]) / (TG[j + 1] - TG[j]), 0.0), 1.0)
    return (1-a)*((1-b)*F[i, j] + b*F[i, j+1]) + a*((1-b)*F[i+1, j] + b*F[i+1, j+1])

@njit(cache=True)
def run(F, VG, TG, m, f, v0, T0hot, N, ihot, dt, nsteps, nrec, fsmode, fsT0, vlo, fc, dx_only=0):
    c = 1.0 / np.sqrt(m)
    u = np.zeros(N); v = np.zeros(N); ep = np.zeros(N); em = np.zeros(N)
    v[:ihot] = v0; ep[:ihot] = 0.5*T0hot; em[:ihot] = 0.5*T0hot
    stuck = np.ones(N, np.bool_); stuck[:ihot] = False
    nr = nsteps // nrec
    front = np.zeros(nr); Tfront = np.zeros(nr); tt = np.zeros(nr)
    Phi = np.zeros(N); g = np.zeros(N)
    cf = c * dt
    for n in range(nsteps):
        for i in range(N):
            T = ep[i] + em[i]
            l = u[i-1] if i > 0 else u[i]
            r = u[i+1] if i < N-1 else u[i]
            g[i] = l - 2*u[i] + r + f
            if stuck[i]:
                if fsmode == 0: fs = fc * np.exp(-T / fsT0)
                elif fsmode == 1: fs = Fint(F, np.log(vlo), T, VG, TG)
                else: fs = fc
                if abs(g[i]) > fs:
                    stuck[i] = False
                    v[i] = 0.0
                    Phi[i] = fs * np.sign(g[i])
                else:
                    Phi[i] = g[i]
            if not stuck[i]:
                av = abs(v[i])
                if av < vlo:
                    if fsmode == 0: fs = fc * np.exp(-T / fsT0)
                    elif fsmode == 1: fs = Fint(F, np.log(vlo), T, VG, TG)
                    else: fs = fc
                    if abs(g[i]) < fs:
                        stuck[i] = True; v[i] = 0.0; Phi[i] = g[i]
                    else:
                        Phi[i] = np.sign(g[i]) * Fint(F, np.log(vlo), T, VG, TG)
                else:
                    Phi[i] = np.sign(v[i]) * Fint(F, np.log(av), T, VG, TG)
        # update
        for i in range(N):
            if stuck[i]:
                v[i] = 0.0
            else:
                v[i] = (v[i] + dt/m * (g[i] - Phi[i])) / (1 + dt/m)
            u[i] += dt * v[i]
        # temperature: heating from sliding sites, damping, upwind advection
        q = np.zeros(N)
        for i in range(N):
            if not stuck[i]:
                q[i] = 0.5 * abs(v[i]) * abs(Phi[i])
        nep = np.empty(N); nem = np.empty(N)
        for i in range(N):
            lp = ep[i-1] if i > 0 else 0.0
            rm = em[i+1] if i < N-1 else 0.0
            nep[i] = ep[i] - cf*(ep[i]-lp) + dt*(-ep[i]/m + q[i])
            nem[i] = em[i] - cf*(em[i]-rm) + dt*(-em[i]/m + q[i])
            # exact-ish decay
        ep = nep; em = nem
        if (n+1) % nrec == 0:
            k = (n+1)//nrec - 1
            ifr = -1
            for i in range(N-1, -1, -1):
                if (not stuck[i]) and abs(v[i]) > 0.0:
                    ifr = i; break
            front[k] = ifr; Tfront[k] = ep[max(ifr, 0)] + em[max(ifr, 0)]; tt[k] = (n+1)*dt
    return tt, front, Tfront, u, v, ep+em, stuck

def speed(m, f, fsmode=0, fsT0=1.0, vlo=1e-3, dt=None, N=None, tmax=None, tab=None, hot_below=True, nrec=20):
    F, fmin, vmin, vs, Tv = tab if tab is not None else table(m)
    c = 1/np.sqrt(m)
    if f >= fmin:  # closure velocity on the upper (stable) branch for f
        fl = np.load(f'{D}/front_table_m{int(m)}.npz')['fl']
        i0 = np.argmin(fl); j = i0 + np.argmax(fl[i0:] >= f) if fl[i0:].max() >= f else len(fl)-1
        v0 = np.interp(f, fl[i0:j+1], vs[i0:j+1]) if j > i0 else vmin
        T0 = np.interp(v0, vs, Tv)
    else:
        v0, T0 = vmin, np.interp(vmin, vs, Tv)
    dt = dt or 0.5*np.sqrt(m)/2
    tmax = tmax or 40*m
    N = N or int(1.2*(1.2*c*tmax) + 400)
    nsteps = int(tmax/dt)//nrec*nrec
    ihot = 300
    tt, fr, Tf, u, v, T, st = run(F, VG, TG, float(m), f, v0, T0, N, ihot, dt, nsteps, nrec, fsmode, fsT0, vlo, 0.308)
    # fit front slope over middle 20%-90% of the run until edges
    ok = (fr > 0) & (fr < N-60)
    sel = ok & (tt > 0.2*tmax)
    if sel.sum() < 5: return np.nan, tt, fr
    p = np.polyfit(tt[sel], fr[sel], 1)
    return p[0], tt, fr

if __name__ == '__main__':
    for m in (256, 1024):
        table(m)
        print('table', m, table(m)[1:3])

def scan():
    rs = np.array([0.9,0.95,1.0,1.02,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.15,1.2,1.25])
    var = {'base': {}, 'vlo1e-4': dict(vlo=1e-4), 'vlo1e-2': dict(vlo=1e-2), 'fsexp0.3': dict(fsT0=0.3),
           'fs_closure': dict(fsmode=1), 'fs_const': dict(fsmode=2)}
    for m in (256, 1024):
        tab = table(m); out = {'r': rs, 'fmin': tab[1], 'c': m**-.5}
        for k, kw in var.items():
            out[k] = np.array([speed(m, r*tab[1], tab=tab, tmax=15*m, **kw)[0] for r in rs]) / m**-.5
        np.savez(f'{D}/front_scan_m{m}.npz', **out)
