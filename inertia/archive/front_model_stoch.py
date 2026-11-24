"""Stochastic, spatially resolved front model (blocks of size ell = c*m = sqrt(m)).
Per block i: phonon energies e_pm (T = e_+ + e_-), no displacement field.
  closure:  v_i = largest root of  f = v + F(v,T_i)   (stuck, v=0, if no root: T_i < T_c(f))
            heating  q_i = v_i F(v_i,T_i)   (power balance T = m q)
  d e_pm = [ -/+ (c/ell)(e_pm,i - e_pm,i-+1) - e_pm/m + q_i/2 ] dt  (upwind)
           + (T_i/2) sqrt(2 dt/(nu m)) xi_i      (common noise on both: multiplicative, from T itself)
Release of stuck blocks by active neighbours = their noise is proportional to the (invading) local T.
Quenched block force offsets: f_i = f + d_i, d_i ~ N(0, sq^2) (discretised on 13 levels), the local closure uses f_i.
Free parameters: nu (modes per block; noise variance T^2/nu), sq (std of quenched block force offsets);
ell=sqrt(m) (heat-attenuation length c*m) and dt are fixed.
Usage: python front_model_stoch.py {table|fmin|edge|bulk|all}
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np, os, sys, json
from numba import njit
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

def table(m, nv=56, nT=56):
    fn = f'{D}/fs_table_m{m}.npz'
    if os.path.exists(fn):
        z = np.load(fn); return z['F'], z['VG'], z['TG']
    import sc
    L = sc.Line(m)
    VG = np.geomspace(1e-3, 1.5, nv); TG = np.concatenate([[0], np.geomspace(2e-3, 4, nT)])
    F = np.array([[L.F(v, T) for T in TG] for v in VG])
    np.savez(fn, F=F, VG=VG, TG=TG); return F, VG, TG

def local_tables(m, f, nT=1500, nv=1500):  # single f
    """For given f: on T grid, v_root(T) (0 if none) and q(T)=v F.  Returns TT, v, q, Tc."""
    F, VG, TG = table(m)
    from scipy.interpolate import RegularGridInterpolator as RGI
    lv = np.log(VG)
    itp = RGI((lv, TG), F, bounds_error=False, fill_value=None)
    TT = np.concatenate([[0], np.geomspace(1e-3, 4, nT)])
    vv = np.geomspace(1e-3, 1.4, nv)
    v = np.zeros(len(TT)); q = np.zeros(len(TT))
    for j, T in enumerate(TT):
        Fv = itp(np.column_stack([np.log(vv), np.full(nv, T)]))
        g = vv + Fv - f
        idx = np.where((g[:-1] < 0) & (g[1:] >= 0))[0]   # upward crossing = stable (increasing) branch
        if len(idx):
            i = idx[-1]; w = -g[i] / (g[i+1] - g[i])
            v[j] = vv[i] + w*(vv[i+1]-vv[i]); q[j] = v[j] * (Fv[i] + w*(Fv[i+1]-Fv[i]))
    sl = np.where(v > 0)[0]
    Tc = TT[sl[0]] if len(sl) else np.inf
    return TT, v, q, Tc

def fixed_points(m, f):
    TT, v, q, Tc = local_tables(m, f)
    G = q - TT/m
    s = np.where((G[:-1] > 0) & (G[1:] <= 0) & (v[:-1] > 0))[0]
    return TT, G, Tc, TT[s]

def fmin_model(m):
    lo, hi = 0.12*(64/m)**.5, 0.3*(64/m)**.5
    for _ in range(14):
        mid = 0.5*(lo+hi)
        if len(fixed_points(m, mid)[3]): hi = mid
        else: lo = mid
    return 0.5*(lo+hi)

def level_tables(m, f, sq, nl=13, nT=700, nv=800):
    """tables for f+d_k, d_k = linspace(-3,3,nl)*sq; returns TT, V[k,:], Q[k,:], weights w_k"""
    ds = np.linspace(-3, 3, nl)*sq
    w = np.exp(-0.5*(ds/max(sq,1e-9))**2) if sq > 0 else np.ones(nl)
    if sq == 0: ds = np.zeros(1); w = np.ones(1)
    V = []; Q = []
    for d in ds:
        TT, v, q, Tc = local_tables(m, f+d, nT, nv); V.append(v); Q.append(q)
    return TT, np.array(V), np.array(Q), w/w.sum()

def make_kidx(w, N, seed):
    return np.random.RandomState(1000+seed).choice(len(w), size=N, p=w).astype(np.int64)

@njit(cache=True)
def lut(TT, a, T):
    # nearest-lower index in a log grid (TT[0]=0, TT[1]=1e-3, geometric to 4)
    if T <= TT[1]: return a[0] if T < 0.5*TT[1] else a[1]
    n = len(TT) - 1
    x = np.log(T/TT[1]) / np.log(TT[-1]/TT[1]) * (n - 1)
    i = min(int(x), n - 2) + 1
    w = x - (i - 1)
    return a[i]*(1-w) + a[i+1]*w if a[i] > 0 and a[i+1] > 0 else a[i]

@njit(cache=True)
def sim(TT, vtab, qtab, kidx, m, nu, N, nsteps, nrec, T0, init, seed, dtf, bc):
    """init: 0 = all hot at T0; 1 = left half hot (T0), right cold; bc 0 periodic, 1 open.  Returns time, active count, activity trace, mean v, final T."""
    np.random.seed(seed)
    ell = np.sqrt(m); c = 1.0/np.sqrt(m); dt = dtf*m; r = c*dt/ell
    ep = np.zeros(N); em = np.zeros(N)
    nh = N if init == 0 else N//2
    ep[:nh] = 0.5*T0; em[:nh] = 0.5*T0
    nr = nsteps // nrec
    act = np.zeros(nr); vm = np.zeros(nr); fr = np.zeros(nr)
    sig = np.sqrt(2*dt/(nu*m))
    q = np.zeros(N); nep = np.zeros(N); nem = np.zeros(N)
    for n in range(nsteps):
        for i in range(N):
            q[i] = lut(TT, qtab[kidx[i]], ep[i]+em[i])
        for i in range(N):
            if bc == 0:
                lp = ep[i-1]; rm = em[(i+1) % N]
            else:
                lp = ep[i-1] if i > 0 else 0.0
                rm = em[i+1] if i < N-1 else 0.0
            T = ep[i]+em[i]
            xi = np.random.randn()*sig*0.5*T
            nep[i] = max(ep[i] - r*(ep[i]-lp) + dt*(-ep[i]/m + 0.5*q[i]) + xi, 0.0)
            nem[i] = max(em[i] - r*(em[i]-rm) + dt*(-em[i]/m + 0.5*q[i]) + xi, 0.0)
        for i in range(N):
            ep[i] = nep[i]; em[i] = nem[i]
        if (n+1) % nrec == 0:
            k = (n+1)//nrec - 1
            a = 0; s = 0.0; last = -1
            for i in range(N):
                if q[i] > 0:
                    a += 1; s += lut(TT, vtab[kidx[i]], ep[i]+em[i]); last = i
            act[k] = a/N; vm[k] = s/N; fr[k] = last
    return act, vm, fr, ep+em


@njit(cache=True)
def Fval(F, VG, TG, v, T):
    """bilinear in (log v, T) of the closure table; returns F and dF/dv."""
    lv0 = np.log(VG[0]); lv1 = np.log(VG[-1]); n = len(VG)
    x = (np.log(max(v, VG[0])) - lv0)/(lv1 - lv0)*(n - 1)
    x = min(max(x, 0.0), n - 1.000001); i = int(x); a = x - i
    T = min(T, TG[-1]*0.999999)
    j = np.searchsorted(TG, T) - 1
    j = min(max(j, 0), len(TG) - 2)
    b = (T - TG[j])/(TG[j+1] - TG[j])
    F0 = (1-b)*F[i, j] + b*F[i, j+1]; F1 = (1-b)*F[i+1, j] + b*F[i+1, j+1]
    dlv = (lv1 - lv0)/(n - 1)
    return (1-a)*F0 + a*F1, (F1 - F0)/dlv/max(v, VG[0])

@njit(cache=True)
def sim2(F, VG, TG, m, f, nu, N, nsteps, nrec, T0, v0, init, seed, dtf, K):
    """Blocks i=0..N-1 (periodic), size ell=sqrt(m).  State u,v,e+,e-.
    f_i = f + K (u_{i-1}+u_{i+1}-2u_i),  K = kappa/ell^2 = 1/m
    m dv = (f_i - v - F(v,T)) dt + sqrt(2 T dt/nu) xi   (v>=0)       <- noise from the local temperature
    dT/dt = v F - T/m (+ upwind phonon transport at c)               <- power balance
    init 0: all hot (v0,T0); 1: first half hot, rest stuck cold."""
    np.random.seed(seed)
    ell = np.sqrt(m); c = 1.0/np.sqrt(m); dt = dtf*m; r = c*dt/ell
    u = np.zeros(N); v = np.zeros(N); ep = np.zeros(N); em = np.zeros(N)
    nh = N if init == 0 else N//2
    v[:nh] = v0; ep[:nh] = 0.5*T0; em[:nh] = 0.5*T0
    for i in range(nh):
        u[i] = v0*0.0
    nr = nsteps // nrec
    act = np.zeros(nr); vm = np.zeros(nr); Tm = np.zeros(nr)
    q = np.zeros(N); nep = np.zeros(N); nem = np.zeros(N); fl = np.zeros(N)
    vthr = 0.3*v0 if v0 > 0 else 0.02
    for n in range(nsteps):
        for i in range(N):
            fl[i] = f + K*(u[(i-1) % N] + u[(i+1) % N] - 2*u[i])
        for i in range(N):
            T = ep[i] + em[i]
            Fv, dF = Fval(F, VG, TG, v[i], T)
            q[i] = v[i]*Fv
            den = 1.0 + dt/m*max(1.0 + dF, 0.0)
            vn = v[i] + dt/m*(fl[i] - v[i] - Fv)/den + np.sqrt(2*T*dt/nu)/m*np.random.randn()
            v[i] = max(vn, 0.0)
        for i in range(N):
            u[i] += dt*v[i]
        for i in range(N):
            lp = ep[(i-1) % N]; rm = em[(i+1) % N]
            T = ep[i] + em[i]
            nep[i] = ep[i] - r*(ep[i]-lp) + dt*(-ep[i]/m + 0.5*q[i])
            nem[i] = em[i] - r*(em[i]-rm) + dt*(-em[i]/m + 0.5*q[i])
        for i in range(N):
            ep[i] = nep[i]; em[i] = nem[i]
        if (n+1) % nrec == 0:
            k = (n+1)//nrec - 1
            a = 0; s = 0.0; ts = 0.0
            for i in range(N):
                if v[i] > vthr: a += 1
                s += v[i]; ts += ep[i] + em[i]
            act[k] = a/N; vm[k] = s/N; Tm[k] = ts/N
    return act, vm, Tm, v, ep+em


FMIN = {64: 0.20076, 256: 0.139577}

def edge_speed(m, nu, r, nseed=3, N=3000, dtf=0.02, tmax=70, seed0=0):
    """s/c of a hot|cold front (periodic: two fronts, hot half initially), from the active-block count slope, t in [8,60]m."""
    F, VG, TG = table(m); f = r*FMIN[m]; c = m**-.5; ell = m**.5
    fp = fixed_points(m, f)[3]; T0 = float(fp[-1]) if len(fp) else 0.58*(64/m)**-0.1
    v0 = 0.07*(64/m)**0.5 if r < 1.03 else 0.1*(64/m)**0.5
    nrec = 50; ss = []
    for sd in range(seed0, seed0+nseed):
        act, vm, Tm, v, T = sim2(F, VG, TG, float(m), f, nu, N, int(tmax/dtf), nrec, T0, v0, 1, sd, dtf, 1/m)
        t = (np.arange(len(act))+1)*nrec*dtf*m; n = act*N
        sel = (t > 8*m) & (t < 60*m) & (n > 0) & (n < 0.95*N)
        ss.append(np.polyfit(t[sel], n[sel], 1)[0]*ell/2/c if sel.sum() > 5 else (-30.0 if n[-1] == 0 else np.nan))
    return np.array(ss)

def bulk(m, nu, r, N=4096, tmax=3000, dtf=0.02, seed=0, nrec=None):
    F, VG, TG = table(m); f = r*FMIN[m]
    fp = fixed_points(m, f)[3]; T0 = float(fp[-1]) if len(fp) else 0.6
    v0 = 0.08*(64/m)**0.5
    nrec = nrec or int(25/dtf)
    act, vm, Tm, v, T = sim2(F, VG, TG, float(m), f, nu, N, int(tmax/dtf), nrec, T0, v0, 0, seed, dtf, 1/m)
    return (np.arange(len(act))+1)*nrec*dtf*m, act

if __name__ == '__main__' and len(sys.argv) > 1:
    job = sys.argv[1]
    m = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    nu = float(sys.argv[3]) if len(sys.argv) > 3 else 10.
    if job == 'edge':
        rs = [0.96, 1.0, 1.03, 1.05, 1.07, 1.09, 1.12, 1.2]
        out = np.array([edge_speed(m, nu, r) for r in rs])
        for r, o in zip(rs, out): print(m, nu, r, np.round(o, 2), flush=True)
        np.savez(f'{D}/stoch_edge_m{m}_nu{nu:g}.npz', r=rs, s=out)
    if job == 'bulk':
        rs = [float(x) for x in sys.argv[4:]]
        for r in rs:
            t, a = bulk(m, nu, r)
            np.savez(f'{D}/stoch_bulk_m{m}_nu{nu:g}_r{r:g}.npz', t=t, a=a)
            print(m, nu, r, np.round(a[::10], 3), flush=True)

def make_fig():
    import glob, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for m, nu, mk in ((64, 4, 'o-'), (64, 10, 's-'), (64, 60, '^-'), (256, 4, 'd--')):
        z = np.load(f'{D}/stoch_edge_m{m}_nu{nu:g}.npz'); s = np.clip(np.nanmedian(z['s'], 1), -3, 1)
        ax[0].plot(z['r'], s, mk, label=f'm={m}, nu={nu}')
    ax[0].axhline(0, c='k', lw=.5); ax[0].axvline(1, c='gray', ls=':'); ax[0].axvspan(0.975, 0.981, color='g', alpha=.2, label='measured zero')
    ax[0].set_xlabel('f/f_min'); ax[0].set_ylabel('s/c (clipped at -3)'); ax[0].legend(fontsize=8); ax[0].set_title('edge speed')
    rs = []; st = []
    for fn in sorted(glob.glob(f'{D}/stoch_bulk_m64_nu4_r*.npz')):
        r = float(fn.split('_r')[1][:-4]); z = np.load(fn); rs.append(r)
        ax[1].plot(z['t']/64, z['a'], label=f'{r:g}'); st.append(z['a'][-40:].mean())
    ax[1].set_xlabel('t/m'); ax[1].set_ylabel('active fraction'); ax[1].legend(title='f/f_min', fontsize=7); ax[1].set_title('m=64, nu=4: rho(t) from all-hot start')
    ax[1].set_xscale('log')
    ax[2].plot(rs, st, 'o-', label='m=64'); r2 = []; s2 = []
    for fn in sorted(glob.glob(f'{D}/stoch_bulk_m256_nu4_r*.npz')):
        z = np.load(fn); r2.append(float(fn.split('_r')[1][:-4])); s2.append(z['a'][-40:].mean())
    ax[2].plot(r2, s2, 's--', label='m=256'); ax[2].axhline(0.4, c='r', ls=':', label='measured 0.4 (0.989 f_min)')
    ax[2].set_xlabel('f/f_min'); ax[2].set_ylabel('rho at t = 3000m'); ax[2].legend(fontsize=8); ax[2].set_title('late-time activity (discontinuous)')
    plt.tight_layout(); plt.savefig(os.path.join(os.path.dirname(D), 'fig_front_stoch.png'), dpi=130)

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'fig':
    make_fig()
