"""Non-equilibrium Gaussian closure at fixed v (see notes_noneq_closure.md).

Iterate: B -> C_F(tau)=Delta_{B(tau)}(v tau) -> S_F -> S_du = S_F Im chi/(eta w) -> B.
All transforms use DCT-IV on a midpoint tau grid and a midpoint omega grid
(w_k = pi (k+1/2)/(P N dtau)), P = zero-padding factor (resolves the w^{-3/2} tail).
"""
import numpy as np
from scipy.fft import dct
from scipy.special import ndtr
import sc, linresp

SQ2PI = np.sqrt(2 * np.pi)


def _p3(x, sd):
    """E[(x + sd z)_+^3], z ~ N(0,1)."""
    u = x / sd
    return x * (x * x + 3 * sd * sd) * ndtr(u) + sd * (x * x + 2 * sd * sd) * np.exp(-0.5 * u * u) / SQ2PI


def Delta_s(x, s, A2=1.0):
    """(A^2/3) M4 smeared by Gaussian of variance s, x >= 0 array, s array (same shape)."""
    sd = np.sqrt(np.maximum(s, 1e-12))
    out = np.zeros_like(x)
    m = x < 2.0 + 8.0 * sd
    xx, ss = x[m], sd[m]
    tot = 0.0
    for j, c in enumerate([1, -4, 6, -4, 1]):
        tot = tot + c * _p3(xx + 2.0 - j, ss)
    out[m] = A2 / 3.0 * tot / 6.0
    return out


class NonEq:
    def __init__(self, m, A=1.0, eta=1.0, kappa=1.0, dtau=0.1, pad=4, nq=2048, tau_max=None, vmin=0.005):
        self.m, self.A, self.eta, self.kappa, self.dtau, self.pad = m, A, eta, kappa, dtau, pad
        if tau_max is None:
            tau_max = max(40.0 * m, 50.0 / vmin)
        self.N = int(tau_max / dtau)
        self.taus = (np.arange(self.N) + 0.5) * dtau
        self.R, _ = sc.response_and_b(self.taus, float(m), eta, kappa, nq)
        M = pad * self.N
        self.M = M
        self.dw = np.pi / (M * dtau)
        self.w = (np.arange(M) + 0.5) * self.dw
        self.imchi = np.imag(linresp.chi_loc(self.w, m, eta, kappa))
        self.Hdu = self.imchi / (eta * self.w)      # S_du = S_F * Hdu

    def _dct4(self, x):
        return 0.5 * dct(x, type=4)

    def S_F(self, C):
        c = np.zeros(self.M); c[:self.N] = C
        return 2.0 * self.dtau * self._dct4(c)

    def B_from_S(self, Sdu):
        c = self._dct4(Sdu)[:self.N]
        return 2.0 / np.pi * self.dw * (Sdu.sum() - c)

    def CF(self, B, v):
        return Delta_s(v * self.taus, B, self.A ** 2)

    def step(self, B, v):
        S = self.S_F(self.CF(B, v))
        Sdu = S * self.Hdu
        return self.B_from_S(Sdu), S, Sdu

    def P2_from(self, Sdu):
        return self.eta / np.pi * self.dw * np.sum(self.w ** 2 * Sdu)

    def P_second_order(self, v):
        """P2 with B=0 (C_F = Delta(v tau)): must equal linresp.power."""
        B0 = np.zeros(self.N)
        _, S, Sdu = self.step(B0, v)
        return self.P2_from(Sdu)

    def solve(self, v, B0=None, alpha=0.5, tol=1e-6, maxit=300, verbose=False):
        B = np.zeros(self.N) if B0 is None else B0.copy()
        sel = None
        for it in range(maxit):
            Bn, S, Sdu = self.step(B, v)
            # relative change measured where B matters (v*tau < ~4 + sqrt(B))
            den = np.maximum(np.abs(Bn), 1e-3)
            msk = v * self.taus < 40.0
            err = np.max(np.abs(Bn - B)[msk] / den[msk])
            B = (1 - alpha) * B + alpha * Bn
            if verbose:
                print(it, err)
            if err < tol:
                break
        Bn, S, Sdu = self.step(B, v)
        Fp = sc.fpin(self.dtau, self.R, B, v, 1.0, self.A ** 2, sc._PIECES, 0)
        P1 = v * Fp
        P2 = self.P2_from(Sdu)
        return dict(B=B, Fp=Fp, P1=P1, P2=P2, it=it + 1, err=err)

    def flow(self, vs, **kw):
        res = [self.solve(v, **kw) for v in vs]
        g = lambda k: np.array([r[k] for r in res])
        return dict(f1=self.eta * vs + g('Fp'), Fp=g('Fp'), P1=g('P1'), P2=g('P2'),
                    f2=self.eta * vs + g('P2') / vs, it=g('it'), err=g('err'))


def main():
    import glob, sys, time
    vs = np.geomspace(0.005, 0.3, 30)
    ms = [int(a) for a in sys.argv[1:]] or [256, 1024]
    out = {'vs': vs}
    for m in ms:
        kw = dict(dtau=0.2, pad=2) if m >= 4096 else {}
        t = time.time()
        L = NonEq(m, **kw)
        fl = L.flow(vs)
        for k, a in fl.items():
            out[f'{k}_m{m}'] = a
        for k in ('f1', 'f2'):
            i = np.argmin(fl[k])
            print(m, k, 'fmin', fl[k][i], 'vmin', vs[i], flush=True)
        print(m, 'max err', fl['err'].max(), 'its', fl['it'].max(), 'P2/P1', (fl['P2'] / fl['P1'])[::5], time.time() - t, flush=True)
        for fn in sorted(glob.glob(f'data/closure_check_m{m}_f*.npz')):
            d = np.load(fn)
            r = L.solve(float(d['v']))
            Bn = np.interp(d['taus'], L.taus, r['B'])
            tag = f"m{m}"
            out[f'cmp_taus_{tag}'] = d['taus']; out[f'cmp_Bmeas_{tag}'] = d['B']
            out[f'cmp_Bnoneq_{tag}'] = Bn; out[f'cmp_v_{tag}'] = float(d['v']); out[f'cmp_Fp_{tag}'] = r['Fp']
            out[f'cmp_Fmeas_{tag}'] = float(d['F_meas']); out[f'cmp_P2P1_{tag}'] = r['P2'] / r['P1']
            out[f'cmp_bcl_T_{tag}'] = d['bcl'] * float(d['Tcl'])
            print(m, 'v', float(d['v']), 'Fp', r['Fp'], 'meas', float(d['F_meas']), 'P2/P1', r['P2'] / r['P1'], 'ratio', Bn / d['B'], flush=True)
    np.savez(f'data/sc_noneq_part_{"_".join(map(str, ms))}.npz', **out)


if __name__ == '__main__':
    main()
