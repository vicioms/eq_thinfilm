"""Hybrid closure: driven-spectrum SHAPE of B(tau), amplitude fixed by exact power balance.
See notes_hybrid_closure.md.  Fixed point at fixed v of
   S_F = FT[Delta_B(v tau)],  S0 = S_F Im chi/(eta w),  mu = v F_pin(B)/P2[S0],
   B = (2/pi) int dw mu S0 (1 - cos w tau),
with P2[S] = eta int dw/pi w^2 S.  Then eta<du'^2> = v F_pin exactly.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np, sys, time
import sc, sc_noneq


class Hybrid(sc_noneq.NonEq):
    def solve(self, v, B0=None, alpha=0.5, tol=1e-6, maxit=400):
        B = np.zeros(self.N) if B0 is None else B0.copy()
        msk = v * self.taus < 40.0
        for it in range(maxit):
            _, S, S0 = self.step(B, v)
            Fp = sc.fpin(self.dtau, self.R, B, v, 1.0, self.A ** 2, sc._PIECES, 0)
            P2 = self.P2_from(S0)
            mu = v * Fp / P2
            Bn = self.B_from_S(mu * S0)
            err = np.max(np.abs(Bn - B)[msk] / np.maximum(np.abs(Bn), 1e-3)[msk])
            B = (1 - alpha) * B + alpha * Bn
            if err < tol:
                break
        _, S, S0 = self.step(B, v)
        Fp = sc.fpin(self.dtau, self.R, B, v, 1.0, self.A ** 2, sc._PIECES, 0)
        mu = v * Fp / self.P2_from(S0)
        Bf = self.B_from_S(mu * S0)    # consistent B
        Fp2 = sc.fpin(self.dtau, self.R, Bf, v, 1.0, self.A ** 2, sc._PIECES, 0)
        P2 = self.P2_from(mu * S0)
        return dict(B=B, Fp=Fp, P1=v * Fp, P2=P2, mu=mu, it=it + 1, err=err, Bmax=B.max())

    def flow(self, vs):
        """continuation from the largest v downward"""
        res = {}
        B0 = None
        for v in vs[::-1]:
            r = self.solve(v, B0=B0)
            B0 = r['B']
            res[v] = r
        out = [res[v] for v in vs]
        g = lambda k: np.array([r[k] for r in out])
        return dict(Fp=g('Fp'), f=self.eta * vs + g('Fp'), P1=g('P1'), P2=g('P2'), mu=g('mu'),
                    it=g('it'), err=g('err'), Bmax=g('Bmax'), Bs=[r['B'] for r in out])


def main():
    import glob
    ms = [int(a) for a in sys.argv[1:]] or [256, 1024]
    vs = np.geomspace(0.005, 0.3, 25)
    out = {'vs': vs}
    for m in ms:
        t = time.time()
        kw = dict(dtau=0.2, pad=2) if m >= 4096 else {}
        H = Hybrid(m, **kw)
        fl = H.flow(vs)
        for k in ('Fp', 'f', 'P1', 'P2', 'mu', 'it', 'err', 'Bmax'):
            out[f'{k}_m{m}'] = fl[k]
        i = np.argmin(fl['f'])
        # parabolic refinement in log v
        lv = np.log(vs)
        if 0 < i < len(vs) - 1:
            c = np.polyfit(lv[i-1:i+2], fl['f'][i-1:i+2], 2)
            x0 = -c[1] / (2 * c[0]); fm = np.polyval(c, x0); vm = np.exp(x0)
        else:
            fm, vm = fl['f'][i], vs[i]
        out[f'fmin_m{m}'] = fm; out[f'vmin_m{m}'] = vm
        L = sc.Line(m)
        fe, Te, Fe = L.flow(vs)
        out[f'f_eq_m{m}'] = fe
        ie = np.argmin(fe)
        print(m, 'hybrid fmin', fm, 'vmin', vm, 'grid', fl['f'][i], vs[i], '| eq fmin', fe[ie], vs[ie],
              'maxerr', fl['err'].max(), 'maxit', fl['it'].max(), 'time', time.time() - t, flush=True)
        print(' f_hyb', np.round(fl['f'][::3], 4), '\n f_eq ', np.round(fe[::3], 4),
              '\n mu', np.round(fl['mu'][::3], 3), '\n Bmax', np.round(fl['Bmax'][::3], 1), flush=True)
        # measured flow
        d = np.load(f'data/flow_m{m}.npz')
        vm_, fm_ = d['v'].mean(1), d['f'].mean(1)
        o = np.argsort(vm_); vm_, fm_ = vm_[o], fm_[o]
        sel = vm_ > vs[0]
        fh = np.interp(vm_, vs, fl['f']); fq = np.interp(vm_, vs, fe)
        out[f'meas_v_m{m}'] = vm_; out[f'meas_f_m{m}'] = fm_
        out[f'hyb_at_meas_m{m}'] = fh; out[f'eq_at_meas_m{m}'] = fq
        for a, b, c in zip(vm_, fm_, zip(fh, fq)):
            print(f'  v={a:.4f} meas f={b:.4f} hyb={c[0]:.4f} eq={c[1]:.4f}')
        # B(tau) comparison
        for fn in sorted(glob.glob(f'data/closure_check_m{m}_f*.npz')):
            c = np.load(fn); v = float(c['v'])
            r = H.solve(v)
            Bh = np.interp(c['taus'], H.taus, r['B'])
            out[f'cmp_taus_m{m}'] = c['taus']; out[f'cmp_Bmeas_m{m}'] = c['B']; out[f'cmp_Bhyb_m{m}'] = Bh
            out[f'cmp_Beq_m{m}'] = c['bcl'] * float(c['Tcl']); out[f'cmp_v_m{m}'] = v
            out[f'cmp_Fhyb_m{m}'] = r['Fp']; out[f'cmp_Fmeas_m{m}'] = float(c['F_meas'])
            vt = v * c['taus']
            print(m, 'v', v, 'Fp hyb', r['Fp'], 'meas', float(c['F_meas']), 'mu', r['mu'])
            for lo, hi in [(0, 1), (1, 3), (3, 10), (10, 30), (30, 100), (100, 1e9)]:
                s = (vt >= lo) & (vt < hi)
                if s.any():
                    print(f'   v tau [{lo},{hi}) hyb/meas {np.mean(Bh[s]/c["B"][s]):.2f}  eq/meas {np.mean(out[f"cmp_Beq_m{m}"][s]/c["B"][s]):.2f}')
    np.savez('data/sc_hybrid.npz', **out)


if __name__ == '__main__':
    main()
