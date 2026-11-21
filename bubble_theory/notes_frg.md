# Functional RG for the wall, and the threshold anisotropy

Continues [notes_phi4_interface.md](notes_phi4_interface.md). Builds on the TKD derivation in [`bubble_sim/tkd_derivation.md`](../bubble_sim/tkd_derivation.md), whose notation is used here: α = angle between the wall **normal** and **b** (edge α = 0, `--burgers-angle 90`; screw α = π/2, `--burgers-angle 0`), a = noise anisotropy, $A(\alpha)=1+a\cos2\alpha$.

All numbers in sections 3–8 are produced by `frg_anisotropy.py`; reports for the two aspect ratios are in `frg_anisotropy_A2.txt` and `frg_anisotropy_A4.txt`.

## 1. The interface action

Projection (notes_phi4_interface §4) gives, in force units ($c=\varepsilon_0=1$, $\mu=2/\eta\approx2.12$):

$$\mu^{-1}\partial_t u=(c/\mu)\,\partial_x^2u+h+F(x,u),\qquad
\overline{F(x,u)F(x',u')}=\Delta(x-x',u-u'),$$

$$\hat\Delta(k_t,u)=\int\frac{dk_n}{2\pi}\,e^{ik_nu}\,S(k_t\hat t+k_n\hat n)\,|\hat\rho(k_n)|^2,\qquad
S(\mathbf k)\propto|\mathbf k|\,[1+a\cos2(\chi_{\mathbf k}-\beta)]\,e^{-k^2w^2}.$$

MSR action: $S=\int\hat u\,[\mu^{-1}\partial_tu-(c/\mu)\partial_x^2u-h]-\tfrac12\int_{x,x',t,t'}\hat u\hat u'\,\Delta(x-x',u-u')$.

Two properties of the dislocation noise matter:

- **$\int du\,\hat\Delta(0,u)=S(0)=0$**: random-bond-like. At depinning RB and RF flow to the same fixed point, so the class is qEW for every orientation as long as the KPZ term is absent.
- **$\hat\Delta(k_t,u)$ is analytic in $k_t$** on the scale $1/w$. Only $\hat\Delta(0,u)\equiv\Delta_\alpha(u)$ survives at large scales; $k_t^2$ corrections carry two extra powers of the rescaling factor and are irrelevant at the qEW fixed point. At large scales $\Delta_\alpha(u)=A(\alpha)\,\Delta_0(u)$ (tkd_derivation §3, checked numerically there).

## 2. One loop

Rescale $x\to e^lx$, $u\to e^{\zeta l}u$, $t\to e^{zl}t$. In $d=1$, $\varepsilon=4-d=3$. With $\tilde\Delta$ the dimensionless correlator,

$$\partial_l\tilde\Delta(u)=(\varepsilon-2\zeta)\,\tilde\Delta(u)+\zeta\,u\,\tilde\Delta'(u)-\tfrac12\,\frac{d^2}{du^2}\big[\tilde\Delta(u)-\tilde\Delta(0)\big]^2 .$$

- The fixed point $\tilde\Delta^*$ has a cusp, $\tilde\Delta^{*\prime}(0^+)<0$, which is what gives a finite threshold.
- At depinning the non-potential flow generates $\int\Delta\neq0$, so RB-type disorder reaches the RF fixed point: $\zeta=\varepsilon/3=1$ at one loop, about 1.25 at two loops and in numerics. $z=2-\tfrac29\varepsilon$, $\nu=1/(2-\zeta)$.
- $\varepsilon=3$ is not small. Exponents come from simulations. What the FRG controls are scaling forms, symmetries and amplitude relations, which is what the rest of this note uses.

## 3. Exact amplitude scaling, and why it fails here

If $\Delta_\alpha(x,u)=A(\alpha)\,\delta(x)\,\Delta_0(u)$ at all scales that matter, rescaling $x\to\lambda x$ maps the disorder strength $A\to A\lambda^{-1}$ at fixed elastic constant. Equating the two fixes $\lambda=L_c\propto A^{-1/3}$, and every threshold scales with it:

$$h_c(L,\alpha)=A(\alpha)^{2/3}\,\mathcal H\big(L\,A(\alpha)^{1/3}\big),\qquad
h_c(\infty,\alpha)=A(\alpha)^{2/3}\,h_c(\infty,A=1).$$

This is exact within the elastic description (not a Larkin estimate), and it predicts large effects. Data, A2, $h_0$ = isotropic threshold at the same L:

| L | edge/h₀ (a = 1) | prediction 2^{2/3} | screw/h₀ (a = 1) | prediction | 45°/h₀ | prediction |
|---|---|---|---|---|---|---|
| 256 | 1.185 | 1.587 | 0.700 | 0 | 0.788 | 1 |
| 2048 | 1.058 | 1.587 | 0.911 | 0 | 0.958 | 1 |

It fails badly, and it fails worse as L grows. The assumption that breaks is the first one: the δ(x) form needs $L_c\gg w$.

## 4. The anisotropy builds up only beyond a few w

The force on a straight segment of length ℓ has variance $\int d^2k\,S(\mathbf k)|\hat\rho(k_n)|^2\,\mathrm{sinc}^2(k_t\ell/2)$. Its anisotropy, normalized to its large-ℓ value,

$$R(\ell)=\frac{1}{a}\,\frac{\mathrm{Var}_\ell(\text{edge})-\mathrm{Var}_\ell(\text{screw})}{\mathrm{Var}_\ell(\text{edge})+\mathrm{Var}_\ell(\text{screw})},$$

is:

| ℓ | 1 | 3 | 6 | 12 | 24 | 48 | 96 |
|---|---|---|---|---|---|---|---|
| R(ℓ) | −0.06 | −0.03 | 0.05 | 0.28 | 0.57 | 0.77 | 0.88 |

Below about 2w = 6 the force is isotropic. This is by construction: the noise is normalized to unit variance at every a, and a single correlation area does not know which way the wall points. The Larkin length is $L_c\approx6$ (tkd_derivation §4), so the threshold, which is set at the Larkin scale, sees $A_{\rm eff}=1+aR(L_c)\cos2\alpha$ with $R(L_c)\approx0.05$:

$$h_c(\infty,\alpha)\simeq h_c^{\rm iso}\,[1+0.05\,a\cos2\alpha]^{2/3}:\quad\text{edge } +3.4\%,\ \text{screw }-3.4\%\ \text{at }a=1.$$

Extrapolations of the data with $L^{-3/4}$ (qEW ν):

| | a = 0 | edge a = 1 | screw a = 1 | prediction (edge / screw) |
|---|---|---|---|---|
| A2 | 0.02544 | 0.02604 (+2.4%) | 0.02450 (−3.7%) | 0.0263 / 0.0246 |
| A4 | 0.02544 | 0.02613 (+2.7%) | 0.02532 (−0.5%) | 0.0263 / 0.0246 |

The infinite-system threshold is isotropic to within a few percent, as predicted. The extrapolation is the weak point. The flat-analysis report finds that fits with the qKPZ ν send all of them to about 0.0254, and the screw A2 fit has χ²/dof = 5.7. The data allow $R(L_c)$ between 0 and about 0.05, so the effective pinning scale is ℓ ≲ 6.

## 5. The anisotropy is a finite-size effect

The finite-size correction is set at scale L, where $R(L)\to1$, so it is fully anisotropic while $h_c(\infty)$ is not:

$$h_c(L,\alpha)\simeq h_c^{\rm iso}(\infty)+B(\alpha,a)\,L^{-1/\nu}.$$

**Exact symmetries.** Rotating **b** by 90° is $a\to-a$, and reflecting the box about the wall normal maps a 45° wall at $a$ onto one at $-a$. Hence:

$$h_{\rm edge}(a)=h_{\rm screw}(-a),\qquad h_{45}(a)=h_{45}(-a).$$

So the odd part $\tfrac12(h_{\rm edge}-h_{\rm screw})$ is odd in a, and the even part $\tfrac12(h_{\rm edge}+h_{\rm screw})-h_0$ and $h_{45}-h_0$ start at $a^2$.

Data, A2 (A4 similar):

| L | odd, a = 0.5 | odd, a = 1 | even, a = 0.5 | even, a = 1 | h₄₅ − h₀, a = 1 |
|---|---|---|---|---|---|
| 256 | +0.00276 | +0.00633 | −0.00031 | −0.00150 | −0.00553 |
| 512 | +0.00178 | +0.00441 | +0.00002 | −0.00089 | −0.00390 |
| 1024 | +0.00111 | +0.00288 | −0.00002 | −0.00068 | −0.00208 |
| 2048 | +0.00071 | +0.00188 | −0.00001 | −0.00039 | −0.00107 |

- The odd part dominates and is close to linear in a: the ratio of a = 0.5 to a = 1 is 0.38–0.44 (linear would give 0.5; the difference is the $a^3$ term).
- The even part at a = 0.5 is zero within errors, as $a^2$ scaling of the a = 1 value requires (−0.0004/4 = −0.0001 at L = 2048).
- Everything decays: the odd part like $L^{-0.6}$.

**One-parameter form.** Writing $h=h_0(L)\,[1+a\,C(L)\cos2\alpha]^{2/3}$, the value of C extracted from edge and screw, at a = 0.5 and a = 1, agrees to a few percent (A2 screw at a = 1 is 30–50% higher):

| L | 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|
| C (A2, mean of four) | 0.33 | 0.23 | 0.15 | 0.10 |

$C\propto L^{-0.60}$ (A2), $L^{-0.65}$ (A4). This is between $1/\nu_{\rm qEW}=0.75$ and $1/\nu_{\rm qKPZ}=0.58$, as the finite-size form requires.

**FSS amplitudes** ($B-B_{a=0}$, fit with $L^{-3/4}$): edge +0.151 → +0.285 and screw −0.185 → −0.447 for a = 0.5 → 1 (A2). They are approximately odd and linear in $a\cos2\alpha$.

## 6. What it is not: orientation smearing

The obvious alternative is that the walls are rough, so a wall of mean orientation α samples many local orientations. Pinning variances add along the wall, so that would give $C=\langle\cos2\theta\rangle$ over the tangent directions (θ measured from the mean wall line). The critical walls give, at chord length 6:

| | L = 256 | 512 | 1024 | 2048 |
|---|---|---|---|---|
| edge a = 1: ⟨cos 2θ⟩ | +0.69 | +0.66 | +0.66 | +0.65 |
| screw a = 1: −⟨cos 2θ⟩ | +0.11 | +0.35 | +0.51 | +0.56 |
| C from thresholds | 0.33 | 0.23 | 0.15 | 0.10 |

Edge walls are nearly straight and stay so: their local orientations would give $C\approx0.65$ at all L. Faceted screw walls go the other way. Smearing is excluded. The small anisotropy comes from the noise (section 4), not from the wall geometry.

## 7. Screw walls at a = 1

At exact screw orientation $\hat\Delta(0,u)=0$. The force is a pure derivative along the wall, $F=\partial_xG$, with $\hat\Delta\propto k_t^2$.

- **It still pins.** The force averaged over ℓ is $[G(\ell)-G(0)]/\ell\sim\ell^{-1}$, against an elastic force $\sim\ell^{-2}$. The Larkin argument goes through, with $L_c\propto\sigma r_f/\sqrt{\Delta_G}$ and $h_c\propto\Delta_G$, linear in the variance instead of the 2/3 power.
- **It generates the missing term.** For a tilted segment, $\int dx\,\partial_xG(x,u(x))=-\int dx\,u'\,\partial_uG$. A slope s produces a $k_t=0$ correlator $A(\pi/2+s)=1-a\cos2s\simeq(1-a)+2as^2$. A screw wall is pinned by its own roughness, and it lowers its energy by rotating segments towards edge-like orientations. This is consistent with the measured tangents: $-\langle\cos2\theta\rangle$ grows from 0.11 to 0.56 with L, which is the faceting of tkd_derivation §6.3.

## 8. KPZ coefficient and crossover length

The TKD coefficient is the tilt curvature of the threshold of a **segment** of size ℓ, not of the infinite system. With the one-parameter form,

$$\frac{h_c''}{h_c}\Big|_{\rm edge}=-\frac{8}{3}\,\frac{aC(\ell)}{1+aC(\ell)},\qquad
\frac{h_c''}{h_c}\Big|_{\rm screw}=+\frac{8}{3}\,\frac{aC(\ell)}{1-aC(\ell)}.$$

- **Signs.** These are fixed by the symmetry and the sign of the odd part: edge is qKPZ with $\lambda v>0$ (directed-percolation depinning, ζ = 0.63); screw has $\lambda v<0$ (faceting). This is the conclusion of tkd_derivation §5, now with a measured, scale-dependent magnitude.
- **Coupling.** $g(\ell)=|h_c''(\ell)|\,r_f(\ell/L_c)^\zeta/(c/\mu)$, using $C(\ell)=\min[R(\ell),\,C_{\rm fit}(\ell)]$. A segment's threshold cannot be more anisotropic than the force it feels; the power law is extrapolated below ℓ = 256. This gives:

| ℓ | 12 | 24 | 48 | 96 | 192 |
|---|---|---|---|---|---|
| g, edge a = 1 | 0.23 | 0.89 | 2.5 | 5.3 | 9.5 |

so $\ell_\times\approx26$ sites. The flat analysis sees the edge walls cross over to ζ ≈ 0.6 at about 100–150 sites. The estimate is short by a factor of 4–6. Since $g\propto r_f L_c^{-\zeta}$, this needs $L_c$ larger, or $r_f$ smaller, by about 2 relative to the Larkin estimates. Those are O(1) prefactors that are not controlled.
- **At L = ∞** the anisotropic part of $h_c''$ shrinks to the $R(L_c)$ level, but it stays finite. Edge walls remain qKPZ asymptotically, at a crossover scale that grows as the anisotropy fades.

## 9. Predictions to test

| test | prediction |
|---|---|
| 45° walls at a = 0.5 | $h_{45}-h_0$ is $\tfrac14$ of its a = 1 value (even in a): −0.0014 at L = 256, −0.0003 at L = 2048 |
| edge and screw at a = 0.25 | odd part about half the a = 0.5 value; even part zero |
| varying the noise width w at fixed Δ | the anisotropy of $h_c(\infty)$ grows with $R(L_c)$, i.e. as $L_c/w$ grows (smaller w or weaker Δ) |
| screw walls at a = 1, varying Δ | the pinning in the derivative channel scales as $h_c\propto\Delta^2$, edge as $\Delta^{4/3}$, but only before the self-generated $k_t=0$ term takes over |
| $C(L)$ at larger L | $C\propto L^{-1/\nu}$ down to a floor $C_\infty\approx R(L_c)\lesssim0.05$ |

## References

- T. Nattermann, S. Stepanow, L.-H. Tang, H. Leschhorn, J. Phys. II 2, 1483 (1992); O. Narayan, D. S. Fisher, PRB 48, 7030 (1993): one-loop FRG for depinning.
- P. Chauve, P. Le Doussal, K. J. Wiese, PRL 86, 1785 (2001); P. Le Doussal, K. J. Wiese, P. Chauve, PRB 66, 174201 (2002): two-loop FRG at depinning.
- P. Le Doussal, K. J. Wiese, PRE 67, 016121 (2003): FRG with the KPZ term (anisotropic depinning), *to check*.
- L.-H. Tang, M. Kardar, D. Dhar, PRL 74, 920 (1995): KPZ coefficient from anisotropic thresholds.
