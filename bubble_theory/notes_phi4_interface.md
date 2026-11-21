# From the disordered φ⁴ field theory to the bubble interface

Model (as in `bubble_sim/bubble_depinning.py`, unit mobility, T = 0, φ = +1 inside the bubble):

$$\partial_t\phi = c\nabla^2\phi + \varepsilon_0[(1+r)\phi - \phi^3] + h + g .$$

Every number quoted as "checked" is verified by `check_projection.py` (defaults ε₀ = c = 1).

## 0. Energy

$\partial_t\phi = -\delta E/\delta\phi$ with

$$E[\phi]=\int d^2x\Big[\tfrac c2(\nabla\phi)^2+\tfrac{\varepsilon_0}{4}(1-\phi^2)^2-\tfrac{\varepsilon_0}{2}r\phi^2-(h+g)\phi\Big].$$

Disorder: random bond $\langle r(x)r(x')\rangle=\Delta_r\,\delta(x-x')$, random field $\langle g(x)g(x')\rangle=G(x-x')$, with $G=\Delta_g\delta$ for white noise. For the dislocation noise $G$ is the anisotropic kernel of section 4.3. (For $r$ uniform on $[-\Delta,\Delta]$, $\Delta_r=\Delta^2/3$.)

## 1. MSR action and disorder average

$$Z=\int\mathcal D\phi\,\mathcal D\hat\phi\;e^{-S},\qquad
S_0=\int d^2x\,dt\;\hat\phi\,\big[\partial_t\phi-c\nabla^2\phi-\varepsilon_0(\phi-\phi^3)-h\big].$$

The dynamics is deterministic, so $Z=1$ for every disorder realization and $\overline{Z}$ can be taken directly: no replicas. Gaussian averaging of the linear couplings $-\hat\phi(g+\varepsilon_0 r\phi)$ gives

$$S_{\rm dis}=-\frac12\int dt\,dt'\!\int d^2x\,d^2x'\;\hat\phi(x,t)\hat\phi(x',t')\Big[G(x-x')+\varepsilon_0^2\,\Delta_r\delta(x-x')\,\phi(x,t)\phi(x',t')\Big].$$

The term is non-local in time ($t,t'$ unrestricted): this is how quenched disorder appears in a dynamical action, and it is the only place disorder enters.

## 2. Saddle points and the fluctuation spectrum

The causal saddle point has $\hat\phi=0$ and $\phi$ solving the clean equation.

**Kink.** $\phi_0(s)=\tanh(s/w)$, $w^2=2c/\varepsilon_0$. Constants (checked):

| quantity | value |
|---|---|
| $N_1=\int\phi_0'^2$ | $4/(3w)$ |
| $\int\phi_0'$ | 2 |
| $\sigma=cN_1$ = kink energy | $(2\sqrt2/3)\sqrt{c\varepsilon_0}$ |
| $\eta$ (unit mobility) | $N_1$ |

**Fluctuation operator.** $\mathcal L=-c\partial_s^2+\varepsilon_0(3\phi_0^2-1)=\tfrac{\varepsilon_0}{2}\big[-\partial_z^2+4-6\,\mathrm{sech}^2z\big]$, $z=s/w$ (Pöschl–Teller, ℓ = 2). Spectrum (checked):

| mode | eigenvalue | meaning |
|---|---|---|
| $\phi_0'$ | 0 | translation → the interface field $u$ |
| $\mathrm{sech}\,z\tanh z$ | $3\varepsilon_0/2$ | shape mode: wall breathing (width, amplitude) |
| continuum | $\ge 2\varepsilon_0$ | bulk magnons |

Everything except the zero mode is gapped by at least $3\varepsilon_0/2$. That gap is what justifies an interface theory on times $\gg 1/\varepsilon_0$ and lengths $\gg w$.

**Bubble.** $\phi_B=\phi_0(R-|x-X|)$. In the clean case it is Langer's critical nucleus at $R^*=\sigma/2h$: the $m=0$ mode is the single unstable direction (→ collective coordinate $R$), $m=\pm1$ are the two translation zero modes (→ $X$), and $m\ge2$ has stiffness $\sigma(m^2-1)/R^2$.

## 3. Collective coordinates (Faddeev–Popov)

Parametrize the wall as the curve $\Gamma(t):\ X+(R+u(\theta,t))\,\hat e_\rho(\theta)$ and split both fields:

$$\phi=\phi_0\big(R+u(\theta,t)-\rho\big)+\chi,\qquad
\hat\phi=\hat u(\theta,t)\,\phi_0'\big(R+u-\rho\big)+\hat\chi,$$

with the gauge condition $\int d\rho\,\phi_0'\,\chi=0$ (and the same for $\hat\chi$) at every θ: $\chi$ lives entirely in the gapped sector. Inserting $1=\int\mathcal Du\,\prod_\theta\delta\big(\int\phi_0'(\phi-\Phi_u)\big)\det[\int\phi_0'\partial_u\Phi_u]$, the Faddeev–Popov determinant is $\prod_\theta N_1$ up to $O(w/R)$ corrections, i.e. a constant at this order.

The collective coordinates are the low Fourier components of the same fields:

$$\hat u(\theta)=\hat R+\hat{\mathbf X}\cdot\hat e(\theta)+\sum_{|m|\ge2}\hat u_m e^{im\theta},$$

and similarly for $u$, with $R$, $\mathbf X$ removed from $u$.

## 4. Tree-level projection

### 4.1 Deterministic part

With $\zeta=R+u-\rho$,

$$\nabla^2\phi=\phi_0''\Big(1+\frac{u_\theta^2}{\rho^2}\Big)-\frac{\phi_0'}{\rho}+\frac{\phi_0'\,u_{\theta\theta}}{\rho^2}.$$

Insert into $S_0$, use $c\phi_0''=-\varepsilon_0(\phi_0-\phi_0^3)$, and project with $\int\rho\,d\rho\,\phi_0'(\cdot)\simeq R\int d\zeta\,\phi_0'(\cdot)$ (using $\int\phi_0'\phi_0''=0$):

$$S_0\to\int d\theta\,dt\;R\,\hat u\Big[\eta(\dot R+\partial_t u)-\sigma\partial_x^2u+\frac{\sigma}{R+u}-2h\Big]+O(\hat u\,u_\theta^2,\chi),\qquad \partial_x=\partial_\theta/R.$$

- The $m=0$ component gives $\eta\dot R=2h-\sigma/R$. Checked: the projected residual of the clean equation on $\tanh((R-\rho)/w)$ equals $-\sigma/R$ to six digits for R = 5…40.
- Expanding $\sigma/(R+u)$ gives the destabilizing $-\sigma u/R^2$.
- The $\tfrac12\sigma(\partial_xu)^2$ and $u\,\partial_x^2u$ terms of the second-order equation come from using the true normal distance instead of the radial one.

### 4.2 Disorder part

Each $\hat\phi$ in $S_{\rm dis}$ becomes $\hat u\,\phi_0'$, and the integral over the normal coordinate turns the bulk correlator into a function of $u(t)-u(t')$:

$$S_{\rm dis}\to-\frac12\int d\ell\,d\ell'\,dt\,dt'\;\hat u(\ell,t)\hat u(\ell',t')\;\Delta\big(\ell-\ell';\,u(\ell,t)-u(\ell',t')\big).$$

This is the starting action of the functional RG for depinning.

**Random field** (white $g$): $\Delta(\ell,v)=\delta(\ell)\,\Delta_{\rm RF}(v)$,

$$\Delta_{\rm RF}(v)=\Delta_g\!\int\!dy\,\phi_0'(y)\phi_0'(y-v)=\frac{4\Delta_g}{w}\,\frac{b\coth b-1}{\sinh^2 b},\qquad b=v/w$$

(checked), from $\int\mathrm{sech}^2x\,\mathrm{sech}^2(x-b)\,dx=4(b\coth b-1)/\sinh^2b$. Short-ranged in $v$ but $\int\Delta\,dv=4\Delta_g\neq0$: random-field class.

**Random bond.** $\hat\phi\phi\to\hat u\,\phi_0\phi_0'=-\tfrac12\hat u\,(\mathrm{sech}^2)'$, so the force correlator is $\Delta_{\rm RB}=-R''$ with

$$R(v)=\varepsilon_0^2\Delta_r\,w\,\frac{b\coth b-1}{\sinh^2b}$$

(checked; $\int\Delta_{\rm RB}\,dv=0$). With $\Delta_r=\Delta^2/3$ this is the Caballero et al. prefactor. RF and RB are different classes in statics, but flow to the same fixed point at depinning.

### 4.3 Dislocation disorder

For an anisotropic kernel $G$ the same projection gives, in Fourier space,

$$\Delta_\theta(q_\ell,v)=\int\frac{dq_u}{2\pi}\,e^{iq_uv}\,|\hat\phi_0'(q_u)|^2\,\hat G(q_\ell\hat t+q_u\hat n),\qquad
\hat\phi_0'(q)=\frac{\pi qw}{\sinh(\pi qw/2)}$$

(Fourier transform checked). Two things are new compared with white disorder:

1. It depends on the wall angle θ through $\hat t$, $\hat n$, and vanishes at screw orientation to leading order.
2. It can be long-ranged along the wall, through $q_\ell$.

On the bubble θ = θ(ℓ) runs around the loop, so the disorder strength is modulated along the perimeter. This is the source of $f_c(\theta)$ and of the disorder-generated KPZ term $-\tfrac12 f_c''(\theta)(\partial_xu)^2$.

## 5. Collective pinning of the whole bubble

Keeping only $m=0$, the $\hat R\hat R$ term is

$$-\frac12\int dt\,dt'\,\hat R(t)\hat R(t')\,(2\pi R)\,\bar\Delta\big(R(t)-R(t')\big),$$

with $\bar\Delta$ the perimeter average of $\Delta_\theta$ (per unit length, $\hat R$ normalized as $\int d\ell\,\hat u$). A rigid bubble is a single particle in a random force of variance $\propto\bar\Delta/(2\pi R)$ per unit length.

- **$2\pi R<L_c=(\sigma w)^{2/3}\bar\Delta^{-1/3}$ (Larkin length):** the bubble moves as one object, $h_d\simeq\sigma/2R+O\big((\bar\Delta/Rw)^{1/2}\big)$; the pinning part decays like $R^{-1/2}$.
- **$2\pi R>L_c$:** the $m\ge2$ modes relax and the elastic regime applies, $h_d\simeq f_c/2+\sigma/2R-A'(2\pi R)^{-1/\nu}$.

So the finite-size-scaling form has a lower bound of validity, $R\gtrsim L_c/2\pi$.

## 6. What integrating out χ adds

Corrections are in two small parameters, $w/R$ and $\Delta/\varepsilon_0^2$.

- **σ and η renormalize** at order $\Delta/\varepsilon_0^2$, through exchange of the shape mode ($3\varepsilon_0/2$).
- **Width fluctuations.** The shape mode is the breathing of the wall, so the scatter of the fitted width and amplitude in `bubble_sim/rayfit.py` is $\langle\chi^2\rangle$ projected on that mode. Its variance should scale as $\Delta/\varepsilon_0^2$: a direct check on the data.
- **Velocity dependence.** $\phi_0\to\phi_v(s)$ gives $O(v^2)$ corrections to η, irrelevant near threshold.
- **Lattice.** The discrete Laplacian makes σ(θ) anisotropic (consistent with the measured width peak at 45°). Peierls–Nabarro pinning scales as $|\hat\phi_0'(2\pi)|^2\approx2\times10^{-9}$ for $w=\sqrt2$ (checked): negligible.
- **Breakdown.** The local projection fails where the curvature radius approaches $w$: near overhangs, and for bubbles with $w/R$ not small. There χ must be kept.

## 7. Next step

The interface action is now derived rather than assumed: σ, η, 2h from the zero-mode projection; Δ(u) from projecting the disorder correlator; R and X as the m = 0, 1 components of the same field. The next step is the functional RG on $u$ with an angle-dependent $\Delta_\theta$: Larkin scaling $f_c(\theta)\propto\Delta_\theta^{2/3}(\sigma w)^{-1/3}$, the sign of $f_c''$ at screw, and the finite-size form of the anisotropy.
