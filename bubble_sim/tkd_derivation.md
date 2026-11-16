# Tang–Kardar–Dhar nel nostro modello: derivazione da zero

Obiettivo: partire dall'equazione di Ginzburg–Landau che simuliamo
(`bubble_depinning.py`, accoppiamento `field`, rumore di dislocazione) e
arrivare, senza parametri liberi nel ragionamento, a:

1. l'equazione di una parete con orientazione qualsiasi;
2. il disordine che quella parete sente, in funzione dell'angolo tra normale e **b**;
3. la soglia $h_c(\alpha)$ e il termine KPZ che l'anisotropia lascia sopravvivere
   alla soglia (argomento di Tang, Kardar e Dhar, PRL 74, 920 (1995));
4. il segno di quel termine per pareti edge e screw, e quindi la classe attesa;
5. il confronto con quello che abbiamo già nei dati.

Notazione: $\alpha$ = angolo tra la **normale** della parete e **b**. Parete edge:
$\alpha = 0$ (normale ∥ b, linea ⊥ b, `--burgers-angle 90` nella geometria flat).
Parete screw: $\alpha = \pi/2$ (linea ∥ b, `--burgers-angle 0`).

---

## 1. Il modello

Su reticolo quadrato, con $J = 1$, `coupling = field` ($r = 0$):

$$
\partial_t \phi = c\,\nabla^2\phi + \varepsilon_0(\phi - \phi^3) + h + g(\mathbf r),
\qquad g = \Delta\, f(\mathbf r),
$$

con $c = J = 1$, $\varepsilon_0 = 1$, $\Delta = 0.04$ e $f$ gaussiano a varianza 1 con spettro

$$
S(\mathbf k) \;\propto\; \frac{(1+a)(\mathbf k\cdot\hat{\mathbf b})^2 + (1-a)(\mathbf k\cdot\hat{\mathbf b}_\perp)^2}{|\mathbf k|}\,e^{-k^2w^2}
\;=\; |\mathbf k|\,\big[1 + a\cos 2(\chi_{\mathbf k}-\beta)\big]\,e^{-k^2w^2},
$$

dove $\chi_{\mathbf k}$ è l'angolo di $\mathbf k$ e $\beta$ quello di **b**. Normalizzato a
varianza 1 (vedi `dislocation_spectrum`), il prefattore è $A = 8\sqrt\pi\,w^3$:
$\int\!\frac{d^2k}{(2\pi)^2}\,A|\mathbf k|e^{-k^2w^2} = 1$.

Il campo $h$ favorisce $\phi = +1$ ed è **scalare**: agisce sulla parete come
una pressione, uguale per ogni orientazione. Questo è il punto da tenere a mente
per KPZ.

## 2. Dalla GL all'equazione della parete

**Profilo.** La parete piatta stazionaria risolve $c\Phi'' + \varepsilon_0(\Phi-\Phi^3) = 0$:

$$
\Phi(\zeta) = -\tanh\frac{\zeta}{\sqrt2\,\delta},\qquad \delta = \sqrt{c/\varepsilon_0} = 1,
$$

dove $\zeta$ è la distanza lungo la normale, verso la fase $-1$. È il `bubble()`
del codice, con $\phi = +1$ sotto la parete.

**Coordinata collettiva.** Si pone $\phi(x,y,t) = \Phi(y - u(x,t))$ per pendenze
piccole. Allora

$$
\partial_t\phi = -\Phi'\,\partial_t u,\qquad
\nabla^2\phi = \Phi''\,(1+u_x^2) - \Phi'\,u_{xx}.
$$

Si proietta sul modo di traslazione: si moltiplica per $\Phi'$ e si integra in
$\zeta$. Valgono $\int\Phi'\Phi'' = 0$ e $\int\Phi'[c\Phi''+\varepsilon_0(\Phi-\Phi^3)] = 0$,
e inoltre $\int\Phi' = \Phi(\infty)-\Phi(-\infty) = -2$. Si ottiene

$$
-\eta\,\partial_t u = -c\,\eta\,u_{xx} - 2h - \!\int\!\Phi'(\zeta)\,g(x,u+\zeta)\,d\zeta,
\qquad \eta = \int\Phi'^2 d\zeta = \frac{2\sqrt2}{3\delta},
$$

cioè

$$
\boxed{\;\partial_t u = c\,\partial_x^2 u + \mu\,\big[h + F(x,u)\big]\;},\qquad
\mu = \frac{2}{\eta} = \frac{3\delta}{\sqrt2}\approx 2.12,
$$

$$
F(x,u) = \int\rho(\zeta)\,g(x,u+\zeta)\,d\zeta,\qquad \rho = \tfrac12|\Phi'|,\quad \int\rho = 1 .
$$

È una qEW. $F$ è il campo casuale mediato sullo spessore della parete
($\sqrt2\,\delta \approx 1.4 < w = 3$, quindi è un effetto piccolo).

**Forma geometrica.** Senza l'ipotesi di pendenza piccola, la stessa proiezione
fatta lungo la normale locale dà, per qualunque orientazione,

$$
v_n = c\,\kappa + \mu\,\big[h + F(\mathbf r)\big],
$$

dove $v_n$ è la velocità normale e $\kappa$ la curvatura. Né $c\kappa$ né $\mu h$
dipendono dall'orientazione: tutta l'anisotropia sta in $F$.

## 3. Il disordine visto da una parete di orientazione α

Nel riferimento della parete ($x$ lungo la tangente $\hat{\mathbf t}$, $u$ lungo la
normale $\hat{\mathbf n}$) il correlatore della forza è

$$
\langle F(x,u)F(0,0)\rangle = \int\!\frac{d^2k}{(2\pi)^2}\,S(\mathbf k)\,|\hat\rho(k_n)|^2\,
e^{i(k_t x + k_n u)} .
$$

Su scale più lunghe della correlazione lungo la parete ($\sim w$) conta solo
l'integrale lungo $x$, che seleziona $k_t = 0$, cioè $\mathbf k \parallel \hat{\mathbf n}$:

$$
\Delta_\alpha(u) \equiv \int dx\,\langle F(x,u)F(0,0)\rangle
= \int\!\frac{dk_n}{2\pi}\,S(k_n\hat{\mathbf n})\,|\hat\rho(k_n)|^2 e^{ik_n u} .
$$

Per $\mathbf k\parallel\hat{\mathbf n}$ l'angolo tra $\mathbf k$ e **b** è proprio $\alpha$, quindi

$$
\boxed{\;\Delta_\alpha(u) = (1 + a\cos2\alpha)\,\Delta_0(u)\;},\qquad
D_0 \equiv \Delta_0(0) \simeq \frac{4w}{\sqrt\pi}\,\Delta^2 .
$$

La forma in $u$ è la stessa per ogni $\alpha$: cambia solo l'ampiezza.

* edge ($\alpha = 0$): fattore $1 + a$, il disordine più forte;
* screw ($\alpha = \pi/2$): fattore $1 - a$. **A $a = 1$ si annulla**: una parete
  screw è ancorata solo dalla parte a $k_t \neq 0$, cioè da correlazioni finite
  lungo la parete.

Inoltre $\int du\,\Delta_\alpha(u) = S(\mathbf k = 0) = 0$, perché $S \propto |\mathbf k|$.
Il correlatore ha integrale nullo, come per il disordine random-bond. Per la
depinning random-bond e random-field portano allo stesso punto fisso qEW
(Nattermann et al. 1992; Narayan & Fisher 1993), quindi a livello di
proiezione la classe è qEW per **ogni** orientazione. Quello che cambia con
$\alpha$ è solo l'ampiezza, e da qui parte TKD.

**Verifica numerica** sullo spettro del codice ($L = 1024$, $w = 3$; somma di
$C(x, u=0)$ lungo la parete):

| a | edge (α = 0) | α = 45° | screw (α = 90°) | previsione $(1+a\cos2\alpha)\,4w/\sqrt\pi$ |
|---|---|---|---|---|
| 0 | 6.770 | 6.770 | 6.770 | 6.770 |
| 0.5 | 10.155 | 6.770 | 3.385 | 10.155 / 6.770 / 3.385 |
| 1 | 13.540 | 6.770 | 0.000 | 13.541 / 6.770 / 0 |

Inoltre $\sum_u\Delta(u) = 0$ alla precisione di macchina, e $\Delta(u)$ cambia
segno per la prima volta a $u = 2w = 6$.

## 4. La soglia dipende dall'orientazione: $h_c(\alpha)$

**Simmetria.** $S$ è invariante per $\mathbf k\to-\mathbf k$ e per riflessione rispetto a
**b**, quindi $h_c(\alpha) = h_c(-\alpha) = h_c(\pi-\alpha)$. Edge e screw sono
estremi: $h_c'(0) = h_c'(\pi/2) = 0$.

**Massimo o minimo? Stima di Larkin.** Si riscrive l'equazione in unità di
campo, $\mu^{-1}\partial_t u = (c/\mu)\,\partial_x^2u + h + F$. Si confronta un
tratto di lunghezza $\ell$ spostato di $r_f$ (scala di correlazione in $u$):

* forza elastica $\sim (c/\mu)\,r_f/\ell^2$;
* forza casuale media $\sim (D/\ell)^{1/2}$.

Ne segue

$$
L_c = \Big[\frac{(c/\mu)^2 r_f^2}{D}\Big]^{1/3},\qquad
h_c \sim \frac{(c/\mu)\,r_f}{L_c^2} \;\propto\; D^{2/3},
\qquad\Rightarrow\qquad
\boxed{\;h_c(\alpha)\propto(1+a\cos2\alpha)^{2/3}\;}
$$

Numeri: $c/\mu = 0.47$, $r_f\approx w = 3$, $D_0 = 6.77\,\Delta^2 = 0.0108$ danno
$L_c \approx 6$ e $h_c \approx 0.04$. Misuriamo $0.025$–$0.029$: stesso ordine di
grandezza, con prefattori $O(1)$ non controllati.

Quindi **edge = massimo** di $h_c$ (direzione dura) e **screw = minimo**
(direzione facile). Le curvature sono

$$
\frac{h_c''}{h_c}\Big|_{\text{edge}} = -\frac{8a}{3(1+a)}\;\;(=-\tfrac43\text{ a }a=1),\qquad
\frac{h_c''}{h_c}\Big|_{\text{screw}} = +\frac{8a}{3(1-a)} .
$$

**Limiti di questa stima.**
* $L_c \approx 2w$: il regime di Larkin non è ben separato dalla scala del rumore.
* Su scala $L_c$ le pendenze sono già $O(r_f/L_c)\sim 0.5$, quindi la parete
  campiona orientazioni sparse di decine di gradi. Questo media l'anisotropia, e
  l'ampiezza $(1\pm a)^{2/3}$ è sovrastimata (a $a = 1$ la screw non ha
  $h_c = 0$).
* Il **segno** invece è robusto: $h_c$ cresce in modo monotono con l'ampiezza del
  disordine proiettato, e questa è massima per $\alpha = 0$ e minima per
  $\alpha = \pi/2$.

$h_c''$ va comunque misurato (§7).

## 5. L'argomento di Tang–Kardar–Dhar

**Velocità di una parete inclinata.** Si prende una parete con orientazione media
$\alpha_0$ e pendenza locale $s = \partial_x u$; la sua normale è ruotata di
$\arctan s$. Su scale più grandi della lunghezza di correlazione, ogni tratto si
muove come una parete piatta della sua orientazione, con velocità normale
$V(h,\alpha)$. La velocità lungo $y$ è

$$
v_y(h,s) = \sqrt{1+s^2}\;V\big(h,\alpha_0+\arctan s\big)
= V_0 + V_0'\,s + \tfrac12\big(V_0 + V_0''\big)s^2 + O(s^3),
$$

dove $'$ indica $\partial_\alpha$. Il coefficiente di $s^2$ è il $\lambda$ di KPZ:

$$
\partial_t u = c\,\partial_x^2u + \frac{\lambda}{2}(\partial_x u)^2 + \dots,\qquad
\lambda = V + \partial_\alpha^2 V .
$$

**Caso isotropo** ($a = 0$). $V$ non dipende da $\alpha$, quindi $\lambda = V = v$.
È il termine KPZ puramente geometrico, e **si annulla alla soglia**. Resta qEW:
è quello che misuriamo ($\zeta = 1.24\pm0.04$, $\nu = 1.31\pm0.09$).

**Caso anisotropo.** Vicino alla soglia $V(h,\alpha) = \mathcal V\big(h-h_c(\alpha)\big)$,
con $\mathcal V(x)\sim x^\beta$ e $\beta < 1$. Sia $\chi = \partial V/\partial h$ la
mobilità differenziale. In un estremo ($h_c' = 0$):

$$
\boxed{\;\lambda = v - h_c''(\alpha_0)\,\chi(h)\;},\qquad \chi\sim(h-h_c)^{\beta-1}\to\infty .
$$

In unità di velocità $\lambda$ **diverge** alla soglia invece di annullarsi.
Questo è il risultato di TKD.

È più chiaro in unità di forza. La soglia locale dipende dalla pendenza,
$h_c(s) = h_c + \tfrac12 h_c''s^2$, e quindi

$$
\mu^{-1}\partial_t u = (c/\mu)\,\partial_x^2u \;-\; \tfrac12\,h_c''(\alpha_0)\,(\partial_x u)^2 + h + F(x,u),
$$

con un coefficiente **finito** alla soglia, fissato da $h_c''$.

**Rilevanza.** Al punto fisso qEW, con $u\sim\ell^\zeta$, il rapporto tra termine
KPZ e termine elastico vale

$$
\frac{|h_c''|\,u^2/\ell^2}{(c/\mu)\,u/\ell^2} = \frac{|h_c''|\,u(\ell)}{c/\mu}\;\sim\;\ell^{\zeta},
$$

quindi il termine è fortemente rilevante ($\zeta = 1.25$). Con $u(\ell)\approx r_f(\ell/L_c)^\zeta$
il crossover sta a

$$
\boxed{\;\ell_\times \approx L_c\Big[\frac{c/\mu}{|h_c''|\,r_f}\Big]^{1/\zeta}\;}
$$

Per la parete edge a $a = 1$, con $|h_c''| \approx \tfrac43 h_c \approx 0.033$, si
ottiene $\ell_\times\approx 3.4\,L_c\approx 20$ siti. Se l'anisotropia efficace è
3 volte più piccola (§4), $\ell_\times\approx 50$–60. A $a = 0.5$ è circa 1.4 volte
più grande.

Va confrontato con la scala a cui le pendenze arrivano a 1 e compaiono i
ripiegamenti: $(r_f/L_c)(\ell/L_c)^{\zeta-1}\sim1$ dà $\ell\approx 80$–100. Per
le pareti edge a $a = 1$ il crossover KPZ dovrebbe quindi arrivare **prima** dei
ripiegamenti.

**Il segno decide la classe.**

* **Edge** ($h_c'' < 0$, massimo). Nell'equazione in unità di forza il termine
  $+\tfrac12|h_c''|(\partial_xu)^2$ facilita le parti inclinate, quindi
  $\lambda v > 0$. È qKPZ, cioè depinning per percolazione diretta (Tang &
  Leschhorn 1992; Buldyrev et al. 1992):

  $$\zeta = 0.633,\quad \nu_\parallel = 1.733,\quad z = 1,\quad \beta = \nu_\parallel(1-\zeta)\approx0.64 .$$

  Al punto fisso la dipendenza dalla pendenza diventa non analitica (TKD): un
  tratto inclinato di $s$ si scorrela dopo $\xi_s\sim s^{-1/(1-\zeta)}$, quindi

  $$
  h_c(s)-h_c(0)\;\propto\;-|s|^{1/[\nu_\parallel(1-\zeta)]} = -|s|^{1.57},
  \qquad v\big(h_c(0),s\big)\propto|s| .
  $$

  Siccome ruotare **b** equivale a inclinare la parete, si misura con
  `--burgers-angle 90±δ`.

* **Screw** ($h_c'' > 0$, minimo). I tratti inclinati sono più duri, quindi
  $\lambda v < 0$ ("qKPZ negativa"). L'orientazione piatta non è stabile: la
  parete ancorata si spezza in faccette lungo orientazioni dure (Jeong, Kahng &
  Kim, PRL 77, 5094 (1996), *da verificare*). Su scala grande la soglia è
  fissata dall'orientazione più dura, quindi ci si aspetta
  $h_c^{\text{screw}}(L)\to h_c^{\text{edge}}$ al crescere di $L$.

* **Orientazioni generiche** ($h_c'\neq0$). C'è anche un termine lineare
  $V_0's$ (una deriva laterale, eliminabile) e il termine $h_c'^2\,\partial_h\chi$,
  che è negativo e diverge ancora più forte. Qui KPZ non si spegne mai.

## 6. Cosa dicono già i dati

1. **a = 0 è qEW** ($\zeta = 1.24\pm0.04$, $\nu = 1.31\pm0.09$, $\nu = 1/(2-\zeta)$),
   come richiesto da $\lambda = v\to0$.
2. **Segno di $h_c''$**: a $L = 128$, $a = 1$, abbiamo $h_d^{\text{screw}}\approx0.011$
   contro $h_d^{\text{edge}}\approx0.024$–0.029. Edge è la direzione dura, come
   previsto.
3. **Le pareti screw si sfaccettano.** Distribuzione dell'angolo della tangente
   rispetto a **b** (pareti critiche `flat_D0.04_b0_A2`, $L = 1024$, 4 pareti per valore di $a$
   (a $L = 512$ e con lisciatura su 4 siti il quadro è lo stesso), lisciate su 16
   siti di arco; 0° = orientazione screw, ±90° = orientazione edge):

   | a | 0–15° | 15–30° | 30–45° | 45–60° | 60–75° | 75–90° | ⟨\|θ\|⟩ |
   |---|---|---|---|---|---|---|---|
   | 0 | 0.22 | 0.22 | 0.18 | 0.16 | 0.12 | 0.10 | 38° |
   | 0.5 | 0.06 | 0.09 | 0.13 | 0.18 | 0.24 | 0.31 | 58° |
   | 1 | 0.02 | 0.03 | 0.07 | 0.18 | 0.32 | 0.39 | 67° |

   A $a = 1$ una parete "screw" è fatta per circa il 70% di tratti entro 30°
   dall'orientazione edge, e solo per il 2% di tratti vicini a quella screw.
   È la sfaccettatura prevista per $\lambda v < 0$.

   Spiega anche le altre anomalie viste:
   * l'estensione verticale di $1.3\,L$;
   * $P/L\approx3.6$;
   * il 17–21% di colonne multivalore;
   * $h_d^{\text{screw}}$ che cresce con $L$ (0.0193 → 0.0235 da $L = 256$ a 2048)
     verso il valore edge.
4. **Domanda aperta.** Lo $\zeta$ intrinseco delle pareti screw è $\approx 1.2$–1.3,
   tipo qEW. Per una parete fatta di faccette edge, se le faccette fossero
   qKPZ, ci aspetteremmo qualcos'altro. O le faccette sono più corte di
   $\ell_\times$, o lo spettro di $\psi$ della parete sfaccettata misura la
   geometria delle faccette e non la rugosità di una singola parete. Da
   guardare con più statistica e con aspect 4.

## 7. Previsioni da verificare

| test | previsione | come |
|---|---|---|
| ζ, ν pareti edge (b90) | 0.633, 1.733 (qKPZ) per $\ell\gtrsim\ell_\times\approx20$–60 a $a = 1$ | `slurm_fss_flat.sh 90`, `analyze_flat.py` |
| ζ locale pareti edge | da ≈1.25 a ≈0.63 attorno a $q\approx2\pi/\ell_\times$ | pannello "local exponent" |
| ripiegamenti | pareti edge molto meno multivalore delle screw | `frac_multivalued`, $P/L$ |
| $h_c(\alpha)$ vicino all'edge | massimo; $h_c(0)-h_c(\delta)\propto|\delta|^{1.57}$, non $\delta^2$ | `slurm_fss_flat.sh 85`, `80`, `70` ($a = 1$, $L$ fisso) |
| $h_c(\alpha)$ vicino alla screw | minimo analitico; $h_c^{\text{screw}}(L)\to h_c^{\text{edge}}$ | angoli 5, 10, 20; serie in $L$ |
| sfaccettatura screw | istogramma di θ concentrato su orientazioni dure, a ogni $L$ | come in §6.3, anche con aspect 4 |
| controllo a = 0 | nessuna dipendenza da α, ζ = 1.25 per ogni angolo | già fatto |

## Riferimenti

* L.-H. Tang, M. Kardar, D. Dhar, *Driven depinning in anisotropic media*, PRL 74, 920 (1995), cond-mat/9409002.
* L. A. N. Amaral, A.-L. Barabási, H. E. Stanley, *Universality classes for interface growth with quenched disorder*, PRL 73, 62 (1994), cond-mat/9403089 (test dell'inclinazione).
* L.-H. Tang, H. Leschhorn, PRA 45, R8309 (1992); S. V. Buldyrev et al., PRA 45, R8313 (1992) (depinning per percolazione diretta).
* H. Jeong, B. Kahng, D. Kim, PRL 77, 5094 (1996) (qKPZ con λ negativo, faccette) — *riferimento da verificare*.
* T. Nattermann, S. Stepanow, L.-H. Tang, H. Leschhorn, J. Phys. II 2, 1483 (1992); O. Narayan, D. S. Fisher, PRB 48, 7030 (1993) (qEW, RB e RF alla depinning).
* E. E. Ferrero, S. Bustingorry, A. B. Kolton, PRE 87, 032122 (2013) (esponenti qEW in 1+1).
* A. B. Kolton, E. E. Ferrero, A. Rosso, arXiv:2306.13415 (GL, ripiegamenti, $\ell_o$).
