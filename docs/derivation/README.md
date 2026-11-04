# Derivation of the DrillPyze coupled axial–torsional model

This document derives the model implemented in `DrillPyze` (kernel:
`src/DrillPyze/kernel.py`), from the governing equations to the discrete
scheme.

* The torsional model follows Aarsnes & Shor (2018) and the control papers of
  Aarsnes, di Meglio & Shor (2018) and Auriol et al. (2022).
* The axial model uses an identical formulation, with Young's modulus in place
  of the shear modulus.
* In the base model the two directions are coupled **only at the drill bit**,
  through a Coulomb-type bit–rock interaction.
* An opt-in extension that also couples them along the string is described in
  §5.4.

All quantities are SI. Notation is collected in §11.

---

## 1. Assumptions

1. The drillstring is a slender, linearly elastic rod of length $L$, with
   curvilinear coordinate $x \in [0, L]$ measured from the top drive ($x=0$)
   to the bit ($x=L$).
2. Only axial and torsional motion is modelled; lateral (bending/whirl)
   dynamics are not. The string lies on the low or high side of the borehole
   as determined by the soft-string normal force (§5.1).
3. Material properties are piecewise constant over sections (drill pipe, HWDP,
   collars, …). Tool-joint mass is included through an effective density
   $\rho = m/A$, where $m$ is the mass per unit length and $A$ the body area.
4. Small strains, so axial and torsional wave speeds are independent of load.
5. Drilling fluid enters through buoyancy (weight) and a linear viscous damping
   term. Hydraulics do not otherwise enter the mechanics.
6. Contact with the borehole is Coulomb friction with distinct static and
   kinetic coefficients, separately for cased and open hole.

---

## 2. Kinematics and constitutive relations

### 2.1 Torsion

Let $\phi(t,x)$ be the angular displacement and
$\omega = \partial\phi/\partial t$ the angular velocity. Following Aarsnes &
Shor (2018), the torque transmitted downward through a cross-section is
proportional to the twist per unit length:

$$
\tau(t,x) = J G \,\frac{\phi(t,x) - \phi(t,x+\mathrm{d}x)}{\mathrm{d}x} = -JG\,\frac{\partial \phi}{\partial x},
$$

where $J = \tfrac{\pi}{32}(D_o^4 - D_i^4)$ is the polar moment of area and $G$
the shear modulus. Differentiating in time gives the **compatibility
equation**

$$
\frac{\partial \tau}{\partial t} + JG\,\frac{\partial \omega}{\partial x} = 0. \tag{2.1}
$$

### 2.2 Axial (identical structure, $G \to E$, $J \to A$)

Let $u(t,x)$ be the axial displacement (positive downhole) and
$v = \partial u/\partial t$ the axial velocity. The axial force transmitted
downward is the **compressive** force

$$
q(t,x) = -EA\,\frac{\partial u}{\partial x},
$$

with $A = \tfrac{\pi}{4}(D_o^2 - D_i^2)$, so $q > 0$ is compression and the
effective tension is $T = -q$. Differentiating in time:

$$
\frac{\partial q}{\partial t} + EA\,\frac{\partial v}{\partial x} = 0. \tag{2.2}
$$

Choosing $q$ (compression) rather than $T$ makes the axial equations term for
term identical to the torsional ones under $(\omega, \tau, J, G) \to (v, q, A, E)$.

---

## 3. Balance of momentum

Consider an element $[x, x + \mathrm{d}x]$.

### 3.1 Torsion

The rotational inertia of the element is $\rho J\,\mathrm{d}x$. The net
transmitted torque is $\tau(x) - \tau(x+\mathrm{d}x)$, and the distributed
torques are viscous damping and Coulomb wall friction:

$$
\rho J\,\frac{\partial \omega}{\partial t} + \frac{\partial \tau}{\partial x} = -k_t\,\rho J\,\omega - \mathcal{F}_t(t,x), \tag{3.1}
$$

where $k_t$ (1/s) is the viscous damping rate due to the drilling fluid, and
$\mathcal{F}_t$ is the frictional torque per unit length (an inclusion, §5).

### 3.2 Axial

The mass of the element is $\rho A\,\mathrm{d}x$. In addition to the
transmitted force, the element carries its buoyed weight per length $w_b$
projected on the pipe axis, $w_b \cos\theta$, where $\theta(x)$ is the
inclination:

$$
\rho A\,\frac{\partial v}{\partial t} + \frac{\partial q}{\partial x} = w_b\cos\theta - k_a\,\rho A\, v - \mathcal{F}_a(t,x), \tag{3.2}
$$

with

$$
w_b = m\,g\left(1 - \frac{\rho_\text{mud}}{\rho_\text{steel}}\right),
$$

$k_a$ the axial damping rate, and $\mathcal{F}_a$ the axial friction force per
unit length.

### 3.3 Wave equations

Eliminating $\tau$ from (2.1) and (3.1) in a uniform section without sources
gives $\partial_t^2\phi = c_t^2\,\partial_x^2\phi$. Likewise, (2.2) and (3.2)
give $\partial_t^2 u = c_a^2\,\partial_x^2 u$. The wave speeds are

$$
c_t = \sqrt{G/\rho}, \qquad c_a = \sqrt{E/\rho}.
$$

For steel, $c_t \approx 3{,}180$ m/s and $c_a \approx 5{,}140$ m/s.

---

## 4. Characteristic form

### 4.1 Hyperbolic system

Write the torsional system as $U_t + \mathbf{A}\,U_x = S$ with
$U = (\tau, \omega)^\top$:

$$
\mathbf{A} = \begin{pmatrix} 0 & JG \\ 1/(\rho J) & 0 \end{pmatrix},
\qquad
S = \begin{pmatrix} 0 \\ S_t/(\rho J) \end{pmatrix},
\qquad
S_t = -k_t\rho J\omega - \mathcal{F}_t .
$$

The eigenvalues are $\pm\sqrt{G/\rho} = \pm c_t$.

### 4.2 Riemann invariants and impedance

Define the **characteristic impedance**

$$
Z_t = \frac{JG}{c_t} = \rho J\,c_t = J\sqrt{G\rho}
\qquad\text{(axial: } Z_a = A\sqrt{E\rho}\text{)},
$$

and the Riemann invariants (Aarsnes & Shor 2018, eq. 11)

$$
\alpha = \omega + \frac{\tau}{Z_t}, \qquad \beta = \omega - \frac{\tau}{Z_t}. \tag{4.1}
$$

Using (2.1), (3.1) and $JG/Z_t = c_t$, $1/(\rho J) = c_t/Z_t$:

$$
\partial_t\alpha
= \partial_t\omega + \tfrac{1}{Z_t}\partial_t\tau
= \frac{S_t - \partial_x\tau}{\rho J} - \frac{JG}{Z_t}\partial_x\omega
= -c_t\,\partial_x\!\left(\omega + \tfrac{\tau}{Z_t}\right) + \frac{S_t}{\rho J},
$$

so

$$
\frac{\partial \alpha}{\partial t} + c_t\,\frac{\partial \alpha}{\partial x} = \frac{S_t}{\rho J},
\qquad
\frac{\partial \beta}{\partial t} - c_t\,\frac{\partial \beta}{\partial x} = \frac{S_t}{\rho J}. \tag{4.2}
$$

$\alpha$ is a wave travelling **down** the string and $\beta$ a wave
travelling **up**. The physical variables are recovered as
$\omega = (\alpha+\beta)/2$ and $\tau = Z_t(\alpha-\beta)/2$.

The axial invariants are
$\alpha_a = v + q/Z_a$, $\beta_a = v - q/Z_a$, which travel at $\pm c_a$
with source $S_a/(\rho A)$.

### 4.3 Junctions between sections

At a change of section at $x_i$ (e.g. pipe to collar), continuity requires
$\omega^- = \omega^+$ and $\tau^- = \tau^+$. Given the incoming invariants,
$\alpha^-$ from above and $\beta^+$ from below, these two conditions give the
interface state

$$
\tau^\ast = \frac{\alpha^- - \beta^+}{1/Z^- + 1/Z^+},
\qquad
\omega^\ast = \alpha^- - \frac{\tau^\ast}{Z^-}. \tag{4.3}
$$

The outgoing waves then follow from (4.1). For
$\bar Z = Z^+/Z^-$ (for one material, $\bar Z = J_c/J_p$), this reproduces the
transmission and reflection relations of Aarsnes & Shor (2018), eqs. 16–17.
Equation (4.3) is also the exact Riemann solution used at **every** cell face
in the numerical scheme (§8).

### 4.4 Gravity: well-balanced split

The weight term in (3.2) is time-invariant. Write $q = q_s(x) + q_d(t,x)$,
where $q_s$ is the compression of the frictionless string hanging freely with
the bit off bottom:

$$
\frac{\mathrm{d}q_s}{\mathrm{d}x} = w_b\cos\theta, \qquad q_s(L) = 0
\;\;\Rightarrow\;\;
q_s(x) = -\int_x^L w_b\cos\theta\,\mathrm{d}\xi \;(\le 0,\ \text{i.e. tension}).
$$

Then $q_d$ satisfies (2.2) and (3.2) **without** the weight term. Only $q_d$ is
transported. This makes the discrete scheme well balanced: a hanging string at
rest is an exact discrete equilibrium. Hookload and WOB are reconstructed as
$\text{HKLD} = -(q_s(0) + q_d(0))$ and $\text{WOB} = q(L) = q_d(L)$.

---

## 5. Along-string friction

### 5.1 Normal force (soft string)

For a cable element with unit tangent $\mathbf{t}(x)$ and effective tension
$T$, equilibrium normal to the pipe gives the contact force per unit length
(Johancsik et al. 1984, in vector form):

$$
F_N(x) = \left\| T\,\frac{\mathrm{d}\mathbf{t}}{\mathrm{d}x} + w_b\,\mathbf{z}_\perp \right\|,
\qquad
\mathbf{z}_\perp = \hat{\mathbf{z}} - (\hat{\mathbf{z}}\cdot\mathbf{t})\,\mathbf{t}, \tag{5.1}
$$

where $\hat{\mathbf z}$ points vertically down. The tangent is interpolated
along the minimum-curvature arcs of the survey. In a 2-D build section,
(5.1) reduces to $|T\,\mathrm{d}\theta/\mathrm{d}x - w_b\sin\theta|$: tension
pulls the pipe against the high side and gravity against the low side. In a
straight inclined hole it reduces to $w_b\sin\theta$, the normal force of
Aarsnes & Shor (2018), eq. 6.

In the base model the normal force is evaluated with the tension of the static
hanging string, $T_s = -q_s$, so $F_N = F_N(x)$ does not depend on the dynamic
state. It is re-evaluated whenever the string has moved more than
`refresh_distance` along the hole.

### 5.2 Coulomb inclusion: torsion (Aarsnes & Shor 2018, eqs. 4–5)

With friction factors $\mu_s \ge \mu_k$ (static, kinetic), contact radius
$r_c$ and threshold $\omega_c$, define the friction torque limits

$$
\mathcal{F}_{c,s}(x) = \mu_s\,F_N(x)\,r_c(x),
\qquad
\mathcal{F}_{c,k}(x) = \mu_k\,F_N(x)\,r_c(x).
$$

The friction torque per length is the set-valued map

$$
\mathcal{F}_t(\omega, x) \in
\begin{cases}
\{\,\mathcal{F}_{c,k}(x)\,\}, & \omega > \omega_c,\\[2pt]
[-\mathcal{F}_{c,s}(x),\ \mathcal{F}_{c,s}(x)], & |\omega| \le \omega_c,\\[2pt]
\{-\mathcal{F}_{c,k}(x)\}, & \omega < -\omega_c .
\end{cases} \tag{5.2}
$$

When $|\omega| \le \omega_c$, $\mathcal{F}_t$ takes the value that keeps the
element at rest,

$$
\mathcal F_t = -\frac{\partial\tau}{\partial x} - k_t\rho J\omega,
$$

provided this lies within the static limit; otherwise it takes the bound.
With $f_\text{rat} = \mu_k/\mu_s$ this is exactly the published model.

### 5.3 Coulomb inclusion: axial (identical)

$$
\mathcal{F}_a(v, x) \in
\begin{cases}
\{\,\mu_k F_N\,\}, & v > v_c,\\[2pt]
[-\mu_s F_N,\ \mu_s F_N], & |v| \le v_c,\\[2pt]
\{-\mu_k F_N\}, & v < -v_c .
\end{cases} \tag{5.3}
$$

In the base model (5.2) and (5.3) are **independent**: neither depends on the
other direction's velocity, and $F_N$ does not depend on the dynamic tension.
The axial and torsional systems are therefore uncoupled along the string, and
**the only coupling is at the bit (§6.3)**.

### 5.4 Optional along-string coupling (`along_string_coupling=True`)

For later studies, an opt-in model treats friction as a single 2-D Coulomb
force opposing the sliding velocity of the pipe surface,
$\mathbf{V} = (v,\ \omega r_c)$:

* static if $(v/v_c)^2 + (\omega/\omega_c)^2 \le 1$, with
  $\|\boldsymbol{\mathcal F}\| \le \mu_s F_N$;
* otherwise $\boldsymbol{\mathcal F} = \mu_k F_N\,\mathbf{V}/\|\mathbf V\|$.

The normal force is evaluated with the live tension $T = -(q_s + q_d)$. This
reproduces velocity-split torque and drag (rotation reduces drag), and lets
axial loads change torque along the string.

---

## 6. Boundary conditions

At each boundary one invariant arrives from the interior. The boundary model
supplies one relation, which closes the system.

### 6.1 Top drive (torsion)

At $x = 0$ the incoming invariant is $\beta_1$ (travelling up). The
characteristic relation $\beta = \omega - \tau/Z_t$ gives the pipe torque at
the top drive:

$$
\tau_0 = Z_t\,(\omega_{TD} - \beta_1). \tag{6.1}
$$

The top drive is a rigid inertia $I_{TD}$ driven by motor torque $\tau_m$:

$$
I_{TD}\,\dot\omega_{TD} = \tau_m - \tau_0 = \tau_m - Z_t\,(\omega_{TD} - \beta_1). \tag{6.2}
$$

Seen from the string, (6.2) is a damper of impedance $Z_t$ in parallel with
the inertia. It is integrated implicitly in $\omega_{TD}$, which is
unconditionally stable:

$$
\omega_{TD}^{k+1} = \frac{\dfrac{I_{TD}}{\Delta t}\,\omega_{TD}^{k} + \tau_m + Z_t\,\beta_1}{\dfrac{I_{TD}}{\Delta t} + Z_t}.
$$

The motor torque comes from a controller (§7), from an external source, or is
implied by a prescribed $\omega_{TD}(t)$ via
$\tau_m = \tau_0 + I_{TD}\,\dot\omega_{TD}$.

### 6.2 Drawworks (axial)

At $x = 0$ the block velocity $v_\text{blk}$ (positive lowering) is
prescribed. The incoming axial invariant $\beta_{a,1}$ gives

$$
q_d(0) = Z_a\,(v_\text{blk} - \beta_{a,1}),
\qquad
\text{HKLD} = -\bigl(q_s(0) + q_d(0)\bigr).
$$

In autodriller mode, $v_\text{blk}$ is the output of a PI controller acting on
the surface WOB estimate $\text{HKLD}_\text{off-bottom} - \text{HKLD}$.

### 6.3 Drill bit: the axial–torsional coupling

At $x = L$ the incoming invariants are $\alpha_t$ (torsion) and $\alpha_a$
(axial). The characteristic relations are

$$
\omega_b = \alpha_t - \frac{\text{TOB}}{Z_t}, \qquad v_b = \alpha_a - \frac{\text{WOB}}{Z_a}. \tag{6.3}
$$

Both directions are closed by the bit–rock interaction law below. Because TOB
depends on WOB and the penetration rate depends on $\omega_b$, this law is the
coupling between the two systems.

**Off bottom** (bit depth $<$ hole depth): $\text{WOB} = \text{TOB} = 0$, so
$\omega_b = \alpha_t$ and $v_b = \alpha_a$. Each free end reflects its incident
wave with the same sign.

**Unilateral contact.** On bottom, the bit can only push on the rock:
$\text{WOB} \ge 0$. If $\alpha_a \le 0$ (the string is pulling the bit up),
the bit lifts off: $\text{WOB} = 0$ and $v_b = \alpha_a$.

**Torque: Coulomb-type law with sticking.** With bit radius $r_b$ (default
$D_\text{bit}/3$), static and kinetic bit friction $\mu_{bs} \ge \mu_{bk}$,
sticking torque $T_{s0}$ and threshold $\omega_c$:

$$
\text{TOB} \in
\begin{cases}
[-(T_{s0} + \mu_{bs}\,r_b\,\text{WOB}),\ T_{s0} + \mu_{bs}\,r_b\,\text{WOB}], & \text{stuck } (|\omega_b| \le \omega_c),\\[3pt]
\mu_{bk}\,r_b\,\text{WOB}\;\operatorname{sign}(\omega_b), & \text{sliding}.
\end{cases} \tag{6.4}
$$

* **Stuck.** The bit stays at rest ($\omega_b = 0$, so
  $\text{TOB} = Z_t\alpha_t$ from (6.3)) as long as
  $|Z_t\alpha_t| \le T_{s0} + \mu_{bs} r_b \text{WOB}$. Breakout occurs when the
  torque delivered by the string exceeds this limit.
* **Sliding.** The reactive torque is linear in WOB. TOB is limited to
  $|Z_t\alpha_t|$, so the bit is never driven backwards within one step.

**Penetration law.** On bottom, the rock advances at

$$
v_b = \text{ROP} = k_r\,\max(\text{WOB} - W_0,\ 0)\;|\omega_b|. \tag{6.5}
$$

Combining with (6.3): if $Z_a\alpha_a \le W_0$ or $\omega_b = 0$, the bit does
not advance ($v_b = 0$, $\text{WOB} = Z_a\alpha_a$). Otherwise

$$
\text{WOB} = \frac{\alpha_a + k_r W_0 |\omega_b|}{1/Z_a + k_r|\omega_b|}.
$$

**Solution.** Equations (6.3)–(6.5) are solved by fixed-point iteration on
$\omega_b$, which converges in a few iterations because each step changes the
state little. The hole depth is updated as
$\max(\text{hole depth},\ \text{bit depth})$.

**Coupling mechanisms.**
* Torsion depends on axial: TOB $\propto$ WOB in (6.4), so axial waves that
  change WOB modulate the torsional boundary and can trigger or suppress
  stick-slip.
* Axial depends on torsion: ROP $\propto |\omega_b|$ in (6.5). A stalled bit
  stops penetrating, WOB rises, and the stick limit rises with it.
* Steady drilling at block speed $V$ and bit speed $\Omega$ therefore gives

$$
\text{WOB} = W_0 + \frac{V}{k_r\Omega}, \qquad \text{TOB} = \mu_{bk} r_b\,\text{WOB}.
$$

A Detournay-type cutting/contact law with depth-of-cut delay (the regenerative
effect) is a planned refinement that would replace (6.4)–(6.5).

---

## 7. Top-drive controllers (Aarsnes, di Meglio & Shor 2018)

All three industrial controllers fit one block diagram. A rate-limited set
point $\omega_{SP}$, a filtered measured speed $\hat\omega_0$ (time constant
$t_{en}$) and an optional impedance feedback of band-pass filtered pipe torque
give

$$
e = \omega_{SP} - \hat\omega_0 - Z\,\frac{s}{(s + 1/t_{hp})(1 + t_{lp}s)}\,\tau_0,
\qquad
\tau_m = K_p\,e + K_i\!\int e\,\mathrm{d}t,
$$

with saturation $|\tau_m| \le \tau_\text{max}$ and anti-windup. With the pipe
impedance $\zeta_p = Z_t(0)$:

| controller | $K_p$ | $K_i$ | $Z$ |
|---|---|---|---|
| stiff | $100\,\zeta_p$ | $5\,I_{TD}$ | 0 |
| SoftTorque | $4\,\zeta_p$ | $(2\pi f_c)^2 I_{TD}$ | 0 |
| ZTorque | $100\,\zeta_p$ | $5\,I_{TD}$ | $1/\zeta_p$ |

For SoftTorque, $f_c$ is the frequency of minimum reflection, usually the
first torsional mode. ZTorque requires a stiff inner loop so the top drive
actually follows the impedance-corrected reference; combining impedance
feedback with soft gains is unstable in this model.

The top-side reflection coefficient, ignoring the impedance path, is

$$
R(\omega) = \left|\frac{\bar C(j\omega) - \zeta_p}{\bar C(j\omega) + \zeta_p}\right|,
\qquad
\bar C(s) = K_p + \frac{K_i}{s} + I_{TD}\,s .
$$

---

## 8. Numerical scheme

### 8.1 Finite volumes with Godunov fluxes

The string is divided into cells $j = 1..N$ of length $h_j \le \Delta x$,
aligned with section boundaries so each cell has uniform properties. Cell
averages $(\omega_j, \tau_j)$ and $(v_j, q_{d,j})$ are stored.

At each face $j+\tfrac12$ the exact Riemann solution (4.3) is evaluated from
$\alpha_j$ (cell above) and $\beta_{j+1}$ (cell below). The boundary faces use
§6. The conservative update is

$$
\tau_j^{*} = \tau_j^k - \frac{\Delta t\,J_jG_j}{h_j}\left(\omega^\ast_{j+1/2} - \omega^\ast_{j-1/2}\right),
\qquad
\omega_j^{*} = \omega_j^k - \frac{\Delta t}{\rho_j J_j h_j}\left(\tau^\ast_{j+1/2} - \tau^\ast_{j-1/2}\right),
$$

and identically for $(q_d, v)$ with $(EA, \rho A)$.

Inside a uniform section this is exactly the first-order upwind scheme on the
Riemann invariants (Aarsnes & Shor 2018, eqs. 23–24). At junctions it enforces
(4.3) without special cases. First-order upwinding is used deliberately: the
set-valued friction introduces temporal discontinuities that make higher-order
schemes oscillate.

### 8.2 Time step

Stability requires $\Delta t \le h_j / c$ for both waves. The axial wave is
faster, so

$$
\Delta t = \text{CFL}\cdot\min_j \frac{h_j}{c_{a,j}}, \qquad \text{CFL} = 0.95,
$$

and the torsional Courant number is $\text{CFL}\,c_t/c_a \approx 0.59$. At
this Courant number upwinding is slightly diffusive, but the tests show
negligible amplitude loss for the low modes that dominate stick-slip (§10).

### 8.3 Source terms (operator splitting)

After the flux update, each cell applies its sources:

1. **Damping, implicitly:**
   $\omega^{**} = \omega^*/(1 + k_t\Delta t)$ and
   $v^{**} = v^*/(1 + k_a\Delta t)$.
2. **Friction, by projection** (Aarsnes & Shor 2018, eqs. 20–22). The friction
   that would bring the cell exactly to rest in this step is
   $\tilde{\mathcal F}_t = \rho J\,\omega^{**}/\Delta t$
   (and $\tilde{\mathcal F}_a = \rho A\,v^{**}/\Delta t$).
   The applied friction is

$$
\mathcal F = \operatorname{sign}(\tilde{\mathcal F})\,\min\!\left(|\tilde{\mathcal F}|,\ \mathcal F_c\right),
\qquad
\mathcal F_c = \begin{cases}\mathcal F_{c,s}, & |\omega^k| \le \omega_c\\ \mathcal F_{c,k}, & \text{otherwise,}\end{cases}
$$

   followed by $\omega^{k+1} = \omega^{**} - \Delta t\,\mathcal F/(\rho J)$.

This realises the inclusion (5.2) exactly in discrete form:
* a cell sticks exactly ($\omega^{k+1} = 0$) when the static limit suffices;
* friction never reverses the motion within a step;
* the regime is selected from the previous step's velocity, as in eq. 22 of
  the paper.

### 8.4 Moving string

The grid is attached to the pipe. As the bit advances by $s$, cell $j$ moves
to $x_j + s$. Every `refresh_distance` of travel, the cell geometry, the
cased/open friction zones, $q_s$ and $F_N$ are re-evaluated. Total compression
is kept continuous across a refresh by setting
$q_d \leftarrow q_d + q_s^\text{old} - q_s^\text{new}$.

### 8.5 Pausing

The state consists of the cell arrays, the top-drive, controller and filter
states, and the bit and hole depths. All of it persists between calls, and
inputs are held constant within a call. A sequence of short `step` calls is
therefore bit-identical to one long call; the steps are counted in integers
against the accumulated requested time.

---

## 9. Steady-state limits and analytic results

These limits connect the model to classical results and are used in the tests.

* **Torque and drag.** Setting $\partial_t = 0$ in (3.1)–(3.2) with kinetic
  friction gives the soft-string torque and drag equations, integrated from
  the bit up:

$$
\frac{\mathrm{d}T}{\mathrm{d}x} = \mathcal F_a - w_b\cos\theta,
\qquad
\frac{\mathrm{d}\tau}{\mathrm{d}x} = -\mathcal F_t,
\qquad
T(L) = -\text{WOB},\ \tau(L) = \text{TOB}.
$$

  This is `statics.torque_drag`. The decoupled model is compared with it
  using the static normal force (`normal_force=sim.fn0`).

* **Break-off torque** (bit off bottom; Aarsnes & Shor 2018, eqs. 26–27):

$$
\tau_\text{max}(x) = \int_x^L \mathcal F_{c,s}\,\mathrm{d}\xi .
$$

* **Stick-slip amplitude and period** (Aarsnes & Shor 2018, eqs. 39–40):

$$
\Delta\tau_m \approx 2\,\tau_\text{max}(0)\,(1 - f_\text{rat}),
\qquad
t_p \approx \frac{\Delta\phi(0)}{\omega_{TD}} + t_{fp},
\qquad
\Delta\phi(0) = 2(1-f_\text{rat})\int_0^L\frac{\tau_\text{max}}{JG}\,\mathrm{d}x .
$$

* **Free vibration** of a uniform string with a fixed top and free bit: wave
  arrival at the bit after $L/c$, and fundamental frequency $f_1 = c/(4L)$.

* **Bit.**
  * Stuck-bit breakout at $\text{TOB} = T_{s0} + \mu_{bs} r_b\,\text{WOB}$.
  * While sliding, $\text{TOB}/\text{WOB} = \mu_{bk} r_b$.
  * Steady drilling at $\text{WOB} = W_0 + V/(k_r\Omega)$.

---

## 10. Verification summary (`tests/`)

| check | result |
|---|---|
| torsional / axial wave arrival at $L/c$ | within 1% |
| fundamental period $4L/c_t$ | within 0.2%; amplitude loss < 1% over 8 round trips |
| hanging string at rest | exactly at rest (well balanced) |
| steady rotation, pick-up, slack-off vs torque & drag | within 2% |
| decoupled model: rotation does not change pick-up drag | within 0.2% |
| opt-in coupling vs velocity-split torque & drag | within 2–3% |
| off-bottom stick-slip only when $\mu_s > \mu_k$ | yes |
| ZTorque suppresses the stiff-controller limit cycle | yes |
| bit breakout torque, TOB/WOB slope, drill-off at constant block | exact / 1% |
| pause–resume, external Python PI vs internal PI | exact / 2% |

---

## 11. Notation

| symbol | meaning | unit |
|---|---|---|
| $x$ | measured depth from surface | m |
| $\omega, \tau$ | angular velocity, torque (transmitted downward) | rad/s, N m |
| $v, q$ | axial velocity (+ down), compressive force; $T=-q$ | m/s, N |
| $J, A$ | polar moment of area, cross-sectional area | m⁴, m² |
| $G, E, \rho$ | shear modulus, Young's modulus, effective density | Pa, Pa, kg/m³ |
| $c_t, c_a$ | torsional, axial wave speed | m/s |
| $Z_t, Z_a$ | torsional, axial characteristic impedance | N m s, N s/m |
| $\alpha, \beta$ | down- and up-going Riemann invariants | rad/s or m/s |
| $w_b$ | buoyed weight per length | N/m |
| $\theta, \mathbf{t}$ | inclination, unit tangent | rad, – |
| $F_N$ | normal force per length | N/m |
| $\mu_s, \mu_k$ | static, kinetic friction factors (cased/open) | – |
| $\omega_c, v_c$ | static → kinetic thresholds | rad/s, m/s |
| $k_t, k_a$ | viscous damping rates | 1/s |
| $I_{TD}, \tau_m, \tau_0$ | top-drive inertia, motor torque, pipe torque at top | kg m², N m, N m |
| $r_b, \mu_{bs}, \mu_{bk}, T_{s0}$ | bit radius, bit static/kinetic friction, sticking torque | m, –, –, N m |
| $k_r, W_0$ | penetration coefficient, threshold WOB | m/(N rad), N |

## References

* Aarsnes, U.J.F. and Shor, R.J., 2018. Torsional vibrations with bit off bottom: Modeling, characterization and field data validation. *J. Petroleum Science and Engineering*, 163, 712–721.
* Aarsnes, U.J.F., di Meglio, F. and Shor, R.J., 2018. Benchmarking of industrial stick-slip mitigation controllers. *IFAC-PapersOnLine*, 51(8), 233–238.
* Auriol, J., Boussaada, I., Shor, R.J., Mounier, H. and Niculescu, S.-I., 2022. Comparing advanced control strategies to eliminate stick-slip oscillations in drillstrings. *IEEE Access*, 10, 10949–10969.
* Johancsik, C.A., Friesen, D.B. and Dawson, R., 1984. Torque and drag in directional wells — prediction and measurement. *JPT*, 36, 987–992.
* Shor, R.J. and Aarsnes, U.J.F. Coupled axial-torsional vibrations in a drillstring: modelling and validation. Draft manuscript.
