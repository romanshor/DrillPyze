# Derivation of the DrillPyze coupled axial–torsional model

This document derives the model implemented in `DrillPyze` (kernel:
`src/DrillPyze/kernel.py`), from the governing equations to the discrete
scheme.

* The torsional model follows Aarsnes & Shor (2018) and the control papers of
  Aarsnes, di Meglio & Shor (2018) and Auriol et al. (2022).
* The axial model uses an identical formulation, with Young's modulus in place
  of the shear modulus.
* The two directions are cross-coupled in two places:
  1. **along the string**, through a single 2-D Coulomb friction acting on the
     combined axial and rotational sliding velocity, with a normal force that
     follows the live tension (§5);
  2. **at the drill bit**, through a Coulomb-type bit–rock interaction (§6.3).
* A decoupled variant (`along_string_coupling=False`), in which the only
  coupling is at the bit, is kept for comparison (§5.4).

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
5. Drilling fluid enters through buoyancy (weight), linear viscous damping
   (steady Couette drag plus an oscillatory Stokes boundary layer) and, with
   the pumps off, a Bingham yield-stress drag (§5.6). Steel damping enters
   through a Rayleigh / Kelvin–Voigt term (§5.6). Hydraulics do not otherwise
   enter the mechanics.
6. Contact with the borehole is Coulomb friction with distinct static and
   kinetic coefficients, separately for cased and open hole.
7. **Frequency range.** The damping models are fitted at the first mode and at
   $f_\max = 10$ Hz. Together with the numerical dissipation of the
   first-order scheme, this limits the model to vibrations below 10 Hz
   (§5.6.6). High-frequency torsional oscillations (HFTO, 50–500 Hz) and bit
   bounce harmonics are out of scope.

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

where $k_t$ (1/s) stands for the linear damping operators of §5.6 (fluid and
material; in the original model a single constant rate), and
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

In the coupled model (default) the normal force is evaluated with the **live**
tension $T = -(q_s + q_d)$. Axial loading (pick-up, slack-off, WOB transfer)
therefore changes the normal force, and through it the rotary friction torque.

In the decoupled variant (§5.4) it is evaluated with the tension of the static
hanging string, $T_s = -q_s$, so $F_N = F_N(x)$ does not depend on the dynamic
state. In both cases it is re-evaluated whenever the string has moved more
than `refresh_distance` along the hole.

### 5.2 Cross-coupled Coulomb inclusion (default)

The pipe surface at $x$ slides against the borehole with the velocity

$$
\mathbf V = \begin{pmatrix} v \\ \omega\, r_c \end{pmatrix}
\qquad\text{(axial, tangential)},
$$

where $r_c$ is the contact radius (tool joint or stabiliser). Coulomb's law
for a contact with normal force $F_N$ is a **single** friction force in the
plane of sliding, $\boldsymbol{\mathcal F} = (\mathcal F_a,\ \mathcal F_\theta)$,
with

$$
\boldsymbol{\mathcal F} \in
\begin{cases}
\left\{\,\mu_k F_N\,\dfrac{\mathbf V}{\|\mathbf V\|}\,\right\}, & \text{sliding},\\[8pt]
\left\{\,\boldsymbol{\mathcal F} : \|\boldsymbol{\mathcal F}\| \le \mu_s F_N\,\right\}, & \text{sticking}.
\end{cases} \tag{5.2}
$$

The regime is set by the normalised sliding speed: the contact is **sticking**
when

$$
\left(\frac{v}{v_c}\right)^2 + \left(\frac{\omega}{\omega_c}\right)^2 \le 1,
$$

and sliding otherwise. The axial friction force per length is $\mathcal F_a$,
and the frictional torque per length in (3.1) is

$$
\mathcal F_t = r_c\,\mathcal F_\theta .
$$

Two cross-coupling mechanisms follow from (5.2):

* **Velocity split.** Kinetic friction has fixed magnitude $\mu_k F_N$ but is
  shared between the two directions in proportion to the velocity components:

$$
\mathcal F_a = \mu_k F_N\,\frac{v}{\sqrt{v^2 + \omega^2 r_c^2}},
\qquad
\mathcal F_t = \mu_k F_N\,r_c\,\frac{\omega r_c}{\sqrt{v^2 + \omega^2 r_c^2}} .
$$

  Rotating while tripping therefore reduces axial drag. Conversely, axial
  motion reduces the rotary torque and can change the stick-slip behaviour.

* **Normal-force coupling.** $F_N$ depends on the live tension (5.1), so axial
  loads change the torque along the string. This is the effect that Aarsnes &
  Shor (2018, §5.5) observed as a phase shift once axial motion started.

A sticking contact also breaks out under **combined** load. The friction capacity
$\mu_s F_N$ is shared, so torque and axial pull together can free the string
when neither alone would.

### 5.3 Special case: no axial motion

For $v \equiv 0$, (5.2) reduces exactly to the torsional inclusion of Aarsnes &
Shor (2018, eqs. 4–5), with limits
$\mathcal F_{c,s} = \mu_s F_N r_c$ and $\mathcal F_{c,k} = \mu_k F_N r_c$ and
$f_\text{rat} = \mu_k/\mu_s$:

$$
\mathcal{F}_t(\omega, x) \in
\begin{cases}
\{\,\mathcal{F}_{c,k}(x)\,\}, & \omega > \omega_c,\\[2pt]
[-\mathcal{F}_{c,s}(x),\ \mathcal{F}_{c,s}(x)], & |\omega| \le \omega_c,\\[2pt]
\{-\mathcal{F}_{c,k}(x)\}, & \omega < -\omega_c .
\end{cases}
$$

The published off-bottom results (break-off torque, stick-slip amplitude and
period, §9) therefore carry over unchanged.

### 5.4 Decoupled variant (`along_string_coupling=False`)

For comparison, the axial and torsional inclusions can be applied
independently, each in the 1-D form above, with thresholds $v_c$ and
$\omega_c$. The normal force is then taken from the static hanging string
(§5.1). There is no velocity split and no normal-force coupling, so the only
axial–torsional coupling is at the bit.

---

### 5.5 Axial static friction

The static and kinetic factors $\mu_s > \mu_k$ were identified from
**torsional** start-ups (Aarsnes & Shor 2018). Axial friction factors from
pick-up and slack-off are kinetic values (Johancsik et al. 1984).

**What goes wrong if the torsional contrast is also applied axially.** Applying
the full contrast axially, with the step law and $v_c = 0.01$ m/s, gives
sustained axial stick-slip at 0 rpm. In slack-off at 200 m/h on Well A:
* the bit alternates between rest and 600 m/h;
* the cycle is about 3 s;
* hookload swings by 40 kN.

The mechanism is the classic one for velocity-weakening friction on a long,
compliant string driven at constant speed: a step from $\mu_s$ to $\mu_k$ is
infinitely steep weakening.

**Options evaluated.** All were run on Well A at 0 rpm, using the bit axial
speed after 40 s; the bit is off bottom.

| option | slack-off 200 m/h | slack-off 20 m/h | pick-up 200 m/h | torsional stick-slip kept |
|---|---|---|---|---|
| full contrast (previous default) | 5–658 m/h | 5–28 | −1128 to −2 | yes |
| $v_c$ = 10⁻³ m/s | 154–245 | 4–367 (stick-slip moves to lower speed) | −236 to −165 | yes |
| $v_c$ = 10⁻⁴ m/s | 154–245 | 4–658 | −806 to −3 | yes |
| axial contrast 0.25 | 4–605 | 6–32 | −805 to 0 | yes |
| exponential weakening $\mu_k + (\mu_s-\mu_k)e^{-s}$ | 5–644 | 5–28 | −1267 to 1 | yes |
| block acceleration ≤ 0.5 m/s² | 5–658 | 5–28 | −1115 to −2 | yes |
| raise axial damping $k_a$ to 1–5 /s | 190–206 | – | – | yes, but not physical: steady Couette mud drag is about 10⁻² /s and the oscillatory Stokes layer about 0.1 /s (§5.6) |
| **axial contrast 0** (new default) | **170–225**, settling to ±5 within 100 s | **6–33** | **−233 to −172** | **yes** |

**Default.** Static friction now acts on an **ellipse** in force space:

$$
\left(\frac{\mathcal F_a}{\mu_{s,a}F_N}\right)^2 + \left(\frac{\mathcal F_\theta}{\mu_s F_N}\right)^2 \le 1,
\qquad \mu_{s,a} = \mu_k + c_a(\mu_s - \mu_k),
$$

with $c_a$ = `axial_static_contrast` = 0 by default. A cell that breaks out
slides with the ellipse radius in the direction of motion for that step, then
with $\mu_k$. The torsional behaviour, and the cross-coupling of §5.2, are
unchanged. The full contrast ($c_a = 1$) and the exponential weakening law
remain available for research, for example to study axial breakover after
connections.

The ±10% start-up ringing that remains is the first axial mode (0.42 Hz),
excited by the step change in block speed. It decays with the axial damping,
which §5.6 now derives from the mud and the steel rather than setting by hand.

### 5.6 Distributed damping: fluid and material

The lumped rates $k_t, k_a$ of (3.1)–(3.2) had no physical basis and damped
every mode at the same rate. They are replaced by three mechanisms with
parameters from the mud report and the string. The constant rates remain as
optional calibration terms (`DampingParameters.k_t`, `k_a`, default 0).

**Choice of mechanisms.** Five candidates were considered:
* steady laminar (Couette) annular shear;
* the oscillatory Stokes boundary layer;
* the yield stress of a Bingham or Herschel–Bulkley mud;
* surge/swab piston flow;
* material damping of steel and tool joints.

For 5" drill pipe in an 8½" hole with a 20 cP mud:
* **Couette** drag is $c_0 \approx 0.24$ N s/m², so $k \approx 0.008$ /s.
* **Stokes layer** at the 0.4 Hz first axial mode ($\delta \approx 3.6$ mm) gives $c_1 \approx 4$ N s/m², so $k \approx 0.1$ /s and $\zeta \approx 0.02$–0.03.
* **Surge/swab** through the DP annulus is $\approx 0.06$ N s/m² per metre. It matters only across tight BHA annuli and is left to a future lumped BHA term.
* **Steel:** the measured loss factor of a steel member is $\eta \approx 0.014$, and jointed metal structures are higher (Orban 2011). This gives $\zeta = \eta/2 \approx 0.005$ for the string.

The oscillatory fluid layer therefore dominates, as Moharrami et al. (2024) and
Jansen & van den Steen (1995) assume with lumped coefficients. The
confined-annulus added mass of Chen et al. (1976) and Dolfo et al. (2020)
applies to lateral motion and is not used here.

#### 5.6.1 F1 — fluid: steady Couette drag

Laminar flow between the pipe (radius $r_o$) and the hole or casing (radius
$R$), with viscosity $\mu$ (the Bingham plastic viscosity if given), gives
drag per unit length

$$
f_{a,0} = \frac{2\pi\mu}{\ln(R/r_o)}\, v,
\qquad
m_{t,0} = \frac{4\pi\mu\, r_o^2 R^2}{R^2 - r_o^2}\,\omega. \tag{5.6}
$$

These act on the full velocity and are the only viscous loads in steady
rotation or tripping, so they enter torque-and-drag directly.

#### 5.6.2 F1 — fluid: oscillatory Stokes layer

A wall oscillating tangentially at frequency $\omega_f$ in a viscous fluid
drives a boundary layer of thickness $\delta = \sqrt{2\nu/\omega_f}$, with
$\nu = \mu/\rho_f$ (Stokes' second problem). The in-phase (dissipative) part
of the wall shear stress is

$$
\tau_w = \sqrt{\rho_f\mu\,\omega_f/2}\; u_w .
$$

For $\delta \ll r_i$ and $\delta \ll R - r_o$ (both hold for $f \gtrsim 0.1$
Hz), each wall acts as a flat plate. Summing the inner and outer walls of the
pipe, and integrating the moment $\tau_w r$ for torsion, gives the drag per
unit length

$$
f_{a,1} = \underbrace{\sqrt{\tfrac{\rho_f\mu\omega_a}{2}}\;2\pi(r_o + r_i)}_{c_{a,1}}\; \tilde v,
\qquad
m_{t,1} = \underbrace{\sqrt{\tfrac{\rho_f\mu\omega_t}{2}}\;2\pi(r_o^3 + r_i^3)}_{c_{t,1}}\; \tilde\omega . \tag{5.7}
$$

The mud inside the pipe is not entrained, consistent with the steel-only mass
$\rho A$. The coefficient grows as $\sqrt{\omega_f}$, which is a fractional
(half-order) time derivative. It is evaluated once at a reference frequency:
by default the first (quarter-wave) mode of the string held at the top,

$$
f_{a,1} = \frac{1}{4\sum_j h_j/c_{a,j}}, \qquad f_{t,1} = \frac{1}{4\sum_j h_j/c_{t,j}} . \tag{5.8}
$$

This is the mode that dominates stick-slip and the tripping transients.

$\tilde v, \tilde\omega$ are the **oscillatory** parts of the velocities. A
uniformly moving pipe has no Stokes layer, only (5.6). Each cell therefore
carries a low-pass state

$$
\dot{\bar v} = \omega_c(v - \bar v), \qquad \tilde v = v - \bar v, \qquad \omega_c = \phi\,\omega_{a,1},
$$

with $\phi$ = `high_pass_fraction` = 0.25. At the reference frequency the
filter passes $1/(1+\phi^2) = 94\%$ of the Stokes damping; at steady state it
passes none.

#### 5.6.3 F2 — fluid: yield-stress drag

A Bingham mud with yield point $\tau_y$ resists the start of motion of a gelled
column with a velocity-independent wall stress. Following the
Herschel–Bulkley drag form $F_d = A_d\tau_y + k_f A_d(V/\Delta)^n$
(Lyons & Plisga 2005, as used in the Swansea thesis), the yield part per unit
length,

$$
f_y = \tau_y\,\pi d_o, \tag{5.9}
$$

is added to the Coulomb capacity of §5.2 for both directions:
$\mu F_N \to \mu F_N + f_y$, in the static ellipse and while sliding. It is
applied with the sliding velocity at the contact radius, which is a
simplification when $r_c \ne r_o$. While the mud circulates, it is already
yielded in the annulus, so by default (5.9) applies only with the pumps off
(`yield_drag='pumps_off'`; also `'always'`, `'off'`).

#### 5.6.4 M1 — material: Rayleigh / Kelvin–Voigt

For a linear viscoelastic (Kelvin–Voigt) rod, $\sigma = E(\varepsilon + \beta\dot\varepsilon)$.
The strain rate is $\dot\varepsilon = \partial_x v$, which adds

$$
\rho A\,\partial_t v \mathrel{+}= \partial_x\!\left(\beta EA\,\partial_x v\right),
\qquad
\rho J\,\partial_t \omega \mathrel{+}= \partial_x\!\left(\beta GJ\,\partial_x \omega\right) \tag{5.10}
$$

to the momentum equations, with modal damping ratio $\zeta_n = \beta\omega_n/2$.
Steel and its connections dissipate a nearly frequency-independent fraction
of energy per cycle (hysteretic damping, $\zeta = \eta/2$). That law is not
causal in the time domain, so it is approximated by a Rayleigh pair, adding a
mass-proportional rate $\alpha$ that acts on $\tilde v, \tilde\omega$ like the
Stokes term:

$$
\zeta(\omega) = \frac{\alpha}{2\omega} + \frac{\beta\omega}{2},
\qquad
\alpha = \frac{\eta\,\omega_1\omega_2}{\omega_1+\omega_2},
\quad
\beta = \frac{\eta}{\omega_1+\omega_2}, \tag{5.11}
$$

which gives $\zeta = \eta/2$ exactly at the first mode $\omega_1$ and at
$\omega_2 = 2\pi f_\max$, and $\zeta < \eta/2$ in between. Moharrami et al.
(2024) fit Rayleigh damping over modes 1–10 in the same way. A Kelvin–Voigt
term tuned only at the first mode ($\beta = \eta/\omega_1$) was rejected: its
$\zeta \propto f$ over-damps 10 Hz content by roughly 40×. The default is
$\eta = 0.01$.

#### 5.6.5 Discretisation

Each step, after the Godunov flux update (§8.1):

1. **Linear fluid and $\alpha$ terms, implicitly per cell:**

$$
\omega^{**} = \frac{\omega^* + \Delta t\,k_{t,1}\bar\omega}{1 + \Delta t(k_{t,0} + k_{t,1})},
\qquad
k_{t,0} = \frac{c_{t,0}}{\rho J} + k_t,\quad k_{t,1} = \frac{c_{t,1}}{\rho J} + \alpha_t,
$$

   and likewise for $v$.

2. **Kelvin–Voigt, implicitly, by one tridiagonal (Thomas) solve per
   direction:**

$$
\frac{m_j h_j}{\Delta t}(x_j' - x_j) = \kappa_{j+\frac12}(x'_{j+1} - x'_j) - \kappa_{j-\frac12}(x'_j - x'_{j-1}),
\qquad
\kappa_{j+\frac12} = \beta\,\frac{\overline{EA}}{\bar h},
$$

   with $\kappa = 0$ at the top-drive and bit faces. An explicit step would be
   unstable: on a 10 m grid the diffusion limit $h^2/(2\beta c^2)$ is far below
   the CFL step. The implicit step:
   * conserves momentum $\sum m_j h_j x_j$;
   * dissipates $\sum\kappa\,(\Delta x')^2 \ge 0$;
   * costs $O(N)$.

3. **Friction** (§8.3), with capacity $\mu F_N + f_y$.

4. **Low-pass update:** $\bar x \leftarrow \bar x + (1 - e^{-\omega_c\Delta t})(x - \bar x)$.

The low-pass states are part of the simulation state, so pausing stays exact
(§8.5).

#### 5.6.6 Damping ratios, numerical dissipation and the 10 Hz limit

For a mode of frequency $f$ ($\omega = 2\pi f$), the model gives

$$
\zeta_\text{fluid} = \frac{k_0 + k_{1,S}H}{2\omega},\quad
\zeta_\text{mat} = \frac{\alpha H}{2\omega} + \frac{\beta\omega}{2},\quad
\zeta_\text{num} = \frac{\pi f h (1-\nu)}{c},\quad H = \frac{1}{1 + (\phi f_1/f)^2}, \tag{5.12}
$$

where $\zeta_\text{num}$ is the amplitude loss of first-order upwinding.
Upwinding has modified-equation diffusion $c h(1-\nu)/2$, with $\nu$ the
Courant number of that wave. `Simulation.damping_summary(f)` tabulates these
per component.

For Well A drill pipe ($f_{a,1}$ = 0.51 Hz, $f_{t,1}$ = 0.32 Hz), the first
modes have $\zeta_\text{fluid}$ = 0.027 (axial) and 0.037 (torsional), with
$\zeta_\text{mat}$ = 0.0047. Previously the default $k$ = 0.1 /s gave 0.016
and 0.025.

**Limitations (stated validity: $f \le 10$ Hz).**

* The Stokes coefficient is fixed at $f_1$. Above it, fluid damping is
  under-predicted by $\sqrt{f/f_1}$: 4–6× at 10 Hz, an absolute shortfall of
  $\zeta \approx 0.005$.
* Frequency-dependent Stokes damping would need a fractional-derivative
  (Prony / auxiliary-state) form; it is deferred.
* Material damping follows (5.11) within the band and over-damps above $f_\max$.
* Numerical dissipation grows as $f\,h$. With $h$ = 10 m, the torsional
  $\zeta_\text{num}$ reaches 0.04 at 10 Hz and exceeds the physical damping
  above about 2 Hz; the axial $\zeta_\text{num}$ stays below 0.003. For
  torsional content up to 10 Hz, use $h \lesssim 2.5$ m ($\zeta_\text{num}
  \approx 0.01$), at about 4× the run time.
* Herschel–Bulkley shear thinning in the Stokes layer, which depends on
  amplitude, and the gel-strength build-up with time are not modelled.

**Verification** (`tests/test_damping.py`, `docs/derivation/damping_verification.py`).
Free decay of the first mode of 3000 m of 5" DP gives measured
vs (5.12)-predicted damping ratios:

| model | axial measured / predicted | torsional measured / predicted |
|---|---|---|
| lumped $k$ = 0.1 /s | 0.0187 / 0.0187 | 0.0310 / 0.0311 |
| F1 fluid | 0.0290 / 0.0293 | 0.0414 / 0.0415 |
| M1 material, $\eta$ = 0.02 | 0.0095 / 0.0096 | 0.0105 / 0.0105 |
| F1 + M1, $\eta$ = 0.01 | 0.0335 / 0.0340 | 0.0461 / 0.0462 |

For Well A tripping in at 200 m/h and 0 rpm, the bit-speed range after 30 s
narrows from 154–248 m/h (lumped $k$ = 0.1 /s) to 172–230 m/h (F1 + M1).

Figures: `docs/figures/damping_free_decay.png`, `damping_vs_frequency.png`,
`damping_trip_well_a.png`.


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

**The hook cannot push.** If the prescribed block speed would make
$\text{HKLD} < 0$, the block stacks out. Hookload is then held at zero, so
$q_d(0) = -q_s(0)$, and the top of the string moves freely with
$v(0) = \beta_{a,1} + q_d(0)/Z_a$. For example, with the bit on bottom at
0 rpm, WOB levels off at the buoyed string weight less friction, instead of
growing without limit.

**Set-point ramps.** Every input approaches a new value linearly over its
ramp time (`ModelParameters.ramps`). The ramp is evaluated inside the time
stepping from a table of (start value, start time, ramp time, target). This
holds for the speed set point, the prescribed top-drive speed, the block
velocity, the WOB set point, the pump rate and the inlet temperature. A set
point that takes over from another mode starts from the present physical
value. See [`hydraulics.md`](hydraulics.md) §6.

**Hydraulics.** Standpipe pressure, bottomhole pressure and ECD come from
[`hydraulics.md`](hydraulics.md): chapter 7 steady losses, plus single-phase
drift-flux transients after Wang et al. (2023).

An optional block-acceleration limit (`DrawworksParameters.block_acceleration`)
rate-limits changes in $v_\text{blk}$.

### 6.3 Drill bit: coupling at the bit

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

**Rock strength.** The drilling strength $S$ is either a constant or, with
`strength_source='formation'`, $S = \zeta\,\text{CCS}_\text{eff}/\eta$. In the
second case $S$ is computed from the formation, pore pressure, bottomhole
pressure and temperatures at the bit. See
[`rock-strength-and-temperature.md`](rock-strength-and-temperature.md), which
also derives the wellbore temperature model (Wang et al. 2023,
SPE-214836-MS).

**Bit–rock law.** Three reduced models are available
(`BitParameters.model`). Their derivation from the Bit Model Library is in
[`bit-models.md`](bit-models.md). With $a = D/2$ and the drilling strength
$S = \zeta\varepsilon$:

| model | parameters | torque | penetration |
|---|---|---|---|
| `coulomb` | $\mu_s, \mu_k, S$ | sticking $\lvert\text{TOB}\rvert \le \mu_s\,\text{WOB}\,D/3$; sliding $\mu_k\,\text{WOB}\,D/3$ | $\text{ROP} = \omega\,\text{WOB}/(2\pi a S)$ |
| `cutting` | $S, \zeta, W_f$ | $\text{WOB}\,a/(2\zeta)$, both sticking and sliding | $\text{ROP} = \omega\,(\text{WOB}-W_f)^+/(2\pi a S)$ |
| `regenerative` | $S, \zeta$ + blades $n$ | $\text{WOB}\,a/(2\zeta)$ | from the hole-bottom profile (below) |

**Quasi-static models (`coulomb`, `cutting`).** Sticking applies when
$|\omega_b| \le \omega_c$: the bit stays at rest, and $\text{TOB} = Z_t\alpha_t$
from (6.3), as long as that is within the sticking limit. When sliding,
$\text{TOB}$ is limited to $|Z_t\alpha_t|$, so the bit is never driven backwards
within a step. Combining the penetration law with (6.3): if
$Z_a\alpha_a \le W_0$ (with $W_0 = W_f$, or 0 for Coulomb) or $\omega_b = 0$, then
$v_b = 0$ and $\text{WOB} = Z_a\alpha_a$. Otherwise

$$
\text{WOB} = \frac{\alpha_a + k_r W_0 |\omega_b|}{1/Z_a + k_r|\omega_b|},
\qquad k_r = \frac{1}{2\pi a S}.
$$

These are solved by fixed-point iteration on $\omega_b$.

**Regenerative model.**
* **Profile.** The hole bottom is stored as a profile $s(\theta)$ over 720
  angular bins.
* **Depth of cut.** Blade $j$, at angle $\varphi_b + 2\pi j/n$, cuts to depth
  $d_j = (u_b - s(\theta_j))^+$, measured against the surface left by the
  previous blade:

$$
\text{WOB} = S a \sum_j d_j, \qquad \text{TOB} = \frac{a}{2\zeta}\,\text{WOB}.
$$

* **Axial contact.** The contact is solved implicitly for the new bit depth,
  $u_b^{k+1} = u_b^k + \Delta t\,(\alpha_a - \text{WOB}/Z_a)$, by iterating on
  the set of blades in contact. If the bit rises above the profile, WOB is
  zero (lift-off).
* **Sweep.** After each step the blades sweep from $\varphi_b^k$ to
  $\varphi_b^{k+1}$, and the bins they leave are set to the bit depth. The
  depth of cut therefore carries the delay of one blade passage. That delay is
  the regenerative effect, a coupling between axial and torsional motion that
  can self-excite.
* **Hole depth.** The hole depth is the deepest point of the profile.

**Coupling mechanisms at the bit.** These act in addition to the along-string
friction coupling of §5.2.
* Torsion depends on axial: TOB is proportional to WOB in every model, so
  axial waves that change WOB modulate the torsional boundary and can trigger
  or suppress stick-slip.
* Axial depends on torsion: penetration needs rotation. A stalled bit stops
  penetrating, WOB rises, and the sticking limit rises with it.
* In the regenerative model the depth of cut also remembers the previous
  blade passage. That delay can destabilise steady drilling, producing bit
  bounce and stick-slip.

---

## 7. Top-drive controllers (Aarsnes, di Meglio & Shor 2018)

All three industrial controllers fit one block diagram. A ramped set
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

1. **Damping, implicitly:** the fluid and Rayleigh-$\alpha$ rates per cell,
   then the Kelvin–Voigt tridiagonal solve (§5.6.5). With only the
   calibration rates this reduces to $\omega^{**} = \omega^*/(1 + k_t\Delta t)$.
2. **Friction, implicitly (coupled model).** Let
   $\mathbf V^{*} = (v^{**},\ \omega^{**} r_c)$ be the predicted sliding
   velocity and $M = \operatorname{diag}(m_a, m_\theta)$ the mass per length in
   each direction, with $m_a = \rho A$ and $m_\theta = \rho J / r_c^2$. Coulomb
   friction (5.2) is applied with a backward-Euler step:

$$
\mathbf V^{k+1} = \mathbf V^{*} - \Delta t\, M^{-1}\,\mu F_N\,\frac{\mathbf V^{k+1}}{\|\mathbf V^{k+1}\|},
\qquad
c \equiv \Delta t\,\mu F_N,
$$

   with $\mu = \mu_s$ if the cell was sticking at step $k$ and $\mu = \mu_k$
   otherwise. Two cases follow:

   * **Stick.** If $\|M\mathbf V^{*}\| \le c$, the friction impulse can absorb
     the cell's momentum, and $\mathbf V^{k+1} = 0$ exactly. The impulse needed
     to stop the cell, $M\mathbf V^*/\Delta t$, is the 2-D form of
     $\tilde{\mathcal F}$ in Aarsnes & Shor (2018, eqs. 20–21).
   * **Slide.** Otherwise each component satisfies

$$
V_i^{k+1} = \frac{m_i V_i^{*}\, s}{m_i s + c},
\qquad
g(s) = \sum_i \left(\frac{m_i V_i^{*}}{m_i s + c}\right)^2 - 1 = 0,
$$

     where $s = \|\mathbf V^{k+1}\|$. The function $g$ is convex and
     decreasing, with $g(0) > 0 \ge g(\|\mathbf V^*\|)$, so Newton's method
     from $s = \|\mathbf V^*\|$ converges monotonically after the first step.
     For equal masses this reduces to $s = \|\mathbf V^*\| - c/m$.

   The implicit form guarantees three properties:
   * friction opposes the **new** sliding direction exactly;
   * friction never reverses the motion within a step;
   * sticking is exact.

   An explicit update with each component clipped separately would distort the
   friction direction when one component is small.

   **Decoupled variant.** Each direction is projected separately:

$$
\mathcal F = \operatorname{sign}(\tilde{\mathcal F})\,\min\!\left(|\tilde{\mathcal F}|,\ \mathcal F_c\right),
\qquad
\tilde{\mathcal F}_t = \rho J\,\omega^{**}/\Delta t,
\qquad
\tilde{\mathcal F}_a = \rho A\,v^{**}/\Delta t,
$$

   which is the published scheme (eq. 22). For $v \equiv 0$ the coupled update
   reduces to exactly this.

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
  friction (5.2) gives the soft-string torque and drag equations, integrated from
  the bit up:

$$
\frac{\mathrm{d}T}{\mathrm{d}x} = \mathcal F_a - w_b\cos\theta,
\qquad
\frac{\mathrm{d}\tau}{\mathrm{d}x} = -\mathcal F_t,
\qquad
T(L) = -\text{WOB},\ \tau(L) = \text{TOB}.
$$

  This is `statics.torque_drag`, with friction split by the velocity direction
  exactly as in (5.2). The coupled dynamic model converges to it in steady
  combined motion. The decoupled variant is compared with it using the static
  normal force (`normal_force=sim.fn0`).

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
  * Coulomb: breakout at $\mu_s\,\text{WOB}\,D/3$; while sliding,
    $\text{TOB}/\text{WOB} = \mu_k D/3$.
  * Cutting: steady drilling at block speed $V$ and bit speed $\Omega$ gives
    $\text{WOB} = W_f + 2\pi a S V/\Omega$ and
    $\text{TOB} = a\,\text{WOB}/(2\zeta)$.
  * Regenerative: the same steady state as cutting with $W_f = 0$ when it is
    stable.

---

## 10. Verification summary (`tests/`)

| check | result |
|---|---|
| torsional / axial wave arrival at $L/c$ | within 1% |
| fundamental period $4L/c_t$ | within 0.2%; amplitude loss < 1% over 8 round trips |
| hanging string at rest | exactly at rest (well balanced) |
| steady rotation, pick-up, slack-off vs torque & drag (both friction models) | within 2% |
| tripping at 0 rpm (±200, 20 m/h): no axial stick-slip; the full axial contrast reproduces it | yes |
| on bottom at 0 rpm: hookload ≥ 0, WOB levels off | yes |
| coupled friction vs velocity-split torque & drag (rotate while picking up / slacking off) | within 2% |
| rotation reduces pick-up drag (coupled); has no effect (decoupled) | yes / within 0.2% |
| combined axial + torsional breakout, implicit friction exact | to $10^{-10}$ |
| off-bottom stick-slip only when $\mu_s > \mu_k$ | yes |
| ZTorque suppresses the stiff-controller limit cycle | yes |
| bit breakout torque, TOB/WOB slope, drill-off at constant block | exact / 1% |
| pause–resume, external Python PI vs internal PI | exact / 2% |
| free decay of first axial / torsional mode vs (5.12): lumped, F1, M1, F1+M1 | within 2% (8 cases) |
| steady rotation with F1+M1: surface torque = Couette (5.6) only | within 2% |
| yield drag (5.9) with pumps off; absent with pumps on | within 5% |
| Kelvin–Voigt solve conserves momentum, dissipates energy | to $10^{-12}$ |

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
| $k_t, k_a$ | additional (calibration) viscous damping rates | 1/s |
| $c_{\cdot,0}, c_{\cdot,1}$ | Couette, Stokes-layer damping coefficients per length | N s/m², N m s |
| $\delta = \sqrt{2\nu/\omega}$ | Stokes layer thickness | m |
| $\eta;\ \alpha, \beta$ | material loss factor; Rayleigh coefficients | –; 1/s, s |
| $\tau_y, f_y$ | mud yield point, yield drag per length | Pa, N/m |
| $f_1, f_\max$ | first-mode and maximum modelled frequency | Hz |
| $I_{TD}, \tau_m, \tau_0$ | top-drive inertia, motor torque, pipe torque at top | kg m², N m, N m |
| $a = D/2$ | bit radius | m |
| $S = \zeta\varepsilon$ | drilling strength (axial cutting force per unit cut area) | Pa |
| $\zeta, \varepsilon$ | normal/tangential cutting force ratio, intrinsic specific energy | –, Pa |
| $\mu_s, \mu_k$ (bit) | Coulomb bit static/kinetic friction (lever $D/3$) | – |
| $W_f$ | wear-flat contact WOB (cutting model) | N |
| $n$ | blades (regenerative model) | – |

## References

* Aarsnes, U.J.F. and Shor, R.J., 2018. Torsional vibrations with bit off bottom: Modeling, characterization and field data validation. *J. Petroleum Science and Engineering*, 163, 712–721.
* Aarsnes, U.J.F., di Meglio, F. and Shor, R.J., 2018. Benchmarking of industrial stick-slip mitigation controllers. *IFAC-PapersOnLine*, 51(8), 233–238.
* Auriol, J., Boussaada, I., Shor, R.J., Mounier, H. and Niculescu, S.-I., 2022. Comparing advanced control strategies to eliminate stick-slip oscillations in drillstrings. *IEEE Access*, 10, 10949–10969.
* Chen, S.S., Wambsganss, M.W. and Jendrzejczyk, J.A., 1976. Added mass and damping of a vibrating rod in confined viscous fluids. *J. Applied Mechanics*, 43(2), 325–329.
* Dolfo, G., Vigué, J. and Lhuillier, D., 2020. Motion of a cylinder oscillating in a confined viscous fluid. arXiv:2011.12000.
* Jansen, J.D. and van den Steen, L., 1995. Active damping of self-excited torsional vibrations in oil well drillstrings. *J. Sound and Vibration*, 179(4), 647–668.
* Lyons, W.C. and Plisga, G.J., 2005. *Standard Handbook of Petroleum and Natural Gas Engineering*, 2nd ed. Gulf Professional Publishing.
* Moharrami, M.J., Shiri, H. and Martins, C.A., 2024. Vibration (MDPI), 7(4), 56. https://www.mdpi.com/2571-631X/7/4/56
* Orban, F., 2011. Damping of materials and members in structures. *J. Physics: Conference Series*, 268, 012022.
* Stokes, G.G., 1851. On the effect of the internal friction of fluids on the motion of pendulums. *Trans. Cambridge Phil. Soc.*, 9, 8–106.
* Johancsik, C.A., Friesen, D.B. and Dawson, R., 1984. Torque and drag in directional wells — prediction and measurement. *JPT*, 36, 987–992.
* Shor, R.J. and Aarsnes, U.J.F. Coupled axial-torsional vibrations in a drillstring: modelling and validation. Draft manuscript.
