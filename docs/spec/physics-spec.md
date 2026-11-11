# DrillPyze physics specification

**Spec revision: 1** (matches the Python model on branch `browser-lite`, November 2026)

This document specifies, at implementation level, every piece of physics and
numerics in the Python model (`src/DrillPyze`). It is precise enough to build
an independent implementation; the browser version in `lite/` is built from
it. The derivations behind the equations are in `docs/derivation/` and are
not repeated here.

> **Maintenance rule.** Any change to the physics or numerics of the Python
> model must update this specification in the same commit, increase the
> **spec revision** above and refresh `docs/spec/manifest.json` with
> `python docs/spec/update_manifest.py`. `tests/test_spec_sync.py` fails
> when a physics source file changes without a manifest update, and when the
> browser version (`lite/engine/spec.js`) implements an older revision.
> A change that does not alter the physics (comments, refactoring) only needs
> the manifest refresh.

## Contents

0. [Conventions and coverage](#0-conventions-and-coverage)
1. [Case description and defaults](#1-case-description-and-defaults)
2. [Geometry](#2-geometry)
3. [String model: state, properties and time step](#3-string-model-state-properties-and-time-step)
4. [Parameter set-up](#4-parameter-set-up)
5. [The time step (kernel)](#5-the-time-step-kernel)
6. [Inputs, ramps and modes](#6-inputs-ramps-and-modes)
7. [Slow tick](#7-slow-tick)
8. [Hydraulics](#8-hydraulics)
9. [Rock strength at the bit](#9-rock-strength-at-the-bit)
10. [Wellbore temperature (full model only)](#10-wellbore-temperature-full-model-only)
11. [Steady torque and drag (reference profiles)](#11-steady-torque-and-drag-reference-profiles)
12. [Outputs](#12-outputs)
13. [Conformance](#13-conformance)

---

## 0. Conventions and coverage

**Units.** SI throughout: m, s, kg, N, Pa, rad, rad/s; temperatures in °C.
Rotary speeds are converted with 1 rpm = $`2\pi/60`$ rad/s. Constants:
$`g = 9.80665`$ m/s², steel density $`\rho_s = 7850`$ kg/m³,
$`E_s = 207`$ GPa, $`G_s = 79.3`$ GPa.

**Signs.**
- Measured depth (MD) and axial velocity are positive downhole.
- The axial force variable $`q`$ is compression positive; effective tension is $`T = -q`$.
- Rotation and torque are positive in the drilling direction.
- Hookload is the effective tension at the top of the string, without block weight.
- WOB is compression at the bit, at least 0.

**Indexing.** String cells $`j = 0 \dots n-1`$ run from surface to bit; faces
$`j = 0 \dots n`$, where face $`j`$ is the top of cell $`j`$ and face $`n`$ is the bit.

**Coverage of the browser version (lite).**

| Section | Lite |
|---|---|
| 2 Geometry, 3–7 string model, kernel, inputs, slow tick | yes |
| 8 Hydraulics: steady and transient (constant mud temperature) | yes |
| 9 Rock strength (P1 confinement, T1 weakening; no T2 thermal stress) | yes |
| 10 Wellbore temperature | no (temperature = reference temperature) |
| 11 Steady torque and drag (reference curves) | yes |

The lite version is the Python model with `parameters.thermal.enabled = False`.

---

## 1. Case description and defaults

A case is a JSON object (see `src/DrillPyze/config.py`):

| Key | Content |
|---|---|
| `well.survey` | stations `{md, inc, azi}` (degrees) |
| `well.casing_shoe` | MD of the deepest shoe (m) |
| `well.hole_diameter` | bit diameter $`D`$ (m) |
| `fluid` | density $`\rho_m`$, plastic viscosity PV (or null), yield point $`\tau_y`$ (or null), viscosity $`\mu`$, heat capacity, conductivity, thermal expansion $`\beta_T`$, sound speed $`c_m`$, reference temperature $`T_\text{ref}`$ |
| `formation` | layers `{top_md, ucs, friction_angle, porosity}` and the scalars of section 9 |
| `drillstring` | components bit-first `{name, length, od, id, linear_mass, contact_od}`; the last row is stretched to surface |
| `parameters` | friction, damping, bit, top_drive, controller, autodriller, drawworks, ramps, hydraulics, thermal, `cfl` |
| `simulation` | `bit_depth`, `hole_depth`, `dx`, `record_rate`, `start_time`, `block_height` |

**Defaults** (dataclasses in `parameters.py`, `wellbore.py`, `formation.py`):

| Group | Parameter | Default |
|---|---|---|
| friction | $`\mu_s, \mu_k`$ cased / open | 0.30, 0.20 / 0.35, 0.25 |
| | $`\omega_c`$, $`v_c`$ | 1 rad/s, 0.01 m/s |
| | `along_string_coupling` | true |
| | `axial_static_contrast` $`c_{ax}`$ | 0 |
| | `weakening` | `step` (or `exponential`) |
| damping | $`k_t, k_a`$ | 0, 0 (1/s) |
| | `fluid` | true |
| | `loss_factor` $`\eta`$ | 0.01 |
| | reference frequencies | first modes (null) |
| | `high_pass_fraction` $`f_c`$ | 0.25 |
| | `yield_drag` | `pumps_off` |
| | `max_frequency` | 10 Hz |
| bit | `model` | `coulomb` (`cutting`, `regenerative`) |
| | $`\mu_{bs}, \mu_{bk}`$ | 0.6, 0.4 |
| | drilling strength $`S`$ | 50 MPa |
| | $`\zeta`$ | 0.8 |
| | contact WOB $`W_f`$ | 0 |
| | blades | 5 |
| | radius | null (= $`D/3`$) |
| | `strength_source` | `constant` (or `formation`) |
| | cutting efficiency | 1 |
| | nozzle area $`A_n`$ | 3.5628e-4 m² |
| | $`C_d`$ | 0.95 |
| top drive | inertia $`I`$ | 2000 kg m² |
| | max torque | 0 (unlimited) |
| autodriller | $`K_p, K_i, v_{max}`$ | 1e-7 m/(s N), 1e-8 m/(s² N), 60 m/h |
| drawworks | block acceleration | 0 (unlimited) |
| ramps (s) | ω_sp, ω_td, τ_m, v_block, WOB_sp, flow, T_in | 0, 1, 0, 2, 0, 3, 3 |
| hydraulics | `enabled`, dx, cut factor, CFL, Picard iterations | true, 30 m, 1.5, 0.8, 2 |
| | back pressure, surface loss $`c`$, exponent $`m`$ | 0, 0, 1.8 |
| | clinging constant $`K_c`$ | 0.45 |
| thermal | flow rate, inlet temperature | 0.0347 m³/s, 30 °C (+ section 10) |
| model | `cfl` | 0.95 |
| | `thermal_update_period` | 2 (slow ticks) |
| fluid | $`\rho_m`$, PV, $`\tau_y`$ | 1000, null, null |
| | $`\mu`$ | 0.020 Pa s |
| | $`\beta_T`$, $`c_m`$, $`T_\text{ref}`$ | 1e-4 /K, 1500 m/s, 20 °C |

**Well A** (`config.well_a()`, the default case):
- **Survey:** vertical to 1500 m, building to 60° at 2100 m, then tangent to 3000 m (azimuth 0).
- **Wellbore:** shoe at 1000 m; hole 8.5 in.
- **Mud:** $`\rho_m = 1200`$ kg/m³.
- **Formation:** the default `Formation()`.
- **Drillstring:**
  - 230 m of collars, 5.75 in × 2.8125 in;
  - drill pipe, 5 in × 4.276 in (no tool joints).
- **Friction:** cased 0.40 / 0.22, open 0.45 / 0.25.
- **Bit:** `cutting` with `strength_source = formation`.
- **Controller:** `stiff`.
- **Circulation:** flow 2000 L/min (0.03333 m³/s), inlet 30 °C.
- **Run settings:**
  - bit at 2500 m, hole at 2505 m;
  - dx 10 m;
  - record rate 50 Hz.

---

## 2. Geometry

### 2.1 Survey (minimum curvature)

Stations are sorted by MD and duplicates (ΔMD ≤ 1e-9) dropped. If the first
station is below 0, a vertical station is prepended at MD 0 with the first
azimuth. Angles are converted to radians. The unit tangent (north, east, down) is

```math
\mathbf t(I, A) = (\sin I \cos A,\; \sin I \sin A,\; \cos I).
```

For each segment $`k`$, the dogleg is
$`\beta_k = \operatorname{atan2}(|\mathbf t_k \times \mathbf t_{k+1}|, \mathbf t_k\cdot\mathbf t_{k+1})`$
(dot product clipped to [−1, 1]).

**Tangent at MD $`s`$.**
- Segment: $`k`$ = index of the last station with MD ≤ s, clipped to $`[0, N-2]`$.
- Fraction: $`f = \operatorname{clip}((s - s_k)/(s_{k+1}-s_k), 0, 1)`$. Depths beyond the last station therefore use the last tangent.
- If $`\beta_k < 10^{-9}`$, the tangent is $`\mathbf t_k`$. Otherwise it is the spherical interpolation, normalised:

```math
\mathbf t(s) = \frac{\sin((1-f)\beta)}{\sin\beta}\mathbf t_k + \frac{\sin(f\beta)}{\sin\beta}\mathbf t_{k+1}.
```

**Positions** (north, east, TVD) at arbitrary MD are interpolated linearly in a table.
- The table grid is the union of the station MDs and 6000 evenly spaced points on $`[0, 1.25 L + 100]`$, where $`L`$ is the larger of the last station MD and the deepest MD requested.
- The table is built by trapezoidal integration of $`\mathbf t(s)`$ from position 0 at MD 0.
- It is rebuilt when a deeper MD is requested.
- TVD is the third component.

### 2.2 Wellbore

A point at MD $`s`$ is **cased** if $`s \le`$ `casing_shoe`. The wall diameter
is the casing ID where cased (if `casing_id` is given; it is not in the case
format, so the wall is always $`D`$) and $`D`$ in open hole.

### 2.3 Drillstring and grid

Each component has:
- body area $`A = \frac{\pi}{4}(OD^2 - ID^2)`$;
- polar moment $`J = \frac{\pi}{32}(OD^4 - ID^4)`$;
- mass per length $`m' =`$ `linear_mass` if given, else $`\rho_s A`$;
- effective density $`\rho = m'/A`$;
- contact radius $`r_c = \tfrac12`$ `contact_od` (default OD);
- $`E = E_s`$, $`G = G_s`$.

Discretising with the bit at MD $`s_b`$ and maximum cell length $`\Delta x`$:
- The top component's length is $`s_b - L_{BHA}`$, where $`L_{BHA}`$ is the sum of the other components. It must be positive.
- In surface-down order, each component of length $`L`$ is split into $`n_c = \max(1, \lceil L/\Delta x - 10^{-9}\rceil)`$ equal cells of length $`h = L/n_c`$.
- Cell tops $`x_j`$ are cumulative from 0.
- Each cell stores $`h_j, OD_j, ID_j, r_{c,j}, A_j, J_j, \rho_j, E_j, G_j`$ and its component index.
- BHA cells are those not in the top component.

### 2.4 Cell geometry, static load and friction factors (refresh)

The grid moves with the string. With the string shifted by
$`\delta = s_{bit} - s_{bit,0}`$ from its initial position, cell tops are
$`x_j + \delta`$.

On a **refresh**, all of the following are re-evaluated:
- **When:** at initialisation, at every slow tick (section 7) once the bit has moved at least 1 m (`refresh_distance`) since the last refresh, and when parameters are updated.
- **Edges and tangents:** edges $`e_0 \dots e_n`$ are clipped at ≥ 0. Edge tangents are $`\mathbf t(e_j)`$ and centre tangents $`\mathbf t_c = \mathbf t(\max(e_j + h_j/2, 0))`$.
- **Curvature vector:** $`\boldsymbol\kappa_j = (\mathbf t(e_{j+1}) - \mathbf t(e_j))/h_j`$.
- **Normal component of the downward unit vector:** $`\mathbf z_{\perp,j} = \hat{\mathbf z} - t_{c,z}\,\mathbf t_c`$, where $`\hat{\mathbf z} = (0,0,1)`$.
- **Inclination cosine:** $`t_{z,j} = t_{c,z}`$. Cell-centre MD is $`e_j + h_j/2`$.
- **Friction factors per cell:** $`\mu_{s,j}, \mu_{k,j}`$ are the cased values where the cell centre is cased, else the open-hole values.
- **Buoyed weight per length** (fixed at initialisation): $`w_{b,j} = \rho_j A_j g (1 - \rho_m/\rho_s)`$.
- **Static compression** of the frictionless hanging string, at the edges and the cell centres:

  ```math
  q^s_{e,j} = -\sum_{i \ge j} w_{b,i}\, t_{z,i}\, h_i \quad (q^s_{e,n} = 0),\qquad
  q^s_{c,j} = \tfrac12(q^s_{e,j} + q^s_{e,j+1}),\qquad q^s_{top} = q^s_{e,0}.
  ```

  The total axial force is $`q^s_c + q^d`$. On a refresh after the first, the
  dynamic part absorbs the change so that the total is continuous:
  $`q^d \mathrel{+}= q^s_{c,\text{old}} - q^s_{c,\text{new}}`$. The kernel
  parameter $`q_{s,top}`$ is updated.
- **Static normal force** (decoupled friction only): $`F_{N0,j} = | -q^s_{c,j}\boldsymbol\kappa_j + w_{b,j}\mathbf z_{\perp,j}|`$.

The **soft-string normal force per length** for tension $`T`$ is

```math
F_N = \left| T\,\boldsymbol\kappa + w_b\,\mathbf z_\perp \right|.
```

---

## 3. String model: state, properties and time step

**Per-cell wave properties.** Each cell has:
- stiffnesses $`K_t = GJ`$ and $`K_a = EA`$;
- inertias per length $`m_t = \rho J`$ and $`m_a = \rho A`$;
- impedances $`Z_t = J\sqrt{G\rho}`$ and $`Z_a = A\sqrt{E\rho}`$;
- wave speeds $`c_a = \sqrt{E/\rho}`$ and $`c_t = \sqrt{G/\rho}`$.

**Time step.** $`\Delta t = \text{cfl}\cdot\min_j(h_j / c_{a,j})`$, fixed for the whole run.

**Field state per cell:**
- angular velocity $`\omega_j`$;
- torque $`\tau_j`$;
- axial velocity $`v_j`$;
- dynamic compression $`q^d_j`$;
- low-passed velocities $`\bar\omega_j, \bar v_j`$ (for damping).

All start at 0.

**Scalar state.**

| Symbol | Meaning |
|---|---|
| $`t`$ | time |
| $`\omega_{td}`$ | top-drive speed |
| $`\tau_m`$ | motor torque |
| $`\tau_0`$ | pipe torque at the top |
| $`I_{PI}`$ | PI integrator |
| $`y_{lp}`$, $`y_{hp}`$ | torque filter states |
| $`\omega_f`$ | filtered top-drive speed |
| $`\omega_{sp}`$ | active set point |
| $`s_{bit}`$ | bit depth |
| $`s_{hole}`$ | hole depth |
| $`\omega_b`$ | bit speed |
| TOB | torque on bit |
| $`v_b`$ | bit axial velocity |
| WOB | weight on bit |
| $`v_{blk}`$ | block velocity |
| $`H`$ | hookload |
| $`I_{AD}`$ | autodriller integrator |
| $`x_{blk}`$ | block travel |
| stuck | bit sticking flag |
| $`\phi_b`$ | bit angle |
| $`W = \int|\omega_b|dt`$ | accumulated bit rotation |

**Initial state:**
- $`s_{bit}`$ is the bit depth; $`s_{hole} = \max(\text{hole depth}, s_{bit})`$.
- $`H = -q^s_{top}`$, and the off-bottom reference hookload is $`H_{ref} = -q^s_{top}`$.
- Everything else is 0.
- The regenerative hole-bottom profile is 720 angular bins, all at $`s_{hole}`$.

---

## 4. Parameter set-up

### 4.1 Controller

The pipe impedance is $`\zeta_p = Z_{t,0}`$ (the top cell), and the first modes are

```math
f_{a1} = \frac{1}{4\sum_j h_j/c_{a,j}},\qquad f_{t1} = \frac{1}{4\sum_j h_j/c_{t,j}}.
```

Controller gains (Aarsnes, di Meglio & Shor 2018), selected by `controller.type`:

| Type | $`K_p`$ | $`K_i`$ | $`Z`$ | filters $`(t_{hp}, t_{lp}, t_{en})`$ |
|---|---|---|---|---|
| `stiff` | $`100\zeta_p`$ | $`5I`$ | 0 | (10, 0.005, 0.001), not configurable |
| `soft_torque` | $`4\zeta_p`$ | $`(2\pi f_c)^2 I`$ | 0 | (10, 0.005, 0.001) |
| `z_torque` | $`100\zeta_p`$ | $`5I`$ | $`1/\zeta_p`$ | from the case (defaults 10, 0.001, 0.001) |
| `custom` | case `kp` (default $`100\zeta_p`$) | case `ki` (default $`5I`$) | case `impedance` (0) | case (10, 0.005, 0.001) |

For `soft_torque`, $`f_c`$ is the case value, else $`f_{t1}`$ (else 0.2 Hz).
For the types other than `custom`, a non-empty case `kp`, `ki` or
`impedance` overrides the factory value.

### 4.2 Bit

$`a = D/2`$. With the drilling strength $`S`$ (section 9 when `strength_source = formation`):

| Model | sticking lever $`\ell_s`$ | sliding lever $`\ell_k`$ | $`k_r`$ | $`W_0`$ |
|---|---|---|---|---|
| `coulomb` | $`\mu_{bs} r`$ | $`\mu_{bk} r`$, $`r`$ = radius or $`D/3`$ | $`1/(2\pi a S)`$ | 0 |
| `cutting` | $`a/(2\zeta)`$ | $`a/(2\zeta)`$ | $`1/(2\pi a S)`$ | $`W_f`$ |
| `regenerative` | $`a/(2\zeta)`$ | $`a/(2\zeta)`$ | $`1/(2\pi a S)`$ | 0 |

The regenerative model also uses the cutting stiffness $`k_c = S a`$ and the blade count $`n_b \ge 1`$.

### 4.3 Damping (derivation sec. 5.6)

**Reference frequencies:** $`f_a`$ and $`f_t`$ are the case values if given, else $`f_{a1}`$ and $`f_{t1}`$.
Angular frequencies are $`\omega_a = 2\pi f_a`$ and $`\omega_t = 2\pi f_t`$.

**F1, fluid.**
- **Inputs:** fluid viscosity $`\mu_f`$ = PV if given, else $`\mu`$. Wall radius $`R = \max(D/2, 1.05 r_o)`$ per cell, with $`r_o = OD/2`$ and $`r_i = ID/2`$.
- **Couette drag:**

  ```math
  c_{0a} = \frac{2\pi\mu_f}{\ln(R/r_o)},\qquad c_{0t} = \frac{4\pi\mu_f r_o^2 R^2}{R^2 - r_o^2}.
  ```

- **Stokes layer:**

  ```math
  c_{1a} = \sqrt{\rho_m\mu_f\omega_a/2}\;2\pi(r_o + r_i),\qquad c_{1t} = \sqrt{\rho_m\mu_f\omega_t/2}\;2\pi(r_o^3 + r_i^3).
  ```

- **Rates:** $`k_{0a} = k_a + c_{0a}/m_a`$ and $`k_{0t} = k_t + c_{0t}/m_t`$; $`k_{1a} = c_{1a}/m_a`$ and $`k_{1t} = c_{1t}/m_t`$.
- **Fluid damping off:** when `fluid` is false or $`\mu_f \le 0`$, $`k_0 = k_a, k_t`$ and $`k_1 = 0`$.
- **High-pass coefficients:**
  $`a_{hp,t} = 1 - e^{-\Delta t\, f_c\, \omega_t}`$ and $`a_{hp,a} = 1 - e^{-\Delta t\, f_c\, \omega_a}`$
  (1 if $`f_c \le 0`$).

**M1, material** (when $`\eta > 0`$).
- $`z = \eta/2`$ and $`\omega_2 = 2\pi f_{max}`$.
- For each direction with reference $`\omega_1`$: if $`\omega_2 > \omega_1`$, then
  $`\alpha = 2z\omega_1\omega_2/(\omega_1+\omega_2)`$ and $`\beta = 2z/(\omega_1+\omega_2)`$; otherwise $`\alpha = 0`$ and $`\beta = 2z/\omega_1`$.
- $`k_1 \mathrel{+}= \alpha`$.
- **Kelvin–Voigt face conductances**, for interior faces $`j = 1\dots n-1`$:

  ```math
  k_{v,j} = \beta\,\tfrac12(K_{j-1} + K_j)\,/\,\tfrac12(h_{j-1} + h_j),
  ```

  with $`k_{v,0} = k_{v,n} = 0`$, separately for torsion ($`K_t`$) and axial ($`K_a`$).
- The KV solve is active only if $`\eta > 0`$.

**F2, yield drag.** $`f_{y,j} = \tau_y\,\pi\,OD_j`$ (N/m) if $`\tau_y > 0`$ and `yield_drag` ≠ `off`. It applies:
- always, for `always`;
- only while the current pump rate is ≤ 0, for `pumps_off`.

---

## 5. The time step (kernel)

One step advances all states by $`\Delta t`$. Notation: $`t^+ = t + \Delta t`$.

### 5.1 Ramped inputs

Five ramped inputs, each a row (start value $`v_0`$, start time $`t_0`$, ramp time $`T_r`$, target $`v_1`$):
- speed set point $`\omega_{sp}`$;
- block velocity;
- motor torque;
- prescribed top-drive speed;
- WOB set point.

Each is evaluated at $`t^+`$:

```math
r(t) = \begin{cases} v_1 & T_r \le 0 \text{ or } t \ge t_0 + T_r\\ v_0 & t \le t_0\\ v_0 + (v_1 - v_0)(t - t_0)/T_r & \text{otherwise.}\end{cases}
```

### 5.2 Top-drive control (mode `ctrl`)

Filter factors are $`a_x = e^{-\Delta t/t_x}`$ (0 if $`t_x \le 0`$).

**`ctrl = 0`, internal controller:**

```math
\begin{aligned}
\omega_{sp} &\leftarrow r_{\omega sp}(t^+),\\
\omega_f &\leftarrow a_{en}\omega_f + (1-a_{en})\,\omega_{td},\\
y_{lp} &\leftarrow a_{lp}y_{lp} + (1-a_{lp})\,\tau_0,\qquad y_{hp} \leftarrow a_{hp}y_{hp} + (1-a_{hp})\,y_{lp},\\
e &= \omega_{sp} - \omega_f - Z\,(y_{lp} - y_{hp}),\\
u &= K_p e + K_i (I_{PI} + e\Delta t).
\end{aligned}
```

If a torque limit $`\tau_{max} > 0`$ applies and $`|u| > \tau_{max}`$, then
$`u \leftarrow \tau_{max}\operatorname{sign}u`$ and the integrator is frozen;
otherwise $`I_{PI} \mathrel{+}= e\Delta t`$. Finally $`\tau_m \leftarrow u`$.

**`ctrl = 1`, external torque:** $`\tau_m \leftarrow r_{\tau m}(t^+)`$ and $`\omega_{sp} \leftarrow r_{\omega sp}(t^+)`$.

**Top-drive dynamics.** Let the incoming invariant be $`\beta_1 = \omega_0 - \tau_0^{cell}/Z_{t,0}`$
(cell 0 values). Then

```math
\omega_{td}^+ = \begin{cases} r_{\omega td}(t^+) & \text{ctrl} = 2\\[2pt]
\dfrac{I\omega_{td}/\Delta t + \tau_m + Z_{t,0}\beta_1}{I/\Delta t + Z_{t,0}} & \text{otherwise,}\end{cases}
\qquad \tau_0 = Z_{t,0}(\omega_{td}^+ - \beta_1).
```

For `ctrl = 2`: $`\tau_m = \tau_0 + I(\omega_{td}^+ - \omega_{td})/\Delta t`$ and $`\omega_{sp} = \omega_{td}^+`$.
Face 0 gets $`(\omega_{td}^+, \tau_0)`$.

### 5.3 Drawworks (mode `axial`)

- **Autodriller** (axial = 1):

  ```math
  W_s = H_{ref} - H,\quad \epsilon = r_{WOB}(t^+) - W_s,\quad v = K_{p,AD}\epsilon + K_{i,AD}(I_{AD} + \epsilon\Delta t),
  ```

  clipped to $`\pm v_{max}`$. When clipped the integrator is frozen; otherwise $`I_{AD} \mathrel{+}= \epsilon\Delta t`$. Then $`v_{blk} = v`$.
- **Prescribed block velocity** (axial = 0): $`v_{blk}`$ follows $`r_{vblk}(t^+)`$. With an acceleration limit $`a > 0`$, it changes by at most $`a\Delta t`$ per step.
- **Top face.** With $`\beta_{a1} = v_0 - q^d_0/Z_{a,0}`$:

  ```math
  q_0 = Z_{a,0}(v_{blk} - \beta_{a1}),\quad v_{top} = v_{blk}.
  ```

  The hook can only pull up: if $`-(q_0 + q^s_{top}) < 0`$, then
  $`q_0 = -q^s_{top}`$ and $`v_{top} = \beta_{a1} + q_0/Z_{a,0}`$.
  Face 0 gets $`(v_{top}, q_0)`$, and the hookload is $`H = -(q_0 + q^s_{top})`$.

### 5.4 Interior faces (Godunov)

For faces $`j = 1\dots n-1`$, with the outgoing invariant from the left $`\alpha = \omega_{j-1} + \tau_{j-1}/Z_{j-1}`$,
the incoming invariant from the right $`\beta = \omega_j - \tau_j/Z_j`$, and the same with $`(v, q^d, Z_a)`$:

```math
F = \frac{\alpha - \beta}{1/Z_{j-1} + 1/Z_j},\qquad \text{velocity} = \alpha - F/Z_{j-1},\qquad \text{force} = F.
```

### 5.5 Bit face

Outgoing invariants of the bottom cell:
$`\alpha_t = \omega_{n-1} + \tau_{n-1}/Z_{t,n-1}`$ and $`\alpha_a = v_{n-1} + q^d_{n-1}/Z_{a,n-1}`$.

**Quasi-static models (coulomb, cutting).**
- **Off bottom:** if $`s_{bit} < s_{hole} - 10^{-9}`$, or if $`\alpha_a \le 0`$ (the string is lifting the bit), the bit is free:
  $`(\omega_b, \text{TOB}, v_b, \text{WOB}, \text{stuck}) = (\alpha_t, 0, \alpha_a, 0, 0)`$.
- **On bottom:** start with $`w = \max(\alpha_t, 0)`$ and iterate at most 8 times:
  1. **Axial.** $`q_{lin} = Z_a\alpha_a`$.
     - If $`q_{lin} \le W_0`$ or $`k_r|w| = 0`$: WOB = $`q_{lin}`$ and $`v_b = 0`$.
     - Otherwise:

       ```math
       \text{WOB} = \frac{\alpha_a + k_r W_0 |w|}{1/Z_a + k_r|w|},\qquad v_b = k_r(\text{WOB} - W_0)|w|.
       ```

  2. **Torsional.** $`\tau_{stop} = Z_t\alpha_t`$.
     - **Sticking:** if ($`|\omega_b^{prev}| \le \omega_c`$ or stuck$`^{prev}`$) and $`|\tau_{stop}| \le \ell_s\,\text{WOB}`$, then TOB = $`\tau_{stop}`$, $`w' = 0`$, stuck = 1.
     - **Sliding:** otherwise TOB $`= \operatorname{sign}(\alpha_t)\min(\ell_k\text{WOB}, |\tau_{stop}|)`$, $`w' = \alpha_t - \text{TOB}/Z_t`$, stuck = 0.
  3. **Convergence:** stop when $`|w' - w| < 10^{-12}`$; otherwise set $`w \leftarrow w'`$.

  Result: $`\omega_b = w`$.

**Regenerative model.**
- **Blades:** blade $`i`$ ($`i = 0\dots n_b-1`$) sits at $`\phi_b + 2\pi i/n_b`$. Bin width is $`2\pi/M`$ with $`M = 720`$, and the surface under blade $`i`$ is $`s_i`$ = surf[bin].
- **Axial (implicit).** With $`u_f = s_{bit} + \Delta t\,\alpha_a`$, solve the following for $`u_1`$ by fixed-point iteration (at most $`n_b+2`$ times, starting from $`u_1 = u_f`$, until the change is below 1e-15):

  ```math
  u_1 = \frac{u_f + g\sum_{s_i<u_1}s_i}{1 + g\,m},\qquad g = \Delta t\,k_c/Z_a,
  ```

  where $`m`$ is the number of blades with $`s_i < u_1`$; if $`m = 0`$, $`u_1 = u_f`$.
- **Contact:** WOB $`= k_c\sum_i \max(u_1 - s_i, 0)`$ and $`v_b = (u_1 - s_{bit})/\Delta t`$.
- **Torque:** the same sticking and sliding rule as above, with $`\ell = \ell_k`$ for both.

**Bit face values.** Face $`n`$ gets $`(\omega_b, \text{TOB})`$ and $`(v_b, \text{WOB})`$. The scalar bit state is updated.

### 5.6 Cell update

For every cell (all with old values on the right):

```math
\begin{aligned}
\tau_j &\mathrel{-}= \Delta t\,\frac{K_{t,j}}{h_j}(\omega^f_{j+1} - \omega^f_j), &
\omega^*_j &= \omega_j - \frac{\Delta t}{m_{t,j}h_j}(\tau^f_{j+1} - \tau^f_j),\\
q^d_j &\mathrel{-}= \Delta t\,\frac{K_{a,j}}{h_j}(v^f_{j+1} - v^f_j), &
v^*_j &= v_j - \frac{\Delta t}{m_{a,j}h_j}(q^f_{j+1} - q^f_j).
\end{aligned}
```

**Linear fluid damping (implicit):**

```math
\omega^{**}_j = \frac{\omega^*_j + \Delta t\,k_{1t,j}\bar\omega_j}{1 + \Delta t(k_{0t,j} + k_{1t,j})},\qquad
v^{**}_j = \frac{v^*_j + \Delta t\,k_{1a,j}\bar v_j}{1 + \Delta t(k_{0a,j} + k_{1a,j})}.
```

**Kelvin–Voigt (implicit, if active).** For $`x = \omega^{**}`$ with $`(m_t, k_{vt})`$,
and $`x = v^{**}`$ with $`(m_a, k_{va})`$, solve the tridiagonal system with the Thomas algorithm:

```math
\frac{m_j h_j}{\Delta t}(x'_j - x_j) = k_{v,j+1}(x'_{j+1} - x'_j) - k_{v,j}(x'_j - x'_{j-1}).
```

### 5.7 Along-string friction

Let $`\omega^{old}, v^{old}`$ be the cell velocities before this step, $`(\omega^s, v^s)`$ the
damped predictions, $`r = r_{c,j}`$, and $`f_y = f_{y,j}`$.

**Coupled** (`along_string_coupling`, default).
- **Normal force:** $`T = -(q^s_{c,j} + q^d_j)`$ (with $`q^d`$ already updated) gives $`F_N = |T\boldsymbol\kappa_j + w_{b,j}\mathbf z_{\perp,j}|`$.
- **Friction factors:** $`\mu_{sa} = \mu_k + c_{ax}(\mu_s - \mu_k)`$ and $`\mu_{st} = \mu_s`$.
- **Sliding direction:** $`v_t = \omega^s r`$, speed $`\sigma = \sqrt{(v^s)^2 + v_t^2}`$, direction $`(e_a, e_t) = (v^s, v_t)/\sigma`$ (if $`\sigma = 0`$: $`(0, 1)`$).
- **Normalised old speed:** $`s_{old} = \sqrt{(v^{old}/v_c)^2 + (\omega^{old}/\omega_c)^2}`$.
- **Static capacities:** $`f_{as} = \mu_{sa}F_N + f_y`$ and $`f_{ts} = \mu_{st}F_N + f_y`$.
  - The static resistance along the sliding direction is $`f_{sd} = 1/\sqrt{(e_a/f_{as})^2 + (e_t/f_{ts})^2}`$, if $`\sigma > 0`$ and both capacities are positive.
  - Otherwise $`f_{sd} = f_{ts}`$.
- **Stick check** when $`s_{old} \le 1`$:

  ```math
  p_a = \frac{m_a v^s}{\Delta t f_{as}},\quad p_t = \frac{m_t\omega^s}{r\Delta t f_{ts}}
  ```

  (each $`10^{30}`$ if its capacity is ≤ 0). If $`p_a^2 + p_t^2 \le 1`$, set $`v = \omega = 0`$ and skip to the next cell.
- **Resistance per length** $`f_{res}`$:
  - $`f_{sd}`$ if $`s_{old} \le 1`$;
  - for `exponential` weakening, $`\mu_kF_N + f_y + (f_{sd} - \mu_kF_N - f_y)e^{-s_{old}}`$;
  - otherwise $`\mu_kF_N + f_y`$.
- **Implicit 2-D Coulomb** on $`(v^s, v_t)`$ with masses $`m_a`$ and $`m_t/r^2`$ and impulse $`c = \Delta t f_{res}`$:
  - With momenta $`p_a = m_a v^s`$ and $`p_t = (m_t/r^2)v_t`$, if $`p_a^2 + p_t^2 \le c^2`$ the cell stops.
  - If $`c \le 0`$, the velocities are unchanged.
  - Otherwise find $`\sigma' \ge 0`$ by Newton iteration from $`\sigma`$ (at most 50 iterations; stop when the change is ≤ $`10^{-14}\sigma`$; clip at 0):

    ```math
    g(\sigma') = \sum_{i\in\{a,t\}}\left(\frac{p_i}{m_i\sigma' + c}\right)^2 - 1 = 0.
    ```

  - Then $`v = v^s m_a\sigma'/(m_a\sigma' + c)`$, $`v_t' = v_t m_t'\sigma'/(m_t'\sigma' + c)`$ with $`m_t' = m_t/r^2`$, and $`\omega = v_t'/r`$.

**Decoupled.** The normal force is $`F_N = F_{N0,j}`$.
- **Torsion.** $`t_{stop} = m_t\omega^s/\Delta t`$ is clipped to $`\pm C_t`$, and then $`\omega = \omega^s - \Delta t\,t_{stop}/m_t`$. The capacity $`C_t`$ is:
  - $`(\mu_sF_N + f_y)r`$ if $`|\omega^{old}| \le \omega_c`$;
  - for exponential weakening, $`((\mu_k + (\mu_s-\mu_k)e^{-|\omega^{old}|/\omega_c})F_N + f_y)r`$;
  - otherwise $`(\mu_kF_N + f_y)r`$.
- **Axial.** The same rule with $`f_{stop} = m_a v^s/\Delta t`$, $`\mu_{sa}`$ in place of $`\mu_s`$, $`v_c`$ in place of $`\omega_c`$, and no factor $`r`$.

### 5.8 Low-pass states and bookkeeping

```math
\bar\omega_j \mathrel{+}= a_{hp,t}(\omega_j - \bar\omega_j),\qquad \bar v_j \mathrel{+}= a_{hp,a}(v_j - \bar v_j).
```

Then:
- $`t \mathrel{+}= \Delta t`$;
- $`s_{bit} \mathrel{+}= v_b\Delta t`$;
- $`x_{blk} \mathrel{+}= v_{blk}\Delta t`$;
- $`W \mathrel{+}= |\omega_b|\Delta t`$.

**Regenerative model only:** the blades sweep from $`\phi_b`$ to $`\phi_b' = \phi_b + \omega_b\Delta t`$.
- Only forward rotation sweeps.
- For each blade, every bin from $`\lfloor(\phi_b + \text{off})/\Delta\theta\rfloor`$ up to (excluding) $`\min(\lfloor(\phi_b' + \text{off})/\Delta\theta\rfloor, \text{start} + M)`$ has its surface raised to $`s_{bit}`$ if it is below.
- $`\phi_b \leftarrow \phi_b'`$, reduced modulo $`2\pi`$ once it exceeds $`10^6`$.

In all models, $`s_{hole} \leftarrow \max(s_{hole}, s_{bit})`$.

**Recording.** Every `rec_every` = $`\max(1, \operatorname{round}(1/(f_{rec}\Delta t)))`$
steps, one row is recorded with these columns:
- time;
- $`\omega_{td}`$, $`\tau_m`$, $`\tau_0`$;
- $`\omega_b`$, TOB;
- $`v_{blk}`$, $`H`$;
- $`v_b`$, WOB;
- $`s_{bit}`$, $`s_{hole}`$;
- $`\omega_{sp}`$.

---

## 6. Inputs, ramps and modes

`step(duration, ...)` takes optional inputs; an input that is not given keeps its value.
Giving an input selects its mode:

| Input | Mode set | Ramp row | Start value when taking over from another mode |
|---|---|---|---|
| `omega_sp` | ctrl 0 | ω_sp | $`\omega_{td}`$ |
| `motor_torque` | ctrl 1 | τ_m | $`\tau_m`$ |
| `omega_td` | ctrl 2 | ω_td | $`\omega_{td}`$ |
| `v_block` | axial 0 | v_blk | $`v_{blk}`$ |
| `wob_sp` | axial 1 | WOB | $`H_{ref} - H`$ |

**Ramp start.** Within the same mode, a ramp to a new target starts from the
ramp's current value at the present time. A target equal to the current
target is ignored. The ramp row becomes (start, $`t`$, ramp time, target),
with the ramp times of section 1.

**Slow ramps.** Pump rate and inlet temperature ramp the same way, evaluated
at the time they are read. At initialisation the pump ramps from 0 at
$`t = 0`$ to the case flow rate over the flow ramp time; the inlet
temperature starts at its value.

**Zeroing the weight indicator.** This sets $`H_{ref} \leftarrow H`$.

**Splitting and pauses.** Time stepping uses integer step counts:
$`N_{target} = \operatorname{round}(t_{requested}/\Delta t)`$ (accumulated
over all calls). Steps run in chunks that end on multiples of the slow-tick
interval $`K_s = \max(1, \lfloor 1/\Delta t\rfloor)`$ steps, and a slow tick
follows each completed multiple. Splitting a run into calls of any length
therefore gives identical results.

---

## 7. Slow tick

At every slow tick (about every 1 s of simulated time):

1. **Geometry refresh** (section 2.4), if the bit has moved ≥ 1 m.
2. **Hydraulics:** advance from the last hydraulics time to now (section 8.3).
   - The pump rate ramps linearly from its value at the last tick to the current ramped value.
   - The bit and hole depths are the current ones.
3. **Every `thermal_update_period` ticks** (default 2):
   - average ROP = $`\max(s_{hole} - s_{hole}^{last}, 0)/\Delta t_{slow}`$;
   - average bit speed = $`(W - W^{last})/\Delta t_{slow}`$;
   - advance the temperature model (section 10; full model only);
   - update the rock strength (section 9) with these averages.

   The rock strength is also computed at initialisation and on parameter updates, with ROP = ω = 0.

**Initialisation order:**
1. Grid and properties.
2. Geometry refresh.
3. Parameters, with the `stiff` controller.
4. $`H_{ref}`$.
5. Hydraulics in the static state (zero flow).
6. Rock strength.

Then the case controller is installed and the parameters are packed again,
which repeats the geometry refresh with no change and the rock strength.

---

## 8. Hydraulics

### 8.1 Rheology and friction gradients

Bingham if a yield point is given (PV = `plastic_viscosity` or $`\mu`$, $`\tau_y`$ = `yield_point`),
else Newtonian (PV = $`\mu`$, $`\tau_y = 0`$).

The friction gradient at speed $`|u| > 0`$ is $`dp/dL = a|u| + b`$, signed with $`u`$:
- **Pipe** (diameter $`D`$): $`a_\ell = 32\,PV/D^2`$, $`b_\ell = 16\tau_y/(3D)`$, $`D_e = D`$, $`\mu_a = PV + \tau_yD/(6|u|)`$.
- **Annulus** (outer $`D_2`$, inner $`D_1`$, $`h = D_2 - D_1`$): $`a_\ell = 48\,PV/h^2`$, $`b_\ell = 6\tau_y/h`$, $`D_e = 0.816h`$, $`\mu_a = PV + \tau_yh/(8|u|)`$.
- **Regime:** $`Re = \rho|u|D_e/\mu_a`$.
  - Below 2100: $`(a, b) = (a_\ell, b_\ell)`$.
  - Otherwise $`f = 0.0791\,(\rho|u|D_e/PV)^{-0.25}`$ and $`(a, b) = (2f\rho|u|/D_e, 0)`$.

**Bit nozzles:** $`\Delta p_{bit} = \rho Q|Q|/(2C_d^2A_n^2)`$. **Surface loss:** $`c|Q|^m`$.

**Steady textbook system** (`steady_state`, constant density):
- $`SPP = \Delta p_{surf} + \sum\Delta p_{pipe} + \Delta p_{bit} + \sum\Delta p_{ann}`$;
- $`BHP = \rho g\,TVD + \sum\Delta p_{ann}`$;
- $`ECD = \rho + \sum\Delta p_{ann}/(g\,TVD)`$.

### 8.2 Transient model: grid and geometry

There are two flow paths on fixed cells:
- **Drillstring bore:** surface to bit; the edges are $`\operatorname{linspace}(0, s_{bit}, n_d+1)`$ with $`n_d = \max(1, \operatorname{round}(s_{bit}/\Delta x_h))`$.
- **Annulus:** surface to the hole bottom, built the same way with $`s_{hole}`$.

Cell centres have TVD $`z`$, and $`TVD_{bottom}`$ is the TVD at $`s_{hole}`$.

**Volume-averaged cross sections.** The string occupies MD $`[s_{bit} - L_s, s_{bit}]`$, with $`L_s = \sum h_j`$.
Integrals of $`\frac\pi4 OD^2`$ and $`\frac\pi4 ID^2`$ along the string (piecewise
constant per string cell, linear interpolation of the cumulative integral)
give, per hydraulic cell:

| Quantity | Annulus | Bore |
|---|---|---|
| wall / pipe diameter | $`D_o`$ = wall diameter (2.2) | $`D_d = \sqrt{4A_{bore}/\pi}`$ |
| inner diameter | $`D_i = \sqrt{4A_{od}/\pi}`$ (0 if < 1 mm) | 0 |
| area $`A`$ | $`\frac\pi4D_o^2 - A_{od}`$ | $`\max(A_{bore}, 10^{-8})`$ |
| volume | $`A L`$ | $`A L`$ |

The bit enters annulus cell $`k_{ab}`$, the cell containing $`s_{bit}`$ (clipped to the last cell).

**Density** (Wang et al. 2023, eq. 2):
$`\rho = \rho_T + p/c_m^2`$, with $`\rho_T = (1 - (T - T_\text{ref})\beta_T)\rho_m`$.
In the lite version $`T = T_\text{ref}`$, so $`\rho_T = \rho_m`$. The cell
pressure (gauge) is $`p = c_m^2(m/V - \rho_T)`$.

**State:**
- mass per cell, $`m_d`$ and $`m_a`$;
- mass flux per face, in kg/s, positive down:
  - bore faces $`G_d[k]`$ (top of cell $`k`$; $`G_d[0]`$ is the pump inlet);
  - annulus faces $`G_a[k]`$ (top of cell $`k`$; $`G_a[0]`$ is the flowline outlet);
  - the bit face $`G_b`$.

### 8.3 Transient model: advance

`advance(duration, q0, q1, s_bit, s_hole)`:

1. **Regrid** the moving bottom cells on the current state (grow–cut–grow).
   - **Annulus,** when the hole deepens: if the new bottom cell length exceeds $`1.5\Delta x_h`$, split it at the midpoint of its new extent (if that is inside the old cell). Mass and volume split by fraction, the new face flux is half the bottom face flux, and the edges are updated.
   - **Bore:**
     - Let $`L`$ be the new bottom length.
     - If $`L > 1.5\Delta x_h`$, split at the midpoint (kept 1e-6 m inside the old cell); the new face flux is the mean of the bottom face flux and $`G_b`$.
     - If $`L < 0.5\Delta x_h`$ and there is more than one cell, merge the bottom two cells.
2. **Pipe velocity** $`v_p = (s_{bit} - s_{bit}^{old})/\text{duration}`$. Set the last edges to the new depths and recompute the geometry (8.2) at the end of the interval. Volumes interpolate linearly from the old geometry (or the new one, if the cell count changed) to the new one over the interval.
3. **Substeps:** $`n_{sub} = \lceil \text{duration}/(\text{CFL}\,L_{min}/c_m)\rceil`$, each of length $`\Delta t_h`$, with $`w = (s+1)/n_{sub}`$ at substep $`s`$:
   1. **Pump face:** $`G_d[0] = \rho_{T,d0}\,(q_0 + (q_1 - q_0)w)`$.
   2. **Cell values:** for each cell, $`V`$ by interpolation, $`A = V/L`$, $`\rho = m/V`$, $`p = c_m^2(\rho - \rho_T)`$.
   3. **Face updates** with the face rule below:
      - interior bore faces (slip factor 1 relative to $`v_p`$);
      - the bit face, bore bottom cell to annulus cell $`k_{ab}`$: no slip, with the nozzle term $`R_n|G|`$, where $`R_n = 1/(2C_d^2A_n^2\bar\rho)`$;
      - interior annulus faces (slip factor $`K_c`$ if both cells contain pipe, else 0);
      - the outlet face. Its drive is $`p_{back} - p_{a0} + \rho_{a0} g z_{a0}`$; there is no upstream half cell, and the slip factor is $`K_c`$ if the cell contains pipe.
   4. **Mass:**
      - $`m_{d,j} \mathrel{+}= \Delta t_h(G_{d,j} - G_{d,j+1})`$, with $`G_b`$ leaving the bottom bore cell;
      - $`m_{a,j} \mathrel{+}= \Delta t_h(G_{a,j} - G_{a,j+1})`$, with no flux through the bottom of the annulus;
      - $`m_{a,k_{ab}} \mathrel{+}= \Delta t_h G_b`$;
      - accumulators: pumped $`\mathrel{+}= \Delta t_h G_d[0]`$, returned $`\mathrel{-}= \Delta t_h G_a[0]`$.

**Face rule** (backward Euler, Picard iterations).
- **Inputs:** upstream (u) and downstream (d) half cells, drive $`P`$ = pressure difference plus $`\bar\rho g\Delta z`$, and inertance $`1/\mathcal I = 1/\big(\tfrac12(L_u/A_u + L_d/A_d)\big)`$.
- **Pipe-dragged flux:** $`G_s = k_s v_p\bar\rho\bar A`$.
- **Iterations** (a set number, default 2): from the current flux $`G_k`$ accumulate over the half cells
  $`R \mathrel{+}= \tfrac12 L\,a(|G_k/(\rho A) - k_sv_p|)/(\rho A)`$ and $`B \mathrel{+}= \tfrac12 L\,b`$, plus $`R \mathrel{+}= R_n|G_k|`$ at the bit. Then:

  ```math
  c = \Delta t_h/\mathcal I,\quad F = G + c(P + RG_s),\quad \delta = 1 + cR,\quad r = F - \delta G_s,
  ```

  - **Yield holds:** if $`|r| \le cB`$, then $`G_k = G_s`$.
  - **Otherwise:** $`G_k = (F - cB\operatorname{sign} r)/\delta`$.

**Static and steady states** (`set_steady(Q)`, $`Q = 0`$ for static).
This is the exact discrete steady state with the pipe at rest:
- **Fluxes:** $`G_d = \rho_{T,d0}Q`$ on every bore face and at the bit; $`G_a = -G_d`$ on annulus faces $`\le k_{ab}`$, 0 below.
- **Annulus pressures,** from the outlet down, by fixed point (50 iterations) on the face densities:
  - outlet: $`p_{a0} = p_{back} + \rho_0 g z_0 - \tfrac12L_0\,(dp/dL)_0`$;
  - each face $`k`$: $`p_k = p_{k-1} + \bar\rho g(z_k - z_{k-1}) - \text{friction of the two half cells}`$.
- **Bit:** $`p_{d,last} = p_{a,k_{ab}} - \bar\rho g(z_{a,k_{ab}} - z_{d,last}) + \Delta p_{bit} + \text{friction}`$.
- **Bore pressures:** upward from the bit, $`p_{k-1} = p_k - \bar\rho g(z_k - z_{k-1}) + \text{friction}`$.
- **Masses:** $`m = (\rho_T + p/c_m^2)V`$.

**Outputs** (`state`):
- **Standpipe pressure:**
  - $`p_{top} = p_{d0} - \rho_{d0}gz_{d0} + \tfrac12L_{d0}\,(dp/dL)(u_0 - v_p)`$, with $`u_0 = G_d[0]/(\rho_{d0}A_{d0})`$;
  - $`SPP = p_{top} + c|q|^m`$, with $`q = G_d[0]/\rho_{T,d0}`$.
- **Bottomhole pressure and ECD:**
  - $`BHP = p_{a,last} + \rho_{a,last}\,g\,(TVD_{bottom} - z_{a,last})`$;
  - $`ECD = BHP/(g\,TVD_{bottom})`$.
- **Flows:**
  - $`q_{in} = q`$;
  - $`q_{out} = -G_a[0]/\rho_{a0}`$;
  - $`q_{bit} = G_b/\tfrac12(\rho_{d0} + \rho_{a0})`$.
- **Nozzle loss:** $`\Delta p_{bit}`$ at the mean of the bit-face densities.
- **Pit gain:** $`-(\text{pumped} - \text{returned})/\rho_m`$.

---

## 9. Rock strength at the bit

For the layer containing the hole depth (the deepest layer with top ≤ MD),
at the hole TVD $`z`$:

- **Rock temperature:** $`T_r = T_{surf} + \gamma z`$ (defaults 40 °C, 0.0691 °C/m).
- **T1, thermal weakening** (defaults $`T^* = 193`$ °C, $`b = 1.41\times10^{-3}`$ /K):

  ```math
  UCS(T) = UCS_0\max\big(1 - b\max(T_r - T^*, 0),\, 0.05\big).
  ```

- **P1, confinement.**
  - Pressures: pore $`p_p = 9.79\text{ kPa/m}\cdot z`$; overburden $`\sigma_v = 25\text{ kPa/m}\cdot z`$.
  - Overbalances:
    - drained, $`\Delta p_{hp} = p_{bh} - p_p`$;
    - undrained, $`\Delta p_{lp} = p_{bh} - (p_p - B(\sigma_v - p_{bh})/3)`$, with $`B = 1`$.
  - Blend weight $`w`$ from the effective porosity $`\phi`$: 0 at or below 0.05, 1 at or above 0.20, linear between.
  - Confined strength:

    ```math
    CCS = UCS(T) + \big((1-w)\Delta p_{lp} + w\Delta p_{hp}\big)\frac{1+\sin\varphi}{1-\sin\varphi}.
    ```

- **T2, thermal stress** (full model only, when a mud temperature at the bit is known):
  - Stress: $`\sigma_{th} = \alpha E(T_r - T_{mud})/(1-\nu)`$.
  - Skin ratio: if ω > 1e-6 and ROP > 1e-9, it is $`\min(1, \delta/d_b)`$, with:
    - time between blade passes $`t_p = 2\pi/(n_b\omega)`$;
    - thermal skin $`\delta = \sqrt{\pi\kappa_r t_p}`$, where $`\kappa_r = k_r/(\rho_r c_r)`$;
    - depth of cut per blade $`d_b = \text{ROP}\,t_p`$.

    Otherwise the ratio is 1.
  - Effective strength: $`CCS_{eff} = \max(CCS - \kappa_{th}\sigma_{th}\cdot\text{ratio},\; 0.1\,UCS(T))`$.

  Without T2, $`CCS_{eff} = \max(CCS, 0.1\,UCS(T))`$.
- **Bottomhole pressure** $`p_{bh}`$: from the hydraulics if enabled; else from the temperature model if enabled; else $`\rho_m g z`$.
- **Drilling strength:** $`S = \zeta\,CCS_{eff}/\eta_c`$. With `strength_source = formation`, $`k_r = 1/(2\pi aS)`$ and $`k_c = Sa`$ are updated.

---

## 10. Wellbore temperature (full model only)

The Wang et al. (2023) single-phase energy balance (`thermal.py`).

**Grid.** Axial cells of about 30 m cover the hole; the bottom cell grows with
the hole and is cut at 1.5 dx. Each cell has:
- drillstring and annulus fluid temperatures;
- 12 radial formation cells, log-spaced to 5 m beyond the wall;
- a far-field geothermal temperature.

**Heat transfer.**
- Conductances drillstring ↔ annulus and annulus ↔ formation are series resistances:
  - convection: $`Nu = 0.023Re^{0.8}Pr^{0.4}`$ above Re 2300, else 3.66;
  - pipe-wall conduction ($`k_s`$ = 45 W/m K);
  - casing (0.012 m) where cased.
- Both are volume-averaged over the string in each cell.

**Time stepping.** Backward Euler with substeps ≤ 2 s and 3 coupling iterations:
- drillstring cells upwind from the inlet down to the bit;
- annulus cells upwind from the bit up (below the bit: stagnant);
- cuttings enter the bottom cell at the formation temperature, at a rate of $`\rho_r c_r A_{hole}\,\text{ROP}`$;
- each formation layer is a tridiagonal radial conduction solve.

On growth, the bottom formation cells relax toward the far-field temperature.

**Annulus pressure.** It is integrated with $`\rho(T, p)`$ and Newtonian
friction to give $`p_{bh}`$ and ECD. When hydraulics is enabled, the
hydraulics instead take the mud temperatures (interpolated to their cells)
for $`\rho_T`$.

Outputs: mud temperature at the bit, outlet temperature and profiles. The
lite version omits this section.

---

## 11. Steady torque and drag (reference profiles)

For one operating state (axial velocity $`v`$, speed $`\omega`$, WOB, TOB),
integrate from the bit up.

**At the bit:** $`T_n = -\text{WOB}`$ and $`M_n = \text{TOB}`$.

**For each cell, bottom to top:**
- $`F_N = |T_{j+1}\boldsymbol\kappa_j + w_b\mathbf z_\perp|`$, or a fixed profile;
- direction cosines $`(c_a, c_t) = (v, \omega r_c)/\sqrt{v^2 + \omega^2r_c^2}`$ (0 if both are 0);
- $`f = \mu_jF_N`$;
- tension: $`T_j = T_{j+1} + (w_bt_z - fc_a)h_j`$;
- torque: $`M_j = M_{j+1} + fc_tr_ch_j`$.

The geometry here uses the unshifted grid (initial position).

**UI reference profiles:**
- rotating ($`\omega = 1`$) with the kinetic factors;
- rotating with the static factors;
- hanging (μ = 0).

With decoupled friction, the static normal force is used.

---

## 12. Outputs

- **`latest`:**
  - top drive: time, $`\omega_{td}`$, $`\tau_m`$, $`\tau_0`$, $`\omega_{sp}`$;
  - bit: $`\omega_b`$, TOB, $`v_b`$, WOB;
  - axial: $`v_{blk}`$, $`H`$;
  - depths: $`s_{bit}`$, $`s_{hole}`$;
  - circulation: the current pump rate, plus hydraulics SPP, $`p_{bh}`$, ECD, $`q_{out}`$.
- **Recorded channels:** section 5.8.
- **Profiles per cell:**
  - MD, ω, τ, $`v`$;
  - tension $`-(q^s_c + q^d)`$;
  - normal force (live for coupled friction, static otherwise);
  - OD and inclination $`\arccos t_z`$.
- **Surface ROP / downhole ROP** in the UI:
  - block speed $`v_{blk}`$;
  - bit ROP $`\max(v_b, 0)`$ when $`s_{bit} \ge s_{hole} - 10^{-6}`$, else 0.
- **Surface WOB:** $`H_{ref} - H`$.

The UI frame format (`src/DrillPyze/ui/server.py`, `Session.frame`) is part of
the browser contract: the lite engine reproduces it field for field.

---

## 13. Conformance

`tests/test_lite_engine.py` runs the browser engine under Node.js against the
Python model on the same case (Well A with the temperature model off) and the
same input script. It requires:
- agreement to rounding level for the first seconds, before stick-slip chaos amplifies round-off;
- agreement of the statistics over the whole run (mean torque and hookload, stick-slip amplitude, hydraulic pressures).

A new spec revision must keep this test passing for the lite version.
