# DrillPyze coupled axial–torsional model

This note records the model equations and numerical choices implemented in
`DrillPyze.kernel`. The torsional part follows Aarsnes & Shor (2018) and the
benchmarking and control papers (Aarsnes et al. 2018; Auriol et al. 2022). The
axial part, the 2-D Coulomb inclusion and the bit boundary are extensions that
the earlier papers left as future work. They are marked **[extension]** below so
they can be revisited when the coupled-model paper is finished.

All quantities are SI. `x` is measured depth from surface (0) to the bit (L).

## 1. Field equations

Torsion (Aarsnes & Shor 2018, eqs. 1–3):

    dτ/dt + J G dω/dx = 0
    ρJ dω/dt + dτ/dx = −k_t ρJ ω − F_t(x,t)

Axial **[extension]**: the same structure, with `q` the axial *compressive*
force (q = −T, where T is effective tension) and `v` the axial velocity
(positive downhole):

    dq/dt + E A dv/dx = 0
    ρA dv/dt + dq/dx = w_b cos(inc) − k_a ρA v − F_a(x,t)

`w_b` is the buoyed weight per length, `k_t` and `k_a` are viscous damping
rates (1/s), and `F_t` (torque per length) and `F_a` (force per length) are
the along-string Coulomb friction.

Gravity is balanced exactly by splitting `q = q_s(x) + q_d(x,t)`. Here `q_s` is
the frictionless static compression of the hanging string, with `q_s = 0` at the
bit and `dq_s/dx = w_b cos(inc)`. Only `q_d` is transported. This makes the scheme
well balanced: a hanging string at rest stays exactly at rest.

## 2. Coupling through friction **[extension]**

At each cell, the pipe surface slides with velocity `V = (v, ω r_c)`, where `r_c`
is the contact radius. The Coulomb inclusion of Aarsnes & Shor (eq. 4) is
generalised to this 2-D velocity:

* **Normal force.** Uses the soft-string vector form, with the live tension
  `T = −(q_s + q_d)`:

      F_N = | T dt/ds + w_b z_perp |

  Torque therefore couples into friction through the axial stress state, the
  effect identified as missing in Section 5.5 of Aarsnes & Shor.
* **Regime.** Cells with `sqrt((v/v_c)² + (ω/ω_c)²) ≤ 1` at the previous step are
  *static*, with friction capped at `μ_s F_N`. Otherwise they are *kinetic*.
* **Static.** The friction required to hold the cell at rest this step is
  computed (the 2-D form of F̃ in eq. 20–21). It is applied in full if its
  magnitude is ≤ `μ_s F_N`, and otherwise scaled to that limit.
* **Kinetic.** Friction of magnitude `μ_k F_N` acts opposite to `V`. This is the
  velocity-split model, so rotation reduces axial drag and vice versa. Each
  component is clipped so friction never reverses the motion within a step,
  which is the `min(|F̃|, f_rat F_c)` limiter of eq. 22.

For pure torsion (`v ≡ 0`) this reduces exactly to the published scheme, with
`μ_k = f_rat μ_s`.

Cased and open hole have separate `μ_s` and `μ_k` (shoe depth from the
wellbore).

## 3. Numerical scheme

* **Discretisation.** Cell-centred finite volumes with Godunov (exact Riemann)
  fluxes at every cell face. Inside a uniform section this is the first-order
  upwind scheme on the Riemann invariants `α = ω + τ/Z_t` and `β = ω − τ/Z_t`
  (Aarsnes & Shor eqs. 11–24). At a change in impedance (collar–pipe, eqs.
  15–19) the face solution enforces continuity of `ω, τ` and of `v, q`
  automatically, so no special junction code is needed.
* **Impedances.** `Z_t = J sqrt(Gρ)` and `Z_a = A sqrt(Eρ)`.
* **Time step.** `Δt = CFL · min(Δx / c_a)`. The axial wave is about 1.6× faster
  than the torsional wave, so torsion runs at CFL ≈ 0.6. Upwinding at CFL < 1 is
  slightly diffusive. The analytic tests quantify this, and the error falls with
  `Δx`.
* **Source terms.** Applied by operator splitting after the flux update:
  gravity is already balanced, damping is implicit, and friction is applied as
  above.

## 4. Boundary conditions

### Top drive (torsion)

    I_TD dω_TD/dt = τ_m − τ_0,   τ_0 = Z_t1 (ω_TD − β_1)

This is integrated implicitly in `ω_TD`, so it is stable for any `I_TD`.
`τ_0` is the pipe torque (what a torque sub measures). `τ_m` is the motor
torque (what a VFD reports).

Control modes:

0. **Internal controller** (Aarsnes et al. 2018, Fig. 4). The motor torque is a
   PI controller acting on the speed error, with optional impedance feedback of
   filtered pipe torque:

       e = ω_SP − ω_0,filt − Z · HP(LP(τ_0))
       τ_m = K_p e + K_i ∫e,   saturated at ±τ_max, with anti-windup

   * Stiff controller: `Z = 0`, `K_p = 100 ζ_p`, `K_i = 5 I_TD`
   * SoftTorque: `Z = 0`, `K_p = 4 ζ_p`, `K_i = (2π f_c)² I_TD`
   * ZTorque: stiff PI gains with `Z = 1/ζ_p`, and high-pass and low-pass filters on `τ_0`.
     The inner loop must be stiff; impedance feedback on soft gains is unstable.
   * `ω_SP` is rate limited.
1. **External motor torque.** `τ_m` is supplied by the caller. This is for
   controllers written in Python and called at their own sample rate.
2. **Prescribed top-drive speed.** `ω_TD` is supplied directly. This is for
   replaying measured surface RPM and comparing the simulated torque against
   the measured torque.

### Drawworks (axial)

A prescribed block velocity `v_blk` (positive = lowering):

    q_0 = Z_a1 (v_blk − β_a1),   hookload = −q_0   (block weight excluded)

An autodriller mode adjusts `v_blk` with a PI controller on surface WOB
(reference off-bottom hookload − hookload).

### Bit **[extension]**

Incoming invariants `α_t = ω_N + τ_N/Z_t` and `α_a = v_N + q_N/Z_a`, with
`ω_b = α_t − TOB/Z_t` and `v_b = α_a − WOB/Z_a`.

* **Off bottom** (bit depth < hole depth): `WOB = TOB = 0`.
* **On bottom:** unilateral contact (`WOB ≥ 0`, so the bit lifts off if the
  string pulls it up) and a rock penetration law:

      ROP = k_r · max(WOB − W_0, 0) · |ω_b|

  with bit friction torque `TOB = μ_b(ω) · WOB · r_b`, where
  `μ_b(ω) = μ_bk + (μ_bs − μ_bk) exp(−|ω|/ω_b)`.
  The bit sticks when the torque needed to stop it is within `μ_bs · WOB · r_b`.
  The coupled bit equations are solved by fixed-point iteration each step.

Coupling between WOB, ROP and RPM enters in two places: TOB ∝ WOB, and
ROP ∝ WOB·RPM. This is a deliberately simple, calibratable surrogate. A
Detournay-style cutting and contact law with depth-of-cut delay (the
regenerative effect) is the planned next refinement.

## 5. Moving string

The grid is attached to the pipe. As the bit advances, cell positions in the
wellbore shift. `Simulation` periodically re-evaluates the cell geometry, the
cased/open friction and `q_s`. Total `q` is kept continuous by adjusting `q_d`.

## 6. Parameters (defaults)

| symbol | meaning | default |
|---|---|---|
| μ_s, μ_k (cased/open) | static / kinetic friction factors | 0.30/0.20, 0.35/0.25 |
| ω_c, v_c | static → kinetic thresholds | 1 rad/s, 0.01 m/s |
| k_t, k_a | viscous damping | 0.05 1/s, 0.05 1/s |
| I_TD | top-drive inertia | 2000 kg m² |
| r_b, μ_bs, μ_bk, ω_b | bit radius, bit friction | hole/2, 0.6, 0.4, 2 rad/s |
| k_r, W_0 | ROP law | 1e-9 m/(N rad), 0 |

Damping and the bit parameters are placeholders until they are calibrated
against field data.
