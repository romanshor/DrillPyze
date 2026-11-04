# DrillPyze — coupled drillstring dynamics

DrillPyze simulates the coupled axial and torsional dynamics of a drillstring
in a 3-D wellbore, together with steady-state torque and drag. It is built to:

* calibrate along-string friction factors from field data,
* study stick-slip in its off-bottom, on-bottom and axially driven forms, and
* develop and benchmark top-drive RPM controllers.

All quantities are **SI** internally. Conversion helpers in `DrillPyze.units`
handle oilfield units at input and output.

> **Branch `coupled-model`:** a rewrite of the 2023 prototype. The full
> mathematical derivation is in [`docs/derivation/`](docs/derivation/README.md).

## Mathematical model

The string is a slender elastic rod with $x \in [0, L]$ measured from the top
drive to the bit. Torsion follows Aarsnes & Shor (2018). The axial direction
uses the **identical** formulation, with Young's modulus $E$ and area $A$ in
place of $G$ and $J$:

$$
\begin{aligned}
&\partial_t \tau + JG\,\partial_x\omega = 0, &\qquad
&\rho J\,\partial_t\omega + \partial_x\tau = -k_t\rho J\omega - \mathcal F_t,\\
&\partial_t q + EA\,\partial_x v = 0, &\qquad
&\rho A\,\partial_t v + \partial_x q = w_b\cos\theta - k_a\rho A v - \mathcal F_a,
\end{aligned}
$$

where:
* $\omega, \tau$ are angular velocity and torque;
* $v, q$ are axial velocity and compressive force (tension $T = -q$);
* $w_b$ is the buoyed weight per length;
* $\mathcal F_t, \mathcal F_a$ are Coulomb friction inclusions with static and
  kinetic factors, separately for cased and open hole.

In Riemann invariants $\alpha = \omega + \tau/Z_t$ and
$\beta = \omega - \tau/Z_t$, with impedance $Z_t = J\sqrt{G\rho}$ (axially
$Z_a = A\sqrt{E\rho}$), each direction is a pair of transport equations at
$\pm c$, with $c_t = \sqrt{G/\rho}$ and $c_a = \sqrt{E/\rho}$.

**The two directions are coupled at the drill bit.** The characteristic
relations there are

$$
\omega_b = \alpha_t - \text{TOB}/Z_t, \qquad v_b = \alpha_a - \text{WOB}/Z_a,
$$

and they are closed by a Coulomb-type bit–rock law:

$$
|\text{TOB}| \le T_{s0} + \mu_{bs}\,r_b\,\text{WOB}\ \ \text{(stuck)},
\qquad
\text{TOB} = \mu_{bk}\,r_b\,\text{WOB}\ \ \text{(sliding)},
\qquad
\text{ROP} = k_r\,(\text{WOB} - W_0)^+\,|\omega_b| .
$$

A sticking torque must be overcome first; after that, the reactive torque is
linear in WOB. TOB therefore depends on the axial state, and penetration on
the torsional state. Along the string, axial and torsional friction are
independent, with the normal force of the static soft string. An opt-in 2-D
friction model (`along_string_coupling=True`) adds along-string coupling for
later studies.

Other model components:
* **Top drive:** an inertia driven by a stiff, SoftTorque or ZTorque
  controller, by an external torque, or by a prescribed speed.
* **Drawworks:** a prescribed block velocity, or an autodriller.
* **Gravity:** handled by a well-balanced static/dynamic split.
* **Numerics:** a Godunov finite-volume scheme (upwind on the invariants,
  exact at section changes) with a CFL-limited time step and projected
  friction (Aarsnes & Shor eqs. 20–22).

The derivation, boundary conditions, discrete scheme and analytic limits are
in [`docs/derivation/README.md`](docs/derivation/README.md).

## Installation

```bash
pip install -e ".[test]"
pytest -q
```

## Workflow: set up, initialise, step

```python
from DrillPyze import (Survey, Component, Drillstring, Wellbore, Fluid, Simulation,
                       ModelParameters, FrictionParameters, units, controllers, statics)

# 1. describe the well and string (SI; convert with units.FT, units.IN, ...)
survey = Survey(md=[0, 1500, 2100, 3000], inc=[0, 0, 60, 60], azi=[0, 0, 0, 0])
well = Wellbore(survey, casing_shoe=1000.0, hole_diameter=8.5 * units.IN)
string = Drillstring([                       # listed from the bit up
    Component('Collars', 230.0, od=5.75 * units.IN, id=2.8125 * units.IN),
    Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN),  # stretched to surface
])
params = ModelParameters(friction=FrictionParameters(0.40, 0.22, 0.45, 0.25))

# 2. set up and initialise
sim = Simulation(well, string, Fluid(1200.0), params)
sim.initialize(bit_depth=2500.0, hole_depth=2505.0, dx=10.0)

# 3. step: pause, change inputs, resume at any time
sim.step(30.0, omega_sp=60 * units.RPM)      # internal controller
sim.step(10.0, v_block=-0.05)                # pick up while rotating
df = sim.results()                           # recorded channels (DataFrame, SI)
prof = sim.profiles()                        # along-string state now

# steady-state torque and drag on the same grid
td = statics.torque_drag(well, sim.grid, 1200.0, mu_cased=0.22, mu_open=0.25, omega=1.0)
```

* Calling `step` repeatedly with short durations gives exactly the same result
  as one long call. That allows live streaming, and controllers in Python that
  run at their own sample rate via `step(dt, motor_torque=u)`.
* `snapshot()` returns an independent copy of the whole simulation, for
  ensembles.

See `examples/well_a_controllers.py`, which compares the stiff, SoftTorque and
ZTorque controllers on an off-bottom start-up.

## Performance

The time-stepping kernel is compiled with Numba. Its time step is set by the
axial wave speed (CFL 0.95), so cost scales with the square of the number of
cells per unit length. Measured on a 3,337 m string in a 2-core sandbox:

| cell size | cells | time step | speed vs real time | 60 min of drilling |
|---|---|---|---|---|
| 10 m  | 334   | 1.85 ms | 155× | ~23 s |
| 5 m   | 668   | 0.92 ms | 38×  | ~95 s |
| 2.5 m | 1,335 | 0.46 ms | 10×  | ~6 min |

The 2023 prototype ran at about 2× real time on a shorter 2,100 m string.

## Verification (`tests/`)

* Survey: minimum-curvature arc geometry, and azimuth crossing north.
* Torque and drag: analytic hanging weight; analytic straight inclined-hole
  pick-up, slack-off and rotating torque; rotation reduces drag.
* Dynamics:
  * torsional and axial wave arrival times (within 1%);
  * fundamental torsional period c/4L (within 0.2%), with less than 1%
    numerical decay;
  * a string at rest stays exactly at rest;
  * steady rotation, pick-up and slack-off match torque and drag (within 2%);
  * in the default decoupled model, rotation does not change pick-up drag;
  * the opt-in along-string coupling matches velocity-split torque and drag;
  * off-bottom stick-slip appears only when μ_s > μ_k;
  * ZTorque suppresses the limit cycle that a stiff controller sustains;
  * on-bottom steady state satisfies the ROP law and the torque balance;
  * the bit breaks out at the sticking torque, sliding TOB is linear in WOB,
    and the bit drills off WOB at a fixed block position;
  * the bit lifts off when picked up;
  * pause and resume are exact;
  * an external Python PI reproduces the internal controller.

## References

* Aarsnes, U.J.F. and Shor, R.J., 2018. Torsional vibrations with bit off bottom: Modeling, characterization and field data validation. *J. Petroleum Science and Engineering*, 163, 712–721.
* Aarsnes, U.J.F., di Meglio, F. and Shor, R.J., 2018. Benchmarking of industrial stick-slip mitigation controllers. *IFAC-PapersOnLine*, 51(8), 233–238.
* Aarsnes, U.J.F., Auriol, J., Di Meglio, F. and Shor, R.J., 2019. Estimating friction factors while drilling. *J. Petroleum Science and Engineering*, 179, 80–91.
* Auriol, J., Aarsnes, U.J.F. and Shor, R.J., 2020. Self-tuning torsional drilling model for real-time applications. *ACC 2020*, 3091–3096.
* Auriol, J., Boussaada, I., Shor, R.J., Mounier, H. and Niculescu, S.-I., 2022. Comparing advanced control strategies to eliminate stick-slip oscillations in drillstrings. *IEEE Access*, 10, 10949–10969.
* Johancsik, C.A., Friesen, D.B. and Dawson, R., 1984. Torque and drag in directional wells — prediction and measurement. *JPT*, 36, 987–992.
