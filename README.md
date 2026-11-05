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

**The two directions are cross-coupled in two places.**

*Along the string*, friction is a single 2-D Coulomb force on the sliding
velocity $\mathbf V = (v,\ \omega r_c)$:

$$
\boldsymbol{\mathcal F} = \mu_k F_N\,\frac{\mathbf V}{\|\mathbf V\|}\ \ \text{(sliding)},
\qquad
\|\boldsymbol{\mathcal F}\| \le \mu_s F_N\ \ \text{(sticking)},
\qquad
\mathcal F_t = r_c\,\mathcal F_\theta,
$$

and the soft-string normal force $F_N$ follows the live tension. As a result:
* rotation reduces drag, and axial motion reduces rotary torque (velocity
  split);
* axial load changes torque through $F_N$;
* combined axial and torsional loads overcome static friction sooner than either alone.

*At the bit*, the characteristic relations

$$
\omega_b = \alpha_t - \text{TOB}/Z_t, \qquad v_b = \alpha_a - \text{WOB}/Z_a,
$$

are closed by one of three reduced bit–rock laws. Each has two or three
parameters and was derived from the
[Bit Model Library](https://github.com/Kartikay100/Bit-Model-Library); see
[`docs/derivation/bit-models.md`](docs/derivation/bit-models.md). With
$a = D/2$ and the drilling strength $S = \zeta\varepsilon$:

| model | parameters | torque on bit | penetration |
|---|---|---|---|
| Coulomb (friction-dominated) | $\mu_s, \mu_k, S$ | sticking $\le \mu_s\,\text{WOB}\,D/3$, sliding $\mu_k\,\text{WOB}\,D/3$ | $\text{ROP} = \omega\,\text{WOB}/(2\pi a S)$ |
| Cutting (Detournay & Defourny) | $S, \zeta, W_f$ | $\text{WOB}\,a/(2\zeta)$ | $\text{ROP} = \omega\,(\text{WOB}-W_f)^+/(2\pi a S)$ |
| Regenerative cutting | $S, \zeta$ + blades | $\text{WOB}\,a/(2\zeta)$ | depth of cut from the hole bottom left by the previous blade |

In every model TOB is proportional to WOB, and penetration needs rotation.

**Rock strength, pore pressure and temperature.** With
`strength_source='formation'` the drilling strength is taken from the rock at
the bit: $S = \zeta\,\text{CCS}_\text{eff}/\eta$. The effective CCS combines
three terms:
* Mohr–Coulomb confinement by the overbalance. Low-porosity,
  low-permeability rock such as granite uses an undrained (Skempton) pore
  pressure (Caicedo, Calhoun & Ewy 2005).
* Thermal weakening above a threshold, fitted to granite data.
* Weakening of the cooled skin at the hole bottom by thermoelastic stress.

The mud temperature at the bit, and the bottomhole pressure, come from the
growing-wellbore drift-flux model of **Wang, Khaled, Luu, Ashok & van Oort
(2023), SPE-214836-MS**, implemented for single-phase mud. It reproduces their
Utah FORGE case: downhole temperature falls from about 310 °F to 210 °F in
30 min at 550 gpm, and drilling adds 1.3–2.2 °F, linear in ROP. See
[`docs/derivation/rock-strength-and-temperature.md`](docs/derivation/rock-strength-and-temperature.md).
The regenerative model adds the delay of one blade passage, which can
self-excite bit bounce and stick-slip.

Static friction is anisotropic. The torsional static/kinetic contrast of
Aarsnes & Shor applies to rotation. Axially, the static factor defaults to the
kinetic one (`axial_static_contrast=0`). Applying the torsional contrast
axially produces artificial axial stick-slip when tripping; see the derivation,
§5.5.

For comparison, `along_string_coupling=False` uses independent 1-D axial and
torsional friction with a static normal force, leaving the bit as the only
coupling. With no axial motion, both friction models reduce exactly to
Aarsnes & Shor (2018).

Other model components:
* **Top drive:** an inertia driven by a stiff, SoftTorque or ZTorque
  controller, by an external torque, or by a prescribed speed.
* **Drawworks:** a prescribed block velocity, or an autodriller. The hook can
  only hold the string up, so hookload is never negative.
* **Wellbore temperature and pressure:** the drift-flux model of Wang et al.
  (2023), coupled every 2 s of simulated time.
* **Gravity:** handled by a well-balanced static/dynamic split.
* **Numerics:** a Godunov finite-volume scheme (upwind on the invariants,
  exact at section changes) with a CFL-limited time step. Coulomb friction is
  applied implicitly (backward Euler), so it opposes the sliding direction
  exactly and sticking is exact.

The derivation, boundary conditions, discrete scheme and analytic limits are
in [`docs/derivation/README.md`](docs/derivation/README.md).

## Installation

```bash
pip install -e ".[test]"
pytest -q
```

## Interactive simulator

```bash
pip install -e ".[ui]"
drillpyze-ui            # opens http://127.0.0.1:8050 in your browser
```

The UI starts with the Well A example case preloaded: the Aarsnes & Shor
geometry in FORGE-like granite, with the wellbore temperature model on. The **SI / Field**
switch at the top of the input panel changes every input, readout, plot axis
and CSV import between SI (m, mm, kN, kN·m, kg/m³, m/h) and field units (ft,
in, klbf, kft·lbf, ppg, ft/h). The case itself is always stored in SI, so
switching never changes it, and the choice is remembered in the browser.

* **Edit the case** in the left panel: survey, drillstring, bit position,
  string and bit friction, damping, top-drive controller and numerics.
  * *Apply and reset* restarts the simulation with the new inputs.
  * The survey and drillstring tables can be imported from CSV files.
  * The whole case can be saved and reopened as JSON (`DrillPyze.config`).
* **Play, pause and step** through the simulation.
  * Choose the frame length and the playback speed (from ¼× real time to as
    fast as possible).
  * *Run ahead* jumps forward without animating each frame. Space toggles
    play; the right arrow steps one frame.
* **Change rig inputs while running** in the transport bar:
  * top drive: controller set point, or a prescribed speed;
  * drawworks: block speed, or an autodriller WOB set point;
  * circulation: pump rate and mud inlet temperature.
* **Watch the string**:
  * the well path, with the casing, drill pipe, BHA and bit positions;
  * rpm, torque and tension against depth, compared with torque-and-drag
    references;
  * drillstring, annulus and formation temperature against depth;
  * time histories with 60 s, 10 min, 60 min and full-run windows, including
    mud and rock temperature at the bit and the effective CCS.
* **Scrub back** through earlier frames with the replay slider.
* **Download** the recorded channels as CSV.

The server runs the same `Simulation` objects as the Python API (see
`DrillPyze/ui/server.py`). Plotly is served from the installed `plotly`
package, so the UI also works offline.

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

| cell size | cells | time step | coupled friction (default) | decoupled option |
|---|---|---|---|---|
| 10 m  | 334   | 1.85 ms | 85× real time (~42 s per 60 min) | 200× (~18 s) |
| 5 m   | 668   | 0.92 ms | 22× (~2.7 min) | 49× (~1.2 min) |
| 2.5 m | 1,335 | 0.46 ms | 6× (~10 min) | 12× (~5 min) |

The coupled model costs about 2× more because it solves an implicit 2-D
friction update (a short Newton iteration) in every cell, every step.

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
  * rotating while picking up or slacking off matches velocity-split torque
    and drag (within 2%);
  * rotation reduces pick-up drag with coupled friction, and has no effect with
    the decoupled option;
  * the implicit 2-D friction is exact, including combined-load breakout;
  * off-bottom stick-slip appears only when μ_s > μ_k;
  * ZTorque suppresses the limit cycle that a stiff controller sustains;
  * on-bottom steady state satisfies the ROP law and the torque balance;
  * the bit breaks out at each model's sticking limit, and the Coulomb bit's
    sliding TOB is linear in WOB;
  * steady drilling matches each bit law; the regenerative model converges
    to the cutting state when damped and self-excites when lightly damped;
  * every bit model drills off WOB at a fixed block position;
* Rock strength and temperature:
  * drained and undrained CCS, the porosity blend, the thermal-weakening fit,
    and the thermal stress and skin ratio;
  * a static well stays at the geothermal profile;
  * the energy balance closes within 3% while circulating;
  * the SPE-214836 FORGE cooling curve and the linear rise of downhole
    temperature with ROP are reproduced;
  * ROP → 0 converges to the static wellbore;
  * the formation-driven drilling strength reproduces the cutting law;
  * pause and resume stay exact with the slow models running;
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
