# DrillPyze — coupled drillstring dynamics

DrillPyze simulates the coupled axial and torsional dynamics of a drillstring
in a 3-D wellbore, together with steady-state torque and drag. It is built to:

* calibrate along-string friction factors from field data,
* study stick-slip in its off-bottom, on-bottom and axially driven forms, and
* develop and benchmark top-drive RPM controllers.

All quantities are **SI** internally. Conversion helpers in `DrillPyze.units`
handle oilfield units at input and output.

> **Branch `coupled-model`:** a rewrite of the 2023 prototype. Model equations
> and numerical choices are documented in [`docs/model.md`](docs/model.md).

## Model at a glance

* **Torsional and axial wave equations.** Solved with a Godunov (exact Riemann)
  finite-volume scheme, which is first-order upwind on the Riemann invariants.
  Impedance changes between pipe, HWDP and collars are handled exactly.
  Torsion follows Aarsnes & Shor (2018); the same scheme is extended to the
  axial direction.
* **Along-string Coulomb friction.** Static and kinetic factors for cased and
  open hole. Friction acts opposite to the combined axial and tangential
  sliding velocity, and the normal force follows the live tension (soft
  string). This is the source of the axial–torsional coupling: rotation reduces
  drag, and axial loads change torque.
* **Top drive.** Motor inertia with three ways to drive it:
  * an internal stiff, SoftTorque or ZTorque controller
    (Aarsnes, di Meglio & Shor 2018);
  * an external motor torque, for controllers written in Python;
  * a prescribed speed, for replaying measured RPM.
* **Drawworks.** Prescribed block velocity, or an autodriller that holds a WOB
  set point.
* **Bit.** Off bottom, or on bottom with unilateral contact. Torque on bit is
  WOB-dependent and velocity-weakening, with sticking. A ROP law links WOB, RPM
  and ROP.
* **Geometry.** Minimum-curvature surveys. Cell geometry follows the string as
  it moves. Gravity is handled with a well-balanced split, so a hanging string
  at rest stays exactly at rest.

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
  * steady rotation, pick-up, slack-off and combined rotate-and-pick-up all
    match torque and drag (within 2–3%);
  * off-bottom stick-slip appears only when μ_s > μ_k;
  * ZTorque suppresses the limit cycle that a stiff controller sustains;
  * on-bottom steady state satisfies the ROP law and the torque balance;
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
