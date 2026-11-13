# DrillPyze — coupled drillstring dynamics

[![Tests](https://github.com/romanshor/DrillPyze/actions/workflows/tests.yml/badge.svg)](https://github.com/romanshor/DrillPyze/actions/workflows/tests.yml)
[![Documentation](https://github.com/romanshor/DrillPyze/actions/workflows/main.yml/badge.svg)](https://romanshor.github.io/DrillPyze/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Version 1.0 alpha (1.0.0a1, released 2026-11-13).**
[Documentation](https://romanshor.github.io/DrillPyze/) ·
[Try it in your browser (DrillPyze Lite)](https://romanshor.github.io/DrillPyze/lite/) ·
[Changelog](changelog.md)

DrillPyze simulates the coupled axial and torsional dynamics of a drillstring
in a 3-D wellbore, together with steady-state torque and drag. It is built to:

* calibrate along-string friction factors from field data,
* study stick-slip in its off-bottom, on-bottom and axially driven forms, and
* develop and benchmark top-drive RPM controllers.

All quantities are **SI** internally. Conversion helpers in `DrillPyze.units`
handle oilfield units at input and output.

DrillPyze 1.0 is a rewrite of the 2023 prototype. The full mathematical
derivation is in [`docs/derivation/`](derivation-readme.md).

## Mathematical model

The string is a slender elastic rod with $x \in [0, L]$ measured from the top
drive to the bit. Torsion follows Aarsnes & Shor (2018). The axial direction
uses the **identical** formulation, with Young's modulus $E$ and area $A$ in
place of $G$ and $J$:

```{math}
\begin{aligned}
&\partial_t \tau + JG\,\partial_x\omega = 0, &\qquad
&\rho J\,\partial_t\omega + \partial_x\tau = -k_t\rho J\omega - \mathcal F_t,\\
&\partial_t q + EA\,\partial_x v = 0, &\qquad
&\rho A\,\partial_t v + \partial_x q = w_b\cos\theta - k_a\rho A v - \mathcal F_a,
\end{aligned}
```

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

```{math}
\boldsymbol{\mathcal F} = \mu_k F_N\,\frac{\mathbf V}{\|\mathbf V\|}\ \ \text{(sliding)},
\qquad
\|\boldsymbol{\mathcal F}\| \le \mu_s F_N\ \ \text{(sticking)},
\qquad
\mathcal F_t = r_c\,\mathcal F_\theta,
```

and the soft-string normal force $F_N$ follows the live tension. As a result:
* rotation reduces drag, and axial motion reduces rotary torque (velocity
  split);
* axial load changes torque through $F_N$;
* combined axial and torsional loads overcome static friction sooner than either alone.

*At the bit*, the characteristic relations

```{math}
\omega_b = \alpha_t - \text{TOB}/Z_t, \qquad v_b = \alpha_a - \text{WOB}/Z_a,
```

are closed by one of three reduced bit–rock laws. Each has two or three
parameters and was derived from the
[Bit Model Library](https://github.com/Kartikay100/Bit-Model-Library); see
[`docs/derivation/bit-models.md`](derivation-bit-models.md). With
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
[`docs/derivation/rock-strength-and-temperature.md`](derivation-rock-strength-and-temperature.md).
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
* **Wellbore temperature:** the drift-flux model of Wang et al. (2023),
  coupled every 2 s of simulated time.
* **Hydraulics** ([`docs/derivation/hydraulics.md`](derivation-hydraulics.md)):
  standpipe pressure, bottomhole pressure and ECD.
  * **Steady losses:** Newtonian or Bingham pipe and annular losses, bit
    nozzle loss and the SPP sum, from the borehole-hydraulics chapter of
    *Introduction to Drilling Engineering* (ch. 7).
  * **Transients:** pump start and stop, water hammer, surge and swab, from
    the single-phase drift-flux model of Wang et al. (2023). It uses fixed
    cells with volume-averaged sections, ρ(T, p), a pump inlet, a pressure
    outlet and the bit nozzles between the two flow paths.
  * **Coupling:** updated every ~1 s. Its BHP feeds the rock-strength model.
* **Set-point ramps:** inputs without a controller approach a new value
  linearly over a short ramp time (`ModelParameters.ramps`). Defaults:
  prescribed top-drive speed 1 s, block speed 2 s, pump rate and inlet
  temperature 3 s. Set points that feed a controller (RPM, autodriller WOB)
  step, and the controller shapes the response.
* **Damping** (derivation §5.6): no tuning knob is needed; the coefficients
  come from the mud report and the string.
  * **F1, mud:** steady laminar Couette drag on the full velocity, plus a
    Stokes boundary layer $\sqrt{\rho_f\mu\omega_1/2}$ on the inner and outer
    pipe walls. The Stokes term acts on the oscillatory part of the velocity
    and is evaluated at the first axial and torsional modes.
  * **F2, mud yield stress:** a Bingham yield drag $\tau_y\pi d_o$ added to
    the Coulomb capacity, by default only with the pumps off.
  * **M1, steel:** Rayleigh / Kelvin–Voigt damping fitted to $\zeta = \eta/2$
    at the first mode and at 10 Hz, solved implicitly.
  * **Calibration:** $k_t, k_a$ remain as optional extra rates.
  * **Limitation:** the model resolves and damps vibrations **up to 10 Hz**.
    Above the first mode, fluid damping is under-predicted by
    $\sqrt{f/f_1}$. For torsional content above about 2 Hz, the grid
    spacing should be reduced (about 2.5 m for 10 Hz) so numerical
    dissipation stays below the physical damping.
* **Gravity:** handled by a well-balanced static/dynamic split.
* **Numerics:** a Godunov finite-volume scheme (upwind on the invariants,
  exact at section changes) with a CFL-limited time step. Coulomb friction is
  applied implicitly (backward Euler), so it opposes the sliding direction
  exactly and sticking is exact.

The derivation, boundary conditions, discrete scheme and analytic limits are
in [`docs/derivation/README.md`](derivation-readme.md).

## Installation

```bash
git clone https://github.com/romanshor/DrillPyze.git
cd DrillPyze
pip install -e ".[ui,test]"     # ui: the interactive simulator; test: pytest
pytest -q                       # optional: run the verification suite
python -m DrillPyze.ui          # start the simulator in your browser
```

DrillPyze needs Python 3.10 or newer. To run it without installing anything,
use [DrillPyze Lite](https://romanshor.github.io/DrillPyze/lite/).

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
* **Drive the rig from the driller's display** at the top of the workspace.
  It has six dials (WOB, rotary speed, torque, ROP, standpipe pressure,
  bottomhole pressure) and two bars (pump rate, mud weight / ECD).
  * Each dial's needle and large number are the surface measurement. The
    model's downhole estimate is an orange arrow outside the dial, and a band
    shows its range over the last 10 s, so stick-slip is visible as a wide
    band on the rotary dial.
  * Bottomhole pressure (a model value, orange needle) shows the static
    mud-weight pressure as a grey tick. The mud-weight bar shows static mud
    weight as the bar and ECD as an arrow. Both show the pore- and
    fracture-pressure window.
  * Set points are cyan diamonds on the dials themselves: drag one, or focus
    it and use the arrow keys. They cover the top-drive set point (or
    prescribed speed), the autodriller WOB or block speed (the chips choose
    which), and the pump rate. The mud inlet temperature sits under the pump bar.
  * Surface WOB comes from the weight indicator: off-bottom hookload minus
    hookload. *Zero* takes the present hookload as off-bottom. Choosing the
    autodriller while off bottom zeroes it automatically
    (`Simulation.tare_hookload`).
  * Limits are red markers and zones, set in the *Driller's display* section
    of the inputs (`display.*` in the case). They change the display only.
  * The status strip shows bit and hole depth, hookload, mud temperature at
    the bit, rock strength, the drawworks mode and a stick-slip alarm.
* **Watch the string**:
  * the well path, with the casing, drill pipe, BHA and bit positions. The open
    hole is drawn to the current hole depth; any survey beyond it is a light
    dashed line;
  * rpm, torque and tension against depth, compared with torque-and-drag
    references;
  * drillstring, annulus and formation temperature against depth;
  * time histories as a drilling-log (EDR-style) display, with seven tracks
    side by side:
    * drilling: block height, ROP, bit and hole depth;
    * rotary speed;
    * downhole: stick-slip index over 10 s and bit angular acceleration;
    * torque;
    * axial loads;
    * hydraulics: SPP, ECD, flow in and out;
    * bit and rock: temperatures, effective CCS.

    The tracks share a time axis that increases downward, with the newest
    data at the bottom; 60 s, 10 min, 60 min and full-run windows. The axis
    shows clock time when the case has a start date/time
    (`simulation.start_time`), and elapsed seconds otherwise;
  * per-curve scales: each curve header shows its name, current value and
    scale ends. Click a header to fix that curve's limits or return it to
    automatic. Fixed limits are marked with a lock and remembered in the
    browser, separately for SI and field units. Curves of the same quantity
    share one scale, and setting it on either sets both: bit and hole depth,
    RPM set point and top drive, motor and surface-pipe torque, flow in and
    out, and rock and mud temperature at the bit.
* **Scrub back** through earlier frames with the replay slider.
* **Download** the recorded channels as CSV.

The server runs the same `Simulation` objects as the Python API (see
`DrillPyze/ui/server.py`). Plotly is served from the installed `plotly`
package, so the UI also works offline.

## Physics specification

`docs/spec/physics-spec.md` specifies every equation, algorithm and default
of the model at implementation level: geometry, the time step, inputs and
ramps, hydraulics, rock strength, temperature and outputs. It is the
reference for independent implementations, including DrillPyze Lite.

**Keep it current.** A change to the physics or numerics updates the spec in
the same commit, raises its revision and refreshes the manifest:

```bash
python docs/spec/update_manifest.py
```

`tests/test_spec_sync.py` fails when a physics source file changes without
this step. It also fails when the browser engine implements an older spec
revision.

## DrillPyze Lite (browser version)

`lite/` runs the same simulator entirely in a web browser, with no install.
- **Engine:** a JavaScript engine written from the specification runs in a
  Web Worker behind the unchanged front end.
- **Physics:** the whole model, including the wellbore temperature model.
- **Speed:** about 20–30× real time for Well A.
- **Conformance:** `tests/test_lite_engine.py` (Node.js) checks that it
  reproduces the Python model to round-off, with and without the
  temperature model.

```bash
python lite/build_site.py && python -m http.server -d _site 8000   # http://localhost:8000
```

**Publishing on GitHub Pages.** The workflow
`.github/workflows/pages-lite.yml` runs the conformance tests, builds the
site and publishes it to the `lite/` folder of the `gh-pages` branch. It
runs on pushes to `main` that touch `lite/` or the front
end, or by hand from the Actions tab.

The site appears at `https://<owner>.github.io/<repo>/lite/`, next to the
documentation. The one-time setting is *Settings → Pages → Build and
deployment → Source: Deploy from a branch, `gh-pages` / `(root)`*. The
documentation workflow keeps the `lite/` folder when it publishes.

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
sim.step(60.0, flow_rate=0.03)               # pumps ramp to 0.03 m3/s over params.ramps.flow_rate
df = sim.results()                           # recorded channels (DataFrame, SI)
hyd = sim.hydraulics_results()               # SPP, BHP, ECD, flow in / out, once per ~1 s
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
* Damping:
  * measured free-decay damping ratios of the first axial and torsional modes
    match the predictions within 0.5% (lumped, fluid, material, combined);
  * steady rotation sees only the Couette drag;
  * the yield drag applies with the pumps off only;
  * the Kelvin–Voigt solve conserves momentum;
  * pause and resume are exact.
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
* Hydraulics:
  * chapter 7 Example 7.1 and Practice Problems 7.1 and 7.2: velocities,
    apparent viscosities, Reynolds numbers, regimes, section losses, bit loss
    and ECD reproduced (three printed turbulent losses in the chapter are
    arithmetic errata; see the hydraulics derivation);
  * the transient model's discrete steady state matches the textbook SPP and
    BHP;
  * a pump start converges to it, and mass is conserved to round-off;
  * the Joukowsky surge is within 1% and the acoustic travel time exact;
  * surge and swab have the right sign and are symmetric;
  * set-point ramps are linear in time;
  * pause and resume stay exact with hydraulics running;
* General:
  * pause and resume are exact;
  * an external Python PI reproduces the internal controller.

## Documentation

The API reference (from the docstrings), this README and the model
derivations are built with Sphinx:

```bash
./docs_generator.sh          # or: sphinx-build -b html docs/source docs/build
```

The workflow `.github/workflows/docs-coupled.yml` builds them on every push to
`coupled-model` that changes the code or the docs. It treats Sphinx warnings as
errors and publishes the pages to
`https://<owner>.github.io/<repo>/coupled-model/`, next to the `main`
documentation.

## Citing DrillPyze

If you use DrillPyze in your work, please cite it. The citation details are
in [`CITATION.cff`](../../../CITATION.cff); GitHub shows them under *Cite this
repository*.

## Contributing and licence

Contributions are welcome: see [`CONTRIBUTING.md`](../../../CONTRIBUTING.md). DrillPyze
is released under the [MIT License](../../../LICENSE).

## References

* Aarsnes, U.J.F. and Shor, R.J., 2018. Torsional vibrations with bit off bottom: Modeling, characterization and field data validation. *J. Petroleum Science and Engineering*, 163, 712–721.
* Aarsnes, U.J.F., di Meglio, F. and Shor, R.J., 2018. Benchmarking of industrial stick-slip mitigation controllers. *IFAC-PapersOnLine*, 51(8), 233–238.
* Aarsnes, U.J.F., Auriol, J., Di Meglio, F. and Shor, R.J., 2019. Estimating friction factors while drilling. *J. Petroleum Science and Engineering*, 179, 80–91.
* Auriol, J., Aarsnes, U.J.F. and Shor, R.J., 2020. Self-tuning torsional drilling model for real-time applications. *ACC 2020*, 3091–3096.
* Auriol, J., Boussaada, I., Shor, R.J., Mounier, H. and Niculescu, S.-I., 2022. Comparing advanced control strategies to eliminate stick-slip oscillations in drillstrings. *IEEE Access*, 10, 10949–10969.
* Johancsik, C.A., Friesen, D.B. and Dawson, R., 1984. Torque and drag in directional wells — prediction and measurement. *JPT*, 36, 987–992.
