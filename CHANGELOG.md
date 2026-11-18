# Changelog

## 1.1.0 (2026-11-18)

Adds energy accounting to the Python model and DrillPyze Lite (physics spec revision 4; derivation §11): where the power goes, in real time, and how much a hybrid rig could recover by regenerative braking. The accumulators change no state, so every case gives the same results as in 1.0.0.

**Energy accounting.**
- Every time step accumulates the work of the top-drive motor (driving and braking), the power into the string at surface, the hook (hoisting and braking), gravity on the moving string, the bit (torsional and axial), the mud motor, and the losses to wall friction and damping (from the kinetic energy each implicit sub-step removes). The balances close to the scheme's numerical dissipation, about 0.6 % for Well A.
- The hydraulics split the standpipe pressure into surface lines, friction inside the string, the mud motor, the bit nozzles, friction in the annulus and the remainder (`WellboreHydraulics.loss_split`); each times the flow is a power.
- `Simulation.energy()` returns the cumulative energies (J) and the energy stored in the string and the top drive; `DrillPyze.energy.flows(e0, e1)` gives mean powers between two snapshots, with rock destruction, the motor's losses, the pump drive, braking and the recoverable part.
- New case section `parameters.energy` (`EnergyParameters`): pump efficiency and regeneration efficiency, 0.9 each by default, editable under *Energy* in the case panel.
- `examples/tripping_energy.py`: tripping Well A to 2500 m brakes 255 kWh, of which 230 kWh is recoverable at 90 %; the trip out lifts 280 kWh.

**Simulator** (Python UI and DrillPyze Lite).
- An energy panel under the drilling log: a Sankey diagram from the Mud Pumps, top drive, drawworks and gravity to Rock Destruction, Frictional Losses inside the string, Frictional Losses in the Annulus, the nozzles, wall friction, damping and recoverable braking, with the mud motor branching off the mud circulation; summary tiles (MSE, recoverable braking, energy since the start in kWh); and power over the same span as the log.
- 1, 10 or 60 s means; *Everything* or the drillstring alone; kW in SI, hp in field units.
- On a phone the losses are grouped so every label fits, and demo mode has a fourth page, *Energy*.

## 1.0.0 (2026-11-16)

First stable release, archived on Zenodo: https://doi.org/10.5281/zenodo.24289880. Since 1.0.0a1 it adds a mud motor (physics spec revision 3), a phone layout and a full-screen demo mode for the simulator, anonymous visit counts on the hosted DrillPyze Lite, and documentation fixes. Cases without a motor give the same results as in 1.0.0a1.

**Mud motor** (physics spec revision 3; derivation §6.4; Python model and DrillPyze Lite).
- An optional positive-displacement motor between the string and the bit, after Nguyen et al. (2017, AADE-17-NTCE-026). Well A defaults to no motor; ticking *Mud motor in the BHA* adds a generic 6¾" 7/8-lobe motor (0.28 rev/gal) whose values can be edited.
- **Speed:** motor constant (rev/gal) × flow, with optional slip at the rated differential and bypass flow.
- **Torque:** T = K_T (Δp − Δp₀) with K_T = η / (2π k_n), so the differential reads the bit torque.
- **Stall:** above the stall differential the bit stops, the motor passes its stall torque, and SPP jumps; picking up clears it.
- **Reactive torque and toolface:** the bit torque acts back on the string, so while sliding the toolface turns left as WOB rises.
- **Hydraulics:** the motor differential is in series with the bit nozzles, so it adds to SPP but not to BHP.
- **Driller's display:**
  - a *BHA: Mud motor / No motor* flag in the status strip and a *Motor* badge on the rotary dial;
  - a toolface dial with motor rpm, differential pressure, a stall warning, and ⟲ / ⟳ buttons that orient the toolface 15° while sliding;
  - an off-bottom pressure marker on the standpipe dial.
- **Drilling log:** motor Δp on the hydraulics track and toolface on the downhole track.
- New outputs: `omega_motor`, `dp_motor`, `motor_stalled` and `toolface` in `Simulation.latest()`; `omega_motor`, `dp_motor` and `stator_angle` in the recorded channels.

**Documentation.**
- README references now include Wang et al. (2023) (drift-flux temperature and hydraulic transients), the *Introduction to Drilling Engineering* hydraulics chapter, Detournay & Defourny (1992), Richard, Germay & Detournay (2007) and Caicedo, Calhoun & Ewy (2005).
- Fixed the typesetting of equation (5.2) in the model derivation ("\tag not allowed in split environment").

**Simulator UI on phones** (shared by the Python UI and DrillPyze Lite).
- The inputs fold away behind an *Inputs* button under a compact header, and *Apply* returns to the simulator.
- Transport, rig status and the driller's display are laid out for a narrow screen; the page scrolls over the dials, and set points still drag.
- The along-string views show one at a time, picked with tabs.
- The drilling log shows two tracks at a time; change them with the arrows or by swiping.
- No plot spills past the screen edge, and plots no longer capture the page's scrolling.

**Demo mode.** *Demo* in the transport bar, or a link ending in `#demo`, goes full screen (where the browser allows):
- The driller's display fills the first page, scaled to fit any screen.
- The along-string views and the log are further pages: swipe up or left, scroll, use the arrow keys, or tap *Dials · String · Log*.
- A slim bar keeps Play, the clock and the simulation speed in reach.
- Esc or × exits.

**Visit counts for DrillPyze Lite.**
- The hosted site counts anonymous visits with GoatCounter, which uses no cookies, and says so under the title.
- It also counts a few events: first Play, run-ahead, demo mode, a new case and a results download.
- The Python UI is unchanged. `lite/build_site.py --no-analytics` leaves the counter out.

## 1.0.0a1 — 1.0 alpha (2026-11-13)

First alpha of DrillPyze 1.0.

**Coupled model.** A real-time model of coupled axial–torsional drillstring dynamics:
- soft-string friction along the well path (coupled or decoupled);
- three bit models (Coulomb, cutting, regenerative);
- fluid and material damping up to 10 Hz;
- top-drive controllers (stiff, SoftTorque, Z-Torque, custom) and a PI autodriller;
- ramped set points;
- circulating hydraulics (steady losses and a drift-flux transient model);
- a wellbore temperature model;
- rock strength at the bit.

**Simulator UI.**
- A driller's display with draggable set points.
- Pason-style log tracks with clock time, shared scales and per-curve limits.
- Along-string views.
- A case editor.
- SI and field units.

**DrillPyze Lite.** The same simulator running entirely in a web browser
(`lite/`), built from the physics specification in `docs/spec`.

**Documentation.** Sphinx API reference, the model derivations and the
verification cases.
