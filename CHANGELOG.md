# Changelog

## Unreleased

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
