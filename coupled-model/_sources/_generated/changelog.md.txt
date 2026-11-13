# Changelog

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
(branch `browser-lite`), built from the physics specification in `docs/spec`.

**Documentation.** Sphinx API reference, the model derivations and the
verification cases.
