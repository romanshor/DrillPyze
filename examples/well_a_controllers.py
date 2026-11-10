"""
Off-bottom start-up in the synthetic Well A of Aarsnes & Shor (2018) with the
three industrial top-drive controllers of Aarsnes, di Meglio & Shor (2018).

Run:  python examples/well_a_controllers.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from DrillPyze import (Survey, Component, Drillstring, Wellbore, Fluid, Simulation,
                       ModelParameters, FrictionParameters, DampingParameters, RampParameters,
                       units, controllers)

# --- 1. import / describe the well and string (SI) ------------------------
survey = Survey(md=[0, 1500, 2100, 3000], inc=[0, 0, 60, 60], azi=[0, 0, 0, 0])
well = Wellbore(survey, casing_shoe=1000.0, hole_diameter=8.5 * units.IN)
string = Drillstring([                                   # bit first
    Component('Collars', 230.0, od=5.75 * units.IN, id=2.8125 * units.IN),
    Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN),
])
mud = Fluid(density=1200.0)

params = ModelParameters(
    friction=FrictionParameters(mu_static_cased=0.40, mu_kinetic_cased=0.22,
                                mu_static_open=0.45, mu_kinetic_open=0.25, omega_c=1.0),
    damping=DampingParameters(k_t=0.1, k_a=0.1),
    # a fast (1 s) start-up ramp excites off-bottom stick-slip for the controllers to act on;
    # with the default 5 s ramp the string breaks free gently and the stiff controller stays at 60 rpm
    ramps=RampParameters(omega_sp=1.0),
)

# first torsional mode for SoftTorque/ZTorque tuning: f ~ c / (4 L)
f1 = 3178.0 / (4 * 2500.0)

results = {}
for name in ('stiff', 'soft_torque', 'z_torque'):
    sim = Simulation(well, string, mud, params.copy())
    sim.initialize(bit_depth=2500.0, hole_depth=2505.0, dx=10.0)   # --- 2. set up / initialise
    I = sim.params.top_drive.inertia
    sim.params.controller = {
        'stiff': controllers.stiff(sim.zeta_p, I),
        'soft_torque': controllers.soft_torque(sim.zeta_p, I, f1),
        'z_torque': controllers.z_torque(sim.zeta_p, I),
    }[name]
    sim.update_parameters()
    sim.step(2.0, omega_sp=0.0)                                      # --- 3. step
    sim.step(88.0, omega_sp=60 * units.RPM)
    results[name] = sim.results()

fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
for name, df in results.items():
    ax[0].plot(df.time, units.rads_to_rpm(df.omega_bit), lw=1, label=name)
    ax[1].plot(df.time, df.motor_torque / 1e3, lw=1, label=name)
ax[0].set_ylabel('Bit speed (rpm)')
ax[1].set_ylabel('Motor torque (kN m)')
ax[1].set_xlabel('Time (s)')
ax[0].legend()
fig.suptitle('Well A off-bottom start-up, 60 rpm set point')
fig.tight_layout()
fig.savefig('well_a_controllers.png', dpi=130)
for name, df in results.items():
    late = df[df.time > 60]
    rpm = units.rads_to_rpm(late.omega_bit)
    print(f'{name:12s} bit rpm last 30 s: min {rpm.min():6.1f}  max {rpm.max():6.1f}')
