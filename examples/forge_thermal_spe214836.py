"""
Reproduce the growing-wellbore temperature study of Wang, Khaled, Luu, Ashok &
van Oort (2023), SPE-214836-MS (Figs. 6-7): Utah FORGE 16(78)-32 near 5,900 ft
MD, 550 gpm, geothermal 40 C + 69.1 C/km, mud and rock properties from their
Table 2.  The mud inlet temperature is not stated in the paper; 35 C is
assumed here.

Run:  python examples/forge_thermal_spe214836.py
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from DrillPyze import Survey, Component, Drillstring, Wellbore, Fluid, units
from DrillPyze.formation import Formation
from DrillPyze.thermal import WellboreThermal, ThermalParameters

FT = units.FT
F = lambda c: np.asarray(c) * 9 / 5 + 32

survey = Survey([0, 5000 * FT, 5900 * FT, 6200 * FT], [0, 0.8, 18, 18], [0, 0, 0, 0])
well = Wellbore(survey, casing_shoe=4000 * FT, hole_diameter=0.22, casing_id=0.22)
string = Drillstring([Component('Bit+BHA', 0.8, od=0.20, id=0.07), Component('DP', 100.0, od=0.127, id=0.11)])
mud = Fluid(density=1078, viscosity=0.020, heat_capacity=3750, thermal_conductivity=0.75, thermal_expansion=1e-4)
rock = Formation()      # FORGE defaults: 40 C, 69.1 C/km, 2800 kg/m3, 930 J/kg K, 2.31 W/m K


def run(rop_fth, minutes=30.0, dt=10.0):
    th = WellboreThermal(well, string, mud, rock, ThermalParameters(flow_rate=550 * units.GPM, inlet_temperature=35.0),
                         bit_depth=5900 * FT, hole_depth=5900 * FT)
    hole, rop = 5900 * FT, rop_fth * FT / 3600
    t, dht, out = [0.0], [th.state()['t_bit']], [th.state()['t_out']]
    for _ in range(int(minutes * 60 / dt)):
        hole += rop * dt
        s = th.advance(dt, bit_depth=hole, hole_depth=hole, rop=rop)
        t.append(th.time / 60); dht.append(s['t_bit']); out.append(s['t_out'])
    return np.array(t), np.array(dht), np.array(out)


runs = {r: run(r) for r in (0, 40, 60, 72)}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
t, d0, o0 = runs[0]
ax[0].plot(t, F(d0), label='Downhole, not drilling')
ax[0].plot(t, F(runs[72][1]), '--', label='Downhole, ROP 72 ft/h')
ax[0].plot(t, F(o0), label='Surface outlet, not drilling')
ax[0].set(xlabel='Time (min)', ylabel='Temperature (°F)', title='Mud temperature (cf. SPE-214836 Fig. 6)')
ax[0].legend()
for r, paper in ((40, 1.3), (60, 2.2), (72, 2.7)):
    ax[1].plot(t, (F(runs[r][1]) - F(d0)), label=f'ROP {r} ft/h (paper ≈ {paper} °F at 30 min)')
ax[1].set(xlabel='Time (min)', ylabel='DHT increase from drilling (°F)', title='Effect of ROP (cf. Fig. 7)')
ax[1].legend()
fig.tight_layout()
fig.savefig('forge_thermal_spe214836.png', dpi=130)
for r in (40, 60, 72):
    print(f'ROP {r} ft/h: DHT increase after 30 min = {F(runs[r][1][-1]) - F(d0[-1]):.2f} °F')
print(f'DHT: {F(d0[0]):.0f} °F at pump start, {F(d0[-1]):.0f} °F after 30 min')
