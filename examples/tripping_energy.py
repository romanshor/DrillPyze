"""
Energy a regenerative drawworks could recover while tripping Well A
(docs/derivation/README.md sec. 11).

One 90 ft (27.4 m) stand is run in and pulled out at 0.5 m/s, pumps off and no
rotation, at several bit depths.  Going in, gravity does the work and the drawworks
brakes it: a hybrid rig could store that energy.  Coming out, the drawworks lifts the
string against gravity and wall friction.  The per-stand energies are integrated over
depth for a full trip, and the top drive's braking when 80 rpm rotation is stopped is
shown for comparison.

Run:  python examples/tripping_energy.py
"""
from DrillPyze import config, energy, units

STAND = 27.432            # m, a 90 ft triple
V_TRIP = 0.5              # m/s block speed, about 100 ft/min
DEPTHS = (300.0, 700.0, 1100.0, 1500.0, 1900.0, 2300.0, 2500.0)
KWH = 3.6e6


def stand(depth, direction):
    """Hook energy (J) for one stand ending (in) or starting (out) at `depth`."""
    cfg = config.well_a()
    p = cfg['parameters']
    p['hydraulics']['enabled'] = False
    p['thermal']['enabled'] = False
    p['thermal']['flow_rate'] = 0.0
    start = depth - STAND if direction == 'in' else depth
    cfg['simulation'].update(bit_depth=start, hole_depth=2505.0, block_height=40.0)
    sim = config.build(cfg)
    sim.step(5.0, flow_rate=0.0, omega_sp=0.0)              # settle, the string hanging still
    e0 = sim.energy()
    sim.step(STAND / V_TRIP, v_block=V_TRIP if direction == 'in' else -V_TRIP)
    sim.step(4.0, v_block=0.0)                               # stop and let it settle
    f = energy.flows(e0, sim.energy(), sim.params.energy)
    dt = f['interval']
    return dict(braked=f['hook_neg'] * dt, lifted=f['hook_pos'] * dt, friction=max(0.0, f['friction_axial']) * dt,
                damping=f['damping_axial'] * dt,
                eta=sim.params.energy.regen_efficiency)


def trip_total(per_stand):
    """Integrate the per-metre energy over depth from surface (linear to the first depth)."""
    s, d0, e0 = 0.0, 0.0, 0.0
    for d, e in zip(DEPTHS, per_stand):
        e /= STAND
        s += 0.5 * (e + e0) * (d - d0)
        d0, e0 = d, e
    return s


def stop_rotation(rpm=80.0):
    cfg = config.well_a()
    cfg['parameters']['thermal']['enabled'] = False
    sim = config.build(cfg)
    sim.step(30.0, omega_sp=rpm * units.RPM)
    e0 = sim.energy()
    sim.step(10.0, omega_sp=0.0)
    return energy.flows(e0, sim.energy())['td_neg'] * 10.0


if __name__ == '__main__':
    ins, outs = [], []
    for d in DEPTHS:
        a, b = stand(d, 'in'), stand(d, 'out')
        ins.append(a['braked'])
        outs.append(b['lifted'])
        print(f'{d:6.0f} m: in {a["braked"] / KWH:6.2f} kWh braked, out {b["lifted"] / KWH:6.2f} kWh lifted '
              f'(wall friction {b["friction"] / KWH:4.2f}, '
              f'damping {b["damping"] / KWH:4.2f} kWh)')
    eta = a['eta']
    e_in, e_out = trip_total(ins), trip_total(outs)
    print(f'\nTrip in to {DEPTHS[-1]:.0f} m ({DEPTHS[-1] / STAND:.0f} stands): {e_in / KWH:.0f} kWh braked, '
          f'{eta * e_in / KWH:.0f} kWh recoverable at {eta:.0%}')
    print(f'Trip out from {DEPTHS[-1]:.0f} m: {e_out / KWH:.0f} kWh lifted; round trip costs {(e_out - e_in) / KWH:.0f} kWh '
          '(mostly wall friction)')
    e_stop = stop_rotation()
    print(f'Stopping 80 rpm rotation: {e_stop / KWH:.3f} kWh braked by the top drive')
