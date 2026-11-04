"""
Verification of the coupled dynamic model against analytic results and the
steady-state torque & drag model.
"""
import numpy as np
import pytest

from DrillPyze import (Simulation, ModelParameters, FrictionParameters, DampingParameters,
                       BitParameters, Fluid, units, statics, controllers)
from fixtures import well_a, vertical_well, uniform_string, standard_string


def frictionless():
    return ModelParameters(
        friction=FrictionParameters(0.0, 0.0, 0.0, 0.0),
        damping=DampingParameters(0.0, 0.0))


def friction(mus_c, muk_c, mus_o, muk_o, k_t=0.05, k_a=0.05, **kw):
    return ModelParameters(friction=FrictionParameters(mus_c, muk_c, mus_o, muk_o, **kw),
                           damping=DampingParameters(k_t, k_a))


def crossing_time(t, y, level):
    i = np.argmax(y >= level)
    return np.interp(level, y[i - 1:i + 1], t[i - 1:i + 1])


# --------------------------------------------------------------- wave physics
def test_torsional_wave_arrival_time():
    L = 3000.0
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), frictionless())
    sim.initialize(L, hole_depth=L + 10, dx=5.0, record_rate=2000)
    w0 = 1.0
    sim.step(2.5, omega_td=w0)
    df = sim.results()
    c = sim.grid.c_torsional[0]
    # free end doubles the incident wave; half-height crossing marks arrival
    t_arr = crossing_time(df.time.values, df.omega_bit.values, w0)
    assert t_arr == pytest.approx(L / c, rel=0.01)


def test_axial_wave_arrival_time():
    L = 3000.0
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), frictionless())
    sim.initialize(L, hole_depth=L + 10, dx=5.0, record_rate=2000)
    v0 = 0.01
    sim.step(1.0, v_block=v0)
    df = sim.results()
    c = sim.grid.c_axial[0]
    t_arr = crossing_time(df.time.values, df.v_bit.values, v0)
    assert t_arr == pytest.approx(L / c, rel=0.01)


def test_torsional_fundamental_frequency():
    """Fixed top (stiff top drive), free bottom: f1 = c / (4L)."""
    L = 2000.0
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), frictionless())
    sim.initialize(L, hole_depth=L + 10, dx=5.0, record_rate=500)
    sim.step(0.3, omega_td=1.0)          # kick
    sim.step(20.0, omega_td=0.0)         # hold top fixed
    df = sim.results()
    df = df[df.time > 1.0]
    y, t = df.omega_bit.values, df.time.values
    # the kick travels as a pulse; successive bit peaks are one period apart
    peaks = t[1:-1][(y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > 0.5 * y.max())]
    period = np.diff(peaks).mean()
    c = sim.grid.c_torsional[0]
    assert 1 / period == pytest.approx(c / (4 * L), rel=0.002)
    # first-order upwind diffusion is negligible at this resolution
    assert y.max() == pytest.approx(2.0, rel=0.01)


# ----------------------------------------------------------------- statics
def test_string_at_rest_stays_at_rest():
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0), friction(0.3, 0.2, 0.35, 0.25))
    sim.initialize(2500.0, hole_depth=2510.0)
    hk0 = sim.latest()['hookload']
    sim.step(5.0)
    assert np.abs(sim.w).max() == 0.0
    assert np.abs(sim.v).max() == 0.0
    assert sim.latest()['hookload'] == pytest.approx(hk0, rel=1e-12)


def test_steady_rotation_matches_torque_and_drag():
    mu_c, mu_o = 0.2, 0.25
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0),
                     friction(mu_c, mu_c, mu_o, mu_o, k_t=0.0, k_a=0.0))
    sim.initialize(2500.0, hole_depth=2510.0, dx=5.0)
    sim.step(40.0, omega_td=0.0)
    for k in range(60):                     # ramp to 80 rpm
        sim.step(0.25, omega_td=80 * units.RPM * (k + 1) / 60)
    sim.step(40.0)
    df = sim.results()
    tq = df.pipe_torque[df.time > df.time.max() - 10].mean()
    td = statics.torque_drag(sim.wellbore, sim.grid, 1200.0, mu_c, mu_o, omega=1.0)
    assert tq == pytest.approx(td.surface_torque, rel=0.02)
    # rotating off bottom with no axial motion: hookload is the free-hanging weight
    hk = df.hookload[df.time > df.time.max() - 10].mean()
    free = statics.torque_drag(sim.wellbore, sim.grid, 1200.0, 0.0).hookload
    assert hk == pytest.approx(free, rel=0.01)


@pytest.mark.parametrize('direction', [-1, 1])
def test_pickup_and_slackoff_match_torque_and_drag(direction):
    mu_c, mu_o = 0.2, 0.25
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0),
                     friction(mu_c, mu_c, mu_o, mu_o))
    sim.initialize(2500.0, hole_depth=2600.0, dx=5.0, refresh_distance=1e9)
    sim.step(25.0, v_block=direction * 0.1)
    df = sim.results()
    hk = df.hookload[df.time > 20].mean()
    td = statics.torque_drag(sim.wellbore, sim.grid, 1200.0, mu_c, mu_o, v=direction)
    assert hk == pytest.approx(td.hookload, rel=0.02)


def test_rotation_while_picking_up_matches_velocity_split():
    mu = 0.25
    v, rpm = -0.05, 100.0
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0), friction(mu, mu, mu, mu, k_t=0.0, k_a=0.0))
    sim.initialize(2500.0, hole_depth=2600.0, dx=5.0, refresh_distance=1e9)
    for k in range(40):
        sim.step(0.25, omega_td=rpm * units.RPM * (k + 1) / 40)
    sim.step(30.0, v_block=v)
    df = sim.results()
    late = df.time > df.time.max() - 8
    td = statics.torque_drag(sim.wellbore, sim.grid, 1200.0, mu, v=v, omega=rpm * units.RPM)
    assert df.hookload[late].mean() == pytest.approx(td.hookload, rel=0.02)
    assert df.pipe_torque[late].mean() == pytest.approx(td.surface_torque, rel=0.03)


# ------------------------------------------------------------- stick-slip
def bit_rpm_swing(sim, last=20.0):
    df = sim.results()
    late = df[df.time > df.time.max() - last]
    rpm = units.rads_to_rpm(late.omega_bit.values)
    return rpm.max() - rpm.min(), rpm.mean()


def test_off_bottom_stick_slip_needs_static_kinetic_difference():
    """Aarsnes & Shor (2018): off-bottom stick-slip arises from mu_s > mu_k."""
    swings = {}
    for label, mus in (('equal', 0.25), ('drop', 0.45)):
        sim = Simulation(well_a(), standard_string(), Fluid(1200.0),
                         friction(mus, 0.25, mus, 0.25, k_t=0.2))
        sim.initialize(2500.0, hole_depth=2510.0, dx=10.0)
        sim.step(90.0, omega_sp=60 * units.RPM)
        swings[label] = bit_rpm_swing(sim)
    assert swings['equal'][0] < 10.0                 # start-up transient has decayed
    assert swings['drop'][0] > 100.0                 # sustained stick-slip limit cycle


# ------------------------------------------------------------------ drilling
def test_on_bottom_drilling_steady_state():
    p = frictionless()
    p.damping = DampingParameters(0.5, 0.5)
    p.bit = BitParameters(radius=0.1, mu_static=0.4, mu_kinetic=0.4, omega_b=0.0, k_rop=1e-7)
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), p)
    sim.initialize(2000.0, hole_depth=2000.0, dx=10.0)
    w = 100 * units.RPM
    for k in range(20):
        sim.step(0.25, omega_td=w * (k + 1) / 20)
    rop = 20.0 / 3600.0
    sim.step(60.0, v_block=rop)
    s = sim.latest()
    expected_wob = rop / (1e-7 * w)
    assert s['v_bit'] == pytest.approx(rop, rel=0.02)
    assert s['wob'] == pytest.approx(expected_wob, rel=0.02)
    assert s['tob'] == pytest.approx(0.4 * expected_wob * 0.1, rel=0.02)
    # no wall friction: surface torque = TOB + viscous damping torque k_t rho J omega L
    g = sim.grid
    viscous = 0.5 * np.sum(g.rho * g.J * g.dx) * w
    assert s['pipe_torque'] == pytest.approx(s['tob'] + viscous, rel=0.01)
    assert s['hole_depth'] > 2000.0 + 0.9 * 60 * rop


def test_bit_lifts_off_when_picked_up():
    p = frictionless()
    p.damping = DampingParameters(0.5, 0.5)
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), p)
    sim.initialize(2000.0, hole_depth=2000.0)
    sim.step(5.0, v_block=0.01)
    assert sim.latest()['wob'] > 0
    sim.step(10.0, v_block=-0.05)
    s = sim.latest()
    assert s['wob'] == 0.0
    assert s['bit_depth'] < s['hole_depth']


# ------------------------------------------------------------- API behaviour
def test_pause_resume_is_exact():
    def make():
        sim = Simulation(well_a(), standard_string(), Fluid(1200.0), friction(0.4, 0.25, 0.4, 0.25))
        sim.initialize(2500.0, hole_depth=2510.0)
        return sim
    a, b = make(), make()
    a.step(10.0, omega_sp=60 * units.RPM)
    for _ in range(40):
        b.step(0.25, omega_sp=60 * units.RPM)
    np.testing.assert_array_equal(a.w, b.w)
    np.testing.assert_array_equal(a.tau, b.tau)
    np.testing.assert_allclose(a.results().values, b.results().values)


def test_external_controller_matches_internal_pi():
    """A Python PI run through motor_torque reproduces the internal controller."""
    sim_i = Simulation(well_a(), standard_string(), Fluid(1200.0), friction(0.25, 0.25, 0.25, 0.25))
    sim_i.initialize(2500.0, hole_depth=2510.0)
    c = sim_i.controller
    sim_e = Simulation(well_a(), standard_string(), Fluid(1200.0), friction(0.25, 0.25, 0.25, 0.25))
    sim_e.initialize(2500.0, hole_depth=2510.0)
    sp = 30 * units.RPM
    sim_i.params.top_drive.setpoint_ramp = 0.0
    sim_i.update_parameters()
    integ = 0.0
    for _ in range(int(5.0 / sim_e.dt)):
        e = sp - sim_e.latest()['omega_td']
        u = c.kp * e + c.ki * (integ + e * sim_e.dt)
        integ += e * sim_e.dt
        sim_e.step(sim_e.dt, motor_torque=u)
    sim_i.step(5.0, omega_sp=sp)
    assert sim_e.latest()['omega_td'] == pytest.approx(sim_i.latest()['omega_td'], rel=0.02)


def test_soft_torque_reflection_minimum_at_fc():
    zeta, I, fc = 300.0, 2000.0, 0.25
    c = controllers.soft_torque(zeta, I, fc)
    f = np.linspace(0.05, 1.0, 2000)
    r = controllers.reflection_coefficient(c, zeta, I, f)
    assert f[np.argmin(r)] == pytest.approx(fc, rel=0.02)


def test_ztorque_suppresses_off_bottom_stick_slip():
    """Aarsnes, di Meglio & Shor (2018): ZTorque avoids the limit cycle a stiff controller sustains."""
    swing = {}
    for name in ('stiff', 'z_torque'):
        sim = Simulation(well_a(), standard_string(), Fluid(1200.0), friction(0.4, 0.22, 0.45, 0.25, k_t=0.1))
        sim.initialize(2500.0, hole_depth=2505.0, dx=10.0)
        I = sim.params.top_drive.inertia
        sim.params.controller = (controllers.stiff(sim.zeta_p, I) if name == 'stiff'
                                 else controllers.z_torque(sim.zeta_p, I, t_hp=2.0))
        sim.update_parameters()
        sim.step(60.0, omega_sp=60 * units.RPM)
        swing[name] = bit_rpm_swing(sim)[0]
    assert swing['stiff'] > 100.0
    assert swing['z_torque'] < 5.0
