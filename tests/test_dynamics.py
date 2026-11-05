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


@pytest.mark.parametrize('coupled', [True, False])
@pytest.mark.parametrize('direction', [-1, 1])
def test_pickup_and_slackoff_match_torque_and_drag(direction, coupled):
    mu_c, mu_o = 0.2, 0.25
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0),
                     friction(mu_c, mu_c, mu_o, mu_o, along_string_coupling=coupled))
    sim.initialize(2500.0, hole_depth=2600.0, dx=5.0, refresh_distance=1e9)
    sim.step(25.0, v_block=direction * 0.1)
    df = sim.results()
    hk = df.hookload[df.time > 20].mean()
    # coupled: normal force follows the live tension (full soft string);
    # decoupled: normal force of the static hanging string
    td = statics.torque_drag(sim.wellbore, sim.grid, 1200.0, mu_c, mu_o, v=direction,
                             normal_force=None if coupled else sim.fn0)
    assert hk == pytest.approx(td.hookload, rel=0.02)


@pytest.mark.parametrize('v, rpm', [(-0.05, 100.0), (-0.2, 40.0), (0.1, 60.0)])
def test_rotating_while_moving_matches_velocity_split(v, rpm):
    """Cross-coupled friction: steady combined motion matches velocity-split T&D."""
    mu = 0.25
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0),
                     friction(mu, mu, mu, mu, k_t=0.0, k_a=0.0))
    sim.initialize(2500.0, hole_depth=2700.0, dx=5.0, refresh_distance=1e9)
    for k in range(40):
        sim.step(0.25, omega_td=rpm * units.RPM * (k + 1) / 40)
    sim.step(30.0, v_block=v)
    df = sim.results()
    late = df.time > df.time.max() - 8
    td = statics.torque_drag(sim.wellbore, sim.grid, 1200.0, mu, v=v, omega=rpm * units.RPM)
    assert df.hookload[late].mean() == pytest.approx(td.hookload, rel=0.02)
    assert df.pipe_torque[late].mean() == pytest.approx(td.surface_torque, rel=0.02)


def pickup_hookload(rpm, coupled):
    mu = 0.25
    sim = Simulation(well_a(), standard_string(), Fluid(1200.0),
                     friction(mu, mu, mu, mu, k_t=0.0, k_a=0.0, along_string_coupling=coupled))
    sim.initialize(2500.0, hole_depth=2600.0, dx=5.0, refresh_distance=1e9)
    for k in range(40):
        sim.step(0.25, omega_td=rpm * units.RPM * (k + 1) / 40)
    sim.step(30.0, v_block=-0.05)
    df = sim.results()
    late = df.time > df.time.max() - 8
    return df.hookload[late].mean(), df.pipe_torque[late].mean()


def test_cross_coupling_rotation_reduces_drag():
    hk_still, _ = pickup_hookload(0.0, coupled=True)
    hk_rot, _ = pickup_hookload(100.0, coupled=True)
    free = statics.torque_drag(well_a(), standard_string().discretize(2500.0, 5.0), 1200.0, 0.0).hookload
    assert hk_rot - free < 0.2 * (hk_still - free)


def test_decoupled_option_has_no_along_string_coupling():
    """With along_string_coupling=False, rotating does not change pick-up drag."""
    hk_still, _ = pickup_hookload(0.0, coupled=False)
    hk_rot, _ = pickup_hookload(100.0, coupled=False)
    assert hk_rot == pytest.approx(hk_still, rel=0.002)


def test_combined_breakout_uses_vector_friction_limit():
    """
    A stuck string breaks free when the combined axial + tangential load reaches
    mu_s F_N, i.e. earlier than either direction alone would.
    """
    from DrillPyze.kernel import _coulomb_2d
    m_a, m_t, c = 20.0, 5.0, 10.0
    # each component alone is within the static impulse ...
    va, vt = 0.8 * c / m_a, 0.8 * c / m_t
    assert abs(m_a * va) < c and abs(m_t * vt) < c
    # ... but together they exceed it, so the cell slides
    a, t = _coulomb_2d(va, vt, m_a, m_t, c)
    assert a > 0 and t > 0
    # and friction acts exactly against the new sliding direction
    V = np.hypot(a, t)
    assert a == pytest.approx(va - c / m_a * a / V, rel=1e-10)
    assert t == pytest.approx(vt - c / m_t * t / V, rel=1e-10)


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
def drilling_sim(model, rpm=100.0, rop_mh=20.0, kt=0.5, T=60.0, **bit):
    p = frictionless()
    p.damping = DampingParameters(kt, kt)
    p.bit = BitParameters(model=model, **bit)
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), p)
    sim.initialize(2000.0, hole_depth=2000.0, dx=10.0, record_rate=200)
    w = rpm * units.RPM
    for k in range(20):
        sim.step(0.25, omega_td=w * (k + 1) / 20)
    sim.step(T, v_block=rop_mh / 3600.0)
    df = sim.results()
    return sim, df[df.time > df.time.max() - 10]


def test_coulomb_bit_drilling_steady_state():
    S, mu, rpm, rop = 50e6, 0.4, 100.0, 20.0
    sim, late = drilling_sim('coulomb', rpm, rop, mu_static=mu, mu_kinetic=mu, drilling_strength=S)
    D = sim.wellbore.hole_diameter
    a, w, v = D / 2, rpm * units.RPM, rop / 3600
    wob = 2 * np.pi * a * S * v / w                     # ROP = omega WOB / (2 pi a S)
    assert late.v_bit.mean() == pytest.approx(v, rel=0.01)
    assert late.wob.mean() == pytest.approx(wob, rel=0.01)
    assert late.tob.mean() == pytest.approx(mu * wob * D / 3, rel=0.01)
    # no wall friction: surface torque = TOB + viscous damping torque
    g = sim.grid
    viscous = 0.5 * np.sum(g.rho * g.J * g.dx) * w
    assert late.pipe_torque.mean() == pytest.approx(late.tob.mean() + viscous, rel=0.01)


@pytest.mark.parametrize('contact_wob', [0.0, 10e3])
def test_cutting_bit_drilling_steady_state(contact_wob):
    """Detournay & Defourny phase II: WOB = W_f + 2 pi a S v / omega, TOB = a WOB / (2 zeta)."""
    S, zeta, rpm, rop = 50e6, 0.8, 100.0, 20.0
    sim, late = drilling_sim('cutting', rpm, rop, drilling_strength=S, zeta=zeta, contact_wob=contact_wob)
    a, w, v = sim.wellbore.hole_diameter / 2, rpm * units.RPM, rop / 3600
    wob = contact_wob + 2 * np.pi * a * S * v / w
    assert late.v_bit.mean() == pytest.approx(v, rel=0.01)
    assert late.wob.mean() == pytest.approx(wob, rel=0.01)
    assert late.tob.mean() == pytest.approx(a * wob / (2 * zeta), rel=0.01)


def test_regenerative_bit_converges_to_cutting_when_stable():
    """With enough damping the regenerative model drills at the quasi-static cutting state."""
    S, zeta, rpm, rop = 50e6, 0.8, 100.0, 20.0
    sim, late = drilling_sim('regenerative', rpm, rop, kt=5.0, drilling_strength=S, zeta=zeta, blades=5)
    a, w, v = sim.wellbore.hole_diameter / 2, rpm * units.RPM, rop / 3600
    wob = 2 * np.pi * a * S * v / w
    assert late.wob.mean() == pytest.approx(wob, rel=0.01)
    assert late.tob.mean() == pytest.approx(a * wob / (2 * zeta), rel=0.01)
    assert late.v_bit.mean() == pytest.approx(v, rel=0.01)
    assert late.wob.std() < 0.01 * wob


def test_regenerative_bit_self_excites_when_lightly_damped():
    """The depth-of-cut delay destabilises steady drilling: bit bounce and torsional oscillation."""
    sim, late = drilling_sim('regenerative', 100.0, 20.0, kt=0.5, drilling_strength=50e6, zeta=0.8, blades=5)
    rpm = units.rads_to_rpm(late.omega_bit)
    assert (late.wob == 0).mean() > 0.05            # the bit leaves bottom (bit bounce)
    assert rpm.max() - rpm.min() > 100.0           # large torsional oscillation


def test_regenerative_hole_bottom_follows_the_bit():
    sim, late = drilling_sim('regenerative', 100.0, 20.0, kt=5.0, T=30.0, drilling_strength=50e6, zeta=0.8, blades=4)
    s = sim.latest()
    # hole bottom profile: deepest point equals the reported hole depth, and every bin has been cut
    assert sim.surf.max() == pytest.approx(s['hole_depth'], abs=1e-9)
    assert sim.surf.min() > 2000.0 + 0.5 * (s['hole_depth'] - 2000.0)


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


# ------------------------------------------------------------ bit coupling
def loaded_bit_sim(model='coulomb', **bit):
    """Vertical, frictionless string standing on bottom with WOB applied, not rotating."""
    p = frictionless()
    p.damping = DampingParameters(0.5, 0.5)
    p.bit = BitParameters(model=model, **bit)
    sim = Simulation(vertical_well(), uniform_string(), Fluid(1000.0), p)
    sim.initialize(2000.0, hole_depth=2000.0, dx=10.0, record_rate=1000)
    sim.step(4.0, v_block=0.02, omega_td=0.0)       # set down ~80 mm: WOB builds, bit cannot drill
    sim.step(6.0, v_block=0.0)                      # hold the block, let axial waves damp out
    return sim


@pytest.mark.parametrize('model, lever', [('coulomb', 'mu_static'), ('cutting', 'zeta'), ('regenerative', 'zeta')])
def test_bit_sticks_until_breakout_torque(model, lever):
    mus, zeta, r = 0.6, 0.8, 0.08
    sim = loaded_bit_sim(model, mu_static=mus, mu_kinetic=0.4, radius=r, zeta=zeta)
    a = sim.wellbore.hole_diameter / 2
    wob = sim.latest()['wob']
    assert wob > 20e3                                # ~80 mm set-down x EA/L ~ 28 kN
    n0 = len(sim.results())
    sim.step(20.0, omega_td=10 * units.RPM)         # twist the string from surface
    df = sim.results().iloc[n0:]
    moving = np.abs(df.omega_bit.values) > 1e-6
    i = np.argmax(moving)
    assert i > 0                                    # bit was sticking at first
    breakout = np.abs(df.tob.values[:i]).max()
    limit = mus * r * df.wob.values[i - 1] if model == 'coulomb' else a / (2 * zeta) * df.wob.values[i - 1]
    assert breakout == pytest.approx(limit, rel=0.02)


def test_sliding_coulomb_bit_torque_is_linear_in_wob():
    muk, r = 0.4, 0.08
    sim = loaded_bit_sim('coulomb', mu_static=muk, mu_kinetic=muk, radius=r, drilling_strength=500e6)
    n0 = len(sim.results())
    sim.step(30.0, omega_td=60 * units.RPM)
    df = sim.results().iloc[n0:]
    sliding = (df.omega_bit > 1.0) & (df.wob > 1e3)
    np.testing.assert_allclose(df.tob[sliding] / df.wob[sliding], muk * r, rtol=1e-9)


@pytest.mark.parametrize('model', ['coulomb', 'cutting', 'regenerative'])
def test_bit_is_the_coupling_point(model):
    """Rotation changes WOB only through the bit: drilling off at constant block position."""
    sim = loaded_bit_sim(model, mu_static=0.4, mu_kinetic=0.4, zeta=0.8, drilling_strength=50e6)
    wob0 = sim.latest()['wob']
    sim.step(20.0, omega_td=100 * units.RPM)        # block held: rotating drills off the WOB
    s = sim.latest()
    assert s['wob'] < 0.7 * wob0
    assert s['hole_depth'] > 2000.0
