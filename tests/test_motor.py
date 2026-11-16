"""
Mud motor (docs/derivation/README.md sec. 6.4): speed from flow, differential pressure
from torque, stall, reactive torque and toolface, the motor loss in the hydraulics, and
exact pause/resume.
"""
import numpy as np
import pytest

from DrillPyze import config, units, MotorParameters


def _cfg(**motor):
    cfg = config.well_a()
    cfg['parameters']['thermal']['enabled'] = False
    cfg['simulation']['hole_depth'] = 2500.1          # bit tags bottom early
    cfg['parameters']['motor'].update(enabled=True, **motor)
    return cfg


def _drill(sim, wob=40e3, t_rot=20.0, t_drill=40.0, rpm=40.0):
    sim.step(t_rot, omega_sp=rpm * units.RPM)
    sim.tare_hookload()
    return sim.step(t_drill, wob_sp=wob)


@pytest.fixture(scope='module')
def drilling():
    sim = config.build(_cfg())
    L = _drill(sim)
    return sim, L


def test_well_a_defaults_to_no_motor():
    cfg = config.well_a()
    assert cfg['parameters']['motor']['enabled'] is False
    sim = config.build(cfg, initialize=False)
    assert not sim.params.motor.enabled
    # the generic motor is ready to switch on: 0.28 rev/gal, K_T and stall torque consistent
    m = MotorParameters()
    assert m.speed_constant * units.GAL == pytest.approx(0.28)
    assert m.torque_constant() == pytest.approx(m.efficiency / (2 * np.pi * m.speed_constant))
    assert m.stall_torque() == pytest.approx(m.torque_constant() * (m.dp_stall - m.dp_no_load))


def test_no_motor_reports_no_motor_channels():
    sim = config.build(config.well_a())
    L = sim.step(1.0, omega_sp=30 * units.RPM)
    assert 'dp_motor' not in L and 'toolface' not in L
    r = sim.results()
    assert np.all(r['dp_motor'] == 0) and np.all(r['omega_motor'] == 0)


def test_motor_speed_from_flow():
    """Off bottom: rotor speed = 2 pi k_n Q, differential = off-bottom differential."""
    sim = config.build(_cfg())
    L = sim.step(20.0, omega_sp=0.0)
    m = sim.params.motor
    q = sim.flow_rate_now()
    assert L['omega_motor'] == pytest.approx(2 * np.pi * m.speed_constant * q, rel=1e-12)
    assert L['dp_motor'] == pytest.approx(m.dp_no_load)
    assert not L['motor_stalled']
    # no surface rotation and the string at rest: the bit turns at the rotor speed
    assert L['omega_bit'] == pytest.approx(L['omega_motor'], abs=0.05)


def test_slip_and_bypass_reduce_speed():
    q_by = 0.002
    sim = config.build(_cfg(slip=0.2, bypass_flow=q_by))
    L = sim.step(20.0, omega_sp=0.0)
    m = sim.params.motor
    expect = 2 * np.pi * m.speed_constant * (sim.flow_rate_now() - q_by) * (1 - 0.2 * m.dp_no_load / m.dp_max)
    assert L['omega_motor'] == pytest.approx(expect, rel=1e-12)


def test_differential_from_bit_torque(drilling):
    """On bottom: dp = dp_0 + TOB / K_T in every recorded row."""
    sim, L = drilling
    m = sim.params.motor
    r = sim.results()
    on = r['tob'] > 100.0
    assert on.sum() > 100
    assert np.allclose(r['dp_motor'][on], m.dp_no_load + r['tob'][on] / m.torque_constant(), rtol=1e-12)
    assert L['dp_motor'] > m.dp_no_load


def test_motor_differential_adds_to_spp():
    """
    The motor differential is a loss in series with the nozzles: once the circulation
    settles (off bottom, steady pump rate) SPP is higher by exactly that differential,
    and the bottomhole pressure is unchanged.
    """
    cfg = _cfg()
    cfg['simulation']['hole_depth'] = 2510.0
    sims = []
    for on in (True, False):
        c = config.merge(cfg, {'parameters': {'motor': {'enabled': on}}})
        s = config.build(c)
        s.step(90.0, omega_sp=0.0)
        sims.append(s.hydraulics_results().iloc[-1])
    h, hr = sims
    # to within the compressibility of the mud: the higher bore pressure packs the mud a
    # little denser, so it flows slower and loses ~1 % less to friction above the motor
    assert h['spp'] - hr['spp'] == pytest.approx(h['dp_motor'], rel=0.02)
    assert h['p_bh'] == pytest.approx(hr['p_bh'], rel=1e-6)


def test_stall():
    """A low stall differential: TOB is capped at the stall torque, the bit stops, dp = dp_stall."""
    sim = config.build(_cfg(dp_max=230 * units.PSI, dp_stall=250 * units.PSI))
    m = sim.params.motor
    _drill(sim, wob=60e3, t_drill=30.0)
    r = sim.results()
    st = r['dp_motor'] >= m.dp_stall * (1 - 1e-12)
    assert st.sum() > 10
    assert np.all(r['tob'] <= m.stall_torque() * (1 + 1e-12))
    assert np.allclose(r['omega_bit'][st], 0.0)
    assert np.allclose(r['tob'][st], m.stall_torque())
    # picking up clears the stall
    L = sim.step(20.0, v_block=-60 / 3600)
    assert not L['motor_stalled']


def test_slide_reactive_torque_turns_toolface_left():
    """Sliding (no surface rotation): more WOB, more reactive torque, the toolface turns left."""
    sim = config.build(_cfg())
    sim.step(20.0, omega_sp=0.0)          # never rotated: the string holds its toolface
    sim.tare_hookload()
    sim.step(60.0, wob_sp=20e3)
    a1 = sim.sc[sim_k().S_PHIS]
    tob1 = sim.latest()['tob']
    sim.step(60.0, wob_sp=60e3)
    a2 = sim.sc[sim_k().S_PHIS]
    tob2 = sim.latest()['tob']
    assert tob2 > tob1
    assert a2 < a1                         # counter-clockwise looking downhole
    tf = sim.latest()['toolface']
    assert 0.0 <= tf < 2 * np.pi
    assert tf == pytest.approx((sim.params.motor.toolface_offset + a2) % (2 * np.pi))


def sim_k():
    from DrillPyze import kernel
    return kernel


def test_pause_resume_exact():
    a = config.build(_cfg())
    b = config.build(_cfg())
    _drill(a, t_rot=5.0, t_drill=10.0)
    b.step(5.0, omega_sp=40 * units.RPM)
    b.tare_hookload()
    for _ in range(20):
        b.step(0.5, wob_sp=40e3)
    ra, rb = a.results(), b.results()
    assert ra.shape == rb.shape
    assert np.array_equal(ra.to_numpy(), rb.to_numpy())


@pytest.mark.parametrize('change, message', [
    (dict(speed_constant=0.0), 'speed constant'),
    (dict(efficiency=1.5), 'efficiency'),
    (dict(dp_stall=50 * units.PSI), 'stall differential'),
    (dict(slip=0.9), 'slip'),
])
def test_validation(change, message):
    errs = config.validate(_cfg(**change))
    assert any(message in e for e in errs), errs
