"""
Energy accounting (docs/derivation/README.md sec. 11): the work at the string's ends,
the losses to wall friction and damping, the mud's pressure drops times the flow, and
the braking a regenerative drive could recover.  The balances close to the scheme's
numerical dissipation, which is a fraction of a per cent.
"""
import numpy as np
import pytest

from DrillPyze import config, units, energy, EnergyParameters


def _cfg(motor=False, hole=2500.2, thermal=False):
    cfg = config.well_a()
    cfg['parameters']['thermal']['enabled'] = thermal
    cfg['simulation']['hole_depth'] = hole            # bit tags bottom soon after lowering starts
    if motor:
        cfg['parameters']['motor']['enabled'] = True
    return cfg


def _run(cfg, rpm=80.0, t_rot=40.0, t_drill=40.0, wob=60e3):
    """Rotate off bottom, then drill on the autodriller; energy snapshots at the phase ends."""
    sim = config.build(cfg)
    sim.step(t_rot - 20.0, omega_sp=rpm * units.RPM)
    e0 = sim.energy()
    sim.step(20.0)
    e1 = sim.energy()
    sim.tare_hookload()
    sim.step(t_drill - 20.0, wob_sp=wob)
    e2 = sim.energy()
    sim.step(20.0)
    e3 = sim.energy()
    return sim, energy.flows(e0, e1, sim.params.energy), energy.flows(e2, e3, sim.params.energy)


@pytest.fixture(scope='module')
def rotary():
    return _run(_cfg())


@pytest.fixture(scope='module')
def motor():
    return _run(_cfg(motor=True), rpm=40.0, wob=40e3)


@pytest.mark.parametrize('phase', [1, 2])
@pytest.mark.parametrize('case', ['rotary', 'motor'])
def test_balances_close(request, case, phase):
    """Rotation, axial and top drive balances close to the numerical dissipation of the scheme."""
    p = request.getfixturevalue(case)[phase]
    through_t = p['top'] + 1.0
    assert abs(p['numerical_torsion']) <= 0.02 * through_t
    # the top drive's backward-Euler inertia loses only I(Δω)²/2 per step
    assert abs(p['numerical_top_drive']) <= 1e-4 * p['td_pos']
    through_a = abs(p['gravity']) + p['hook_pos'] + p['hook_neg'] + p['bit_axial'] + 1.0
    assert abs(p['numerical_axial']) <= 0.02 * through_a + 50.0


def test_losses_are_dissipative(rotary):
    for p in rotary[1:]:
        for k in ('friction_torsion', 'damping_torsion', 'damping_axial', 'td_pos', 'td_neg', 'hook_pos', 'hook_neg'):
            assert p[k] >= -1e-6, k


def test_rotary_drilling_split(rotary):
    """Off bottom the top drive's work goes to wall friction; on bottom part of it reaches the rock."""
    _, off, on = rotary
    assert off['bit_torsion'] == pytest.approx(0.0, abs=1.0)
    assert off['friction_torsion'] > 0.8 * off['top']
    assert on['bit_torsion'] > 0.1 * on['top']
    assert on['rock'] == pytest.approx(on['bit_torsion'] + on['bit_axial'])
    assert on['motor'] == 0.0


def test_hydraulic_split(rotary):
    """The pump's SPP·Q is shared by the string, nozzles and annulus; it closes by construction."""
    p = rotary[1]
    assert p['h_pump'] > 1e5
    parts = sum(p[f'h_{k}'] for k in ('surface', 'string', 'motor', 'nozzles', 'annulus', 'other'))
    assert parts == pytest.approx(p['h_pump'], rel=1e-9)
    assert p['h_string'] > p['h_annulus'] > 0 and p['h_nozzles'] > 0
    assert abs(p['h_other']) < 0.05 * p['h_pump']
    assert p['pump_drive'] == pytest.approx(p['h_pump'] / 0.9)


def test_motor_drilling(motor):
    """The pumps drive the motor: part of its hydraulic power reaches the bit, the rest is lost."""
    sim, off, on = motor
    assert on['h_motor'] > on['motor'] > 0
    assert on['motor_loss'] > 0
    assert on['rock'] == pytest.approx(on['bit_torsion'] + on['bit_axial'] + on['motor'])
    assert on['h_motor'] > off['h_motor']          # the differential rises with bit torque
    assert off['h_motor'] > 0 and off['motor'] == pytest.approx(0.0, abs=1.0)


def test_no_friction_no_friction_energy():
    cfg = _cfg()
    f = cfg['parameters']['friction']
    for k in ('mu_static_cased', 'mu_kinetic_cased', 'mu_static_open', 'mu_kinetic_open'):
        f[k] = 0.0
    sim = config.build(cfg)
    sim.step(10.0, omega_sp=60 * units.RPM)
    e = sim.energy()
    assert e['friction_torsion'] == pytest.approx(0.0, abs=1e-9 * e['top'])
    assert e['friction_axial'] == pytest.approx(0.0, abs=1e-9 * e['top'])
    assert e['top'] > 0


def _trip(v):
    cfg = _cfg(hole=2505.0)
    cfg['parameters']['hydraulics']['enabled'] = False
    cfg['parameters']['thermal']['flow_rate'] = 0.0
    cfg['simulation'].update(bit_depth=2400.0, block_height=40.0)
    sim = config.build(cfg)
    sim.step(5.0, flow_rate=0.0, omega_sp=0.0)
    e0 = sim.energy()
    sim.step(10.0, v_block=v)
    sim.step(4.0, v_block=0.0)
    return energy.flows(e0, sim.energy())


def test_tripping_in_brakes():
    """Lowering: gravity does the work and the drawworks brakes it (recoverable)."""
    p = _trip(0.5)
    assert p['gravity'] > 0 and p['hook_neg'] > 0.8 * p['gravity']
    assert p['hook_pos'] < 0.01 * p['hook_neg']
    assert p['recoverable'] == pytest.approx(0.9 * p['braking'])


def test_tripping_out_lifts():
    """Hoisting: the drawworks lifts the string against gravity and wall friction."""
    p = _trip(-0.5)
    assert p['gravity'] < 0 and p['hook_pos'] > -p['gravity']
    assert p['hook_pos'] - (-p['gravity']) == pytest.approx(p['friction_axial'] + p['damping_axial'] + p['d_store_axial']
                                                              + p['numerical_axial'], abs=1.0)


def test_stopping_rotation_brakes_the_top_drive():
    sim = config.build(_cfg())
    sim.step(20.0, omega_sp=80 * units.RPM)
    e0 = sim.energy()
    sim.step(5.0, omega_sp=0.0)
    p = energy.flows(e0, sim.energy())
    assert p['td_neg'] > 0
    assert p['d_store_top_drive'] < 0


def test_efficiencies():
    cfg = _cfg()
    cfg['parameters']['energy'] = dict(pump_efficiency=0.8, regen_efficiency=0.5)
    sim = config.build(cfg)
    assert sim.params.energy == EnergyParameters(pump_efficiency=0.8, regen_efficiency=0.5)
    sim.step(5.0, omega_sp=60 * units.RPM)
    e0 = sim.energy()
    sim.step(5.0)
    p = energy.flows(e0, sim.energy(), sim.params.energy)
    assert p['pump_drive'] == pytest.approx(p['h_pump'] / 0.8)
    assert p['pump_loss'] == pytest.approx(p['h_pump'] * 0.25)
    for bad in (0.0, 1.5, -0.1):
        cfg['parameters']['energy']['pump_efficiency'] = bad
        assert any('Pump efficiency' in e for e in config.validate(cfg))
    cfg['parameters']['energy'] = dict(pump_efficiency=0.9, regen_efficiency=2)
    assert any('Regeneration efficiency' in e for e in config.validate(cfg))
    del cfg['parameters']['energy']                 # older cases: the defaults
    assert config.validate(cfg) == [] and config.build(cfg, initialize=False).params.energy == EnergyParameters()


def test_flows_needs_a_later_snapshot():
    sim = config.build(_cfg())
    e = sim.energy()
    with pytest.raises(ValueError):
        energy.flows(e, e)


def test_pause_resume_exact():
    a = config.build(_cfg())
    a.step(0.0, omega_sp=60 * units.RPM)
    a.step(12.0)
    b = config.build(_cfg())
    b.step(0.0, omega_sp=60 * units.RPM)
    for _ in range(4):
        b.step(3.0)
    ea, eb = a.energy(), b.energy()
    for k in ea:
        assert ea[k] == eb[k], k
    assert np.isfinite(list(ea.values())).all()
