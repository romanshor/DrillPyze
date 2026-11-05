"""Verification of the distributed damping models (F1 fluid, M1 material, F2 yield)."""
import numpy as np
import pytest
from scipy.signal import hilbert, butter, sosfiltfilt

from DrillPyze import (Simulation, ModelParameters, FrictionParameters, DampingParameters, Fluid)
from DrillPyze.kernel import _kv_solve
from fixtures import vertical_well, uniform_string


def _sim(damping, fluid=None, hole=3500.0):
    p = ModelParameters(friction=FrictionParameters(0, 0, 0, 0), damping=damping)
    sim = Simulation(vertical_well(4000), uniform_string(), fluid or Fluid(1200.0, viscosity=0.02), p)
    sim.initialize(3000.0, hole_depth=hole, dx=10.0, record_rate=200)
    return sim


def _decay_ratio(sim, mode):
    """Free decay of the first mode after a short kick; damping ratio from the envelope."""
    if mode == 'axial':
        sim.step(0.5, v_block=0.2, omega_td=0.0)
        sim.step(60.0, v_block=0.0)
        y, f = sim.results()['v_bit'].values, sim.modal_frequencies()[0]
    else:
        sim.step(0.5, omega_td=2.0, v_block=0.0)
        sim.step(60.0, omega_td=0.0)
        y, f = sim.results()['omega_bit'].values, sim.modal_frequencies()[1]
    t = sim.results()['time'].values
    m = t > 1.0
    fs = 1.0 / np.median(np.diff(t))
    yb = sosfiltfilt(butter(2, [0.6 * f, 1.4 * f], btype='band', fs=fs, output='sos'), y[m])
    env, tt = np.abs(hilbert(yb)), t[m]
    span = tt[-1] - tt[0]
    k = (tt > tt[0] + 0.15 * span) & (tt < tt[-1] - 0.15 * span)
    return -np.polyfit(tt[k], np.log(env[k]), 1)[0] / (2 * np.pi * f)


def test_kv_solve_conserves_momentum_and_dissipates():
    rng = np.random.default_rng(0)
    n = 50
    x = rng.normal(size=n)
    m = rng.uniform(1, 3, n)
    h = rng.uniform(5, 10, n)
    k = np.concatenate([[0.0], rng.uniform(10, 100, n - 1), [0.0]])
    x0 = x.copy()
    _kv_solve(x, m, h, k, 0.01, np.empty(n), np.empty(n))
    assert np.sum(m * h * x) == pytest.approx(np.sum(m * h * x0), rel=1e-12)
    assert np.sum(m * h * x ** 2) < np.sum(m * h * x0 ** 2)


def test_first_mode_frequencies_quarter_wave():
    sim = _sim(DampingParameters())
    fa, ft = sim.modal_frequencies()
    g = sim.grid
    assert fa == pytest.approx(g.c_axial[0] / (4 * 3000.0), rel=1e-9)
    assert ft == pytest.approx(g.c_torsional[0] / (4 * 3000.0), rel=1e-9)


@pytest.mark.parametrize('label, damping', [
    ('lumped', DampingParameters(0.1, 0.1, fluid=False, loss_factor=0.0)),
    ('fluid', DampingParameters(fluid=True, loss_factor=0.0)),
    ('material', DampingParameters(fluid=False, loss_factor=0.02)),
    ('fluid+material', DampingParameters(fluid=True, loss_factor=0.02)),
])
@pytest.mark.parametrize('mode', ['axial', 'torsional'])
def test_free_decay_matches_predicted_damping_ratio(label, damping, mode):
    sim = _sim(damping)
    row = sim.damping_summary().set_index('mode').loc[mode]
    if label == 'lumped':
        w = 2 * np.pi * row['f']
        predicted = 0.1 / (2 * w)       # k0 only (no fluid, no material)
    else:
        predicted = row['zeta_fluid'] + row['zeta_material']
    predicted += row['zeta_numerical']
    measured = _decay_ratio(sim, mode)
    assert measured == pytest.approx(predicted, rel=0.08, abs=3e-4), (label, mode)


def test_steady_rotation_sees_only_couette_drag():
    """The Stokes-layer and Rayleigh-alpha terms act on oscillations only."""
    sim = _sim(DampingParameters(fluid=True, loss_factor=0.01))
    omega = 2 * np.pi * 120 / 60
    sim.step(60.0, omega_td=omega, v_block=0.0)
    r = sim.results()
    tau0 = r['pipe_torque'][r['time'] > 60.0 - 4 / sim.modal_frequencies()[1]].mean()   # 4 periods
    g = sim.grid
    expected = np.sum(sim.k0t * sim.mt * omega * g.dx)        # Couette only
    with_stokes = expected + np.sum(sim.k1t * sim.mt * omega * g.dx)
    assert tau0 == pytest.approx(expected, rel=0.02)
    assert with_stokes > 3 * expected                          # would have been far larger


def test_yield_drag_applies_with_pumps_off_only():
    mud = Fluid(1200.0, plastic_viscosity=0.02, yield_point=10.0)
    v = 0.05
    hk = {}
    for label, q in (('off', 0.0), ('on', 0.03)):
        sim = _sim(DampingParameters(fluid=True, loss_factor=0.01), fluid=mud)
        sim.step(1.0, omega_td=0.0, v_block=0.0, flow_rate=q)
        h_static = sim.latest()['hookload']
        sim.step(60.0, v_block=v)
        r = sim.results()
        hk[label] = h_static - r['hookload'][r['time'] > 61.0 - 4 / sim.modal_frequencies()[0]].mean()
    g = sim.grid
    fy_total = np.sum(10.0 * np.pi * g.od * g.dx)
    couette = np.sum(sim.k0a * sim.ma * v * g.dx)
    assert hk['off'] == pytest.approx(fy_total + couette, rel=0.05)
    assert hk['on'] == pytest.approx(couette, rel=0.05, abs=50.0)


def test_pause_resume_exact_with_damping():
    a = _sim(DampingParameters(fluid=True, loss_factor=0.02))
    b = _sim(DampingParameters(fluid=True, loss_factor=0.02))
    a.step(4.0, omega_sp=6.0, v_block=0.02)
    for _ in range(8):
        b.step(0.5, omega_sp=6.0, v_block=0.02)
    assert np.array_equal(a.w, b.w) and np.array_equal(a.v, b.v)
    assert np.array_equal(a.wbar, b.wbar)
