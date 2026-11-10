"""
Verification of the hydraulics model:

* steady state against Shor, Introduction to Drilling Engineering, ch. 7
  (Example 7.1, Practice Problems 7.1 and 7.2);
* the transient (single-phase drift flux, Wang et al. 2023) against its own
  exact discrete steady state, hydrostatics, mass conservation, the
  Joukowsky pressure surge and the acoustic travel time;
* surge and swab from pipe motion;
* set-point ramps.
"""
import copy
import numpy as np
import pytest

from DrillPyze import (Simulation, ModelParameters, FrictionParameters, DampingParameters, Fluid,
                       Survey, Wellbore, Drillstring, Component, RampParameters, HydraulicsParameters, units)
from DrillPyze.hydraulics import (steady_state, nozzle_area, section, WellboreHydraulics, G)

PPG = 119.826427
IN = 0.0254
FT = 0.3048
GPM = 6.30901964e-5
CP = 1e-3
YP = 0.4788026          # Pa per lbf/100 ft^2
PSI = 6894.7573


# ------------------------------------------------------------ textbook cases
def _textbook(mw, pv, yp, q, dp_od, dp_id, dp_len, dc, hole, nozzles, tvd):
    pipe = [('DP', dp_id * IN, dp_len * FT)]
    ann = []
    if dc:
        pipe.append(('DC', dc[1] * IN, dc[2] * FT))
        ann.append(('DC/OH', hole * IN, dc[0] * IN, dc[2] * FT))
    ann.append(('DP/OH', hole * IN, dp_od * IN, dp_len * FT))
    return steady_state(q * GPM, mw * PPG, pv * CP, yp * YP, pipe, ann, nozzle_area(nozzles), tvd * FT)


EX71 = (11.5, 25, 14, 400, 5.0, 4.276, 9800, (6.75, 2.875, 500), 8.75, [13] * 3, 10300)
PP71 = (9.5, 20, 0, 320, 4.5, 3.826, 7500, (6.25, 2.5, 600), 8.5, [11] * 3, 8100)
PP72 = (12.2, 30, 18, 450, 5.5, 4.778, 11500, None, 9.875, [14] * 3, 11500)

# Values as printed in the chapter.  Sections marked * do not follow from the
# chapter's own turbulent formula (its intermediate powers are mis-evaluated);
# for those the reference is the formula re-evaluated: 191, 248 and 532 psi.
CASES = [
    ('Example 7.1', EX71,
     dict(v=[8.94, 19.77, 5.27, 3.17], mua=[69.6, 38.5, 51.6, 107.8], re=[5861, 15776, 1777, 957],
          regime=['turbulent', 'turbulent', 'laminar', 'laminar'], dp=[571, 191, 34, 238]),   # DC* (printed 172)
     dict(bit=1119, ecd=12.01)),
    ('Problem 7.1', PP71,
     dict(v=[8.93, 20.92, 3.94, 2.51], mua=[20, 20, 20, 20], re=[15052, 23038, 3187, 3611],
          regime=['turbulent'] * 4, dp=[408, 248, 18.8, 55.9]),                                # DC* (printed 228)
     dict(bit=1157, ecd=9.68)),
    ('Problem 7.2', PP72,
     dict(v=[8.05, 2.73], mua=[101.2, 174.2], re=[4306, 633],
          regime=['turbulent', 'laminar'], dp=[532, 286]),                                     # DP* (printed 599)
     dict(bit=1120, ecd=12.68)),
]


@pytest.mark.parametrize('name, data, sec, tot', CASES, ids=[c[0] for c in CASES])
def test_textbook_steady_state(name, data, sec, tot):
    r = _textbook(*data)
    df = r['sections']
    # pipe sections first, then annulus from the bit up, as in the chapter tables
    order = [0, 1, 2, 3] if len(df) == 4 else [0, 1]
    df = df.iloc[order]
    np.testing.assert_allclose(df.velocity / FT, sec['v'], rtol=0.01)
    np.testing.assert_allclose(df.apparent_viscosity / CP, sec['mua'], rtol=0.01)
    np.testing.assert_allclose(df.reynolds, sec['re'], rtol=0.01)
    assert list(df.regime) == sec['regime']
    np.testing.assert_allclose(df.dp / PSI, sec['dp'], rtol=0.05)   # printed values are rounded
    assert r['bit'] / PSI == pytest.approx(tot['bit'], rel=0.005)
    assert r['ecd'] / PPG == pytest.approx(tot['ecd'], abs=0.01)
    assert r['spp'] == pytest.approx(r['string'] + r['bit'] + r['annulus'] + r['surface'])
    assert r['bhp'] == pytest.approx(r['hydrostatic'] + r['annulus'])


def test_newtonian_is_bingham_with_zero_yield():
    a = section(0.03, 1200.0, 0.02, 0.0, 1000.0, 0.2159, 0.127)
    b = section(0.03, 1200.0, 0.02, 1e-12, 1000.0, 0.2159, 0.127)
    assert a['dp'] == pytest.approx(b['dp'], rel=1e-9)


# --------------------------------------------- model geometry = textbook well
def _example71_model(sound_speed=1500.0, expansion=0.0, dx=30.48, newtonian=False):
    """Example 7.1 as a model well: vertical, 10,300 ft, 8.75 in hole, DC 500 ft, DP 9,800 ft."""
    L = 10300 * FT
    sv = Survey([0.0, L + 200.0], [0.0, 0.0], [0.0, 0.0])
    wb = Wellbore(sv, casing_shoe=0.0, hole_diameter=8.75 * IN)
    ds = Drillstring([Component('DC', 500 * FT, 6.75 * IN, 2.875 * IN),
                      Component('DP', 9800 * FT, 5.0 * IN, 4.276 * IN)])
    if newtonian:
        fluid = Fluid(11.5 * PPG, viscosity=1 * CP, sound_speed=sound_speed, thermal_expansion=expansion)
    else:
        fluid = Fluid(11.5 * PPG, plastic_viscosity=25 * CP, yield_point=14 * YP,
                      sound_speed=sound_speed, thermal_expansion=expansion)
    grid = ds.discretize(L, dx)
    hp = HydraulicsParameters(dx=dx)
    return WellboreHydraulics(wb, grid, fluid, hp, nozzle_area([13] * 3), 0.95, L, L, L)


def test_model_steady_state_matches_textbook():
    """The transient model's exact discrete steady state reproduces Example 7.1 (near-incompressible mud)."""
    h = _example71_model(sound_speed=1e6)
    h.set_steady(400 * GPM)
    st = h.state()
    ref = _textbook(*EX71)
    assert st['spp'] == pytest.approx(ref['spp'], rel=0.005)
    assert st['p_bh'] == pytest.approx(ref['bhp'], rel=0.002)
    assert st['ecd'] == pytest.approx(ref['ecd'], rel=0.002)
    assert st['dp_bit'] == pytest.approx(ref['bit'], rel=0.02)


def test_compressible_mud_raises_bhp():
    """With rho(p) of Wang et al. eq. 2 the column is heavier than the textbook constant density."""
    a, b = _example71_model(sound_speed=1e6), _example71_model(sound_speed=1500.0)
    a.set_static()
    b.set_static()
    pa, pb = a.state()['p_bh'], b.state()['p_bh']
    rho, tvd = 11.5 * PPG, 10300 * FT
    # first-order: p = rho g z (1 + rho g z / (2 rho c^2))
    expected = rho * G * tvd * (1 + rho * G * tvd / (2 * rho * 1500.0 ** 2))
    assert pa == pytest.approx(rho * G * tvd, rel=1e-6)
    assert pb == pytest.approx(expected, rel=2e-4)


def test_static_column_stays_static():
    h = _example71_model()
    h.set_static()
    p0 = h.profiles()['p'].values.copy()
    for _ in range(20):
        h.advance(1.0, 0.0, 0.0, h.bit_depth, h.hole_depth)
    np.testing.assert_allclose(h.profiles()['p'].values, p0, rtol=1e-9, atol=1e-3)
    assert abs(h.state()['q_out']) < 1e-9


def test_transient_converges_to_steady_state_and_conserves_mass():
    h = _example71_model()
    h.set_static()
    m0 = h.mass()
    Q = 400 * GPM
    h.advance(20.0, 0.0, Q, h.bit_depth, h.hole_depth)       # 20 s pump ramp
    for _ in range(280):
        h.advance(1.0, Q, Q, h.bit_depth, h.hole_depth)
    ref = copy.deepcopy(h)
    ref.set_steady(Q)
    st, rs = h.state(), ref.state()
    assert st['spp'] == pytest.approx(rs['spp'], rel=1e-3)
    assert st['p_bh'] == pytest.approx(rs['p_bh'], rel=1e-4)
    assert st['q_out'] == pytest.approx(Q, rel=1e-3)
    # mass balance to round-off: change in the well = pumped in - returned
    assert h.mass() - m0 == pytest.approx(h.acc[0] - h.acc[1], abs=1e-6 * m0)


def test_joukowsky_surge_and_acoustic_travel_time():
    """A fast pump start raises the standpipe pressure by rho c dQ / A and reaches the bit after L / c."""
    h = _example71_model(newtonian=True)     # 1 cP: friction small against the surge
    h.set_static()
    p_static = h.state()
    bore0 = h.profiles().query("path == 'drillstring'")['p'].values.copy()
    Q = 100 * GPM
    A = np.pi / 4 * (4.276 * IN) ** 2
    rho, c = 11.5 * PPG, 1500.0
    h.advance(0.02, 0.0, Q, h.bit_depth, h.hole_depth)
    t, rise_bit, spp = [0.02], [], []
    for _ in range(250):
        h.advance(0.01, Q, Q, h.bit_depth, h.hole_depth)
        t.append(t[-1] + 0.01)
        bore = h.profiles().query("path == 'drillstring'")['p'].values
        rise_bit.append(bore[-2] - bore0[-2])
        spp.append(h.state()['spp'] - p_static['spp'])
    t = np.array(t[1:])
    # Joukowsky at the standpipe before the reflection returns (t < 2 L / c)
    dp_j = rho * c * Q / A
    early = (t > 0.2) & (t < 0.6)
    assert np.mean(np.array(spp)[early]) == pytest.approx(dp_j, rel=0.08)
    # arrival at the bit: first time the bore pressure there rises by half of dp_j
    L_bit = h.e_d[-2] + 0.5 * h.L_d[-2]
    t_arr = t[np.argmax(np.array(rise_bit) > 0.5 * dp_j)]
    assert t_arr == pytest.approx(L_bit / c, rel=0.05)


def test_surge_and_swab():
    """Running in raises the bottomhole pressure, pulling out lowers it, by similar amounts."""
    out = {}
    for label, v in (('surge', 0.5), ('swab', -0.5)):
        h = _example71_model()
        L = h.bit_depth
        # the bit starts 200 m off bottom so it can move both ways
        h.bit_depth = L - 200.0
        h.e_d = h._edges(h.bit_depth)
        h.set_static()
        p0 = h.state()['p_bh']
        bit = h.bit_depth
        for _ in range(60):
            bit += 0.5 * v
            h.advance(0.5, 0.0, 0.0, bit, h.hole_depth)
        out[label] = h.state()['p_bh'] - p0
        out[label + '_qout'] = h.state()['q_out']
    assert out['surge'] > 0 > out['swab']
    assert abs(out['surge']) == pytest.approx(abs(out['swab']), rel=0.25)
    # running in at 0.5 m/s displaces at least the steel volume at the flowline
    A_steel = np.pi / 4 * ((5.0 * IN) ** 2 - (4.276 * IN) ** 2)
    assert out['surge_qout'] > 0.9 * 0.5 * A_steel


# ------------------------------------------------------- simulation coupling
def _sim(ramps=None, fluid=None, **hp):
    from fixtures import well_a, standard_string
    p = ModelParameters(friction=FrictionParameters(0.25, 0.25, 0.25, 0.25),
                        ramps=ramps or RampParameters(), hydraulics=HydraulicsParameters(**hp))
    s = Simulation(well_a(), standard_string(), fluid or Fluid(1200.0), p)
    return s.initialize(2500.0, hole_depth=2510.0)


def test_pump_start_in_simulation_reaches_steady_pressures():
    s = _sim()
    Q = 0.03
    s.step(150.0, omega_sp=0.0, v_block=0.0, flow_rate=Q)
    ref = copy.deepcopy(s.hydraulics)
    ref.set_steady(Q)
    st = s.hydraulics.state()
    assert st['spp'] == pytest.approx(ref.state()['spp'], rel=2e-3)
    assert st['p_bh'] == pytest.approx(ref.state()['p_bh'], rel=2e-4)
    r = s.hydraulics_results()
    assert len(r) >= 149 and r['spp'].iloc[-1] > 5e6
    # the rock-strength model sees the circulating bottomhole pressure
    assert s.rock is None or s.rock['p_bh'] == pytest.approx(st['p_bh'], rel=1e-3)


def test_setpoint_ramps_are_linear_in_time():
    s = _sim(ramps=RampParameters(omega_sp=5.0, v_block=2.0, flow_rate=20.0))
    w = 60 * units.RPM
    s.step(2.5, omega_sp=w, v_block=0.01, flow_rate=0.02)
    assert s.latest()['omega_sp'] == pytest.approx(0.5 * w, rel=1e-3)
    assert s.latest()['flow_rate'] == pytest.approx(0.02 * 2.5 / 20.0, rel=1e-3)
    assert s.latest()['v_block'] == pytest.approx(0.01, rel=1e-3)     # 2 s ramp already complete
    s.step(2.5)
    assert s.latest()['omega_sp'] == pytest.approx(w, rel=1e-6)
    assert s.latest()['v_block'] == pytest.approx(0.01, rel=1e-6)
    # a new target ramps from the current value
    s.step(1.0, omega_sp=0.0)
    assert s.latest()['omega_sp'] == pytest.approx(0.8 * w, rel=1e-3)


def test_mode_switch_ramps_from_present_speed():
    s = _sim(ramps=RampParameters(omega_sp=4.0, omega_td=4.0))
    s.step(6.0, omega_td=60 * units.RPM)
    w_now = s.latest()['omega_td']
    s.step(1.0, omega_sp=120 * units.RPM)
    # one second into a 4 s ramp from the present speed
    assert s.latest()['omega_sp'] == pytest.approx(w_now + 0.25 * (120 * units.RPM - w_now), rel=1e-3)


def test_pause_resume_exact_with_hydraulics():
    a, b = _sim(), _sim()
    a.step(6.0, omega_sp=6.0, v_block=0.02, flow_rate=0.03)
    for _ in range(12):
        b.step(0.5, omega_sp=6.0, v_block=0.02, flow_rate=0.03)
    assert np.array_equal(a.hydraulics.m_a, b.hydraulics.m_a)
    assert a.hydraulics.state()['spp'] == b.hydraulics.state()['spp']
