"""
Rock strength at the bit (P1 + T1 + T2) and the wellbore temperature model
(Wang et al. 2023, SPE-214836-MS).
"""
import numpy as np
import pytest

from DrillPyze import (Survey, Component, Drillstring, Wellbore, Fluid, Simulation, ModelParameters,
                       FrictionParameters, DampingParameters, BitParameters, units)
from DrillPyze.formation import Formation, FormationLayer
from DrillPyze.thermal import WellboreThermal, ThermalParameters

FT = units.FT


# ------------------------------------------------------------------ strength
def test_drained_ccs_is_mohr_coulomb():
    f = Formation(layers=[FormationLayer(0, ucs=100e6, friction_angle=30.0, porosity=0.25)],
                  weakening_threshold=1e9)
    tvd, pbh = 2000.0, 25e6
    r = f.strength(2000.0, tvd, pbh)
    n_phi = (1 + np.sin(np.radians(30))) / (1 - np.sin(np.radians(30)))      # = 3
    assert r['w_drained'] == 1.0
    assert r['ccs'] == pytest.approx(100e6 + (pbh - f.pore_pressure(tvd)) * n_phi)
    # identical to the form in Caicedo et al.: UCS + DP + 2 DP sin(phi) / (1 - sin(phi))
    dp = pbh - f.pore_pressure(tvd)
    assert r['ccs'] == pytest.approx(100e6 + dp + 2 * dp * 0.5 / 0.5)


def test_low_porosity_uses_undrained_skempton_pore_pressure():
    lay = FormationLayer(0, ucs=100e6, friction_angle=36.0, porosity=0.005)
    f = Formation(layers=[lay], weakening_threshold=1e9)
    tvd, pbh = 2000.0, 22e6
    r = f.strength(2000.0, tvd, pbh)
    pp_lp = f.pore_pressure(tvd) - f.skempton_b * (f.overburden(tvd) - pbh) / 3
    assert r['w_drained'] == 0.0
    assert r['dp_undrained'] == pytest.approx(pbh - pp_lp)
    assert r['ccs'] > f.strength(2000.0, tvd, pbh)['ucs_t'] + (pbh - f.pore_pressure(tvd)) * r['n_phi']
    # porosity between the cutoffs blends linearly
    lay.porosity = 0.125
    assert f.strength(2000.0, tvd, pbh)['w_drained'] == pytest.approx(0.5)


def test_thermal_weakening_matches_fit():
    f = Formation(layers=[FormationLayer(0, ucs=125e6)])
    # TVD that puts the rock at 450 C: Xu et al. (2019) air-cooled loss 38 %, fit gives 36 %
    tvd = (450.0 - f.surface_temperature) / f.geothermal_gradient
    r = f.strength(tvd, tvd, 0.0)
    assert r['ucs_t'] / 125e6 == pytest.approx(1 - 1.41e-3 * (450 - 193), rel=1e-6)
    tvd = (150.0 - f.surface_temperature) / f.geothermal_gradient
    assert f.strength(tvd, tvd, 0.0)['ucs_t'] == 125e6            # below the threshold


def test_cooling_thermal_stress_and_skin_ratio():
    f = Formation()
    tvd = 2500.0
    base = f.strength(tvd, tvd, 25e6)
    cool = f.strength(tvd, tvd, 25e6, t_mud=base['t_rock'] - 50.0, rop=0.0, omega=0.0)
    sig = f.thermal_expansion * f.youngs_modulus * 50.0 / (1 - f.poisson_ratio)
    assert cool['sigma_th'] == pytest.approx(sig)
    assert cool['ccs_eff'] == pytest.approx(base['ccs'] - f.thermal_stress_efficiency * sig)
    # fast drilling cuts deeper than the cooled skin: smaller effect
    fast = f.strength(tvd, tvd, 25e6, t_mud=base['t_rock'] - 50.0, rop=100 / 3600, omega=60 * units.RPM, blades=5)
    t_pass = 2 * np.pi / (5 * 60 * units.RPM)
    expected = np.sqrt(np.pi * f.thermal_diffusivity * t_pass) / (100 / 3600 * t_pass)
    assert fast['skin_ratio'] == pytest.approx(min(1, expected))
    assert cool['ccs_eff'] < fast['ccs_eff'] < base['ccs']
    # strength floor
    f.thermal_stress_efficiency = 1e3
    assert f.strength(tvd, tvd, 25e6, t_mud=0.0)['ccs_eff'] == pytest.approx(f.strength_floor * base['ucs_t'])


# ------------------------------------------------------------------- thermal
def forge_case():
    """SPE-214836-MS validation case: FORGE 16(78)-32, ~5900 ft MD, 550 gpm (Tables 1-2)."""
    sv = Survey([0, 5000 * FT, 5900 * FT, 6200 * FT], [0, 0.8, 18, 18], [0, 0, 0, 0])
    wb = Wellbore(sv, casing_shoe=4000 * FT, hole_diameter=0.22, casing_id=0.22)
    ds = Drillstring([Component('Bit+BHA', 0.8, od=0.20, id=0.07), Component('DP', 100.0, od=0.127, id=0.11)])
    fl = Fluid(density=1078, viscosity=0.020, heat_capacity=3750, thermal_conductivity=0.75, thermal_expansion=1e-4)
    return wb, ds, fl, Formation()


def run_forge(rop_fth=0.0, minutes=30, flow=550 * units.GPM, t_in=35.0):
    wb, ds, fl, fm = forge_case()
    th = WellboreThermal(wb, ds, fl, fm, ThermalParameters(flow_rate=flow, inlet_temperature=t_in),
                         bit_depth=5900 * FT, hole_depth=5900 * FT)
    hole, rop = 5900 * FT, rop_fth * FT / 3600
    out = []
    for _ in range(int(minutes * 6)):
        hole += rop * 10
        s = th.advance(10.0, bit_depth=hole, hole_depth=hole, rop=rop)
        out.append(s['t_bit'])
    return th, np.array(out)


def test_static_well_stays_at_geothermal():
    th, dht = run_forge(flow=0.0, minutes=20)
    np.testing.assert_allclose(th.T_ann, th.T_far, atol=1e-9)
    np.testing.assert_allclose(th.T_f, np.broadcast_to(th.T_far[:, None], th.T_f.shape), atol=1e-9)


def test_energy_balance_while_circulating():
    wb, ds, fl, fm = forge_case()
    th = WellboreThermal(wb, ds, fl, fm, ThermalParameters(flow_rate=550 * units.GPM, inlet_temperature=35.0),
                         bit_depth=5900 * FT, hole_depth=5900 * FT)
    th.advance(600.0)
    # outer boundary heat flow is negligible over this time (5 m of formation); the energy
    # removed by the mud equals the loss of stored energy in fluid + formation
    E0 = th.stored_energy()
    removed = 0.0
    F = fl.density * fl.heat_capacity * 550 * units.GPM
    for _ in range(60):
        th.advance(5.0)
        removed += F * (th.T_ann[0] - th.inlet_temperature) * 5.0
    loss = E0 - th.stored_energy()
    assert removed == pytest.approx(loss, rel=0.03)


def test_circulation_cools_bottom_like_spe214836():
    """Wang et al. (2023) Fig. 6: DHT falls from ~310 F to ~210 F in 30 min at 550 gpm."""
    th, dht = run_forge(minutes=30)
    to_f = lambda c: c * 9 / 5 + 32
    assert to_f(th.state()['t_rock']) == pytest.approx(327, abs=2)
    assert 300 < to_f(dht[0]) < 330
    assert 195 < to_f(dht[-1]) < 225
    assert np.all(np.diff(dht) < 1e-9)                      # monotone cooling


def test_drilling_raises_dht_linearly_with_rop():
    """Wang et al. (2023) Fig. 7: +1.3, +2.2, +2.7 F at 40, 60, 72 ft/h after 30 min."""
    _, base = run_forge(0.0)
    rise = {r: (run_forge(r)[1][-1] - base[-1]) * 9 / 5 for r in (40, 60, 72)}
    paper = {40: 1.3, 60: 2.2, 72: 2.7}
    for r in rise:
        assert rise[r] == pytest.approx(paper[r], rel=0.3)
    # approximately linear in ROP
    assert rise[72] / rise[40] == pytest.approx(72 / 40, rel=0.1)


def test_rop_to_zero_converges_to_static_wellbore():
    """The paper's convergence check: as ROP -> 0 the growing-wellbore result tends to the static one."""
    _, base = run_forge(0.0, minutes=10)
    diffs = [abs(run_forge(r, minutes=10)[1][-1] - base[-1]) for r in (20, 5, 1)]
    assert diffs[0] > diffs[1] > diffs[2]
    assert diffs[2] < 0.02


def test_grid_grows_with_the_hole():
    th, _ = run_forge(rop_fth=60, minutes=30)                # 30 ft drilled; cells cut at 1.5 dx
    edges = th.x_top[-1] + th.dx[-1]
    assert edges == pytest.approx(5930 * FT, abs=0.05)
    assert th.dx.max() <= 1.5 * th.p.dx + 1e-9


def test_ecd_exceeds_static_density_when_circulating():
    th, _ = run_forge(minutes=5)
    assert th.state()['ecd'] > 1078 * 0.98                  # hot mud is lighter, friction adds back
    th0, _ = run_forge(minutes=5, flow=0.0)
    assert th.state()['p_bh'] > th0.state()['p_bh'] - 0.5e6


# --------------------------------------------------------------- integration
def test_formation_strength_sets_bit_drilling_strength():
    """Cutting bit with S = zeta * CCS_eff / eta: steady drilling WOB follows the rock strength."""
    sv = Survey([0, 3000], [0, 0], [0, 0])
    wb = Wellbore(sv, casing_shoe=0.0, hole_diameter=8.5 * units.IN)
    ds = Drillstring([Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN)])
    fm = Formation(layers=[FormationLayer(0, ucs=80e6, friction_angle=36.0, porosity=0.005)],
                   thermal_stress_efficiency=0.0)
    p = ModelParameters(friction=FrictionParameters(0, 0, 0, 0), damping=DampingParameters(0.5, 0.5, fluid=False, loss_factor=0.0),
                        bit=BitParameters(model='cutting', zeta=0.8, strength_source='formation', cutting_efficiency=1.0))
    p.thermal.enabled = False
    sim = Simulation(wb, ds, Fluid(1000.0), p, formation=fm)
    sim.initialize(2000.0, hole_depth=2000.0, dx=10.0)
    S = sim.rock['S']
    assert S == pytest.approx(0.8 * sim.rock['ccs_eff'])
    assert sim.rock['ccs_eff'] > 80e6                        # overbalance confines the rock
    w = 100 * units.RPM
    for k in range(20):
        sim.step(0.25, omega_td=w * (k + 1) / 20)
    sim.step(60.0, v_block=20 / 3600)
    df = sim.results()
    late = df[df.time > df.time.max() - 10]
    a = wb.hole_diameter / 2
    # quasi-static cutting law with the formation-derived S: WOB = 2 pi a S v / omega
    assert late.wob.mean() == pytest.approx(2 * np.pi * a * sim.rock['S'] * late.v_bit.mean() / late.omega_bit.mean(), rel=0.01)
    assert late.v_bit.mean() == pytest.approx(20 / 3600, rel=0.1)


def test_simulation_with_thermal_model_records_slow_channels():
    wb, ds, fl, fm = forge_case()
    p = ModelParameters(bit=BitParameters(model='cutting', strength_source='formation'))
    sim = Simulation(wb, ds, fl, p, formation=fm)
    sim.initialize(5900 * FT, hole_depth=5900 * FT, dx=20.0)
    sim.step(20.0, omega_sp=60 * units.RPM, flow_rate=550 * units.GPM, inlet_temperature=35.0)
    slow = sim.slow_results()
    assert len(slow) == 10                                   # every 2 s of simulated time
    assert slow.t_bit.iloc[-1] < slow.t_rock.iloc[-1]        # circulation cools the bit
    assert (slow.sigma_th > 0).all() and (slow.ccs_eff < slow.ccs).all()


def test_pause_resume_exact_with_slow_models():
    def make():
        wb, ds, fl, fm = forge_case()
        p = ModelParameters(bit=BitParameters(model='cutting', strength_source='formation'))
        sim = Simulation(wb, ds, fl, p, formation=fm)
        sim.initialize(5900 * FT, hole_depth=5900 * FT, dx=20.0)
        return sim
    a, b = make(), make()
    a.step(9.0, omega_sp=60 * units.RPM, v_block=0.005, flow_rate=0.03)
    for _ in range(36):
        b.step(0.25, omega_sp=60 * units.RPM, v_block=0.005, flow_rate=0.03)
    np.testing.assert_array_equal(a.w, b.w)
    np.testing.assert_array_equal(a.thermal.T_ann, b.thermal.T_ann)
    assert a.slow_results().equals(b.slow_results())
