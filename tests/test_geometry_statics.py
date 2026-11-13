import numpy as np
import pytest

from DrillPyze import Survey, Component, Drillstring, Wellbore, units, statics
from DrillPyze.geometry import buoyed_weight
from fixtures import well_a_survey, uniform_string, vertical_well, standard_string, well_a


def test_vertical_survey_tvd():
    s = Survey([0, 1000, 2000], [0, 0, 0], [0, 0, 0])
    np.testing.assert_allclose(s.tvd_at([500, 1500]), [500, 1500], atol=1e-6)


def test_constant_build_arc():
    s = well_a_survey()
    R = 600.0 / np.radians(60.0)
    assert s.tvd[2] == pytest.approx(1500 + R * np.sin(np.radians(60)), rel=1e-9)
    inc, _ = s.inc_azi([1800.0])
    assert np.degrees(inc[0]) == pytest.approx(30.0, abs=1e-6)


def test_azimuth_wraps_through_north():
    s = Survey([0, 500, 600, 700], [0, 0, 30, 30], [350, 350, 350, 10])
    _, azi = s.inc_azi([650.0])
    d = np.degrees(azi[0])
    assert d > 355 or d < 5          # interpolates across north, not through south


def test_hanging_vertical_hookload_is_buoyed_weight():
    w = vertical_well()
    g = uniform_string().discretize(3000.0, 10.0)
    r = statics.torque_drag(w, g, 1200.0, 0.3)
    expected = np.sum(buoyed_weight(g, 1200.0) * g.dx)
    assert r.hookload == pytest.approx(expected, rel=1e-9)
    # rotation, pick-up and slack-off are all friction free in a vertical hole
    for kw in (dict(omega=1.0), dict(v=-1.0), dict(v=1.0)):
        assert statics.torque_drag(w, g, 1200.0, 0.3, **kw).hookload == pytest.approx(expected, rel=1e-9)


def test_inclined_straight_hole_torque_and_drag():
    inc = 60.0
    w = Wellbore(Survey([0, 3000], [inc, inc], [45, 45]), 0.0, 8.5 * units.IN)
    g = uniform_string().discretize(2000.0, 10.0)
    wb = buoyed_weight(g, 1000.0)
    L = g.dx.sum()
    mu = 0.3
    s, c = np.sin(np.radians(inc)), np.cos(np.radians(inc))
    rot = statics.torque_drag(w, g, 1000.0, mu, omega=1.0)
    assert rot.surface_torque == pytest.approx(mu * wb[0] * s * g.r_contact[0] * L, rel=1e-9)
    assert rot.hookload == pytest.approx(wb[0] * c * L, rel=1e-9)
    pu = statics.torque_drag(w, g, 1000.0, mu, v=-1.0)
    so = statics.torque_drag(w, g, 1000.0, mu, v=1.0)
    assert pu.hookload == pytest.approx(wb[0] * (c + mu * s) * L, rel=1e-9)
    assert so.hookload == pytest.approx(wb[0] * (c - mu * s) * L, rel=1e-9)


def test_rotation_reduces_drag():
    w = well_a()
    g = standard_string().discretize(2500.0, 10.0)
    free = statics.torque_drag(w, g, 1200.0, 0.3).hookload
    pu = statics.torque_drag(w, g, 1200.0, 0.3, v=-1.0).hookload
    pu_rot = statics.torque_drag(w, g, 1200.0, 0.3, v=-0.05, omega=10.0).hookload
    assert free < pu_rot < pu
    assert (pu_rot - free) < 0.2 * (pu - free)


def test_discretization_respects_component_boundaries():
    g = standard_string().discretize(2500.0, 7.0)
    assert g.dx.sum() == pytest.approx(2500.0)
    assert np.all(g.dx <= 7.0 + 1e-9)
    # last 230 m are collars
    collars = np.isclose(g.od, 5.75 * units.IN)
    assert g.dx[collars].sum() == pytest.approx(230.0)
    assert np.all(np.diff(collars.astype(int)) >= 0)
