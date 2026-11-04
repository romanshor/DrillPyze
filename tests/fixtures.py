"""Shared synthetic cases (SI)."""
import numpy as np
from DrillPyze import Survey, Component, Drillstring, Wellbore, Fluid, units


def well_a_survey():
    """Aarsnes & Shor (2018) Well A: vertical to 1500 m, build 3 deg/30 m to 60 deg, hold."""
    md = [0, 1500, 2100, 3000]
    inc = [0, 0, 60, 60]
    azi = [0, 0, 0, 0]
    return Survey(md, inc, azi)


def standard_string(dp_mass=None):
    return Drillstring([
        Component('Collars', 230.0, od=5.75 * units.IN, id=2.8125 * units.IN),
        Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN, linear_mass=dp_mass),
    ])


def well_a(casing_shoe=1000.0):
    return Wellbore(well_a_survey(), casing_shoe=casing_shoe, hole_diameter=8.5 * units.IN)


def vertical_well(td=3000.0):
    return Wellbore(Survey([0, td], [0, 0], [0, 0]), casing_shoe=0.0, hole_diameter=8.5 * units.IN)


def uniform_string():
    return Drillstring([Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN)])
