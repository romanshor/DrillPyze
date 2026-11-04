"""
Steady-state Newtonian hydraulics for pipe and annular flow (SI units).

Fanning friction factor: 16/Re (laminar, Re < 2100) and Blasius
0.0791 Re^-0.25 (turbulent).  Annular flow uses the hydraulic diameter
(Do - Di).  Non-Newtonian rheologies are planned.
"""

import numpy as np


def area(Do, Di=None):
    """Flow area (m^2) of a pipe (``Di`` is None) or annulus."""
    return np.pi / 4 * (Do ** 2 - (0.0 if Di is None else Di ** 2))


def fluid_velocity(Q, Do, Di=None):
    """Mean velocity (m/s) for flow rate ``Q`` (m^3/s)."""
    return Q / area(Do, Di)


def hydraulic_diameter(Do, Di=None):
    return Do if Di is None else Do - Di


def reynolds_number(rho, v, mu, Do, Di=None):
    """Reynolds number; ``rho`` kg/m^3, ``v`` m/s, ``mu`` Pa s, diameters m."""
    return rho * abs(v) * hydraulic_diameter(Do, Di) / mu


def friction_factor(Re):
    """Fanning friction factor."""
    Re = max(Re, 1e-12)
    return 16.0 / Re if Re < 2100 else 0.0791 / Re ** 0.25


def frictional_pressure_drop(rho, Q, mu, L, Do, Di=None):
    """
    Frictional pressure drop (Pa) over length ``L`` (m) for flow rate ``Q``
    (m^3/s), density ``rho`` (kg/m^3) and viscosity ``mu`` (Pa s).
    """
    v = fluid_velocity(Q, Do, Di)
    Dh = hydraulic_diameter(Do, Di)
    f = friction_factor(reynolds_number(rho, v, mu, Do, Di))
    return 2.0 * f * rho * v * abs(v) * L / Dh
