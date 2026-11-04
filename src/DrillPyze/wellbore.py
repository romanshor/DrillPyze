"""
Wellbore (survey + casing/open-hole geometry) and drilling fluid.  SI units.
"""

from dataclasses import dataclass
import numpy as np

from .survey import Survey


@dataclass
class Fluid:
    """
    Drilling fluid.

    Parameters
    ----------
    density : float
        Mud density (kg/m^3).
    plastic_viscosity : float, optional
        Bingham plastic viscosity (Pa s); not used by the mechanics model yet.
    yield_point : float, optional
        Bingham yield point (Pa); not used by the mechanics model yet.
    """
    density: float = 1000.0
    plastic_viscosity: float = None
    yield_point: float = None


class Wellbore:
    """
    A wellbore: survey, casing shoe depth and hole sizes.

    Parameters
    ----------
    survey : Survey
    casing_shoe : float
        Measured depth of the deepest casing shoe (m).  Friction above the shoe
        uses the cased-hole friction factors.
    hole_diameter : float
        Open-hole (bit) diameter (m).
    casing_id : float, optional
        Inner diameter of the casing (m); informational.
    """

    def __init__(self, survey: Survey, casing_shoe: float, hole_diameter: float, casing_id: float = None):
        self.survey = survey
        self.casing_shoe = float(casing_shoe)
        self.hole_diameter = float(hole_diameter)
        self.casing_id = casing_id

    def is_cased(self, md):
        """Boolean mask: True where ``md`` is inside casing."""
        return np.asarray(md) <= self.casing_shoe
