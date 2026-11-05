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

    Mechanical:
        density            mud density at the reference temperature and surface (kg/m^3)
        plastic_viscosity  Bingham PV (Pa s); not used by the mechanics yet
        yield_point        Bingham YP (Pa); not used yet
    Thermal and hydraulic (Wang et al. 2023, SPE-214836-MS, Table 2 defaults):
        viscosity              effective viscosity for hydraulics and heat transfer (Pa s)
        heat_capacity          J/(kg K)
        thermal_conductivity   W/(m K)
        thermal_expansion      volumetric expansion coefficient (1/K), eq. 2
        sound_speed            m/s, sets compressibility in eq. 2
        reference_temperature  C, temperature at which ``density`` applies
    """
    density: float = 1000.0
    plastic_viscosity: float = None
    yield_point: float = None
    viscosity: float = 0.020
    heat_capacity: float = 3750.0
    thermal_conductivity: float = 0.75
    thermal_expansion: float = 1e-4
    sound_speed: float = 1500.0
    reference_temperature: float = 20.0

    def density_at(self, T, p):
        """Mud density at temperature T (C) and gauge pressure p (Pa), Wang et al. (2023) eq. 2."""
        return (1.0 - (T - self.reference_temperature) * self.thermal_expansion) * self.density + p / self.sound_speed ** 2


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
