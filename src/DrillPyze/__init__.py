"""
DrillPyze: a coupled axial-torsional drillstring dynamics simulator,
torque & drag model and controller test bench.  SI units throughout.
"""

from . import units, survey, drillstring, wellbore, geometry, statics, kernel
from . import parameters, controllers, simulation, hydraulics, formation, thermal, config
from .survey import Survey
from .drillstring import Component, Drillstring
from .wellbore import Wellbore, Fluid
from .parameters import (ModelParameters, FrictionParameters, DampingParameters,
                         BitParameters, MotorParameters, TopDriveParameters, ControllerParameters,
                         AutodrillerParameters, RampParameters, HydraulicsParameters)
from .simulation import Simulation
from .statics import torque_drag
from .formation import Formation, FormationLayer
from .thermal import WellboreThermal, ThermalParameters

__version__ = "1.0.0"             # PEP 440 (keep pyproject.toml in step)
__release_date__ = "2026-11-16"
__release_name__ = "1.0"
__author__ = 'Roman J. Shor'
__credits__ = 'Texas A&M University; University of Calgary'
