"""
Unit conversions.

DrillPyze works in SI units internally (m, s, kg, N, rad, Pa).  Field data are
often in oilfield units; these helpers convert at the input/output boundary.
"""

import numpy as np

# length
FT = 0.3048                 # m per ft
IN = 0.0254                 # m per in

# mass / force
LBM = 0.45359237            # kg per lbm
LBF = 4.4482216152605       # N per lbf
KLBF = 1000.0 * LBF         # N per klbf
TONNE_F = 9806.65           # N per metric tonne-force

# torque
FTLBF = FT * LBF            # N m per ft lbf
KFTLBF = 1000.0 * FTLBF     # N m per kft lbf

# pressure
PSI = 6894.757293168        # Pa per psi

# density
PPG = 119.8264273           # kg/m^3 per lbm/gal
LBM_FT3 = LBM / FT**3       # kg/m^3 per lbm/ft^3
SG = 1000.0                 # kg/m^3 per specific gravity

# flow
GPM = 6.30901964e-5         # m^3/s per US gal/min
GAL = 3.785411784e-3        # m^3 per US gal

# linear weight
LB_FT = LBM / FT            # kg/m per lb/ft

# angular velocity
RPM = 2.0 * np.pi / 60.0    # rad/s per rev/min

# velocity
FT_HR = FT / 3600.0         # m/s per ft/hr
M_HR = 1.0 / 3600.0         # m/s per m/hr

# gravity
G = 9.80665                 # m/s^2

# steel defaults
STEEL_DENSITY = 7850.0      # kg/m^3
STEEL_E = 207e9             # Pa
STEEL_G = 79.3e9            # Pa


def rpm_to_rads(rpm):
    """Convert rev/min to rad/s."""
    return np.asarray(rpm) * RPM


def rads_to_rpm(w):
    """Convert rad/s to rev/min."""
    return np.asarray(w) / RPM
