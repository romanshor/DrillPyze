"""
Model parameters (SI).  Defaults are order-of-magnitude placeholders; friction
and damping are meant to be calibrated against field data.
"""

from dataclasses import dataclass, field, asdict, replace


@dataclass
class FrictionParameters:
    """Along-string Coulomb friction (static/kinetic, cased/open hole)."""
    mu_static_cased: float = 0.30
    mu_kinetic_cased: float = 0.20
    mu_static_open: float = 0.35
    mu_kinetic_open: float = 0.25
    omega_c: float = 1.0        # rad/s, torsional static->kinetic threshold
    v_c: float = 0.01           # m/s, axial static->kinetic threshold
    # True (default): one 2-D Coulomb inclusion on the combined sliding velocity
    #   (v, omega r) with the normal force from the live tension, so axial and
    #   torsional motion are coupled along the string as well as at the bit.
    # False: independent 1-D axial and torsional inclusions with the normal force
    #   of the static hanging string; the only coupling is then at the bit.
    along_string_coupling: bool = True


@dataclass
class DampingParameters:
    k_t: float = 0.05           # torsional viscous damping rate (1/s)
    k_a: float = 0.05           # axial viscous damping rate (1/s)


@dataclass
class BitParameters:
    """
    Bit-rock interaction: the axial-torsional coupling point at the bit.
    Three reduced models (see docs/derivation/bit-models.md); D = hole
    diameter, a = D/2, omega = bit speed, d = depth of cut per revolution.

    ``coulomb`` -- friction-dominated (worn bit / contact regime); 3 parameters
        TOB = mu_kinetic * WOB * D/3 sliding;  |TOB| <= mu_static * WOB * D/3 sticking
        ROP = omega * WOB / (2 pi a S)

    ``cutting`` -- Detournay & Defourny (1992), quasi-static; 3 parameters
        WOB = W_f + S a d                    (contact threshold + cutting)
        TOB = WOB a / (2 zeta)               (rate independent, no weakening)
        ROP = omega * (WOB - W_f) / (2 pi a S)

    ``regenerative`` -- cutting with the depth of cut set by the hole-bottom
        profile left by the previous blade (Richard, Germay & Detournay 2007);
        2 rock parameters + the blade count
        WOB = S a sum_j d_j,  TOB = WOB a / (2 zeta)

    S = zeta * epsilon is the drilling strength (Pa); epsilon the intrinsic
    specific energy; zeta the ratio of normal to tangential cutting force.
    """
    model: str = 'coulomb'
    mu_static: float = 0.6          # coulomb
    mu_kinetic: float = 0.4         # coulomb
    drilling_strength: float = 50e6 # all models, S = zeta * epsilon (Pa)
    zeta: float = 0.8               # cutting, regenerative
    contact_wob: float = 0.0        # cutting, W_f (N)
    blades: int = 5                 # regenerative
    radius: float = None            # coulomb torque lever; None = D/3


@dataclass
class TopDriveParameters:
    inertia: float = 2000.0     # kg m^2 (including gearbox, reflected to the quill)
    max_torque: float = 0.0     # N m, 0 = unlimited
    setpoint_ramp: float = 60.0 * 2 * 3.141592653589793 / 60.0  # rad/s per s (60 rpm/s)


@dataclass
class ControllerParameters:
    """
    Top-drive speed controller (Aarsnes et al. 2018, Fig. 4).  Use the
    factories in :mod:`DrillPyze.controllers` rather than setting these by hand.
    """
    kp: float = 0.0             # N m s / rad
    ki: float = 0.0             # N m / rad
    impedance: float = 0.0      # Z (rad/s per N m); 0 = no torque feedback
    t_hp: float = 10.0          # s, high-pass time constant on pipe torque
    t_lp: float = 0.005         # s, low-pass time constant on pipe torque
    t_en: float = 0.001         # s, encoder / speed filter time constant
    name: str = 'unset'


@dataclass
class AutodrillerParameters:
    kp: float = 1e-7            # (m/s) per N of WOB error
    ki: float = 1e-8            # (m/s) per N s
    v_max: float = 60.0 / 3600  # m/s


@dataclass
class ModelParameters:
    friction: FrictionParameters = field(default_factory=FrictionParameters)
    damping: DampingParameters = field(default_factory=DampingParameters)
    bit: BitParameters = field(default_factory=BitParameters)
    top_drive: TopDriveParameters = field(default_factory=TopDriveParameters)
    controller: ControllerParameters = None
    autodriller: AutodrillerParameters = field(default_factory=AutodrillerParameters)
    cfl: float = 0.95           # Courant number on the axial wave speed

    def to_dict(self):
        return asdict(self)

    def copy(self, **changes):
        return replace(self, **changes)
