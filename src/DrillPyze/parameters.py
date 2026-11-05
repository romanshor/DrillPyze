"""
Model parameters (SI).  Defaults are order-of-magnitude placeholders; friction
and damping are meant to be calibrated against field data.
"""

from dataclasses import dataclass, field, asdict, replace

from .thermal import ThermalParameters


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
    # fraction of the static-kinetic difference that applies axially:
    # mu_s,axial = mu_k + axial_static_contrast * (mu_s - mu_k).
    # Default 0: the mu_s / mu_k contrast was identified from torsional start-ups
    # (Aarsnes & Shor 2018); with a step law along a multi-km string any axial
    # contrast produces slack-off / pick-up stick-slip (docs/derivation, sec. 5.5).
    axial_static_contrast: float = 0.0
    # 'step': mu_s below the threshold, mu_k above (Aarsnes & Shor 2018)
    # 'exponential': mu_k + (mu_s - mu_k) exp(-s) with s the normalised sliding speed
    weakening: str = 'step'


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
    blades: int = 5                 # regenerative (also sets the time between blade passes for T2)
    radius: float = None            # coulomb torque lever; None = D/3
    # 'constant': S = drilling_strength.  'formation': S = zeta * CCS_eff / cutting_efficiency,
    # with CCS_eff from the formation at the bit (P1 + T1 + T2, see DrillPyze.formation).
    strength_source: str = 'constant'
    cutting_efficiency: float = 1.0 # eta: epsilon = CCS_eff / eta (1 for a sharp cutter)


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
class DrawworksParameters:
    block_acceleration: float = 0.0   # m/s^2, 0 = block speed changes instantly


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
    drawworks: DrawworksParameters = field(default_factory=DrawworksParameters)
    thermal: ThermalParameters = field(default_factory=ThermalParameters)
    thermal_update_period: float = 2.0   # s between thermal / rock-strength updates
    cfl: float = 0.95           # Courant number on the axial wave speed

    def to_dict(self):
        return asdict(self)

    def copy(self, **changes):
        return replace(self, **changes)
