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
    # False (default): axial and torsional friction are independent 1-D Coulomb
    #   inclusions with the normal force of the static hanging string, so the
    #   only axial-torsional coupling is at the bit.
    # True: 2-D Coulomb inclusion on the combined sliding velocity with the
    #   normal force from the live tension (along-string coupling).
    along_string_coupling: bool = False


@dataclass
class DampingParameters:
    k_t: float = 0.05           # torsional viscous damping rate (1/s)
    k_a: float = 0.05           # axial viscous damping rate (1/s)


@dataclass
class BitParameters:
    """
    Coulomb-type bit-rock interaction (the axial-torsional coupling point).

    stuck   (|omega_b| <= omega_c):  |TOB| <= stick_torque + mu_static * radius * WOB
    sliding:                         TOB  =  mu_kinetic * radius * WOB

    ``radius=None`` uses D/3 (hole diameter / 3), so ``mu_kinetic`` is the
    conventional bit friction coefficient in TOB = mu * WOB * D/3.

    Penetration law: ROP = k_rop * max(WOB - wob_threshold, 0) * |omega_bit|.
    """
    radius: float = None
    mu_static: float = 0.6
    mu_kinetic: float = 0.4
    stick_torque: float = 0.0   # N m, sticking torque at zero WOB
    k_rop: float = 1e-8         # m / (N rad)
    wob_threshold: float = 0.0  # N


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
