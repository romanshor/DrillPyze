"""
Model parameters (SI).  Defaults are order-of-magnitude placeholders; friction
and damping are meant to be calibrated against field data.
"""

from dataclasses import dataclass, field, asdict, replace

from .thermal import ThermalParameters
from .units import GAL, PSI


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
    """
    Distributed damping (docs/derivation/README.md, sec. 5.6).

    Fluid (F1, from the mud and the string/hole geometry, no free parameter):
        steady annular Couette drag on the full velocity, plus Stokes-layer
        drag sqrt(rho mu omega_ref / 2) on the inner and outer walls acting on
        the oscillatory part of the velocity (high-passed below
        ``high_pass_fraction * omega_ref``).  omega_ref defaults to the first
        (quarter-wave) mode of the string, axial and torsional separately.
    Fluid (F2): Bingham yield-stress drag tau_y * pi * OD per length, added to the
        Coulomb wall friction; applied per ``yield_drag``.
    Material (M1): Rayleigh damping zeta(f) = alpha/(2 omega) + beta omega/2 fitted so
        zeta = loss_factor/2 at f_ref and at max_frequency (and below it in
        between); beta is Kelvin-Voigt strain-rate damping sigma = E (eps + beta deps/dt),
        alpha acts on the oscillatory velocity.
    Calibration: k_t, k_a are additional mass-proportional rates (1/s) on the
        full velocity (the original lumped model), default 0.

    Validity: the fluid terms are evaluated at a single reference frequency, so
    modes above it are under-damped by sqrt(f / f_ref); together with the
    numerical dissipation of the first-order scheme this limits the model to
    vibrations below ``max_frequency`` (10 Hz).
    """
    k_t: float = 0.0            # additional torsional viscous damping rate (1/s)
    k_a: float = 0.0            # additional axial viscous damping rate (1/s)
    fluid: bool = True          # F1: Couette + Stokes-layer mud damping
    loss_factor: float = 0.01   # M1: material loss factor eta (zeta = eta/2 at f_ref and max_frequency); 0 = off
    reference_frequency_axial: float = None       # Hz; None = first axial mode
    reference_frequency_torsional: float = None   # Hz; None = first torsional mode
    high_pass_fraction: float = 0.25   # Stokes term acts above this fraction of f_ref
    yield_drag: str = 'pumps_off'      # F2: 'always' | 'pumps_off' | 'off'
    max_frequency: float = 10.0        # Hz, validity limit; upper fit point of the material damping


@dataclass
class BitParameters:
    """
    Bit-rock interaction: the axial-torsional coupling point at the bit.
    Three reduced models (see docs/derivation/bit-models.md); D = hole
    diameter, a = D/2, omega = bit speed, d = depth of cut per revolution.

    ``coulomb`` -- friction-dominated (worn bit / contact regime); 3 parameters::

        TOB = mu_kinetic * WOB * D/3 sliding;  |TOB| <= mu_static * WOB * D/3 sticking
        ROP = omega * WOB / (2 pi a S)

    ``cutting`` -- Detournay & Defourny (1992), quasi-static; 3 parameters::

        WOB = W_f + S a d                    (contact threshold + cutting)
        TOB = WOB a / (2 zeta)               (rate independent, no weakening)
        ROP = omega * (WOB - W_f) / (2 pi a S)

    ``regenerative`` -- cutting with the depth of cut set by the hole-bottom
    profile left by the previous blade (Richard, Germay & Detournay 2007);
    2 rock parameters + the blade count::

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
    # hydraulics: total flow area of the nozzles (m^2; default 5 x 12/32 in) and discharge coefficient
    nozzle_area: float = 3.5628e-4
    discharge_coefficient: float = 0.95


@dataclass
class MotorParameters:
    """
    Positive-displacement mud motor between the string and the bit
    (docs/derivation/README.md, sec. 6.4; Nguyen et al. 2017).  Defaults: a generic
    6-3/4" 7/8-lobe motor.  ``enabled=False`` means no motor in the BHA.

    Speed:  N_r = k_n (Q - Q_bypass) (1 - s dp/dp_max)        (rev/s; k_n in rev/m^3)
    Torque: T = K_T (dp - dp_0),  K_T = eta / (2 pi k_n)       (N m)
    Stall:  T <= K_T (dp_stall - dp_0); beyond it the rotor stops and dp = dp_stall.
    """
    enabled: bool = False
    speed_constant: float = 0.28 / GAL   # k_n, rev per m^3 of flow (0.28 rev/gal)
    efficiency: float = 0.75             # eta, mechanical efficiency
    dp_no_load: float = 100 * PSI        # dp_0, off-bottom differential pressure (Pa)
    dp_max: float = 750 * PSI            # rated (maximum recommended) differential (Pa)
    dp_stall: float = 1125 * PSI         # stall differential (Pa)
    slip: float = 0.0                    # s, fractional speed drop at dp_max (leakage)
    bypass_flow: float = 0.0             # flow not passing through the power section (m^3/s)
    toolface_offset: float = 0.0         # toolface (rad) when the string-bottom angle is zero

    def torque_constant(self):
        """K_T = eta / (2 pi k_n), N m per Pa of differential pressure."""
        import math
        return self.efficiency / (2.0 * math.pi * self.speed_constant)

    def stall_torque(self):
        """Largest torque the motor delivers, K_T (dp_stall - dp_0) (N m)."""
        return self.torque_constant() * (self.dp_stall - self.dp_no_load)


@dataclass
class TopDriveParameters:
    """Top-drive inertia and optional torque limit."""
    inertia: float = 2000.0     # kg m^2 (including gearbox, reflected to the quill)
    max_torque: float = 0.0     # N m, 0 = unlimited


@dataclass
class RampParameters:
    """
    Ramp times (s) for set-point changes on inputs that have no controller.
    Such a set point approaches a new value linearly from its current value
    over its ramp time (0 = instant); ramps are evaluated inside the time
    stepping, so pausing and resuming stays exact.

    Set points that feed a controller -- the top-drive speed set point
    (``omega_sp``, PI / SoftTorque / ZTorque) and the autodriller WOB set
    point (``wob_sp``) -- are not ramped by default: the controller shapes the
    response.  Their entries exist only for studies of set-point profiling.
    """
    omega_sp: float = 0.0           # controlled (top-drive speed controller): no ramp
    omega_td: float = 1.0           # prescribed top-drive speed (1 s = linear between 1 s samples)
    motor_torque: float = 0.0       # torque command from an external controller: no ramp
    v_block: float = 2.0            # block velocity (open loop)
    wob_sp: float = 0.0             # controlled (autodriller): no ramp
    flow_rate: float = 3.0          # pump rate (open loop)
    inlet_temperature: float = 3.0  # mud inlet temperature (open loop)

    @classmethod
    def instant(cls):
        """All ramps off (inputs change in one time step)."""
        return cls(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


@dataclass
class HydraulicsParameters:
    """
    Circulating-system hydraulics (DrillPyze.hydraulics; docs/derivation/hydraulics.md).
    Steady losses per Shor, Introduction to Drilling Engineering, ch. 7; transients from
    the single-phase drift-flux model of Wang et al. (2023, SPE-214836-MS).
    """
    enabled: bool = True
    dx: float = 30.0                   # m, cell length (both flow paths)
    cut_factor: float = 1.5            # bottom cells split above this multiple of dx
    cfl: float = 0.8                   # acoustic Courant number
    iterations: int = 2                # Picard iterations on the implicit face friction
    back_pressure: float = 0.0         # Pa (gauge) at the flowline / choke
    surface_loss_coefficient: float = 0.0   # c in dp_surf = c Q^m (Pa per (m^3/s)^m)
    surface_loss_exponent: float = 1.8      # m
    clinging_constant: float = 0.45    # fraction of pipe velocity dragging the annular mud (surge/swab)


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
    """Drawworks limits for prescribed block velocity."""
    block_acceleration: float = 0.0   # m/s^2, optional extra limit on top of the ramp, 0 = none


@dataclass
class AutodrillerParameters:
    """PI autodriller on surface WOB: gains and maximum block speed."""
    kp: float = 1e-7            # (m/s) per N of WOB error
    ki: float = 1e-8            # (m/s) per N s
    v_max: float = 60.0 / 3600  # m/s


@dataclass
class ModelParameters:
    """
    All model parameters: friction, damping, bit, mud motor, top drive, controller,
    autodriller, drawworks, ramps, hydraulics, thermal and the CFL number.
    """
    friction: FrictionParameters = field(default_factory=FrictionParameters)
    damping: DampingParameters = field(default_factory=DampingParameters)
    bit: BitParameters = field(default_factory=BitParameters)
    motor: MotorParameters = field(default_factory=MotorParameters)
    top_drive: TopDriveParameters = field(default_factory=TopDriveParameters)
    controller: ControllerParameters = None
    autodriller: AutodrillerParameters = field(default_factory=AutodrillerParameters)
    drawworks: DrawworksParameters = field(default_factory=DrawworksParameters)
    ramps: RampParameters = field(default_factory=RampParameters)
    hydraulics: HydraulicsParameters = field(default_factory=HydraulicsParameters)
    thermal: ThermalParameters = field(default_factory=ThermalParameters)
    thermal_update_period: float = 2.0   # s between thermal / rock-strength updates
    cfl: float = 0.95           # Courant number on the axial wave speed

    def to_dict(self):
        """Plain-dict form."""
        return asdict(self)

    def copy(self, **changes):
        """Copy with some fields replaced (dataclasses.replace)."""
        return replace(self, **changes)
