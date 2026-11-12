"""
Formation description and rock strength at the bit.

Rock strength enters the bit models through the drilling strength
S = zeta * epsilon, with the intrinsic specific energy taken as the confined
compressive strength at the bit divided by a cutting efficiency eta::

    epsilon = CCS_eff / eta

CCS_eff combines three effects (see docs/derivation/rock-strength-and-temperature.md)::

    P1  Confinement by the overbalance, Mohr-Coulomb
            CCS = UCS + dp * (1 + sin phi) / (1 - sin phi)
        with dp = p_bh - p_p for permeable rock, and for low-permeability /
        low-porosity rock the undrained ("Skempton") pore pressure
            p_p* = p_p - B (sigma_v - p_bh) / 3
        blended by effective porosity (Caicedo, Calhoun & Ewy 2005; US 7,412,331
        and US 7,555,414).
    T1  Thermal weakening of the rock
            UCS(T) = UCS_0 [1 - beta max(T - T*, 0)]
        fitted to granite data of Xu et al. (2019): T* = 193 C, beta = 1.41e-3 /K.
    T2  Cooling-induced thermal stress at the hole bottom
            sigma_th = alpha E (T_rock - T_mud) / (1 - nu)
            CCS_eff  = CCS - kappa sigma_th min(1, delta / d_b)
        where delta = sqrt(pi D_th t_pass) is the thermal skin cooled between
        blade passes and d_b the depth of cut per blade.

Temperatures are in deg C, pressures in Pa, depths in m.
"""

from dataclasses import dataclass, field, asdict
from typing import List
import numpy as np

from . import units


@dataclass
class FormationLayer:
    """A formation interval starting at ``top_md`` (m MD) and extending to the next layer."""
    top_md: float
    ucs: float = 110e6               # Pa
    friction_angle: float = 36.0     # degrees
    porosity: float = 0.005          # effective porosity (fraction)
    name: str = ''


@dataclass
class Formation:
    """
    Formation, pore-pressure and geothermal description.

    Defaults are for the Utah FORGE granitoid:
      UCS 91-126 MPa (Ye et al. 2022, ARMA-2022-0272), friction angle 36 deg
      (Finnila & Podgorney 2020), matrix porosity < 0.5 %, E ~ 55 GPa,
      nu ~ 0.26, alpha = 6e-6 /K (Podgorney et al. 2021; compiled in the Utah
      FORGE decade review, Processes 2026, 14, 512); geothermal 40 C surface
      and 69.1 C/km, rock density 2800 kg/m3, heat capacity 930 J/kg K,
      conductivity 2.31 W/m K (Wang et al. 2023, SPE-214836-MS, Table 2).
    """
    layers: List[FormationLayer] = field(default_factory=lambda: [FormationLayer(0.0, name='Granitoid')])
    # pressures
    pore_pressure_gradient: float = 9.79e3      # Pa/m (fresh-water hydrostatic)
    overburden_gradient: float = 25.0e3         # Pa/m
    skempton_b: float = 1.0                     # undrained pore-pressure response (1 in Caicedo et al.)
    porosity_low: float = 0.05                  # at or below: low-permeability (undrained) CCS
    porosity_high: float = 0.20                 # at or above: permeable (drained) CCS
    # geothermal and rock thermal properties
    surface_temperature: float = 40.0           # C
    geothermal_gradient: float = 0.0691         # C/m
    rock_density: float = 2800.0                # kg/m3
    rock_heat_capacity: float = 930.0           # J/(kg K)
    rock_conductivity: float = 2.31             # W/(m K)
    # thermoelastic properties
    thermal_expansion: float = 6.0e-6           # 1/K (linear)
    youngs_modulus: float = 55e9                # Pa
    poisson_ratio: float = 0.26
    # T1: thermal weakening (fit to Xu et al. 2019 granite data)
    weakening_threshold: float = 193.0          # C
    weakening_rate: float = 1.41e-3             # 1/K
    # T2: cooling thermal stress
    thermal_stress_efficiency: float = 0.5      # kappa, 0..1 (1 = full thermoelastic stress)
    strength_floor: float = 0.1                 # CCS_eff >= floor * UCS(T)

    # ----------------------------------------------------------- profiles
    def layer_at(self, md):
        """The layer containing measured depth `md` (the deepest layer with top <= md)."""
        lay = sorted(self.layers, key=lambda l: l.top_md)
        out = lay[0]
        for l in lay:
            if md >= l.top_md:
                out = l
        return out

    def temperature(self, tvd):
        """Undisturbed formation temperature (C) at TVD (m)."""
        return self.surface_temperature + self.geothermal_gradient * np.asarray(tvd)

    def pore_pressure(self, tvd):
        """Hydrostatic pore pressure at true vertical depth `tvd` (Pa)."""
        return self.pore_pressure_gradient * np.asarray(tvd)

    def overburden(self, tvd):
        """Overburden stress at true vertical depth `tvd` (Pa)."""
        return self.overburden_gradient * np.asarray(tvd)

    @property
    def thermal_diffusivity(self):
        """Rock thermal diffusivity k / (rho c) (m^2/s)."""
        return self.rock_conductivity / (self.rock_density * self.rock_heat_capacity)

    # ------------------------------------------------------- rock strength
    def strength(self, md, tvd, p_bh, t_mud=None, rop=0.0, omega=0.0, blades=5):
        """
        Rock strength at the bit.

        Parameters
        ----------
        md, tvd : float
            Hole depth (m MD) and its TVD (m).
        p_bh : float
            Bottomhole pressure (Pa), e.g. from the circulating thermal/hydraulic model.
        t_mud : float, optional
            Mud temperature at the bit (C).  None disables the thermal-stress term.
        rop, omega : float
            Penetration rate (m/s) and bit speed (rad/s), for the thermal skin ratio.
        blades : int
            Number of blades (sets the time between blade passes).

        Returns
        -------
        dict with ucs0, t_rock, ucs_t, dp_drained, dp_undrained, w_drained,
        ccs (P1 + T1), sigma_th, skin_ratio, ccs_eff (P1 + T1 + T2).
        """
        lay = self.layer_at(md)
        t_rock = float(self.temperature(tvd))
        # T1
        ucs_t = lay.ucs * max(1.0 - self.weakening_rate * max(t_rock - self.weakening_threshold, 0.0), 0.05)
        # P1: drained and undrained (Skempton) overbalance, blended by porosity
        pp = float(self.pore_pressure(tvd))
        sv = float(self.overburden(tvd))
        dp_hp = p_bh - pp
        pp_lp = pp - self.skempton_b * (sv - p_bh) / 3.0
        dp_lp = p_bh - pp_lp
        phi_e = lay.porosity
        lo, hi = self.porosity_low, self.porosity_high
        w = 0.0 if phi_e <= lo else 1.0 if phi_e >= hi else (phi_e - lo) / (hi - lo)
        s = np.sin(np.radians(lay.friction_angle))
        n_phi = (1 + s) / (1 - s)
        ccs = ucs_t + ((1 - w) * dp_lp + w * dp_hp) * n_phi
        # T2
        sigma_th = 0.0
        ratio = 0.0
        if t_mud is not None:
            sigma_th = self.thermal_expansion * self.youngs_modulus * (t_rock - t_mud) / (1.0 - self.poisson_ratio)
            ratio = 1.0
            if omega > 1e-6 and rop > 1e-9:
                t_pass = 2.0 * np.pi / (max(int(blades), 1) * omega)
                delta = np.sqrt(np.pi * self.thermal_diffusivity * t_pass)
                d_b = rop * t_pass
                ratio = min(1.0, delta / d_b)
        ccs_eff = max(ccs - self.thermal_stress_efficiency * sigma_th * ratio, self.strength_floor * ucs_t)
        return dict(ucs0=lay.ucs, t_rock=t_rock, ucs_t=ucs_t, n_phi=n_phi, dp_drained=dp_hp,
                    dp_undrained=dp_lp, w_drained=w, ccs=ccs, sigma_th=sigma_th, skin_ratio=ratio,
                    ccs_eff=ccs_eff)

    # ------------------------------------------------------------------ io
    def to_dict(self):
        """Plain-dict form (for JSON cases)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        """Formation from its dict form; missing values take the defaults."""
        d = dict(d or {})
        layers = [FormationLayer(**l) for l in d.pop('layers', [])] or None
        f = cls(**d)
        if layers:
            f.layers = layers
        return f
