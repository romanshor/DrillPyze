"""
Transient wellbore temperature while drilling, after the modified drift-flux
model of Wang, Khaled, Luu, Ashok & van Oort (2023), "Downhole Temperature
Estimation of a Growing High-Temperature Wellbore Using a Modified Drift
Flux Modeling Approach", SPE-214836-MS, doi:10.2118/214836-MS.

Implemented here in its single-phase (liquid mud) form:

* drillstring and annulus fluid energy balances on a common axial grid,
  upwind in the flow direction (paper eqs. 3-5, 13);
* radial heat flow drillstring <-> annulus and annulus <-> formation through
  series conduction/convection resistances (eqs. 6-9; Nusselt number from
  Dittus-Boelter for turbulent flow, 3.66 for laminar flow);
* the formation as radial conduction cells per axial cell with the
  undisturbed geothermal temperature in the far field (Khaled et al. 2023c);
* a growing bottom cell, grown at the ROP and cut in two when it reaches a
  threshold (the "grow-cut-grow" strategy, Fig. 2), with the cuttings
  entering at the formation temperature (eqs. 14-20) and the bottom formation
  cells relaxed toward the undisturbed temperature as the hole deepens
  (eq. 10);
* mud density rho(T, p) (eq. 2) integrated up the annulus with the frictional
  pressure loss to give the bottomhole pressure and ECD.

Simplifications relative to the paper, all consistent with its assumptions
1-3: single liquid phase (gas volume fraction zero), velocity from the pump
rate (ROP << mud velocity), axial conduction in the fluid neglected,
transient heat storage in the drillstring steel neglected, cuttings carry
their heat into the mud at the bottom cell.  Moving geometry discontinuities
are handled by volume-averaging the string geometry over each cell at every
update (Fig. 1); the associated ROP-scaled area-change terms of eqs. 11-13
are O(ROP/v) and omitted.

Units: SI, temperatures in deg C.
"""

from dataclasses import dataclass
import numpy as np
from numba import njit

from . import units
from .hydraulics import frictional_pressure_drop


@dataclass
class ThermalParameters:
    enabled: bool = True
    dx: float = 30.0                # nominal axial cell length (m)
    cut_factor: float = 1.5         # bottom cell is split when longer than cut_factor * dx
    n_radial: int = 12              # formation cells per axial cell
    far_field: float = 5.0          # radial extent of the modelled formation beyond the wall (m)
    steel_conductivity: float = 45.0
    casing_thickness: float = 0.012 # m, used for conduction through casing
    flow_rate: float = 0.0347       # default pump rate (m^3/s), 550 gpm
    inlet_temperature: float = 30.0 # default mud inlet temperature (C)
    max_dt: float = 2.0             # s, largest thermal time step
    gs_iterations: int = 3          # fluid / formation coupling iterations per step


def _h(Q, A, Dh, fluid):
    """Convection coefficient (W/m^2 K): Dittus-Boelter (turbulent) or Nu = 3.66 (laminar)."""
    k, mu, cp, rho = fluid.thermal_conductivity, fluid.viscosity, fluid.heat_capacity, fluid.density
    if A <= 0 or Dh <= 0:
        return 0.0
    v = abs(Q) / A
    Re = rho * v * Dh / mu
    Pr = mu * cp / k
    Nu = 0.023 * Re ** 0.8 * Pr ** 0.4 if Re > 2300 else 3.66
    return max(Nu, 3.66) * k / Dh


@njit(cache=True)
def _solve_layer(T, rc_vol, g, g_in, T_in_side, g_out, T_far, dt):
    """Backward-Euler radial conduction in one formation layer (tridiagonal)."""
    m = T.shape[0]
    a = np.zeros(m)
    b = np.zeros(m)
    c = np.zeros(m)
    d = np.zeros(m)
    for i in range(m):
        C = rc_vol[i] / dt
        b[i] = C
        d[i] = C * T[i]
        if i == 0:
            b[i] += g_in
            d[i] += g_in * T_in_side
        else:
            a[i] = -g[i - 1]
            b[i] += g[i - 1]
        if i == m - 1:
            b[i] += g_out
            d[i] += g_out * T_far
        else:
            c[i] = -g[i]
            b[i] += g[i]
    # Thomas algorithm
    for i in range(1, m):
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    T[m - 1] = d[m - 1] / b[m - 1]
    for i in range(m - 2, -1, -1):
        T[i] = (d[i] - c[i] * T[i + 1]) / b[i]


@njit(cache=True)
def _thermal_step(dt, F, T_in, jb, nb, q_cut_g, T_cut,
                  T_ds, T_ann, T_f, C_ds, C_ann, G_ds, G_w, rc_vol, g_rad, g_out, T_far, n_iter):
    """
    One backward-Euler step of the coupled fluid / formation energy balance.

    F       mud heat-capacity rate rho c Q (W/K)
    jb      cell containing the bit (flow turns from drillstring to annulus there)
    nb      bottom (growing) cell, receives the cuttings
    q_cut_g cuttings heat-capacity rate rho_r c_r A_hole ROP (W/K), entering at T_cut
    Returns the heat into the fluid from the formation (W) at the end of the step.
    """
    n = T_ann.shape[0]
    T_ds_old = T_ds.copy()
    T_ann_old = T_ann.copy()
    T_f_old = T_f.copy()
    for it in range(n_iter):
        # drillstring, flowing down (cells 0..jb)
        up = T_in
        for j in range(jb + 1):
            Cd = C_ds[j] / dt
            num = Cd * T_ds_old[j] + F * up + G_ds[j] * T_ann[j]
            den = Cd + F + G_ds[j]
            T_ds[j] = num / den if den > 0 else T_ds_old[j]
            up = T_ds[j]
        # annulus, flowing up from the bit (cells below the bit are stagnant)
        below = T_ds[jb]
        for j in range(n - 1, -1, -1):
            Ca = C_ann[j] / dt
            Fj = F if j <= jb else 0.0
            gds = G_ds[j] if j <= jb else 0.0
            num = Ca * T_ann_old[j] + gds * T_ds[j] + G_w[j] * T_f[j, 0]
            den = Ca + gds + G_w[j]
            if Fj > 0.0:
                num += Fj * below
                den += Fj
            if j == nb:
                num += q_cut_g * T_cut
                den += q_cut_g
            T_ann[j] = num / den
            if j <= jb:
                below = T_ann[j]
        # formation layers, coupled to the annulus through the wall conductance
        for j in range(n):
            T_f[j, :] = T_f_old[j, :]
            _solve_layer(T_f[j], rc_vol[j], g_rad[j], G_w[j], T_ann[j], g_out[j], T_far[j], dt)
    q = 0.0
    for j in range(n):
        q += G_w[j] * (T_f[j, 0] - T_ann[j])
    return q


class WellboreThermal:
    """
    Wellbore temperature model coupled to a :class:`Simulation`.

    Call :meth:`advance` with the elapsed time, pump rate, inlet temperature,
    bit and hole depth and ROP.  The grid follows the hole as it is drilled.
    """

    def __init__(self, wellbore, drillstring, fluid, formation, params: ThermalParameters,
                 bit_depth, hole_depth, initial='geothermal'):
        self.wellbore = wellbore
        self.drillstring = drillstring
        self.fluid = fluid
        self.formation = formation
        self.p = params
        self.time = 0.0
        self.bit_depth = bit_depth
        self.hole_depth = hole_depth
        self.flow_rate = params.flow_rate
        self.inlet_temperature = params.inlet_temperature
        self._build_grid(hole_depth)
        self._geometry(bit_depth, force=True)
        tf = self.T_far
        self.T_ds = tf.copy()
        self.T_ann = tf.copy()
        self.T_f = np.repeat(tf[:, None], params.n_radial, axis=1)
        if initial not in ('geothermal',):
            raise ValueError('initial must be "geothermal"')
        self.q_formation = 0.0
        self._pressure()

    # ----------------------------------------------------------------- grid
    def _build_grid(self, hole_depth):
        dx = self.p.dx
        n = max(1, int(np.floor(hole_depth / dx)))
        edges = np.linspace(0.0, n * dx, n + 1)
        if hole_depth - edges[-1] > 1e-9:
            if hole_depth - edges[-1] < 0.5 * dx and n > 1:
                edges[-1] = hole_depth
            else:
                edges = np.append(edges, hole_depth)
        self.x_top = edges[:-1].copy()
        self.dx = np.diff(edges)

    def _radial(self, r_w):
        m = self.p.n_radial
        f = r_w * ((r_w + self.p.far_field) / r_w) ** (np.arange(m + 1) / m)   # face radii
        rc = np.sqrt(f[:-1] * f[1:])
        k = self.formation.rock_conductivity
        vol = np.pi * (f[1:] ** 2 - f[:-1] ** 2) * self.formation.rock_density * self.formation.rock_heat_capacity
        g = 2 * np.pi * k / np.log(rc[1:] / rc[:-1])          # between nodes, per length
        g_out = 2 * np.pi * k / np.log(f[-1] / rc[-1])
        r_in = np.log(rc[0] / r_w) / (2 * np.pi * k)          # wall to first node, per length
        return rc, vol, g, g_out, r_in

    def _geometry(self, bit_depth, force=False):
        """Volume-averaged geometry and heat-transfer conductances for every cell."""
        n = len(self.dx)
        key = (n, round(float(self.dx[-1]), 2), self.flow_rate)
        if (not force and getattr(self, '_geo_key', None) == key
                and abs(bit_depth - self._geo_bit) < 0.25):
            return
        self._geo_key, self._geo_bit = key, bit_depth
        wb, ds = self.wellbore, self.drillstring
        fl = self.fluid
        mid = self.x_top + 0.5 * self.dx
        self.tvd = wb.survey.tvd_at(mid)
        self.T_far = np.asarray(self.formation.temperature(self.tvd), dtype=float)
        # string geometry, sampled along each cell (volume averaged, Fig. 1)
        top_len = bit_depth - ds.bha_length
        comps = list(reversed(ds.components))
        lens = [top_len] + [c.length for c in comps[1:]]
        bnds = np.concatenate([[0.0], np.cumsum(lens)])
        k_s = self.p.steel_conductivity
        Q = self.flow_rate
        A_ds = np.zeros(n); A_ann = np.zeros(n); G_ds = np.zeros(n); G_w = np.zeros(n)
        rc_vol = np.zeros((n, self.p.n_radial)); g_rad = np.zeros((n, self.p.n_radial - 1)); g_out = np.zeros(n)
        self.D_wall = np.zeros(n)
        ns = 8
        xs = self.x_top[:, None] + (np.arange(ns) + 0.5)[None, :] / ns * self.dx[:, None]
        ic = np.clip(np.searchsorted(bnds, xs, side='right') - 1, 0, len(comps) - 1)
        od_c = np.array([c.od for c in comps]); id_c = np.array([c.id for c in comps])
        inside = xs <= bit_depth
        od = np.where(inside, od_c[ic], 0.0); idd = np.where(inside, id_c[ic], 0.0)
        cased = mid <= wb.casing_shoe
        Dw = np.where(cased, (wb.casing_id or wb.hole_diameter), wb.hole_diameter)
        self.D_wall = Dw
        A_ds = (np.pi / 4 * idd ** 2).mean(axis=1)
        A_ann = np.pi / 4 * Dw ** 2 - (np.pi / 4 * od ** 2).mean(axis=1)
        Di, Do = idd.mean(axis=1), od.mean(axis=1)
        radial = {}
        for j in range(n):
            h_ann = _h(Q, A_ann[j], Dw[j] - Do[j], fl)
            if Do[j] > 0 and Di[j] > 0:
                h_i = _h(Q, A_ds[j], Di[j], fl)
                R_ds = 1 / (h_i * np.pi * Di[j]) + np.log(Do[j] / Di[j]) / (2 * np.pi * k_s) + 1 / (h_ann * np.pi * Do[j])
                G_ds[j] = self.dx[j] / R_ds
            r_w = Dw[j] / 2 + (self.p.casing_thickness if cased[j] else 0.0)
            if r_w not in radial:
                radial[r_w] = self._radial(r_w)
            rc, vol, g, go, r_in = radial[r_w]
            R_w = 1 / (h_ann * np.pi * Dw[j]) + r_in
            if cased[j]:
                R_w += np.log((Dw[j] + 2 * self.p.casing_thickness) / Dw[j]) / (2 * np.pi * k_s)
            G_w[j] = self.dx[j] / R_w
            rc_vol[j] = vol * self.dx[j]
            g_rad[j] = g * self.dx[j]
            g_out[j] = go * self.dx[j]
        rc_ = fl.density * fl.heat_capacity
        self.A_ds, self.A_ann = A_ds, A_ann
        self.C_ds, self.C_ann = rc_ * A_ds * self.dx, rc_ * A_ann * self.dx
        self.G_ds, self.G_w = G_ds, G_w
        self.rc_vol, self.g_rad, self.g_out = rc_vol, g_rad, g_out
        self.jb = int(min(np.searchsorted(self.x_top, bit_depth, side='right') - 1, n - 1))
        self.nb = n - 1

    # ------------------------------------------------------------ growth
    def _grow(self, hole_depth, rop, dt):
        dh = hole_depth - (self.x_top[-1] + self.dx[-1])
        if dh <= 0:
            return
        self.dx[-1] += dh
        # eq. 10: bottom formation cells relax toward the undisturbed temperature
        tvd_b = self.wellbore.survey.tvd_at([hole_depth])[0]
        T_far_b = float(self.formation.temperature(tvd_b))
        w = min(1.0, dh / self.dx[-1])
        self.T_f[-1] += w * (T_far_b - self.T_f[-1])
        # grow-cut-grow (Fig. 2): split the bottom cell, both halves keep its state
        if self.dx[-1] > self.p.cut_factor * self.p.dx:
            half = 0.5 * self.dx[-1]
            self.dx[-1] = half
            self.x_top = np.append(self.x_top, self.x_top[-1] + half)
            self.dx = np.append(self.dx, half)
            self.T_ds = np.append(self.T_ds, self.T_ds[-1])
            self.T_ann = np.append(self.T_ann, self.T_ann[-1])
            self.T_f = np.vstack([self.T_f, self.T_f[-1:]])

    # ---------------------------------------------------------- pressure
    def _pressure(self):
        """Annular pressure profile from rho(T, p) (eq. 2) plus frictional loss."""
        fl = self.fluid
        n = len(self.dx)
        tvd_edges = self.wellbore.survey.tvd_at(np.concatenate([self.x_top, [self.x_top[-1] + self.dx[-1]]]))
        p = 0.0
        self.p_ann = np.zeros(n + 1)
        Q = abs(self.flow_rate)
        for j in range(n):
            rho = fl.density_at(self.T_ann[j], p)
            dp = rho * units.G * (tvd_edges[j + 1] - tvd_edges[j])
            if Q > 0 and j <= self.jb and self.A_ann[j] > 0:
                Do = np.sqrt(max(self.D_wall[j] ** 2 - 4 * self.A_ann[j] / np.pi, 0.0))
                dp += frictional_pressure_drop(rho, Q, fl.viscosity, self.dx[j], self.D_wall[j], Do if Do > 0 else None)
            p += dp
            self.p_ann[j + 1] = p
        self.p_bh = p
        tvd_b = tvd_edges[-1]
        self.ecd = p / (units.G * tvd_b) if tvd_b > 0 else fl.density

    # ------------------------------------------------------------ advance
    def advance(self, duration, flow_rate=None, inlet_temperature=None, bit_depth=None,
                hole_depth=None, rop=0.0):
        """Advance the temperature field by ``duration`` seconds."""
        if flow_rate is not None:
            self.flow_rate = float(flow_rate)
        if inlet_temperature is not None:
            self.inlet_temperature = float(inlet_temperature)
        if bit_depth is not None:
            self.bit_depth = bit_depth
        hole = self.hole_depth if hole_depth is None else max(hole_depth, self.hole_depth)
        nsub = max(1, int(np.ceil(duration / self.p.max_dt)))
        dt = duration / nsub
        fl, fm = self.fluid, self.formation
        hole0 = self.hole_depth
        for k in range(nsub):
            h_k = hole0 + (hole - hole0) * (k + 1) / nsub
            self._grow(h_k, rop, dt)
            self._geometry(self.bit_depth, force=(rop > 0))
            F = fl.density * fl.heat_capacity * abs(self.flow_rate)
            on_bottom = self.bit_depth >= h_k - 0.5
            A_hole = np.pi / 4 * self.wellbore.hole_diameter ** 2
            q_cut_g = fm.rock_density * fm.rock_heat_capacity * A_hole * max(rop, 0.0) if on_bottom else 0.0
            T_cut = float(fm.temperature(self.wellbore.survey.tvd_at([h_k])[0]))
            self.q_formation = _thermal_step(
                dt, F, self.inlet_temperature, self.jb, self.nb, q_cut_g, T_cut,
                self.T_ds, self.T_ann, self.T_f, self.C_ds, self.C_ann, self.G_ds, self.G_w,
                self.rc_vol, self.g_rad, self.g_out, self.T_far, self.p.gs_iterations)
            self.q_cuttings = q_cut_g * (T_cut - self.T_ann[self.nb])
            self.time += dt
        self.hole_depth = hole
        self._pressure()
        return self.state()

    # ------------------------------------------------------------- output
    def state(self):
        return dict(time=self.time, t_bit=float(self.T_ann[self.jb]), t_ds_bit=float(self.T_ds[self.jb]),
                    t_out=float(self.T_ann[0]), t_in=self.inlet_temperature,
                    t_rock=float(self.formation.temperature(self.wellbore.survey.tvd_at([self.hole_depth])[0])),
                    p_bh=float(self.p_bh), ecd=float(self.ecd), flow_rate=self.flow_rate)

    def profiles(self):
        import pandas as pd
        return pd.DataFrame(dict(md=self.x_top + 0.5 * self.dx, tvd=self.tvd, t_ds=self.T_ds,
                                 t_ann=self.T_ann, t_wall=self.T_f[:, 0], t_far=self.T_far,
                                 has_pipe=(self.A_ds > 0)))

    def stored_energy(self):
        """Thermal energy relative to 0 C in fluid and formation (J), for balance checks."""
        return float(np.sum(self.C_ds * self.T_ds) + np.sum(self.C_ann * self.T_ann) + np.sum(self.rc_vol * self.T_f))
