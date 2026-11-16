"""
Borehole hydraulics (SI units).

Steady state
    Newtonian and Bingham-plastic pressure losses in pipe and annular flow,
    bit nozzle loss, standpipe pressure, bottomhole pressure and ECD, as in
    Shor, *Introduction to Drilling Engineering*, ch. 7 (Borehole Hydraulics)::

        * apparent viscosity: pipe  mu_a = PV + tau_y D / (6 v)
                              annulus mu_a = PV + tau_y (D2 - D1) / (8 v)
        * Reynolds number:    pipe  Re = rho v D / mu_a
                              annulus Re = rho v 0.816 (D2 - D1) / mu_a
          laminar below 2100, turbulent (transition treated as turbulent) above
        * laminar:   pipe     dp/dL = 32 PV v / D^2  + 16 tau_y / (3 D)
                     annulus  dp/dL = 48 PV v / h^2  + 6 tau_y / h,       h = D2 - D1
                     (slot flow, 12 PV v / g^2 + 3 tau_y / g with the gap g = h / 2)
        * turbulent: Fanning, Blasius f = 0.0791 Re^-0.25 with Re on PV,
                     dp/dL = 2 f rho v^2 / D_e,  D_e = D (pipe), 0.816 h (annulus)
        * bit:       dp = rho Q^2 / (2 C_d^2 A_n^2)
        * SPP = dp_surface + dp_string + dp_bit + dp_annulus;  BHP = p_hydrostatic + dp_annulus

    (These are the chapter's field-unit equations written in SI.)  A Newtonian
    fluid is the special case tau_y = 0, PV = mu.

Transient
    :class:`WellboreHydraulics` -- the single-phase limit of the drift-flux
    model of Wang, Khaled, Luu, Ashok & van Oort (2023, SPE-214836-MS; after
    Evje & Fjelde 2002): mass and momentum along the drillstring and annulus
    with rho(T, p) of their eq. 2, a volume-rate inlet, a constant-pressure
    outlet, fixed cells with volume-averaged cross sections, and the bit
    nozzles linking the two flow paths.  See docs/derivation/hydraulics.md.
"""

import numpy as np
from numba import njit

G = 9.80665
RE_LAMINAR = 2100.0
ANNULUS_FACTOR = 0.816      # equivalent diameter factor for annular Reynolds number / turbulent loss


# ====================================================================== legacy
def area(Do, Di=None):
    """Flow area (m^2) of a pipe (``Di`` is None) or annulus."""
    return np.pi / 4 * (Do ** 2 - (0.0 if Di is None else Di ** 2))


def fluid_velocity(Q, Do, Di=None):
    """Mean velocity (m/s) for flow rate ``Q`` (m^3/s)."""
    return Q / area(Do, Di)


def hydraulic_diameter(Do, Di=None):
    """Hydraulic diameter: Do for a pipe, Do - Di for an annulus."""
    return Do if Di is None else Do - Di


def reynolds_number(rho, v, mu, Do, Di=None):
    """Newtonian Reynolds number on the hydraulic diameter."""
    return rho * abs(v) * hydraulic_diameter(Do, Di) / mu


def friction_factor(Re):
    """Fanning friction factor (16/Re laminar, Blasius turbulent)."""
    Re = max(Re, 1e-12)
    return 16.0 / Re if Re < RE_LAMINAR else 0.0791 / Re ** 0.25


def frictional_pressure_drop(rho, Q, mu, L, Do, Di=None):
    """Frictional pressure drop (Pa) over ``L`` (m) for a Newtonian fluid (chapter 7 forms)."""
    v = fluid_velocity(Q, Do, Di)
    return gradient(v, rho, mu, 0.0, Do, 0.0 if Di is None else Di) * L


# ================================================================ rheology
def rheology(fluid):
    """(PV, tau_y) of a :class:`~DrillPyze.wellbore.Fluid`: Bingham if a yield point is given."""
    if fluid.yield_point:
        pv = fluid.plastic_viscosity if fluid.plastic_viscosity is not None else fluid.viscosity
        return float(pv), float(fluid.yield_point)
    return float(fluid.viscosity), 0.0


@njit(cache=True)
def _coeffs(au, rho, pv, ty, d_out, d_in):
    """
    Friction gradient dp/dL = a u + b sign(u) at speed |u| = au.

    Pipe when d_in <= 0 (diameter d_out), else annulus (d_out = hole, d_in = pipe).
    Returns (a, b, turbulent).
    """
    if d_in <= 0.0:
        D = d_out
        a_lam = 32.0 * pv / (D * D)
        b_lam = 16.0 * ty / (3.0 * D)
        De = D
        mu_a = pv + (ty * D / (6.0 * au) if au > 0.0 else 1e30)
    else:
        h = d_out - d_in
        # slot flow with gap (D2 - D1)/2: 12 PV v / gap^2 + 3 tau_y / gap
        a_lam = 48.0 * pv / (h * h)
        b_lam = 6.0 * ty / h
        De = ANNULUS_FACTOR * h
        mu_a = pv + (ty * h / (8.0 * au) if au > 0.0 else 1e30)
    Re = rho * au * De / mu_a
    if Re < RE_LAMINAR:
        return a_lam, b_lam, 0
    Re_t = rho * au * De / pv
    f = 0.0791 / Re_t ** 0.25
    return 2.0 * f * rho * au / De, 0.0, 1


@njit(cache=True)
def gradient(u, rho, pv, ty, d_out, d_in):
    """Frictional pressure gradient (Pa/m), signed with the velocity ``u`` (m/s)."""
    au = abs(u)
    if au == 0.0:
        return 0.0
    a, b, _ = _coeffs(au, rho, pv, ty, d_out, d_in)
    g = a * au + b
    return g if u > 0.0 else -g


def section(Q, rho, pv, ty, length, d_out, d_in=0.0):
    """
    One textbook section: velocity, apparent viscosity, Reynolds number,
    regime and pressure loss.  Pipe if ``d_in`` is 0, annulus otherwise.
    """
    A = area(d_out, d_in if d_in > 0 else None)
    v = Q / A
    if d_in <= 0:
        De, mu_a = d_out, pv + ty * d_out / (6.0 * v)
    else:
        h = d_out - d_in
        De, mu_a = ANNULUS_FACTOR * h, pv + ty * h / (8.0 * v)
    Re = rho * v * De / mu_a
    a, b, turb = _coeffs(v, rho, pv, ty, d_out, d_in)
    return dict(velocity=v, apparent_viscosity=mu_a, reynolds=Re,
                regime='turbulent' if turb else 'laminar', dp=(a * v + b) * length)


def nozzle_area(sizes_32nds):
    """Total flow area (m^2) of nozzles given in 32nds of an inch."""
    d = np.asarray(sizes_32nds, dtype=float) / 32.0 * 0.0254
    return float(np.sum(np.pi / 4 * d ** 2))


def bit_pressure_drop(Q, rho, tfa, cd=0.95):
    """Bit nozzle pressure drop (Pa): rho Q^2 / (2 C_d^2 A^2)."""
    return rho * Q * abs(Q) / (2.0 * cd ** 2 * tfa ** 2)


def surface_loss(Q, coefficient=0.0, exponent=1.8):
    """Surface-equipment loss c Q^m (Pa); the chapter's empirical form."""
    return coefficient * abs(Q) ** exponent


def steady_state(Q, density, pv, ty, pipe, annulus, tfa, tvd, cd=0.95, surface=0.0):
    """
    Textbook steady-state circulating system (constant density).

    pipe     list of (name, inner diameter, length)
    annulus  list of (name, hole diameter, pipe outer diameter, length)
    tfa      total nozzle area (m^2); tvd (m); surface loss (Pa)

    Returns dict with a per-section table, bit loss, SPP, annular loss,
    hydrostatic pressure, BHP and ECD (kg/m^3).
    """
    import pandas as pd
    rows = []
    for name, D, L in pipe:
        rows.append(dict(section=name, kind='pipe', **section(Q, density, pv, ty, L, D)))
    for name, D2, D1, L in annulus:
        rows.append(dict(section=name, kind='annulus', **section(Q, density, pv, ty, L, D2, D1)))
    df = pd.DataFrame(rows)
    dp_bit = bit_pressure_drop(Q, density, tfa, cd)
    dp_string = float(df.loc[df.kind == 'pipe', 'dp'].sum())
    dp_ann = float(df.loc[df.kind == 'annulus', 'dp'].sum())
    p_hyd = density * G * tvd
    return dict(sections=df, bit=dp_bit, string=dp_string, annulus=dp_ann, surface=surface,
                spp=surface + dp_string + dp_bit + dp_ann, hydrostatic=p_hyd,
                bhp=p_hyd + dp_ann, ecd=density + dp_ann / (G * tvd))


# ============================================================ transient kernel
@njit(cache=True)
def _face(G, dt, inv_I, drive, rho_u, A_u, L_u, du_o, du_i, rho_d, A_d, L_d, dd_o, dd_i,
          pv, ty, v_slip, k_slip, extra_R, iters):
    """
    Backward-Euler momentum update of one face (mass flux G, kg/s, + down):

        dG/dt = (drive - dp_fric(G)) / I,   dp_fric = R (G - G_s) + B sign(G - G_s)

    R, B from the two half cells (and an optional quadratic nozzle term
    extra_R * |G|), relative to the pipe-dragged flux G_s = k_slip v_slip rho A.
    The yield term is a Coulomb inclusion: the face holds (G = G_s) while the
    drive cannot break the gel.
    """
    rho_f = 0.5 * (rho_u + rho_d)
    A_f = 0.5 * (A_u + A_d)
    G_s = k_slip * v_slip * rho_f * A_f
    Gk = G
    for _ in range(iters):
        R = 0.0
        B = 0.0
        if L_u > 0.0:
            ur = abs(Gk / (rho_u * A_u) - k_slip * v_slip)
            a, b, t = _coeffs(ur, rho_u, pv, ty, du_o, du_i)
            R += 0.5 * L_u * a / (rho_u * A_u)
            B += 0.5 * L_u * b
        if L_d > 0.0:
            ur = abs(Gk / (rho_d * A_d) - k_slip * v_slip)
            a, b, t = _coeffs(ur, rho_d, pv, ty, dd_o, dd_i)
            R += 0.5 * L_d * a / (rho_d * A_d)
            B += 0.5 * L_d * b
        R += extra_R * abs(Gk)
        c = dt * inv_I
        free = G + c * (drive + R * G_s)
        den = 1.0 + c * R
        # stick (relative to the pipe) if the drive cannot overcome the yield stress
        rel = free - den * G_s
        if abs(rel) <= c * B:
            Gk = G_s
        else:
            s = 1.0 if rel > 0.0 else -1.0
            Gk = (free - c * B * s) / den
    return Gk


@njit(cache=True)
def _advance(nsub, dt, c2, pv, ty, k_ann, v_pipe, q_in_a, q_in_b, rho_in, p_back, nozzle_R, dp_motor, iters,
             # drillstring path (n_d cells), faces G_d[0..n_d-1] (face k = top of cell k), bit face Gb
             m_d, V_d0, V_d1, L_d, z_d, rT_d, D_d, G_d,
             # annulus path (n_a cells), faces G_a[0..n_a-1] (face 0 = outlet), bit enters cell kab
             m_a, V_a0, V_a1, L_a, z_a, rT_a, Do_a, Di_a, G_a, kab, Gb_arr, acc):
    """
    Explicit substeps of the transient model (numba kernel).

    Advances the cell masses and face mass fluxes of the drillstring and annulus
    paths over `nsub` substeps of length `dt`: pump inlet face, interior bore
    faces, the bit (nozzle) face -- with the mud motor differential `dp_motor` in series --,
    interior annulus faces and the outlet face, each
    updated with `_face`, then the mass balance. Arrays are updated in place;
    `acc` accumulates the pumped and returned mass.
    """
    n_d = m_d.shape[0]
    n_a = m_a.shape[0]
    p_d = np.empty(n_d)
    p_a = np.empty(n_a)
    r_d = np.empty(n_d)
    r_a = np.empty(n_a)
    A_d = np.empty(n_d)
    A_a = np.empty(n_a)
    Gb = Gb_arr[0]
    for s in range(nsub):
        w1 = (s + 1.0) / nsub
        # pump rate ramps linearly over the interval
        G_d[0] = rho_in * (q_in_a + (q_in_b - q_in_a) * w1)
        for j in range(n_d):
            V = V_d0[j] + (V_d1[j] - V_d0[j]) * w1
            A_d[j] = V / L_d[j]
            r_d[j] = m_d[j] / V
            p_d[j] = c2 * (r_d[j] - rT_d[j])
        for j in range(n_a):
            V = V_a0[j] + (V_a1[j] - V_a0[j]) * w1
            A_a[j] = V / L_a[j]
            r_a[j] = m_a[j] / V
            p_a[j] = c2 * (r_a[j] - rT_a[j])
        # ---- faces inside the drillstring (k_slip = 1: fluid moves with the pipe wall)
        for k in range(1, n_d):
            drive = p_d[k - 1] - p_d[k] + 0.5 * (r_d[k - 1] + r_d[k]) * G * (z_d[k] - z_d[k - 1])
            inv_I = 1.0 / (0.5 * (L_d[k - 1] / A_d[k - 1] + L_d[k] / A_d[k]))
            G_d[k] = _face(G_d[k], dt, inv_I, drive,
                           r_d[k - 1], A_d[k - 1], L_d[k - 1], D_d[k - 1], 0.0,
                           r_d[k], A_d[k], L_d[k], D_d[k], 0.0,
                           pv, ty, v_pipe, 1.0, 0.0, iters)
        # ---- bit nozzles: DS bottom cell -> annulus cell kab
        jd = n_d - 1
        drive = p_d[jd] - p_a[kab] + 0.5 * (r_d[jd] + r_a[kab]) * G * (z_a[kab] - z_d[jd]) - dp_motor
        inv_I = 1.0 / (0.5 * (L_d[jd] / A_d[jd] + L_a[kab] / A_a[kab]))
        Gb = _face(Gb, dt, inv_I, drive,
                   r_d[jd], A_d[jd], L_d[jd], D_d[jd], 0.0,
                   r_a[kab], A_a[kab], L_a[kab], Do_a[kab], Di_a[kab],
                   pv, ty, 0.0, 0.0, nozzle_R / (0.5 * (r_d[jd] + r_a[kab])), iters)
        # ---- annulus faces (k_slip = clinging constant)
        for k in range(1, n_a):
            drive = p_a[k - 1] - p_a[k] + 0.5 * (r_a[k - 1] + r_a[k]) * G * (z_a[k] - z_a[k - 1])
            inv_I = 1.0 / (0.5 * (L_a[k - 1] / A_a[k - 1] + L_a[k] / A_a[k]))
            ks = k_ann if Di_a[k] > 0.0 and Di_a[k - 1] > 0.0 else 0.0
            G_a[k] = _face(G_a[k], dt, inv_I, drive,
                           r_a[k - 1], A_a[k - 1], L_a[k - 1], Do_a[k - 1], Di_a[k - 1],
                           r_a[k], A_a[k], L_a[k], Do_a[k], Di_a[k],
                           pv, ty, v_pipe, ks, 0.0, iters)
        # ---- outlet (flowline at p_back, depth 0)
        drive = p_back - p_a[0] + r_a[0] * G * z_a[0]
        inv_I = 1.0 / (0.5 * L_a[0] / A_a[0])
        G_a[0] = _face(G_a[0], dt, inv_I, drive,
                       r_a[0], A_a[0], 0.0, Do_a[0], Di_a[0],
                       r_a[0], A_a[0], L_a[0], Do_a[0], Di_a[0],
                       pv, ty, v_pipe, k_ann if Di_a[0] > 0.0 else 0.0, 0.0, iters)
        # ---- mass (acc: cumulative mass pumped in, and returned out of the flowline)
        acc[0] += dt * G_d[0]
        acc[1] -= dt * G_a[0]
        for j in range(n_d - 1):
            m_d[j] += dt * (G_d[j] - G_d[j + 1])
        m_d[jd] += dt * (G_d[jd] - Gb)
        for j in range(n_a - 1):
            m_a[j] += dt * (G_a[j] - G_a[j + 1])
        m_a[n_a - 1] += dt * G_a[n_a - 1]
        m_a[kab] += dt * Gb
    Gb_arr[0] = Gb


# ============================================================== the model
class WellboreHydraulics:
    """
    Transient single-phase drift-flux hydraulics of the circulating system.

    Two flow paths on fixed cells (cell edges fixed in measured depth): the
    drillstring bore from surface to the bit, and the annulus (plus rat hole)
    from surface to the hole bottom.  Cross sections are volume-averaged over
    the string occupying each cell (Wang et al. 2023, Fig. 1); the bottom
    annulus cell grows with the hole and the bottom drillstring cell follows
    the bit (grow-cut-grow, Fig. 2).

    State: mass per cell (both paths) and mass flux per face.  Pressure from
    rho(T, p) of Wang et al. eq. 2: p = c^2 (m/V - rho_T).  Inlet: pump volume
    rate.  Outlet: flowline at ``back_pressure``.  The bit nozzles are a face
    between the bottom drillstring cell and the annulus cell at the bit.

    ``advance(t0, t1, ...)`` integrates with an explicit acoustic step on the
    pressures and a backward-Euler, Picard-iterated step on the face fluxes
    (friction, nozzle loss and Bingham yield stress implicit).
    """

    def __init__(self, wellbore, grid, fluid, params, tfa, cd, bit_depth, hole_depth, bit_depth0):
        self.wellbore = wellbore
        self.fluid = fluid
        self.p = params
        self.tfa = float(tfa)
        self.dp_motor = 0.0                 # mud motor differential in series with the nozzles (Pa)
        self.cd = float(cd)
        self.pv, self.ty = rheology(fluid)
        self.c = float(fluid.sound_speed)
        self.dx = float(params.dx)
        # string profile relative to its top (cumulative areas for volume averaging)
        self._s_edges = np.concatenate([[0.0], np.cumsum(grid.dx)])
        self._L_string = float(self._s_edges[-1])
        self._cum_od = np.concatenate([[0.0], np.cumsum(np.pi / 4 * grid.od ** 2 * grid.dx)])
        self._cum_bore = np.concatenate([[0.0], np.cumsum(np.pi / 4 * grid.id ** 2 * grid.dx)])
        self.bit_depth = float(bit_depth)
        self.hole_depth = float(hole_depth)
        self.e_a = self._edges(self.hole_depth)
        self.e_d = self._edges(self.bit_depth)
        self.T_d = None
        self.T_a = None
        self._geometry()
        self.V_d = self._Vd.copy()
        self.V_a = self._Va.copy()
        self.time = 0.0
        self.q_pump = 0.0
        self.v_pipe = 0.0
        self.acc = np.zeros(2)       # cumulative mass pumped in / returned at the flowline (kg)
        self.set_static()

    # ------------------------------------------------------------ geometry
    def _edges(self, depth):
        """Evenly spaced cell edges from surface to `depth`, about `dx` apart."""
        n = max(1, int(round(depth / self.dx)))
        return np.linspace(0.0, depth, n + 1)

    def _wall_diameter(self, md):
        """Wall diameter at `md`: casing ID where cased (if known), else the hole diameter."""
        wb = self.wellbore
        D = np.full(np.shape(md), wb.hole_diameter, dtype=float)
        if wb.casing_id is not None:
            D = np.where(wb.is_cased(md), wb.casing_id, D)
        return D

    def _string_integral(self, cum, md):
        """Integral of a string area from surface to ``md`` with the bit at self.bit_depth."""
        top = self.bit_depth - self._L_string
        s = np.clip(np.asarray(md) - top, 0.0, self._L_string)
        return np.interp(s, self._s_edges, cum) - np.interp(np.clip(-top, 0.0, self._L_string), self._s_edges, cum)

    def _geometry(self):
        """
        Cell lengths, TVDs, diameters, areas and volumes of both flow paths for the
        present bit and hole depth; the string is volume-averaged over each cell.
        """
        sv = self.wellbore.survey
        e_a, e_d = self.e_a, self.e_d
        self.L_a = np.diff(e_a)
        self.L_d = np.diff(e_d)
        c_a = 0.5 * (e_a[:-1] + e_a[1:])
        c_d = 0.5 * (e_d[:-1] + e_d[1:])
        self.z_a = sv.tvd_at(c_a)
        self.z_d = sv.tvd_at(c_d)
        self.tvd_bottom = float(sv.tvd_at([self.hole_depth])[0])
        Dw = self._wall_diameter(c_a)
        A_od = np.diff(self._string_integral(self._cum_od, e_a)) / self.L_a
        A_bore = np.diff(self._string_integral(self._cum_bore, e_d)) / self.L_d
        self.Do_a = Dw
        self.Di_a = np.sqrt(np.maximum(A_od, 0.0) * 4 / np.pi)
        self.Di_a[self.Di_a < 1e-3] = 0.0
        self.D_d = np.sqrt(np.maximum(A_bore, 1e-8) * 4 / np.pi)
        self._Va = (np.pi / 4 * Dw ** 2 - A_od) * self.L_a
        self._Vd = np.maximum(A_bore, 1e-8) * self.L_d
        self.kab = int(min(np.searchsorted(e_a, self.bit_depth, side='right') - 1, len(self.L_a) - 1))
        self.kab = max(self.kab, 0)
        self._rho_T()

    def _rho_T(self):
        """
        Temperature part of the density per cell (kg/m^3), at the reference temperature
        when no mud temperatures have been set.
        """
        fl = self.fluid
        Td = np.full(len(self.L_d), fl.reference_temperature) if self.T_d is None else self.T_d
        Ta = np.full(len(self.L_a), fl.reference_temperature) if self.T_a is None else self.T_a
        self.rT_d = (1.0 - (Td - fl.reference_temperature) * fl.thermal_expansion) * fl.density
        self.rT_a = (1.0 - (Ta - fl.reference_temperature) * fl.thermal_expansion) * fl.density

    def set_temperatures(self, md_ds, t_ds, md_ann, t_ann):
        """Mud temperatures (C) from the thermal model, interpolated to the cells."""
        self._T_src = (np.asarray(md_ds), np.asarray(t_ds), np.asarray(md_ann), np.asarray(t_ann))
        self._interp_T()

    def _interp_T(self):
        """Interpolate the stored mud-temperature profiles onto the present cells."""
        if getattr(self, '_T_src', None) is None:
            return
        md_ds, t_ds, md_ann, t_ann = self._T_src
        self.T_d = np.interp(0.5 * (self.e_d[:-1] + self.e_d[1:]), md_ds, t_ds)
        self.T_a = np.interp(0.5 * (self.e_a[:-1] + self.e_a[1:]), md_ann, t_ann)
        self._rho_T()

    # --------------------------------------------------- topology changes
    def _regrid(self, bit_new, hole_new):
        """Split/merge the moving bottom cells (on the current state) before advancing."""
        dx = self.dx
        # annulus: the bottom cell grows with the hole
        if hole_new > self.e_a[-1]:
            L_new = hole_new - self.e_a[-2]
            if L_new > self.p.cut_factor * dx:
                mid = self.e_a[-2] + 0.5 * L_new
                if mid < self.e_a[-1]:
                    f = (mid - self.e_a[-2]) / (self.e_a[-1] - self.e_a[-2])
                    m = self.m_a[-1]
                    self.m_a = np.concatenate([self.m_a[:-1], [f * m, (1 - f) * m]])
                    self.V_a = np.concatenate([self.V_a[:-1], [f * self.V_a[-1], (1 - f) * self.V_a[-1]]])
                    self.G_a = np.concatenate([self.G_a, [0.5 * self.G_a[-1]]])
                    self.e_a = np.concatenate([self.e_a[:-1], [mid, self.e_a[-1]]])
                    if self.T_a is not None:
                        self.T_a = np.append(self.T_a, self.T_a[-1])
        # drillstring: the bottom cell follows the bit
        L_new = bit_new - self.e_d[-2]
        if L_new > self.p.cut_factor * dx:
            mid = self.e_d[-2] + 0.5 * L_new
            mid = min(mid, self.e_d[-1] - 1e-6)
            f = (mid - self.e_d[-2]) / (self.e_d[-1] - self.e_d[-2])
            m = self.m_d[-1]
            self.m_d = np.concatenate([self.m_d[:-1], [f * m, (1 - f) * m]])
            self.V_d = np.concatenate([self.V_d[:-1], [f * self.V_d[-1], (1 - f) * self.V_d[-1]]])
            self.G_d = np.concatenate([self.G_d, [0.5 * (self.G_d[-1] + self.Gb[0])]])
            self.e_d = np.concatenate([self.e_d[:-1], [mid, self.e_d[-1]]])
            if self.T_d is not None:
                self.T_d = np.append(self.T_d, self.T_d[-1])
        elif L_new < 0.5 * dx and len(self.m_d) > 1:
            self.m_d = np.concatenate([self.m_d[:-2], [self.m_d[-2] + self.m_d[-1]]])
            self.V_d = np.concatenate([self.V_d[:-2], [self.V_d[-2] + self.V_d[-1]]])
            self.G_d = self.G_d[:-1].copy()
            self.e_d = np.concatenate([self.e_d[:-2], [self.e_d[-1]]])
            if self.T_d is not None:
                self.T_d = self.T_d[:-1].copy()

    # ------------------------------------------------------------- states
    def _pressure_profiles(self):
        """Gauge pressure per cell (Pa) of the drillstring and annulus, from mass and volume."""
        c2 = self.c ** 2
        p_d = c2 * (self.m_d / self.V_d - self.rT_d)
        p_a = c2 * (self.m_a / self.V_a - self.rT_a)
        return p_d, p_a

    def _rho_of_p(self, p, rT):
        """Density at pressure `p` for temperature density `rT`."""
        return rT + p / self.c ** 2

    def set_static(self):
        """Hydrostatic state with no flow (exact discrete equilibrium)."""
        self.set_steady(0.0)

    def set_steady(self, Q):
        """
        Exact discrete steady state of the transient scheme at pump rate ``Q``
        (m^3/s) with the pipe at rest: constant mass flux on every face, and on
        each face drive = friction (fixed-point on the face densities).
        """
        self._geometry()
        rho_in = self.rT_d[0]
        Gm = rho_in * Q
        pv, ty, iters = self.pv, self.ty, 50
        n_a, n_d = len(self.L_a), len(self.L_d)
        p_a = np.zeros(n_a)
        # annulus from the outlet down to the bit: return flux -Gm (upward)
        A_a = self._Va / self.L_a
        A_d = self._Vd / self.L_d
        # outlet face: p_back - p0 + rho0 g z0 = dp_fric(half cell 0)
        p = self.p.back_pressure + self.rT_a[0] * G * self.z_a[0]
        for _ in range(iters):
            rho = self._rho_of_p(p, self.rT_a[0])
            fr = 0.5 * self.L_a[0] * gradient(-Gm / (rho * A_a[0]), rho, pv, ty, self.Do_a[0], self.Di_a[0])
            p = self.p.back_pressure + rho * G * self.z_a[0] - fr
        p_a[0] = p
        for k in range(1, n_a):
            Gk = -Gm if k <= self.kab else 0.0
            pk = p_a[k - 1] + self.rT_a[k] * G * (self.z_a[k] - self.z_a[k - 1])
            for _ in range(iters):
                ru, rd = self._rho_of_p(p_a[k - 1], self.rT_a[k - 1]), self._rho_of_p(pk, self.rT_a[k])
                fr = 0.5 * self.L_a[k - 1] * gradient(Gk / (ru * A_a[k - 1]), ru, pv, ty, self.Do_a[k - 1], self.Di_a[k - 1]) \
                    + 0.5 * self.L_a[k] * gradient(Gk / (rd * A_a[k]), rd, pv, ty, self.Do_a[k], self.Di_a[k])
                pk = p_a[k - 1] + 0.5 * (ru + rd) * G * (self.z_a[k] - self.z_a[k - 1]) - fr
            p_a[k] = pk
        # bit face: p_d,last - p_a,kab + rho g dz = nozzle loss
        p_d = np.zeros(n_d)
        jd = n_d - 1
        kab = self.kab
        pk = p_a[kab]
        for _ in range(iters):
            rd_, ra_ = self._rho_of_p(pk, self.rT_d[jd]), self._rho_of_p(p_a[kab], self.rT_a[kab])
            rf = 0.5 * (rd_ + ra_)
            fr = 0.5 * self.L_d[jd] * gradient(Gm / (rd_ * A_d[jd]), rd_, pv, ty, self.D_d[jd], 0.0) \
                + 0.5 * self.L_a[kab] * gradient(Gm / (ra_ * A_a[kab]), ra_, pv, ty, self.Do_a[kab], self.Di_a[kab])
            pk = p_a[kab] - rf * G * (self.z_a[kab] - self.z_d[jd]) + bit_pressure_drop(Gm / rf, rf, self.tfa, self.cd) + fr \
                + self.dp_motor
        p_d[jd] = pk
        # drillstring upward from the bit
        for k in range(jd, 0, -1):
            pu = p_d[k] - self.rT_d[k] * G * (self.z_d[k] - self.z_d[k - 1])
            for _ in range(iters):
                ru, rd = self._rho_of_p(pu, self.rT_d[k - 1]), self._rho_of_p(p_d[k], self.rT_d[k])
                fr = 0.5 * self.L_d[k - 1] * gradient(Gm / (ru * A_d[k - 1]), ru, pv, ty, self.D_d[k - 1], 0.0) \
                    + 0.5 * self.L_d[k] * gradient(Gm / (rd * A_d[k]), rd, pv, ty, self.D_d[k], 0.0)
                pu = p_d[k] - 0.5 * (ru + rd) * G * (self.z_d[k] - self.z_d[k - 1]) + fr
            p_d[k - 1] = pu
        self.V_d, self.V_a = self._Vd.copy(), self._Va.copy()
        self.m_d = self._rho_of_p(p_d, self.rT_d) * self.V_d
        self.m_a = self._rho_of_p(p_a, self.rT_a) * self.V_a
        self.G_d = np.full(n_d, Gm)
        self.G_a = np.where(np.arange(n_a) <= kab, -Gm, 0.0).astype(float)
        self.Gb = np.array([Gm])
        self.q_pump = float(Q)

    # ------------------------------------------------------------ advance
    def advance(self, duration, q_start, q_end, bit_depth, hole_depth, dp_motor=0.0):
        """
        Advance by ``duration`` s; the pump rate ramps linearly from
        ``q_start`` to ``q_end`` (m^3/s) and the bit and hole move linearly to
        the given depths (string motion displaces mud: surge / swab).  ``dp_motor``
        (Pa) is the mud motor differential over the interval, a loss in series with
        the bit nozzles.
        """
        self.dp_motor = float(dp_motor)
        if duration <= 0:
            return self.state()
        bit_depth = float(bit_depth)
        hole_depth = max(float(hole_depth), self.hole_depth)
        self._regrid(bit_depth, hole_depth)
        self.v_pipe = (bit_depth - self.bit_depth) / duration
        # geometry at the end of the interval
        self.bit_depth, self.hole_depth = bit_depth, hole_depth
        self.e_a = np.concatenate([self.e_a[:-1], [hole_depth]])
        self.e_d = np.concatenate([self.e_d[:-1], [bit_depth]])
        self._interp_T()
        self._geometry()
        V_d0, V_a0 = self.V_d, self.V_a
        if len(V_a0) != len(self._Va):
            V_a0 = self._Va.copy()
        if len(V_d0) != len(self._Vd):
            V_d0 = self._Vd.copy()
        L_min = min(self.L_a.min(), self.L_d.min())
        dt_max = self.p.cfl * L_min / self.c
        nsub = max(1, int(np.ceil(duration / dt_max)))
        dt = duration / nsub
        nozzle_R = 1.0 / (2.0 * self.cd ** 2 * self.tfa ** 2)
        _advance(nsub, dt, self.c ** 2, self.pv, self.ty, self.p.clinging_constant, self.v_pipe,
                 q_start, q_end, float(self.rT_d[0]), self.p.back_pressure, nozzle_R, self.dp_motor, self.p.iterations,
                 self.m_d, V_d0, self._Vd, self.L_d, self.z_d, self.rT_d, self.D_d, self.G_d,
                 self.m_a, V_a0, self._Va, self.L_a, self.z_a, self.rT_a, self.Do_a, self.Di_a, self.G_a,
                 self.kab, self.Gb, self.acc)
        self.V_d, self.V_a = self._Vd.copy(), self._Va.copy()
        self.q_pump = float(q_end)
        self.time += duration
        return self.state()

    # ------------------------------------------------------------- output
    def state(self):
        """Pump (standpipe) pressure, bottomhole pressure, ECD, flows (SI, gauge)."""
        p_d, p_a = self._pressure_profiles()
        r_d0 = self.m_d[0] / self.V_d[0]
        A0 = self.V_d[0] / self.L_d[0]
        u0 = self.G_d[0] / (r_d0 * A0)
        fr0 = 0.5 * self.L_d[0] * gradient(u0 - self.v_pipe, r_d0, self.pv, self.ty, self.D_d[0], 0.0)
        p_top = p_d[0] - r_d0 * G * self.z_d[0] + fr0
        q = self.G_d[0] / self.rT_d[0]      # pump volume rate at inlet conditions
        spp = p_top + surface_loss(q, self.p.surface_loss_coefficient, self.p.surface_loss_exponent)
        rb = self.m_a[-1] / self.V_a[-1]
        bhp = p_a[-1] + rb * G * (self.tvd_bottom - self.z_a[-1])
        # pressure at the bit (annulus side)
        p_bit = p_a[self.kab]
        rb_f = 0.5 * (self.m_d[-1] / self.V_d[-1] + self.m_a[self.kab] / self.V_a[self.kab])
        ra0 = self.m_a[0] / self.V_a[0]
        return dict(spp=float(spp), p_bh=float(bhp), ecd=float(bhp / (G * self.tvd_bottom)) if self.tvd_bottom > 0 else float(self.fluid.density),
                    p_bit_annulus=float(p_bit), p_bit_bore=float(p_d[-1]),
                    dp_bit=float(bit_pressure_drop(self.Gb[0] / rb_f, rb_f, self.tfa, self.cd)),
                    q_in=float(q), q_out=float(-self.G_a[0] / ra0), q_bit=float(self.Gb[0] / (0.5 * (r_d0 + ra0))),
                    v_pipe=float(self.v_pipe), surface_loss=float(spp - p_top), dp_motor=self.dp_motor,
                    # pit gain (+) / loss (-) since the start, as mud volume at inlet density
                    pit_gain=float((self.acc[0] - self.acc[1]) / self.fluid.density) * -1.0)

    def profiles(self):
        """Pressure and density per cell of both paths (DataFrame with a `path` column)."""
        import pandas as pd
        p_d, p_a = self._pressure_profiles()
        ann = pd.DataFrame(dict(md=0.5 * (self.e_a[:-1] + self.e_a[1:]), tvd=self.z_a, p=p_a,
                                rho=self.m_a / self.V_a, path='annulus'))
        ds = pd.DataFrame(dict(md=0.5 * (self.e_d[:-1] + self.e_d[1:]), tvd=self.z_d, p=p_d,
                               rho=self.m_d / self.V_d, path='drillstring'))
        return pd.concat([ds, ann], ignore_index=True)

    def mass(self):
        """Total mud mass in the well (kg), for balance checks."""
        return float(self.m_d.sum() + self.m_a.sum())
