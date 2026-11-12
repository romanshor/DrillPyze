"""
User-facing simulator for the coupled axial-torsional drillstring model.

Workflow::

    sim = Simulation(wellbore, drillstring, fluid, params)   # set up
    sim.initialize(bit_depth=3000.0, hole_depth=3005.0)        # initialise
    sim.step(5.0, omega_sp=60 * units.RPM)                     # step ...
    sim.step(1.0, v_block=0.01)                                # ... and step again
    df = sim.results()                                         # recorded data

The simulation can be paused between ``step`` calls for as long as needed,
inputs can change between calls, and calls of any length give identical
results to one long call.  This is the hook for live data streaming and for
controllers written in Python (call ``step`` at the controller sample rate
with ``motor_torque=...``).
"""

import copy
import numpy as np
import pandas as pd

from . import kernel as K
from . import controllers
from .geometry import buoyed_weight
from .parameters import ModelParameters
from .wellbore import Fluid


class Simulation:
    """
    Parameters
    ----------
    wellbore : Wellbore
    drillstring : Drillstring
    fluid : Fluid
    params : ModelParameters, optional
    name : str
    """

    def __init__(self, wellbore, drillstring, fluid=None, params=None, name='simulation', formation=None):
        self.name = name
        self.wellbore = wellbore
        self.drillstring = drillstring
        self.fluid = fluid if fluid is not None else Fluid()
        self.params = params if params is not None else ModelParameters()
        self.formation = formation
        self.thermal = None
        self.rock = None
        self._initialized = False

    # ================================================================ setup
    def initialize(self, bit_depth, hole_depth=None, dx=10.0, record_rate=100.0,
                   refresh_distance=1.0):
        """
        Discretise the string with the bit at ``bit_depth`` (m) and start from
        rest, hanging freely from the block (statically balanced).

        Parameters
        ----------
        bit_depth : float
            Bit measured depth (m).
        hole_depth : float, optional
            Hole measured depth (m).  Defaults to ``bit_depth`` (bit just
            touching bottom with zero WOB).
        dx : float
            Maximum cell length (m).  Smaller is more accurate and slower; the
            time step scales with ``dx``.
        record_rate : float
            Rate (Hz) at which scalar channels are recorded.
        refresh_distance : float
            Axial travel of the string (m) after which the cell geometry and
            static tension are re-evaluated.
        """
        p = self.params
        self.grid = g = self.drillstring.discretize(bit_depth, dx)
        n = g.n

        # material and wave properties per cell
        self.Kt = g.J * g.G
        self.mt = g.rho * g.J
        self.Zt = g.J * np.sqrt(g.G * g.rho)
        self.Ka = g.E * g.A
        self.ma = g.rho * g.A
        self.Za = g.A * np.sqrt(g.E * g.rho)
        self.dt = p.cfl * float(np.min(g.dx / g.c_axial))
        self.zeta_p = controllers.pipe_impedance(g)
        self.w_b = buoyed_weight(g, self.fluid.density)

        # state
        self.w = np.zeros(n)
        self.tau = np.zeros(n)
        self.v = np.zeros(n)
        self.qd = np.zeros(n)
        self.wbar = np.zeros(n)      # low-passed velocities for the Stokes-layer damping
        self.vbar = np.zeros(n)
        self._scratch = [np.zeros(n + 1) for _ in range(4)]

        self.sc = np.zeros(K.N_SCALARS)
        self.sc[K.S_BITD] = bit_depth
        self.sc[K.S_HOLED] = bit_depth if hole_depth is None else max(hole_depth, bit_depth)
        # hole-bottom profile for the regenerative bit model (flat at the hole depth)
        self.surf = np.full(720, self.sc[K.S_HOLED])
        self._bit_depth0 = bit_depth
        self._geom_bitd = None
        self._qs_c = None
        self.refresh_distance = refresh_distance
        self._refresh_geometry()
        self.sc[K.S_HKLD] = -self._qs_top

        # parameters
        self.pf = np.zeros(K.N_PF)
        self.pi = np.zeros(K.N_PI, dtype=np.int64)
        self._pack_parameters()
        self.pf[K.P_HKLD_REF] = -self._qs_top

        # inputs (targets, held between calls); every set point ramps over its ramp time
        self.inputs = dict(omega_sp=0.0, v_block=0.0, motor_torque=0.0, omega_td=0.0, wob_sp=0.0,
                           flow_rate=p.thermal.flow_rate, inlet_temperature=p.thermal.inlet_temperature)
        self.ramps = np.zeros((K.N_RAMPS, 4))        # (start value, start time, ramp time, target)
        # pumps start at t = 0 and ramp to the initial pump rate; inlet temperature starts at its value
        self._slow_ramps = dict(flow_rate=[0.0, 0.0, p.ramps.flow_rate, p.thermal.flow_rate],
                                inlet_temperature=[p.thermal.inlet_temperature, 0.0, 0.0,
                                                   p.thermal.inlet_temperature])

        # slow models: wellbore temperature (Wang et al. 2023) and rock strength at the bit
        self._slow_every = max(1, int(1.0 / self.dt))                         # one slow tick per <= 1 s
        self._thermal_every = max(1, int(round(p.thermal_update_period)))    # ticks between thermal updates
        self._slow_ticks = 0
        self._slow_records = []
        self._last_slow = dict(t=0.0, hole=self.sc[K.S_HOLED], wint=0.0)
        if self.formation is not None and p.thermal.enabled:
            from .thermal import WellboreThermal
            self.thermal = WellboreThermal(self.wellbore, self.drillstring, self.fluid, self.formation,
                                           p.thermal, bit_depth, self.sc[K.S_HOLED])
        # circulating-system hydraulics (ch. 7 losses + single-phase drift flux, Wang et al. 2023)
        self.hydraulics = None
        self._hyd_records = []
        self._hyd_last = dict(t=0.0, q=0.0)
        if p.hydraulics.enabled:
            from .hydraulics import WellboreHydraulics
            self.hydraulics = WellboreHydraulics(self.wellbore, g, self.fluid, p.hydraulics,
                                                 p.bit.nozzle_area, p.bit.discharge_coefficient,
                                                 bit_depth, self.sc[K.S_HOLED], bit_depth)
            self._sync_hydraulics_temperature()
            self.hydraulics.set_static()
        self._update_rock(rop=0.0, omega=0.0)

        # recording
        self.record_rate = record_rate
        self.rec_every = max(1, int(round(1.0 / (record_rate * self.dt))))
        self._rec_phase = 0
        self._records = []
        self._steps_total = 0
        self._t_requested = 0.0
        self._initialized = True
        return self

    # ----------------------------------------------------------- internals
    def _pack_parameters(self):
        """
        Write the parameters (controller, bit, friction, autodriller, damping) into the
        kernel parameter vector and arrays.
        """
        p = self.params
        f, d, b, td = p.friction, p.damping, p.bit, p.top_drive
        c = p.controller if p.controller is not None else controllers.stiff(self.zeta_p, td.inertia)
        self.controller = c
        pf = self.pf
        pf[K.P_DT] = self.dt
        pf[K.P_KT] = d.k_t
        pf[K.P_KA] = d.k_a
        pf[K.P_WC] = f.omega_c
        pf[K.P_VC] = f.v_c
        pf[K.P_ITD] = td.inertia
        pf[K.P_KP] = c.kp
        pf[K.P_KI] = c.ki
        pf[K.P_ZC] = c.impedance
        pf[K.P_THP] = c.t_hp
        pf[K.P_TLP] = c.t_lp
        pf[K.P_TEN] = c.t_en
        pf[K.P_TAUMAX] = td.max_torque
        D = self.wellbore.hole_diameter
        a = D / 2.0
        S = self._S if (b.strength_source == 'formation' and getattr(self, '_S', None)) else float(b.drilling_strength)
        if b.model not in K.BIT_MODELS:
            raise ValueError(f'unknown bit model {b.model!r}; expected one of {K.BIT_MODELS}')
        self.pi[K.I_BIT] = K.BIT_MODELS.index(b.model)
        if b.model == 'coulomb':
            r = b.radius if b.radius is not None else D / 3.0
            pf[K.P_LEVS] = b.mu_static * r
            pf[K.P_LEVK] = b.mu_kinetic * r
            pf[K.P_KROP] = 1.0 / (2.0 * np.pi * a * S)
            pf[K.P_W0] = 0.0
        else:
            lev = a / (2.0 * b.zeta)
            pf[K.P_LEVS] = lev
            pf[K.P_LEVK] = lev
            pf[K.P_KROP] = 1.0 / (2.0 * np.pi * a * S)
            pf[K.P_W0] = b.contact_wob if b.model == 'cutting' else 0.0
        pf[K.P_KCUT] = S * a
        pf[K.P_BLADES] = max(1, int(b.blades))
        pf[K.P_QS_TOP] = self._qs_top
        pf[K.P_AD_KP] = p.autodriller.kp
        pf[K.P_AD_KI] = p.autodriller.ki
        pf[K.P_AD_VMAX] = p.autodriller.v_max
        self.pi[K.I_COUPLE] = 1 if f.along_string_coupling else 0
        pf[K.P_AXC] = f.axial_static_contrast
        if f.weakening not in ('step', 'exponential'):
            raise ValueError("friction weakening must be 'step' or 'exponential'")
        pf[K.P_WEAK] = 1.0 if f.weakening == 'exponential' else 0.0
        pf[K.P_BLK_ACC] = p.drawworks.block_acceleration
        self._pack_damping()

    # ------------------------------------------------------------- damping
    def modal_frequencies(self):
        """
        First axial and torsional natural frequencies (Hz) of the string held
        at the top and free at the bit: f_1 = 1 / (4 T), T the wave travel time.
        """
        g = self.grid
        return (1.0 / (4.0 * float(np.sum(g.dx / g.c_axial))),
                1.0 / (4.0 * float(np.sum(g.dx / g.c_torsional))))

    def _pack_damping(self):
        """
        Per-cell damping coefficients (docs/derivation/README.md sec. 5.6).

        F1  Couette (steady, laminar annulus) + Stokes layer at f_ref on the
            inner and outer walls, as rates per unit cell inertia (1/s);
        M1  Kelvin-Voigt face conductances eta_v E A / dx_face, eta_v = eta / omega_ref;
        F2  Bingham yield drag tau_y pi OD per length (N/m).
        """
        d = self.params.damping
        g = self.grid
        fl = self.fluid
        n = g.n
        f_a1, f_t1 = self.modal_frequencies()
        f_a = d.reference_frequency_axial or f_a1
        f_t = d.reference_frequency_torsional or f_t1
        w_a, w_t = 2 * np.pi * f_a, 2 * np.pi * f_t
        self.damping_reference = dict(f_axial=f_a, f_torsional=f_t, f_axial_1=f_a1, f_torsional_1=f_t1)

        ro, ri = g.od / 2.0, g.id / 2.0
        k0t = np.full(n, float(d.k_t))
        k0a = np.full(n, float(d.k_a))
        k1t = np.zeros(n)
        k1a = np.zeros(n)
        mu = fl.plastic_viscosity if fl.plastic_viscosity is not None else fl.viscosity
        rho = fl.density
        if d.fluid and mu and mu > 0.0:
            # hole (or casing) radius at each cell
            md = getattr(self, 'cell_md', g.x)
            R = np.full(n, self.wellbore.hole_diameter / 2.0)
            if self.wellbore.casing_id is not None:
                R = np.where(self.wellbore.is_cased(md), self.wellbore.casing_id / 2.0, R)
            R = np.maximum(R, 1.05 * ro)
            # steady laminar Couette in the annulus
            c0a = 2.0 * np.pi * mu / np.log(R / ro)                       # N s / m^2
            c0t = 4.0 * np.pi * mu * ro ** 2 * R ** 2 / (R ** 2 - ro ** 2)  # N m s / m
            # oscillatory Stokes layer on both walls
            sa = np.sqrt(rho * mu * w_a / 2.0)
            st = np.sqrt(rho * mu * w_t / 2.0)
            c1a = sa * 2.0 * np.pi * (ro + ri)
            c1t = st * 2.0 * np.pi * (ro ** 3 + ri ** 3)
            k0a = k0a + c0a / self.ma
            k0t = k0t + c0t / self.mt
            k1a = c1a / self.ma
            k1t = c1t / self.mt
        self.k0t, self.k0a, self.k1t, self.k1a = k0t, k0a, k1t, k1a
        fc = d.high_pass_fraction
        dt = self.dt
        self.pf[K.P_HP_T] = 1.0 - np.exp(-dt * fc * w_t) if fc > 0 else 1.0
        self.pf[K.P_HP_A] = 1.0 - np.exp(-dt * fc * w_a) if fc > 0 else 1.0

        # material damping: Rayleigh pair fitted to zeta = eta/2 at f_ref and at
        # max_frequency (zeta <= eta/2 in between).  The mass-proportional part acts
        # on the oscillatory velocity (like the Stokes term), the stiffness-
        # proportional part is Kelvin-Voigt with face conductances (zero at the
        # top drive and bit faces).
        kvt = np.zeros(n + 1)
        kva = np.zeros(n + 1)
        self.rayleigh = dict(alpha_t=0.0, beta_t=0.0, alpha_a=0.0, beta_a=0.0)
        if d.loss_factor and d.loss_factor > 0.0:
            z = 0.5 * d.loss_factor
            w2 = 2 * np.pi * d.max_frequency
            for key, w1 in (('t', w_t), ('a', w_a)):
                if w2 > w1:
                    alpha, beta = 2 * z * w1 * w2 / (w1 + w2), 2 * z / (w1 + w2)
                else:
                    alpha, beta = 0.0, 2 * z / w1
                self.rayleigh['alpha_' + key] = alpha
                self.rayleigh['beta_' + key] = beta
            k1t = k1t + self.rayleigh['alpha_t']
            k1a = k1a + self.rayleigh['alpha_a']
            self.k1t, self.k1a = k1t, k1a
            if n > 1:
                hf = 0.5 * (g.dx[:-1] + g.dx[1:])
                kvt[1:-1] = self.rayleigh['beta_t'] * 0.5 * (self.Kt[:-1] + self.Kt[1:]) / hf
                kva[1:-1] = self.rayleigh['beta_a'] * 0.5 * (self.Ka[:-1] + self.Ka[1:]) / hf
        self.kvt, self.kva = kvt, kva
        self.pf[K.P_KV] = 1.0 if d.loss_factor and d.loss_factor > 0.0 else 0.0

        # Bingham yield-stress drag on the outer wall
        if d.yield_drag not in ('always', 'pumps_off', 'off'):
            raise ValueError("damping.yield_drag must be 'always', 'pumps_off' or 'off'")
        ty = fl.yield_point or 0.0
        self._fy_full = ty * np.pi * g.od if (ty > 0 and d.yield_drag != 'off') else np.zeros(n)
        self._fy_zero = np.zeros(n)

    def _yield_drag_now(self):
        """Yield-stress drag per cell to use now (zero unless enabled for the present pump state)."""
        d = self.params.damping
        if d.yield_drag == 'always':
            return self._fy_full
        if d.yield_drag == 'pumps_off' and self.flow_rate_now() <= 0.0:
            return self._fy_full
        return self._fy_zero

    def damping_summary(self, frequency=None):
        """
        Predicted modal damping ratios per component and direction, at the
        reference frequency (default) or at ``frequency`` (Hz):

            zeta_fluid     = (k0 + k1_stokes H) / (2 omega) * s,   H = 1/(1 + (f_c/f)^2)
            zeta_material  = alpha H / (2 omega) + beta omega / 2   (Rayleigh / Kelvin-Voigt)
            zeta_numerical = pi f dx (1 - nu) / c                   (first-order Godunov)

        with s = 1 at f_ref (the model's frequency-independent Stokes
        coefficient).  ``zeta_fluid_true`` scales the Stokes part by
        sqrt(f/f_ref) as the physical layer would.
        """
        g = self.grid
        d = self.params.damping
        r = self.damping_reference
        rows = []
        for c in np.unique(g.component):
            j = int(np.where(g.component == c)[0][0])
            for mode, fr, k0, k1, cw, key in (('axial', r['f_axial'], self.k0a, self.k1a, g.c_axial, 'a'),
                                               ('torsional', r['f_torsional'], self.k0t, self.k1t, g.c_torsional, 't')):
                f = frequency or fr
                w_ = 2 * np.pi * f
                H = 1.0 / (1.0 + (d.high_pass_fraction * fr / f) ** 2)
                alpha, beta = self.rayleigh['alpha_' + key], self.rayleigh['beta_' + key]
                k1s = k1[j] - alpha
                nu = self.dt * cw[j] / g.dx[j]
                rows.append(dict(component=self.drillstring.components[c].name, mode=mode, f=f,
                                 zeta_fluid=(k0[j] + k1s * H) / (2 * w_),
                                 zeta_fluid_true=(k0[j] + k1s * H * np.sqrt(f / fr)) / (2 * w_),
                                 zeta_material=alpha * H / (2 * w_) + beta * w_ / 2,
                                 zeta_numerical=np.pi * f * g.dx[j] * (1 - nu) / cw[j]))
        return pd.DataFrame(rows)

    def _set_drilling_strength(self, S):
        """Update the drilling strength S used by the bit model (Pa)."""
        self._S = float(S)
        a = self.wellbore.hole_diameter / 2.0
        self.pf[K.P_KROP] = 1.0 / (2.0 * np.pi * a * self._S)
        self.pf[K.P_KCUT] = self._S * a

    def _update_rock(self, rop, omega):
        """Rock strength at the bit from the formation, pressure and temperature (slow update)."""
        if self.formation is None:
            return
        b = self.params.bit
        hole = self.sc[K.S_HOLED]
        tvd = float(self.wellbore.survey.tvd_at([hole])[0])
        t_mud = None
        if self.thermal is not None:
            ts = self.thermal.state()
            p_bh, t_mud = ts['p_bh'], ts['t_bit']
        else:
            p_bh = self.fluid.density * 9.80665 * tvd
        if self.hydraulics is not None:
            p_bh = self.hydraulics.state()['p_bh']
        self.rock = self.formation.strength(hole, tvd, p_bh, t_mud=t_mud, rop=rop, omega=omega, blades=b.blades)
        self.rock['p_bh'] = p_bh
        self.rock['S'] = b.zeta * self.rock['ccs_eff'] / b.cutting_efficiency
        if b.strength_source == 'formation':
            self._set_drilling_strength(self.rock['S'])

    # ------------------------------------------------------------ ramps
    def _slow_value(self, name, t=None):
        """Ramped value of a slow input (pump rate, inlet temperature) at time `t`."""
        v0, t0, T, tgt = self._slow_ramps[name]
        t = self.sc[K.S_T] if t is None else t
        if T <= 0 or t >= t0 + T:
            return tgt
        return v0 + (tgt - v0) * max(0.0, (t - t0) / T)

    def flow_rate_now(self):
        """Current (ramped) pump rate, m^3/s."""
        return self._slow_value('flow_rate')

    def _ramp_value(self, k):
        """Value of kernel ramp row `k` at the present time."""
        v0, t0, T, tgt = self.ramps[k]
        t = self.sc[K.S_T]
        if T <= 0 or t >= t0 + T:
            return tgt
        return v0 + (tgt - v0) * max(0.0, (t - t0) / T)

    def _set_ramp(self, k, target, ramp_time, start=None):
        """Start a ramp of kernel input k to ``target`` from ``start`` (default: its current value)."""
        if start is None:
            if self.ramps[k, 3] == target:
                return
            start = self._ramp_value(k)
        self.ramps[k] = (start, self.sc[K.S_T], ramp_time, target)

    def _set_slow_ramp(self, name, target, ramp_time):
        """Start a ramp of a slow input to `target` from its present value."""
        if self._slow_ramps[name][3] == target:
            return
        self._slow_ramps[name] = [self._slow_value(name), self.sc[K.S_T], ramp_time, target]

    def _sync_hydraulics_temperature(self):
        """Pass the thermal model mud temperatures to the hydraulics model (density)."""
        if self.hydraulics is None or self.thermal is None:
            return
        pr = self.thermal.profiles()
        ds = pr[pr['has_pipe']]
        if len(ds) == 0:
            ds = pr
        self.hydraulics.set_temperatures(ds['md'].values, ds['t_ds'].values, pr['md'].values, pr['t_ann'].values)

    def _hydraulics_tick(self):
        """Advance the hydraulics to the present time and record its state."""
        t = self.sc[K.S_T]
        dt = t - self._hyd_last['t']
        if self.hydraulics is None or dt <= 0:
            return
        q1 = self.flow_rate_now()
        hs = self.hydraulics.advance(dt, self._hyd_last['q'], q1, self.sc[K.S_BITD], self.sc[K.S_HOLED])
        self._hyd_last = dict(t=t, q=q1)
        self._hyd_records.append(dict(time=t, **hs))

    def _slow_tick(self):
        """Called every ~1 s of simulated time at fixed step counts (pause/resume exact)."""
        self._slow_ticks += 1
        self._refresh_geometry()
        self._hydraulics_tick()
        if self._slow_ticks % self._thermal_every:
            return
        t, hole, wint = self.sc[K.S_T], self.sc[K.S_HOLED], self.sc[K.S_WINT]
        dt = t - self._last_slow['t']
        if dt <= 0:
            return
        rop = max(hole - self._last_slow['hole'], 0.0) / dt
        omega = (wint - self._last_slow['wint']) / dt
        self._last_slow = dict(t=t, hole=hole, wint=wint)
        ts = None
        if self.thermal is not None:
            ts = self.thermal.advance(dt, flow_rate=self.flow_rate_now(),
                                      inlet_temperature=self._slow_value('inlet_temperature'),
                                      bit_depth=self.sc[K.S_BITD], hole_depth=hole, rop=rop)
            self._sync_hydraulics_temperature()
        self._update_rock(rop, omega)
        rec = dict(time=t, rop=rop, omega_bit_avg=omega, flow_rate=self.flow_rate_now())
        if ts is not None:
            rec.update({k: ts[k] for k in ('t_bit', 't_ds_bit', 't_out', 't_in', 'ecd')})
        if self.hydraulics is not None:
            hs = self.hydraulics.state()
            rec.update(spp=hs['spp'], ecd=hs['ecd'])
        if self.rock is not None:
            rec.update({k: self.rock[k] for k in ('t_rock', 'p_bh', 'ucs_t', 'ccs', 'ccs_eff', 'sigma_th',
                                                    'skin_ratio', 'w_drained', 'S')})
        self._slow_records.append(rec)

    def update_parameters(self, params=None, **_):
        """
        Replace model parameters mid-simulation (e.g. friction factors from an
        estimator).  The state is untouched.
        """
        if params is not None:
            self.params = params
        self._refresh_geometry(force=True)
        self._pack_parameters()
        if self.hydraulics is not None:
            from .hydraulics import rheology
            h, b = self.hydraulics, self.params.bit
            h.p = self.params.hydraulics
            h.tfa, h.cd = b.nozzle_area, b.discharge_coefficient
            h.pv, h.ty = rheology(self.fluid)
        self._update_rock(0.0, 0.0)

    def _refresh_geometry(self, force=False):
        """Re-evaluate cell geometry, friction zones and static tension."""
        bitd = self.sc[K.S_BITD]
        if not force and self._geom_bitd is not None and abs(bitd - self._geom_bitd) < self.refresh_distance:
            return
        g = self.grid
        shift = bitd - self._bit_depth0
        top = g.x_top + shift
        edges = np.concatenate([top, [top[-1] + g.dx[-1]]])
        sv = self.wellbore.survey
        t_e = sv.tangent(np.clip(edges, 0.0, None))
        t_c = sv.tangent(np.clip(top + 0.5 * g.dx, 0.0, None))
        kappa = (t_e[1:] - t_e[:-1]) / g.dx[:, None]
        zperp = np.array([0.0, 0.0, 1.0])[None, :] - t_c[:, 2:3] * t_c
        self.kappa = kappa
        self.zperp = zperp
        self._geo_cols = tuple(np.ascontiguousarray(a) for a in
                               (kappa[:, 0], kappa[:, 1], kappa[:, 2], zperp[:, 0], zperp[:, 1], zperp[:, 2]))
        self.t_z = t_c[:, 2].copy()
        self.cell_md = top + 0.5 * g.dx

        f = self.params.friction
        cased = self.wellbore.is_cased(self.cell_md)
        self.mus = np.where(cased, f.mu_static_cased, f.mu_static_open).astype(float)
        self.muk = np.where(cased, f.mu_kinetic_cased, f.mu_kinetic_open).astype(float)

        # static compression (frictionless hanging string), well-balanced split
        load = self.w_b * self.t_z * g.dx
        qs_edge = -np.concatenate([np.cumsum(load[::-1])[::-1], [0.0]])
        qs_c_new = 0.5 * (qs_edge[:-1] + qs_edge[1:])
        if self._qs_c is not None:
            self.qd += self._qs_c - qs_c_new      # keep total q continuous
        self._qs_c = qs_c_new
        self._qs_top = qs_edge[0]
        # normal force of the static hanging string (used by the decoupled model)
        Ts = -self._qs_c
        self.fn0 = np.linalg.norm(Ts[:, None] * kappa + self.w_b[:, None] * zperp, axis=1)
        if hasattr(self, 'pf'):
            self.pf[K.P_QS_TOP] = self._qs_top
        self._geom_bitd = bitd

    def _normal_force_now(self):
        """Soft-string normal force per length (N/m) from the present tension."""
        T = -(self._qs_c + self.qd)
        f = T[:, None] * self.kappa + self.w_b[:, None] * self.zperp
        return np.linalg.norm(f, axis=1)

    # ================================================================ step
    def step(self, duration, omega_sp=None, v_block=None, motor_torque=None,
             omega_td=None, wob_sp=None, flow_rate=None, inlet_temperature=None):
        """
        Advance the simulation by ``duration`` seconds with inputs held
        constant.  Inputs not given keep their previous values.

        Top drive (choose one; the most recently given one sets the mode):
          ``omega_sp``      speed set point (rad/s) for the internal controller
          ``motor_torque``  motor torque (N m) from an external controller
          ``omega_td``      prescribed top-drive speed (rad/s), e.g. measured RPM

        Drawworks (choose one):
          ``v_block``       block velocity (m/s, positive lowering)
          ``wob_sp``        WOB set point (N) for the autodriller

        Circulation (hydraulics and wellbore temperature):
          ``flow_rate``          pump rate (m^3/s)
          ``inlet_temperature``  mud temperature entering the drillstring (C)

        Inputs without a controller (prescribed speed, block velocity, pump
        rate, inlet temperature) approach a new value linearly over a short
        ramp time (``params.ramps``; 0 = instant), starting from the current
        value -- or, when taking over from another mode, from the present
        physical value.  Controlled set points (``omega_sp``, ``wob_sp``) are
        not ramped by default: the controller shapes the response.

        Returns a dict of the latest scalar values.
        """
        if not self._initialized:
            raise RuntimeError('Call initialize() first')
        r = self.params.ramps
        ctrl, axial = self.pi[K.I_CTRL], self.pi[K.I_AXIAL]
        # a set point taken over from another mode starts from the present physical value
        if omega_sp is not None:
            self.inputs['omega_sp'] = float(omega_sp)
            self._set_ramp(K.R_WSP, float(omega_sp), r.omega_sp,
                           start=None if ctrl == 0 else self.sc[K.S_WTD])
            self.pi[K.I_CTRL] = 0
        if motor_torque is not None:
            self.inputs['motor_torque'] = float(motor_torque)
            self._set_ramp(K.R_TAUM, float(motor_torque), r.motor_torque,
                           start=None if ctrl == 1 else self.sc[K.S_TAUM])
            self.pi[K.I_CTRL] = 1
        if omega_td is not None:
            self.inputs['omega_td'] = float(omega_td)
            self._set_ramp(K.R_WTD, float(omega_td), r.omega_td,
                           start=None if ctrl == 2 else self.sc[K.S_WTD])
            self.pi[K.I_CTRL] = 2
        if v_block is not None:
            self.inputs['v_block'] = float(v_block)
            self._set_ramp(K.R_VBLK, float(v_block), r.v_block,
                           start=None if axial == 0 else self.sc[K.S_VBLK])
            self.pi[K.I_AXIAL] = 0
        if wob_sp is not None:
            self.inputs['wob_sp'] = float(wob_sp)
            self._set_ramp(K.R_WOB, float(wob_sp), r.wob_sp,
                           start=None if axial == 1 else self.pf[K.P_HKLD_REF] - self.sc[K.S_HKLD])
            self.pi[K.I_AXIAL] = 1
        if flow_rate is not None:
            self.inputs['flow_rate'] = float(flow_rate)
            self._set_slow_ramp('flow_rate', float(flow_rate), r.flow_rate)
        if inlet_temperature is not None:
            self.inputs['inlet_temperature'] = float(inlet_temperature)
            self._set_slow_ramp('inlet_temperature', float(inlet_temperature), r.inlet_temperature)

        # integer step bookkeeping avoids drift over many short calls
        self._t_requested += duration
        target = int(round(self._t_requested / self.dt))
        remaining = target - self._steps_total
        # chunks end on fixed multiples of the slow-update interval, so slow updates
        # (geometry, temperature, rock strength) happen at the same simulated times
        # however the run is split into step() calls
        K_ = self._slow_every
        while remaining > 0:
            to_boundary = K_ - (self._steps_total % K_)
            nsteps = min(remaining, to_boundary)
            self._run(nsteps)
            remaining -= nsteps
            if self._steps_total % K_ == 0:
                self._slow_tick()
        return self.latest()

    def _run(self, nsteps):
        """Run `nsteps` kernel steps and keep the recorded rows."""
        nrec = (self._rec_phase + nsteps) // self.rec_every
        rec = np.empty((nrec, K.N_REC))
        wf, tf, vf, qf = self._scratch
        i = self.inputs
        idx, self._rec_phase = K.advance(
            nsteps, self.sc, self.pf, self.pi,
            self.grid.dx, self.Kt, self.mt, self.Zt, self.Ka, self.ma, self.Za,
            self.grid.r_contact, self.w_b,
            *self._geo_cols,
            self.mus, self.muk, self._qs_c, self.fn0,
            self.w, self.tau, self.v, self.qd,
            self.ramps,
            rec, 0, self.rec_every, self._rec_phase,
            wf, tf, vf, qf, self.surf,
            self.k0t, self.k0a, self.k1t, self.k1a, self.wbar, self.vbar,
            self.kvt, self.kva, self._yield_drag_now())
        if idx:
            self._records.append(rec[:idx])
        self._steps_total += nsteps

    # ============================================================= outputs
    @property
    def time(self):
        """Simulated time (s)."""
        return self.sc[K.S_T]

    def latest(self):
        """Latest scalar values as a dict (SI units)."""
        s = self.sc
        return dict(time=s[K.S_T], omega_td=s[K.S_WTD], motor_torque=s[K.S_TAUM],
                    pipe_torque=s[K.S_TAU0], omega_bit=s[K.S_WBIT], tob=s[K.S_TOB],
                    v_block=s[K.S_VBLK], hookload=s[K.S_HKLD], v_bit=s[K.S_VBIT],
                    wob=s[K.S_WOB], bit_depth=s[K.S_BITD], hole_depth=s[K.S_HOLED],
                    omega_sp=s[K.S_WSP], flow_rate=self.flow_rate_now(),
                    **({k: v for k, v in self.hydraulics.state().items() if k in ('spp', 'p_bh', 'ecd', 'q_out')}
                       if self.hydraulics is not None else {}))

    def tare_hookload(self):
        """
        Zero the weight indicator: take the present hookload as the off-bottom
        reference, so surface WOB = reference - hookload (used by the autodriller).
        """
        self.pf[K.P_HKLD_REF] = self.sc[K.S_HKLD]
        return float(self.pf[K.P_HKLD_REF])

    @property
    def hookload_reference(self):
        """Off-bottom reference hookload (N) for surface WOB."""
        return float(self.pf[K.P_HKLD_REF])

    def hydraulics_results(self):
        """
        Hydraulics channels, one row per ~1 s slow tick: standpipe pressure
        ``spp``, bottomhole pressure ``p_bh`` (Pa gauge), ``ecd`` (kg/m^3), pump
        rate ``q_in``, return flow ``q_out``, bit flow ``q_bit`` (m^3/s), bit
        pressure drop ``dp_bit``, surface loss and the pipe velocity used for
        surge / swab.
        """
        return pd.DataFrame(self._hyd_records)

    def slow_results(self):
        """
        Channels from the slow models (one row per thermal update): mud
        temperatures (C), bottomhole pressure (Pa), ECD (kg/m3), rock
        temperature, UCS(T), CCS, effective CCS, thermal stress (Pa), drilling
        strength S (Pa), average ROP (m/s) and bit speed (rad/s).
        """
        return pd.DataFrame(self._slow_records)

    def results(self):
        """All recorded scalar channels as a DataFrame (SI units)."""
        if not self._records:
            return pd.DataFrame(columns=K.REC_CHANNELS)
        return pd.DataFrame(np.vstack(self._records), columns=K.REC_CHANNELS)

    def profiles(self):
        """Current along-string state as a DataFrame (one row per cell)."""
        return pd.DataFrame(dict(
            md=self.cell_md, omega=self.w.copy(), torque=self.tau.copy(), v=self.v.copy(),
            tension=-(self._qs_c + self.qd),
            normal_force=self._normal_force_now() if self.params.friction.along_string_coupling else self.fn0,
            od=self.grid.od, inc=np.arccos(np.clip(self.t_z, -1, 1))))

    def clear_results(self):
        """Discard recorded data (the state is kept)."""
        self._records = []

    # ======================================================== state copies
    def snapshot(self):
        """A deep copy of the full simulation (for ensembles or rewinding)."""
        return copy.deepcopy(self)

    def info(self):
        """
        Grid and numerics summary: cells, time step, cell sizes, string length,
        pipe impedance, controller name and torsional CFL number.
        """
        g = self.grid
        return dict(cells=g.n, dt=self.dt, dx_min=float(g.dx.min()), dx_max=float(g.dx.max()),
                    length=float(g.dx.sum()), zeta_p=self.zeta_p,
                    controller=self.controller.name,
                    torsional_cfl=float(np.max(self.dt * g.c_torsional / g.dx)))
