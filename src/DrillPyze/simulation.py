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

    def __init__(self, wellbore, drillstring, fluid=None, params=None, name='simulation'):
        self.name = name
        self.wellbore = wellbore
        self.drillstring = drillstring
        self.fluid = fluid if fluid is not None else Fluid()
        self.params = params if params is not None else ModelParameters()
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
        self._scratch = [np.zeros(n + 1) for _ in range(4)]

        self.sc = np.zeros(K.N_SCALARS)
        self.sc[K.S_BITD] = bit_depth
        self.sc[K.S_HOLED] = bit_depth if hole_depth is None else max(hole_depth, bit_depth)
        self._bit_depth0 = bit_depth
        self._geom_bitd = None
        self._qs_c = None
        self.refresh_distance = refresh_distance
        self._refresh_geometry()
        self.fn0 = self._normal_force_now()
        self.sc[K.S_HKLD] = -self._qs_top

        # parameters
        self.pf = np.zeros(K.N_PF)
        self.pi = np.zeros(K.N_PI, dtype=np.int64)
        self._pack_parameters()
        self.pf[K.P_HKLD_REF] = -self._qs_top

        # inputs (held between calls)
        self.inputs = dict(omega_sp=0.0, v_block=0.0, motor_torque=0.0, omega_td=0.0, wob_sp=0.0)

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
        pf[K.P_RAMP] = td.setpoint_ramp
        pf[K.P_RBIT] = b.radius if b.radius is not None else self.wellbore.hole_diameter / 2
        pf[K.P_MUBS] = b.mu_static
        pf[K.P_MUBK] = b.mu_kinetic
        pf[K.P_WB] = b.omega_b
        pf[K.P_KROP] = b.k_rop
        pf[K.P_W0] = b.wob_threshold
        pf[K.P_QS_TOP] = self._qs_top
        pf[K.P_AD_KP] = p.autodriller.kp
        pf[K.P_AD_KI] = p.autodriller.ki
        pf[K.P_AD_VMAX] = p.autodriller.v_max
        self.pi[K.I_FN] = 0 if f.dynamic_normal_force else 1

    def update_parameters(self, params=None, **_):
        """
        Replace model parameters mid-simulation (e.g. friction factors from an
        estimator).  The state is untouched.
        """
        if params is not None:
            self.params = params
        self._refresh_geometry(force=True)
        self._pack_parameters()

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
        if hasattr(self, 'pf'):
            self.pf[K.P_QS_TOP] = self._qs_top
        self._geom_bitd = bitd

    def _normal_force_now(self):
        T = -(self._qs_c + self.qd)
        f = T[:, None] * self.kappa + self.w_b[:, None] * self.zperp
        return np.linalg.norm(f, axis=1)

    # ================================================================ step
    def step(self, duration, omega_sp=None, v_block=None, motor_torque=None,
             omega_td=None, wob_sp=None):
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

        Returns a dict of the latest scalar values.
        """
        if not self._initialized:
            raise RuntimeError('Call initialize() first')
        if omega_sp is not None:
            self.inputs['omega_sp'] = float(omega_sp)
            self.pi[K.I_CTRL] = 0
        if motor_torque is not None:
            self.inputs['motor_torque'] = float(motor_torque)
            self.pi[K.I_CTRL] = 1
        if omega_td is not None:
            self.inputs['omega_td'] = float(omega_td)
            self.pi[K.I_CTRL] = 2
        if v_block is not None:
            self.inputs['v_block'] = float(v_block)
            self.pi[K.I_AXIAL] = 0
        if wob_sp is not None:
            self.inputs['wob_sp'] = float(wob_sp)
            self.pi[K.I_AXIAL] = 1

        # integer step bookkeeping avoids drift over many short calls
        self._t_requested += duration
        target = int(round(self._t_requested / self.dt))
        remaining = target - self._steps_total
        # chunk so the geometry can follow the moving string
        max_chunk = max(1, int(1.0 / self.dt))
        while remaining > 0:
            nsteps = min(remaining, max_chunk)
            self._run(nsteps)
            remaining -= nsteps
            self._refresh_geometry()
        return self.latest()

    def _run(self, nsteps):
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
            i['omega_sp'], i['v_block'], i['motor_torque'], i['omega_td'], i['wob_sp'],
            rec, 0, self.rec_every, self._rec_phase,
            wf, tf, vf, qf)
        if idx:
            self._records.append(rec[:idx])
        self._steps_total += nsteps

    # ============================================================= outputs
    @property
    def time(self):
        return self.sc[K.S_T]

    def latest(self):
        """Latest scalar values as a dict (SI units)."""
        s = self.sc
        return dict(time=s[K.S_T], omega_td=s[K.S_WTD], motor_torque=s[K.S_TAUM],
                    pipe_torque=s[K.S_TAU0], omega_bit=s[K.S_WBIT], tob=s[K.S_TOB],
                    v_block=s[K.S_VBLK], hookload=s[K.S_HKLD], v_bit=s[K.S_VBIT],
                    wob=s[K.S_WOB], bit_depth=s[K.S_BITD], hole_depth=s[K.S_HOLED],
                    omega_sp=s[K.S_WSP])

    def results(self):
        """All recorded scalar channels as a DataFrame (SI units)."""
        if not self._records:
            return pd.DataFrame(columns=K.REC_CHANNELS)
        return pd.DataFrame(np.vstack(self._records), columns=K.REC_CHANNELS)

    def profiles(self):
        """Current along-string state as a DataFrame (one row per cell)."""
        return pd.DataFrame(dict(
            md=self.cell_md, omega=self.w.copy(), torque=self.tau.copy(), v=self.v.copy(),
            tension=-(self._qs_c + self.qd), normal_force=self._normal_force_now(),
            od=self.grid.od, inc=np.arccos(np.clip(self.t_z, -1, 1))))

    def clear_results(self):
        """Discard recorded data (the state is kept)."""
        self._records = []

    # ======================================================== state copies
    def snapshot(self):
        """A deep copy of the full simulation (for ensembles or rewinding)."""
        return copy.deepcopy(self)

    def info(self):
        g = self.grid
        return dict(cells=g.n, dt=self.dt, dx_min=float(g.dx.min()), dx_max=float(g.dx.max()),
                    length=float(g.dx.sum()), zeta_p=self.zeta_p,
                    controller=self.controller.name,
                    torsional_cfl=float(np.max(self.dt * g.c_torsional / g.dx)))
