"""
Steady-state soft-string torque and drag (SI units).

Friction on each element opposes the sliding velocity of the pipe surface,
which has an axial component ``v`` (positive downhole) and a tangential
component ``omega * r``.  Rotation therefore reduces drag, the classic
velocity-split model.  Normal force uses the vector soft-string relation in
:mod:`DrillPyze.geometry` and is updated with the local tension while
integrating from the bit up (Johancsik et al. 1984; Sheppard et al. 1987).

Sign conventions
----------------
* ``tension`` is effective tension, positive in tension.  At the bit it is
  ``-WOB``.
* ``torque`` is the torque transmitted downward, positive in the rotation
  direction.  At the bit it is ``TOB``.
* ``v > 0`` is moving downhole (slack-off / RIH), ``v < 0`` is pick-up.
"""

from dataclasses import dataclass
import numpy as np
from numba import njit

from .geometry import cell_geometry, buoyed_weight, normal_force


@dataclass
class TorqueDragResult:
    md: np.ndarray            # MD of cell edges, surface first (n+1)
    tension: np.ndarray       # effective tension at cell edges (N)
    torque: np.ndarray        # torque at cell edges (N m)
    normal_force: np.ndarray  # normal force per length on each cell (N/m)

    @property
    def hookload(self):
        """Effective tension at surface (N), excluding block weight."""
        return self.tension[0]

    @property
    def surface_torque(self):
        return self.torque[0]


@njit(cache=True)
def _integrate(dx, w_b, t_z, kappa, z_perp, r_c, mu, v, omega, wob, tob):
    n = len(dx)
    T = np.empty(n + 1)
    M = np.empty(n + 1)
    FN = np.empty(n)
    T[n] = -wob
    M[n] = tob
    for j in range(n - 1, -1, -1):
        fn = normal_force(T[j + 1], w_b[j], kappa[j, 0], kappa[j, 1], kappa[j, 2],
                          z_perp[j, 0], z_perp[j, 1], z_perp[j, 2])
        vt = omega * r_c[j]
        speed = np.sqrt(v * v + vt * vt)
        if speed > 0.0:
            ca = v / speed
            ct = vt / speed
        else:
            ca = 0.0
            ct = 0.0
        f = mu[j] * fn
        # dT/ds = f*ca - w_b*t_z  ->  integrate upward
        T[j] = T[j + 1] + (w_b[j] * t_z[j] - f * ca) * dx[j]
        M[j] = M[j + 1] + f * ct * r_c[j] * dx[j]
        FN[j] = fn
    return T, M, FN


def friction_profile(wellbore, grid, mu_cased, mu_open):
    """Per-cell friction factor: ``mu_cased`` above the shoe, ``mu_open`` below."""
    return np.where(wellbore.is_cased(grid.x), mu_cased, mu_open).astype(float)


def torque_drag(wellbore, grid, fluid_density, mu_cased, mu_open=None,
                v=0.0, omega=0.0, wob=0.0, tob=0.0, geometry=None):
    """
    Soft-string torque and drag for one operating state.

    Parameters
    ----------
    wellbore : Wellbore
    grid : drillstring.Grid
        Discretised string (from ``Drillstring.discretize``).
    fluid_density : float
        Mud density (kg/m^3) for buoyancy.
    mu_cased, mu_open : float
        Friction factors inside casing and in open hole.  ``mu_open`` defaults
        to ``mu_cased``.
    v : float
        Axial pipe velocity (m/s), positive downhole.  Only its direction
        relative to ``omega * r`` matters.
    omega : float
        Rotary speed (rad/s).
    wob, tob : float
        Weight on bit (N) and torque on bit (N m).

    Common states
    -------------
    rotating off bottom: ``omega>0, v=0``;  pick-up: ``v<0, omega=0``;
    slack-off: ``v>0, omega=0``;  drilling: ``omega>0, v=ROP, wob, tob``.
    """
    if mu_open is None:
        mu_open = mu_cased
    geo = geometry if geometry is not None else cell_geometry(wellbore.survey, grid)
    w_b = buoyed_weight(grid, fluid_density)
    mu = friction_profile(wellbore, grid, mu_cased, mu_open)
    T, M, FN = _integrate(grid.dx, w_b, geo.t_z, geo.kappa, geo.z_perp, grid.r_contact,
                          mu, float(v), float(omega), float(wob), float(tob))
    md = np.concatenate([grid.x_top, [grid.x_top[-1] + grid.dx[-1]]])
    return TorqueDragResult(md=md, tension=T, torque=M, normal_force=FN)


def broomstick(wellbore, grid, fluid_density, mu_values, mode='rotating', **kwargs):
    """
    Hookload or surface torque over a list of friction factors (the classic
    broomstick plot).  ``mode`` is 'rotating', 'pickup' or 'slackoff'.
    Returns an array with one value per friction factor.
    """
    out = []
    for mu in np.atleast_1d(mu_values):
        if mode == 'rotating':
            r = torque_drag(wellbore, grid, fluid_density, mu, omega=1.0, **kwargs)
            out.append(r.surface_torque)
        elif mode == 'pickup':
            r = torque_drag(wellbore, grid, fluid_density, mu, v=-1.0, **kwargs)
            out.append(r.hookload)
        elif mode == 'slackoff':
            r = torque_drag(wellbore, grid, fluid_density, mu, v=1.0, **kwargs)
            out.append(r.hookload)
        else:
            raise ValueError(mode)
    return np.asarray(out)
