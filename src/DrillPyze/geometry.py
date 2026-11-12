"""
Per-cell wellbore geometry used by both the torque & drag model and the
dynamic simulator, so the two compute identical normal forces.

The soft-string normal force per unit length on a cell is

    F_N = | T * dt/ds + w_b * z_perp |

where ``t`` is the unit tangent, ``T`` the effective tension, ``w_b`` the
buoyed weight per length and ``z_perp`` the component of the downward unit
vector perpendicular to the pipe.  This vector form is equivalent to the
Johancsik et al. (1984) expression and is valid for any change in azimuth
(including crossing north) without special handling.
"""

from dataclasses import dataclass
import numpy as np
from numba import njit

from . import units


@dataclass
class CellGeometry:
    """Per-cell well-path geometry seen by the string at its present position."""
    t_z: np.ndarray        # cos(inclination) at cell centre
    kappa: np.ndarray      # (n, 3) curvature vector dt/ds (1/m)
    z_perp: np.ndarray     # (n, 3) downward unit vector component normal to pipe
    inc: np.ndarray        # inclination at cell centre (rad)
    azi: np.ndarray        # azimuth at cell centre (rad)
    tvd: np.ndarray        # TVD at cell centre (m)


def cell_geometry(survey, grid):
    """Compute :class:`CellGeometry` for a discretised string in a survey."""
    edges = np.concatenate([grid.x_top, [grid.x_top[-1] + grid.dx[-1]]])
    t_edges = survey.tangent(edges)
    t_c = survey.tangent(grid.x)
    kappa = (t_edges[1:] - t_edges[:-1]) / grid.dx[:, None]
    zhat = np.array([0.0, 0.0, 1.0])
    z_perp = zhat[None, :] - t_c[:, 2:3] * t_c
    inc, azi = survey.inc_azi(grid.x)
    return CellGeometry(t_z=t_c[:, 2].copy(), kappa=kappa, z_perp=z_perp,
                        inc=inc, azi=azi, tvd=survey.tvd_at(grid.x))


def buoyed_weight(grid, fluid_density, g=units.G):
    """Buoyed weight per unit length of each cell (N/m)."""
    mass = grid.rho * grid.A
    # displaced volume per length = steel mass / steel density (tool joints included)
    return mass * g * (1.0 - fluid_density / units.STEEL_DENSITY)


@njit(cache=True, fastmath=True)
def normal_force(T, w_b, kx, ky, kz, zx, zy, zz):
    """Soft-string normal force per length for one cell."""
    fx = T * kx + w_b * zx
    fy = T * ky + w_b * zy
    fz = T * kz + w_b * zz
    return np.sqrt(fx * fx + fy * fy + fz * fz)
