"""
Wellbore survey with minimum-curvature geometry.

All quantities are SI: measured depth in m, angles in radians internally.
Inclination and azimuth may be supplied in degrees with ``degrees=True``
(the default, since surveys are almost always reported in degrees).
"""

import numpy as np
import pandas as pd


def _tangent(inc, azi):
    """Unit tangent (north, east, down) for inclination/azimuth in radians."""
    return np.stack([np.sin(inc) * np.cos(azi),
                     np.sin(inc) * np.sin(azi),
                     np.cos(inc)], axis=-1)


def _dogleg(t1, t2):
    """Dogleg angle between unit tangents (robust for small and large angles)."""
    d = np.clip(np.sum(t1 * t2, axis=-1), -1.0, 1.0)
    cross = np.linalg.norm(np.cross(t1, t2), axis=-1)
    return np.arctan2(cross, d)


class Survey:
    """
    A directional survey with minimum-curvature interpolation.

    Parameters
    ----------
    md : array_like
        Measured depth of each station (m), strictly increasing.
    inc, azi : array_like
        Inclination and azimuth of each station.
    degrees : bool
        If True (default) ``inc`` and ``azi`` are in degrees.
    """

    def __init__(self, md, inc, azi, degrees=True):
        md = np.asarray(md, dtype=float)
        inc = np.asarray(inc, dtype=float)
        azi = np.asarray(azi, dtype=float)
        if degrees:
            inc = np.radians(inc)
            azi = np.radians(azi)
        order = np.argsort(md)
        md, inc, azi = md[order], inc[order], azi[order]
        # drop duplicate stations
        keep = np.concatenate([[True], np.diff(md) > 1e-9])
        md, inc, azi = md[keep], inc[keep], azi[keep]
        # prepend a vertical tie-in at surface if the survey starts below 0
        if md[0] > 0:
            md = np.concatenate([[0.0], md])
            inc = np.concatenate([[0.0], inc])
            azi = np.concatenate([[azi[0]], azi])
        if len(md) < 2:
            raise ValueError('A survey needs at least two stations')

        self.md = md
        self.inc = inc
        self.azi = azi
        self.t = _tangent(inc, azi)
        self._compute_positions()

    # ------------------------------------------------------------------ io
    @classmethod
    def from_dataframe(cls, df, md='md', inc='inc', azi='azi', md_scale=1.0, degrees=True):
        """
        Build a survey from a DataFrame.  ``md_scale`` converts the depth
        column to metres (e.g. ``units.FT`` for a survey in feet).
        """
        return cls(df[md].to_numpy() * md_scale, df[inc].to_numpy(), df[azi].to_numpy(), degrees=degrees)

    @classmethod
    def from_file(cls, path, md='md', inc='inc', azi='azi', md_scale=1.0, degrees=True, **read_kwargs):
        """Load a survey from a .csv or .xlsx file."""
        if str(path).lower().endswith(('.xls', '.xlsx')):
            df = pd.read_excel(path, **read_kwargs)
        else:
            df = pd.read_csv(path, **read_kwargs)
        return cls.from_dataframe(df, md, inc, azi, md_scale, degrees)

    # ------------------------------------------------------------- geometry
    def _compute_positions(self):
        """Station positions (north, east, TVD) by minimum curvature."""
        t1, t2 = self.t[:-1], self.t[1:]
        dmd = np.diff(self.md)
        dl = _dogleg(t1, t2)
        rf = np.ones_like(dl)
        nz = dl > 1e-12
        rf[nz] = 2.0 / dl[nz] * np.tan(dl[nz] / 2.0)
        dpos = (dmd * rf / 2.0)[:, None] * (t1 + t2)
        pos = np.vstack([np.zeros(3), np.cumsum(dpos, axis=0)])
        self.north, self.east, self.tvd = pos[:, 0], pos[:, 1], pos[:, 2]
        self.dogleg = dl
        self.dls = np.zeros_like(dl)
        self.dls = dl / dmd             # rad/m

    @property
    def td(self):
        """Deepest surveyed measured depth (m)."""
        return self.md[-1]

    def tangent(self, md):
        """
        Unit tangent vectors (north, east, down) at arbitrary measured depths,
        interpolated along the minimum-curvature arcs.  Depths beyond the last
        station are extrapolated tangentially.
        """
        md = np.atleast_1d(np.asarray(md, dtype=float))
        k = np.clip(np.searchsorted(self.md, md, side='right') - 1, 0, len(self.md) - 2)
        L = self.md[k + 1] - self.md[k]
        f = np.clip((md - self.md[k]) / L, 0.0, 1.0)
        t1, t2 = self.t[k], self.t[k + 1]
        dl = self.dogleg[k]
        out = np.empty_like(t1)
        small = dl < 1e-9
        # straight segment
        out[small] = t1[small]
        # slerp along the arc
        b = dl[~small]
        ff = f[~small]
        s = np.sin(b)
        out[~small] = (np.sin((1 - ff) * b) / s)[:, None] * t1[~small] + (np.sin(ff * b) / s)[:, None] * t2[~small]
        out /= np.linalg.norm(out, axis=1, keepdims=True)
        return out

    def inc_azi(self, md):
        """Interpolated inclination and azimuth (radians) at measured depths."""
        t = self.tangent(md)
        inc = np.arccos(np.clip(t[:, 2], -1, 1))
        azi = np.mod(np.arctan2(t[:, 1], t[:, 0]), 2 * np.pi)
        return inc, azi

    def positions(self, md):
        """
        (north, east, tvd) in metres at measured depths, by integrating the
        interpolated tangent along the minimum-curvature path.
        """
        md = np.atleast_1d(np.asarray(md, dtype=float))
        need = max(self.md[-1], md.max())
        if getattr(self, '_pos_md', None) is None or self._pos_md[-1] < need:
            # cached position table (extended tangentially beyond the last station)
            grid = np.union1d(self.md, np.linspace(0, need * 1.25 + 100.0, 6000))
            t = self.tangent(grid)
            d = 0.5 * (t[1:] + t[:-1]) * np.diff(grid)[:, None]
            self._pos_md = grid
            self._pos = np.vstack([np.zeros(3), np.cumsum(d, axis=0)])
        return tuple(np.interp(md, self._pos_md, self._pos[:, k]) for k in range(3))

    def tvd_at(self, md):
        """True vertical depth (m) at measured depths (numerical arc integration)."""
        return self.positions(md)[2]
