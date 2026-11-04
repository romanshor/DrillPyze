"""
Drillstring description and spatial discretisation.

Components are listed from the bit upward (the usual BHA-report order).  The
top component (normally drill pipe) is stretched or shortened so that the
string reaches the requested bit depth.  All values are SI.
"""

from dataclasses import dataclass, field
import numpy as np
import pandas as pd

from . import units


@dataclass
class Component:
    """
    One drillstring component.

    Parameters
    ----------
    name : str
        Description, e.g. 'Bit', 'Motor', 'DC', 'HWDP', 'DP'.
    length : float
        Length (m).
    od, id : float
        Body outer and inner diameter (m).
    linear_mass : float, optional
        Mass per length including tool joints (kg/m).  If omitted, the mass of
        a plain steel tube is used.
    contact_od : float, optional
        Diameter that contacts the borehole (tool joint or stabiliser OD, m).
        Defaults to ``od``.
    E, G : float
        Young's and shear modulus (Pa).  Steel by default.
    density : float
        Material density (kg/m^3), used when ``linear_mass`` is not given.
    """
    name: str
    length: float
    od: float
    id: float
    linear_mass: float = None
    contact_od: float = None
    E: float = units.STEEL_E
    G: float = units.STEEL_G
    density: float = units.STEEL_DENSITY

    @property
    def area(self):
        return np.pi / 4 * (self.od ** 2 - self.id ** 2)

    @property
    def polar_moment(self):
        return np.pi / 32 * (self.od ** 4 - self.id ** 4)

    @property
    def mass_per_length(self):
        return self.linear_mass if self.linear_mass is not None else self.density * self.area

    @property
    def effective_density(self):
        """Density that gives the stated mass per length over the body area."""
        return self.mass_per_length / self.area


@dataclass
class Grid:
    """Cell-wise properties of a discretised drillstring (index 0 at surface)."""
    x_top: np.ndarray       # MD of the top of each cell (m)
    dx: np.ndarray          # cell length (m)
    od: np.ndarray
    id: np.ndarray
    r_contact: np.ndarray   # contact radius (m)
    A: np.ndarray           # cross-sectional area (m^2)
    J: np.ndarray           # polar moment (m^4)
    rho: np.ndarray         # effective density (kg/m^3)
    E: np.ndarray
    G: np.ndarray
    component: np.ndarray   # index into the component list

    @property
    def n(self):
        return len(self.dx)

    @property
    def x(self):
        """MD of cell centres (m)."""
        return self.x_top + 0.5 * self.dx

    @property
    def x_bottom(self):
        return self.x_top + self.dx

    @property
    def c_axial(self):
        return np.sqrt(self.E / self.rho)

    @property
    def c_torsional(self):
        return np.sqrt(self.G / self.rho)


class Drillstring:
    """
    A drillstring built from components listed bit-first.

    The last component in the list (top of the string) is the one that is
    lengthened to reach surface; normally drill pipe.
    """

    def __init__(self, components, name='Drillstring'):
        if len(components) == 0:
            raise ValueError('At least one component is required')
        self.name = name
        self.components = list(components)

    # ------------------------------------------------------------------ io
    @classmethod
    def from_dataframe(cls, df, name='Drillstring', length_scale=1.0, diameter_scale=1.0,
                       mass_scale=1.0, bit_first=True):
        """
        Build from a table with columns ``name, length, od, id`` and optional
        ``linear_mass, contact_od``.  The scale arguments convert the table's
        units to SI (e.g. ``length_scale=units.FT, diameter_scale=units.IN,
        mass_scale=units.LB_FT``).
        """
        cols = {c.lower(): c for c in df.columns}
        rows = df.iterrows() if bit_first else df.iloc[::-1].iterrows()
        comps = []
        for _, r in rows:
            def get(key, scale):
                c = cols.get(key)
                if c is None or pd.isna(r[c]):
                    return None
                return float(r[c]) * scale
            comps.append(Component(
                name=str(r[cols.get('name', cols.get('component'))]),
                length=get('length', length_scale),
                od=get('od', diameter_scale),
                id=get('id', diameter_scale),
                linear_mass=get('linear_mass', mass_scale),
                contact_od=get('contact_od', diameter_scale),
            ))
        return cls(comps, name)

    @property
    def bha_length(self):
        """Length of all components except the top (filler) component."""
        return sum(c.length for c in self.components[:-1])

    # ------------------------------------------------------- discretisation
    def discretize(self, bit_depth, dx=10.0):
        """
        Discretise the string with the bit at ``bit_depth`` (m MD).  Each
        component is split into equal cells no longer than ``dx`` so cell edges
        align with component boundaries.

        Returns
        -------
        Grid
        """
        top_len = bit_depth - self.bha_length
        if top_len <= 0:
            raise ValueError('Bit depth is shallower than the BHA length')

        # surface-down order
        comps = list(reversed(self.components))
        lengths = [top_len] + [c.length for c in comps[1:]]

        arrays = {k: [] for k in ('dx', 'od', 'id', 'rc', 'A', 'J', 'rho', 'E', 'G', 'comp')}
        for k, (c, L) in enumerate(zip(comps, lengths)):
            n = max(1, int(np.ceil(L / dx - 1e-9)))
            h = L / n
            cod = c.contact_od if c.contact_od is not None else c.od
            arrays['dx'] += [h] * n
            arrays['od'] += [c.od] * n
            arrays['id'] += [c.id] * n
            arrays['rc'] += [cod / 2] * n
            arrays['A'] += [c.area] * n
            arrays['J'] += [c.polar_moment] * n
            arrays['rho'] += [c.effective_density] * n
            arrays['E'] += [c.E] * n
            arrays['G'] += [c.G] * n
            arrays['comp'] += [len(comps) - 1 - k] * n   # index in the bit-first list

        dxs = np.asarray(arrays['dx'])
        x_top = np.concatenate([[0.0], np.cumsum(dxs)[:-1]])
        return Grid(x_top=x_top, dx=dxs,
                    od=np.asarray(arrays['od']), id=np.asarray(arrays['id']),
                    r_contact=np.asarray(arrays['rc']), A=np.asarray(arrays['A']),
                    J=np.asarray(arrays['J']), rho=np.asarray(arrays['rho']),
                    E=np.asarray(arrays['E']), G=np.asarray(arrays['G']),
                    component=np.asarray(arrays['comp'], dtype=int))
