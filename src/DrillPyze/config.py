"""
Case configuration: a plain, JSON-serialisable description of a simulation
case in SI units, and a builder that turns it into a :class:`Simulation`.

The same structure is used by the user interface, scripts and (later) the
data-driven agent, so a case can be saved, shared and reproduced::

    from DrillPyze import config
    cfg = config.well_a()                 # the example case
    cfg['parameters']['friction']['mu_static_open'] = 0.5
    sim = config.build(cfg)               # initialised Simulation
    config.save(cfg, 'case.json')
"""

import copy
import json
from dataclasses import asdict

import numpy as np

from . import units, controllers
from .survey import Survey
from .drillstring import Component, Drillstring
from .wellbore import Wellbore, Fluid
from .parameters import (ModelParameters, FrictionParameters, DampingParameters, BitParameters,
                         TopDriveParameters, AutodrillerParameters)
from .simulation import Simulation

CONTROLLER_TYPES = ('stiff', 'soft_torque', 'z_torque', 'custom')


# --------------------------------------------------------------- example case
def well_a():
    """
    Synthetic Well A of Aarsnes & Shor (2018): vertical to 1500 m, build at
    3 deg/30 m to 60 deg, then tangent; 230 m of 5-3/4" collars below 5" drill
    pipe; bit off bottom at 2500 m.  Friction gives off-bottom stick-slip with a
    stiff top-drive controller.
    """
    return {
        'name': 'Well A (Aarsnes & Shor 2018)',
        'well': {
            'survey': [
                {'md': 0.0, 'inc': 0.0, 'azi': 0.0},
                {'md': 1500.0, 'inc': 0.0, 'azi': 0.0},
                {'md': 2100.0, 'inc': 60.0, 'azi': 0.0},
                {'md': 3000.0, 'inc': 60.0, 'azi': 0.0},
            ],
            'casing_shoe': 1000.0,
            'hole_diameter': 8.5 * units.IN,
        },
        'fluid': {'density': 1200.0},
        'drillstring': [                       # bit first; last row is stretched to surface
            {'name': 'Collars', 'length': 230.0, 'od': 5.75 * units.IN, 'id': 2.8125 * units.IN,
             'linear_mass': None, 'contact_od': None},
            {'name': 'DP', 'length': 100.0, 'od': 5.0 * units.IN, 'id': 4.276 * units.IN,
             'linear_mass': None, 'contact_od': None},
        ],
        'parameters': {
            'friction': asdict(FrictionParameters(mu_static_cased=0.40, mu_kinetic_cased=0.22,
                                                  mu_static_open=0.45, mu_kinetic_open=0.25)),
            'damping': asdict(DampingParameters(k_t=0.1, k_a=0.1)),
            'bit': asdict(BitParameters()),
            'top_drive': asdict(TopDriveParameters()),
            'controller': {'type': 'stiff', 'f_c': None, 't_hp': 10.0, 't_lp': 0.001,
                           't_en': 0.001, 'kp': None, 'ki': None, 'impedance': None},
            'autodriller': asdict(AutodrillerParameters()),
            'cfl': 0.95,
        },
        'simulation': {'bit_depth': 2500.0, 'hole_depth': 2505.0, 'dx': 10.0, 'record_rate': 50.0},
    }


def default():
    """The default case (currently Well A)."""
    return well_a()


# ------------------------------------------------------------------ builder
def _num(x, default=None):
    if x is None or x == '':
        return default
    return float(x)


def build_well(cfg):
    w = cfg['well']
    s = w['survey']
    survey = Survey([_num(r['md']) for r in s], [_num(r['inc']) for r in s],
                    [_num(r.get('azi', 0.0), 0.0) for r in s])
    return Wellbore(survey, _num(w['casing_shoe'], 0.0), _num(w['hole_diameter']))


def build_drillstring(cfg):
    comps = []
    for r in cfg['drillstring']:
        comps.append(Component(name=str(r.get('name', '')), length=_num(r['length']),
                               od=_num(r['od']), id=_num(r['id']),
                               linear_mass=_num(r.get('linear_mass')),
                               contact_od=_num(r.get('contact_od'))))
    return Drillstring(comps)


LEGACY_BIT_KEYS = ('stick_torque', 'k_rop', 'wob_threshold', 'omega_b')


def bit_parameters(b, hole_diameter=None):
    """BitParameters from a config entry, translating pre-0.1 keys."""
    b = dict(b)
    k_rop = b.pop('k_rop', None)
    w0 = b.pop('wob_threshold', None)
    for k in LEGACY_BIT_KEYS:
        b.pop(k, None)
    if 'drilling_strength' not in b and k_rop and hole_diameter:
        b['drilling_strength'] = 1.0 / (2 * np.pi * (hole_diameter / 2) * k_rop)
    if w0 and 'contact_wob' not in b:
        b['contact_wob'] = w0
    if b.get('blades') is not None:
        b['blades'] = int(b['blades'])
    return BitParameters(**b)


def build_parameters(cfg, zeta_p=None, first_mode=None):
    p = cfg['parameters']
    params = ModelParameters(
        friction=FrictionParameters(**p['friction']),
        damping=DampingParameters(**p['damping']),
        bit=bit_parameters(p['bit'], _num(cfg['well'].get('hole_diameter'))),
        top_drive=TopDriveParameters(**p['top_drive']),
        autodriller=AutodrillerParameters(**p.get('autodriller', {})),
        cfl=float(p.get('cfl', 0.95)),
    )
    c = p.get('controller', {'type': 'stiff'})
    if zeta_p is not None:
        params.controller = make_controller(c, zeta_p, params.top_drive.inertia, first_mode)
    return params


def make_controller(c, zeta_p, inertia, first_mode=None):
    """Controller parameters from a controller config entry."""
    kind = c.get('type', 'stiff')
    f_c = _num(c.get('f_c')) or first_mode or 0.2
    if kind == 'stiff':
        ctrl = controllers.stiff(zeta_p, inertia)
    elif kind == 'soft_torque':
        ctrl = controllers.soft_torque(zeta_p, inertia, f_c)
    elif kind == 'z_torque':
        ctrl = controllers.z_torque(zeta_p, inertia, t_hp=_num(c.get('t_hp'), 10.0),
                                    t_lp=_num(c.get('t_lp'), 0.001), t_en=_num(c.get('t_en'), 0.001))
    elif kind == 'custom':
        ctrl = controllers.custom(_num(c.get('kp'), 100 * zeta_p), _num(c.get('ki'), 5 * inertia),
                                  impedance=_num(c.get('impedance'), 0.0),
                                  t_hp=_num(c.get('t_hp'), 10.0), t_lp=_num(c.get('t_lp'), 0.005),
                                  t_en=_num(c.get('t_en'), 0.001))
    else:
        raise ValueError(f'unknown controller type {kind!r}; expected one of {CONTROLLER_TYPES}')
    # explicit gain overrides
    for k in ('kp', 'ki', 'impedance'):
        if kind != 'custom' and c.get(k) not in (None, ''):
            setattr(ctrl, k, float(c[k]))
    return ctrl


def first_torsional_mode(sim):
    """Estimate of the first torsional mode (Hz): 1 / (4 * one-way travel time)."""
    g = sim.grid
    return float(1.0 / (4.0 * np.sum(g.dx / g.c_torsional)))


def build(cfg, initialize=True):
    """Build (and by default initialise) a :class:`Simulation` from a config dict."""
    well = build_well(cfg)
    string = build_drillstring(cfg)
    fluid = Fluid(density=_num(cfg['fluid']['density'], 1000.0))
    params = build_parameters(cfg)
    sim = Simulation(well, string, fluid, params, name=cfg.get('name', 'simulation'))
    if not initialize:
        return sim
    s = cfg['simulation']
    sim.initialize(_num(s['bit_depth']), hole_depth=_num(s.get('hole_depth')),
                   dx=_num(s.get('dx'), 10.0), record_rate=_num(s.get('record_rate'), 50.0))
    sim.params.controller = make_controller(cfg['parameters'].get('controller', {'type': 'stiff'}),
                                            sim.zeta_p, sim.params.top_drive.inertia,
                                            first_torsional_mode(sim))
    sim.update_parameters()
    return sim


def validate(cfg):
    """Return a list of human-readable problems with a config (empty if fine)."""
    errs = []
    try:
        s = cfg['well']['survey']
        md = [float(r['md']) for r in s]
        if len(md) < 2:
            errs.append('The survey needs at least two stations.')
        if any(b <= a for a, b in zip(md, md[1:])):
            errs.append('Survey measured depths must increase.')
        bit = float(cfg['simulation']['bit_depth'])
        bha = sum(float(r['length']) for r in cfg['drillstring'][:-1])
        if bit <= bha:
            errs.append(f'Bit depth ({bit:.0f} m) must exceed the BHA length ({bha:.0f} m).')
        hole = cfg['simulation'].get('hole_depth')
        if hole not in (None, '') and float(hole) < bit:
            errs.append('Hole depth must be at or below the bit depth.')
        for r in cfg['drillstring']:
            if float(r['id']) >= float(r['od']):
                errs.append(f"Component {r.get('name', '')!r}: ID must be smaller than OD.")
        if float(cfg['simulation'].get('dx', 10)) <= 0:
            errs.append('Cell size must be positive.')
        b = cfg['parameters'].get('bit', {})
        from .kernel import BIT_MODELS
        if b.get('model', 'coulomb') not in BIT_MODELS:
            errs.append(f"Bit model must be one of {', '.join(BIT_MODELS)}.")
        if b.get('model') == 'coulomb' and float(b.get('mu_kinetic', 0)) > float(b.get('mu_static', 1)):
            errs.append('Bit kinetic friction exceeds static friction.')
        if float(b.get('drilling_strength', 1) or 0) <= 0:
            errs.append('Bit drilling strength must be positive.')
        if b.get('model') in ('cutting', 'regenerative') and float(b.get('zeta', 1) or 0) <= 0:
            errs.append('Bit ζ must be positive.')
        if b.get('model') == 'regenerative' and int(b.get('blades', 1) or 0) < 1:
            errs.append('The regenerative bit needs at least one blade.')
        f = cfg['parameters']['friction']
        for loc in ('cased', 'open'):
            if float(f[f'mu_kinetic_{loc}']) > float(f[f'mu_static_{loc}']):
                errs.append(f'{loc.capitalize()}-hole kinetic friction exceeds static friction.')
    except (KeyError, TypeError, ValueError) as e:
        errs.append(f'Malformed configuration: {e}')
    return errs


# ----------------------------------------------------------------------- io
def save(cfg, path):
    with open(path, 'w') as f:
        json.dump(cfg, f, indent=2)


def load(path):
    with open(path) as f:
        return json.load(f)


def merge(base, override):
    """Deep-merge ``override`` into a copy of ``base``."""
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out
