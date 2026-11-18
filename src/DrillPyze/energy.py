"""
Energy accounting: mean power flows between two ``Simulation.energy()`` snapshots
(docs/derivation/README.md sec. 11).

Every term is a mean power (W) over the interval.  The three balances close up to the
numerical dissipation of the scheme, which is reported as ``numerical_*``:

* rotation:  top = bit_torsion + friction_torsion + damping_torsion + d(store_torsion)/dt
* axial:     gravity + hook_pos - hook_neg = bit_axial + friction_axial + damping_axial
             + d(store_axial)/dt
* top drive: td_pos - td_neg = top + d(store_top_drive)/dt
* mud:       h_pump = h_surface + h_string + h_motor + h_nozzles + h_annulus + h_other,
             and h_motor = motor (to the bit) + motor_loss

``pump_drive`` (h_pump / pump efficiency) and ``recoverable`` ((td_neg + hook_neg) x
regeneration efficiency) use the case's EnergyParameters.
"""

from .kernel import EN_NAMES
from .simulation import HYD_ENERGY

CUMULATIVE = EN_NAMES + tuple(f'h_{k}' for k in HYD_ENERGY)
STORES = ('store_torsion', 'store_axial', 'store_top_drive')


def flows(e0, e1, energy_params=None):
    """
    Mean power flows (W) between two energy snapshots ``e0`` and ``e1`` (dicts from
    ``Simulation.energy()``, ``e1`` later).  Returns a dict with every cumulative term,
    the storage rates ``d_store_*``, derived terms (``rock`` at the bit, ``motor_loss``,
    ``pump_drive``, ``pump_loss``, ``braking``, ``recoverable``) and the closure residuals
    ``numerical_torsion``, ``numerical_axial`` and ``numerical_top_drive``.
    """
    dt = e1['time'] - e0['time']
    if dt <= 0:
        raise ValueError('the second snapshot must be later than the first')
    p = {k: (e1[k] - e0[k]) / dt for k in CUMULATIVE}
    for k in STORES:
        p['d_' + k] = (e1[k] - e0[k]) / dt
    eta_p = getattr(energy_params, 'pump_efficiency', 0.9) if energy_params is not None else 0.9
    eta_r = getattr(energy_params, 'regen_efficiency', 0.9) if energy_params is not None else 0.9
    p['rock'] = p['bit_torsion'] + p['bit_axial'] + p['motor']
    p['motor_loss'] = p['h_motor'] - p['motor']
    p['pump_drive'] = p['h_pump'] / eta_p
    p['pump_loss'] = p['pump_drive'] - p['h_pump']
    p['braking'] = p['td_neg'] + p['hook_neg']
    p['recoverable'] = p['braking'] * eta_r
    p['numerical_torsion'] = p['top'] - p['bit_torsion'] - p['friction_torsion'] - p['damping_torsion'] - p['d_store_torsion']
    p['numerical_axial'] = (p['gravity'] + p['hook_pos'] - p['hook_neg'] - p['bit_axial'] - p['friction_axial']
                            - p['damping_axial'] - p['d_store_axial'])
    p['numerical_top_drive'] = p['td_pos'] - p['td_neg'] - p['top'] - p['d_store_top_drive']
    p['interval'] = dt
    return p
