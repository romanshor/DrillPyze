"""
Top-drive controller factories (Aarsnes, di Meglio & Shor 2018, Sec. 3).

All three industrial controllers share one block diagram: a PI speed loop,
optionally with impedance feedback of band-pass filtered pipe torque.  The
factories compute gains from the drill-pipe characteristic impedance
``zeta_p = J_p sqrt(G rho)`` and the top-drive inertia.

Controllers that do not fit this structure can be written in Python and run in
the simulator's external-torque mode (see :meth:`Simulation.step`).
"""

import numpy as np

from .parameters import ControllerParameters


def pipe_impedance(grid):
    """Torsional characteristic impedance of the top cell, zeta_p = J sqrt(G rho) (N m s)."""
    return float(grid.J[0] * np.sqrt(grid.G[0] * grid.rho[0]))


def stiff(zeta_p, inertia):
    """Industry-standard stiff PI: Kp = 100 zeta_p, Ki = 5 I_TD (eq. 8)."""
    return ControllerParameters(kp=100.0 * zeta_p, ki=5.0 * inertia, impedance=0.0, name='stiff')


def soft_torque(zeta_p, inertia, f_c):
    """
    SoftTorque/SoftSpeed: Kp = 4 zeta_p, Ki = (2 pi f_c)^2 I_TD (eqs. 11-12),
    where ``f_c`` (Hz) is the frequency at which the reflection coefficient is
    minimised, normally the first torsional mode of the string.
    """
    return ControllerParameters(kp=4.0 * zeta_p, ki=(2 * np.pi * f_c) ** 2 * inertia,
                                impedance=0.0, name='soft_torque')


def z_torque(zeta_p, inertia, t_hp=10.0, t_lp=0.001, t_en=0.001, kp=None, ki=None):
    """
    ZTorque: a stiff inner speed loop whose reference is corrected by the
    band-pass filtered pipe torque, Z = 1/zeta_p (eqs. 14-15), so the top drive
    behaves like a matched impedance.  The inner loop must be stiff for the top
    drive to actually follow the corrected reference; pairing the impedance
    feedback with soft PI gains is unstable.  Defaults follow the paper (1 ms
    speed and low-pass filters, 10 s high-pass).  Convergence to the set point
    is governed by ``t_hp``.
    """
    c = stiff(zeta_p, inertia)
    if kp is not None:
        c.kp = kp
    if ki is not None:
        c.ki = ki
    c.impedance = 1.0 / zeta_p
    c.t_hp, c.t_lp, c.t_en = t_hp, t_lp, t_en
    c.name = 'z_torque'
    return c


def custom(kp, ki, impedance=0.0, t_hp=10.0, t_lp=0.005, t_en=0.001, name='custom'):
    """Any controller with the same block structure."""
    return ControllerParameters(kp=kp, ki=ki, impedance=impedance, t_hp=t_hp,
                                t_lp=t_lp, t_en=t_en, name=name)


def reflection_coefficient(c: ControllerParameters, zeta_p, inertia, freqs):
    """
    Top-side torsional reflection coefficient ``|(C(s)+I s - zeta)/(C(s)+I s + zeta)|``
    at frequencies ``freqs`` (Hz), ignoring the impedance-feedback path (eq. 10).
    """
    s = 2j * np.pi * np.asarray(freqs)
    C = c.kp + c.ki / s + inertia * s
    return np.abs((C - zeta_p) / (C + zeta_p))
