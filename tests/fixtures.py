"""Shared synthetic cases (SI)."""
import numpy as np
from DrillPyze import Survey, Component, Drillstring, Wellbore, Fluid, units


def well_a_survey():
    """Aarsnes & Shor (2018) Well A: vertical to 1500 m, build 3 deg/30 m to 60 deg, hold."""
    md = [0, 1500, 2100, 3000]
    inc = [0, 0, 60, 60]
    azi = [0, 0, 0, 0]
    return Survey(md, inc, azi)


def standard_string(dp_mass=None):
    return Drillstring([
        Component('Collars', 230.0, od=5.75 * units.IN, id=2.8125 * units.IN),
        Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN, linear_mass=dp_mass),
    ])


def well_a(casing_shoe=1000.0):
    return Wellbore(well_a_survey(), casing_shoe=casing_shoe, hole_diameter=8.5 * units.IN)


def vertical_well(td=3000.0):
    return Wellbore(Survey([0, td], [0, 0], [0, 0]), casing_shoe=0.0, hole_diameter=8.5 * units.IN)


def uniform_string():
    return Drillstring([Component('DP', 100.0, od=5.0 * units.IN, id=4.276 * units.IN)])


def bandpass_envelope(y, fs, f_lo, f_hi, pad=4):
    """
    Zero-phase band-pass and Hilbert envelope with numpy FFTs only (no scipy):
    keep the [f_lo, f_hi] band (cosine-tapered edges) and return the magnitude
    of the analytic signal.  The record is zero-padded ``pad`` times so the
    circular FFT does not wrap the large early amplitude onto the tail.
    """
    y = np.asarray(y, dtype=float)
    y = y - np.mean(y)
    n0 = len(y)
    n = pad * n0
    Y = np.fft.fft(y, n)
    f = np.fft.fftfreq(n, 1.0 / fs)
    af = np.abs(f)
    bw = 0.2 * (f_hi - f_lo)
    w = np.zeros(n)
    w[(af >= f_lo) & (af <= f_hi)] = 1.0
    lo = (af > f_lo - bw) & (af < f_lo)
    hi = (af > f_hi) & (af < f_hi + bw)
    w[lo] = 0.5 * (1 - np.cos(np.pi * (af[lo] - (f_lo - bw)) / bw))
    w[hi] = 0.5 * (1 + np.cos(np.pi * (af[hi] - f_hi) / bw))
    # analytic signal: double the positive frequencies, drop the negative ones
    w = w * np.where(f > 0, 2.0, np.where(f == 0, 1.0, 0.0))
    return np.abs(np.fft.ifft(Y * w))[:n0]


def decay_ratio(t, y, f, t_start=1.0, trim=0.15):
    """Damping ratio of the mode near ``f`` (Hz) from the log-envelope slope of a free decay."""
    t = np.asarray(t)
    m = t > t_start
    fs = 1.0 / np.median(np.diff(t))
    env, tt = bandpass_envelope(np.asarray(y)[m], fs, 0.6 * f, 1.4 * f), t[m]
    span = tt[-1] - tt[0]
    k = (tt > tt[0] + trim * span) & (tt < tt[-1] - trim * span)
    return -np.polyfit(tt[k], np.log(env[k]), 1)[0] / (2 * np.pi * f)
