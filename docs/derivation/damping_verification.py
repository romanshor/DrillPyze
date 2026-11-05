"""
Damping verification figures (docs/derivation/README.md sec. 4.6):

1. free decay of the first axial and torsional modes for each damping model,
   measured vs predicted damping ratio;
2. predicted damping ratio vs frequency (fluid, material, numerical) for Well A DP;
3. Well A tripping in at 200 m/h with zero RPM: lumped k = 0.1 1/s vs F1+M1.

Run:  python docs/derivation/damping_verification.py   (writes docs/figures/damping_*.png)
"""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'tests'))
from fixtures import vertical_well, uniform_string  # noqa: E402
from DrillPyze import (Simulation, ModelParameters, FrictionParameters, DampingParameters, Fluid, config)  # noqa: E402
from scipy.signal import hilbert, butter, sosfiltfilt  # noqa: E402

FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
CASES = [('lumped k=0.1/s', DampingParameters(0.1, 0.1, fluid=False, loss_factor=0.0)),
         ('F1 fluid', DampingParameters(fluid=True, loss_factor=0.0)),
         ('M1 material (eta=0.02)', DampingParameters(fluid=False, loss_factor=0.02)),
         ('F1+M1 (eta=0.01)', DampingParameters(fluid=True, loss_factor=0.01))]


def free_decay(damping, mode):
    p = ModelParameters(friction=FrictionParameters(0, 0, 0, 0), damping=damping)
    sim = Simulation(vertical_well(4000), uniform_string(), Fluid(1200.0, viscosity=0.02), p)
    sim.initialize(3000.0, hole_depth=3500.0, dx=10.0, record_rate=200)
    if mode == 'axial':
        sim.step(0.5, v_block=0.2, omega_td=0.0); sim.step(60.0, v_block=0.0)
        ch, f = 'v_bit', sim.modal_frequencies()[0]
    else:
        sim.step(0.5, omega_td=2.0, v_block=0.0); sim.step(60.0, omega_td=0.0)
        ch, f = 'omega_bit', sim.modal_frequencies()[1]
    r = sim.results()
    t, y = r['time'].values, r[ch].values
    m = t > 1.0
    fs = 1.0 / np.median(np.diff(t))
    yb = sosfiltfilt(butter(2, [0.6 * f, 1.4 * f], btype='band', fs=fs, output='sos'), y[m])
    env, tt = np.abs(hilbert(yb)), t[m]
    span = tt[-1] - tt[0]
    k = (tt > tt[0] + 0.15 * span) & (tt < tt[-1] - 0.15 * span)
    zeta = -np.polyfit(tt[k], np.log(env[k]), 1)[0] / (2 * np.pi * f)
    row = sim.damping_summary().set_index('mode').loc[mode]
    if damping.fluid or damping.loss_factor:
        pred = row['zeta_fluid'] + row['zeta_material'] + row['zeta_numerical']
    else:
        pred = damping.k_a / (4 * np.pi * f) + row['zeta_numerical']
    return dict(t=t, y=y, f=f, zeta=zeta, pred=pred)


def main():
    results = {}
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for name, d in CASES:
        for ax, mode in zip(axs, ('axial', 'torsional')):
            r = free_decay(d, mode)
            results[f'{name} | {mode}'] = dict(f=r['f'], measured=r['zeta'], predicted=r['pred'])
            ax.plot(r['t'], r['y'], lw=0.8, label=f"{name}: zeta={r['zeta']:.4f} (pred. {r['pred']:.4f})")
    axs[0].set_ylabel('bit velocity (m/s)'); axs[1].set_ylabel('bit speed (rad/s)')
    axs[1].set_xlabel('time (s)')
    for ax, mode in zip(axs, ('axial', 'torsional')):
        ax.set_title(f'Free decay, first {mode} mode, 3000 m of 5" DP, vertical, no wall friction', fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'damping_free_decay.png'), dpi=150)

    # damping ratio vs frequency, Well A DP section, default damping (F1+M1, eta=0.01)
    sim = config.build(config.well_a())
    fr = np.logspace(np.log10(0.1), np.log10(30), 80)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, mode in zip(axs, ('axial', 'torsional')):
        rows = [sim.damping_summary(f).query("component == 'DP' and mode == @mode").iloc[0] for f in fr]
        zf = np.array([r['zeta_fluid'] for r in rows]); zt = np.array([r['zeta_fluid_true'] for r in rows])
        zm = np.array([r['zeta_material'] for r in rows]); zn = np.array([r['zeta_numerical'] for r in rows])
        ax.loglog(fr, zf, label='fluid, model (Stokes at f_ref)')
        ax.loglog(fr, zt, '--', label='fluid, Stokes sqrt(f) scaling')
        ax.loglog(fr, zm, label='material, Rayleigh fit (eta=0.01)')
        ax.loglog(fr, 0 * fr + 0.005, ':', c='gray', label='material, hysteretic eta/2')
        ax.loglog(fr, zn, label='numerical (dx = 10 m)')
        ax.loglog(fr, zf + zm + zn, 'k', lw=2, label='model total')
        ax.axvline(10, c='r', lw=0.8); ax.set_title(f'Well A drill pipe, {mode}', fontsize=9)
        ax.set_xlabel('frequency (Hz)')
    axs[0].set_ylabel('damping ratio'); axs[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'damping_vs_frequency.png'), dpi=150)
    results['well_a_reference'] = sim.damping_reference
    results['well_a_summary'] = sim.damping_summary().round(5).to_dict('records')

    # Well A tripping in at 200 m/h, zero RPM
    fig, ax = plt.subplots(figsize=(8, 3.5))
    for name, d in (('lumped k=0.1/s (previous default)', DampingParameters(0.1, 0.1, fluid=False, loss_factor=0.0)),
                    ('F1+M1 (new default)', DampingParameters())):
        cfg = config.well_a()
        cfg['parameters']['damping'] = vars(d).copy()
        cfg['simulation']['hole_depth'] = cfg['simulation']['bit_depth'] + 100.0
        s = config.build(cfg)
        s.step(120.0, omega_td=0.0, v_block=200.0 / 3600.0)
        r = s.results()
        ax.plot(r['time'], r['v_bit'] * 3600, lw=0.8, label=name)
        late = r[r['time'] > 30]
        results[f'trip {name}'] = dict(v_bit_range_after_30s_mph=[float(late['v_bit'].min() * 3600), float(late['v_bit'].max() * 3600)])
    ax.set_xlabel('time (s)'); ax.set_ylabel('bit speed (m/h)'); ax.legend(fontsize=8)
    ax.set_title('Well A, tripping in at 200 m/h block speed, 0 RPM', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'damping_trip_well_a.png'), dpi=150)
    with open(os.path.join(FIG, 'damping_results.json'), 'w') as fh:
        json.dump(results, fh, indent=1, default=float)
    print(json.dumps(results, indent=1, default=float))


if __name__ == '__main__':
    main()
