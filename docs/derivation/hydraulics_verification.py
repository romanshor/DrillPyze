"""
Hydraulics verification figures and numbers (docs/derivation/hydraulics.md, sec. 5).

1. Textbook steady state (ch. 7 Example 7.1, Practice Problems 7.1, 7.2): model vs chapter.
2. Water hammer: Joukowsky surge and acoustic travel time on the Example 7.1 well.
3. Well A pump start-up in the full simulation: SPP, BHP / ECD and flow in / out.
4. Surge and swab while tripping on the Example 7.1 well.

Run:  python docs/derivation/hydraulics_verification.py   (writes docs/figures/hydraulics_*.png)
"""
import os
import sys
import json
import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'tests'))
from test_hydraulics import (_textbook, _example71_model, EX71, PP71, PP72, GPM, PSI, PPG, IN, FT)  # noqa: E402
from DrillPyze import config, units  # noqa: E402

FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
out = {}

# 1. textbook ------------------------------------------------------------------
printed = {'Example 7.1': dict(spp=2134, ecd=12.01, bit=1119),
           'Problem 7.1': dict(spp=1868, ecd=9.68, bit=1157),
           'Problem 7.2': dict(spp=2005, ecd=12.68, bit=1120)}
for name, data in (('Example 7.1', EX71), ('Problem 7.1', PP71), ('Problem 7.2', PP72)):
    r = _textbook(*data)
    out[name] = dict(sections_psi=[round(x, 1) for x in (r['sections'].dp / PSI)],
                     bit_psi=round(r['bit'] / PSI, 1), spp_psi=round(r['spp'] / PSI, 1),
                     ecd_ppg=round(r['ecd'] / PPG, 3), printed=printed[name])
h = _example71_model(sound_speed=1e6)
h.set_steady(400 * GPM)
st = h.state()
out['Example 7.1 model'] = dict(spp_psi=st['spp'] / PSI, bhp_psi=st['p_bh'] / PSI, ecd_ppg=st['ecd'] / PPG,
                                dp_bit_psi=st['dp_bit'] / PSI)
hc = _example71_model(sound_speed=1500.0)
hc.set_steady(400 * GPM)
out['Example 7.1 model, compressible'] = dict(spp_psi=hc.state()['spp'] / PSI, bhp_psi=hc.state()['p_bh'] / PSI,
                                              ecd_ppg=hc.state()['ecd'] / PPG)

# 2. water hammer ----------------------------------------------------------------
h = _example71_model(newtonian=True)
h.set_static()
s0 = h.state()['spp']
bore0 = h.profiles().query("path == 'drillstring'")['p'].values.copy()
Q = 100 * GPM
A = np.pi / 4 * (4.276 * IN) ** 2
rho, c = 11.5 * PPG, 1500.0
h.advance(0.02, 0.0, Q, h.bit_depth, h.hole_depth)
t, spp, pbit = [0.02], [h.state()['spp'] - s0], [0.0]
for _ in range(800):
    h.advance(0.01, Q, Q, h.bit_depth, h.hole_depth)
    t.append(t[-1] + 0.01)
    spp.append(h.state()['spp'] - s0)
    pbit.append(h.profiles().query("path == 'drillstring'")['p'].values[-2] - bore0[-2])
t, spp, pbit = map(np.array, (t, spp, pbit))
dpj = rho * c * Q / A
L_bit = h.e_d[-2] + 0.5 * h.L_d[-2]
early = (t > 0.2) & (t < 0.6)
out['water hammer'] = dict(joukowsky_MPa=dpj / 1e6, spp_early_MPa=float(np.mean(spp[early]) / 1e6),
                           arrival_s=float(t[np.argmax(pbit > 0.5 * dpj)]), L_over_c_s=L_bit / c)
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(t, spp / 1e6, label='standpipe pressure rise')
ax.plot(t, pbit / 1e6, label='bore pressure rise at the bit')
ax.axhline(dpj / 1e6, ls=':', c='gray', label='Joukowsky rho c dQ/A')
ax.axvline(L_bit / c, ls='--', c='gray', lw=0.8)
ax.axvline(2 * L_bit / c, ls='--', c='gray', lw=0.8)
ax.set_xlabel('time (s)'); ax.set_ylabel('pressure rise (MPa)'); ax.legend(fontsize=8)
ax.set_title('Water hammer: 100 gpm pump step, Example 7.1 well, 1 cP mud (dashed: L/c, 2L/c)', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hydraulics_water_hammer.png'), dpi=150)

# 3. Well A pump start-up ------------------------------------------------------
cfg = config.well_a()
s = config.build(cfg)
Q = s.inputs['flow_rate']
s.step(5.0, omega_sp=0.0, v_block=0.0)
s.step(175.0)
r = s.hydraulics_results()
ref = copy.deepcopy(s.hydraulics)
ref.set_steady(Q)
out['Well A pump start'] = dict(Q_m3s=Q, spp_end_psi=float(r.spp.iloc[-1] / PSI),
                                spp_steady_psi=ref.state()['spp'] / PSI,
                                pbh_end_psi=float(r.p_bh.iloc[-1] / PSI), pbh_steady_psi=ref.state()['p_bh'] / PSI,
                                ecd_end_ppg=float(r.ecd.iloc[-1] / PPG))
fig, axs = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
axs[0].plot(r.time, r.q_in / GPM, label='pump'); axs[0].plot(r.time, r.q_out / GPM, label='return')
axs[0].set_ylabel('flow (gpm)'); axs[0].legend(fontsize=8)
axs[1].plot(r.time, r.spp / PSI); axs[1].axhline(ref.state()['spp'] / PSI, ls=':', c='gray')
axs[1].set_ylabel('SPP (psi)')
axs[2].plot(r.time, r.ecd / PPG); axs[2].set_ylabel('ECD (ppg)'); axs[2].set_xlabel('time (s)')
axs[0].set_title('Well A: pumps ramped to %.0f gpm over %.0f s (dotted: steady state)' % (Q / GPM, s.params.ramps.flow_rate), fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hydraulics_pump_start.png'), dpi=150)

# 4. surge and swab -----------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 3.5))
for label, v in (('running in 0.5 m/s', 0.5), ('pulling out 0.5 m/s', -0.5)):
    h = _example71_model()
    h.bit_depth = h.bit_depth - 200.0
    h.e_d = h._edges(h.bit_depth)
    h.set_static()
    p0 = h.state()['p_bh']
    bit, tt, dp = h.bit_depth, [], []
    for k in range(120):
        vv = v if k < 60 else 0.0
        bit += 0.5 * vv
        h.advance(0.5, 0.0, 0.0, bit, h.hole_depth)
        tt.append(0.5 * (k + 1)); dp.append(h.state()['p_bh'] - p0)
    ax.plot(tt, np.array(dp) / PSI, label=label)
    out['surge/swab ' + label] = dict(dp_bh_psi_at_30s=float(dp[59] / PSI), emw_ppg=float(dp[59] / (9.80665 * h.tvd_bottom) / PPG))
ax.axvline(30, ls=':', c='gray')
ax.set_xlabel('time (s)'); ax.set_ylabel('BHP change (psi)'); ax.legend(fontsize=8)
ax.set_title('Surge and swab, Example 7.1 well and mud, pumps off (pipe stops at 30 s)', fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hydraulics_surge_swab.png'), dpi=150)

with open(os.path.join(FIG, 'hydraulics_results.json'), 'w') as fh:
    json.dump(out, fh, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
