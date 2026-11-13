"""
Conformance of the browser engine (lite/engine, built from docs/spec/physics-spec.md)
with the Python model: Well A with and without the wellbore temperature model, one input
script covering rotation off bottom, lowering, autodriller drilling, a pump-rate
change and pick-up (spec section 13).  Needs Node.js; skipped without it.
"""
import json
import pathlib
import shutil
import subprocess

import numpy as np
import pytest

from DrillPyze import config, units

from branches import lite_check, lite_lag_allowed

ROOT = pathlib.Path(__file__).resolve().parents[1]
NODE = shutil.which('node')
pytestmark = pytest.mark.skipif(NODE is None, reason='Node.js not installed')


def _case(thermal):
    cfg = config.well_a()
    cfg['parameters']['thermal']['enabled'] = thermal
    cfg['simulation']['hole_depth'] = 2500.5         # bit reaches bottom within the script
    return cfg


def _script():
    events = {0: dict(omega_sp=100 * units.RPM, v_block=0.0), 10: dict(v_block=60 / 3600),
              25: dict(wob_sp=60e3), 45: dict(flow_rate=0.025), 60: dict(v_block=-120 / 3600)}
    return [dict(duration=1.0, inputs=events.get(t, {})) for t in range(72)]


@pytest.fixture(scope='module', params=[True, False], ids=['thermal', 'no-thermal'])
def runs(request, tmp_path_factory):
    cfg, script = _case(request.param), _script()
    sim = config.build(cfg)
    py = [sim.step(s['duration'], **s['inputs']) for s in script]
    path = tmp_path_factory.mktemp('lite') / 'scenario.json'
    path.write_text(json.dumps(dict(config=cfg, script=script)))
    out = subprocess.run([NODE, str(ROOT / 'lite' / 'tests' / 'run_scenario.mjs'), str(path)],
                         capture_output=True, text=True)
    if out.returncode != 0:
        # the browser engine could not run the case (e.g. a new case key it does not know yet)
        return sim, py, {'error': out.stderr.strip()[-2000:]}
    return sim, py, json.loads(out.stdout)


def _js(runs):
    """Browser-engine results, or an expected failure / error where it could not run."""
    js = runs[2]
    lite_check('error' not in js, f"browser engine failed: {js.get('error')}")
    return js


def test_setup_matches(runs):
    sim, js = runs[0], _js(runs)
    for k, v in sim.info().items():
        same = js['info'][k] == pytest.approx(v, rel=1e-12) if isinstance(v, float) else js['info'][k] == v
        lite_check(same, f'info[{k!r}]: Python {v!r}, Lite {js["info"][k]!r}')


@pytest.mark.parametrize('channel', ['omega_td', 'pipe_torque', 'motor_torque', 'omega_bit', 'tob', 'hookload',
                                     'wob', 'v_bit', 'bit_depth', 'hole_depth', 'spp', 'p_bh', 'ecd', 'q_out'])
def test_channels_match(runs, channel):
    py, js = runs[1], _js(runs)
    a = np.array([r[channel] for r in py])
    b = np.array([r[channel] for r in js['latest']])
    scale = max(np.abs(a).max(), 1e-9)
    # identical arithmetic up to evaluation order: agreement to round-off over the whole script
    lite_check(np.max(np.abs(a - b)) <= 1e-6 * scale, f'{channel}: max |diff| {np.max(np.abs(a - b)):.3g}')


def test_script_drills(runs):
    """The script exercises the bit on bottom (autodriller WOB) and hoisting."""
    _, py, _ = runs
    assert max(r['wob'] for r in py) > 30e3
    assert py[-1]['bit_depth'] < max(r['bit_depth'] for r in py)


def test_recorded_rows_match(runs):
    sim, js = runs[0], _js(runs)
    a = sim.results().to_numpy()
    b = np.array(js['rec'])
    lite_check(a.shape == b.shape, f'recorded rows: Python {a.shape}, Lite {b.shape}')
    lite_check(np.allclose(a, b, rtol=1e-6, atol=1e-6 * np.abs(a).max(axis=0)), 'recorded rows differ')


def test_slow_channels_match(runs):
    """Temperatures, rock strength and hydraulics from the slow models."""
    sim, js = runs[0], _js(runs)
    py = sim.slow_results()
    jsr = js['slow']
    lite_check(len(py) == len(jsr), f'slow rows: Python {len(py)}, Lite {len(jsr)}')
    for k in ('t_bit', 't_out', 'ecd', 'ccs_eff', 'S', 'p_bh'):
        if k not in py:
            continue
        a = py[k].to_numpy(dtype=float)
        b = np.array([r.get(k, np.nan) for r in jsr], dtype=float)
        lite_check(np.allclose(a, b, rtol=1e-6, atol=1e-9), f'slow channel {k}')
