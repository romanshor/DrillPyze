"""
Conformance of the browser engine (lite/engine, built from docs/spec/physics-spec.md)
with the Python model: Well A without the wellbore temperature model, one input
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

ROOT = pathlib.Path(__file__).resolve().parents[1]
NODE = shutil.which('node')
pytestmark = pytest.mark.skipif(NODE is None, reason='Node.js not installed')


def _case():
    cfg = config.well_a()
    cfg['parameters']['thermal']['enabled'] = False
    cfg['simulation']['hole_depth'] = 2500.5         # bit reaches bottom within the script
    return cfg


def _script():
    events = {0: dict(omega_sp=100 * units.RPM, v_block=0.0), 10: dict(v_block=60 / 3600),
              25: dict(wob_sp=60e3), 45: dict(flow_rate=0.025), 60: dict(v_block=-120 / 3600)}
    return [dict(duration=1.0, inputs=events.get(t, {})) for t in range(72)]


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    cfg, script = _case(), _script()
    sim = config.build(cfg)
    py = [sim.step(s['duration'], **s['inputs']) for s in script]
    path = tmp_path_factory.mktemp('lite') / 'scenario.json'
    path.write_text(json.dumps(dict(config=cfg, script=script)))
    out = subprocess.run([NODE, str(ROOT / 'lite' / 'tests' / 'run_scenario.mjs'), str(path)],
                         capture_output=True, text=True, check=True)
    return sim, py, json.loads(out.stdout)


def test_setup_matches(runs):
    sim, _, js = runs
    for k, v in sim.info().items():
        if isinstance(v, float):
            assert js['info'][k] == pytest.approx(v, rel=1e-12), k
        else:
            assert js['info'][k] == v, k


@pytest.mark.parametrize('channel', ['omega_td', 'pipe_torque', 'motor_torque', 'omega_bit', 'tob', 'hookload',
                                     'wob', 'v_bit', 'bit_depth', 'hole_depth', 'spp', 'p_bh', 'ecd', 'q_out'])
def test_channels_match(runs, channel):
    _, py, js = runs
    a = np.array([r[channel] for r in py])
    b = np.array([r[channel] for r in js['latest']])
    scale = max(np.abs(a).max(), 1e-9)
    # identical arithmetic up to evaluation order: agreement to round-off over the whole script
    assert np.max(np.abs(a - b)) <= 1e-6 * scale, f'{channel}: max |diff| {np.max(np.abs(a - b)):.3g}'


def test_script_drills(runs):
    """The script exercises the bit on bottom (autodriller WOB) and hoisting."""
    _, py, _ = runs
    assert max(r['wob'] for r in py) > 30e3
    assert py[-1]['bit_depth'] < max(r['bit_depth'] for r in py)


def test_recorded_rows_match(runs):
    sim, _, js = runs
    a = sim.results().to_numpy()
    b = np.array(js['rec'])
    assert a.shape == b.shape
    assert np.allclose(a, b, rtol=1e-6, atol=1e-6 * np.abs(a).max(axis=0))
