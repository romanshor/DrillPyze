"""Case configuration and the UI server API."""
import json
import threading
import urllib.request

import numpy as np
import pytest

from DrillPyze import config, energy, units
from DrillPyze.ui import server


def test_default_case_is_json_and_valid(tmp_path):
    cfg = config.default()
    assert config.validate(cfg) == []
    p = tmp_path / 'case.json'
    config.save(cfg, p)
    assert config.load(p) == json.loads(json.dumps(cfg))


@pytest.mark.parametrize('kind', config.CONTROLLER_TYPES)
def test_build_each_controller(kind):
    cfg = config.merge(config.default(), {'parameters': {'controller': {'type': kind}}})
    sim = config.build(cfg)
    assert sim.controller.name == kind
    sim.step(0.5, omega_sp=30 * units.RPM)
    assert np.isfinite(sim.w).all()


def test_validation_catches_bad_input():
    cfg = config.default()
    cfg['simulation']['bit_depth'] = 100.0
    cfg['drillstring'][0]['id'] = 1.0
    cfg['parameters']['friction']['mu_kinetic_open'] = 0.9
    errs = ' '.join(config.validate(cfg))
    assert 'BHA length' in errs and 'ID must be smaller' in errs and 'kinetic friction exceeds' in errs


def test_start_time_validation():
    cfg = config.default()
    for ok in (None, '', '2021-01-07 06:30:00', '2021-01-07T06:30:00', '2021-01-07 06:30'):
        cfg['simulation']['start_time'] = ok
        assert not [e for e in config.validate(cfg) if 'Start date' in e], ok
    cfg['simulation']['start_time'] = '7 Jan 2021'
    assert any('Start date' in e for e in config.validate(cfg))


def test_survey_positions_match_tvd():
    s = config.build_well(config.default()).survey
    n, e, z = s.positions([1500.0, 2100.0])
    R = 600.0 / np.radians(60)
    assert z[1] == pytest.approx(1500 + R * np.sin(np.radians(60)), rel=1e-4)
    assert n[1] == pytest.approx(R * (1 - np.cos(np.radians(60))), rel=1e-3)


@pytest.fixture(scope='module')
def ui_url():
    srv = server.make_server('127.0.0.1', 0)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    yield f'http://127.0.0.1:{srv.server_address[1]}'
    srv.shutdown()


def _post(url, body):
    req = urllib.request.Request(url, data=json.dumps(body).encode(), headers={'Content-Type': 'application/json'})
    return json.loads(urllib.request.urlopen(req).read())


def test_ui_server_init_step_and_csv(ui_url):
    d = json.loads(urllib.request.urlopen(ui_url + '/api/defaults').read())
    init = _post(ui_url + '/api/init', {'config': d['config']})
    sid = init['session']
    assert init['info']['cells'] > 0 and len(init['reference']['md']) == init['info']['cells'] + 1
    fr = _post(ui_url + '/api/step', {'session': sid, 'duration': 6.0,
                                      'inputs': {'td_mode': 'setpoint', 'rpm': 60, 'dw_mode': 'velocity', 'v_block_mh': 0}})
    assert fr['time'] == pytest.approx(6.0, abs=0.01)
    assert len(fr['profile']['md']) == init['info']['cells']
    assert 0 < len(fr['trace']['time']) <= server.MAX_TRACE_ROWS
    assert fr['latest']['rpm_td'] > 30
    # block position is integrated from the block velocity at the full record rate
    fr2 = _post(ui_url + '/api/step', {'session': sid, 'duration': 3.0,
                                       'inputs': {'td_mode': 'setpoint', 'rpm': 60, 'dw_mode': 'velocity', 'v_block_mh': 360}})
    bp = fr2['trace']['block_pos']
    assert bp[0] == pytest.approx(25.0, abs=0.05)                   # 0.1 m/s, ramped over 2 s: ~0.2 m in 3 s
    assert 25.0 - bp[-1] == pytest.approx(0.1 * (3.0 - 1.0), rel=0.05)
    csv = urllib.request.urlopen(f'{ui_url}/api/results.csv?session={sid}').read().decode()
    assert csv.splitlines()[0].startswith('time,omega_td')
    html = urllib.request.urlopen(ui_url + '/').read().decode()
    assert 'DrillPyze' in html


def test_ui_server_drillers_display_fields(ui_url):
    """Channels for the driller's display, and zeroing the weight indicator with a zero-length step."""
    cfg = config.default()
    assert set(cfg['display']) >= {'rpm_limit', 'wob_limit', 'torque_limit', 'rop_limit', 'spp_limit',
                                   'pore_density', 'frac_density'}
    sid = _post(ui_url + '/api/init', {'config': cfg})['session']
    inp = {'td_mode': 'setpoint', 'rpm': 60, 'dw_mode': 'velocity', 'v_block_mh': 0}
    L = _post(ui_url + '/api/step', {'session': sid, 'duration': 3.0, 'inputs': inp})['latest']
    for k in ('wob_surf_kn', 'hk_ref_kn', 'block_mh', 'p_bh_hyd', 'tvd_bottom', 'mud_density'):
        assert L[k] is not None, k
    assert L['wob_surf_kn'] == pytest.approx(L['hk_ref_kn'] - L['hookload_kn'])
    # static pressure at bottom lies below the circulating bottomhole pressure
    assert L['mud_density'] * 9.80665 * L['tvd_bottom'] < L['p_bh_hyd']
    fr = _post(ui_url + '/api/step', {'session': sid, 'duration': 0.0, 'inputs': dict(inp, tare=True)})
    assert fr['time'] == pytest.approx(3.0, abs=0.01)
    assert fr['latest']['hk_ref_kn'] == pytest.approx(fr['latest']['hookload_kn'])
    assert fr['latest']['wob_surf_kn'] == pytest.approx(0.0, abs=1e-9)
    html = urllib.request.urlopen(ui_url + '/').read().decode()
    assert 'id="deck"' in html and '/static/deck.js' in html
    js = urllib.request.urlopen(ui_url + '/static/deck.js').read().decode()
    assert 'function deckInputs' in js


def test_ui_server_energy_frame(ui_url):
    """Frames carry every energy accumulator, the efficiencies and the motor flag (energy panel)."""
    cfg = config.default()
    cfg['parameters']['motor']['enabled'] = True
    cfg['parameters']['energy'] = {'pump_efficiency': 0.85, 'regen_efficiency': 0.8}
    sid = _post(ui_url + '/api/init', {'config': cfg})['session']
    inp = {'td_mode': 'setpoint', 'rpm': 40, 'dw_mode': 'velocity', 'v_block_mh': 0}
    e = _post(ui_url + '/api/step', {'session': sid, 'duration': 4.0, 'inputs': inp})['energy']
    for k in energy.CUMULATIVE + energy.STORES:
        assert isinstance(e[k], float), k
    assert e['has_motor'] is True and e['pump_efficiency'] == 0.85 and e['regen_efficiency'] == 0.8
    assert e['h_motor'] > 0 and e['h_pump'] > e['h_motor'] and e['td_pos'] > 0
    assert e['bit_diameter'] == pytest.approx(cfg['well']['hole_diameter'])
    html = urllib.request.urlopen(ui_url + '/').read().decode()
    assert 'id="energy"' in html and '/static/energy.js' in html and 'data-p="3"' in html


def test_tare_hookload_sets_autodriller_reference():
    sim = config.build(config.default())
    sim.step(2.0, omega_sp=60 * units.RPM)
    hk = sim.latest()['hookload']
    assert sim.tare_hookload() == pytest.approx(hk)
    assert sim.hookload_reference == pytest.approx(hk)


def test_ui_server_reports_invalid_case(ui_url):
    cfg = config.default()
    cfg['simulation']['bit_depth'] = 50.0
    with pytest.raises(urllib.error.HTTPError) as e:
        _post(ui_url + '/api/init', {'config': cfg})
    assert e.value.code == 400
    assert 'BHA length' in json.loads(e.value.read())['error']


def test_legacy_bit_keys_are_translated():
    cfg = config.default()
    a = cfg['well']['hole_diameter'] / 2
    cfg['parameters']['bit'] = {'radius': None, 'mu_static': 0.6, 'mu_kinetic': 0.4, 'stick_torque': 0.0,
                                'k_rop': 1e-8, 'wob_threshold': 5e3}
    b = config.bit_parameters(cfg['parameters']['bit'], cfg['well']['hole_diameter'])
    assert b.model == 'coulomb'
    assert b.drilling_strength == pytest.approx(1 / (2 * np.pi * a * 1e-8))
    assert b.contact_wob == 5e3
    config.build(cfg).step(0.1)


@pytest.mark.parametrize('model', ['coulomb', 'cutting', 'regenerative'])
def test_each_bit_model_builds_from_config(model):
    cfg = config.merge(config.default(), {'parameters': {'bit': {'model': model}}, 'simulation': {'hole_depth': 2500.0}})
    assert config.validate(cfg) == []
    sim = config.build(cfg)
    sim.step(1.0, omega_sp=60 * units.RPM, v_block=0.01)
    assert np.isfinite(sim.w).all() and sim.latest()['wob'] >= 0


def test_version_is_consistent_and_served(ui_url):
    """One version everywhere: package, pyproject.toml and the UI (with the release date)."""
    import pathlib
    import re
    import DrillPyze
    root = pathlib.Path(__file__).resolve().parents[1]
    assert re.search(r'^version = "([^"]+)"', (root / 'pyproject.toml').read_text(), re.M).group(1) == DrillPyze.__version__
    d = json.loads(urllib.request.urlopen(ui_url + '/api/defaults').read())
    assert d['version'] == DrillPyze.__version__
    assert d['release_date'] == DrillPyze.__release_date__ and re.match(r'\d{4}-\d{2}-\d{2}$', d['release_date'])
    assert 'id="version"' in urllib.request.urlopen(ui_url + '/').read().decode()
