"""Case configuration and the UI server API."""
import json
import threading
import urllib.request

import numpy as np
import pytest

from DrillPyze import config, units
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
