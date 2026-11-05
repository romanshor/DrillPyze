"""
DrillPyze interactive UI server.

A small standard-library HTTP server that owns :class:`Simulation` objects and
steps them on request; the browser front end (``static/``) edits the case,
drives the transport (play / pause / step) and animates the results.

Run with ``drillpyze-ui`` (or ``python -m DrillPyze.ui``) and open the printed
address.  API (JSON over HTTP):

    GET  /api/defaults            default case config (Well A)
    POST /api/init   {config}     build + initialise a simulation -> session, geometry, first frame
    POST /api/step   {session, duration, inputs}  advance and return one frame
    GET  /api/results.csv?session=...             all recorded channels (SI) as CSV
"""

import argparse
import json
import os
import threading
import time
import uuid
import webbrowser
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

import numpy as np

from .. import config as C
from .. import kernel as K
from .. import statics, units

STATIC = os.path.join(os.path.dirname(__file__), 'static')
PLOTLY_CDN = 'https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.35.2/plotly.min.js'
MAX_TRACE_ROWS = 400          # per frame, after decimation

_sessions = {}
_sessions_lock = threading.Lock()


def _r(a, nd=4):
    """Round an array for compact JSON."""
    a = np.asarray(a, dtype=float)
    out = np.round(a, nd)
    out[~np.isfinite(out)] = 0.0
    return out.tolist()


class Session:
    def __init__(self, cfg):
        errs = C.validate(cfg)
        if errs:
            raise ValueError(' '.join(errs))
        self.cfg = cfg
        self.sim = C.build(cfg)
        self.lock = threading.Lock()
        self.sent_chunks = 0
        self.sent_slow = 0
        self.wall = 0.0
        self.touched = time.time()
        self._reference()

    # ---------------------------------------------------------- reference
    def _reference(self):
        """Static reference profiles for the plots."""
        sim = self.sim
        f = sim.params.friction
        fixed = None if f.along_string_coupling else sim.fn0
        rho = sim.fluid.density
        kin = statics.torque_drag(sim.wellbore, sim.grid, rho, f.mu_kinetic_cased, f.mu_kinetic_open,
                                  omega=1.0, normal_force=fixed)
        stat = statics.torque_drag(sim.wellbore, sim.grid, rho, f.mu_static_cased, f.mu_static_open,
                                   omega=1.0, normal_force=fixed)
        hang = statics.torque_drag(sim.wellbore, sim.grid, rho, 0.0)
        md_e = kin.md
        self.reference = {
            'md': _r(md_e, 2),
            'torque_rotating_knm': _r(kin.torque / 1e3),
            'torque_breakoff_knm': _r(stat.torque / 1e3),
            'tension_hanging_kn': _r(hang.tension / 1e3),
        }
        sv = sim.wellbore.survey
        md_path = np.linspace(0, max(sv.td, sim.sc[K.S_HOLED]), 400)
        n, e, z = sv.positions(md_path)
        az = sv.azi[-1]
        self.vs_azimuth = az
        g = sim.grid
        top = len(sim.drillstring.components) - 1        # the component stretched to surface
        self.path = {'md': _r(md_path, 1), 'vs': _r(n * np.cos(az) + e * np.sin(az), 2), 'tvd': _r(z, 2),
                     'shoe_md': sim.wellbore.casing_shoe,
                     'bha': (g.component != top).astype(int).tolist(),
                     'bha_length': float(sim.drillstring.bha_length)}

    # -------------------------------------------------------------- inputs
    def apply_inputs(self, inp):
        """Translate UI inputs (display units) into Simulation.step keyword arguments."""
        kw = {}
        td_mode = inp.get('td_mode', 'setpoint')
        rpm = float(inp.get('rpm', 0.0))
        if td_mode == 'speed':
            kw['omega_td'] = rpm * units.RPM
        else:
            kw['omega_sp'] = rpm * units.RPM
        if inp.get('dw_mode', 'velocity') == 'wob':
            kw['wob_sp'] = float(inp.get('wob_kn', 0.0)) * 1e3
        else:
            kw['v_block'] = float(inp.get('v_block_mh', 0.0)) / 3600.0
        if inp.get('flow_m3s') is not None:
            kw['flow_rate'] = max(float(inp['flow_m3s']), 0.0)
        if inp.get('inlet_c') is not None:
            kw['inlet_temperature'] = float(inp['inlet_c'])
        return kw

    def step(self, duration, inputs):
        with self.lock:
            self.touched = time.time()
            t0 = time.perf_counter()
            if duration > 0:
                self.sim.step(duration, **self.apply_inputs(inputs or {}))
            self.wall += time.perf_counter() - t0
            return self.frame(time.perf_counter() - t0, duration)

    # --------------------------------------------------------------- frame
    def frame(self, wall=0.0, duration=0.0):
        sim = self.sim
        chunks = sim._records[self.sent_chunks:]
        self.sent_chunks = len(sim._records)
        if chunks:
            rows = np.vstack(chunks)
            if len(rows) > MAX_TRACE_ROWS:
                idx = np.linspace(0, len(rows) - 1, MAX_TRACE_ROWS).round().astype(int)
                rows = rows[idx]
        else:
            rows = np.zeros((0, K.N_REC))
        col = {c: rows[:, i] for i, c in enumerate(K.REC_CHANNELS)}
        trace = {
            'time': _r(col['time'], 3),
            'rpm_td': _r(units.rads_to_rpm(col['omega_td']), 2),
            'rpm_bit': _r(units.rads_to_rpm(col['omega_bit']), 2),
            'rpm_sp': _r(units.rads_to_rpm(col['omega_sp']), 2),
            'motor_torque_knm': _r(col['motor_torque'] / 1e3, 3),
            'pipe_torque_knm': _r(col['pipe_torque'] / 1e3, 3),
            'tob_knm': _r(col['tob'] / 1e3, 3),
            'hookload_kn': _r(col['hookload'] / 1e3, 2),
            'wob_kn': _r(col['wob'] / 1e3, 2),
            'rop_mh': _r(np.where(col['bit_depth'] >= col['hole_depth'] - 1e-6, np.maximum(col['v_bit'], 0.0), 0.0) * 3600.0, 2),
            'block_mh': _r(col['v_block'] * 3600.0, 2),
            'bit_depth': _r(col['bit_depth'], 3),
            'hole_depth': _r(col['hole_depth'], 3),
        }
        prof = sim.profiles()
        sv = sim.wellbore.survey
        n, e, z = sv.positions(np.clip(prof.md.values, 0, None))
        az = self.vs_azimuth
        # slow channels (temperature, rock strength) since the last frame
        slow_rows = sim._slow_records[self.sent_slow:]
        self.sent_slow = len(sim._slow_records)
        slow = {k: [round(float(r.get(k, np.nan)), 4) if r.get(k) is not None else None for r in slow_rows]
                for k in ('time', 't_bit', 't_out', 't_in', 't_rock', 'ccs_eff', 'ccs', 'S', 'ecd', 'p_bh')}
        thermal = None
        if sim.thermal is not None:
            tp = sim.thermal.profiles()
            thermal = {'md': _r(tp.md, 1), 't_ds': [round(float(v), 2) if hp else None for v, hp in zip(tp.t_ds, tp.has_pipe)],
                       't_ann': _r(tp.t_ann, 2), 't_wall': _r(tp.t_wall, 2), 't_far': _r(tp.t_far, 2)}
        rock = sim.rock or {}
        ts = sim.thermal.state() if sim.thermal is not None else {}
        L = sim.latest()
        on_bottom = L['bit_depth'] >= L['hole_depth'] - 1e-6
        rop = max(L['v_bit'], 0.0) if on_bottom else 0.0
        return {
            'time': L['time'],
            'latest': {
                'rpm_td': units.rads_to_rpm(L['omega_td']), 'rpm_bit': units.rads_to_rpm(L['omega_bit']),
                'rpm_sp': units.rads_to_rpm(L['omega_sp']),
                'motor_torque_knm': L['motor_torque'] / 1e3, 'pipe_torque_knm': L['pipe_torque'] / 1e3,
                'tob_knm': L['tob'] / 1e3, 'hookload_kn': L['hookload'] / 1e3, 'wob_kn': L['wob'] / 1e3,
                'rop_mh': rop * 3600.0, 'bit_axial_mh': L['v_bit'] * 3600.0,
                'bit_depth': L['bit_depth'], 'hole_depth': L['hole_depth'], 'on_bottom': bool(on_bottom),
                't_bit': ts.get('t_bit'), 't_out': ts.get('t_out'), 'ecd': ts.get('ecd'),
                't_rock': rock.get('t_rock'), 'ccs_eff': rock.get('ccs_eff'), 'ccs': rock.get('ccs'),
                'S': rock.get('S'), 'p_bh': rock.get('p_bh'),
            },
            'thermal': thermal,
            'slow': slow,
            'profile': {
                'md': _r(prof.md, 2),
                'rpm': _r(units.rads_to_rpm(prof.omega), 2),
                'torque_knm': _r(prof.torque / 1e3, 3),
                'tension_kn': _r(prof.tension / 1e3, 2),
                'v_axial_mms': _r(prof.v * 1e3, 3),
                'normal_knm': _r(prof.normal_force / 1e3, 4),
                'vs': _r(n * np.cos(az) + e * np.sin(az), 2),
                'tvd': _r(z, 2),
            },
            'trace': trace,
            'perf': {'wall': wall, 'sim': duration,
                     'realtime': (duration / wall) if wall > 0 else None},
        }

    def info(self):
        i = self.sim.info()
        i.update(first_mode_hz=C.first_torsional_mode(self.sim),
                 coupling=self.sim.params.friction.along_string_coupling,
                 kp=self.sim.controller.kp, ki=self.sim.controller.ki,
                 impedance=self.sim.controller.impedance)
        return i

    def csv(self):
        df = self.sim.results()
        return df.to_csv(index=False)


def _prune():
    """Drop sessions idle for more than an hour."""
    now = time.time()
    with _sessions_lock:
        for k in [k for k, s in _sessions.items() if now - s.touched > 3600]:
            del _sessions[k]


class Handler(BaseHTTPRequestHandler):
    server_version = 'DrillPyzeUI/0.1'

    def log_message(self, fmt, *args):          # quiet by default
        if getattr(self.server, 'verbose', False):
            super().log_message(fmt, *args)

    # ---------------------------------------------------------------- utils
    def _send(self, code, body, ctype='application/json', headers=None):
        if isinstance(body, (dict, list)):
            body = json.dumps(body).encode()
        elif isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        for k, v in (headers or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self):
        n = int(self.headers.get('Content-Length') or 0)
        return json.loads(self.rfile.read(n) or b'{}')

    def _static(self, name):
        path = os.path.normpath(os.path.join(STATIC, name))
        if not path.startswith(STATIC) or not os.path.isfile(path):
            return self._send(404, {'error': 'not found'})
        ctype = {'.html': 'text/html; charset=utf-8', '.js': 'application/javascript',
                 '.css': 'text/css'}.get(os.path.splitext(path)[1], 'application/octet-stream')
        with open(path, 'rb') as f:
            self._send(200, f.read(), ctype)

    # ---------------------------------------------------------------- GET
    def do_GET(self):
        u = urlparse(self.path)
        if u.path in ('/', '/index.html'):
            return self._static('index.html')
        if u.path == '/static/plotly.min.js':
            try:
                import plotly
                p = os.path.join(os.path.dirname(plotly.__file__), 'package_data', 'plotly.min.js')
                with open(p, 'rb') as f:
                    return self._send(200, f.read(), 'application/javascript')
            except Exception:
                return self._send(302, b'', headers={'Location': PLOTLY_CDN})
        if u.path.startswith('/static/'):
            return self._static(u.path[len('/static/'):])
        if u.path == '/api/defaults':
            return self._send(200, {'config': C.default(), 'controller_types': list(C.CONTROLLER_TYPES)})
        if u.path == '/api/results.csv':
            sid = parse_qs(u.query).get('session', [''])[0]
            s = _sessions.get(sid)
            if s is None:
                return self._send(404, {'error': 'unknown session'})
            with s.lock:
                body = s.csv()
            return self._send(200, body, 'text/csv',
                              {'Content-Disposition': 'attachment; filename="drillpyze_results.csv"'})
        return self._send(404, {'error': 'not found'})

    # --------------------------------------------------------------- POST
    def do_POST(self):
        u = urlparse(self.path)
        try:
            req = self._json()
            if u.path == '/api/init':
                _prune()
                s = Session(req['config'])
                sid = uuid.uuid4().hex
                with _sessions_lock:
                    old = req.get('previous')
                    if old:
                        _sessions.pop(old, None)
                    _sessions[sid] = s
                return self._send(200, {'session': sid, 'info': s.info(), 'reference': s.reference,
                                        'path': s.path, 'frame': s.frame()})
            if u.path == '/api/step':
                s = _sessions.get(req.get('session'))
                if s is None:
                    return self._send(410, {'error': 'Session expired: press Initialise.'})
                dur = max(0.0, min(float(req.get('duration', 0.1)), 3600.0))
                return self._send(200, s.step(dur, req.get('inputs')))
            if u.path == '/api/validate':
                return self._send(200, {'errors': C.validate(req['config'])})
        except (ValueError, KeyError, TypeError) as e:
            return self._send(400, {'error': str(e)})
        except Exception as e:                    # pragma: no cover - surfaced to the UI
            return self._send(500, {'error': f'{type(e).__name__}: {e}'})
        return self._send(404, {'error': 'not found'})


def make_server(host='127.0.0.1', port=8050, verbose=False):
    srv = ThreadingHTTPServer((host, port), Handler)
    srv.verbose = verbose
    return srv


def warm_up():
    """Compile the Numba kernel once so the first Initialise is quick."""
    sim = C.build(C.default())
    sim.step(0.01, omega_sp=0.0)


def main(argv=None):
    ap = argparse.ArgumentParser(description='DrillPyze interactive simulator UI')
    ap.add_argument('--host', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=8050)
    ap.add_argument('--no-browser', action='store_true', help='do not open a browser window')
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args(argv)
    print('DrillPyze UI: compiling the simulation kernel (first run only takes a few seconds)...')
    warm_up()
    srv = make_server(a.host, a.port, a.verbose)
    url = f'http://{a.host}:{srv.server_address[1]}/'
    print(f'DrillPyze UI running at {url}  (Ctrl+C to stop)')
    if not a.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()


if __name__ == '__main__':
    main()
