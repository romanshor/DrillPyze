// Web Worker hosting the browser engine: answers the same API as the Python UI server
// (GET /api/defaults, POST /api/init, POST /api/step, GET /api/results.csv) so the
// simulation runs off the page's main thread.
import { Session } from './session.js';
import { wellA } from './config.js';
import { SPEC_REVISION } from './spec.js';
import { VERSION, RELEASE, RELEASE_DATE } from './version.js';

const sessions = new Map();
let next = 1;

function handle(path, body) {
  if (path === '/api/defaults') return { config: wellA(), controller_types: ['stiff', 'soft_torque', 'z_torque', 'custom'], spec_revision: SPEC_REVISION,
                                         version: VERSION, release: RELEASE, release_date: RELEASE_DATE };
  if (path === '/api/init') {
    const s = new Session(body.config);
    if (body.previous) sessions.delete(body.previous);
    const sid = `lite-${next++}`;
    sessions.set(sid, s);
    return { session: sid, info: s.info(), reference: s.reference, path: s.path, frame: s.frame() };
  }
  if (path === '/api/step') {
    const s = sessions.get(body.session);
    if (!s) throw new Error('Session expired: press Initialise.');
    const dur = Math.max(0, Math.min(Number(body.duration ?? 0.1), 3600));
    return s.step(dur, body.inputs);
  }
  if (path.startsWith('/api/results.csv')) {
    const sid = new URLSearchParams(path.split('?')[1] || '').get('session');
    const s = sessions.get(sid);
    if (!s) throw new Error('unknown session');
    return { csv: s.csv() };
  }
  throw new Error(`not found: ${path}`);
}

self.onmessage = (e) => {
  const { id, path, body } = e.data;
  try { self.postMessage({ id, ok: true, data: handle(path, body) }); }
  catch (err) { self.postMessage({ id, ok: false, error: err.message || String(err) }); }
};
