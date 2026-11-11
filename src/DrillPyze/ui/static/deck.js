/* DrillPyze driller's display ("deck").
 *
 * Arc gauges for WOB, rotary speed, torque, ROP, standpipe and bottomhole
 * pressure, and bars for pump rate and mud weight / ECD.  Each gauge keeps the
 * surface measurement as its needle; the model's downhole estimate is an
 * orange arrow outside the dial with its range over the last 10 s as a band.
 * Set points are cyan diamonds dragged on the dials themselves (or focused and
 * moved with the arrow keys); limits are red markers and zones.
 *
 * Loaded before app.js and uses its globals (state, toDisp, fromDisp, unitOf,
 * srv, SYSTEM, $, get, advance) only from inside functions.
 */
'use strict';

// live control state, SI (rpm in rpm); app.js sends it with every step
const ctl = { tdMode: 'setpoint', rpm: 60, dwMode: 'velocity', wob: 0, vblk: 0, flow: null, tin: null, tare: false };

const DISPLAY_DEFAULTS = { rpm_limit: 220, wob_limit: 150e3, torque_limit: 25e3, rop_limit: 90 / 3600,
                           spp_limit: 4000 * 6894.757293168, pore_density: 1100, frac_density: 1500 };
const DK = { surf: '#f4f7fa', dh: '#f2a33a', sp: '#3fc6e0', lim: '#e5534b', limMark: '#ff8a83', tick: '#7f93a6',
             track: '#26333f', dim: '#93a5b6', off: '#7e8d9c', bg: '#0f1720', ref: '#c3ccd6' };
const G_ACC = 9.80665;

const deck = { key: null, nodes: [], bars: [], drag: null, toastT: null, last: null };

function displayLimit(k) {
  const v = state.cfg && get(state.cfg, `display.${k}`);
  return v == null || !isFinite(v) ? DISPLAY_DEFAULTS[k] : Number(v);
}

// "nice" scale in display units: about `n` major divisions
function niceScale(lo, hi, n = 5) {
  const span = Math.max(hi - lo, 1e-9), raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 3, 5, 10].map(m => m * mag).find(s => s >= raw * 0.999);
  return { min: Math.floor(lo / step + 1e-9) * step, max: Math.ceil(hi / step - 1e-9) * step, step };
}
const decFor = (step) => (step >= 10 ? 0 : step >= 1 ? 1 : 2);
const fx = (v, d) => (v == null || !isFinite(v) ? '–' : Number(v).toFixed(d));
const tickLabel = (v, step) => {
  if (Math.abs(v) < 1e-9) return '0';
  if (Math.abs(v) >= 1000 && step >= 500) return `${+(v / 1000).toFixed(2)}k`;
  return String(+v.toFixed(decFor(step) + 1));
};

// ---------------------------------------------------------------- geometry
const DW = 220, DH = 196, CX = 110, CY = 106, R = 80, A0 = -135, SW = 270;
const BW = 100, BH = 176, BX = 42, BT = 10, BB = 166;
const dAng = (d, v) => A0 + SW * Math.max(0, Math.min(1, (v - d.min) / (d.max - d.min)));
const dVal = (d, a) => d.min + (d.max - d.min) * Math.max(0, Math.min(1, (a - A0) / SW));
const dPt = (a, r) => [CX + r * Math.sin(a * Math.PI / 180), CY - r * Math.cos(a * Math.PI / 180)];
const xy = (p) => `${p[0].toFixed(1)} ${p[1].toFixed(1)}`;
function dArc(a0, a1, r) {
  if (a1 < a0 + 0.3) a1 = a0 + 0.3;
  return `M${xy(dPt(a0, r))} A${r} ${r} 0 ${a1 - a0 > 180 ? 1 : 0} 1 ${xy(dPt(a1, r))}`;
}
function radial(a, r0, r1, col, w) {
  const p = dPt(a, r0), q = dPt(a, r1);
  return `<line x1="${p[0].toFixed(1)}" y1="${p[1].toFixed(1)}" x2="${q[0].toFixed(1)}" y2="${q[1].toFixed(1)}" stroke="${col}" stroke-width="${w}"/>`;
}
const bY = (b, v) => BB - (BB - BT) * Math.max(0, Math.min(1, (v - b.min) / (b.max - b.min)));
const bV = (b, y) => b.min + (b.max - b.min) * Math.max(0, Math.min(1, (BB - y) / (BB - BT)));

// ----------------------------------------------------------- definitions
// Values are read from a frame's `latest` (server units) and returned in display units.
function deckDefs(L) {
  const tvd = L && L.tvd_bottom, rho = (L && L.mud_density) || get(state.cfg, 'fluid.density') || 1200;
  const pd = (dens) => (tvd ? toDisp('pressure', dens * G_ACC * tvd) : null);
  const lim = { rpm: displayLimit('rpm_limit'), wob: toDisp('force', displayLimit('wob_limit')),
                tq: toDisp('torque', displayLimit('torque_limit')), rop: toDisp('rop', displayLimit('rop_limit')),
                spp: toDisp('pressure', displayLimit('spp_limit')) };
  // pump scale from the case's pump rate only: a scale that followed the set point would
  // grow under a drag towards the top and rebuild the bar mid-drag (the handle vanished)
  const flow0 = get(state.cfg, 'parameters.thermal.flow_rate') || 0.03;
  const pore = pd(displayLimit('pore_density')), frac = pd(displayLimit('frac_density'));
  const sc = (lo, hi, n) => niceScale(lo, hi, n);
  const dials = [
    Object.assign(sc(0, lim.wob * 1.3, 4), { id: 'wob', title: 'WOB', kind: 'force', lim: lim.wob,
      surf: (l) => srv('force', l.wob_surf_kn), dh: (l) => srv('force', l.wob_kn), trace: 'wob_kn', traceKind: 'force',
      chip: () => 'Autodriller', chipOn: () => ctl.dwMode === 'wob', chipClick: () => setDwMode('wob'),
      chip2: 'Zero', chip2Title: 'Zero the weight indicator (take the present hookload as off-bottom)', chip2Click: tare,
      sp: () => toDisp('force', ctl.wob), spActive: () => ctl.dwMode === 'wob', spSet: (v) => { ctl.wob = Math.max(0, fromDisp('force', v)); },
      spMsg: () => `Autodriller WOB set point → ${fx(toDisp('force', ctl.wob), 1)} ${unitOf('force')}` }),
    Object.assign(sc(0, lim.rpm * 1.12, 5), { id: 'rpm', title: 'Rotary', kind: 'rpm', lim: lim.rpm,
      surf: (l) => l.rpm_td, dh: (l) => l.rpm_bit, trace: 'rpm_bit', traceKind: 'rpm',
      chip: () => (ctl.tdMode === 'setpoint' ? 'Controller' : 'Prescribed'), chipOn: () => true,
      chipTitle: 'Top drive: speed controller set point, or prescribed top-drive speed (click to switch)',
      chipClick: () => { ctl.tdMode = ctl.tdMode === 'setpoint' ? 'speed' : 'setpoint';
        deckToast(ctl.tdMode === 'setpoint' ? 'Top drive: speed controller' : 'Top drive: prescribed speed'); },
      sp: () => ctl.rpm, spActive: () => true, spSet: (v) => { ctl.rpm = Math.max(0, v); },
      spMsg: () => `Top drive ${ctl.tdMode === 'setpoint' ? 'set point' : 'speed'} → ${fx(ctl.rpm, 0)} rpm` }),
    Object.assign(sc(0, lim.tq * 1.2, 6), { id: 'tq', title: 'Torque', kind: 'torque', lim: lim.tq,
      surf: (l) => srv('torque', l.pipe_torque_knm), dh: (l) => srv('torque', l.tob_knm), trace: 'tob_knm', traceKind: 'torque' }),
    Object.assign(sc(-lim.rop * 0.67, lim.rop * 1.33, 6), { id: 'rop', title: 'ROP', kind: 'rop', lim: lim.rop, zero: true,
      surf: (l) => srv('rop', l.block_mh), dh: (l) => srv('rop', l.rop_mh), trace: 'rop_mh', traceKind: 'rop',
      chip: () => 'Block speed', chipOn: () => ctl.dwMode === 'velocity', chipClick: () => setDwMode('velocity'),
      sp: () => toDisp('rop', ctl.vblk), spActive: () => ctl.dwMode === 'velocity', spSet: (v) => { ctl.vblk = fromDisp('rop', v); },
      spMsg: () => `Block speed → ${fx(toDisp('rop', ctl.vblk), 0)} ${unitOf('rop')}${ctl.vblk < 0 ? ' (hoisting)' : ''}` }),
    Object.assign(sc(0, lim.spp * 1.25, 5), { id: 'spp', title: 'Standpipe', kind: 'pressure', lim: lim.spp,
      surf: (l) => toDisp('pressure', l.spp), surfLabel: 'SPP' }),
  ];
  const bhpLo = pore != null ? pore * 0.85 : 0, bhpHi = frac != null ? frac * 1.08 : lim.spp;
  dials.push(Object.assign(sc(bhpLo, bhpHi, 4), { id: 'bhp', title: 'Bottomhole', kind: 'pressure', lim: frac, limLo: pore,
    needleDh: (l) => toDisp('pressure', l.p_bh_hyd), ref: () => pd(rho), refLabel: 'Static', dhLabel: 'BHP' }));
  const mwLo = Math.min(displayLimit('pore_density'), rho), mwHi = Math.max(displayLimit('frac_density'), rho);
  const bars = [
    Object.assign(sc(0, toDisp('flow', flow0 * 2), 5), { id: 'q', title: 'Pump', kind: 'flow',
      sp: () => toDisp('flow', ctl.flow), spSet: (v) => { ctl.flow = Math.max(0, fromDisp('flow', v)); },
      spMsg: () => `Pump rate → ${fx(toDisp('flow', ctl.flow), 0)} ${unitOf('flow')}` }),
    Object.assign(sc(toDisp('ecd', mwLo * 0.92), toDisp('ecd', mwHi * 1.04), 6), { id: 'mw', title: 'Mud wt', kind: 'ecd',
      pore: toDisp('ecd', displayLimit('pore_density')), frac: toDisp('ecd', displayLimit('frac_density')), mw: toDisp('ecd', rho) }),
  ];
  // decimals from the scale size: 0 for hundreds and up, 1 below 50, 2 below 5
  for (const d of dials.concat(bars)) d.dec = Math.max(Math.abs(d.min), Math.abs(d.max)) > 50 ? 0 : Math.max(Math.abs(d.min), Math.abs(d.max)) > 5 ? 1 : 2;
  return { dials, bars };
}

// ------------------------------------------------------------- dial drawing
function dialFace(d) {
  let g = `<path d="${dArc(A0, A0 + SW, R)}" stroke="${DK.track}" stroke-width="11" fill="none"/>`;
  if (d.lim != null && d.lim < d.max) g += `<path d="${dArc(dAng(d, d.lim), A0 + SW, R)}" stroke="${DK.lim}" stroke-width="11" fill="none" opacity=".7"/>`;
  if (d.limLo != null && d.limLo > d.min) g += `<path d="${dArc(A0, dAng(d, d.limLo), R)}" stroke="${DK.lim}" stroke-width="11" fill="none" opacity=".7"/>`;
  for (const lv of [d.lim, d.limLo]) if (lv != null && lv > d.min && lv < d.max) g += radial(dAng(d, lv), R - 9, R + 9, DK.limMark, 3);
  const n = Math.round((d.max - d.min) / d.step);
  for (let k = 0; k <= 2 * n; k++) {
    const v = d.min + k * d.step / 2, a = dAng(d, v), major = k % 2 === 0;
    g += radial(a, R - 7, R - (major ? 17 : 12), DK.tick, major ? 1.8 : 1);
    if (major) { const t = dPt(a, R - 28); g += `<text x="${t[0].toFixed(1)}" y="${(t[1] + 4).toFixed(1)}" text-anchor="middle" font-size="11.5" fill="${DK.dim}">${tickLabel(v, d.step)}</text>`; }
  }
  if (d.zero && d.min < 0) g += radial(dAng(d, 0), R - 7, R + 7, DK.ref, 2.5);
  g += `<text x="${CX}" y="${CY + 77}" text-anchor="middle" font-size="11.5" fill="${DK.dim}">${unitOf(d.kind)}</text>`;
  return g;
}

function diamondAt(a, active) {
  const p = dPt(a, R);
  return `<g transform="translate(${xy(p)}) rotate(${a.toFixed(1)})"><circle r="15" fill="transparent"/>
    <path class="vis" d="M0 -10 L7.5 0 L0 10 L-7.5 0 Z" fill="${active ? DK.sp : DK.bg}" stroke="${active ? DK.bg : DK.off}" stroke-width="${active ? 1.8 : 2}"/></g>`;
}

// range (display units) of a trace channel over the last 10 s before time t
function traceRange(key, kind, t) {
  const T = state.traces; if (!T || !T.time.length || !T[key]) return null;
  let hi = T.time.length - 1;
  if (t != null) { let lo = 0; while (lo < hi) { const m = (lo + hi + 1) >> 1; if (T.time[m] <= t + 1e-9) lo = m; else hi = m - 1; } hi = lo; }
  let mn = Infinity, mx = -Infinity, s = 0, n = 0;
  for (let k = hi; k >= 0 && T.time[k] >= T.time[hi] - 10; k--) { const v = T[key][k]; if (v == null) continue; mn = Math.min(mn, v); mx = Math.max(mx, v); s += v; n++; }
  if (!n) return null;
  const c = (v) => (kind === 'rpm' ? v : srv(kind, v));
  return [c(mn), c(mx), c(s / n)];
}

function dialDynamic(d, L, t) {
  let g = '';
  const vs = d.surf ? d.surf(L) : null, vd = d.dh ? d.dh(L) : null, vn = d.needleDh ? d.needleDh(L) : null;
  if (d.trace) { const r = traceRange(d.trace, d.traceKind, t);
    if (r) g += `<path d="${dArc(dAng(d, r[0]), dAng(d, r[1]), R + 11)}" stroke="${DK.dh}" stroke-width="6" fill="none" opacity=".38" stroke-linecap="round"/>`; }
  if (d.ref) { const v = d.ref(); if (v != null) g += radial(dAng(d, v), R - 9, R + 12, DK.ref, 2.5); }
  if (vd != null && isFinite(vd)) { const a = dAng(d, vd);
    g += `<path d="M${xy(dPt(a, R + 8))} L${xy(dPt(a - 5.5, R + 23))} L${xy(dPt(a + 5.5, R + 23))} Z" fill="${DK.dh}" stroke="${DK.bg}" stroke-width="1.2"/>`; }
  const needle = (v, col) => { const a = dAng(d, v), tip = dPt(a, R - 10), tail = dPt(a + 180, 12);
    return `<line x1="${tail[0].toFixed(1)}" y1="${tail[1].toFixed(1)}" x2="${tip[0].toFixed(1)}" y2="${tip[1].toFixed(1)}" stroke="${col}" stroke-width="3.5" stroke-linecap="round"/>`; };
  if (vs != null && isFinite(vs)) g += needle(vs, DK.surf);
  if (vn != null && isFinite(vn)) g += needle(vn, DK.dh);
  const big = vs != null ? vs : vn, col = vs != null ? DK.surf : DK.dh;
  g += `<circle cx="${CX}" cy="${CY}" r="7" fill="${DK.bg}" stroke="${col}" stroke-width="2.5"/>`;
  g += `<text x="${CX}" y="${CY + 60}" text-anchor="middle" font-size="25" font-weight="600" fill="${col}">${fx(big, d.dec)}</text>`;
  return g;
}

// ------------------------------------------------------------- bar drawing
function barFace(b) {
  let g = `<rect x="${BX}" y="${BT}" width="22" height="${BB - BT}" rx="3" fill="${DK.track}"/>`;
  if (b.id === 'mw') {
    if (b.frac < b.max) g += `<rect x="${BX}" y="${BT}" width="22" height="${(bY(b, b.frac) - BT).toFixed(1)}" fill="${DK.lim}" opacity=".7"/>`;
    if (b.pore > b.min) g += `<rect x="${BX}" y="${bY(b, b.pore).toFixed(1)}" width="22" height="${(BB - bY(b, b.pore)).toFixed(1)}" fill="${DK.lim}" opacity=".7"/>`;
    for (const v of [b.frac, b.pore]) if (v > b.min && v < b.max) g += `<line x1="${BX - 4}" y1="${bY(b, v).toFixed(1)}" x2="${BX + 26}" y2="${bY(b, v).toFixed(1)}" stroke="${DK.limMark}" stroke-width="2.5"/>`;
  }
  const n = Math.round((b.max - b.min) / b.step);
  for (let k = 0; k <= 2 * n; k++) {
    const v = b.min + k * b.step / 2, y = bY(b, v).toFixed(1), major = k % 2 === 0;
    g += `<line x1="${BX - (major ? 9 : 5)}" y1="${y}" x2="${BX - 1}" y2="${y}" stroke="${DK.tick}" stroke-width="${major ? 1.6 : 1}"/>`;
    if (major) g += `<text x="${BX - 12}" y="${(+y + 4).toFixed(1)}" text-anchor="end" font-size="11" fill="${DK.dim}">${tickLabel(v, b.step)}</text>`;
  }
  return g;
}
function arrowRight(y, label) {
  return `<path d="M${BX + 24} ${y.toFixed(1)} L${BX + 38} ${(y - 6).toFixed(1)} L${BX + 38} ${(y + 6).toFixed(1)} Z" fill="${DK.dh}" stroke="${DK.bg}" stroke-width="1"/>
    <text x="${BX + 41}" y="${(y + 4).toFixed(1)}" font-size="9.5" fill="${DK.dh}">${label}</text>`;
}
function barDynamic(b, L) {
  let g = '';
  const fill = (v) => (v == null || !isFinite(v) ? '' : `<rect x="${BX + 6}" y="${bY(b, v).toFixed(1)}" width="10" height="${(BB - bY(b, v)).toFixed(1)}" rx="2" fill="${DK.surf}"/>`);
  if (b.id === 'q') {
    g += fill(toDisp('flow', L.q_in));
    if (L.q_out != null) g += arrowRight(bY(b, toDisp('flow', L.q_out)), 'out');
  } else {
    g += fill(b.mw);
    if (L.ecd_hyd != null) g += arrowRight(bY(b, toDisp('ecd', L.ecd_hyd)), 'ECD');
  }
  return g;
}
function barHandle(b) {
  const y = bY(b, b.sp()).toFixed(1);
  return `<g class="sp-handle" tabindex="0" role="slider" aria-label="Pump rate set point" aria-valuenow="${fx(b.sp(), 0)}" transform="translate(${BX + 11} ${y})">
    <rect x="-20" y="-12" width="40" height="24" fill="transparent"/>
    <path class="vis" d="M-17 0 L-9 -7 L-9 7 Z M17 0 L9 -7 L9 7 Z" fill="${DK.sp}" stroke="${DK.bg}" stroke-width="1.2"/>
    <line x1="-9" y1="0" x2="9" y2="0" stroke="${DK.sp}" stroke-width="2"/></g>`;
}

// ------------------------------------------------------------------ build
function deckBuild(defs) {
  const host = $('deck'), legend = host.querySelector('.deck-legend');
  for (const el of host.querySelectorAll('.gauge')) el.remove();
  deck.nodes = []; deck.bars = []; deck.drag = null;
  for (const d of defs.dials) {
    const el = document.createElement('div'); el.className = 'gauge';
    const cells = [];
    if (d.surf) cells.push(`<div class="s"><span>${d.surfLabel || 'Surf'}</span><b data-k="s">–</b></div>`);
    if (d.dh || d.needleDh) cells.push(`<div class="d"><span>${d.dhLabel || 'DH'}</span><b data-k="d">–</b></div>`);
    if (d.ref) cells.push(`<div class="h"><span>${d.refLabel}</span><b data-k="r">–</b></div>`);
    if (d.sp) cells.push(`<div class="p"><span>Set</span><b data-k="p">–</b></div>`);
    else cells.push(`<div class="l"><span>Limit</span><b data-k="l">${fx(d.lim, d.dec)}</b></div>`);
    el.innerHTML = `<div class="g-head"><h2>${d.title}</h2>${d.chip ? `<button class="chip" type="button" title="${d.chipTitle || ''}"></button>` : ''}${d.chip2 ? `<button class="chip chip2" type="button" title="${d.chip2Title}">${d.chip2}</button>` : ''}</div>
      <svg viewBox="0 0 ${DW} ${DH}" role="img" aria-label="${d.title} gauge"><g>${dialFace(d)}</g><g class="dyn"></g><g class="spg"></g></svg>
      <div class="reads">${cells.join('')}</div>`;
    host.insertBefore(el, legend);
    const n = { d, el, svg: el.querySelector('svg'), dyn: el.querySelector('.dyn'), spg: el.querySelector('.spg'), chip: el.querySelector('.chip:not(.chip2)') };
    if (n.chip) n.chip.addEventListener('click', () => { d.chipClick(); deckRedraw(); });
    if (d.chip2) el.querySelector('.chip2').addEventListener('click', () => d.chip2Click());
    if (d.sp) wireDialDrag(n);
    deck.nodes.push(n);
  }
  for (const b of defs.bars) {
    const el = document.createElement('div'); el.className = 'gauge bar';
    const reads = b.id === 'q'
      ? '<div class="s"><span>In</span><b data-k="qi">–</b></div><div class="d"><span>Out</span><b data-k="qo">–</b></div><div class="p"><span>Set</span><b data-k="qp">–</b></div>'
      : '<div class="s"><span>Static</span><b data-k="mw">–</b></div><div class="d"><span>ECD</span><b data-k="ecd">–</b></div>';
    const foot = b.id === 'q'
      ? `<label class="tin">In <input type="number" id="t-in" step="any" aria-label="Mud inlet temperature"><span data-kind="temp"></span></label>`
      : `<div class="tin">${unitOf('ecd')}</div>`;
    el.innerHTML = `<div class="g-head"><h2>${b.title}</h2></div>
      <svg viewBox="0 0 ${BW} ${BH}" role="img" aria-label="${b.title} bar"><g>${barFace(b)}</g><g class="dyn"></g><g class="spg"></g></svg>
      <div class="reads">${reads}</div>${foot}`;
    host.insertBefore(el, legend);
    const n = { b, el, svg: el.querySelector('svg'), dyn: el.querySelector('.dyn'), spg: el.querySelector('.spg') };
    if (b.sp) wireBarDrag(n);
    deck.bars.push(n);
  }
  const tin = $('t-in');
  if (tin) {
    tin.value = ctl.tin == null ? '' : +toDisp('temp', ctl.tin).toPrecision(4);
    tin.addEventListener('change', () => { ctl.tin = tin.value === '' ? null : fromDisp('temp', Number(tin.value)); });
    tin.parentElement.querySelector('[data-kind]').textContent = unitOf('temp');
  }
}

// ---------------------------------------------------------------- drawing
function deckUpdate(L, t) {
  if (!L) return;
  deck.last = { L, t };
  const defs = deckDefs(L);
  const key = JSON.stringify([SYSTEM, defs.dials.map(d => [d.min, d.max, d.step, d.lim, d.limLo]), defs.bars.map(b => [b.min, b.max, b.step, b.pore, b.frac, b.mw])]);
  if (deck.drag) {
    // hold the scales still while a set point is being dragged; changes apply on release
  } else if (key !== deck.key) {
    const focus = document.activeElement && $('deck').contains(document.activeElement) && document.activeElement.id === 't-in';
    deckBuild(defs); deck.key = key;
    if (focus) $('t-in').focus();
  } else {
    deck.nodes.forEach((n, i) => { Object.assign(n.d, defs.dials[i]); });
    deck.bars.forEach((n, i) => { Object.assign(n.b, defs.bars[i]); });
  }
  for (const n of deck.nodes) { n.dyn.innerHTML = dialDynamic(n.d, L, t); if (deck.drag !== n) drawDialSp(n); }
  for (const n of deck.bars) { n.dyn.innerHTML = barDynamic(n.b, L); if (deck.drag !== n) drawBarSp(n); }
  deckReadouts(L);
}
function deckRedraw() { if (deck.last) deckUpdate(deck.last.L, deck.last.t); }

function deckReadouts(L) {
  for (const n of deck.nodes) {
    const d = n.d, q = (k) => n.el.querySelector(`[data-k="${k}"]`);
    if (q('s')) q('s').textContent = fx(d.surf(L), d.dec);
    if (q('d')) q('d').textContent = fx(d.dh ? d.dh(L) : d.needleDh(L), d.dec);
    if (q('r')) q('r').textContent = fx(d.ref(), d.dec);
    if (q('p')) q('p').textContent = d.spActive() ? fx(d.sp(), d.dec) : 'off';
    if (q('l')) q('l').textContent = fx(d.lim, d.dec);
    if (n.chip) { n.chip.textContent = d.chip(); n.chip.classList.toggle('on', d.chipOn()); }
  }
  const q = (k) => $('deck').querySelector(`[data-k="${k}"]`);
  const qb = deck.bars.find(n => n.b.id === 'q').b, mb = deck.bars.find(n => n.b.id === 'mw').b;
  q('qi').textContent = fx(toDisp('flow', L.q_in), 0); q('qo').textContent = fx(toDisp('flow', L.q_out), 0);
  q('qp').textContent = fx(qb.sp(), 0);
  const dd = SYSTEM === 'field' ? 2 : 0;
  q('mw').textContent = fx(mb.mw, dd); q('ecd').textContent = fx(toDisp('ecd', L.ecd_hyd), dd);
  // status strip and stick-slip alarm (bit speed range over 10 s, mean above 5 rpm)
  const r = traceRange('rpm_bit', 'rpm', deck.last && deck.last.t);
  $('alarm').hidden = !(r && r[2] > 5 && (r[1] - r[0]) / (2 * r[2]) > 0.5);
  $('st-mode').textContent = ctl.dwMode === 'wob' ? 'Autodriller' : (ctl.vblk === 0 ? 'Block stopped' : ctl.vblk > 0 ? 'Lowering' : 'Hoisting');
}

function drawDialSp(n) {
  const d = n.d; if (!d.sp) return;
  const had = document.activeElement && n.spg.contains(document.activeElement);
  const act = d.spActive();
  n.spg.innerHTML = `<g class="sp-handle${act ? '' : ' inactive'}" tabindex="${act ? 0 : -1}" role="slider" aria-label="${d.title} set point" aria-valuenow="${fx(d.sp(), d.dec)}">${diamondAt(dAng(d, d.sp()), act)}</g>`;
  if (had) n.spg.querySelector('.sp-handle').focus();
}
function drawBarSp(n) {
  if (!n.b.sp) return;
  const had = document.activeElement && n.spg.contains(document.activeElement);
  n.spg.innerHTML = barHandle(n.b);
  if (had) n.spg.querySelector('.sp-handle').focus();
}

// --------------------------------------------------------------- dragging
function svgPoint(svg, e) { const p = svg.createSVGPoint(); p.x = e.clientX; p.y = e.clientY; return p.matrixTransform(svg.getScreenCTM().inverse()); }
const snapTo = (v, s) => Math.round(v / s) * s;
function spStep(d) { return d.step / (d.step >= 10 ? 10 : 5); }

function wireDialDrag(n) {
  n.spg.addEventListener('pointerdown', (e) => {
    if (!n.d.spActive()) return;
    deck.drag = n; n.spg.setPointerCapture(e.pointerId); e.preventDefault();
  });
  n.spg.addEventListener('pointermove', (e) => {
    if (deck.drag !== n) return;
    const p = svgPoint(n.svg, e); let a = Math.atan2(p.x - CX, CY - p.y) * 180 / Math.PI;
    if (a > 157.5) a -= 360;                    // across the bottom gap: snap to the nearer end
    n.d.spSet(snapTo(dVal(n.d, Math.max(A0, Math.min(A0 + SW, a))), spStep(n.d)));
    drawDialSp(n); deckReadouts(deck.last.L);
  });
  const end = () => { if (deck.drag === n) { deck.drag = null; deckToast(n.d.spMsg()); deckRedraw(); } };
  for (const ev of ['pointerup', 'pointercancel', 'lostpointercapture']) n.spg.addEventListener(ev, end);
  n.spg.addEventListener('keydown', (e) => {
    if (!n.d.spActive()) return;
    const k = { ArrowUp: 1, ArrowRight: 1, ArrowDown: -1, ArrowLeft: -1, PageUp: 10, PageDown: -10 }[e.key]; if (!k) return;
    e.preventDefault(); e.stopPropagation();
    n.d.spSet(Math.max(n.d.min, Math.min(n.d.max, snapTo(n.d.sp(), spStep(n.d)) + k * spStep(n.d))));
    drawDialSp(n); deckReadouts(deck.last.L); deckToast(n.d.spMsg());
  });
}
function wireBarDrag(n) {
  n.spg.addEventListener('pointerdown', (e) => { deck.drag = n; n.spg.setPointerCapture(e.pointerId); e.preventDefault(); });
  n.spg.addEventListener('pointermove', (e) => {
    if (deck.drag !== n) return;
    n.b.spSet(snapTo(bV(n.b, svgPoint(n.svg, e).y), spStep(n.b))); drawBarSp(n); deckReadouts(deck.last.L);
  });
  const end = () => { if (deck.drag === n) { deck.drag = null; deckToast(n.b.spMsg()); deckRedraw(); } };
  for (const ev of ['pointerup', 'pointercancel', 'lostpointercapture']) n.spg.addEventListener(ev, end);
  n.spg.addEventListener('keydown', (e) => {
    const k = { ArrowUp: 1, ArrowRight: 1, ArrowDown: -1, ArrowLeft: -1, PageUp: 10, PageDown: -10 }[e.key]; if (!k) return;
    e.preventDefault(); e.stopPropagation();
    n.b.spSet(Math.max(n.b.min, Math.min(n.b.max, snapTo(n.b.sp(), spStep(n.b)) + k * spStep(n.b))));
    drawBarSp(n); deckReadouts(deck.last.L); deckToast(n.b.spMsg());
  });
}

// ---------------------------------------------------------------- actions
function setDwMode(mode) {
  if (ctl.dwMode === mode) return;
  ctl.dwMode = mode;
  // a new drawworks mode starts from rest: zero WOB set point, or a stopped block
  ctl.wob = 0; ctl.vblk = 0;
  const L = deck.last && deck.last.L;
  if (mode === 'wob' && L && !L.on_bottom) { ctl.tare = true; if (!state.playing) advance(0); }
  deckToast(mode === 'wob' ? 'Drawworks: autodriller (WOB set point 0; drag to set)' : 'Drawworks: block speed (stopped; drag to set)');
}
function tare() {
  ctl.tare = true;
  if (!state.playing) advance(0);
  deckToast('Weight indicator zeroed at the present hookload');
}
function deckToast(m) {
  const t = $('deck-toast'); t.textContent = m; t.classList.add('show');
  clearTimeout(deck.toastT); deck.toastT = setTimeout(() => t.classList.remove('show'), 2000);
}
// SI inputs for the server; the zeroing request is sent once
function deckInputs() {
  const out = { td_mode: ctl.tdMode, rpm: ctl.rpm, flow_m3s: ctl.flow, inlet_c: ctl.tin, dw_mode: ctl.dwMode,
                v_block_mh: ctl.dwMode === 'velocity' ? ctl.vblk * 3600 : 0, wob_kn: ctl.dwMode === 'wob' ? ctl.wob / 1e3 : 0,
                tare: ctl.tare };
  ctl.tare = false;
  return out;
}
