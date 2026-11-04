/* DrillPyze simulator UI.
 *
 * The browser edits a case (SI config, shown in convenient display units),
 * asks the Python server to initialise a simulation, and then plays it by
 * requesting one frame at a time.  Inputs in the transport bar are sent with
 * every step, so they can be changed while the simulation runs or is paused.
 */
'use strict';

const RPM = 2 * Math.PI / 60;
const COLORS = { ink: '#1c2733', steel: '#5b6b7c', line: '#d3dae2', rot: '#1f6fb2', trq: '#b7790a',
                 ax: '#2f7d4f', bit: '#b8326e', ref: '#9aa7b4' };
const MAX_FRAMES = 4000;

const $ = (id) => document.getElementById(id);
const state = {
  cfg: null, defaults: null, session: null, info: null, reference: null, path: null,
  frames: [], traces: null, playing: false, busy: false, scrub: null, windowS: 60,
  dirty: false, rpmMax: 200, tqMax: 1, tqMin: 0,
};

// ----------------------------------------------------------------- form spec
// [path in config, label, display unit, factor display->SI, options]
const SECTIONS = [
  { title: 'Well', open: true, note: 'Survey stations in metres and degrees. Casing above the shoe uses the cased-hole friction factors.',
    fields: [
      ['well.casing_shoe', 'Casing shoe', 'm', 1],
      ['well.hole_diameter', 'Hole diameter', 'mm', 1e-3],
      ['fluid.density', 'Mud density', 'kg/m³', 1],
    ], table: 'survey' },
  { title: 'Drillstring', open: true, note: 'Listed from the bit up. The last row is lengthened to reach surface. Leave mass blank for plain steel tube.',
    fields: [], table: 'string' },
  { title: 'Bit position', open: true,
    fields: [
      ['simulation.bit_depth', 'Bit depth', 'm', 1],
      ['simulation.hole_depth', 'Hole depth', 'm', 1, { hint: 'Equal to bit depth puts the bit on bottom' }],
    ] },
  { title: 'String friction',
    fields: [
      ['parameters.friction.mu_static_cased', 'Static, cased', '–', 1],
      ['parameters.friction.mu_kinetic_cased', 'Kinetic, cased', '–', 1],
      ['parameters.friction.mu_static_open', 'Static, open hole', '–', 1],
      ['parameters.friction.mu_kinetic_open', 'Kinetic, open hole', '–', 1],
      ['parameters.friction.omega_c', 'Static friction below, rotation', 'rpm', RPM],
      ['parameters.friction.v_c', 'Static friction below, axial', 'mm/s', 1e-3],
      ['parameters.friction.along_string_coupling', 'Couple axial and rotation in friction', 'bool', 1],
    ] },
  { title: 'Damping',
    fields: [
      ['parameters.damping.k_t', 'Torsional', '1/s', 1],
      ['parameters.damping.k_a', 'Axial', '1/s', 1],
    ] },
  { title: 'Bit',
    note: 'Sticking: |TOB| ≤ T₀ + μs·r·WOB. Sliding: TOB = μk·r·WOB. ROP = k·(WOB − W₀)·|ω|.',
    fields: [
      ['parameters.bit.radius', 'Effective radius r', 'mm', 1e-3, { blank: 'D/3' }],
      ['parameters.bit.mu_static', 'Static friction μs', '–', 1],
      ['parameters.bit.mu_kinetic', 'Kinetic friction μk', '–', 1],
      ['parameters.bit.stick_torque', 'Sticking torque T₀', 'kN·m', 1e3],
      ['parameters.bit.k_rop', 'Penetration coefficient k', 'm/(N·rad)', 1, { step: 1e-9 }],
      ['parameters.bit.wob_threshold', 'Threshold WOB W₀', 'kN', 1e3],
    ] },
  { title: 'Top drive and controller',
    fields: [
      ['parameters.controller.type', 'Controller', 'select', 1,
        { options: [['stiff', 'Stiff PI'], ['soft_torque', 'SoftTorque'], ['z_torque', 'ZTorque'], ['custom', 'Custom PI']] }],
      ['parameters.controller.f_c', 'SoftTorque tuning frequency', 'Hz', 1, { blank: 'auto' }],
      ['parameters.controller.t_hp', 'ZTorque high-pass', 's', 1],
      ['parameters.controller.kp', 'Override Kp', 'N·m·s', 1, { blank: 'auto' }],
      ['parameters.controller.ki', 'Override Ki', 'N·m', 1, { blank: 'auto' }],
      ['parameters.top_drive.inertia', 'Top drive inertia', 'kg·m²', 1],
      ['parameters.top_drive.max_torque', 'Torque limit', 'kN·m', 1e3, { hint: '0 for no limit' }],
      ['parameters.top_drive.setpoint_ramp', 'Set point ramp', 'rpm/s', RPM],
    ] },
  { title: 'Numerics',
    fields: [
      ['simulation.dx', 'Cell size', 'm', 1, { hint: 'Halving the cell size makes runs about 4× slower' }],
      ['simulation.record_rate', 'Recording rate', 'Hz', 1],
      ['parameters.cfl', 'Courant number', '–', 1],
    ] },
];

const SURVEY_COLS = [['md', 'MD', 'm', 1], ['inc', 'Inc', '°', 1], ['azi', 'Azi', '°', 1]];
const STRING_COLS = [['name', 'Name', '', 'text'], ['length', 'Length', 'm', 1], ['od', 'OD', 'mm', 1e-3],
                     ['id', 'ID', 'mm', 1e-3], ['linear_mass', 'Mass', 'kg/m', 1], ['contact_od', 'Contact OD', 'mm', 1e-3]];

// ------------------------------------------------------------------ helpers
function get(obj, path) { return path.split('.').reduce((o, k) => (o == null ? o : o[k]), obj); }
function set(obj, path, v) {
  const ks = path.split('.'); const last = ks.pop();
  const o = ks.reduce((o, k) => (o[k] = o[k] ?? {}), obj); o[last] = v;
}
function fmt(x, nd = 3) {
  if (x == null || x === '' || Number.isNaN(x)) return '';
  if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e6)) return Number(x).toExponential(3);
  return String(+Number(x).toFixed(nd));
}
function fixed(x, nd) { return (x == null || !isFinite(x)) ? '–' : Number(x).toFixed(nd); }
function message(text, kind = 'error') {
  const m = $('message');
  if (!text) { m.hidden = true; return; }
  m.textContent = text; m.className = 'message' + (kind === 'info' ? ' info' : ''); m.hidden = false;
}
async function api(path, body) {
  const r = await fetch(path, body === undefined ? {} :
    { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
  const data = await r.json().catch(() => ({ error: `Server error ${r.status}` }));
  if (!r.ok) throw new Error(data.error || `Server error ${r.status}`);
  return data;
}
function markDirty() { state.dirty = true; $('dirty-hint').hidden = false; }

// ------------------------------------------------------------------- form
function fieldRow(spec) {
  const [path, label, unit, factor, opt = {}] = spec;
  const row = document.createElement('div');
  const id = 'f-' + path.replace(/\./g, '-');
  if (unit === 'bool') {
    row.className = 'field check';
    row.innerHTML = `<input type="checkbox" id="${id}"><label for="${id}">${label}</label>`;
    const el = row.querySelector('input');
    el.checked = !!get(state.cfg, path);
    el.addEventListener('change', () => { set(state.cfg, path, el.checked); markDirty(); });
    return row;
  }
  row.className = 'field';
  if (unit === 'select') {
    row.innerHTML = `<label for="${id}">${label}</label><select id="${id}"></select><span></span>`;
    const el = row.querySelector('select');
    for (const [v, t] of opt.options) el.add(new Option(t, v));
    el.value = get(state.cfg, path);
    el.addEventListener('change', () => { set(state.cfg, path, el.value); markDirty(); });
    return row;
  }
  row.innerHTML = `<label for="${id}">${label}</label><input type="number" id="${id}" step="${opt.step ?? 'any'}"><span class="unit">${unit}</span>`;
  const el = row.querySelector('input');
  const v = get(state.cfg, path);
  el.value = v == null ? '' : fmt(v / factor, 5);
  if (opt.blank) el.placeholder = opt.blank;
  if (opt.hint) el.title = opt.hint;
  el.addEventListener('change', () => {
    set(state.cfg, path, el.value === '' ? null : Number(el.value) * factor); markDirty();
  });
  return row;
}

function tableEditor(kind) {
  const cols = kind === 'survey' ? SURVEY_COLS : STRING_COLS;
  const rows = () => (kind === 'survey' ? state.cfg.well.survey : state.cfg.drillstring);
  const wrap = document.createElement('div');
  const draw = () => {
    wrap.innerHTML = '';
    const t = document.createElement('table'); t.className = 'grid' + (kind === 'string' ? ' wide' : '');
    t.innerHTML = '<thead><tr>' + cols.map(c => `<th>${c[1]}${c[2] ? ` <span>(${c[2]})</span>` : ''}</th>`).join('') + '<th></th></tr></thead>';
    const tb = document.createElement('tbody');
    rows().forEach((r, i) => {
      const tr = document.createElement('tr');
      for (const [key, , , factor] of cols) {
        const td = document.createElement('td');
        const inp = document.createElement('input');
        if (factor === 'text') { inp.type = 'text'; inp.value = r[key] ?? ''; }
        else { inp.type = 'number'; inp.step = 'any'; inp.value = r[key] == null ? '' : fmt(r[key] / factor, 2); }
        inp.setAttribute('aria-label', `${key} row ${i + 1}`);
        inp.addEventListener('change', () => {
          r[key] = factor === 'text' ? inp.value : (inp.value === '' ? null : Number(inp.value) * factor);
          markDirty();
        });
        td.appendChild(inp); tr.appendChild(td);
      }
      const td = document.createElement('td');
      const del = document.createElement('button'); del.className = 'del'; del.textContent = '×';
      del.setAttribute('aria-label', `Remove row ${i + 1}`);
      del.addEventListener('click', () => { rows().splice(i, 1); markDirty(); draw(); });
      td.appendChild(del); tr.appendChild(td); tb.appendChild(tr);
    });
    t.appendChild(tb);
    const tw = document.createElement('div'); tw.className = 'table-wrap'; tw.appendChild(t); wrap.appendChild(tw);
    const act = document.createElement('div'); act.className = 'table-actions';
    const add = document.createElement('button'); add.className = 'btn quiet'; add.textContent = 'Add row';
    add.addEventListener('click', () => {
      const last = rows()[rows().length - 1] || {};
      rows().push(kind === 'survey' ? { md: (last.md || 0) + 100, inc: last.inc || 0, azi: last.azi || 0 }
                                    : { ...last, name: 'New component' });
      markDirty(); draw();
    });
    const imp = document.createElement('label'); imp.className = 'btn quiet file';
    imp.textContent = 'Import CSV';
    const fi = document.createElement('input'); fi.type = 'file'; fi.accept = '.csv,text/csv';
    fi.addEventListener('change', async () => {
      try {
        const parsed = parseCSV(await fi.files[0].text(), cols);
        if (kind === 'survey') state.cfg.well.survey = parsed; else state.cfg.drillstring = parsed;
        markDirty(); draw(); message(`Imported ${parsed.length} rows.`, 'info');
      } catch (e) { message(e.message); }
      fi.value = '';
    });
    imp.appendChild(fi);
    act.append(add, imp); wrap.appendChild(act);
  };
  draw();
  return wrap;
}

function parseCSV(text, cols) {
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
  const head = lines.shift().split(',').map(h => h.trim().toLowerCase().replace(/\s*\(.*\)\s*$/, ''));
  const idx = cols.map(([key]) => head.indexOf(key));
  const need = cols.filter(([, , , f], k) => f !== 'text' && idx[k] < 0 && !['linear_mass', 'contact_od', 'azi'].includes(cols[k][0]));
  if (need.length) throw new Error(`CSV is missing column(s): ${need.map(c => c[0]).join(', ')}. Units: ${cols.map(c => `${c[0]}${c[2] ? ' ' + c[2] : ''}`).join(', ')}.`);
  return lines.map(l => {
    const v = l.split(',').map(s => s.trim()); const r = {};
    cols.forEach(([key, , , f], k) => {
      const s = idx[k] >= 0 ? v[idx[k]] : '';
      r[key] = f === 'text' ? s : (s === '' ? (key === 'azi' ? 0 : null) : Number(s) * f);
    });
    return r;
  });
}

function buildForm() {
  const form = $('form'); form.innerHTML = '';
  for (const sec of SECTIONS) {
    const d = document.createElement('details'); d.className = 'sec'; d.open = !!sec.open;
    d.innerHTML = `<summary>${sec.title}</summary>`;
    const body = document.createElement('div'); body.className = 'sec-body';
    if (sec.note) { const p = document.createElement('p'); p.className = 'sec-note'; p.textContent = sec.note; body.appendChild(p); }
    for (const f of sec.fields) body.appendChild(fieldRow(f));
    if (sec.table) body.appendChild(tableEditor(sec.table));
    d.appendChild(body); form.appendChild(d);
  }
  $('case-name').textContent = state.cfg.name || 'Custom case';
}

// ------------------------------------------------------------------- plots
const PCONF = { displayModeBar: false, responsive: true };
const TITLE = (text) => ({ text, x: 0.02, xanchor: 'left', y: 0.985, yanchor: 'top', yref: 'container', font: { size: 15 } });
const LEGEND_TOP = { orientation: 'h', x: 0, xanchor: 'left', y: 1.0, yanchor: 'bottom', font: { size: 11 }, traceorder: 'reversed' };
const FONT = { family: 'Barlow, Segoe UI, Arial, sans-serif', size: 12, color: COLORS.ink };
function baseLayout(extra) {
  return Object.assign({
    font: FONT, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    margin: { l: 56, r: 14, t: 34, b: 40 }, showlegend: false, uirevision: 'keep',
  }, extra);
}
function axis(title, extra) {
  return Object.assign({ title: { text: title, standoff: 6 }, gridcolor: '#e6ebf0', zerolinecolor: COLORS.line,
                         linecolor: COLORS.line, ticks: 'outside', tickcolor: COLORS.line }, extra);
}
function depthAxis(show = true) {
  const r = state.reference ? [state.reference.md[state.reference.md.length - 1] * 1.01, 0] : undefined;
  return axis(show ? 'Measured depth (m)' : '', { range: r, autorange: false, showticklabels: show });
}

function drawStringViews(fr) {
  const p = fr.profile;
  const ref = state.reference;
  // ---- RPM along the string
  state.rpmMax = Math.max(state.rpmMax, ...p.rpm.map(Math.abs));
  Plotly.react('p-rpm', [
    { x: p.rpm, y: p.md, mode: 'lines', line: { color: COLORS.rot, width: 2 }, hovertemplate: '%{x:.0f} rpm at %{y:.0f} m<extra></extra>' },
  ], baseLayout({ title: TITLE('Rotary speed'), margin: { l: 56, r: 14, t: 74, b: 40 },
    xaxis: axis('rpm', { range: [-0.08 * state.rpmMax, state.rpmMax * 1.05] }), yaxis: depthAxis(true) }), PCONF);

  // ---- torque along the string with T&D references
  const tq = p.torque_knm;
  state.tqMax = Math.max(state.tqMax, ...tq, ...ref.torque_breakoff_knm);
  state.tqMin = Math.min(state.tqMin, ...tq);
  Plotly.react('p-torque', [
    { x: ref.torque_breakoff_knm, y: ref.md, mode: 'lines', line: { color: COLORS.ref, dash: 'dot', width: 1.5 }, name: 'Static μ', hoverinfo: 'skip' },
    { x: ref.torque_rotating_knm, y: ref.md, mode: 'lines', line: { color: COLORS.ref, dash: 'dash', width: 1.5 }, name: 'Kinetic μ', hoverinfo: 'skip' },
    { x: tq, y: p.md, mode: 'lines', line: { color: COLORS.trq, width: 2 }, name: 'Simulated', hovertemplate: '%{x:.2f} kN·m at %{y:.0f} m<extra></extra>' },
  ], baseLayout({ title: TITLE('Torque'), showlegend: true,
    legend: LEGEND_TOP, margin: { l: 14, r: 14, t: 74, b: 40 },
    xaxis: axis('kN·m', { range: [Math.min(0, state.tqMin) * 1.1, state.tqMax * 1.08] }), yaxis: depthAxis(false) }), PCONF);

  // ---- axial force
  Plotly.react('p-tension', [
    { x: ref.tension_hanging_kn, y: ref.md, mode: 'lines', line: { color: COLORS.ref, dash: 'dash', width: 1.5 }, name: 'Hanging', hoverinfo: 'skip' },
    { x: p.tension_kn, y: p.md, mode: 'lines', line: { color: COLORS.ax, width: 2 }, name: 'Simulated', hovertemplate: '%{x:.0f} kN at %{y:.0f} m<extra></extra>' },
  ], baseLayout({ title: TITLE('Effective tension'), showlegend: true,
    legend: LEGEND_TOP, margin: { l: 14, r: 14, t: 74, b: 40 },
    xaxis: axis('kN'), yaxis: depthAxis(false) }), PCONF);

  // ---- well path with the drill pipe, BHA and bit positions
  const path = state.path;
  const cased = path.md.map((m, i) => (m <= path.shoe_md ? i : -1)).filter(i => i >= 0);
  const pick = (arr, want) => arr.map((v, i) => (path.bha[i] === want ? v : null));
  const nb = p.md.length - 1;
  const bitTip = `Bit at ${fr.latest.bit_depth.toFixed(1)} m MD`;
  Plotly.react('p-path', [
    { x: path.vs, y: path.tvd, mode: 'lines', line: { color: '#d9e0e7', width: 16 }, hoverinfo: 'skip', name: 'Open hole' },
    { x: cased.map(i => path.vs[i]), y: cased.map(i => path.tvd[i]), mode: 'lines', line: { color: '#a9b5c1', width: 18 },
      hoverinfo: 'skip', name: `Casing to ${path.shoe_md.toFixed(0)} m` },
    { x: pick(p.vs, 0), y: pick(p.tvd, 0), mode: 'lines', line: { color: COLORS.steel, width: 3 }, name: 'Drill pipe',
      hoverinfo: 'skip' },
    { x: pick(p.vs, 1), y: pick(p.tvd, 1), mode: 'lines', line: { color: COLORS.ink, width: 8 },
      name: `BHA (${path.bha_length.toFixed(0)} m)`, hoverinfo: 'skip' },
    { x: [p.vs[nb]], y: [p.tvd[nb]], mode: 'markers', name: 'Bit', hovertemplate: bitTip + '<extra></extra>',
      marker: { symbol: 'diamond', size: 13, color: COLORS.trq, line: { color: COLORS.ink, width: 1 } } },
  ], baseLayout({ title: TITLE('Well path'), showlegend: true,
    legend: { x: 0.98, xanchor: 'right', y: 0.98, yanchor: 'top', bgcolor: 'rgba(255,255,255,0.85)', font: { size: 11 } },
    xaxis: axis('Horizontal displacement (m)', { scaleanchor: 'y', scaleratio: 1 }),
    yaxis: axis('True vertical depth (m)', { autorange: 'reversed' }), margin: { l: 64, r: 14, t: 40, b: 40 } }), PCONF);
}

function emptyTraces() {
  return { time: [], rpm_td: [], rpm_bit: [], rpm_sp: [], motor_torque_knm: [], pipe_torque_knm: [], tob_knm: [],
           hookload_kn: [], wob_kn: [], rop_mh: [], block_mh: [], bit_depth: [], hole_depth: [] };
}
function appendTrace(t) { for (const k in state.traces) { const a = t[k]; if (a) for (let i = 0; i < a.length; i++) state.traces[k].push(a[i]); } }

function timeRange(tNow) {
  const W = state.windowS;
  if (!W) return [0, Math.max(tNow, 1)];
  return tNow <= W ? [0, W] : [tNow - W, tNow];
}

function drawTraces(cursor) {
  const T = state.traces; const tNow = T.time.length ? T.time[T.time.length - 1] : 0;
  const range = timeRange(tNow);
  const shapes = cursor == null ? [] : [{ type: 'line', x0: cursor, x1: cursor, yref: 'paper', y0: 0, y1: 1, line: { color: COLORS.ink, width: 1, dash: 'dot' } }];
  const lay = (title, yt, extra = {}) => baseLayout(Object.assign({
    title: { text: title, x: 0.01, font: { size: 14 } }, showlegend: true,
    legend: { orientation: 'h', x: 1, xanchor: 'right', y: 1.18, font: { size: 11 } },
    margin: { l: 56, r: 56, t: 30, b: 30 }, xaxis: axis('', { range }), yaxis: axis(yt), shapes,
    datarevision: T.time.length }, extra));
  const sg = (y, name, color, extra = {}) => Object.assign({ x: T.time, y, name, type: 'scattergl', mode: 'lines', line: { color, width: 1.4 } }, extra);
  Plotly.react('t-rpm', [
    sg(T.rpm_sp, 'Set point', COLORS.ref, { line: { color: COLORS.ref, width: 1.2, dash: 'dash' } }),
    sg(T.rpm_td, 'Top drive', COLORS.ink),
    sg(T.rpm_bit, 'Bit', COLORS.rot),
  ], lay('Rotary speed', 'rpm'), PCONF);
  Plotly.react('t-torque', [
    sg(T.motor_torque_knm, 'Motor', COLORS.ink),
    sg(T.pipe_torque_knm, 'Pipe at surface', COLORS.trq),
    sg(T.tob_knm, 'On bit', COLORS.bit),
  ], lay('Torque', 'kN·m'), PCONF);
  Plotly.react('t-axial', [
    sg(T.hookload_kn, 'Hookload', COLORS.ax),
    sg(T.wob_kn, 'WOB', COLORS.ink, { yaxis: 'y2' }),
  ], lay('Axial', 'Hookload (kN)', { xaxis: axis('Time (s)', { range }),
    yaxis2: axis('WOB (kN)', { overlaying: 'y', side: 'right', showgrid: false, rangemode: 'tozero' }) }), PCONF);
}

function showReadouts(L) {
  $('ro-rpm-td').textContent = fixed(L.rpm_td, 0);
  $('ro-rpm-bit').textContent = fixed(L.rpm_bit, 0);
  $('ro-tq').textContent = fixed(L.pipe_torque_knm, 2);
  $('ro-tob').textContent = fixed(L.tob_knm, 2);
  $('ro-hk').textContent = fixed(L.hookload_kn, 0);
  $('ro-wob').textContent = fixed(L.wob_kn, 1);
  $('ro-rop').textContent = fixed(L.rop_mh, 1);
  $('ro-depth').textContent = fixed(L.bit_depth, 2);
}

function render(fr, live = true) {
  drawStringViews(fr);
  showReadouts(fr.latest);
  $('clock').textContent = `${fr.time.toFixed(2)} s`;
  drawTraces(live ? null : fr.time);
}

// -------------------------------------------------------------- simulation
function inputs() {
  return { td_mode: $('td-mode').value, rpm: Number($('rpm').value) || 0,
           dw_mode: $('dw-mode').value,
           v_block_mh: $('dw-mode').value === 'velocity' ? Number($('dw').value) || 0 : 0,
           wob_kn: $('dw-mode').value === 'wob' ? Number($('dw').value) || 0 : 0 };
}

async function initialise() {
  pause();
  message(null);
  $('apply').disabled = true; $('apply').textContent = 'Starting…';
  try {
    const r = await api('/api/init', { config: state.cfg, previous: state.session });
    Object.assign(state, { session: r.session, info: r.info, reference: r.reference, path: r.path,
                           frames: [], traces: emptyTraces(), scrub: null, rpmMax: Math.max(60, Number($('rpm').value) * 1.5),
                           tqMax: 1, tqMin: 0, dirty: false });
    $('dirty-hint').hidden = true;
    pushFrame(r.frame);
    render(r.frame);
    const i = r.info;
    $('speed-out').textContent = `${i.cells} cells, Δt ${(i.dt * 1e3).toFixed(2)} ms, ${i.controller.replace('_', ' ')}`;
    $('download').href = `/api/results.csv?session=${r.session}`;
  } catch (e) {
    message(e.message);
  } finally {
    $('apply').disabled = false; $('apply').textContent = 'Apply and reset';
  }
}

function pushFrame(fr) {
  appendTrace(fr.trace);
  state.frames.push({ time: fr.time, latest: fr.latest, profile: fr.profile });
  if (state.frames.length > MAX_FRAMES) state.frames.splice(0, state.frames.length - MAX_FRAMES);
  const s = $('scrub'); s.max = state.frames.length - 1; if (state.scrub == null) s.value = s.max;
}

async function advance(duration) {
  if (!state.session) return;
  state.busy = true;
  try {
    const fr = await api('/api/step', { session: state.session, duration, inputs: inputs() });
    pushFrame(fr);
    if (state.scrub == null) render(fr);
    const rt = fr.perf.realtime;
    $('speed-out').textContent = rt ? `simulating at ${rt >= 10 ? rt.toFixed(0) : rt.toFixed(1)}× real time` : '';
    return fr;
  } catch (e) {
    message(e.message); pause();
  } finally { state.busy = false; }
}

async function loop() {
  while (state.playing) {
    const dt = Number($('frame-dt').value);
    const speed = Number($('speed').value);
    const t0 = performance.now();
    await advance(dt);
    const wait = speed > 0 ? dt * 1000 / speed - (performance.now() - t0) : 0;
    await new Promise(res => (wait > 1 ? setTimeout(res, wait) : requestAnimationFrame(res)));
  }
}

function play() {
  if (state.playing || !state.session) return;
  if (state.dirty) message('Inputs have changed since the last reset. The run continues with the previous case until you apply.', 'info');
  goLive();
  state.playing = true;
  $('play').textContent = 'Pause'; $('play').classList.add('playing'); $('play').setAttribute('aria-label', 'Pause');
  loop();
}
function pause() {
  state.playing = false;
  $('play').textContent = 'Play'; $('play').classList.remove('playing'); $('play').setAttribute('aria-label', 'Play');
}

async function runAhead() {
  const total = Math.max(0, Number($('ahead').value) || 0);
  if (!total || !state.session) return;
  pause(); goLive();
  const chunk = Math.min(5, Math.max(1, total / 20));
  const btn = $('run-ahead'); btn.disabled = true;
  let done = 0;
  while (done < total - 1e-9) {
    const d = Math.min(chunk, total - done);
    if (!(await advance(d))) break;
    done += d; btn.textContent = `Running… ${Math.round(100 * done / total)}%`;
  }
  btn.disabled = false; btn.textContent = 'Run ahead (s)';
}

// ------------------------------------------------------------------ replay
function scrubTo(i) {
  pause();
  state.scrub = i;
  const fr = state.frames[i];
  render(fr, false);
  $('scrub-label').textContent = `Replay at ${fr.time.toFixed(2)} s`;
  $('go-live').hidden = false;
}
function goLive() {
  if (state.scrub == null) return;
  state.scrub = null;
  $('scrub').value = $('scrub').max;
  $('scrub-label').textContent = 'Live'; $('go-live').hidden = true;
  const fr = state.frames[state.frames.length - 1]; if (fr) render(fr);
}

// ------------------------------------------------------------- live inputs
function syncLive(source, target) { $(target).value = $(source).value; }
function dwModeChanged() {
  const wob = $('dw-mode').value === 'wob';
  const s = $('dw-slider');
  if (wob) { s.min = 0; s.max = 300; s.step = 1; $('dw-unit').textContent = 'kN'; }
  else { s.min = -60; s.max = 60; s.step = 0.5; $('dw-unit').textContent = 'm/h'; }
  $('dw').value = 0; s.value = 0;
}

// -------------------------------------------------------------- case files
function download(name, text, type = 'application/json') {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// ------------------------------------------------------------------ wiring
function wire() {
  $('apply').addEventListener('click', initialise);
  $('reset').addEventListener('click', initialise);
  $('play').addEventListener('click', () => (state.playing ? pause() : play()));
  $('step').addEventListener('click', () => { pause(); goLive(); advance(Number($('frame-dt').value)); });
  $('run-ahead').addEventListener('click', runAhead);
  $('restore').addEventListener('click', () => { state.cfg = structuredClone(state.defaults); buildForm(); initialise(); });
  $('save-case').addEventListener('click', () => download('drillpyze_case.json', JSON.stringify(state.cfg, null, 2)));
  $('open-case').addEventListener('change', async (ev) => {
    try {
      state.cfg = JSON.parse(await ev.target.files[0].text());
      buildForm(); await initialise(); message(`Opened ${ev.target.files[0].name}.`, 'info');
    } catch (e) { message(`Could not open the case file: ${e.message}`); }
    ev.target.value = '';
  });
  $('collapse-case').addEventListener('click', () => { document.querySelector('.app').classList.add('case-hidden'); $('show-case').hidden = false; window.dispatchEvent(new Event('resize')); });
  $('show-case').addEventListener('click', () => { document.querySelector('.app').classList.remove('case-hidden'); $('show-case').hidden = true; window.dispatchEvent(new Event('resize')); });

  $('rpm-slider').addEventListener('input', () => syncLive('rpm-slider', 'rpm'));
  $('rpm').addEventListener('input', () => syncLive('rpm', 'rpm-slider'));
  $('dw-slider').addEventListener('input', () => syncLive('dw-slider', 'dw'));
  $('dw').addEventListener('input', () => syncLive('dw', 'dw-slider'));
  $('dw-mode').addEventListener('change', dwModeChanged);

  $('scrub').addEventListener('input', (e) => scrubTo(Number(e.target.value)));
  $('go-live').addEventListener('click', goLive);

  for (const b of document.querySelectorAll('#window button')) {
    b.addEventListener('click', () => {
      for (const o of document.querySelectorAll('#window button')) o.setAttribute('aria-checked', String(o === b));
      state.windowS = Number(b.dataset.w);
      drawTraces(state.scrub == null ? null : state.frames[state.scrub].time);
    });
  }
  document.addEventListener('keydown', (e) => {
    if (e.target.matches('input, select, textarea')) return;
    if (e.code === 'Space') { e.preventDefault(); state.playing ? pause() : play(); }
    if (e.code === 'ArrowRight' && !state.playing) { advance(Number($('frame-dt').value)); }
  });
}

async function start() {
  wire();
  if (typeof Plotly === 'undefined') {
    message('The plotting library could not be loaded. Install it with "pip install plotly" or connect to the internet, then reload.');
    return;
  }
  try {
    const d = await api('/api/defaults');
    state.defaults = d.config; state.cfg = structuredClone(d.config);
    buildForm();
    await initialise();
    message('Well A is loaded. Press Play (or Space) to start the off-bottom start-up at 60 rpm.', 'info');
  } catch (e) { message(`Could not reach the DrillPyze server: ${e.message}`); }
}
start();
