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


// ------------------------------------------------------------------- units
// Every quantity kind has a display unit in each system and the number of SI
// units per display unit.  The case config and the server stay in SI; the UI
// converts at the edges, so switching systems never changes the case.
const UNITS = {
  si: {
    length: ['m', 1], depth: ['m', 1], diameter: ['mm', 1e-3], density: ['kg/m³', 1],
    linear_mass: ['kg/m', 1], force: ['kN', 1e3], torque: ['kN·m', 1e3], rop: ['m/h', 1 / 3600],
    slow_speed: ['mm/s', 1e-3], inertia: ['kg·m²', 1], k_rop: ['(m/h)/(kN·rpm)', (1 / 3600) / (1e3 * RPM)], stress: ['MPa', 1e6],
    modulus: ['GPa', 1e9], pgrad: ['kPa/m', 1e3], flow: ['L/min', 1 / 60000], tgrad: ['°C/km', 1e-3], per_temp: ['1/°C', 1],
    temp: ['°C', 1], dtemp: ['°C', 1], yp: ['Pa', 1], area: ['mm²', 1e-6], pressure: ['MPa', 1e6], ecd: ['kg/m³', 1],
    kp: ['N·m·s', 1], ki: ['N·m', 1],
  },
  field: {
    length: ['ft', 0.3048], depth: ['ft', 0.3048], diameter: ['in', 0.0254], density: ['ppg', 119.826427],
    linear_mass: ['lb/ft', 1.48816394], force: ['klbf', 4448.2216], torque: ['kft·lbf', 1355.8179],
    rop: ['ft/h', 0.3048 / 3600], slow_speed: ['in/s', 0.0254], inertia: ['lbm·ft²', 0.0421401],
    k_rop: ['(ft/h)/(klbf·rpm)', (0.3048 / 3600) / (4448.2216 * RPM)], stress: ['psi', 6894.7573], modulus: ['Mpsi', 6894757.3], pgrad: ['psi/ft', 6894.7573 / 0.3048],
    flow: ['gpm', 6.30901964e-5], tgrad: ['°F/100 ft', (5 / 9) / 30.48], per_temp: ['1/°F', 1.8], temp: ['°F', 5 / 9], dtemp: ['°F', 5 / 9], yp: ['lbf/100 ft²', 0.4788026], area: ['in²', 6.4516e-4], pressure: ['psi', 6894.7573], ecd: ['ppg', 119.826427], kp: ['ft·lbf·s', 1.3558179], ki: ['ft·lbf', 1.3558179],
  },
};
// kinds whose unit does not change between systems: [label, SI per display unit]
const FIXED = { accel: ['m/s²', 1], cond: ['W/(m·K)', 1], cp: ['J/(kg·K)', 1], visc: ['mPa·s', 1e-3], none: ['–', 1], deg: ['°', 1], rpm: ['rpm', RPM], rpm_s: ['rpm/s', RPM], per_s: ['1/s', 1],
                hz: ['Hz', 1], s: ['s', 1] };
const SERVER = { force: 1e3, torque: 1e3, length: 1, depth: 1, rop: 1 / 3600, temp: 1, stress: 1 };   // SI per unit the server sends

let SYSTEM = 'si';
try { SYSTEM = localStorage.getItem('drillpyze-units') || 'si'; } catch (e) { /* storage unavailable */ }

function unitOf(kind) { return (FIXED[kind] || UNITS[SYSTEM][kind])[0]; }
function perDisp(kind) { return (FIXED[kind] || UNITS[SYSTEM][kind])[1]; }
// temperatures are stored in deg C; the field display (deg F) needs an offset
const OFFSET = { temp: () => (SYSTEM === 'field' ? 32 : 0) };
const toDisp = (kind, si) => (si == null ? null : si / perDisp(kind) + (OFFSET[kind] ? OFFSET[kind]() : 0));
const fromDisp = (kind, d) => (d == null ? null : (d - (OFFSET[kind] ? OFFSET[kind]() : 0)) * perDisp(kind));
// convert an array (or number) the server sent in its fixed units into display units
function srv(kind, a) {
  const k = SERVER[kind] / perDisp(kind), off = OFFSET[kind] ? OFFSET[kind]() : 0;
  return Array.isArray(a) ? a.map(v => (v == null ? null : v * k + off)) : (a == null ? a : a * k + off);
}

const $ = (id) => document.getElementById(id);
const state = {
  cfg: null, defaults: null, session: null, info: null, reference: null, path: null,
  frames: [], traces: null, slow: null, playing: false, busy: false, scrub: null, windowS: 60,
  dirty: false, rpmMax: 200, tqMax: 1, tqMin: 0,
  scales: null, autoScale: {},
};
state.scales = loadScales();

// ----------------------------------------------------------------- form spec
// [path in config, label, quantity kind, options]
const SECTIONS = [
  { title: 'Well', open: true, note: () => `Survey stations in ${unitOf('depth') === 'ft' ? 'feet' : 'metres'} and degrees. Casing above the shoe uses the cased-hole friction factors.`,
    fields: [
      ['well.casing_shoe', 'Casing shoe', 'depth'],
      ['well.hole_diameter', 'Hole diameter', 'diameter'],
      ['fluid.density', 'Mud density', 'density'],
    ], table: 'survey' },
  { title: 'Drillstring', open: true, note: () => 'Listed from the bit up. The last row is lengthened to reach surface. Leave mass blank for plain steel tube.',
    fields: [], table: 'string' },
  { title: 'Bit position', open: true,
    fields: [
      ['simulation.bit_depth', 'Bit depth', 'depth'],
      ['simulation.hole_depth', 'Hole depth', 'depth', { hint: 'Equal to bit depth puts the bit on bottom' }],
    ] },
  { title: 'String friction',
    fields: [
      ['parameters.friction.mu_static_cased', 'Static, cased', 'none'],
      ['parameters.friction.mu_kinetic_cased', 'Kinetic, cased', 'none'],
      ['parameters.friction.mu_static_open', 'Static, open hole', 'none'],
      ['parameters.friction.mu_kinetic_open', 'Kinetic, open hole', 'none'],
      ['parameters.friction.omega_c', 'Static friction below, rotation', 'rpm'],
      ['parameters.friction.v_c', 'Static friction below, axial', 'slow_speed'],
      ['parameters.friction.along_string_coupling', 'Couple axial and rotation in friction', 'bool'],
      ['parameters.friction.axial_static_contrast', 'Axial static contrast', 'none', { hint: '0: axial static = kinetic (default). 1: the full static/kinetic difference applies axially' }],
      ['parameters.friction.weakening', 'Static to kinetic', 'select', { options: [['step', 'Step (Aarsnes & Shor)'], ['exponential', 'Exponential']] }],
    ] },
  { title: 'Damping',
    fields: [
      ['parameters.damping.fluid', 'Mud damping (Couette + Stokes layer)', 'bool'],
      ['parameters.damping.loss_factor', 'Steel loss factor η', 'none', { hint: 'Material damping ratio η/2 at the first mode and at the maximum frequency (Rayleigh). 0 = off' }],
      ['parameters.damping.yield_drag', 'Mud yield-stress drag', 'select', { options: [['pumps_off', 'Pumps off only'], ['always', 'Always'], ['off', 'Off']] }],
      ['parameters.damping.max_frequency', 'Maximum modelled frequency', 'hz', { hint: 'Model valid below this frequency (default 10 Hz)' }],
      ['parameters.damping.k_t', 'Extra torsional damping', 'per_s', { hint: 'Calibration term on top of the physical damping; default 0' }],
      ['parameters.damping.k_a', 'Extra axial damping', 'per_s', { hint: 'Calibration term on top of the physical damping; default 0' }],
    ] },
  { title: 'Bit',
    note: () => ({
      coulomb: 'Friction-dominated bit. Sliding: TOB = μk·WOB·D/3; sticking up to μs·WOB·D/3. ROP = ω·WOB / (2π·a·S).',
      cutting: 'Detournay & Defourny cutting. WOB = W_f + S·a·d, TOB = WOB·a/(2ζ), with d the depth of cut per revolution. No velocity weakening.',
      regenerative: 'Cutting with the depth of cut set by the hole bottom left by the previous blade. WOB = S·a·Σd_j, TOB = WOB·a/(2ζ).',
    })[get(state.cfg, 'parameters.bit.model') || 'coulomb'] + (get(state.cfg, 'parameters.bit.strength_source') === 'formation'
      ? ' Drilling strength S = ζ·CCS/η from the formation at the bit.' : ''),
    fields: [
      ['parameters.bit.model', 'Model', 'select',
        { options: [['coulomb', 'Coulomb (friction)'], ['cutting', 'Cutting (Detournay)'], ['regenerative', 'Regenerative cutting']], rebuild: true }],
      ['parameters.bit.strength_source', 'Rock strength', 'select',
        { options: [['formation', 'From the formation (CCS)'], ['constant', 'Fixed drilling strength']], rebuild: true }],
      ['parameters.bit.mu_static', 'Static friction μs', 'none', { models: ['coulomb'] }],
      ['parameters.bit.mu_kinetic', 'Kinetic friction μk', 'none', { models: ['coulomb'] }],
      ['parameters.bit.drilling_strength', 'Drilling strength S', 'stress', { source: 'constant', hint: 'S = ζ·ε, the axial cutting force per unit cut area' }],
      ['parameters.bit.zeta', 'Force ratio ζ', 'none', { models: ['cutting', 'regenerative'], orSource: 'formation', hint: 'Normal / tangential cutting force' }],
      ['parameters.bit.cutting_efficiency', 'Cutting efficiency η', 'none', { source: 'formation', hint: 'ε = CCS/η; 1 for a sharp cutter' }],
      ['parameters.bit.contact_wob', 'Contact WOB W_f', 'force', { models: ['cutting'], hint: 'WOB carried by wear flats before cutting starts' }],
      ['parameters.bit.blades', 'Blades', 'none', { models: ['regenerative'], orSource: 'formation' }],
      ['parameters.bit.nozzle_area', 'Nozzle flow area (TFA)', 'area', { hint: 'Sum of nozzle areas; one n/32 in nozzle is π(n/64)² in²' }],
      ['parameters.bit.discharge_coefficient', 'Nozzle discharge coefficient', 'none', { hint: '0.95 for jets' }],
    ] },
  { title: 'Hydraulics',
    note: () => 'Steady losses per Introduction to Drilling Engineering ch. 7 (Newtonian, or Bingham when a yield point is given). Transients (pump start, surge and swab) from the single-phase drift-flux model of Wang et al. (2023).',
    fields: [
      ['parameters.hydraulics.enabled', 'Model hydraulics (SPP, BHP, ECD)', 'bool'],
      ['parameters.hydraulics.back_pressure', 'Back pressure at the flowline', 'pressure'],
      ['parameters.hydraulics.surface_loss_coefficient', 'Surface loss c (Δp = c·Q^m, SI)', 'none', { hint: 'Pa per (m³/s)^m; 0 to neglect surface lines' }],
      ['parameters.hydraulics.surface_loss_exponent', 'Surface loss exponent m', 'none'],
      ['parameters.hydraulics.clinging_constant', 'Clinging constant (surge/swab)', 'none', { hint: 'Fraction of pipe speed that drags the annular mud' }],
      ['fluid.sound_speed', 'Speed of sound in mud', 'none', { hint: 'm/s; sets compressibility and the water-hammer speed' }],
      ['parameters.hydraulics.dx', 'Hydraulics cell length', 'length'],
    ] },
  { title: 'Ramp times',
    note: () => 'Inputs without a controller approach a new value linearly over their ramp time (0 = in one time step). Set points that feed a controller (RPM, autodriller WOB) are not ramped by default.',
    fields: [
      ['parameters.ramps.omega_td', 'Prescribed top-drive speed', 's'],
      ['parameters.ramps.v_block', 'Block speed', 's'],
      ['parameters.ramps.flow_rate', 'Pump rate', 's'],
      ['parameters.ramps.inlet_temperature', 'Inlet temperature', 's'],
      ['parameters.ramps.omega_sp', 'RPM set point (controlled; default none)', 's'],
      ['parameters.ramps.wob_sp', 'WOB set point (controlled; default none)', 's'],
    ] },
  { title: 'Formation', when: () => !!state.cfg.formation,
    note: () => 'Layers by top depth. CCS = UCS(T) + Δp·(1+sin φ)/(1−sin φ). Low-porosity rock (≤ 5 %) uses the undrained (Skempton) pore pressure; ≥ 20 % uses the drained one; linear blend between.',
    fields: [
      ['formation.pore_pressure_gradient', 'Pore pressure gradient', 'pgrad'],
      ['formation.overburden_gradient', 'Overburden gradient', 'pgrad'],
      ['formation.skempton_b', 'Skempton B', 'none'],
    ], table: 'layers' },
  { title: 'Rock temperature', when: () => !!state.cfg.formation,
    note: () => 'Thermal weakening above T*, and cooling by the mud at the hole bottom (thermal stress α·E·ΔT/(1−ν), scaled by κ).',
    fields: [
      ['formation.surface_temperature', 'Surface temperature', 'temp'],
      ['formation.geothermal_gradient', 'Geothermal gradient', 'tgrad'],
      ['formation.weakening_threshold', 'Weakening starts at T*', 'temp'],
      ['formation.weakening_rate', 'Weakening rate β', 'per_temp', { hint: 'Fractional UCS loss per degree above T*' }],
      ['formation.thermal_stress_efficiency', 'Thermal-stress efficiency κ', 'none'],
      ['formation.thermal_expansion', 'Thermal expansion α', 'per_temp'],
      ['formation.youngs_modulus', 'Young\'s modulus E', 'modulus'],
      ['formation.poisson_ratio', 'Poisson\'s ratio ν', 'none'],
      ['formation.rock_conductivity', 'Rock conductivity', 'cond'],
      ['formation.rock_heat_capacity', 'Rock heat capacity', 'cp'],
      ['formation.rock_density', 'Rock density', 'density'],
    ] },
  { title: 'Mud and circulation',
    note: () => 'Wellbore temperature: drift-flux model of Wang et al. (2023, SPE-214836-MS). Pump rate and inlet temperature can be changed live in the control bar.',
    fields: [
      ['parameters.thermal.enabled', 'Model wellbore temperature', 'bool'],
      ['parameters.thermal.flow_rate', 'Pump rate at start', 'flow'],
      ['parameters.thermal.inlet_temperature', 'Inlet temperature at start', 'temp'],
      ['fluid.viscosity', 'Mud viscosity', 'visc'],
      ['fluid.plastic_viscosity', 'Plastic viscosity (PV)', 'visc', { blank: 'use mud viscosity', hint: 'Bingham PV; used for the mud damping when given' }],
      ['fluid.yield_point', 'Yield point (YP)', 'yp', { blank: 'none', hint: 'Bingham YP; gives a yield-stress drag on the string (see Damping)' }],
      ['fluid.heat_capacity', 'Mud heat capacity', 'cp'],
      ['fluid.thermal_conductivity', 'Mud conductivity', 'cond'],
      ['fluid.thermal_expansion', 'Mud thermal expansion', 'per_temp'],
      ['parameters.thermal.dx', 'Thermal cell length', 'length'],
    ] },
  { title: 'Top drive and controller',
    fields: [
      ['parameters.controller.type', 'Controller', 'select',
        { options: [['stiff', 'Stiff PI'], ['soft_torque', 'SoftTorque'], ['z_torque', 'ZTorque'], ['custom', 'Custom PI']] }],
      ['parameters.controller.f_c', 'SoftTorque tuning frequency', 'hz', { blank: 'auto' }],
      ['parameters.controller.t_hp', 'ZTorque high-pass', 's'],
      ['parameters.controller.kp', 'Override Kp', 'kp', { blank: 'auto' }],
      ['parameters.controller.ki', 'Override Ki', 'ki', { blank: 'auto' }],
      ['parameters.top_drive.inertia', 'Top drive inertia', 'inertia'],
      ['parameters.top_drive.max_torque', 'Torque limit', 'torque', { hint: '0 for no limit' }],
      ['parameters.drawworks.block_acceleration', 'Block acceleration limit', 'accel', { hint: '0 for instant changes in block speed' }],
    ] },
  { title: 'Numerics',
    fields: [
      ['simulation.dx', 'Cell size', 'length', { hint: 'Halving the cell size makes runs about 4× slower' }],
      ['simulation.record_rate', 'Recording rate', 'hz'],
      ['parameters.cfl', 'Courant number', 'none'],
    ] },
];

// [key, column label, quantity kind ('text' for strings)]
const SURVEY_COLS = [['md', 'MD', 'depth'], ['inc', 'Inc', 'deg'], ['azi', 'Azi', 'deg']];
const LAYER_COLS = [['name', 'Name', 'text'], ['top_md', 'Top MD', 'depth'], ['ucs', 'UCS', 'stress'],
                    ['friction_angle', 'φ', 'deg'], ['porosity', 'Porosity', 'none']];
const STRING_COLS = [['name', 'Name', 'text'], ['length', 'Length', 'length'], ['od', 'OD', 'diameter'],
                     ['id', 'ID', 'diameter'], ['linear_mass', 'Mass', 'linear_mass'], ['contact_od', 'Contact OD', 'diameter']];

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
function sig(x, n = 6) {
  // up to n significant figures, no trailing zeros
  if (x == null || x === '' || Number.isNaN(x)) return '';
  if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e7)) return Number(x).toExponential(3);
  return String(Number(Number(x).toPrecision(n)));
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
  const [path, label, kind, opt = {}] = spec;
  const row = document.createElement('div');
  const id = 'f-' + path.replace(/\./g, '-');
  if (kind === 'bool') {
    row.className = 'field check';
    row.innerHTML = `<input type="checkbox" id="${id}"><label for="${id}">${label}</label>`;
    const el = row.querySelector('input');
    el.checked = !!get(state.cfg, path);
    el.addEventListener('change', () => { set(state.cfg, path, el.checked); markDirty(); });
    return row;
  }
  row.className = 'field';
  if (kind === 'select') {
    row.innerHTML = `<label for="${id}">${label}</label><select id="${id}"></select><span></span>`;
    const el = row.querySelector('select');
    for (const [v, t] of opt.options) el.add(new Option(t, v));
    el.value = get(state.cfg, path);
    el.addEventListener('change', () => { set(state.cfg, path, el.value); markDirty(); if (opt.rebuild) buildForm(); });
    return row;
  }
  row.innerHTML = `<label for="${id}">${label}</label><input type="number" id="${id}" step="any"><span class="unit">${unitOf(kind)}</span>`;
  const el = row.querySelector('input');
  const v = get(state.cfg, path);
  el.value = v == null ? '' : sig(toDisp(kind, v));
  if (opt.blank) el.placeholder = opt.blank;
  if (opt.hint) el.title = opt.hint;
  el.addEventListener('change', () => {
    set(state.cfg, path, el.value === '' ? null : fromDisp(kind, Number(el.value))); markDirty();
  });
  return row;
}

function tableEditor(kind) {
  const cols = { survey: SURVEY_COLS, string: STRING_COLS, layers: LAYER_COLS }[kind];
  const rows = () => ({ survey: state.cfg.well.survey, string: state.cfg.drillstring, layers: state.cfg.formation.layers }[kind]);
  const wrap = document.createElement('div');
  const draw = () => {
    wrap.innerHTML = '';
    const t = document.createElement('table'); t.className = 'grid' + (kind === 'string' ? ' wide' : '');
    t.innerHTML = '<thead><tr>' + cols.map(([, label, kind]) => `<th>${label}${kind !== 'text' ? ` <span>(${unitOf(kind)})</span>` : ''}</th>`).join('') + '<th></th></tr></thead>';
    const tb = document.createElement('tbody');
    rows().forEach((r, i) => {
      const tr = document.createElement('tr');
      for (const [key, , kind] of cols) {
        const td = document.createElement('td');
        const inp = document.createElement('input');
        if (kind === 'text') { inp.type = 'text'; inp.value = r[key] ?? ''; }
        else { inp.type = 'number'; inp.step = 'any'; inp.value = r[key] == null ? '' : sig(toDisp(kind, r[key]), 5); }
        inp.setAttribute('aria-label', `${key} row ${i + 1}`);
        inp.addEventListener('change', () => {
          r[key] = kind === 'text' ? inp.value : (inp.value === '' ? null : fromDisp(kind, Number(inp.value)));
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
      rows().push(kind === 'survey' ? { md: (last.md || 0) + fromDisp('depth', 100), inc: last.inc || 0, azi: last.azi || 0 }
                : kind === 'layers' ? { ...last, name: 'New layer', top_md: (last.top_md || 0) + fromDisp('depth', 100) }
                : { ...last, name: 'New component' });
      markDirty(); draw();
    });
    const imp = document.createElement('label'); imp.className = 'btn quiet file';
    imp.textContent = 'Import CSV';
    const fi = document.createElement('input'); fi.type = 'file'; fi.accept = '.csv,text/csv';
    fi.addEventListener('change', async () => {
      try {
        const parsed = parseCSV(await fi.files[0].text(), cols);
        if (kind === 'survey') state.cfg.well.survey = parsed;
        else if (kind === 'layers') state.cfg.formation.layers = parsed;
        else state.cfg.drillstring = parsed;
        markDirty(); draw(); message(`Imported ${parsed.length} rows (read as ${SYSTEM === 'field' ? 'field' : 'SI'} units).`, 'info');
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
  // Values are read in the units currently shown in the table headers.
  const lines = text.trim().split(/\r?\n/).filter(l => l.trim());
  const head = lines.shift().split(',').map(h => h.trim().toLowerCase().replace(/\s*\(.*\)\s*$/, ''));
  const idx = cols.map(([key]) => head.indexOf(key));
  const optional = ['linear_mass', 'contact_od', 'azi', 'name', 'friction_angle', 'porosity'];
  const need = cols.filter(([key], k) => idx[k] < 0 && !optional.includes(key));
  if (need.length) throw new Error(`The CSV is missing column(s): ${need.map(c => c[0]).join(', ')}. Expected headers: ${cols.map(([k, , kind]) => (kind === 'text' ? k : `${k} (${unitOf(kind)})`)).join(', ')}.`);
  return lines.map(l => {
    const v = l.split(',').map(x => x.trim()); const r = {};
    cols.forEach(([key, , kind], k) => {
      const x = idx[k] >= 0 ? v[idx[k]] : '';
      r[key] = kind === 'text' ? x : (x === '' ? (key === 'azi' ? 0 : null) : fromDisp(kind, Number(x)));
    });
    return r;
  });
}

function buildForm() {
  const form = $('form');
  const wasOpen = [...form.querySelectorAll('details.sec')].map(d => d.open);
  form.innerHTML = '';
  for (const [k, sec] of SECTIONS.entries()) {
    if (sec.when && !sec.when()) continue;
    const d = document.createElement('details'); d.className = 'sec'; d.open = wasOpen.length ? wasOpen[k] : !!sec.open;
    d.innerHTML = `<summary>${sec.title}</summary>`;
    const body = document.createElement('div'); body.className = 'sec-body';
    if (sec.note) { const p = document.createElement('p'); p.className = 'sec-note'; p.textContent = sec.note(); body.appendChild(p); }
    const bitModel = get(state.cfg, 'parameters.bit.model') || 'coulomb';
    const source = get(state.cfg, 'parameters.bit.strength_source') || 'constant';
    for (const f of sec.fields) {
      const opt = f[3] || {};
      const byModel = !opt.models || opt.models.includes(bitModel);
      const bySource = opt.orSource ? (byModel || opt.orSource === source) : byModel && (!opt.source || opt.source === source);
      if (!bySource) continue;
      body.appendChild(fieldRow(f));
    }
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
  return axis(show ? `Measured depth (${unitOf('depth')})` : '', { range: r && srv('depth', r), autorange: false, showticklabels: show });
}

function drawStringViews(fr) {
  const p = fr.profile;
  const ref = state.reference;
  const md = srv('depth', p.md), refMd = srv('depth', ref.md);
  const uD = unitOf('depth'), uT = unitOf('torque'), uF = unitOf('force');
  // ---- RPM along the string
  state.rpmMax = Math.max(state.rpmMax, ...p.rpm.map(Math.abs));
  Plotly.react('p-rpm', [
    { x: p.rpm, y: md, mode: 'lines', line: { color: COLORS.rot, width: 2 }, hovertemplate: `%{x:.0f} rpm at %{y:.0f} ${uD}<extra></extra>` },
  ], baseLayout({ title: TITLE('Rotary speed'), margin: { l: 56, r: 14, t: 74, b: 40 },
    xaxis: axis('rpm', { range: [-0.08 * state.rpmMax, state.rpmMax * 1.05] }), yaxis: depthAxis(true) }), PCONF);

  // ---- torque along the string with T&D references (range tracked in SI)
  state.tqMax = Math.max(state.tqMax, ...p.torque_knm, ...ref.torque_breakoff_knm);
  state.tqMin = Math.min(state.tqMin, ...p.torque_knm);
  Plotly.react('p-torque', [
    { x: srv('torque', ref.torque_breakoff_knm), y: refMd, mode: 'lines', line: { color: COLORS.ref, dash: 'dot', width: 1.5 }, name: 'Static μ', hoverinfo: 'skip' },
    { x: srv('torque', ref.torque_rotating_knm), y: refMd, mode: 'lines', line: { color: COLORS.ref, dash: 'dash', width: 1.5 }, name: 'Kinetic μ', hoverinfo: 'skip' },
    { x: srv('torque', p.torque_knm), y: md, mode: 'lines', line: { color: COLORS.trq, width: 2 }, name: 'Simulated', hovertemplate: `%{x:.2f} ${uT} at %{y:.0f} ${uD}<extra></extra>` },
  ], baseLayout({ title: TITLE('Torque'), showlegend: true,
    legend: LEGEND_TOP, margin: { l: 14, r: 14, t: 74, b: 40 },
    xaxis: axis(uT, { range: srv('torque', [Math.min(0, state.tqMin) * 1.1, state.tqMax * 1.08]) }), yaxis: depthAxis(false) }), PCONF);

  // ---- axial force
  Plotly.react('p-tension', [
    { x: srv('force', ref.tension_hanging_kn), y: refMd, mode: 'lines', line: { color: COLORS.ref, dash: 'dash', width: 1.5 }, name: 'Hanging', hoverinfo: 'skip' },
    { x: srv('force', p.tension_kn), y: md, mode: 'lines', line: { color: COLORS.ax, width: 2 }, name: 'Simulated', hovertemplate: `%{x:.0f} ${uF} at %{y:.0f} ${uD}<extra></extra>` },
  ], baseLayout({ title: TITLE('Effective tension'), showlegend: true,
    legend: LEGEND_TOP, margin: { l: 14, r: 14, t: 74, b: 40 },
    xaxis: axis(uF), yaxis: depthAxis(false) }), PCONF);

  // ---- wellbore temperature (Wang et al. 2023 drift-flux model)
  const th = fr.thermal;
  if (th) {
    const tmd = srv('depth', th.md), uT = unitOf('temp');
    Plotly.react('p-temp', [
      { x: srv('temp', th.t_far), y: tmd, mode: 'lines', line: { color: COLORS.ref, dash: 'dash', width: 1.5 }, name: 'Formation', hoverinfo: 'skip' },
      { x: srv('temp', th.t_ds), y: tmd, mode: 'lines', line: { color: COLORS.rot, width: 2 }, name: 'Drillstring', hovertemplate: `%{x:.0f} ${uT} at %{y:.0f} ${uD}<extra></extra>` },
      { x: srv('temp', th.t_ann), y: tmd, mode: 'lines', line: { color: '#a3402a', width: 2 }, name: 'Annulus', hovertemplate: `%{x:.0f} ${uT} at %{y:.0f} ${uD}<extra></extra>` },
    ], baseLayout({ title: TITLE('Temperature'), showlegend: true, legend: LEGEND_TOP, margin: { l: 14, r: 14, t: 74, b: 40 },
      xaxis: axis(uT), yaxis: depthAxis(false) }), PCONF);
  } else {
    Plotly.react('p-temp', [], baseLayout({ title: TITLE('Temperature'), xaxis: axis(''), yaxis: depthAxis(false),
      annotations: [{ text: 'Wellbore temperature model is off', showarrow: false, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5 }] }), PCONF);
  }

  // ---- well path with the drill pipe, BHA and bit positions
  const path = state.path;
  const cased = path.md.map((m, i) => (m <= path.shoe_md ? i : -1)).filter(i => i >= 0);
  const pvs = srv('length', path.vs), ptvd = srv('depth', path.tvd);
  const vs = srv('length', p.vs), tvd = srv('depth', p.tvd);
  const pick = (arr, want) => arr.map((v, i) => (path.bha[i] === want ? v : null));
  const nb = p.md.length - 1;
  const bitTip = `Bit at ${toDisp('depth', fr.latest.bit_depth).toFixed(1)} ${uD} MD`;
  Plotly.react('p-path', [
    { x: pvs, y: ptvd, mode: 'lines', line: { color: '#d9e0e7', width: 16 }, hoverinfo: 'skip', name: 'Open hole' },
    { x: cased.map(i => pvs[i]), y: cased.map(i => ptvd[i]), mode: 'lines', line: { color: '#a9b5c1', width: 18 },
      hoverinfo: 'skip', name: `Casing to ${toDisp('depth', path.shoe_md).toFixed(0)} ${uD}` },
    { x: pick(vs, 0), y: pick(tvd, 0), mode: 'lines', line: { color: COLORS.steel, width: 3 }, name: 'Drill pipe',
      hoverinfo: 'skip' },
    { x: pick(vs, 1), y: pick(tvd, 1), mode: 'lines', line: { color: COLORS.ink, width: 8 },
      name: `BHA (${toDisp('length', path.bha_length).toFixed(0)} ${unitOf('length')})`, hoverinfo: 'skip' },
    { x: [vs[nb]], y: [tvd[nb]], mode: 'markers', name: 'Bit', hovertemplate: bitTip + '<extra></extra>',
      marker: { symbol: 'diamond', size: 13, color: COLORS.trq, line: { color: COLORS.ink, width: 1 } } },
  ], baseLayout({ title: TITLE('Well path'), showlegend: true,
    legend: { x: 0.98, xanchor: 'right', y: 0.98, yanchor: 'top', bgcolor: 'rgba(255,255,255,0.85)', font: { size: 11 } },
    xaxis: axis(`Horizontal displacement (${unitOf('length')})`, { scaleanchor: 'y', scaleratio: 1 }),
    yaxis: axis(`True vertical depth (${uD})`, { autorange: 'reversed' }), margin: { l: 64, r: 14, t: 40, b: 40 } }), PCONF);
}

function emptyTraces() {
  return { time: [], rpm_td: [], rpm_bit: [], rpm_sp: [], motor_torque_knm: [], pipe_torque_knm: [], tob_knm: [],
           hookload_kn: [], wob_kn: [], rop_mh: [], block_mh: [], bit_depth: [], hole_depth: [] };
}
function appendTrace(t) { for (const k in state.traces) { const a = t[k]; if (a) for (let i = 0; i < a.length; i++) state.traces[k].push(a[i]); } }
function emptyHyd() { return { time: [], spp: [], p_bh: [], ecd: [], q_in: [], q_out: [] }; }
function appendHyd(h) { if (!h) return; for (const k in state.hyd) { const a = h[k]; if (a) for (let i = 0; i < a.length; i++) state.hyd[k].push(a[i]); } }
function emptySlow() { return { time: [], t_bit: [], t_out: [], t_in: [], t_rock: [], ccs_eff: [], ccs: [], S: [], ecd: [], p_bh: [] }; }
function appendSlow(sl) { if (!sl) return; for (const k in state.slow) { const a = sl[k]; if (a) for (let i = 0; i < a.length; i++) state.slow[k].push(a[i]); } }

function timeRange(tNow) {
  const W = state.windowS;
  if (!W) return [0, Math.max(tNow, 1)];
  return tNow <= W ? [0, W] : [tNow - W, tNow];
}

// Time history as a drilling-log (EDR) display: tracks side by side, a shared
// vertical time axis increasing downward (newest data at the bottom), and a
// header per curve with its name, current value and scale ends.
const LOG_GAP = 0.014;          // gap between tracks (fraction of plot width)
const LOG_ROW = 30;             // header height per curve (px)
const LOG_HEIGHT = 780;         // figure height (px)
function logTracks() {
  const T = state.traces, L = state.slow, H = state.hyd;
  const dsp = (k, a) => a.map(v => toDisp(k, v));
  const f = (y, name, color, unit, extra = {}) => Object.assign({ t: T.time, y, name, color, unit }, extra);
  const s = (y, name, color, unit, extra = {}) => Object.assign({ t: L.time, y, name, color, unit }, extra);
  const h = (y, name, color, unit, extra = {}) => Object.assign({ t: H.time, y, name, color, unit }, extra);
  const uT = unitOf('torque'), uF = unitOf('force'), uP = unitOf('pressure'), uE = unitOf('ecd'), uTe = unitOf('temp'), uS = unitOf('stress');
  return [
    { title: 'Rotary speed', curves: [
      f(T.rpm_sp, 'RPM set point', COLORS.ref, 'rpm', { id: 'rpm_sp', dash: 'dash' }),
      f(T.rpm_td, 'Top drive', COLORS.ink, 'rpm', { id: 'rpm_td' }),
      f(T.rpm_bit, 'Bit', COLORS.rot, 'rpm', { id: 'rpm_bit' })] },
    { title: 'Torque', curves: [
      f(srv('torque', T.motor_torque_knm), 'Motor', COLORS.ink, uT, { id: 'tq_motor' }),
      f(srv('torque', T.pipe_torque_knm), 'Surface pipe', COLORS.trq, uT, { id: 'tq_pipe' }),
      f(srv('torque', T.tob_knm), 'On bit', COLORS.bit, uT, { id: 'tq_bit' })] },
    { title: 'Axial', curves: [
      f(srv('force', T.hookload_kn), 'Hookload', COLORS.ax, uF, { id: 'hookload' }),
      f(srv('force', T.wob_kn), 'WOB', COLORS.ink, uF, { id: 'wob', zero: true })] },
    { title: 'Hydraulics', curves: [
      h(dsp('pressure', H.spp), 'SPP', COLORS.ink, uP, { id: 'spp' }),
      h(dsp('ecd', H.ecd), 'ECD', COLORS.ax, uE, { id: 'ecd' })] },
    { title: 'Bit and rock', curves: [
      s(srv('temp', L.t_rock), 'Rock at bit', COLORS.ref, uTe, { id: 't_rock', dash: 'dash' }),
      s(srv('temp', L.t_bit), 'Mud at bit', '#a3402a', uTe, { id: 't_bit' }),
      s(srv('stress', L.ccs_eff), 'Effective CCS', COLORS.ink, uS, { id: 'ccs' })] },
  ];
}

function niceRange(lo, hi, zero) {
  if (!isFinite(lo) || !isFinite(hi)) return [0, 1];
  if (zero) lo = Math.min(lo, 0);
  if (hi - lo < 1e-9) { const d = Math.max(Math.abs(hi) * 0.05, 1); lo -= d; hi += d; }
  const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / 4)));
  const m = [1, 2, 2.5, 5, 10].find(k => (hi - lo) / (k * step) <= 5) * step;
  return [Math.floor(lo / m) * m, Math.ceil(hi / m) * m];
}
function fmtScale(v) { const a = Math.abs(v); return a >= 1000 ? v.toFixed(0) : a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1).replace(/\.0$/, '') : v.toPrecision(2).replace(/\.?0+$/, ''); }

function drawTraces(cursor) {
  const T = state.traces; const tNow = T.time.length ? T.time[T.time.length - 1] : 0;
  const [t0, t1] = timeRange(tNow);
  const tracks = logTracks();
  const n = tracks.length, wTrack = (1 - (n - 1) * LOG_GAP) / n;
  const rows = Math.max(...tracks.map(k => k.curves.length));
  const el = $('t-log');
  const height = LOG_HEIGHT, top = 30 + rows * LOG_ROW, bottom = 14;
  const px = (p) => p / (height - top - bottom);          // pixels -> paper units above the plot
  const data = [], shapes = [], annotations = [], hits = [];
  const layout = baseLayout({ height, margin: { l: 64, r: 12, t: top, b: bottom }, showlegend: false, hovermode: 'y unified',
    datarevision: `${T.time.length}-${SYSTEM}`,
    yaxis: axis('', { range: [t1, t0], autorange: false, title: { text: 'Time (s)', standoff: 4 },
                      gridcolor: '#dfe5eb', zeroline: false }) });
  let ax = 0;
  tracks.forEach((tr, i) => {
    const d0 = i * (wTrack + LOG_GAP), d1 = d0 + wTrack;
    // track frame and title
    shapes.push({ type: 'rect', xref: 'paper', yref: 'paper', x0: d0, x1: d1, y0: 0, y1: 1, line: { color: '#b9c3ce', width: 1 }, layer: 'below' });
    annotations.push({ xref: 'paper', yref: 'paper', x: (d0 + d1) / 2, y: 1 + px(top - 4), xanchor: 'center', yanchor: 'top', showarrow: false,
                       text: `<b>${tr.title}</b>`, font: { size: 13, color: COLORS.ink } });
    // curves fill the header rows from the plot upward, so the last curve sits next to the plot
    const off = rows - tr.curves.length;
    let base = null;
    tr.curves.forEach((c, j) => {
      ax += 1; const xa = ax === 1 ? 'x' : 'x' + ax, key = ax === 1 ? 'xaxis' : 'xaxis' + ax;
      // scale from the data inside the time window
      let lo = Infinity, hi = -Infinity, last = null;
      for (let k = 0; k < c.t.length; k++) { const tk = c.t[k], v = c.y[k]; if (v == null || !isFinite(v)) continue;
        if (tk >= t0 && tk <= t1) { if (v < lo) lo = v; if (v > hi) hi = v; } last = v; }
      const auto = niceRange(lo, hi, c.zero);
      const man = manualScale(c.id);
      const [r0, r1] = man || auto;
      state.autoScale[c.id] = auto;
      layout[key] = { domain: [d0, d1], range: [r0, r1], anchor: 'y', showticklabels: false, ticks: '', zeroline: false,
                      showgrid: j === 0, gridcolor: '#eef1f4', nticks: 5, fixedrange: true };
      if (base) layout[key].overlaying = base; else base = xa;
      data.push({ x: c.y, y: c.t, xaxis: xa, yaxis: 'y', type: c.t === T.time ? 'scattergl' : 'scatter', mode: 'lines', name: c.name,
                  line: { color: c.color, width: 1.5, dash: c.dash || 'solid' },
                  hovertemplate: `${c.name}: %{x:.4g} ${c.unit}<extra></extra>` });
      // curve header row (Pason style): scale minimum | name and current value | scale maximum, coloured bar below
      const rowTop = top - 26 - (off + j) * LOG_ROW;              // px above the plot area
      hits.push({ id: c.id, name: c.name, unit: c.unit, color: c.color, man: !!man, d0, d1, y: top - rowTop - 2, h: 26 });
      const yt = 1 + px(rowTop - 9), yb = 1 + px(rowTop - 21);
      const pad = 0.005;
      annotations.push({ xref: 'paper', yref: 'paper', x: (d0 + d1) / 2, y: yt, xanchor: 'center', yanchor: 'middle', showarrow: false,
                         text: `${man ? '\u{1F512} ' : ''}<b>${c.name}</b> ${last == null ? '–' : fmtScale(last)} ${c.unit}`,
                         font: { size: 11.5, color: c.color } });
      annotations.push({ xref: 'paper', yref: 'paper', x: d0 + pad, y: yt, xanchor: 'left', yanchor: 'middle', showarrow: false,
                         text: man ? `<b>${fmtScale(r0)}</b>` : fmtScale(r0), font: { size: 10.5, color: c.color } });
      annotations.push({ xref: 'paper', yref: 'paper', x: d1 - pad, y: yt, xanchor: 'right', yanchor: 'middle', showarrow: false,
                         text: man ? `<b>${fmtScale(r1)}</b>` : fmtScale(r1), font: { size: 10.5, color: c.color } });
      shapes.push({ type: 'line', xref: 'paper', yref: 'paper', x0: d0 + pad, x1: d1 - pad, y0: yb, y1: yb,
                    line: { color: c.color, width: 2.5, dash: c.dash || 'solid' } });
    });
  });
  if (cursor != null) shapes.push({ type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: cursor, y1: cursor,
                                    line: { color: COLORS.ink, width: 1, dash: 'dot' } });
  layout.shapes = shapes; layout.annotations = annotations;
  Plotly.react('t-log', data, layout, PCONF);
  placeScaleHits(hits, layout.margin);
}

// ------------------------------------------------- manual curve scales
// Limits are kept per unit system, in display units, and remembered in this browser.
function loadScales() {
  try { const s = JSON.parse(localStorage.getItem('drillpyze-scales') || '{}'); return { si: s.si || {}, field: s.field || {} }; }
  catch (e) { return { si: {}, field: {} }; }
}
function saveScales() { try { localStorage.setItem('drillpyze-scales', JSON.stringify(state.scales)); } catch (e) { /* storage unavailable */ } }
function manualScale(id) {
  const m = state.scales[SYSTEM] && state.scales[SYSTEM][id];
  return m && isFinite(m[0]) && isFinite(m[1]) && m[1] > m[0] ? m : null;
}
function setScale(id, lim) {
  if (lim) state.scales[SYSTEM][id] = lim; else delete state.scales[SYSTEM][id];
  saveScales(); drawTraces(state.scrub == null ? null : state.frames[state.scrub].time);
}
// Invisible buttons over each curve header; clicking one opens a small pop-up
// with that curve's limits only.
function placeScaleHits(hits, m) {
  const el = $('t-log'), layer = $('scale-hits');
  const W = el.clientWidth - m.l - m.r;
  const sig = JSON.stringify([W, SYSTEM, hits.map(h => [h.id, h.man, h.y])]);
  if (layer.dataset.sig !== sig) {
    layer.dataset.sig = sig; layer.innerHTML = '';
    for (const h of hits) {
      const btn = document.createElement('button');
      btn.type = 'button'; btn.className = 'scale-hit'; btn.dataset.id = h.id;
      btn.style.left = `${m.l + h.d0 * W}px`; btn.style.width = `${(h.d1 - h.d0) * W}px`;
      btn.style.top = `${h.y}px`; btn.style.height = `${h.h}px`;
      btn.title = `${h.name}: ${h.man ? 'fixed scale' : 'automatic scale'} (click to set limits)`;
      btn.setAttribute('aria-label', `Set scale limits for ${h.name} (${h.unit})`);
      btn.addEventListener('click', (ev) => { ev.stopPropagation(); openScalePop(h, btn); });
      layer.appendChild(btn);
    }
  }
  $('scales-reset').hidden = !Object.keys(state.scales[SYSTEM] || {}).length;
}

function openScalePop(h, anchor) {
  const pop = $('scale-pop'), auto = state.autoScale[h.id] || [0, 1], man = manualScale(h.id);
  pop.dataset.id = h.id;
  $('sp-title').innerHTML = `<span class="swatch" style="border-color:${h.color}"></span>${h.name}`;
  $('sp-unit-lo').textContent = h.unit; $('sp-unit-hi').textContent = h.unit;
  $('sp-state').textContent = man ? 'Fixed scale' : 'Automatic scale';
  $('sp-lo').value = (man || auto)[0]; $('sp-hi').value = (man || auto)[1];
  $('sp-msg').textContent = '';
  const box = $('t-log'), a = anchor.getBoundingClientRect(), b = box.getBoundingClientRect();
  pop.hidden = false;
  const w = pop.offsetWidth;
  pop.style.left = `${Math.min(Math.max(a.left - b.left + a.width / 2 - w / 2, 4), b.width - w - 4)}px`;
  pop.style.top = `${a.bottom - b.top + 6}px`;
  $('sp-lo').focus(); $('sp-lo').select();
}
function closeScalePop() { $('scale-pop').hidden = true; }
function applyScalePop() {
  const id = $('scale-pop').dataset.id, a = parseFloat($('sp-lo').value), b = parseFloat($('sp-hi').value);
  if (!(isFinite(a) && isFinite(b) && b > a)) { $('sp-msg').textContent = 'Max must be greater than min.'; return; }
  setScale(id, [a, b]); closeScalePop();
}

function showReadouts(L) {
  $('ro-rpm-td').textContent = fixed(L.rpm_td, 0);
  $('ro-rpm-bit').textContent = fixed(L.rpm_bit, 0);
  $('ro-tq').textContent = fixed(srv('torque', L.pipe_torque_knm), 2);
  $('ro-tob').textContent = fixed(srv('torque', L.tob_knm), 2);
  $('ro-hk').textContent = fixed(srv('force', L.hookload_kn), 0);
  $('ro-wob').textContent = fixed(srv('force', L.wob_kn), 1);
  $('ro-rop').textContent = fixed(srv('rop', L.rop_mh), 1);
  $('ro-depth').textContent = fixed(toDisp('depth', L.bit_depth), 2);
  $('ro-tbit').textContent = fixed(srv('temp', L.t_bit), 0);
  $('ro-ccs').textContent = fixed(srv('stress', L.ccs_eff), 0);
  $('ro-spp').textContent = fixed(toDisp('pressure', L.spp), SYSTEM === 'field' ? 0 : 2);
  $('ro-ecd').textContent = fixed(toDisp('ecd', L.ecd_hyd), SYSTEM === 'field' ? 2 : 0);
  for (const el of document.querySelectorAll('[data-kind]')) el.textContent = unitOf(el.dataset.kind);
}

function render(fr, live = true) {
  drawStringViews(fr);
  showReadouts(fr.latest);
  $('clock').textContent = `${fr.time.toFixed(2)} s`;
  drawTraces(live ? null : fr.time);
}

// -------------------------------------------------------------- simulation
function inputs() {
  // the server takes block speed in m/h and WOB in kN
  const dw = Number($('dw').value) || 0;
  const vel = $('dw-mode').value === 'velocity';
  const flow = $('flow').value === '' ? null : fromDisp('flow', Number($('flow').value));
  const tin = $('t-in').value === '' ? null : fromDisp('temp', Number($('t-in').value));
  return { td_mode: $('td-mode').value, rpm: Number($('rpm').value) || 0,
           flow_m3s: flow, inlet_c: tin,
           dw_mode: $('dw-mode').value,
           v_block_mh: vel ? fromDisp('rop', dw) * 3600 : 0,
           wob_kn: vel ? 0 : fromDisp('force', dw) / 1e3 };
}

async function initialise() {
  pause();
  message(null);
  $('apply').disabled = true; $('apply').textContent = 'Starting…';
  try {
    const r = await api('/api/init', { config: state.cfg, previous: state.session });
    Object.assign(state, { session: r.session, info: r.info, reference: r.reference, path: r.path,
                           frames: [], traces: emptyTraces(), slow: emptySlow(), hyd: emptyHyd(), scrub: null, rpmMax: Math.max(60, Number($('rpm').value) * 1.5),
                           tqMax: 1, tqMin: 0, dirty: false });
    $('dirty-hint').hidden = true;
    setCirculation(get(state.cfg, 'parameters.thermal.flow_rate'), get(state.cfg, 'parameters.thermal.inlet_temperature'));
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
  appendSlow(fr.slow);
  appendHyd(fr.hyd);
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
function dwLimits() {
  // slider ranges chosen in SI and expressed in display units
  return $('dw-mode').value === 'wob'
    ? { kind: 'force', min: 0, max: toDisp('force', 300e3) }
    : { kind: 'rop', min: toDisp('rop', -60 / 3600), max: toDisp('rop', 60 / 3600) };
}
function applyDwUnits(valueSI) {
  const L = dwLimits(), s = $('dw-slider');
  const step = (L.max - L.min) / 240;
  s.min = L.min; s.max = L.max; s.step = step;
  $('dw-unit').textContent = unitOf(L.kind);
  const v = valueSI == null ? 0 : toDisp(L.kind, valueSI);
  $('dw').value = fmt(v, 1); $('dw').step = Number(step.toPrecision(1)); s.value = v;
}
function dwModeChanged() { applyDwUnits(0); }
function setCirculation(flowSI, tinSI) {
  if (flowSI != null) $('flow').value = sig(toDisp('flow', flowSI), 4);
  if (tinSI != null) $('t-in').value = sig(toDisp('temp', tinSI), 4);
}
function circulationSI() {
  return [$('flow').value === '' ? null : fromDisp('flow', Number($('flow').value)),
          $('t-in').value === '' ? null : fromDisp('temp', Number($('t-in').value))];
}

// ------------------------------------------------------------ unit system
function setSystem(sys) {
  if (sys === SYSTEM) return;
  // keep the live drawworks input physically the same across the switch
  const L = dwLimits();
  const dwSI = fromDisp(L.kind, Number($('dw').value) || 0);
  const circ = circulationSI();
  SYSTEM = sys;
  try { localStorage.setItem('drillpyze-units', sys); } catch (e) { /* storage unavailable */ }
  closeScalePop();
  for (const b of document.querySelectorAll('#units button')) b.setAttribute('aria-checked', String(b.dataset.sys === sys));
  $('dw-mode').options[0].text = `Block speed (+ lowering)`;
  applyDwUnits(dwSI);
  setCirculation(circ[0], circ[1]);
  for (const el of document.querySelectorAll('[data-kind]')) el.textContent = unitOf(el.dataset.kind);
  if (state.cfg) buildForm();
  const fr = state.scrub == null ? state.frames[state.frames.length - 1] : state.frames[state.scrub];
  if (fr) render(fr, state.scrub == null);
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
  for (const b of document.querySelectorAll('#units button')) b.addEventListener('click', () => setSystem(b.dataset.sys));

  $('scrub').addEventListener('input', (e) => scrubTo(Number(e.target.value)));
  $('go-live').addEventListener('click', goLive);

  $('scales-reset').addEventListener('click', () => { state.scales[SYSTEM] = {}; saveScales(); closeScalePop();
    drawTraces(state.scrub == null ? null : state.frames[state.scrub].time); });
  $('sp-apply').addEventListener('click', applyScalePop);
  $('sp-auto').addEventListener('click', () => { setScale($('scale-pop').dataset.id, null); closeScalePop(); });
  $('sp-close').addEventListener('click', closeScalePop);
  $('scale-pop').addEventListener('keydown', (ev) => {
    if (ev.key === 'Enter') { ev.preventDefault(); applyScalePop(); }
    if (ev.key === 'Escape') { ev.preventDefault(); closeScalePop(); }
  });
  document.addEventListener('click', (ev) => { const p = $('scale-pop');
    if (!p.hidden && !p.contains(ev.target) && !ev.target.classList.contains('scale-hit')) closeScalePop(); });
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
  for (const b of document.querySelectorAll('#units button')) b.setAttribute('aria-checked', String(b.dataset.sys === SYSTEM));
  applyDwUnits(0);
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
