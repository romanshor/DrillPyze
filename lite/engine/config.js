// Case description, defaults, validation and builder (spec sections 1, 4.1 and 7)
import { Survey } from './survey.js';
import { Drillstring } from './drillstring.js';
import { Simulation, stiffController } from './simulation.js';
import { FORMATION_DEFAULTS } from './formation.js';
import { BIT_MODELS } from './kernel.js';
import { IN, PSI } from './units.js';

export const FLUID_DEFAULTS = { density: 1000, plastic_viscosity: null, yield_point: null, viscosity: 0.020, heat_capacity: 3750,
                                thermal_conductivity: 0.75, thermal_expansion: 1e-4, sound_speed: 1500, reference_temperature: 20 };
const DEFAULTS = {
  friction: { mu_static_cased: 0.30, mu_kinetic_cased: 0.20, mu_static_open: 0.35, mu_kinetic_open: 0.25, omega_c: 1.0, v_c: 0.01,
              along_string_coupling: true, axial_static_contrast: 0.0, weakening: 'step' },
  damping: { k_t: 0, k_a: 0, fluid: true, loss_factor: 0.01, reference_frequency_axial: null, reference_frequency_torsional: null,
             high_pass_fraction: 0.25, yield_drag: 'pumps_off', max_frequency: 10 },
  bit: { model: 'coulomb', mu_static: 0.6, mu_kinetic: 0.4, drilling_strength: 50e6, zeta: 0.8, contact_wob: 0, blades: 5, radius: null,
         strength_source: 'constant', cutting_efficiency: 1.0, nozzle_area: 3.5628e-4, discharge_coefficient: 0.95 },
  top_drive: { inertia: 2000, max_torque: 0 },
  autodriller: { kp: 1e-7, ki: 1e-8, v_max: 60 / 3600 },
  drawworks: { block_acceleration: 0 },
  ramps: { omega_sp: 0, omega_td: 1, motor_torque: 0, v_block: 2, wob_sp: 0, flow_rate: 3, inlet_temperature: 3 },
  hydraulics: { enabled: true, dx: 30, cut_factor: 1.5, cfl: 0.8, iterations: 2, back_pressure: 0, surface_loss_coefficient: 0,
                surface_loss_exponent: 1.8, clinging_constant: 0.45 },
  thermal: { enabled: false, dx: 30, cut_factor: 1.5, n_radial: 12, far_field: 5, steel_conductivity: 45, casing_thickness: 0.012,
             flow_rate: 0.0347, inlet_temperature: 30, max_dt: 2, gs_iterations: 3 },
};
const clone = (o) => JSON.parse(JSON.stringify(o));

// Well A (Aarsnes & Shor 2018); the lite version runs without the wellbore temperature model
export function wellA() {
  const p = clone(DEFAULTS);
  Object.assign(p.friction, { mu_static_cased: 0.40, mu_kinetic_cased: 0.22, mu_static_open: 0.45, mu_kinetic_open: 0.25 });
  Object.assign(p.bit, { model: 'cutting', strength_source: 'formation' });
  Object.assign(p.thermal, { enabled: false, flow_rate: 2000 / 60000, inlet_temperature: 30.0 });
  return {
    name: 'Well A (Aarsnes & Shor 2018)',
    well: { survey: [{ md: 0, inc: 0, azi: 0 }, { md: 1500, inc: 0, azi: 0 }, { md: 2100, inc: 60, azi: 0 }, { md: 3000, inc: 60, azi: 0 }],
            casing_shoe: 1000.0, hole_diameter: 8.5 * IN },
    fluid: { ...FLUID_DEFAULTS, density: 1200 },
    formation: clone(FORMATION_DEFAULTS),
    drillstring: [
      { name: 'Collars', length: 230, od: 5.75 * IN, id: 2.8125 * IN, linear_mass: null, contact_od: null },
      { name: 'DP', length: 100, od: 5.0 * IN, id: 4.276 * IN, linear_mass: null, contact_od: null },
    ],
    parameters: { ...p, controller: { type: 'stiff', f_c: null, t_hp: 10.0, t_lp: 0.001, t_en: 0.001, kp: null, ki: null, impedance: null }, cfl: 0.95 },
    simulation: { bit_depth: 2500, hole_depth: 2505, dx: 10, record_rate: 50, start_time: null, block_height: 25.0 },
    display: { rpm_limit: 220, wob_limit: 150e3, torque_limit: 25e3, rop_limit: 90 / 3600, spp_limit: 4000 * PSI, pore_density: 1100, frac_density: 1500 },
  };
}

const num = (x, d = null) => (x === null || x === undefined || x === '' ? d : Number(x));

export function validate(cfg) {
  const errs = [];
  try {
    const md = cfg.well.survey.map(r => Number(r.md));
    if (md.length < 2) errs.push('The survey needs at least two stations.');
    if (md.some((m, i) => i > 0 && m <= md[i - 1])) errs.push('Survey measured depths must increase.');
    const bit = Number(cfg.simulation.bit_depth), bha = cfg.drillstring.slice(0, -1).reduce((s, r) => s + Number(r.length), 0);
    if (bit <= bha) errs.push(`Bit depth (${bit.toFixed(0)} m) must exceed the BHA length (${bha.toFixed(0)} m).`);
    const hole = cfg.simulation.hole_depth;
    if (hole !== null && hole !== '' && hole !== undefined && Number(hole) < bit) errs.push('Hole depth must be at or below the bit depth.');
    for (const r of cfg.drillstring) if (Number(r.id) >= Number(r.od)) errs.push(`Component '${r.name || ''}': ID must be smaller than OD.`);
    if (Number(cfg.simulation.dx ?? 10) <= 0) errs.push('Cell size must be positive.');
    const st = cfg.simulation.start_time;
    if (st && isNaN(Date.parse(String(st).trim().replace(' ', 'T')))) errs.push('Start date/time must look like 2021-01-07 06:30:00.');
    const b = cfg.parameters.bit || {};
    if (!BIT_MODELS.includes(b.model || 'coulomb')) errs.push(`Bit model must be one of ${BIT_MODELS.join(', ')}.`);
    if (b.model === 'coulomb' && Number(b.mu_kinetic ?? 0) > Number(b.mu_static ?? 1)) errs.push('Bit kinetic friction exceeds static friction.');
    if (!(Number(b.drilling_strength ?? 1) > 0)) errs.push('Bit drilling strength must be positive.');
    if (['cutting', 'regenerative'].includes(b.model) && !(Number(b.zeta ?? 1) > 0)) errs.push('Bit ζ must be positive.');
    if (b.model === 'regenerative' && !(Number(b.blades ?? 1) >= 1)) errs.push('The regenerative bit needs at least one blade.');
    const h = cfg.parameters.hydraulics || {};
    if (h.enabled !== false) {
      if (!(Number(b.nozzle_area ?? 1) > 0)) errs.push('Bit nozzle area must be positive.');
      const cd = Number(b.discharge_coefficient ?? 0.95);
      if (!(cd > 0 && cd <= 1)) errs.push('Nozzle discharge coefficient must be between 0 and 1.');
      if (!(Number(h.dx ?? 30) > 0)) errs.push('Hydraulics cell size must be positive.');
    }
    for (const [k, v] of Object.entries(cfg.parameters.ramps || {})) if (v !== null && Number(v) < 0) errs.push(`Ramp time for ${k} must not be negative.`);
    if ((b.strength_source || 'constant') === 'formation' && !cfg.formation) errs.push('The bit takes its strength from the formation, but the case has no formation.');
    for (const l of (cfg.formation || {}).layers || []) {
      if (Number(l.ucs ?? 1) <= 0) errs.push('Formation UCS must be positive.');
      if (!(Number(l.porosity ?? 0) >= 0 && Number(l.porosity ?? 0) < 1)) errs.push('Formation porosity must be between 0 and 1.');
      if (!(Number(l.friction_angle ?? 30) > 0 && Number(l.friction_angle ?? 30) < 90)) errs.push('Formation friction angle must be between 0 and 90 degrees.');
    }
  } catch (e) {
    errs.push(`Malformed case: ${e.message}`);
  }
  return errs;
}

export function makeController(c, zetaP, inertia, firstMode) {
  const kind = c.type || 'stiff', fc = num(c.f_c) || firstMode || 0.2;
  let ctrl;
  if (kind === 'stiff') ctrl = stiffController(zetaP, inertia);
  else if (kind === 'soft_torque') ctrl = { kp: 4 * zetaP, ki: (2 * Math.PI * fc) ** 2 * inertia, impedance: 0, t_hp: 10, t_lp: 0.005, t_en: 0.001, name: 'soft_torque' };
  else if (kind === 'z_torque') ctrl = { ...stiffController(zetaP, inertia), impedance: 1 / zetaP, t_hp: num(c.t_hp, 10), t_lp: num(c.t_lp, 0.001), t_en: num(c.t_en, 0.001), name: 'z_torque' };
  else if (kind === 'custom') ctrl = { kp: num(c.kp, 100 * zetaP), ki: num(c.ki, 5 * inertia), impedance: num(c.impedance, 0), t_hp: num(c.t_hp, 10), t_lp: num(c.t_lp, 0.005), t_en: num(c.t_en, 0.001), name: 'custom' };
  else throw new Error(`unknown controller type '${kind}'`);
  if (kind !== 'custom') for (const k of ['kp', 'ki', 'impedance']) if (c[k] !== null && c[k] !== undefined && c[k] !== '') ctrl[k] = Number(c[k]);
  return ctrl;
}

function section(p, key) { return { ...DEFAULTS[key], ...(p[key] || {}) }; }

export function build(cfg) {
  const w = cfg.well;
  const survey = new Survey(w.survey.map(r => num(r.md)), w.survey.map(r => num(r.inc)), w.survey.map(r => num(r.azi, 0)));
  const wellbore = { survey, casing_shoe: num(w.casing_shoe, 0), hole_diameter: num(w.hole_diameter), casing_id: null,
                     isCased(md) { return md <= this.casing_shoe; },
                     wallDiameter(md) { return this.casing_id != null && this.isCased(md) ? this.casing_id : this.hole_diameter; } };
  const string = new Drillstring(cfg.drillstring.map(r => ({ name: r.name, length: num(r.length), od: num(r.od), id: num(r.id),
                                                             linear_mass: num(r.linear_mass), contact_od: num(r.contact_od) })));
  const fluid = { ...FLUID_DEFAULTS };
  for (const [k, v] of Object.entries(cfg.fluid || {})) if (k in FLUID_DEFAULTS) fluid[k] = (k === 'plastic_viscosity' || k === 'yield_point') ? num(v) : num(v, FLUID_DEFAULTS[k]);
  const p = cfg.parameters;
  const params = {};
  for (const k of Object.keys(DEFAULTS)) params[k] = section(p, k);
  params.bit.blades = Math.trunc(params.bit.blades);
  params.thermal.enabled = false;                    // lite: no wellbore temperature model (spec 10)
  params.cfl = num(p.cfl, 0.95);
  params.thermal_update_period = 2.0;
  params.controller = null;
  let formation = null;
  if (cfg.formation) {
    formation = { ...clone(FORMATION_DEFAULTS), ...clone(cfg.formation) };
    if (!formation.layers || !formation.layers.length) formation.layers = clone(FORMATION_DEFAULTS.layers);
  }
  const sim = new Simulation(wellbore, string, fluid, params, formation);
  const s = cfg.simulation;
  sim.initialize(num(s.bit_depth), num(s.hole_depth), num(s.dx, 10), num(s.record_rate, 50));
  params.controller = makeController(p.controller || { type: 'stiff' }, sim.zetaP, params.top_drive.inertia, sim.modalFrequencies()[1]);
  sim.updateParameters();
  return sim;
}
