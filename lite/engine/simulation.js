// Simulation: set-up, inputs and ramps, slow tick and outputs (spec sections 2.4, 3, 4, 6, 7, 12)
import { G, STEEL_DENSITY } from './units.js';
import { advance, rampValue, BIT_MODELS, EN_NAMES, R_WSP, R_VBLK, R_TAUM, R_WTD, R_WOB, R_FLOW } from './kernel.js';
import { WellboreHydraulics } from './hydraulics.js';
import { rockStrength } from './formation.js';
import { WellboreThermal } from './thermal.js';

const F64 = (n, v = 0) => new Float64Array(n).fill(v);
// hydraulic energy terms (spec 8.4), from WellboreHydraulics.lossSplit
export const HYD_ENERGY = ['pump', 'surface', 'string', 'motor', 'nozzles', 'annulus', 'other'];

export function stiffController(zetaP, inertia) {
  return { kp: 100 * zetaP, ki: 5 * inertia, impedance: 0, t_hp: 10, t_lp: 0.005, t_en: 0.001, name: 'stiff' };
}

export class Simulation {
  constructor(wellbore, drillstring, fluid, params, formation = null) {
    Object.assign(this, { wellbore, drillstring, fluid, params, formation });
    this.rock = null; this.hydraulics = null; this.thermal = null;
  }

  // ------------------------------------------------------------- set-up
  initialize(bitDepth, holeDepth = null, dx = 10, recordRate = 100, refreshDistance = 1) {
    const P = this.params, g = this.grid = this.drillstring.discretize(bitDepth, dx), n = this.n = g.n;
    this.dx = g.dx; this.rc = g.rc;
    this.Kt = g.J.map((J, j) => J * g.G[j]); this.mt = g.J.map((J, j) => g.rho[j] * J);
    this.Zt = g.J.map((J, j) => J * Math.sqrt(g.G[j] * g.rho[j]));
    this.Ka = g.A.map((A, j) => g.E[j] * A); this.ma = g.A.map((A, j) => g.rho[j] * A);
    this.Za = g.A.map((A, j) => A * Math.sqrt(g.E[j] * g.rho[j]));
    let m = Infinity; for (let j = 0; j < n; j++) m = Math.min(m, g.dx[j] / g.cA[j]);
    this.dt = P.cfl * m;
    this.zetaP = this.Zt[0];
    this.wb = g.A.map((A, j) => g.rho[j] * A * G * (1 - this.fluid.density / STEEL_DENSITY));
    for (const k of ['w', 'tau', 'v', 'qd', 'wbar', 'vbar']) this[k] = F64(n);
    for (const k of ['wf', 'tf', 'vf', 'qf']) this[k] = F64(n + 1);
    this.tmp = { wst: F64(n), vst: F64(n), wOld: F64(n), vOld: F64(n), cpr: F64(n), dpr: F64(n) };
    const hole = holeDepth == null ? bitDepth : Math.max(holeDepth, bitDepth);
    this.sc = { t: 0, wtd: 0, taum: 0, tau0: 0, int: 0, lptau: 0, hplp: 0, wfilt: 0, wsp: 0, bitd: bitDepth, holed: hole,
                wbit: 0, tob: 0, vbit: 0, wob: 0, vblk: 0, hkld: 0, adInt: 0, blkpos: 0, stuck: 0, phib: 0, wint: 0,
                phis: 0, dpm: 0, dpint: 0, wmot: 0, stall: 0 };
    this.surf = F64(720, hole);
    this.en = F64(EN_NAMES.length);                         // energy accumulators (J), spec 5.9
    this.enHyd = Object.fromEntries(HYD_ENERGY.map(k => [k, 0]));   // hydraulic energy (J), spec 8.4
    this.bitDepth0 = bitDepth; this.geomBitd = null; this.qsC = null; this.refreshDistance = refreshDistance;
    this.p = {};
    this.refreshGeometry();
    this.sc.hkld = -this.qsTop;
    this.ctrl = 0; this.axial = 0;
    this.controller = null;
    this.packParameters();
    this.p.hkldRef = -this.qsTop;
    this.inputs = { omega_sp: 0, v_block: 0, motor_torque: 0, omega_td: 0, wob_sp: 0,
                    flow_rate: P.thermal.flow_rate, inlet_temperature: P.thermal.inlet_temperature };
    this.ramps = [0, 1, 2, 3, 4, 5].map(() => [0, 0, 0, 0]);
    this.slowRamps = { flow_rate: [0, 0, P.ramps.flow_rate, P.thermal.flow_rate],
                       inlet_temperature: [P.thermal.inlet_temperature, 0, 0, P.thermal.inlet_temperature] };
    this.ramps[R_FLOW] = this.slowRamps.flow_rate;     // the kernel sees the pump rate (mud motor)
    this.slowEvery = Math.max(1, Math.trunc(1 / this.dt));
    this.thermalEvery = Math.max(1, Math.round(P.thermal_update_period));
    this.slowTicks = 0; this.slowRecords = []; this.lastSlow = { t: 0, hole, wint: 0 };
    this.hydRecords = []; this.hydLast = { t: 0, q: 0, dpint: 0 };
    // wellbore temperature (spec 10), then hydraulics with its mud temperatures (10.7)
    if (this.formation && P.thermal.enabled)
      this.thermal = new WellboreThermal(this.wellbore, this.drillstring, this.fluid, this.formation, P.thermal, bitDepth, hole);
    if (P.hydraulics.enabled) {
      this.hydraulics = new WellboreHydraulics(this.wellbore, g, this.fluid, P.hydraulics, P.bit.nozzle_area,
                                               P.bit.discharge_coefficient, bitDepth, hole);
      this.syncHydraulicsTemperature();
      this.hydraulics.setStatic();
    }
    this.updateRock(0, 0);
    this.recordRate = recordRate;
    this.recEvery = Math.max(1, Math.round(1 / (recordRate * this.dt)));
    this.recPhase = 0; this.rec = []; this.stepsTotal = 0; this.tRequested = 0;
    return this;
  }

  modalFrequencies() {
    let ta = 0, tt = 0; const g = this.grid;
    for (let j = 0; j < this.n; j++) { ta += g.dx[j] / g.cA[j]; tt += g.dx[j] / g.cT[j]; }
    return [1 / (4 * ta), 1 / (4 * tt)];
  }

  // spec 4: kernel parameters
  packParameters() {
    const P = this.params, f = P.friction, b = P.bit, td = P.top_drive, p = this.p;
    const c = this.controller = P.controller || stiffController(this.zetaP, td.inertia);
    Object.assign(p, { dt: this.dt, wc: f.omega_c, vc: f.v_c, itd: td.inertia, kp: c.kp, ki: c.ki, zc: c.impedance,
                       thp: c.t_hp, tlp: c.t_lp, ten: c.t_en, taumax: td.max_torque });
    const D = this.wellbore.hole_diameter, a = D / 2;
    const S = b.strength_source === 'formation' && this.S ? this.S : b.drilling_strength;
    if (!BIT_MODELS.includes(b.model)) throw new Error(`unknown bit model ${b.model}`);
    p.bitModel = BIT_MODELS.indexOf(b.model);
    if (b.model === 'coulomb') {
      const r = b.radius != null ? b.radius : D / 3;
      p.levS = b.mu_static * r; p.levK = b.mu_kinetic * r; p.w0 = 0;
    } else {
      p.levS = p.levK = a / (2 * b.zeta); p.w0 = b.model === 'cutting' ? b.contact_wob : 0;
    }
    p.krop = 1 / (2 * Math.PI * a * S); p.kcut = S * a;
    p.blades = Math.max(1, Math.trunc(b.blades));
    p.qsTop = this.qsTop;
    p.adKp = P.autodriller.kp; p.adKi = P.autodriller.ki; p.adVmax = P.autodriller.v_max;
    p.coupled = !!f.along_string_coupling; p.axc = f.axial_static_contrast;
    if (!['step', 'exponential'].includes(f.weakening)) throw new Error("friction weakening must be 'step' or 'exponential'");
    p.weakExp = f.weakening === 'exponential';
    p.blkAcc = P.drawworks.block_acceleration;
    const m = P.motor;
    Object.assign(p, { motor: !!m.enabled, mKn: 2 * Math.PI * m.speed_constant, mKt: m.efficiency / (2 * Math.PI * m.speed_constant),
                       mDp0: m.dp_no_load, mDpMax: m.dp_max, mDpStall: m.dp_stall, mSlip: m.slip, mQby: m.bypass_flow });
    this.packDamping();
  }

  // spec 4.3
  packDamping() {
    const d = this.params.damping, g = this.grid, fl = this.fluid, n = this.n, p = this.p;
    const [fa1, ft1] = this.modalFrequencies();
    const fa = d.reference_frequency_axial || fa1, ft = d.reference_frequency_torsional || ft1;
    const wA = 2 * Math.PI * fa, wT = 2 * Math.PI * ft;
    this.k0t = F64(n, d.k_t); this.k0a = F64(n, d.k_a); this.k1t = F64(n); this.k1a = F64(n);
    const mu = fl.plastic_viscosity != null ? fl.plastic_viscosity : fl.viscosity, rho = fl.density;
    if (d.fluid && mu && mu > 0) {
      for (let j = 0; j < n; j++) {
        const ro = g.od[j] / 2, ri = g.id[j] / 2;
        const R = Math.max(this.wellbore.wallDiameter(this.cellMd[j]) / 2, 1.05 * ro);
        const c0a = 2 * Math.PI * mu / Math.log(R / ro), c0t = 4 * Math.PI * mu * ro * ro * R * R / (R * R - ro * ro);
        const c1a = Math.sqrt(rho * mu * wA / 2) * 2 * Math.PI * (ro + ri), c1t = Math.sqrt(rho * mu * wT / 2) * 2 * Math.PI * (ro ** 3 + ri ** 3);
        this.k0a[j] += c0a / this.ma[j]; this.k0t[j] += c0t / this.mt[j];
        this.k1a[j] = c1a / this.ma[j]; this.k1t[j] = c1t / this.mt[j];
      }
    }
    const fc = d.high_pass_fraction, dt = this.dt;
    p.hpT = fc > 0 ? 1 - Math.exp(-dt * fc * wT) : 1;
    p.hpA = fc > 0 ? 1 - Math.exp(-dt * fc * wA) : 1;
    this.kvt = F64(n + 1); this.kva = F64(n + 1);
    p.kv = !!(d.loss_factor && d.loss_factor > 0);
    if (p.kv) {
      const z = 0.5 * d.loss_factor, w2 = 2 * Math.PI * d.max_frequency, ray = {};
      for (const [key, w1] of [['t', wT], ['a', wA]]) {
        if (w2 > w1) { ray['alpha_' + key] = 2 * z * w1 * w2 / (w1 + w2); ray['beta_' + key] = 2 * z / (w1 + w2); }
        else { ray['alpha_' + key] = 0; ray['beta_' + key] = 2 * z / w1; }
      }
      for (let j = 0; j < n; j++) { this.k1t[j] += ray.alpha_t; this.k1a[j] += ray.alpha_a; }
      for (let j = 1; j < n; j++) {
        const hf = 0.5 * (g.dx[j - 1] + g.dx[j]);
        this.kvt[j] = ray.beta_t * 0.5 * (this.Kt[j - 1] + this.Kt[j]) / hf;
        this.kva[j] = ray.beta_a * 0.5 * (this.Ka[j - 1] + this.Ka[j]) / hf;
      }
    }
    if (!['always', 'pumps_off', 'off'].includes(d.yield_drag)) throw new Error("damping.yield_drag must be 'always', 'pumps_off' or 'off'");
    const ty = fl.yield_point || 0;
    this.fyFull = ty > 0 && d.yield_drag !== 'off' ? g.od.map(od => ty * Math.PI * od) : F64(n);
    this.fyZero = F64(n);
  }
  yieldDragNow() {
    const d = this.params.damping;
    if (d.yield_drag === 'always') return this.fyFull;
    if (d.yield_drag === 'pumps_off' && this.flowRateNow() <= 0) return this.fyFull;
    return this.fyZero;
  }

  // spec 2.4
  refreshGeometry(force = false) {
    const bitd = this.sc.bitd;
    if (!force && this.geomBitd != null && Math.abs(bitd - this.geomBitd) < this.refreshDistance) return;
    const g = this.grid, n = this.n, sv = this.wellbore.survey, shift = bitd - this.bitDepth0;
    const top = g.xTop.map(x => x + shift);
    const tE = []; for (let j = 0; j <= n; j++) tE.push(sv.tangent(Math.max(j < n ? top[j] : top[n - 1] + g.dx[n - 1], 0)));
    const geo = { kx: F64(n), ky: F64(n), kz: F64(n), zx: F64(n), zy: F64(n), zz: F64(n) };
    this.tz = F64(n); this.cellMd = F64(n);
    for (let j = 0; j < n; j++) {
      const tc = sv.tangent(Math.max(top[j] + 0.5 * g.dx[j], 0));
      geo.kx[j] = (tE[j + 1][0] - tE[j][0]) / g.dx[j]; geo.ky[j] = (tE[j + 1][1] - tE[j][1]) / g.dx[j]; geo.kz[j] = (tE[j + 1][2] - tE[j][2]) / g.dx[j];
      geo.zx[j] = 0 - tc[2] * tc[0]; geo.zy[j] = 0 - tc[2] * tc[1]; geo.zz[j] = 1 - tc[2] * tc[2];
      this.tz[j] = tc[2]; this.cellMd[j] = top[j] + 0.5 * g.dx[j];
    }
    Object.assign(this, geo);
    const f = this.params.friction;
    this.mus = this.cellMd.map(m => (this.wellbore.isCased(m) ? f.mu_static_cased : f.mu_static_open));
    this.muk = this.cellMd.map(m => (this.wellbore.isCased(m) ? f.mu_kinetic_cased : f.mu_kinetic_open));
    // static compression, well-balanced split
    const qsE = F64(n + 1);
    for (let j = n - 1; j >= 0; j--) qsE[j] = qsE[j + 1] - this.wb[j] * this.tz[j] * g.dx[j];
    const qsNew = F64(n); for (let j = 0; j < n; j++) qsNew[j] = 0.5 * (qsE[j] + qsE[j + 1]);
    if (this.qsC) for (let j = 0; j < n; j++) this.qd[j] += this.qsC[j] - qsNew[j];
    this.qsC = qsNew; this.qsTop = qsE[0];
    this.fn0 = F64(n);
    for (let j = 0; j < n; j++) {
      const T = -qsNew[j];
      this.fn0[j] = Math.hypot(T * this.kx[j] + this.wb[j] * this.zx[j], T * this.ky[j] + this.wb[j] * this.zy[j], T * this.kz[j] + this.wb[j] * this.zz[j]);
    }
    if (this.p) this.p.qsTop = this.qsTop;
    this.geomBitd = bitd;
  }
  normalForceNow() {
    const out = F64(this.n);
    for (let j = 0; j < this.n; j++) {
      const T = -(this.qsC[j] + this.qd[j]);
      out[j] = Math.hypot(T * this.kx[j] + this.wb[j] * this.zx[j], T * this.ky[j] + this.wb[j] * this.zy[j], T * this.kz[j] + this.wb[j] * this.zz[j]);
    }
    return out;
  }

  setDrillingStrength(S) {
    this.S = S;
    const a = this.wellbore.hole_diameter / 2;
    this.p.krop = 1 / (2 * Math.PI * a * S); this.p.kcut = S * a;
  }
  // spec 9 (T2 uses the mud temperature at the bit from the thermal model)
  updateRock(rop, omega) {
    if (!this.formation) return;
    const b = this.params.bit, hole = this.sc.holed, tvd = this.wellbore.survey.tvdAt(hole);
    let tMud = null, pbh = this.fluid.density * G * tvd;
    if (this.thermal) { const ts = this.thermal.state(); pbh = ts.p_bh; tMud = ts.t_bit; }
    if (this.hydraulics) pbh = this.hydraulics.state().p_bh;
    this.rock = rockStrength(this.formation, hole, tvd, pbh, tMud, rop, omega, b.blades);
    this.rock.p_bh = pbh;
    this.rock.S = b.zeta * this.rock.ccs_eff / b.cutting_efficiency;
    if (b.strength_source === 'formation') this.setDrillingStrength(this.rock.S);
  }

  // ------------------------------------------------------------- ramps (spec 6)
  slowValue(name, t = this.sc.t) {
    const [v0, t0, T, tgt] = this.slowRamps[name];
    if (T <= 0 || t >= t0 + T) return tgt;
    return v0 + (tgt - v0) * Math.max(0, (t - t0) / T);
  }
  flowRateNow() { return this.slowValue('flow_rate'); }
  rampNow(k) {
    const [v0, t0, T, tgt] = this.ramps[k], t = this.sc.t;
    if (T <= 0 || t >= t0 + T) return tgt;
    return v0 + (tgt - v0) * Math.max(0, (t - t0) / T);
  }
  setRamp(k, target, rampTime, start = null) {
    if (start == null) { if (this.ramps[k][3] === target) return; start = this.rampNow(k); }
    this.ramps[k] = [start, this.sc.t, rampTime, target];
  }
  setSlowRamp(name, target, rampTime) {
    if (this.slowRamps[name][3] === target) return;
    this.slowRamps[name] = [this.slowValue(name), this.sc.t, rampTime, target];
  }

  updateParameters() {
    this.refreshGeometry(true);
    this.packParameters();
    this.updateRock(0, 0);
  }

  syncHydraulicsTemperature() {
    if (!this.hydraulics || !this.thermal) return;
    const pr = this.thermal.profiles(), idx = pr.md.map((_, i) => i).filter(i => pr.has_pipe[i]);
    const use = idx.length ? idx : pr.md.map((_, i) => i);
    this.hydraulics.setTemperatures(use.map(i => pr.md[i]), use.map(i => pr.t_ds[i]), pr.md, pr.t_ann);
  }

  // ------------------------------------------------------------- slow tick (spec 7)
  hydraulicsTick() {
    const t = this.sc.t, dt = t - this.hydLast.t;
    if (!this.hydraulics || dt <= 0) return;
    const q1 = this.flowRateNow(), dpint = this.sc.dpint;
    // mud motor differential, averaged over the interval, in series with the bit nozzles
    const hs = this.hydraulics.advance(dt, this.hydLast.q, q1, this.sc.bitd, this.sc.holed, (dpint - this.hydLast.dpint) / dt);
    this.hydLast = { t, q: q1, dpint };
    this.hydRecords.push({ time: t, ...hs });
    // hydraulic energy over the interval: each pressure drop times the pump rate
    const ls = this.hydraulics.lossSplit();
    for (const k of HYD_ENERGY) this.enHyd[k] += (k === 'pump' ? ls.spp : ls[k]) * ls.q * dt;
  }
  slowTick() {
    this.slowTicks++;
    this.refreshGeometry();
    this.hydraulicsTick();
    if (this.slowTicks % this.thermalEvery) return;
    const { t, holed: hole, wint } = this.sc, dt = t - this.lastSlow.t;
    if (dt <= 0) return;
    const rop = Math.max(hole - this.lastSlow.hole, 0) / dt, omega = (wint - this.lastSlow.wint) / dt;
    this.lastSlow = { t, hole, wint };
    let ts = null;
    if (this.thermal) {
      ts = this.thermal.advance(dt, this.flowRateNow(), this.slowValue('inlet_temperature'), this.sc.bitd, hole, rop);
      this.syncHydraulicsTemperature();
    }
    this.updateRock(rop, omega);
    const rec = { time: t, rop, omega_bit_avg: omega, flow_rate: this.flowRateNow() };
    if (ts) for (const k of ['t_bit', 't_ds_bit', 't_out', 't_in', 'ecd']) rec[k] = ts[k];
    if (this.hydraulics) { const hs = this.hydraulics.state(); rec.spp = hs.spp; rec.ecd = hs.ecd; }
    if (this.rock) for (const k of ['t_rock', 'p_bh', 'ucs_t', 'ccs', 'ccs_eff', 'sigma_th', 'skin_ratio', 'w_drained', 'S']) rec[k] = this.rock[k];
    this.slowRecords.push(rec);
  }

  // ------------------------------------------------------------- step (spec 6)
  step(duration, inp = {}) {
    const r = this.params.ramps, ctrl = this.ctrl, axial = this.axial, sc = this.sc;
    if (inp.omega_sp != null) { this.inputs.omega_sp = inp.omega_sp; this.setRamp(R_WSP, inp.omega_sp, r.omega_sp, ctrl === 0 ? null : sc.wtd); this.ctrl = 0; }
    if (inp.motor_torque != null) { this.inputs.motor_torque = inp.motor_torque; this.setRamp(R_TAUM, inp.motor_torque, r.motor_torque, ctrl === 1 ? null : sc.taum); this.ctrl = 1; }
    if (inp.omega_td != null) { this.inputs.omega_td = inp.omega_td; this.setRamp(R_WTD, inp.omega_td, r.omega_td, ctrl === 2 ? null : sc.wtd); this.ctrl = 2; }
    if (inp.v_block != null) { this.inputs.v_block = inp.v_block; this.setRamp(R_VBLK, inp.v_block, r.v_block, axial === 0 ? null : sc.vblk); this.axial = 0; }
    if (inp.wob_sp != null) { this.inputs.wob_sp = inp.wob_sp; this.setRamp(R_WOB, inp.wob_sp, r.wob_sp, axial === 1 ? null : this.p.hkldRef - sc.hkld); this.axial = 1; }
    if (inp.flow_rate != null) { this.inputs.flow_rate = inp.flow_rate; this.setSlowRamp('flow_rate', inp.flow_rate, r.flow_rate); this.ramps[R_FLOW] = this.slowRamps.flow_rate; }
    if (inp.inlet_temperature != null) { this.inputs.inlet_temperature = inp.inlet_temperature; this.setSlowRamp('inlet_temperature', inp.inlet_temperature, r.inlet_temperature); }
    this.tRequested += duration;
    let remaining = Math.round(this.tRequested / this.dt) - this.stepsTotal;
    const K = this.slowEvery;
    while (remaining > 0) {
      const nsteps = Math.min(remaining, K - (this.stepsTotal % K));
      advance(this, nsteps);
      this.stepsTotal += nsteps;
      remaining -= nsteps;
      if (this.stepsTotal % K === 0) this.slowTick();
    }
    return this.latest();
  }

  tareHookload() { this.p.hkldRef = this.sc.hkld; return this.p.hkldRef; }
  get hookloadReference() { return this.p.hkldRef; }

  // ------------------------------------------------------------- outputs (spec 12)
  latest() {
    const s = this.sc;
    const out = { time: s.t, omega_td: s.wtd, motor_torque: s.taum, pipe_torque: s.tau0, omega_bit: s.wbit, tob: s.tob,
                  v_block: s.vblk, hookload: s.hkld, v_bit: s.vbit, wob: s.wob, bit_depth: s.bitd, hole_depth: s.holed,
                  omega_sp: s.wsp, flow_rate: this.flowRateNow() };
    if (this.hydraulics) { const h = this.hydraulics.state(); Object.assign(out, { spp: h.spp, p_bh: h.p_bh, ecd: h.ecd, q_out: h.q_out }); }
    if (this.params.motor.enabled) Object.assign(out, this.motorState());
    return out;
  }
  // mud motor channels (spec 12): rotor speed, differential, stall flag, toolface (rad, 0..2 pi)
  // energy accounting (J, cumulative since the start; spec 5.9 and 8.4)
  energy() {
    const h = this.grid.dx, n = this.w.length;
    const out = { time: this.sc.t };
    EN_NAMES.forEach((k, i) => { out[k] = this.en[i]; });
    for (const k of HYD_ENERGY) out['h_' + k] = this.enHyd[k];
    let st = 0, sa = 0;
    for (let j = 0; j < n; j++) {
      st += h[j] * (0.5 * this.mt[j] * this.w[j] ** 2 + 0.5 * this.tau[j] ** 2 / this.Kt[j]);
      sa += h[j] * (0.5 * this.ma[j] * this.v[j] ** 2 + 0.5 * this.qd[j] ** 2 / this.Ka[j]);
    }
    out.store_torsion = st; out.store_axial = sa;
    out.store_top_drive = 0.5 * this.params.top_drive.inertia * this.sc.wtd ** 2;
    return out;
  }

  motorState() {
    const s = this.sc, tf = this.params.motor.toolface_offset + s.phis, tp = 2 * Math.PI;
    return { omega_motor: s.wmot, dp_motor: s.dpm, motor_stalled: s.stall > 0.5, toolface: ((tf % tp) + tp) % tp };
  }
  profiles() {
    const n = this.n;
    return { md: Array.from(this.cellMd), omega: Array.from(this.w), torque: Array.from(this.tau), v: Array.from(this.v),
             tension: Array.from({ length: n }, (_, j) => -(this.qsC[j] + this.qd[j])),
             normal_force: Array.from(this.params.friction.along_string_coupling ? this.normalForceNow() : this.fn0),
             od: Array.from(this.grid.od), inc: Array.from(this.tz, z => Math.acos(Math.max(-1, Math.min(1, z)))) };
  }
  info() {
    const g = this.grid;
    let mx = 0; for (let j = 0; j < this.n; j++) mx = Math.max(mx, this.dt * g.cT[j] / g.dx[j]);
    return { cells: this.n, dt: this.dt, dx_min: Math.min(...g.dx), dx_max: Math.max(...g.dx), length: g.dx.reduce((a, b) => a + b, 0),
             zeta_p: this.zetaP, controller: this.controller.name, torsional_cfl: mx };
  }
}
