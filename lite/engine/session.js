// UI session: the same frames as the Python server (src/DrillPyze/ui/server.py, Session) so the
// unchanged front end can run on the browser engine (spec section 12).
import { build, validate } from './config.js';
import { torqueDrag } from './statics.js';
import { REC_CHANNELS } from './kernel.js';
import { HYD_ENERGY } from './simulation.js';
import { RPM } from './units.js';

const MAX_TRACE_ROWS = 400;
const rnd = (x, nd = 4) => { const f = 10 ** nd; const v = Math.round(x * f) / f; return Number.isFinite(v) ? v : 0; };
const R = (arr, nd = 4) => Array.from(arr, x => rnd(x, nd));
const rpm = (w) => w / RPM;
const C = Object.fromEntries(REC_CHANNELS.map((c, i) => [c, i]));

export class Session {
  constructor(cfg) {
    const errs = validate(cfg);
    if (errs.length) throw new Error(errs.join(' '));
    this.cfg = cfg;
    this.sim = build(cfg);
    this.sentRows = 0; this.sentSlow = 0; this.sentHyd = 0;
    const bh = (cfg.simulation || {}).block_height;
    this.blockPos = bh === null || bh === undefined || bh === '' ? 25.0 : Number(bh);
    this.blockT = null;
    this.buildReference();
  }

  buildReference() {
    const s = this.sim, f = s.params.friction, n = s.n;
    const geo = { n, dx: s.grid.dx, rc: s.grid.rc, wb: s.wb, tz: s.tz, kx: s.kx, ky: s.ky, kz: s.kz, zx: s.zx, zy: s.zy, zz: s.zz };
    const fixed = f.along_string_coupling ? null : s.fn0;
    const kin = torqueDrag(geo, s.muk, 0, 1, 0, 0, fixed), stat = torqueDrag(geo, s.mus, 0, 1, 0, 0, fixed);
    const hang = torqueDrag(geo, new Float64Array(n), 0, 0);
    const mdE = Array.from({ length: n + 1 }, (_, j) => (j < n ? s.grid.xTop[j] : s.grid.xTop[n - 1] + s.grid.dx[n - 1]));
    this.reference = { md: R(mdE, 2), torque_rotating_knm: R(kin.torque.map(x => x / 1e3)),
                       torque_breakoff_knm: R(stat.torque.map(x => x / 1e3)), tension_hanging_kn: R(hang.tension.map(x => x / 1e3)) };
    const sv = s.wellbore.survey, top = Math.max(sv.td, s.sc.holed);
    const mdPath = Array.from({ length: 400 }, (_, i) => top * i / 399);
    const [nn, ee, zz] = sv.positions(mdPath);
    const az = sv.azi[sv.azi.length - 1];
    this.vsAz = az;
    const topComp = s.drillstring.components.length - 1;
    this.path = { md: R(mdPath, 1), vs: R(nn.map((x, i) => x * Math.cos(az) + ee[i] * Math.sin(az)), 2), tvd: R(zz, 2),
                  shoe_md: s.wellbore.casing_shoe, bha: Array.from(s.grid.comp, c => (c !== topComp ? 1 : 0)),
                  bha_length: s.drillstring.bhaLength };
  }

  applyInputs(inp) {
    const kw = {};
    if (inp.tare) this.sim.tareHookload();
    const r = Number(inp.rpm || 0);
    if ((inp.td_mode || 'setpoint') === 'speed') kw.omega_td = r * RPM; else kw.omega_sp = r * RPM;
    if ((inp.dw_mode || 'velocity') === 'wob') kw.wob_sp = Number(inp.wob_kn || 0) * 1e3;
    else kw.v_block = Number(inp.v_block_mh || 0) / 3600;
    if (inp.flow_m3s !== null && inp.flow_m3s !== undefined) kw.flow_rate = Math.max(Number(inp.flow_m3s), 0);
    if (inp.inlet_c !== null && inp.inlet_c !== undefined) kw.inlet_temperature = Number(inp.inlet_c);
    return kw;
  }

  step(duration, inputs) {
    const t0 = performance.now();
    const kw = this.applyInputs(inputs || {});
    if (duration > 0) this.sim.step(duration, kw);
    const wall = (performance.now() - t0) / 1000;
    return this.frame(wall, duration);
  }

  frame(wall = 0, duration = 0) {
    const sim = this.sim;
    let rows = sim.rec.slice(this.sentRows);
    this.sentRows = sim.rec.length;
    if (rows.length) {
      let tPrev = this.blockT === null ? rows[0][C.time] : this.blockT, pos = this.blockPos;
      rows = rows.map(r => { pos -= r[C.v_block] * (r[C.time] - tPrev); tPrev = r[C.time]; return [...r, pos]; });
      this.blockPos = pos; this.blockT = tPrev;
      if (rows.length > MAX_TRACE_ROWS) rows = Array.from({ length: MAX_TRACE_ROWS }, (_, i) => rows[Math.round(i * (rows.length - 1) / (MAX_TRACE_ROWS - 1))]);
    }
    const col = (c) => rows.map(r => r[C[c]]);
    const bp = rows.map(r => r[REC_CHANNELS.length]);
    const bitd = col('bit_depth'), hole = col('hole_depth'), vbit = col('v_bit');
    const trace = {
      time: R(col('time'), 3), rpm_td: R(col('omega_td').map(rpm), 2), rpm_bit: R(col('omega_bit').map(rpm), 2),
      rpm_sp: R(col('omega_sp').map(rpm), 2), motor_torque_knm: R(col('motor_torque').map(x => x / 1e3), 3),
      pipe_torque_knm: R(col('pipe_torque').map(x => x / 1e3), 3), tob_knm: R(col('tob').map(x => x / 1e3), 3),
      hookload_kn: R(col('hookload').map(x => x / 1e3), 2), wob_kn: R(col('wob').map(x => x / 1e3), 2),
      rop_mh: R(vbit.map((v, i) => (bitd[i] >= hole[i] - 1e-6 ? Math.max(v, 0) : 0) * 3600), 2),
      block_mh: R(col('v_block').map(x => x * 3600), 2), bit_depth: R(bitd, 3), hole_depth: R(hole, 3), block_pos: R(bp, 3),
      rpm_motor: R(col('omega_motor').map(rpm), 2), dp_motor: R(col('dp_motor'), 0),
      toolface_deg: R(col('stator_angle').map(a => { const tp = 2 * Math.PI, x = this.sim.params.motor.toolface_offset + a;
                                                        return ((x % tp) + tp) % tp * 180 / Math.PI; }), 2),
    };
    const prof = sim.profiles(), sv = sim.wellbore.survey;
    const [nn, ee, zz] = sv.positions(prof.md.map(m => Math.max(m, 0)));
    const slowRows = sim.slowRecords.slice(this.sentSlow); this.sentSlow = sim.slowRecords.length;
    const hydRows = sim.hydRecords.slice(this.sentHyd); this.sentHyd = sim.hydRecords.length;
    const hyd = Object.fromEntries(['time', 'spp', 'p_bh', 'ecd', 'q_in', 'q_out', 'dp_motor'].map(k => [k, hydRows.map(r => rnd(r[k], 6))]));
    const slow = Object.fromEntries(['time', 't_bit', 't_out', 't_in', 't_rock', 'ccs_eff', 'ccs', 'S', 'ecd', 'p_bh']
      .map(k => [k, slowRows.map(r => (r[k] === undefined || r[k] === null ? null : rnd(r[k], 4)))]));
    let thermal = null;
    if (sim.thermal) {
      const tp = sim.thermal.profiles();
      thermal = { md: R(tp.md, 1), t_ds: tp.t_ds.map((v, i) => (tp.has_pipe[i] ? rnd(v, 2) : null)), t_ann: R(tp.t_ann, 2),
                  t_wall: R(tp.t_wall, 2), t_far: R(tp.t_far, 2) };
    }
    const ts = sim.thermal ? sim.thermal.state() : {};
    const rock = sim.rock || {}, hs = sim.hydraulics ? sim.hydraulics.state() : {}, L = sim.latest();
    const onBottom = L.bit_depth >= L.hole_depth - 1e-6, rop = onBottom ? Math.max(L.v_bit, 0) : 0;
    return {
      time: L.time,
      latest: {
        rpm_td: rpm(L.omega_td), rpm_bit: rpm(L.omega_bit), rpm_sp: rpm(L.omega_sp),
        motor_torque_knm: L.motor_torque / 1e3, pipe_torque_knm: L.pipe_torque / 1e3, tob_knm: L.tob / 1e3,
        hookload_kn: L.hookload / 1e3, wob_kn: L.wob / 1e3, rop_mh: rop * 3600, bit_axial_mh: L.v_bit * 3600,
        bit_depth: L.bit_depth, hole_depth: L.hole_depth, on_bottom: onBottom,
        t_bit: ts.t_bit ?? null, t_out: ts.t_out ?? null, ecd: ts.ecd ?? null,
        t_rock: rock.t_rock ?? null, ccs_eff: rock.ccs_eff ?? null, ccs: rock.ccs ?? null, S: rock.S ?? null, p_bh: rock.p_bh ?? null,
        spp: hs.spp ?? null, ecd_hyd: hs.ecd ?? null, q_out: hs.q_out ?? null, q_in: hs.q_in ?? null, pit_gain: hs.pit_gain ?? null,
        hk_ref_kn: sim.hookloadReference / 1e3, wob_surf_kn: (sim.hookloadReference - L.hookload) / 1e3,
        block_mh: L.v_block * 3600, wob_sp_kn: sim.axial === 1 ? sim.inputs.wob_sp / 1e3 : null,
        p_bh_hyd: hs.p_bh ?? null, tvd_bottom: sim.hydraulics ? sim.hydraulics.tvdBottom : null, mud_density: sim.fluid.density,
        ...this.motorLatest(L),
      },
      hyd, thermal, slow,
      profile: {
        md: R(prof.md, 2), rpm: R(prof.omega.map(rpm), 2), torque_knm: R(prof.torque.map(x => x / 1e3), 3),
        tension_kn: R(prof.tension.map(x => x / 1e3), 2), v_axial_mms: R(prof.v.map(x => x * 1e3), 3),
        normal_knm: R(prof.normal_force.map(x => x / 1e3), 4),
        vs: R(nn.map((x, i) => x * Math.cos(this.vsAz) + ee[i] * Math.sin(this.vsAz)), 2), tvd: R(zz, 2),
      },
      trace,
      energy: this.energy(),
      perf: { wall, sim: duration, realtime: wall > 0 ? duration / wall : null },
    };
  }

  // cumulative energy (J) for the energy panel; the hydraulic terms are carried forward from
  // the last hydraulics tick at the latest split so that short averaging windows stay smooth
  energy() {
    const sim = this.sim, e = sim.energy();
    if (sim.hydraulics) {
      const ls = sim.hydraulics.lossSplit(), lag = e.time - sim.hydLast.t;
      for (const k of HYD_ENERGY) e['h_' + k] += (k === 'pump' ? ls.spp : ls[k]) * ls.q * lag;
    }
    const p = sim.params.energy;
    return { ...e, pump_efficiency: p.pump_efficiency, regen_efficiency: p.regen_efficiency,
             bit_diameter: sim.wellbore.hole_diameter, has_motor: !!sim.params.motor.enabled };
  }

  // driller's-display motor values (null without a motor)
  motorLatest(L) {
    const m = this.sim.params.motor;
    if (!m.enabled) return { motor: false, rpm_motor: null, dp_motor: null, motor_stalled: null, toolface_deg: null,
                             dp_no_load: null, dp_max: null, dp_stall: null };
    return { motor: true, rpm_motor: rpm(L.omega_motor), dp_motor: L.dp_motor, motor_stalled: L.motor_stalled,
             toolface_deg: L.toolface * 180 / Math.PI, dp_no_load: m.dp_no_load, dp_max: m.dp_max, dp_stall: m.dp_stall };
  }

  info() {
    const i = this.sim.info(), c = this.sim.controller;
    return { ...i, first_mode_hz: this.sim.modalFrequencies()[1], coupling: this.sim.params.friction.along_string_coupling,
             kp: c.kp, ki: c.ki, impedance: c.impedance };
  }

  csv() {
    return [REC_CHANNELS.join(',')].concat(this.sim.rec.map(r => r.join(','))).join('\n') + '\n';
  }
}
