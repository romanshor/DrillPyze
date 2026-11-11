// Wellbore temperature (spec section 10): single-phase energy balance of Wang et al. (2023)
// with radial conduction in the formation and a growing bottom cell.
import { G } from './units.js';
import { gradient } from './hydraulics.js';

export const THERMAL_DEFAULTS = { enabled: true, dx: 30, cut_factor: 1.5, n_radial: 12, far_field: 5, steel_conductivity: 45,
                                  casing_thickness: 0.012, flow_rate: 0.0347, inlet_temperature: 30, max_dt: 2, gs_iterations: 3 };

const geoT = (f, z) => f.surface_temperature + f.geothermal_gradient * z;
const densityAt = (fl, T, p) => (1 - (T - fl.reference_temperature) * fl.thermal_expansion) * fl.density + p / fl.sound_speed ** 2;

// convection coefficient (10.3)
function hConv(Q, A, Dh, fl) {
  if (A <= 0 || Dh <= 0) return 0;
  const k = fl.thermal_conductivity, mu = fl.viscosity, cp = fl.heat_capacity, rho = fl.density;
  const Re = rho * (Math.abs(Q) / A) * Dh / mu, Pr = mu * cp / k;
  const Nu = Re > 2300 ? 0.023 * Re ** 0.8 * Pr ** 0.4 : 3.66;
  return Math.max(Nu, 3.66) * k / Dh;
}

// backward-Euler radial conduction in one formation layer (10.5)
function solveLayer(T, rcVol, g, gIn, Tin, gOut, Tfar, dt) {
  const m = T.length, a = new Float64Array(m), b = new Float64Array(m), c = new Float64Array(m), d = new Float64Array(m);
  for (let i = 0; i < m; i++) {
    const C = rcVol[i] / dt;
    b[i] = C; d[i] = C * T[i];
    if (i === 0) { b[i] += gIn; d[i] += gIn * Tin; } else { a[i] = -g[i - 1]; b[i] += g[i - 1]; }
    if (i === m - 1) { b[i] += gOut; d[i] += gOut * Tfar; } else { c[i] = -g[i]; b[i] += g[i]; }
  }
  for (let i = 1; i < m; i++) { const w = a[i] / b[i - 1]; b[i] -= w * c[i - 1]; d[i] -= w * d[i - 1]; }
  T[m - 1] = d[m - 1] / b[m - 1];
  for (let i = m - 2; i >= 0; i--) T[i] = (d[i] - c[i] * T[i + 1]) / b[i];
}

export class WellboreThermal {
  constructor(wellbore, drillstring, fluid, formation, params, bitDepth, holeDepth) {
    Object.assign(this, { wellbore, drillstring, fluid, formation, p: params });
    this.time = 0; this.bitDepth = bitDepth; this.holeDepth = holeDepth;
    this.flowRate = params.flow_rate; this.inletTemperature = params.inlet_temperature;
    this.buildGrid(holeDepth);
    this.geometry(bitDepth, true);
    this.Tds = Float64Array.from(this.Tfar); this.Tann = Float64Array.from(this.Tfar);
    this.Tf = this.Tfar.map(t => new Float64Array(params.n_radial).fill(t));
    this.pressure();
  }
  // 10.2
  buildGrid(hole) {
    const dx = this.p.dx, n = Math.max(1, Math.floor(hole / dx));
    const e = Array.from({ length: n + 1 }, (_, i) => i * dx);
    if (hole - e[n] > 1e-9) { if (hole - e[n] < 0.5 * dx && n > 1) e[n] = hole; else e.push(hole); }
    this.xTop = e.slice(0, -1); this.dx = e.slice(1).map((x, i) => x - e[i]);
  }
  radial(rw) {
    const m = this.p.n_radial, f = this.formation, k = f.rock_conductivity;
    const fr = Array.from({ length: m + 1 }, (_, i) => rw * ((rw + this.p.far_field) / rw) ** (i / m));
    const rc = fr.slice(0, -1).map((x, i) => Math.sqrt(x * fr[i + 1]));
    const vol = rc.map((_, i) => Math.PI * (fr[i + 1] ** 2 - fr[i] ** 2) * f.rock_density * f.rock_heat_capacity);
    const g = rc.slice(0, -1).map((r, i) => 2 * Math.PI * k / Math.log(rc[i + 1] / r));
    return { vol, g, gOut: 2 * Math.PI * k / Math.log(fr[m] / rc[m - 1]), rIn: Math.log(rc[0] / rw) / (2 * Math.PI * k) };
  }
  // 10.3
  geometry(bitDepth, force = false) {
    const n = this.dx.length;
    const key = `${n}|${Math.round(this.dx[n - 1] * 100) / 100}|${this.flowRate}`;
    if (!force && this.geoKey === key && Math.abs(bitDepth - this.geoBit) < 0.25) return;
    this.geoKey = key; this.geoBit = bitDepth;
    const wb = this.wellbore, fl = this.fluid, p = this.p, sv = wb.survey;
    const mid = this.xTop.map((x, j) => x + 0.5 * this.dx[j]);
    this.tvd = sv.positions(mid)[2];
    this.Tfar = this.tvd.map(z => geoT(this.formation, z));
    const comps = this.drillstring.components.slice().reverse();
    const lens = comps.map((c, i) => (i === 0 ? bitDepth - this.drillstring.bhaLength : c.length));
    const bnds = [0]; for (const L of lens) bnds.push(bnds[bnds.length - 1] + L);
    const ks = p.steel_conductivity, Q = this.flowRate, rcp = fl.density * fl.heat_capacity;
    const A = (d) => Math.PI / 4 * d * d;
    Object.assign(this, { Ads: new Float64Array(n), Aann: new Float64Array(n), Gds: new Float64Array(n), Gw: new Float64Array(n),
                          Dwall: new Float64Array(n), rcVol: [], gRad: [], gOut: new Float64Array(n) });
    const cache = new Map();
    for (let j = 0; j < n; j++) {
      let sOd2 = 0, sId2 = 0, sOd = 0, sId = 0;
      for (let k = 0; k < 8; k++) {
        const x = this.xTop[j] + (k + 0.5) / 8 * this.dx[j];
        let ic = 0; for (let i = 0; i < bnds.length; i++) if (bnds[i] <= x) ic = i;   // searchsorted right - 1
        ic = Math.min(Math.max(ic, 0), comps.length - 1);
        const od = x <= bitDepth ? comps[ic].od : 0, id = x <= bitDepth ? comps[ic].id : 0;
        sOd2 += A(od); sId2 += A(id); sOd += od; sId += id;
      }
      const cased = mid[j] <= wb.casing_shoe;
      const Dw = cased ? (wb.casing_id || wb.hole_diameter) : wb.hole_diameter;
      this.Dwall[j] = Dw;
      this.Ads[j] = sId2 / 8; this.Aann[j] = A(Dw) - sOd2 / 8;
      const Di = sId / 8, Do = sOd / 8;
      const hA = hConv(Q, this.Aann[j], Dw - Do, fl);
      if (Do > 0 && Di > 0) {
        const hI = hConv(Q, this.Ads[j], Di, fl);
        const Rds = 1 / (hI * Math.PI * Di) + Math.log(Do / Di) / (2 * Math.PI * ks) + 1 / (hA * Math.PI * Do);
        this.Gds[j] = this.dx[j] / Rds;
      }
      const rw = Dw / 2 + (cased ? p.casing_thickness : 0);
      if (!cache.has(rw)) cache.set(rw, this.radial(rw));
      const r = cache.get(rw);
      let Rw = 1 / (hA * Math.PI * Dw) + r.rIn;
      if (cased) Rw += Math.log((Dw + 2 * p.casing_thickness) / Dw) / (2 * Math.PI * ks);
      this.Gw[j] = this.dx[j] / Rw;
      this.rcVol.push(r.vol.map(v => v * this.dx[j]));
      this.gRad.push(r.g.map(v => v * this.dx[j]));
      this.gOut[j] = r.gOut * this.dx[j];
    }
    this.Cds = this.Ads.map((a, j) => rcp * a * this.dx[j]);
    this.Cann = this.Aann.map((a, j) => rcp * a * this.dx[j]);
    let jb = 0; for (let j = 0; j < n; j++) if (this.xTop[j] <= bitDepth) jb = j;
    this.jb = Math.min(jb, n - 1); this.nb = n - 1;
  }
  // 10.4
  grow(h) {
    const n = this.dx.length, dh = h - (this.xTop[n - 1] + this.dx[n - 1]);
    if (dh <= 0) return;
    this.dx[n - 1] += dh;
    const Tb = geoT(this.formation, this.wellbore.survey.tvdAt(h)), w = Math.min(1, dh / this.dx[n - 1]);
    const tf = this.Tf[n - 1]; for (let i = 0; i < tf.length; i++) tf[i] += w * (Tb - tf[i]);
    if (this.dx[n - 1] > this.p.cut_factor * this.p.dx) {
      const half = 0.5 * this.dx[n - 1];
      this.dx[n - 1] = half;
      this.xTop.push(this.xTop[n - 1] + half); this.dx.push(half);
      this.Tds = Float64Array.from([...this.Tds, this.Tds[n - 1]]);
      this.Tann = Float64Array.from([...this.Tann, this.Tann[n - 1]]);
      this.Tf.push(Float64Array.from(tf));
    }
  }
  // 10.5 (one backward-Euler step with Gauss-Seidel coupling)
  step(dt, F, Tin, qc, Tc) {
    const n = this.Tann.length, jb = this.jb, nb = this.nb, { Tds, Tann, Tf, Cds, Cann, Gds, Gw } = this;
    const TdsO = Tds.slice(), TannO = Tann.slice(), TfO = Tf.map(r => r.slice());
    for (let it = 0; it < this.p.gs_iterations; it++) {
      let up = Tin;
      for (let j = 0; j <= jb; j++) {
        const Cd = Cds[j] / dt, num = Cd * TdsO[j] + F * up + Gds[j] * Tann[j], den = Cd + F + Gds[j];
        Tds[j] = den > 0 ? num / den : TdsO[j];
        up = Tds[j];
      }
      let below = Tds[jb];
      for (let j = n - 1; j >= 0; j--) {
        const Ca = Cann[j] / dt, Fj = j <= jb ? F : 0, gds = j <= jb ? Gds[j] : 0;
        let num = Ca * TannO[j] + gds * Tds[j] + Gw[j] * Tf[j][0], den = Ca + gds + Gw[j];
        if (Fj > 0) { num += Fj * below; den += Fj; }
        if (j === nb) { num += qc * Tc; den += qc; }
        Tann[j] = num / den;
        if (j <= jb) below = Tann[j];
      }
      for (let j = 0; j < n; j++) {
        Tf[j].set(TfO[j]);
        solveLayer(Tf[j], this.rcVol[j], this.gRad[j], Gw[j], Tann[j], this.gOut[j], this.Tfar[j], dt);
      }
    }
  }
  advance(duration, flowRate, inletTemperature, bitDepth, holeDepth, rop) {
    this.flowRate = flowRate; this.inletTemperature = inletTemperature; this.bitDepth = bitDepth;
    const hole = Math.max(holeDepth, this.holeDepth), nsub = Math.max(1, Math.ceil(duration / this.p.max_dt)), dt = duration / nsub;
    const fl = this.fluid, fm = this.formation, h0 = this.holeDepth;
    for (let k = 0; k < nsub; k++) {
      const hk = h0 + (hole - h0) * (k + 1) / nsub;
      this.grow(hk);
      this.geometry(this.bitDepth, rop > 0);
      const F = fl.density * fl.heat_capacity * Math.abs(this.flowRate);
      const onBottom = this.bitDepth >= hk - 0.5, Ahole = Math.PI / 4 * this.wellbore.hole_diameter ** 2;
      const qc = onBottom ? fm.rock_density * fm.rock_heat_capacity * Ahole * Math.max(rop, 0) : 0;
      this.step(dt, F, this.inletTemperature, qc, geoT(fm, this.wellbore.survey.tvdAt(hk)));
      this.time += dt;
    }
    this.holeDepth = hole;
    this.pressure();
    return this.state();
  }
  // 10.6
  pressure() {
    const fl = this.fluid, n = this.dx.length, sv = this.wellbore.survey;
    const edges = [...this.xTop, this.xTop[n - 1] + this.dx[n - 1]], ze = sv.positions(edges)[2];
    let p = 0; const Q = Math.abs(this.flowRate);
    for (let j = 0; j < n; j++) {
      const rho = densityAt(fl, this.Tann[j], p);
      let dp = rho * G * (ze[j + 1] - ze[j]);
      if (Q > 0 && j <= this.jb && this.Aann[j] > 0) {
        const Dw = this.Dwall[j], Do = Math.sqrt(Math.max(Dw * Dw - 4 * this.Aann[j] / Math.PI, 0));
        const area = Math.PI / 4 * (Dw * Dw - (Do > 0 ? Do * Do : 0));
        dp += gradient(Q / area, rho, fl.viscosity, 0, Dw, Do > 0 ? Do : 0) * this.dx[j];
      }
      p += dp;
    }
    this.pBh = p;
    this.ecd = ze[n] > 0 ? p / (G * ze[n]) : fl.density;
  }
  state() {
    return { time: this.time, t_bit: this.Tann[this.jb], t_ds_bit: this.Tds[this.jb], t_out: this.Tann[0], t_in: this.inletTemperature,
             t_rock: geoT(this.formation, this.wellbore.survey.tvdAt(this.holeDepth)), p_bh: this.pBh, ecd: this.ecd, flow_rate: this.flowRate };
  }
  profiles() {
    return { md: this.xTop.map((x, j) => x + 0.5 * this.dx[j]), tvd: Array.from(this.tvd), t_ds: Array.from(this.Tds),
             t_ann: Array.from(this.Tann), t_wall: this.Tf.map(r => r[0]), t_far: Array.from(this.Tfar), has_pipe: Array.from(this.Ads, a => a > 0) };
  }
}
