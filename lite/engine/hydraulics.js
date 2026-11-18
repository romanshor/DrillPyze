// Circulating-system hydraulics (spec section 8): steady losses and the transient
// single-phase drift-flux model, at the reference mud temperature (rho_T = rho_m).
import { G } from './units.js';
import { interp } from './survey.js';

const RE_LAMINAR = 2100, ANNULUS_FACTOR = 0.816;

export function rheology(fluid) {
  if (fluid.yield_point) return [fluid.plastic_viscosity != null ? fluid.plastic_viscosity : fluid.viscosity, fluid.yield_point];
  return [fluid.viscosity, 0];
}

// friction gradient dp/dL = a u + b sign(u) at speed au (spec 8.1); returns [a, b]
export function coeffs(au, rho, pv, ty, dOut, dIn) {
  let aLam, bLam, De, muA;
  if (dIn <= 0) {
    const D = dOut;
    aLam = 32 * pv / (D * D); bLam = 16 * ty / (3 * D); De = D;
    muA = pv + (au > 0 ? ty * D / (6 * au) : 1e30);
  } else {
    const hh = dOut - dIn;
    aLam = 48 * pv / (hh * hh); bLam = 6 * ty / hh; De = ANNULUS_FACTOR * hh;
    muA = pv + (au > 0 ? ty * hh / (8 * au) : 1e30);
  }
  const Re = rho * au * De / muA;
  if (Re < RE_LAMINAR) return [aLam, bLam];
  const f = 0.0791 / Math.pow(rho * au * De / pv, 0.25);
  return [2 * f * rho * au / De, 0];
}
export function gradient(u, rho, pv, ty, dOut, dIn) {
  const au = Math.abs(u);
  if (au === 0) return 0;
  const [a, b] = coeffs(au, rho, pv, ty, dOut, dIn);
  const g = a * au + b;
  return u > 0 ? g : -g;
}
export const bitPressureDrop = (Q, rho, tfa, cd) => rho * Q * Math.abs(Q) / (2 * cd * cd * tfa * tfa);
export const surfaceLoss = (Q, c, m) => c * Math.pow(Math.abs(Q), m);

// backward-Euler face update with Picard iterations (spec 8.3, "face rule")
function face(Gf, dt, invI, drive, ru, Au, Lu, duo, dui, rd, Ad, Ld, ddo, ddi, pv, ty, vSlip, kSlip, extraR, iters) {
  const rhoF = 0.5 * (ru + rd), AF = 0.5 * (Au + Ad), Gs = kSlip * vSlip * rhoF * AF;
  let Gk = Gf;
  for (let it = 0; it < iters; it++) {
    let R = 0, B = 0;
    if (Lu > 0) { const [a, b] = coeffs(Math.abs(Gk / (ru * Au) - kSlip * vSlip), ru, pv, ty, duo, dui); R += 0.5 * Lu * a / (ru * Au); B += 0.5 * Lu * b; }
    if (Ld > 0) { const [a, b] = coeffs(Math.abs(Gk / (rd * Ad) - kSlip * vSlip), rd, pv, ty, ddo, ddi); R += 0.5 * Ld * a / (rd * Ad); B += 0.5 * Ld * b; }
    R += extraR * Math.abs(Gk);
    const c = dt * invI, free = Gf + c * (drive + R * Gs), den = 1 + c * R, rel = free - den * Gs;
    if (Math.abs(rel) <= c * B) Gk = Gs;
    else Gk = (free - c * B * (rel > 0 ? 1 : -1)) / den;
  }
  return Gk;
}

const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => (n === 1 ? a : a + (b - a) * i / (n - 1)));
const diff = (e) => e.slice(1).map((x, i) => x - e[i]);

export class WellboreHydraulics {
  constructor(wellbore, grid, fluid, params, tfa, cd, bitDepth, holeDepth) {
    this.wellbore = wellbore; this.fluid = fluid; this.p = params; this.tfa = tfa; this.cd = cd; this.dpMotor = 0;
    [this.pv, this.ty] = rheology(fluid);
    this.c = fluid.sound_speed; this.dx = params.dx;
    this.sEdges = [0]; this.cumOd = [0]; this.cumBore = [0];
    for (let j = 0; j < grid.n; j++) {
      this.sEdges.push(this.sEdges[j] + grid.dx[j]);
      this.cumOd.push(this.cumOd[j] + Math.PI / 4 * grid.od[j] ** 2 * grid.dx[j]);
      this.cumBore.push(this.cumBore[j] + Math.PI / 4 * grid.id[j] ** 2 * grid.dx[j]);
    }
    this.Ls = this.sEdges[grid.n];
    this.bitDepth = bitDepth; this.holeDepth = holeDepth;
    this.eA = this._edges(holeDepth); this.eD = this._edges(bitDepth);
    this.TD = null; this.TA = null; this.Tsrc = null;      // mud temperatures from the thermal model (spec 10.7)
    this._geometry();
    this.VD = this._Vd.slice(); this.VA = this._Va.slice();
    this.qPump = 0; this.vPipe = 0; this.acc = [0, 0];
    this.setStatic();
  }
  _edges(depth) { const n = Math.max(1, Math.round(depth / this.dx)); return linspace(0, depth, n + 1); }
  _stringIntegral(cum, md) {
    const top = this.bitDepth - this.Ls, clip = (x) => Math.min(Math.max(x, 0), this.Ls);
    return interp(clip(md - top), this.sEdges, cum) - interp(clip(-top), this.sEdges, cum);
  }
  _geometry() {
    const sv = this.wellbore.survey, eA = this.eA, eD = this.eD;
    this.LA = diff(eA); this.LD = diff(eD);
    const cA = eA.slice(0, -1).map((x, i) => 0.5 * (x + eA[i + 1])), cD = eD.slice(0, -1).map((x, i) => 0.5 * (x + eD[i + 1]));
    this.zA = sv.positions(cA)[2]; this.zD = sv.positions(cD)[2];
    this.tvdBottom = sv.tvdAt(this.holeDepth);
    const Dw = cA.map(m => this.wellbore.wallDiameter(m));
    const Aod = this.LA.map((L, i) => (this._stringIntegral(this.cumOd, eA[i + 1]) - this._stringIntegral(this.cumOd, eA[i])) / L);
    const Abore = this.LD.map((L, i) => (this._stringIntegral(this.cumBore, eD[i + 1]) - this._stringIntegral(this.cumBore, eD[i])) / L);
    this.DoA = Dw;
    this.DiA = Aod.map(a => { const d = Math.sqrt(Math.max(a, 0) * 4 / Math.PI); return d < 1e-3 ? 0 : d; });
    this.DD = Abore.map(a => Math.sqrt(Math.max(a, 1e-8) * 4 / Math.PI));
    this._Va = Dw.map((D, i) => (Math.PI / 4 * D * D - Aod[i]) * this.LA[i]);
    this._Vd = Abore.map((a, i) => Math.max(a, 1e-8) * this.LD[i]);
    // annulus cell containing the bit: last edge <= bit depth, clipped
    let k = 0;
    for (let i = 0; i < eA.length; i++) if (eA[i] <= this.bitDepth) k = i;
    this.kab = Math.max(0, Math.min(k, this.LA.length - 1));
    this._rhoT();
  }
  // rho_T = (1 - (T - T_ref) beta) rho_m, T = T_ref without the thermal model (spec 8.2)
  _rhoT() {
    const fl = this.fluid, f = (T) => (1 - (T - fl.reference_temperature) * fl.thermal_expansion) * fl.density;
    this.rTD = this.LD.map((_, i) => f(this.TD ? this.TD[i] : fl.reference_temperature));
    this.rTA = this.LA.map((_, i) => f(this.TA ? this.TA[i] : fl.reference_temperature));
  }
  setTemperatures(mdDs, tDs, mdAnn, tAnn) {
    this.Tsrc = [mdDs, tDs, mdAnn, tAnn];
    this._interpT();
  }
  _interpT() {
    if (!this.Tsrc) return;
    const [mdDs, tDs, mdAnn, tAnn] = this.Tsrc, mid = (e) => e.slice(0, -1).map((x, i) => 0.5 * (x + e[i + 1]));
    this.TD = mid(this.eD).map(m => interp(m, mdDs, tDs));
    this.TA = mid(this.eA).map(m => interp(m, mdAnn, tAnn));
    this._rhoT();
  }
  _rhoOfP(p, rT) { return rT + p / (this.c * this.c); }

  _regrid(bitNew, holeNew) {
    const dx = this.dx, cut = this.p.cut_factor;
    if (holeNew > this.eA[this.eA.length - 1]) {
      const e = this.eA, Lnew = holeNew - e[e.length - 2];
      if (Lnew > cut * dx) {
        const mid = e[e.length - 2] + 0.5 * Lnew;
        if (mid < e[e.length - 1]) {
          const f = (mid - e[e.length - 2]) / (e[e.length - 1] - e[e.length - 2]);
          const m = this.mA.pop(), V = this.VA.pop();
          this.mA.push(f * m, (1 - f) * m); this.VA.push(f * V, (1 - f) * V);
          this.GA.push(0.5 * this.GA[this.GA.length - 1]);
          this.eA = [...e.slice(0, -1), mid, e[e.length - 1]];
          if (this.TA) this.TA.push(this.TA[this.TA.length - 1]);
        }
      }
    }
    const e = this.eD, Lnew = bitNew - e[e.length - 2];
    if (Lnew > cut * dx) {
      const mid = Math.min(e[e.length - 2] + 0.5 * Lnew, e[e.length - 1] - 1e-6);
      const f = (mid - e[e.length - 2]) / (e[e.length - 1] - e[e.length - 2]);
      const m = this.mD.pop(), V = this.VD.pop();
      this.mD.push(f * m, (1 - f) * m); this.VD.push(f * V, (1 - f) * V);
      this.GD.push(0.5 * (this.GD[this.GD.length - 1] + this.Gb));
      this.eD = [...e.slice(0, -1), mid, e[e.length - 1]];
      if (this.TD) this.TD.push(this.TD[this.TD.length - 1]);
    } else if (Lnew < 0.5 * dx && this.mD.length > 1) {
      const m = this.mD.pop(), V = this.VD.pop();
      this.mD[this.mD.length - 1] += m; this.VD[this.VD.length - 1] += V;
      this.GD.pop();
      this.eD = [...e.slice(0, -2), e[e.length - 1]];
      if (this.TD) this.TD.pop();
    }
  }

  setStatic() { this.setSteady(0); }
  setSteady(Q) {
    this._geometry();
    const rhoIn = this.rTD[0], Gm = rhoIn * Q, { pv, ty } = this, iters = 50, bp = this.p.back_pressure;
    const nA = this.LA.length, nD = this.LD.length, kab = this.kab;
    const AA = this._Va.map((V, i) => V / this.LA[i]), AD = this._Vd.map((V, i) => V / this.LD[i]);
    const pA = new Array(nA).fill(0);
    let p = bp + this.rTA[0] * G * this.zA[0];
    for (let it = 0; it < iters; it++) {
      const rho = this._rhoOfP(p, this.rTA[0]);
      const fr = 0.5 * this.LA[0] * gradient(-Gm / (rho * AA[0]), rho, pv, ty, this.DoA[0], this.DiA[0]);
      p = bp + rho * G * this.zA[0] - fr;
    }
    pA[0] = p;
    for (let k = 1; k < nA; k++) {
      const Gk = k <= kab ? -Gm : 0;
      let pk = pA[k - 1] + this.rTA[k] * G * (this.zA[k] - this.zA[k - 1]);
      for (let it = 0; it < iters; it++) {
        const ru = this._rhoOfP(pA[k - 1], this.rTA[k - 1]), rd = this._rhoOfP(pk, this.rTA[k]);
        const fr = 0.5 * this.LA[k - 1] * gradient(Gk / (ru * AA[k - 1]), ru, pv, ty, this.DoA[k - 1], this.DiA[k - 1])
                 + 0.5 * this.LA[k] * gradient(Gk / (rd * AA[k]), rd, pv, ty, this.DoA[k], this.DiA[k]);
        pk = pA[k - 1] + 0.5 * (ru + rd) * G * (this.zA[k] - this.zA[k - 1]) - fr;
      }
      pA[k] = pk;
    }
    const pD = new Array(nD).fill(0), jd = nD - 1;
    let pk = pA[kab];
    for (let it = 0; it < iters; it++) {
      const rd_ = this._rhoOfP(pk, this.rTD[jd]), ra_ = this._rhoOfP(pA[kab], this.rTA[kab]), rf = 0.5 * (rd_ + ra_);
      const fr = 0.5 * this.LD[jd] * gradient(Gm / (rd_ * AD[jd]), rd_, pv, ty, this.DD[jd], 0)
               + 0.5 * this.LA[kab] * gradient(Gm / (ra_ * AA[kab]), ra_, pv, ty, this.DoA[kab], this.DiA[kab]);
      pk = pA[kab] - rf * G * (this.zA[kab] - this.zD[jd]) + bitPressureDrop(Gm / rf, rf, this.tfa, this.cd) + fr + this.dpMotor;
    }
    pD[jd] = pk;
    for (let k = jd; k > 0; k--) {
      let pu = pD[k] - this.rTD[k] * G * (this.zD[k] - this.zD[k - 1]);
      for (let it = 0; it < iters; it++) {
        const ru = this._rhoOfP(pu, this.rTD[k - 1]), rd = this._rhoOfP(pD[k], this.rTD[k]);
        const fr = 0.5 * this.LD[k - 1] * gradient(Gm / (ru * AD[k - 1]), ru, pv, ty, this.DD[k - 1], 0)
                 + 0.5 * this.LD[k] * gradient(Gm / (rd * AD[k]), rd, pv, ty, this.DD[k], 0);
        pu = pD[k] - 0.5 * (ru + rd) * G * (this.zD[k] - this.zD[k - 1]) + fr;
      }
      pD[k - 1] = pu;
    }
    this.VD = this._Vd.slice(); this.VA = this._Va.slice();
    this.mD = pD.map((p_, i) => this._rhoOfP(p_, this.rTD[i]) * this.VD[i]);
    this.mA = pA.map((p_, i) => this._rhoOfP(p_, this.rTA[i]) * this.VA[i]);
    this.GD = new Array(nD).fill(Gm);
    this.GA = pA.map((_, i) => (i <= kab ? -Gm : 0));
    this.Gb = Gm;
    this.qPump = Q;
  }

  advance(duration, q0, q1, bitDepth, holeDepth, dpMotor = 0) {
    if (duration <= 0) return this.state();
    this.dpMotor = dpMotor;      // mud motor differential (Pa), in series with the bit nozzles
    holeDepth = Math.max(holeDepth, this.holeDepth);
    this._regrid(bitDepth, holeDepth);
    this.vPipe = (bitDepth - this.bitDepth) / duration;
    this.bitDepth = bitDepth; this.holeDepth = holeDepth;
    this.eA = [...this.eA.slice(0, -1), holeDepth];
    this.eD = [...this.eD.slice(0, -1), bitDepth];
    this._interpT();
    this._geometry();
    const VD0 = this.VD.length === this._Vd.length ? this.VD : this._Vd.slice();
    const VA0 = this.VA.length === this._Va.length ? this.VA : this._Va.slice();
    const Lmin = Math.min(...this.LA, ...this.LD);
    const nsub = Math.max(1, Math.ceil(duration / (this.p.cfl * Lmin / this.c)));
    const dt = duration / nsub, c2 = this.c * this.c;
    const nozzleR = 1 / (2 * this.cd ** 2 * this.tfa ** 2);
    const { pv, ty } = this, kAnn = this.p.clinging_constant, vp = this.vPipe, iters = this.p.iterations, pBack = this.p.back_pressure;
    const nD = this.mD.length, nA = this.mA.length, kab = this.kab, rhoIn = this.rTD[0];
    const LD = this.LD, LA = this.LA, zD = this.zD, zA = this.zA, DD = this.DD, DoA = this.DoA, DiA = this.DiA;
    const mD = this.mD, mA = this.mA, GD = this.GD, GA = this.GA;
    const pD = new Float64Array(nD), pA = new Float64Array(nA), rD = new Float64Array(nD), rA = new Float64Array(nA),
          AD = new Float64Array(nD), AA = new Float64Array(nA);
    let Gb = this.Gb;
    for (let s = 0; s < nsub; s++) {
      const w1 = (s + 1) / nsub;
      GD[0] = rhoIn * (q0 + (q1 - q0) * w1);
      for (let j = 0; j < nD; j++) { const V = VD0[j] + (this._Vd[j] - VD0[j]) * w1; AD[j] = V / LD[j]; rD[j] = mD[j] / V; pD[j] = c2 * (rD[j] - this.rTD[j]); }
      for (let j = 0; j < nA; j++) { const V = VA0[j] + (this._Va[j] - VA0[j]) * w1; AA[j] = V / LA[j]; rA[j] = mA[j] / V; pA[j] = c2 * (rA[j] - this.rTA[j]); }
      for (let k = 1; k < nD; k++) {
        const drive = pD[k - 1] - pD[k] + 0.5 * (rD[k - 1] + rD[k]) * G * (zD[k] - zD[k - 1]);
        const invI = 1 / (0.5 * (LD[k - 1] / AD[k - 1] + LD[k] / AD[k]));
        GD[k] = face(GD[k], dt, invI, drive, rD[k - 1], AD[k - 1], LD[k - 1], DD[k - 1], 0, rD[k], AD[k], LD[k], DD[k], 0, pv, ty, vp, 1, 0, iters);
      }
      const jd = nD - 1;
      {
        const drive = pD[jd] - pA[kab] + 0.5 * (rD[jd] + rA[kab]) * G * (zA[kab] - zD[jd]) - this.dpMotor;
        const invI = 1 / (0.5 * (LD[jd] / AD[jd] + LA[kab] / AA[kab]));
        Gb = face(Gb, dt, invI, drive, rD[jd], AD[jd], LD[jd], DD[jd], 0, rA[kab], AA[kab], LA[kab], DoA[kab], DiA[kab],
                  pv, ty, 0, 0, nozzleR / (0.5 * (rD[jd] + rA[kab])), iters);
      }
      for (let k = 1; k < nA; k++) {
        const drive = pA[k - 1] - pA[k] + 0.5 * (rA[k - 1] + rA[k]) * G * (zA[k] - zA[k - 1]);
        const invI = 1 / (0.5 * (LA[k - 1] / AA[k - 1] + LA[k] / AA[k]));
        const ks = DiA[k] > 0 && DiA[k - 1] > 0 ? kAnn : 0;
        GA[k] = face(GA[k], dt, invI, drive, rA[k - 1], AA[k - 1], LA[k - 1], DoA[k - 1], DiA[k - 1], rA[k], AA[k], LA[k], DoA[k], DiA[k], pv, ty, vp, ks, 0, iters);
      }
      {
        const drive = pBack - pA[0] + rA[0] * G * zA[0], invI = 1 / (0.5 * LA[0] / AA[0]);
        GA[0] = face(GA[0], dt, invI, drive, rA[0], AA[0], 0, DoA[0], DiA[0], rA[0], AA[0], LA[0], DoA[0], DiA[0], pv, ty, vp, DiA[0] > 0 ? kAnn : 0, 0, iters);
      }
      this.acc[0] += dt * GD[0];
      this.acc[1] -= dt * GA[0];
      for (let j = 0; j < nD - 1; j++) mD[j] += dt * (GD[j] - GD[j + 1]);
      mD[jd] += dt * (GD[jd] - Gb);
      for (let j = 0; j < nA - 1; j++) mA[j] += dt * (GA[j] - GA[j + 1]);
      mA[nA - 1] += dt * GA[nA - 1];
      mA[kab] += dt * Gb;
    }
    this.Gb = Gb;
    this.VD = this._Vd.slice(); this.VA = this._Va.slice();
    this.qPump = q1;
    return this.state();
  }

  state() {
    const c2 = this.c * this.c, G0 = this.GD[0];
    const pD0 = c2 * (this.mD[0] / this.VD[0] - this.rTD[0]);
    const rd0 = this.mD[0] / this.VD[0], A0 = this.VD[0] / this.LD[0], u0 = G0 / (rd0 * A0);
    const fr0 = 0.5 * this.LD[0] * gradient(u0 - this.vPipe, rd0, this.pv, this.ty, this.DD[0], 0);
    const pTop = pD0 - rd0 * G * this.zD[0] + fr0;
    const q = G0 / this.rTD[0];
    const spp = pTop + surfaceLoss(q, this.p.surface_loss_coefficient, this.p.surface_loss_exponent);
    const nA = this.mA.length, rb = this.mA[nA - 1] / this.VA[nA - 1];
    const pAl = c2 * (rb - this.rTA[nA - 1]);
    const bhp = pAl + rb * G * (this.tvdBottom - this.zA[nA - 1]);
    const ra0 = this.mA[0] / this.VA[0], nD = this.mD.length, kab = this.kab;
    const rbF = 0.5 * (this.mD[nD - 1] / this.VD[nD - 1] + this.mA[kab] / this.VA[kab]);
    return { spp, p_bh: bhp, ecd: this.tvdBottom > 0 ? bhp / (G * this.tvdBottom) : this.fluid.density,
             dp_bit: bitPressureDrop(this.Gb / rbF, rbF, this.tfa, this.cd),
             q_in: q, q_out: -this.GA[0] / ra0, q_bit: this.Gb / (0.5 * (rd0 + ra0)),
             v_pipe: this.vPipe, surface_loss: spp - pTop, dp_motor: this.dpMotor,
             pit_gain: -(this.acc[0] - this.acc[1]) / this.fluid.density };
  }

  // where the standpipe pressure goes (Pa, spec 8.4): surface lines, friction inside the
  // string, the mud motor, the bit nozzles, friction in the annulus, and the rest
  lossSplit() {
    const c2 = this.c * this.c, hs = this.state(), nD = this.mD.length, k = this.kab;
    let hydD = 0, hydA = 0, z = 0;
    for (let j = 0; j < nD; j++) { hydD += this.mD[j] / this.VD[j] * G * (this.zD[j] - z); z = this.zD[j]; }
    z = 0;
    for (let j = 0; j <= k; j++) { hydA += this.mA[j] / this.VA[j] * G * (this.zA[j] - z); z = this.zA[j]; }
    const pDb = c2 * (this.mD[nD - 1] / this.VD[nD - 1] - this.rTD[nD - 1]);
    const pAk = c2 * (this.mA[k] / this.VA[k] - this.rTA[k]);
    const pTop = hs.spp - hs.surface_loss;
    const out = { q: hs.q_in, spp: hs.spp, surface: hs.surface_loss, string: pTop + hydD - pDb, motor: this.dpMotor,
                  nozzles: hs.dp_bit, annulus: pAk - hydA - this.p.back_pressure };
    out.other = out.spp - (out.surface + out.string + out.motor + out.nozzles + out.annulus);
    return out;
  }
}
