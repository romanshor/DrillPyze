// Time step of the coupled axial-torsional string model (spec section 5).
// Operates in place on the simulation object `s` (arrays and scalar fields).

export const BIT_MODELS = ['coulomb', 'cutting', 'regenerative'];
export const REC_CHANNELS = ['time', 'omega_td', 'motor_torque', 'pipe_torque', 'omega_bit', 'tob', 'v_block',
                             'hookload', 'v_bit', 'wob', 'bit_depth', 'hole_depth', 'omega_sp'];
// ramp rows (spec 5.1)
export const R_WSP = 0, R_VBLK = 1, R_TAUM = 2, R_WTD = 3, R_WOB = 4;

const sign = (x) => (x > 0 ? 1 : x < 0 ? -1 : 0);

export function rampValue(r, t) {
  const [v0, t0, T, tgt] = r;
  if (T <= 0) return tgt;
  const f = (t - t0) / T;
  if (f >= 1) return tgt;
  if (f <= 0) return v0;
  return v0 + (tgt - v0) * f;
}

// quasi-static bit, coulomb and cutting models (spec 5.5)
function bitQS(alT, alA, Zt, Za, onBottom, wPrev, stuckPrev, p) {
  if (!onBottom || alA <= 0) return [alT, 0, alA, 0, 0];
  const { levS, levK, krop: kr, w0, wc } = p;
  let w = alT > 0 ? alT : 0, wob = 0, v = 0, tob = 0, stuck = 0;
  for (let it = 0; it < 8; it++) {
    const qlin = Za * alA, aw = Math.abs(w);
    if (qlin <= w0 || kr * aw === 0) { wob = qlin; v = 0; }
    else { wob = (alA + kr * w0 * aw) / (1 / Za + kr * aw); v = kr * (wob - w0) * aw; }
    const tStop = Zt * alT;
    let wNew;
    if ((Math.abs(wPrev) <= wc || stuckPrev > 0.5) && Math.abs(tStop) <= levS * wob) { tob = tStop; wNew = 0; stuck = 1; }
    else { tob = sign(alT) * Math.min(levK * wob, Math.abs(tStop)); wNew = alT - tob / Zt; stuck = 0; }
    if (Math.abs(wNew - w) < 1e-12) { w = wNew; break; }
    w = wNew;
  }
  return [w, tob, v, wob, stuck];
}

// regenerative bit (spec 5.5)
function bitRegen(alT, alA, Zt, Za, u, dt, wPrev, stuckPrev, phi, surf, p) {
  const kc = p.kcut, n = p.blades, lev = p.levK, wc = p.wc, M = surf.length, dth = 2 * Math.PI / M;
  const sBl = new Float64Array(n);
  for (let j = 0; j < n; j++) {
    const b = ((Math.floor((phi + 2 * Math.PI * j / n) / dth) % M) + M) % M;
    sBl[j] = surf[b];
  }
  const uFree = u + dt * alA;
  let u1 = uFree;
  for (let it = 0; it < n + 2; it++) {
    let m = 0, ssum = 0;
    for (let j = 0; j < n; j++) if (sBl[j] < u1) { m++; ssum += sBl[j]; }
    let uNew;
    if (m === 0) uNew = uFree;
    else { const g = dt * kc / Za; uNew = (uFree + g * ssum) / (1 + g * m); }
    if (Math.abs(uNew - u1) < 1e-15) { u1 = uNew; break; }
    u1 = uNew;
  }
  let wob = 0;
  for (let j = 0; j < n; j++) if (u1 > sBl[j]) wob += kc * (u1 - sBl[j]);
  const v = (u1 - u) / dt, tStop = Zt * alT;
  if ((Math.abs(wPrev) <= wc || stuckPrev > 0.5) && Math.abs(tStop) <= lev * wob) return [0, tStop, v, wob, 1];
  const tob = sign(alT) * Math.min(lev * wob, Math.abs(tStop));
  return [alT - tob / Zt, tob, v, wob, 0];
}

function sweep(phiOld, phiNew, u, surf, n) {
  if (phiNew <= phiOld) return;
  const M = surf.length, dth = 2 * Math.PI / M;
  for (let j = 0; j < n; j++) {
    const off = 2 * Math.PI * j / n;
    const b0 = Math.floor((phiOld + off) / dth), b1 = Math.floor((phiNew + off) / dth);
    const last = Math.min(b1, b0 + M);
    for (let b = b0; b < last; b++) { const k = ((b % M) + M) % M; if (u > surf[k]) surf[k] = u; }
  }
}

// implicit 2-D Coulomb friction on one cell (spec 5.7)
function coulomb2d(va, vt, ma, mt, c) {
  const pa = ma * va, pt = mt * vt;
  if (pa * pa + pt * pt <= c * c) return [0, 0];
  if (c <= 0) return [va, vt];
  let s = Math.sqrt(va * va + vt * vt);
  for (let it = 0; it < 50; it++) {
    const da = ma * s + c, dt_ = mt * s + c, ra = pa / da, rt = pt / dt_;
    const g = ra * ra + rt * rt - 1, dg = -2 * (ra * ra * ma / da + rt * rt * mt / dt_);
    let sNew = s - g / dg;
    if (sNew < 0) sNew = 0;
    if (Math.abs(sNew - s) <= 1e-14 * (s + 1e-30)) { s = sNew; break; }
    s = sNew;
  }
  return [va * ma * s / (ma * s + c), vt * mt * s / (mt * s + c)];
}

// backward-Euler Kelvin-Voigt step, Thomas algorithm (spec 5.6)
function kvSolve(x, m, h, k, dt, cp, dp) {
  const n = x.length;
  for (let j = 0; j < n; j++) {
    const mh = m[j] * h[j] / dt, a = -k[j], b = mh + k[j] + k[j + 1], c = -k[j + 1], d = mh * x[j];
    if (j === 0) { cp[j] = c / b; dp[j] = d / b; }
    else { const den = b - a * cp[j - 1]; cp[j] = c / den; dp[j] = (d - a * dp[j - 1]) / den; }
  }
  x[n - 1] = dp[n - 1];
  for (let j = n - 2; j >= 0; j--) x[j] = dp[j] - cp[j] * x[j + 1];
}

/** Advance `nsteps` steps; recorded rows are pushed to s.rec (one array per channel). */
export function advance(s, nsteps) {
  const n = s.n, p = s.p, dt = p.dt, h = s.dx;
  const { Kt, mt, Zt, Ka, ma, Za, rc, wb, kx, ky, kz, zx, zy, zz, mus, muk, qsC, fn0,
          w, tau, v, qd, wf, tf, vf, qf, k0t, k0a, k1t, k1a, wbar, vbar, kvt, kva, surf } = s;
  const fy = s.yieldDragNow();
  const { wst, vst, wOld, vOld, cpr, dpr } = s.tmp;
  const wc = p.wc, vc = p.vc, I = p.itd, ctrl = s.ctrl, axialMode = s.axial, coupled = p.coupled,
        bitModel = p.bitModel, nBl = p.blades, axc = p.axc, weakExp = p.weakExp, useKv = p.kv;
  const aEn = p.ten > 0 ? Math.exp(-dt / p.ten) : 0, aLp = p.tlp > 0 ? Math.exp(-dt / p.tlp) : 0,
        aHp = p.thp > 0 ? Math.exp(-dt / p.thp) : 0;
  const R = s.ramps, sc = s.sc, rec = s.rec;

  for (let step = 0; step < nsteps; step++) {
    // ---- top drive input (5.1, 5.2)
    const wtd = sc.wtd, tNow = sc.t + dt;
    const wspT = rampValue(R[R_WSP], tNow), vBlkIn = rampValue(R[R_VBLK], tNow), taumExt = rampValue(R[R_TAUM], tNow),
          wtdExt = rampValue(R[R_WTD], tNow), wobSp = rampValue(R[R_WOB], tNow);
    if (ctrl === 0) {
      sc.wsp = wspT;
      sc.wfilt = aEn * sc.wfilt + (1 - aEn) * wtd;
      sc.lptau = aLp * sc.lptau + (1 - aLp) * sc.tau0;
      sc.hplp = aHp * sc.hplp + (1 - aHp) * sc.lptau;
      const e = sc.wsp - sc.wfilt - p.zc * (sc.lptau - sc.hplp);
      let u = p.kp * e + p.ki * (sc.int + e * dt);
      if (p.taumax > 0 && Math.abs(u) > p.taumax) u = p.taumax * sign(u);
      else sc.int += e * dt;
      sc.taum = u;
    } else if (ctrl === 1) { sc.taum = taumExt; sc.wsp = wspT; }
    const beta1 = w[0] - tau[0] / Zt[0];
    const wNewTd = ctrl === 2 ? wtdExt : (I * wtd / dt + sc.taum + Zt[0] * beta1) / (I / dt + Zt[0]);
    const tau0 = Zt[0] * (wNewTd - beta1);
    if (ctrl === 2) { sc.taum = tau0 + I * (wNewTd - wtd) / dt; sc.wsp = wtdExt; }
    sc.wtd = wNewTd; sc.tau0 = tau0; wf[0] = wNewTd; tf[0] = tau0;

    // ---- drawworks (5.3)
    if (axialMode === 1) {
      const err = wobSp - (p.hkldRef - sc.hkld);
      let vb = p.adKp * err + p.adKi * (sc.adInt + err * dt);
      if (vb > p.adVmax) vb = p.adVmax;
      else if (vb < -p.adVmax) vb = -p.adVmax;
      else sc.adInt += err * dt;
      sc.vblk = vb;
    } else {
      const acc = p.blkAcc, dv = vBlkIn - sc.vblk;
      if (acc > 0 && Math.abs(dv) > acc * dt) sc.vblk += acc * dt * sign(dv);
      else sc.vblk = vBlkIn;
    }
    const betaA1 = v[0] - qd[0] / Za[0];
    let q0 = Za[0] * (sc.vblk - betaA1), vTop = sc.vblk;
    if (-(q0 + p.qsTop) < 0) { q0 = -p.qsTop; vTop = betaA1 + q0 / Za[0]; }
    vf[0] = vTop; qf[0] = q0; sc.hkld = -(q0 + p.qsTop);

    // ---- interior faces (5.4)
    for (let j = 0; j < n - 1; j++) {
      let a = w[j] + tau[j] / Zt[j], b = w[j + 1] - tau[j + 1] / Zt[j + 1];
      let F = (a - b) / (1 / Zt[j] + 1 / Zt[j + 1]);
      wf[j + 1] = a - F / Zt[j]; tf[j + 1] = F;
      a = v[j] + qd[j] / Za[j]; b = v[j + 1] - qd[j + 1] / Za[j + 1];
      F = (a - b) / (1 / Za[j] + 1 / Za[j + 1]);
      vf[j + 1] = a - F / Za[j]; qf[j + 1] = F;
    }

    // ---- bit face (5.5)
    const alT = w[n - 1] + tau[n - 1] / Zt[n - 1], alA = v[n - 1] + qd[n - 1] / Za[n - 1];
    const res = bitModel === 2
      ? bitRegen(alT, alA, Zt[n - 1], Za[n - 1], sc.bitd, dt, sc.wbit, sc.stuck, sc.phib, surf, p)
      : bitQS(alT, alA, Zt[n - 1], Za[n - 1], sc.bitd >= sc.holed - 1e-9, sc.wbit, sc.stuck, p);
    wf[n] = res[0]; tf[n] = res[1]; vf[n] = res[2]; qf[n] = res[3];
    sc.wbit = res[0]; sc.tob = res[1]; sc.vbit = res[2]; sc.wob = res[3]; sc.stuck = res[4];

    // ---- cell update and linear damping (5.6)
    for (let j = 0; j < n; j++) {
      wOld[j] = w[j]; vOld[j] = v[j];
      tau[j] -= dt * Kt[j] / h[j] * (wf[j + 1] - wf[j]);
      const ws = w[j] - dt / (mt[j] * h[j]) * (tf[j + 1] - tf[j]);
      qd[j] -= dt * Ka[j] / h[j] * (vf[j + 1] - vf[j]);
      const vs = v[j] - dt / (ma[j] * h[j]) * (qf[j + 1] - qf[j]);
      wst[j] = (ws + dt * k1t[j] * wbar[j]) / (1 + dt * (k0t[j] + k1t[j]));
      vst[j] = (vs + dt * k1a[j] * vbar[j]) / (1 + dt * (k0a[j] + k1a[j]));
    }
    if (useKv) { kvSolve(wst, mt, h, kvt, dt, cpr, dpr); kvSolve(vst, ma, h, kva, dt, cpr, dpr); }

    // ---- along-string friction (5.7)
    for (let j = 0; j < n; j++) {
      const wo = wOld[j], vo = vOld[j], ws = wst[j], vs = vst[j], fyj = fy[j], r = rc[j];
      if (coupled) {
        const T = -(qsC[j] + qd[j]);
        const fx = T * kx[j] + wb[j] * zx[j], fyy = T * ky[j] + wb[j] * zy[j], fz = T * kz[j] + wb[j] * zz[j];
        const fn = Math.sqrt(fx * fx + fyy * fyy + fz * fz);
        const muSa = muk[j] + axc * (mus[j] - muk[j]), muSt = mus[j];
        const vtS = ws * r, spd = Math.sqrt(vs * vs + vtS * vtS);
        const ea = spd > 0 ? vs / spd : 0, et = spd > 0 ? vtS / spd : 1;
        const sOld = Math.sqrt((vo / vc) ** 2 + (wo / wc) ** 2);
        const faS = muSa * fn + fyj, ftS = muSt * fn + fyj;
        const fSd = (spd > 0 && faS > 0 && ftS > 0) ? 1 / Math.sqrt((ea / faS) ** 2 + (et / ftS) ** 2) : ftS;
        let fRes;
        if (sOld <= 1) {
          const pa = faS > 0 ? ma[j] * vs / (dt * faS) : 1e30, pt = ftS > 0 ? mt[j] * ws / (r * dt * ftS) : 1e30;
          if (pa * pa + pt * pt <= 1) { v[j] = 0; w[j] = 0; continue; }
          fRes = fSd;
        } else if (weakExp) fRes = muk[j] * fn + fyj + (fSd - muk[j] * fn - fyj) * Math.exp(-sOld);
        else fRes = muk[j] * fn + fyj;
        const [va, vt] = coulomb2d(vs, vtS, ma[j], mt[j] / (r * r), dt * fRes);
        v[j] = va; w[j] = vt / r;
      } else {
        const fn = fn0[j];
        let tStop = mt[j] * ws / dt, cap;
        if (Math.abs(wo) <= wc) cap = (mus[j] * fn + fyj) * r;
        else if (weakExp) cap = ((muk[j] + (mus[j] - muk[j]) * Math.exp(-Math.abs(wo) / wc)) * fn + fyj) * r;
        else cap = (muk[j] * fn + fyj) * r;
        tStop = Math.max(-cap, Math.min(cap, tStop));
        w[j] = ws - dt * tStop / mt[j];
        let fStop = ma[j] * vs / dt;
        const muSa = muk[j] + axc * (mus[j] - muk[j]);
        if (Math.abs(vo) <= vc) cap = muSa * fn + fyj;
        else if (weakExp) cap = (muk[j] + (muSa - muk[j]) * Math.exp(-Math.abs(vo) / vc)) * fn + fyj;
        else cap = muk[j] * fn + fyj;
        fStop = Math.max(-cap, Math.min(cap, fStop));
        v[j] = vs - dt * fStop / ma[j];
      }
    }
    // ---- low-pass states and bookkeeping (5.8)
    for (let j = 0; j < n; j++) { wbar[j] += p.hpT * (w[j] - wbar[j]); vbar[j] += p.hpA * (v[j] - vbar[j]); }
    sc.t += dt;
    sc.bitd += sc.vbit * dt;
    sc.blkpos += sc.vblk * dt;
    sc.wint += Math.abs(sc.wbit) * dt;
    if (bitModel === 2) {
      const phiNew = sc.phib + sc.wbit * dt;
      sweep(sc.phib, phiNew, sc.bitd, surf, nBl);
      sc.phib = phiNew > 1e6 ? phiNew - 2 * Math.PI * Math.floor(phiNew / (2 * Math.PI)) : phiNew;
    }
    if (sc.bitd > sc.holed) sc.holed = sc.bitd;

    if (++s.recPhase >= s.recEvery) {
      s.recPhase = 0;
      rec.push([sc.t, sc.wtd, sc.taum, sc.tau0, sc.wbit, sc.tob, sc.vblk, sc.hkld, sc.vbit, sc.wob, sc.bitd, sc.holed, sc.wsp]);
    }
  }
}
