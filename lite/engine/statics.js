// Steady soft-string torque and drag from the bit up (spec 11)
export function torqueDrag(geo, mu, v, omega, wob = 0, tob = 0, fnFixed = null) {
  const n = geo.n, T = new Float64Array(n + 1), M = new Float64Array(n + 1), FN = new Float64Array(n);
  T[n] = -wob; M[n] = tob;
  for (let j = n - 1; j >= 0; j--) {
    let fn;
    if (fnFixed) fn = fnFixed[j];
    else {
      const fx = T[j + 1] * geo.kx[j] + geo.wb[j] * geo.zx[j], fy = T[j + 1] * geo.ky[j] + geo.wb[j] * geo.zy[j],
            fz = T[j + 1] * geo.kz[j] + geo.wb[j] * geo.zz[j];
      fn = Math.sqrt(fx * fx + fy * fy + fz * fz);
    }
    const vt = omega * geo.rc[j], sp = Math.sqrt(v * v + vt * vt);
    const ca = sp > 0 ? v / sp : 0, ct = sp > 0 ? vt / sp : 0, f = mu[j] * fn;
    T[j] = T[j + 1] + (geo.wb[j] * geo.tz[j] - f * ca) * geo.dx[j];
    M[j] = M[j + 1] + f * ct * geo.rc[j] * geo.dx[j];
    FN[j] = fn;
  }
  return { tension: T, torque: M, normal: FN };
}
