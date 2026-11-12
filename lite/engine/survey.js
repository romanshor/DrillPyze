// Wellbore survey with minimum-curvature interpolation (spec 2.1)

function tangent(inc, azi) {
  return [Math.sin(inc) * Math.cos(azi), Math.sin(inc) * Math.sin(azi), Math.cos(inc)];
}
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const crossNorm = (a, b) => Math.hypot(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);

// numpy.interp: linear, clamped at the ends; xs increasing
export function interp(x, xs, ys) {
  const n = xs.length;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[n - 1]) return ys[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
  const f = (x - xs[lo]) / (xs[hi] - xs[lo]);
  return ys[lo] + f * (ys[hi] - ys[lo]);
}

export class Survey {
  constructor(md, inc, azi, degrees = true) {
    let rows = md.map((m, i) => [Number(m), Number(inc[i]) * (degrees ? Math.PI / 180 : 1), Number(azi[i]) * (degrees ? Math.PI / 180 : 1)]);
    rows.sort((a, b) => a[0] - b[0]);
    rows = rows.filter((r, i) => i === 0 || r[0] - rows[i - 1][0] > 1e-9);
    if (rows[0][0] > 0) rows.unshift([0, 0, rows[0][2]]);
    if (rows.length < 2) throw new Error('A survey needs at least two stations');
    this.md = rows.map(r => r[0]); this.inc = rows.map(r => r[1]); this.azi = rows.map(r => r[2]);
    this.t = rows.map(r => tangent(r[1], r[2]));
    this.dogleg = [];
    for (let k = 0; k < this.md.length - 1; k++) {
      const d = Math.min(1, Math.max(-1, dot(this.t[k], this.t[k + 1])));
      this.dogleg.push(Math.atan2(crossNorm(this.t[k], this.t[k + 1]), d));
    }
    this._posMd = null;
  }
  get td() { return this.md[this.md.length - 1]; }

  tangent(s) {
    const md = this.md, N = md.length;
    // last station with md <= s, clipped to [0, N-2]
    let k = 0, lo = 0, hi = N;
    while (lo < hi) { const m = (lo + hi) >> 1; if (md[m] <= s) lo = m + 1; else hi = m; }
    k = Math.min(Math.max(lo - 1, 0), N - 2);
    const L = md[k + 1] - md[k];
    const f = Math.min(1, Math.max(0, (s - md[k]) / L));
    const t1 = this.t[k], t2 = this.t[k + 1], b = this.dogleg[k];
    let out;
    if (b < 1e-9) out = [t1[0], t1[1], t1[2]];
    else {
      const sb = Math.sin(b), c1 = Math.sin((1 - f) * b) / sb, c2 = Math.sin(f * b) / sb;
      out = [c1 * t1[0] + c2 * t2[0], c1 * t1[1] + c2 * t2[1], c1 * t1[2] + c2 * t2[2]];
    }
    const nrm = Math.hypot(out[0], out[1], out[2]);
    return [out[0] / nrm, out[1] / nrm, out[2] / nrm];
  }

  _table(need) {
    if (this._posMd && this._posMd[this._posMd.length - 1] >= need) return;
    const top = need * 1.25 + 100.0, n = 6000;
    const set = new Set(this.md);
    for (let i = 0; i < n; i++) set.add(top * i / (n - 1));
    const grid = Array.from(set).sort((a, b) => a - b);
    const pos = [[0, 0, 0]];
    let prev = this.tangent(grid[0]);
    for (let i = 1; i < grid.length; i++) {
      const t = this.tangent(grid[i]), d = grid[i] - grid[i - 1], p = pos[i - 1];
      pos.push([p[0] + 0.5 * (t[0] + prev[0]) * d, p[1] + 0.5 * (t[1] + prev[1]) * d, p[2] + 0.5 * (t[2] + prev[2]) * d]);
      prev = t;
    }
    this._posMd = grid;
    this._pos = [pos.map(p => p[0]), pos.map(p => p[1]), pos.map(p => p[2])];
  }
  // (north, east, tvd) arrays at measured depths
  positions(mds) {
    const need = Math.max(this.td, ...mds);
    this._table(need);
    return [0, 1, 2].map(c => mds.map(m => interp(m, this._posMd, this._pos[c])));
  }
  tvdAt(md) { return this.positions([md])[2][0]; }
}
