// Drillstring components and grid (spec 2.3)
import { STEEL_DENSITY, STEEL_E, STEEL_G } from './units.js';

export function component(c) {
  const od = c.od, id = c.id;
  const A = Math.PI / 4 * (od * od - id * id), J = Math.PI / 32 * (od ** 4 - id ** 4);
  const m = c.linear_mass != null ? c.linear_mass : STEEL_DENSITY * A;
  return { name: c.name || '', length: c.length, od, id, A, J, rho: m / A,
           rc: (c.contact_od != null ? c.contact_od : od) / 2, E: STEEL_E, G: STEEL_G };
}

export class Drillstring {
  constructor(components) {
    if (!components.length) throw new Error('At least one component is required');
    this.components = components.map(component);        // bit first
  }
  get bhaLength() { return this.components.slice(0, -1).reduce((s, c) => s + c.length, 0); }

  discretize(bitDepth, dx) {
    const topLen = bitDepth - this.bhaLength;
    if (topLen <= 0) throw new Error('Bit depth is shallower than the BHA length');
    const comps = this.components.slice().reverse();          // surface down
    const g = { dx: [], od: [], id: [], rc: [], A: [], J: [], rho: [], E: [], G: [], comp: [] };
    comps.forEach((c, k) => {
      const L = k === 0 ? topLen : c.length;
      const n = Math.max(1, Math.ceil(L / dx - 1e-9)), h = L / n;
      for (let i = 0; i < n; i++) {
        g.dx.push(h); g.od.push(c.od); g.id.push(c.id); g.rc.push(c.rc); g.A.push(c.A); g.J.push(c.J);
        g.rho.push(c.rho); g.E.push(c.E); g.G.push(c.G); g.comp.push(comps.length - 1 - k);
      }
    });
    const n = g.dx.length, F = (a) => Float64Array.from(a);
    const grid = { n, dx: F(g.dx), od: F(g.od), id: F(g.id), rc: F(g.rc), A: F(g.A), J: F(g.J), rho: F(g.rho),
                   E: F(g.E), G: F(g.G), comp: Int32Array.from(g.comp), xTop: new Float64Array(n) };
    for (let j = 1; j < n; j++) grid.xTop[j] = grid.xTop[j - 1] + grid.dx[j - 1];
    grid.cA = grid.E.map((E, j) => Math.sqrt(E / grid.rho[j]));
    grid.cT = grid.G.map((G, j) => Math.sqrt(G / grid.rho[j]));
    return grid;
  }
}
