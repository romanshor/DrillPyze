// Rock strength at the bit (spec 9): P1 confinement and T1 weakening; T2 needs a mud temperature
export const FORMATION_DEFAULTS = {
  layers: [{ top_md: 0, ucs: 110e6, friction_angle: 36, porosity: 0.005, name: 'Granitoid' }],
  pore_pressure_gradient: 9.79e3, overburden_gradient: 25.0e3, skempton_b: 1.0, porosity_low: 0.05, porosity_high: 0.20,
  surface_temperature: 40.0, geothermal_gradient: 0.0691, rock_density: 2800.0, rock_heat_capacity: 930.0,
  rock_conductivity: 2.31, thermal_expansion: 6.0e-6, youngs_modulus: 55e9, poisson_ratio: 0.26,
  weakening_threshold: 193.0, weakening_rate: 1.41e-3, thermal_stress_efficiency: 0.5, strength_floor: 0.1,
};

export function rockStrength(f, md, tvd, pbh, tMud = null, rop = 0, omega = 0, blades = 5) {
  const lay = f.layers.slice().sort((a, b) => a.top_md - b.top_md).reduce((o, l) => (md >= l.top_md ? l : o), f.layers.slice().sort((a, b) => a.top_md - b.top_md)[0]);
  const tRock = f.surface_temperature + f.geothermal_gradient * tvd;
  const ucsT = lay.ucs * Math.max(1 - f.weakening_rate * Math.max(tRock - f.weakening_threshold, 0), 0.05);
  const pp = f.pore_pressure_gradient * tvd, sv = f.overburden_gradient * tvd;
  const dpHp = pbh - pp, dpLp = pbh - (pp - f.skempton_b * (sv - pbh) / 3);
  const phi = lay.porosity, lo = f.porosity_low, hi = f.porosity_high;
  const w = phi <= lo ? 0 : phi >= hi ? 1 : (phi - lo) / (hi - lo);
  const s = Math.sin(lay.friction_angle * Math.PI / 180), nPhi = (1 + s) / (1 - s);
  const ccs = ucsT + ((1 - w) * dpLp + w * dpHp) * nPhi;
  let sigmaTh = 0, ratio = 0;
  if (tMud != null) {
    sigmaTh = f.thermal_expansion * f.youngs_modulus * (tRock - tMud) / (1 - f.poisson_ratio);
    ratio = 1;
    if (omega > 1e-6 && rop > 1e-9) {
      const tp = 2 * Math.PI / (Math.max(Math.trunc(blades), 1) * omega);
      const kap = f.rock_conductivity / (f.rock_density * f.rock_heat_capacity);
      ratio = Math.min(1, Math.sqrt(Math.PI * kap * tp) / (rop * tp));
    }
  }
  const ccsEff = Math.max(ccs - f.thermal_stress_efficiency * sigmaTh * ratio, f.strength_floor * ucsT);
  return { ucs0: lay.ucs, t_rock: tRock, ucs_t: ucsT, n_phi: nPhi, dp_drained: dpHp, dp_undrained: dpLp,
           w_drained: w, ccs, sigma_th: sigmaTh, skin_ratio: ratio, ccs_eff: ccsEff };
}
