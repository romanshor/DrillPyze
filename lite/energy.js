// Energy panel: where the power goes (docs/derivation sec. 11, spec 5.9 and 8.4).
//
// Each frame carries the engine's cumulative energies (J) since the start of the run.  The
// panel keeps one snapshot per frame and differences two of them for a mean power over the
// averaging window (1, 10 or 60 s), which it draws as a Sankey diagram (sources, hubs, the
// mud motor, sinks), a power time series and summary tiles.  On a phone the sinks are
// grouped so every label fits; in demo mode the panel is the fourth page.
/* global window, document, getComputedStyle, IntersectionObserver, requestAnimationFrame */
(function () {
  'use strict';
  const SVGNS = 'http://www.w3.org/2000/svg', MAX_HIST = 40000, KWH = 3.6e6, HP = 745.69987, PSI = 6894.757;
  const S = { hist: [], win: 60, view: 'all', opts: null, visible: true, lastDraw: 0, timer: null, shown: null };
  const $ = (id) => document.getElementById(id);
  const css = (v) => getComputedStyle($('energy')).getPropertyValue(v).trim();
  const el = (tag, attrs = {}, parent) => {
    const e = document.createElementNS(SVGNS, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (parent) parent.appendChild(e);
    return e;
  };
  const field = () => S.opts && S.opts.system() === 'field';
  const unit = () => (field() ? 'hp' : 'kW');
  const fac = () => (field() ? 1 / HP : 1 / 1000);
  const fmt = (w) => { const x = w * fac(), a = Math.abs(x); return a >= 100 ? x.toFixed(0) : a >= 10 ? x.toFixed(1) : x.toFixed(2); };
  const fmtE = (j) => {
    const k = j / KWH, a = Math.abs(k);
    if (a >= 1000) return `${(k / 1000).toFixed(a >= 1e4 ? 1 : 2)} MWh`;
    return `${a >= 100 ? k.toFixed(0) : a >= 10 ? k.toFixed(1) : a >= 0.1 ? k.toFixed(2) : k.toFixed(3)} kWh`;
  };
  const pct = (a, b) => (b > 0 ? `${Math.round(100 * a / b)} %` : '–');

  // efficiencies: the case editor's values (live), else the ones the run started with
  function eff(e) {
    const c = S.opts ? S.opts.efficiencies() : {};
    const ok = (v) => v != null && v !== '' && Number(v) > 0 && Number(v) <= 1;
    return { p: ok(c.pump) ? Number(c.pump) : (e.pump_efficiency || 0.9), r: ok(c.regen) ? Number(c.regen) : (e.regen_efficiency || 0.9) };
  }

  // ------------------------------------------------------------------ history and means
  function push(e, latest) {
    if (!e) return;
    const h = S.hist, last = h[h.length - 1];
    if (last && e.time < last.time - 1e-9) h.length = 0;              // a new run
    if (last && Math.abs(e.time - last.time) < 1e-9) h.pop();
    h.push({ ...e, hole_depth: latest ? latest.hole_depth : null, bit_depth: latest ? latest.bit_depth : null,
             spp: latest ? latest.spp : null, q_in: latest ? latest.q_in : null });
    if (h.length > MAX_HIST) h.splice(0, h.length - MAX_HIST);
  }
  function reset() { S.hist.length = 0; S.shown = null; }
  // index of the last snapshot at or before t
  function at(t) {
    const h = S.hist; let lo = 0, hi = h.length - 1;
    if (hi < 0 || t == null) return hi;
    if (h[0].time > t) return 0;
    while (lo < hi) { const m = (lo + hi + 1) >> 1; if (h[m].time <= t + 1e-9) lo = m; else hi = m - 1; }
    return lo;
  }
  // mean power flows (W) over the window ending at snapshot i (the terms of DrillPyze.energy.flows)
  function meansAt(i, w) {
    const h = S.hist; if (i < 1) return null;
    let j = at(h[i].time - w); if (j >= i) j = i - 1;
    const a = h[j], b = h[i], dt = b.time - a.time;
    if (!(dt > 0)) return null;
    const d = (k) => ((b[k] ?? 0) - (a[k] ?? 0)) / dt, n = eff(b);
    const p = {
      interval: dt, time: b.time, hasMotor: !!b.has_motor, eta_p: n.p, eta_r: n.r, bitDiameter: b.bit_diameter,
      td_pos: d('td_pos'), td_neg: d('td_neg'), top: d('top'), bit_t: d('bit_torsion'), motor: d('motor'), bit_a: d('bit_axial'),
      hook_pos: d('hook_pos'), hook_neg: d('hook_neg'), gravity: d('gravity'),
      fric_t: d('friction_torsion'), fric_a: d('friction_axial'), damp_t: d('damping_torsion'), damp_a: d('damping_axial'),
      h_pump: d('h_pump'), h_surface: d('h_surface'), h_string: d('h_string'), h_motor: d('h_motor'),
      h_nozzles: d('h_nozzles'), h_annulus: d('h_annulus'), h_other: d('h_other'),
      rop: b.hole_depth != null && a.hole_depth != null ? Math.max(0, b.hole_depth - a.hole_depth) / dt : 0,
    };
    p.rock = p.bit_t + p.bit_a + p.motor;
    p.pump_drive = p.h_pump / n.p;
    p.braking = p.td_neg + p.hook_neg;
    p.recoverable = p.braking * n.r;
    p.wall = p.fric_t + p.fric_a;
    p.damping = p.damp_t + p.damp_a;
    p.input = p.pump_drive + p.td_pos + p.hook_pos + Math.max(0, p.gravity);
    return p;
  }

  // ------------------------------------------------------------------ flows
  // sources (column 0) -> hubs (1) -> mud motor (2) -> where the power ends up (3)
  const SHORT = { grav: 'Gravity', hoist: 'Drawworks', release: 'Released' };
  function flows(p, grouped) {
    const C = { mud: css('--e-mud'), rot: css('--e-rot'), ax: css('--e-ax'), heat: css('--e-heat'), store: css('--e-store'), regen: css('--e-regen') };
    const nodes = {}, links = [];
    const DEF = {
      pump: ['Mud Pumps', 0, C.mud, 0], td: ['Top drive', 0, C.rot, 1], hoist: ['Drawworks (hoisting)', 0, C.ax, 2],
      grav: ['Gravity (lowering)', 0, C.ax, 3], release: ['Released by the string', 0, C.store, 4],
      mud: ['Mud circulation', 1, C.mud, 0], ploss: ['Pump losses', 1, C.heat, 0.5], rot: ['String rotation', 1, C.rot, 1],
      ax: ['String axial', 1, C.ax, 2], motor: ['Mud motor', 2, C.mud, 0],
      hpipe: ['Frictional Losses inside the string', 3, C.heat, 0], hnoz: ['Bit nozzles (jet)', 3, C.heat, 1],
      hann: ['Frictional Losses in the Annulus', 3, C.heat, 2], hsurf: ['Surface lines', 3, C.heat, 3],
      hres: ['Hydrostatic and compression', 3, C.store, 4], mloss: ['Motor losses', 3, C.heat, 5],
      rock: ['Rock Destruction', 3, C.heat, 7], tfric: ['Wall friction (torque)', 3, C.heat, 8], afric: ['Wall friction (drag)', 3, C.heat, 9],
      damp: ['Viscous damping', 3, C.heat, 10], tdbrake: ['Top-drive braking', 3, C.regen, 11], dwbrake: ['Drawworks braking', 3, C.regen, 12],
      lift: ['Potential energy gained', 3, C.store, 13], stored: ['Stored in the string', 3, C.store, 14],
    };
    const GROUP = { hpipe: 'mudf', hann: 'mudf', hsurf: 'mudf', hres: 'mudf', tfric: 'wall', afric: 'wall', damp: 'wall',
                    tdbrake: 'regen', dwbrake: 'regen', lift: 'stored' };
    const GDEF = { mudf: ['Mud friction', 3, C.heat, 0], wall: ['Wall friction', 3, C.heat, 8], regen: ['Recoverable braking', 3, C.regen, 11] };
    const node = (id) => {
      const gid = grouped && GROUP[id] ? GROUP[id] : id, d = GDEF[gid] || DEF[gid];
      return (nodes[gid] ??= { id: gid, label: d[0], col: d[1], color: d[2], order: d[3], in: 0, out: 0 });
    };
    const flow = (a, b, v, color, min = 20) => {
      if (!Number.isFinite(v) || v < min) return;
      const s = node(a), t = node(b), ex = links.find(l => l.s === s && l.t === t);
      if (ex) ex.v += v; else links.push({ s, t, v, color });
      s.out += v; t.in += v;
    };
    const sflow = (hub, v, color, min) => (v >= 0 ? flow(hub, 'stored', v, color, min) : flow('release', hub, -v, color, min));
    if (S.view === 'all') {
      flow('pump', 'mud', p.h_pump, C.mud); flow('pump', 'ploss', p.pump_drive - p.h_pump, C.mud);
      flow('mud', 'hpipe', p.h_string, C.mud); flow('mud', 'hnoz', p.h_nozzles, C.mud); flow('mud', 'hann', p.h_annulus, C.mud);
      flow('mud', 'hsurf', p.h_surface, C.mud); flow('mud', 'hres', p.h_other, C.mud);
      if (p.hasMotor) { flow('mud', 'motor', p.h_motor, C.mud); flow('motor', 'rock', p.motor, C.mud); flow('motor', 'mloss', p.h_motor - p.motor, C.mud); }
    }
    // rotation: the top-drive motor's positive work in, its braking out; the rest of the balance
    // (wind-up of the string and spin of the top drive, and the scheme's numerical loss) is storage
    const minT = Math.max(20, 0.01 * p.td_pos);
    flow('td', 'rot', p.td_pos, C.rot); flow('rot', 'tdbrake', p.td_neg, C.regen);
    flow('rot', 'rock', p.bit_t, C.rot); flow('rot', 'tfric', p.fric_t, C.rot); flow('rot', 'damp', p.damp_t, C.rot);
    sflow('rot', p.td_pos - p.td_neg - p.bit_t - p.fric_t - p.damp_t, C.store, minT);
    // axial: gravity does work lowering (or the string gains potential energy hoisting)
    const minA = Math.max(20, 0.01 * Math.abs(p.gravity) + 0.01 * p.hook_pos);
    if (p.gravity >= 0) flow('grav', 'ax', p.gravity, C.ax); else flow('ax', 'lift', -p.gravity, C.ax);
    flow('hoist', 'ax', p.hook_pos, C.ax); flow('ax', 'dwbrake', p.hook_neg, C.regen);
    flow('ax', 'rock', p.bit_a, C.ax); flow('ax', 'afric', p.fric_a, C.ax); flow('ax', 'damp', p.damp_a, C.ax);
    sflow('ax', p.gravity + p.hook_pos - p.hook_neg - p.bit_a - p.fric_a - p.damp_a, C.store, minA);
    return { nodes: Object.values(nodes), links };
  }

  // ------------------------------------------------------------------ Sankey (SVG)
  function drawSankey(svg, p, opt) {
    const { width, height, grouped, labelSize = 12.5, compact = false } = opt;
    svg.innerHTML = '';
    svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
    if (compact) { svg.setAttribute('width', '100%'); svg.removeAttribute('height'); }
    else { svg.setAttribute('width', width); svg.setAttribute('height', height); }
    if (!p) return;
    const { nodes, links } = flows(p, grouped);
    if (!nodes.length) {
      const t = el('text', { x: width / 2, y: height / 2, 'text-anchor': 'middle', class: 'sk-value' }, svg);
      t.textContent = 'No power flowing: start the pumps, the rotation or the drawworks.';
      return;
    }
    for (const n of nodes) n.v = Math.max(n.in, n.out);
    const cols = [...new Set(nodes.map(n => n.col))].sort((a, b) => a - b);
    const colIndex = Object.fromEntries(cols.map((c, k) => [c, k]));
    const nodeW = compact ? 9 : 12, pad = compact ? 5 : 8, minSlot = compact ? 29 : 30;
    const FR = compact ? (cols.length === 4 ? [0, 0.3, 0.62, 1] : [0, 0.3, 1]) : null;
    const xOf = (c) => {
      const k = colIndex[c], nC = cols.length;
      if (k === nC - 1) return width - nodeW;
      return FR ? FR[k] * (width - nodeW) : (k / Math.max(1, nC - 1)) * (width - nodeW) * 0.7;
    };
    // scale: the tallest column fits, with every node at least one label high
    const byCol = cols.map(c => nodes.filter(n => n.col === c).sort((a, b) => a.order - b.order));
    let k = Infinity;
    for (const col of byCol) k = Math.min(k, (height - pad * (col.length - 1)) / col.reduce((s, n) => s + n.v, 0));
    for (let it = 0; it < 40; it++) {
      if (byCol.every(col => col.reduce((s, n) => s + Math.max(n.v * k, minSlot), 0) + pad * (col.length - 1) <= height)) break;
      k *= 0.93;
    }
    for (const col of byCol) {
      let y = 0;
      for (const n of col) {
        n.h = Math.max(1.5, n.v * k); n.slot = Math.max(n.h, minSlot); n.x = xOf(n.col);
        n.y = y + (n.slot - n.h) / 2; y += n.slot + pad;
      }
    }
    // the mud motor branches off the bottom of the mud circulation, level with it
    for (const n of nodes.filter(n => n.id === 'motor')) {
      const hub = links.find(l => l.t === n)?.s;
      if (hub) n.y = hub.y + hub.h - n.h;
    }
    const outs = new Map(), ins = new Map();
    for (const l of links) {
      if (!outs.has(l.s)) outs.set(l.s, []); outs.get(l.s).push(l);
      if (!ins.has(l.t)) ins.set(l.t, []); ins.get(l.t).push(l);
    }
    for (const [n, ls] of outs) {
      ls.sort((a, b) => (a.t.id === 'motor') - (b.t.id === 'motor') || a.t.y - b.t.y);
      let y = n.y; for (const l of ls) { l.w = l.v * k; l.y0 = y; y += l.w; }
    }
    for (const [n, ls] of ins) { ls.sort((a, b) => a.s.y - b.s.y); let y = n.y; for (const l of ls) { l.y1 = y; y += l.v * k; } }
    const gL = el('g', {}, svg), gN = el('g', {}, svg), gT = el('g', {}, svg), panel = css('--panel'), u = unit();
    for (const l of links) {
      const x0 = l.s.x + nodeW, x1 = l.t.x, xm = (x0 + x1) / 2, w = Math.max(l.w, 0.8);
      const path = el('path', { d: `M${x0},${l.y0} C${xm},${l.y0} ${xm},${l.y1} ${x1},${l.y1} L${x1},${l.y1 + w} C${xm},${l.y1 + w} ${xm},${l.y0 + w} ${x0},${l.y0 + w} Z`,
                                fill: l.color, 'fill-opacity': 0.32 }, gL);
      el('title', {}, path).textContent = `${l.s.label} → ${l.t.label}: ${fmt(l.v)} ${u}`;
    }
    for (const n of nodes) {
      const r = el('rect', { x: n.x, y: n.y, width: nodeW, height: n.h, fill: n.color, rx: 1.5 }, gN);
      el('title', {}, r).textContent = `${n.label}: ${fmt(n.v)} ${u}`;
      const last = colIndex[n.col] === cols.length - 1;
      // sinks are labelled on the left of their bar; on a phone the motor's label sits on its
      // left too (clear of the sink labels) and a tall hub is named at its top
      const left = last || (compact && n.id === 'motor');
      let cy = n.y + n.h / 2;
      if (compact && !last && n.col > 0 && n.h > 44) cy = n.y + 14;
      const tx = left ? n.x - 6 : n.x + nodeW + 6, anchor = left ? 'end' : 'start';
      const t = el('text', { x: tx, y: cy - (compact ? 1 : 2), 'text-anchor': anchor, class: 'sk-label', 'font-size': labelSize }, gT);
      t.textContent = compact && SHORT[n.id] ? SHORT[n.id] : n.label;
      const v = el('text', { x: tx, y: cy + (compact ? 10 : 12), 'text-anchor': anchor, class: 'sk-value', 'font-size': labelSize - 1.5 }, gT);
      v.textContent = `${fmt(n.v)} ${u}`;
      if (!last) for (const e of [t, v]) { e.setAttribute('paint-order', 'stroke'); e.setAttribute('stroke', panel); e.setAttribute('stroke-width', 3); e.setAttribute('stroke-linejoin', 'round'); }
    }
  }

  // ------------------------------------------------------------------ tiles and state line
  function drawTiles(p, i, phone) {
    const u = unit(), h = S.hist[i], e0 = S.hist[0];
    let tiles;
    if (!p) tiles = [];
    else if (phone) {
      tiles = [['heat', 'Rock Destruction', `${fmt(p.rock)} ${u}`], ['mud', 'Mud Pumps', `${fmt(p.pump_drive)} ${u}`],
               ['regen', 'Recoverable', `${fmt(p.recoverable)} ${u}`]];
    } else {
      const area = Math.PI / 4 * (p.bitDiameter || 0) ** 2, mse = p.rop * 3600 > 1 && area > 0 ? p.rock / (area * p.rop) : null;
      const n = eff(h), cum = (k) => (h[k] ?? 0) - (e0[k] ?? 0);
      const inE = cum('h_pump') / n.p + cum('td_pos') + cum('hook_pos') + Math.max(0, cum('gravity'));
      const brE = cum('td_neg') + cum('hook_neg');
      tiles = [
        ['', 'Total input', `${fmt(p.input)} ${u}`, `mud pumps ${pct(p.pump_drive, p.input)}`],
        ['heat', 'Rock Destruction', `${fmt(p.rock)} ${u}`,
         mse ? `MSE ${field() ? `${(mse / PSI / 1000).toFixed(1)} ksi` : `${(mse / 1e6).toFixed(0)} MPa`}` : 'no rock work off bottom'],
        p.hasMotor ? ['mud', 'Mud motor', `${fmt(p.h_motor)} ${u}`, p.h_motor > 1 ? `from the mud, ${pct(p.motor, p.h_motor)} reaches the bit` : 'no flow']
                   : ['rot', 'Top drive', `${fmt(p.td_pos)} ${u}`, `${pct(p.bit_t, p.td_pos)} reaches the bit`],
        ['ax', 'Wall friction', `${fmt(p.wall)} ${u}`, `damping ${fmt(p.damping)} ${u}`],
        ['mud', 'Mud Pumps', `${fmt(p.pump_drive)} ${u}`, `SPP·Q ${fmt(p.h_pump)} ${u} at ${Math.round(p.eta_p * 100)} %`],
        ['regen', 'Recoverable braking', `${fmt(p.recoverable)} ${u}`, `of ${fmt(p.braking)} ${u} braked, at ${Math.round(p.eta_r * 100)} %`],
        ['', 'Since the start', fmtE(inE), `in; ${fmtE(brE * n.r)} recoverable`],
      ];
    }
    $('e-tiles').innerHTML = tiles.map(([c, k, v, s]) => `<div class="e-tile ${c}"><div class="k">${k}</div><div class="v">${v}</div>${s ? `<div class="s">${s}</div>` : ''}</div>`).join('');
    if (!p) { $('e-state').textContent = 'Energy flows appear after the first frame.'; return; }
    const dep = field() ? `${(h.bit_depth / 0.3048).toFixed(0)} ft` : `${h.bit_depth.toFixed(1)} m`;
    const extra = h.spp != null && !phone ? ` · SPP ${field() ? `${(h.spp / PSI).toFixed(0)} psi` : `${(h.spp / 1e5).toFixed(1)} bar`}` : '';
    $('e-state').textContent = `Mean over ${p.interval >= 10 ? p.interval.toFixed(0) : p.interval.toFixed(1)} s ending at ${p.time.toFixed(1)} s · bit at ${dep}${extra}${p.hasMotor ? ' · mud motor' : ''}`;
  }

  // ------------------------------------------------------------------ power over time (SVG)
  function series(p) {
    const s = [['Mud Pumps', '--e-mud', p.pump_drive, true], ['Top drive', '--e-rot', p.td_pos], ['Rock Destruction', '--e-heat', p.rock]];
    if (p.hasMotor) s.push(['Into the motor', '--e-motor', p.h_motor, true]);
    s.push(['Wall friction', '--e-ax', p.wall], ['Braking', '--e-regen', p.braking]);
    return s.filter(x => S.view === 'all' || !x[3]);
  }
  const niceStep = (max) => { const r = max / 5, q = 10 ** Math.floor(Math.log10(r)); return [1, 2, 5, 10].map(x => x * q).find(x => x >= r); };
  function drawTimeline(i, phone) {
    const svg = $('e-timeline'), h = S.hist, W = Math.max(phone ? 300 : 560, svg.parentElement.clientWidth), H = phone ? 220 : 250;
    const m = { l: 48, r: 10, t: 10, b: 34 };
    svg.innerHTML = ''; svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.setAttribute('width', W); svg.setAttribute('height', H);
    if (i < 1) { $('e-legend-t').innerHTML = ''; return; }
    const tEnd = h[i].time, span = S.opts ? S.opts.logWindow() : 0, tStart = span > 0 ? Math.max(h[0].time, tEnd - span) : h[0].time;
    let j0 = at(tStart); if (h[j0].time < tStart - 1e-9 && j0 < i) j0++;
    j0 = Math.max(1, j0);
    const n = i - j0 + 1, stride = Math.max(1, Math.ceil(n / 300)), idx = [];
    for (let k = i; k >= j0; k -= stride) idx.unshift(k);
    const pts = idx.map(k => meansAt(k, S.win)).filter(Boolean);
    if (!pts.length) return;
    const names = series(pts[pts.length - 1]).map(x => x[0]);
    const ser = names.map((name) => ({ name, color: css(series(pts[pts.length - 1]).find(x => x[0] === name)[1]),
                                       y: pts.map(q => (series(q).find(x => x[0] === name) || [0, 0, 0])[2] * fac()) }));
    const ymax = Math.max(1, ...ser.flatMap(s => s.y)) * 1.06, t0 = Math.min(tStart, pts[0].time), t1 = Math.max(tEnd, t0 + 1);
    const X = (t) => m.l + ((t - t0) / (t1 - t0)) * (W - m.l - m.r), Y = (v) => H - m.b - (Math.max(0, v) / ymax) * (H - m.t - m.b);
    const ink = css('--steel'), grid = css('--line');
    const step = niceStep(ymax);
    for (let v = 0; v <= ymax; v += step) {
      el('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), stroke: grid }, svg);
      el('text', { x: m.l - 6, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 11, fill: ink }, svg).textContent = v >= 10 ? v.toFixed(0) : +v.toFixed(2);
    }
    const ts = niceStep(t1 - t0);
    for (let t = Math.ceil(t0 / ts) * ts; t <= t1 + 1e-9; t += ts) {
      el('text', { x: X(t), y: H - m.b + 15, 'text-anchor': 'middle', 'font-size': 11, fill: ink }, svg).textContent = +t.toFixed(1);
    }
    el('text', { x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', 'font-size': 11.5, fill: ink }, svg).textContent = 'Time (s)';
    el('text', { x: 12, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', 'font-size': 11.5, fill: ink,
                 transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})` }, svg).textContent = `Power (${unit()})`;
    for (const s of ser) {
      el('path', { d: s.y.map((v, k) => `${k ? 'L' : 'M'}${X(pts[k].time).toFixed(1)},${Y(v).toFixed(1)}`).join(''),
                   fill: 'none', stroke: s.color, 'stroke-width': 2, 'stroke-linejoin': 'round' }, svg);
    }
    $('e-legend-t').innerHTML = ser.map(s => `<span><i style="background:${s.color}"></i>${s.name}</span>`).join('');
  }

  // ------------------------------------------------------------------ draw
  function draw(tShown = null, force = false) {
    S.shown = tShown;
    if (!$('energy') || (!S.visible && !force)) return;
    const now = performance.now();
    // at most ~5 redraws a second while playing; the last frame is always drawn
    if (!force && now - S.lastDraw < 200) {
      clearTimeout(S.timer); S.timer = setTimeout(() => draw(S.shown, true), 200 - (now - S.lastDraw));
      return;
    }
    clearTimeout(S.timer); S.lastDraw = now;
    const phone = S.opts ? S.opts.phone() : false, i = at(tShown == null ? Infinity : tShown);
    const p = i >= 1 ? meansAt(i, S.win) : null;
    drawTiles(p, i, phone);
    const svg = $('e-sankey');
    if (phone) drawSankey(svg, p, { width: 320, height: 600, grouped: true, labelSize: 10.5, compact: true });
    else drawSankey(svg, p, { width: Math.max(680, svg.parentElement.clientWidth), height: S.view === 'all' ? 540 : 420, grouped: false });
    $('e-sk-title').textContent = p ? `Power flow, ${S.view === 'mech' ? 'drillstring only, ' : ''}${unit()}` : 'Power flow';
    drawTimeline(i, phone);
  }

  function radio(id, key, attr, parse) {
    for (const b of document.querySelectorAll(`#${id} button`)) {
      b.addEventListener('click', () => {
        S[key] = parse(b.dataset[attr]);
        for (const o of document.querySelectorAll(`#${id} button`)) o.setAttribute('aria-checked', String(o === b));
        draw(S.shown, true);
      });
    }
  }
  function init(opts) {
    S.opts = opts;
    radio('e-win', 'win', 'w', Number);
    radio('e-view', 'view', 'v', String);
    if (typeof IntersectionObserver !== 'undefined') {
      new IntersectionObserver((es) => {
        const v = es.some(e => e.isIntersecting);
        if (v && !S.visible) { S.visible = true; requestAnimationFrame(() => draw(S.shown, true)); } else S.visible = v;
      }).observe($('energy'));
    }
  }
  window.EnergyPanel = { init, push, reset, draw, history: () => S.hist, means: (t, w) => meansAt(at(t ?? Infinity), w ?? S.win) };
}());
