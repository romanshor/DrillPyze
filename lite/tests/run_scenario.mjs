// Run an input script through the lite engine; print JSON results (used by tests/test_lite_engine.py).
// usage: node run_scenario.mjs scenario.json   where scenario = {config, script: [{duration, inputs}]}
import { readFileSync } from 'node:fs';
import { build } from '../engine/config.js';

const sc = JSON.parse(readFileSync(process.argv[2], 'utf8'));
const t0 = performance.now();
const sim = build(sc.config);
const out = { info: sim.info(), rock0: sim.rock, latest: [] };
for (const seg of sc.script) out.latest.push(sim.step(seg.duration, seg.inputs || {}));
out.wall = (performance.now() - t0) / 1000;
out.rec = sim.rec;
out.hyd = sim.hydRecords;
out.slow = sim.slowRecords;
process.stdout.write(JSON.stringify(out));
