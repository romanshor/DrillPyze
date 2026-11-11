# DrillPyze Lite

DrillPyze running entirely in a web browser, with no server or Python install.
It is the regular DrillPyze front end; the simulation runs in a Web Worker using
a JavaScript engine (`engine/`) written from the physics specification
(`docs/spec/physics-spec.md`).

**Physics:** the full Python model:
- the coupled axial–torsional string with soft-string friction and damping;
- the bit models, top-drive controllers and autodriller;
- steady and transient hydraulics;
- the wellbore temperature model;
- rock strength from the formation.

`tests/test_lite_engine.py` checks that the engine reproduces the Python model
to round-off on Well A, with and without the temperature model.

**Speed:** about 20–30× real time for Well A on a laptop browser.

## Files

| Path | Purpose |
|---|---|
| `engine/*.js` | the simulation engine (ES modules, no dependencies) |
| `engine/spec.js` | the spec revision the engine implements |
| `engine/worker.js` | the Web Worker answering the UI's API calls |
| `lite-api.js` | routes the front end's API calls to the worker |
| `build_site.py` | assembles the static site |
| `tests/run_scenario.mjs` | runs an input script under Node.js (used by the conformance test) |

## Build and try locally

```bash
python lite/build_site.py                 # writes _site/
python -m http.server -d _site 8000       # open http://localhost:8000
```

A web server is needed: browsers do not run module workers from `file://`.

## Keeping it in step with the Python model

When the Python physics changes, `tests/test_spec_sync.py` fails until:

1. the specification is updated and its revision raised;
2. `python docs/spec/update_manifest.py` is run.

Then:

3. update the engine from the changed spec sections;
4. set `SPEC_REVISION` in `engine/spec.js`;
5. make `tests/test_lite_engine.py` pass.

## Publishing on GitHub Pages

`.github/workflows/pages-lite.yml` runs the conformance tests, builds the site
and publishes it to the `lite/` folder of the `gh-pages` branch, next to the
documentation. See the main README for the one-time set-up.
