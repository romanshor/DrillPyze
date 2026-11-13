# Contributing to DrillPyze

Contributions are welcome. Bug reports, feature requests and pull requests can
be opened on GitHub. Please use the issue templates for bugs and features.

## Development set-up

```bash
git clone https://github.com/romanshor/DrillPyze.git
cd DrillPyze
pip install -e ".[test,ui,docs]"
pytest                       # the full test suite (about a minute)
python -m DrillPyze.ui       # the interactive simulator
./docs_generator.sh          # the documentation, into docs/build
```

## Guidelines

- **Units:** keep everything SI inside the model. Convert only at input and
  output (`DrillPyze.units`).
- **Tests:** add or update tests with every change to the physics. The
  verification cases in `tests/` compare against analytical solutions and
  published results; keep them passing.
- **Docstrings:** document new functions and classes (numpy style). The API
  reference is built from them, and the documentation build treats warnings
  as errors.
- **Maths in Markdown:** use GitHub's protected forms (```` ```math ```` blocks
  and `` $`...`$ ``). `tests/test_docs.py` checks this.
- **Changelog:** add user-visible changes to `CHANGELOG.md`.

## Physics specification and the browser version

DrillPyze Lite (the browser version, `lite/`) is built from the
implementation-level specification in `docs/spec/physics-spec.md`. Both live on
the `browser-lite` branch until it is merged into `main`. A change to the physics or numerics of the Python model must:

1. update the specification in the same pull request and raise its revision;
2. refresh the manifest with `python docs/spec/update_manifest.py`;
3. bring the browser engine up to the new revision, so that
   `tests/test_lite_engine.py` still passes.

`tests/test_spec_sync.py` enforces steps 1–3.

## Releases

The version lives in `pyproject.toml`, `DrillPyze/__init__.py` (with the
release date) and, for DrillPyze Lite, `lite/engine/version.js`. Tests check that they agree.
Update all three, add the `CHANGELOG.md` entry, and tag the release
`vX.Y.Z`.
