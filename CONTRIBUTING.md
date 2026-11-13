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
implementation-level specification in `docs/spec/physics-spec.md`. A change to the physics or numerics of the Python model must:

1. update the specification in the same pull request and raise its revision;
2. refresh the manifest with `python docs/spec/update_manifest.py`;
3. bring the browser engine up to the new revision, so that
   `tests/test_lite_engine.py` still passes.

`tests/test_spec_sync.py` enforces steps 1–3.

## Branches and releases

- **`coupled-model`:** development of the Python model, the simulator UI and
  the documentation. Work here, or in feature branches merged into it. Each
  push republishes the development docs at `/coupled-model/` on the
  documentation site.
- **`browser-lite`:** merge `coupled-model` into it, then bring DrillPyze Lite
  up to date. Update the specification, run the manifest script and update the
  browser engine, as described above.
- **`main`:** releases only. Release by pull request from `browser-lite`
  into `main`. Each push to `main` publishes the release docs and DrillPyze Lite.
- **After a release:** merge `main` back into `coupled-model` and
  `browser-lite`.

**Branch-aware tests.** The Lite checks are strict on `browser-lite`, on
`main` and in pull requests into them. These are `tests/test_lite_engine.py`
and the spec-revision check in `tests/test_spec_sync.py`.

On `coupled-model` the browser engine may lag the physics, so a Lite mismatch
is reported as an expected failure (`xfail`, listed in the test summary)
instead of an error. Physics regressions still fail on every branch, and so
does a physics change without a spec review.

The branch is taken from the CI environment or from `git`. To override it,
set `DRILLPYZE_LITE_LAG=strict` or `allow`. For example, run
`DRILLPYZE_LITE_LAG=strict pytest` on `coupled-model` to see how far Lite has
fallen behind. The rules are in `tests/branches.py`.

**Development versions.** After a release, set the development branches to
the next pre-release with a `.dev0` suffix (for example `1.0.0a2.dev0`). The
UI and the docs then show that the build is not a release.

## Releases

The version lives in `pyproject.toml`, `DrillPyze/__init__.py` (with the
release date) and, for DrillPyze Lite, `lite/engine/version.js`. Tests check that they agree.
Update all three, add the `CHANGELOG.md` entry, and tag the release
`vX.Y.Z`.
