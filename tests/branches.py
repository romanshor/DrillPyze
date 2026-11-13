"""
Branch-aware strictness for the browser-version checks.

Development happens on two branches and releases on ``main``:

* ``coupled-model`` -- Python physics; the browser engine may lag behind here,
* ``browser-lite``  -- the browser engine is brought up to date with the physics,
* ``main``          -- releases; everything must agree.

On a branch where the browser engine may lag, a Lite conformance failure is
reported as an expected failure (``xfail``, shown in the summary) instead of an
error.  Pull requests into ``browser-lite`` or ``main`` are always strict.  Everywhere else it fails the run.  Physics regressions and an unreviewed
physics change (``tests/test_spec_sync.py::test_physics_sources_unchanged``) fail
on every branch.

Override with the environment variable ``DRILLPYZE_LITE_LAG``: ``allow`` or
``strict``.
"""
import os
import subprocess

import pytest

LAG_ALLOWED_BRANCHES = {'coupled-model'}


def current_branch():
    """Branch under test: CI variables first (pull requests use the source branch), then git."""
    for var in ('DRILLPYZE_BRANCH', 'GITHUB_HEAD_REF', 'GITHUB_REF_NAME'):
        v = os.environ.get(var)
        if v:
            return v
    try:
        out = subprocess.run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], capture_output=True, text=True,
                             cwd=os.path.dirname(__file__), timeout=10)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def lite_lag_allowed():
    """True where the browser engine may lag the Python physics (see module docstring)."""
    mode = os.environ.get('DRILLPYZE_LITE_LAG', '').lower()
    if mode in ('allow', 'strict'):
        return mode == 'allow'
    target = os.environ.get('GITHUB_BASE_REF')        # pull requests: strict unless merging into a lag branch
    if target and target not in LAG_ALLOWED_BRANCHES:
        return False
    return current_branch() in LAG_ALLOWED_BRANCHES


def lite_check(ok, message):
    """Assert a browser-version check: a failure is an expected failure where Lite may lag."""
    if ok:
        return
    if lite_lag_allowed():
        pytest.xfail(f'DrillPyze Lite lags the Python model on branch {current_branch()!r} '
                     f'(update it on browser-lite): {message}')
    pytest.fail(message)
