"""The static DrillPyze Lite site assembles from the front end and the browser engine."""
import importlib.util
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('build_site', ROOT / 'lite' / 'build_site.py')
bs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bs)


import pytest


def test_site_builds(tmp_path):
    pytest.importorskip('plotly')
    out = bs.build(tmp_path / 'site')
    html = (out / 'index.html').read_text()
    for f in ('app.js', 'deck.js', 'app.css', 'lite-api.js', 'plotly.min.js', 'engine/worker.js', 'engine/kernel.js',
              'engine/hydraulics.js', 'engine/session.js', '.nojekyll'):
        assert (out / f).exists(), f
    assert '/static/' not in html
    assert '<script src="plotly.min.js"></script>' in html          # bundled, no CDN
    assert (html.index('src="lite-api.js"') < html.index('src="deck.js"') < html.index('src="app.js"'))
    assert 'window.DRILLPYZE_LITE' in (out / 'app.js').read_text()


def test_lite_version_matches_package():
    """DrillPyze Lite carries the same release as the Python package."""
    import re
    import DrillPyze
    js = (ROOT / 'lite' / 'engine' / 'version.js').read_text()
    get = lambda k: re.search(rf"{k} = '([^']+)'", js).group(1)
    assert get('VERSION') == DrillPyze.__version__
    assert get('RELEASE') == DrillPyze.__release_name__
    assert get('RELEASE_DATE') == DrillPyze.__release_date__
