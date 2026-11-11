"""The static DrillPyze Lite site assembles from the front end and the browser engine."""
import importlib.util
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('build_site', ROOT / 'lite' / 'build_site.py')
bs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bs)


def test_site_builds(tmp_path):
    out = bs.build(tmp_path / 'site')
    html = (out / 'index.html').read_text()
    for f in ('app.js', 'deck.js', 'app.css', 'lite-api.js', 'engine/worker.js', 'engine/kernel.js',
              'engine/hydraulics.js', 'engine/session.js', '.nojekyll'):
        assert (out / f).exists(), f
    assert '/static/' not in html
    assert (html.index('src="lite-api.js"') < html.index('src="deck.js"') < html.index('src="app.js"'))
    assert 'window.DRILLPYZE_LITE' in (out / 'app.js').read_text()
