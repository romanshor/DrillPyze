"""
Assemble DrillPyze Lite as a static site (for GitHub Pages or any web server).

The site is the regular front end (src/DrillPyze/ui/static) with its API answered
by the browser engine (lite/engine, built from docs/spec/physics-spec.md) in a
Web Worker instead of the Python server.

    python lite/build_site.py            # writes _site/
    python lite/build_site.py --out dist
    python -m http.server -d _site 8000  # try it locally at http://localhost:8000
"""
import argparse
import pathlib
import re
import shutil

ROOT = pathlib.Path(__file__).resolve().parents[1]
STATIC = ROOT / 'src' / 'DrillPyze' / 'ui' / 'static'
LITE = ROOT / 'lite'
PLOTLY_CDN = 'https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.35.2/plotly.min.js'


def index_html():
    html = (STATIC / 'index.html').read_text()
    reps = [
        ('<title>DrillPyze simulator</title>', '<title>DrillPyze Lite</title>'),
        ('<script src="/static/plotly.min.js"></script>', f'<script src="{PLOTLY_CDN}"></script>'),
        ('<h1>DrillPyze</h1>', '<h1>DrillPyze <span class="lite-badge" title="Runs entirely in your browser">Lite</span></h1>'),
        ('<script src="/static/deck.js"></script>', '<script src="lite-api.js"></script>\n<script src="deck.js"></script>'),
        ('<script src="/static/app.js"></script>', '<script src="app.js"></script>'),
        ('<link rel="stylesheet" href="/static/app.css">',
         '<link rel="stylesheet" href="app.css">\n<style>.lite-badge { font: 600 12px/1 var(--font); vertical-align: middle; '
         'padding: 3px 7px; border-radius: 10px; background: #e8f1f9; color: #164e7e; border: 1px solid #c6dbee; }</style>'),
    ]
    for a, b in reps:
        if a not in html:
            raise SystemExit(f'index.html changed: cannot find {a!r}')
        html = html.replace(a, b)
    if re.search(r'["\']/static/', html):
        raise SystemExit('index.html still references /static/ paths')
    return html


def build(out):
    out = pathlib.Path(out)
    if out.exists():
        shutil.rmtree(out)
    (out / 'engine').mkdir(parents=True)
    (out / 'index.html').write_text(index_html())
    for f in ('app.js', 'deck.js', 'app.css'):
        shutil.copy(STATIC / f, out / f)
    shutil.copy(LITE / 'lite-api.js', out / 'lite-api.js')
    for f in sorted((LITE / 'engine').glob('*.js')):
        shutil.copy(f, out / 'engine' / f.name)
    (out / '.nojekyll').write_text('')          # serve files as they are (no Jekyll processing)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=str(ROOT / '_site'))
    out = build(ap.parse_args().out)
    print(f'DrillPyze Lite written to {out}')
