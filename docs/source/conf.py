# Sphinx configuration for the DrillPyze documentation (coupled model).
#
# The narrative pages are the repository's Markdown files (README, derivations and,
# when present, the physics specification).  They use GitHub's maths syntax
# (```math fences and $`...`$), so at build time they are copied into
# docs/source/_generated with the maths converted to MyST and image paths rebased.
# The API reference comes from the docstrings (autodoc + napoleon).
import os
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'src'))

# -- Project information -----------------------------------------------------
project = 'DrillPyze'
author = 'Roman J. Shor'
copyright = '2023-2026, Roman J. Shor'
try:
    import DrillPyze
    release = DrillPyze.__version__
except Exception:          # documentation still builds if the package cannot import
    release = '0.1'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'myst_parser',
]
myst_enable_extensions = ['dollarmath', 'amsmath']
myst_heading_anchors = 3
source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}
exclude_patterns = ['_build']
templates_path = []
autodoc_member_order = 'bysource'
autodoc_default_options = {'members': True, 'undoc-members': True, 'show-inheritance': True}
autodoc_typehints = 'description'
# numba is optional at documentation time: mock it so the modules import anywhere
autodoc_mock_imports = []
napoleon_numpy_docstring = True
napoleon_google_docstring = False
suppress_warnings = ['myst.xref_missing', 'myst.header']

# -- HTML ----------------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_title = f'DrillPyze {release}'
html_static_path = []
html_theme_options = {'navigation_depth': 3}

# -- Markdown pages from the repository -----------------------------------------
GEN = HERE / '_generated'
PAGES = [('readme', ROOT / 'README.md', 'Overview (README)')]
PAGES += [('changelog', ROOT / 'CHANGELOG.md', 'Changelog')]
PAGES += [(f'derivation-{p.stem.lower()}', p, None) for p in sorted((ROOT / 'docs' / 'derivation').glob('*.md'))]
PAGES += [(f'spec-{p.stem.lower()}', p, None) for p in sorted((ROOT / 'docs' / 'spec').glob('*.md'))]


def _github_to_myst(text, src, dst):
    """GitHub maths -> MyST maths; relative image/file links rebased from src to dst."""
    # Sphinx wraps a display equation that contains line breaks (\\) in a split
    # environment, where \tag is not allowed (a cases block with a tag, say); such
    # equations go out unwrapped, inside an equation environment of their own.
    def fence(m):
        indent, body = m.group(1), m.group(2)
        if '\\tag' in body and '\\\\' in body:
            return f'{indent}```{{math}}\n{indent}:nowrap:\n\n{indent}\\begin{{equation}}\n{body}{indent}\\end{{equation}}\n{indent}```'
        return f'{indent}```{{math}}\n{body}{indent}```'
    text = re.sub(r'^([ \t]*)```math[ \t]*\n(.*?)^[ \t]*```', fence, text, flags=re.M | re.S)
    text = re.sub(r'\$`([^`]+)`\$', r'$\1$', text)

    def rebase(m):
        target = m.group(2)
        if re.match(r'^(https?:|mailto:|#|/)', target):
            return m.group(0)
        path, _, anchor = target.partition('#')
        p = (src.parent / path).resolve()
        if p.suffix.lower() == '.md':          # links to other Markdown pages -> their generated page
            for name, f, _ in PAGES:
                if f.resolve() == p:
                    return f'{m.group(1)}({name}.md{"#" + anchor if anchor else ""})'
        return f'{m.group(1)}({os.path.relpath(p, dst.parent)}{"#" + anchor if anchor else ""})'
    return re.sub(r'(!?\[[^\]]*\])\(([^)\s]+)\)', rebase, text)


def generate_pages(app=None):
    GEN.mkdir(exist_ok=True)
    for old in GEN.glob('*.*'):                    # drop pages left from another branch's build
        old.unlink()
    for name, src, _ in PAGES:
        dst = GEN / f'{name}.md'
        dst.write_text(_github_to_myst(src.read_text(), src, dst))
    # physics specification: listed when the branch has one (docs/spec), else a pointer
    spec = [name for name, src, _ in PAGES if name.startswith('spec-')]
    body = ('.. toctree::\n   :maxdepth: 1\n\n' + ''.join(f'   {n}\n' for n in spec)) if spec else (
        'The implementation-level physics specification (``docs/spec``) is not part of this branch.\n')
    (GEN / 'specification.rst').write_text('Physics specification\n=====================\n\n' + body)


generate_pages()


def setup(app):
    app.connect('builder-inited', generate_pages)
