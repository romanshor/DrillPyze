"""
Documentation hygiene: maths must be written in GitHub's protected forms so
Markdown cannot strip backslashes or braces before it is typeset
(``\\{`` -> ``{``, ``\\,`` -> ``,``, ``\\\\`` -> ``\\``, ``_..._`` -> italics):

* display maths in fenced  ```math  blocks (no ``$$``);
* inline maths as ``$`...`$``;
* code fences balanced.
"""
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOCS = ([ROOT / 'README.md'] + sorted((ROOT / 'docs' / 'derivation').glob('*.md'))
        + sorted((ROOT / 'docs' / 'spec').glob('*.md')))


@pytest.mark.parametrize('path', DOCS, ids=[str(p.relative_to(ROOT)) for p in DOCS])
def test_maths_is_github_safe(path):
    fence = False
    problems = []
    for i, line in enumerate(path.read_text().split('\n'), 1):
        if line.strip().startswith('```'):
            fence = not fence
            continue
        if fence:
            continue
        rest = re.sub(r'\$`[^`]*`\$', '', line)     # protected inline maths
        rest = re.sub(r'`[^`]*`', '', rest)         # code spans
        if '$' in rest:
            problems.append(f'{i}: {line.strip()[:80]}')
    assert not fence, 'unbalanced code fence'
    assert not problems, 'unprotected $ maths (use ```math blocks or $`...`$):\n' + '\n'.join(problems)
