"""
Refresh docs/spec/manifest.json: the spec revision and a hash of every
physics source file, so tests/test_spec_sync.py can tell when the Python
physics changed without the specification being reviewed.

Run after updating docs/spec/physics-spec.md (or after confirming that a
source change does not alter the physics):

    python docs/spec/update_manifest.py
"""
import hashlib
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
SPEC = ROOT / 'docs' / 'spec' / 'physics-spec.md'
MANIFEST = ROOT / 'docs' / 'spec' / 'manifest.json'
# files whose content defines the physics and numerics (config.py is excluded: it changes for
# UI reasons; its build order is specified in sections 4 and 7)
PHYSICS = ['kernel.py', 'simulation.py', 'geometry.py', 'survey.py', 'drillstring.py', 'controllers.py',
           'hydraulics.py', 'thermal.py', 'formation.py', 'statics.py', 'parameters.py', 'wellbore.py',
           'units.py']


def spec_revision(text=None):
    text = SPEC.read_text() if text is None else text
    m = re.search(r'\*\*Spec revision: (\d+)\*\*', text)
    if not m:
        raise ValueError('physics-spec.md has no "**Spec revision: N**" line')
    return int(m.group(1))


def hashes():
    src = ROOT / 'src' / 'DrillPyze'
    # normalise line endings so checkouts on any platform hash the same
    return {f: hashlib.sha256((src / f).read_bytes().replace(b'\r\n', b'\n')).hexdigest() for f in PHYSICS}


def build():
    return {'spec_revision': spec_revision(), 'files': hashes()}


if __name__ == '__main__':
    MANIFEST.write_text(json.dumps(build(), indent=2) + '\n')
    print(f'wrote {MANIFEST.relative_to(ROOT)} (spec revision {spec_revision()})')
