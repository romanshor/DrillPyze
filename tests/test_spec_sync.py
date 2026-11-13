"""
The physics specification (docs/spec/physics-spec.md) must follow the Python model,
and the browser engine must implement the current spec revision.

If test_physics_sources_unchanged fails: review the change against the spec, update
physics-spec.md (and raise its revision) if the physics or numerics changed, then
run  python docs/spec/update_manifest.py .
"""
import importlib.util
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('update_manifest', ROOT / 'docs' / 'spec' / 'update_manifest.py')
um = importlib.util.module_from_spec(spec)
spec.loader.exec_module(um)


def test_physics_sources_unchanged():
    manifest = json.loads(um.MANIFEST.read_text())
    changed = [f for f, h in um.hashes().items() if manifest['files'].get(f) != h]
    assert not changed, ('physics source changed since the spec was last reviewed: ' + ', '.join(changed)
                         + '. Update docs/spec/physics-spec.md if the physics changed, then run '
                           'python docs/spec/update_manifest.py')


def test_manifest_revision_matches_spec():
    assert json.loads(um.MANIFEST.read_text())['spec_revision'] == um.spec_revision()


def test_lite_engine_implements_current_revision():
    src = (ROOT / 'lite' / 'engine' / 'spec.js').read_text()
    rev = int(re.search(r'SPEC_REVISION\s*=\s*(\d+)', src).group(1))
    assert rev == um.spec_revision(), (f'lite/engine implements spec revision {rev}, the spec is at '
                                       f'{um.spec_revision()}: update the browser engine')
