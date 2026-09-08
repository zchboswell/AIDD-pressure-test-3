import importlib
import json
from pathlib import Path

environment = json.loads(Path('admin/ENVIRONMENT.json').read_text())
checks = {}
for name, expected in environment['packages'].items():
    try:
        module = importlib.import_module(name)
        version = module.cmd.get_version()[0] if name == 'pymol' else getattr(module, '__version__', None)
        checks[name] = dict(import_ok=True, version=version, path=module.__file__, expected_version=expected['version'], version_matches=version == expected['version'])
    except Exception as error:
        checks[name] = dict(import_ok=False, error=str(error))
Path('admin/RUNTIME_CHECKS.json').write_text(json.dumps(checks, indent=2) + '\n')
print(json.dumps(checks, indent=2))
