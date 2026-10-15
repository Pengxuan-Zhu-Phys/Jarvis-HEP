#!/usr/bin/env bash
# Complete HEP release: wheel + sdist, public dependencies, no extras/editable.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${JARVIS_RELEASE_PYTHON:-python3}"
output=""
operas_dist=""
upload=0
dry=0
requested_version=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --upload) upload=1; shift ;;
    --dry) dry=1; shift ;;
    --output|--operas-dist|--python)
      [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
      case "$1" in
        --output) output="$2" ;;
        --operas-dist) operas_dist="$2" ;;
        --python) python_bin="$2" ;;
      esac
      shift 2 ;;
    -h|--help)
      echo 'Usage: release [X.Y.Z] [--output DIR] [--operas-dist DIR] [--python PYTHON] [--upload]'
      echo 'Default: build wheel/sdist, install all dependencies in a clean venv, verify; no upload.'
      exit 0 ;;
    *)
      [[ -z "$requested_version" && "$1" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "Unknown argument: $1" >&2; exit 2; }
      requested_version="$1"; shift ;;
  esac
done
[[ "$dry" -eq 0 || "$upload" -eq 0 ]] || { echo '--dry and --upload cannot be combined' >&2; exit 2; }
# Resolve user paths before changing to the repository directory.
if [[ -n "$output" && "$output" != /* ]]; then output="$PWD/$output"; fi
if [[ -n "$operas_dist" && "$operas_dist" != /* ]]; then operas_dist="$PWD/$operas_dist"; fi
cd "$root"
work="$(mktemp -d "${TMPDIR:-/tmp}/jarvis-full-release.XXXXXXXX")"
trap 'rm -rf "$work"' EXIT
# Release tooling is separate from the distribution; every user installs the same complete package.
"$python_bin" -m venv "$work/tools"
build_python="$work/tools/bin/python"
"$build_python" -m pip install --disable-pip-version-check build twine tomli
if [[ -n "$requested_version" ]]; then
  "$build_python" - "$requested_version" <<'PY'
from pathlib import Path
import re, sys, tomli as tomllib
p = Path('pyproject.toml')
project = tomllib.loads(p.read_text())['project']
version = sys.argv[1]
assert project['name'] == 'Jarvis-HEP', 'This helper only releases Jarvis-HEP'
f = Path('jarvishep2/__init__.py')
updated, n = re.subn(r'(?m)^version\s*=\s*"[^"]+"$', f'version = "{version}"', p.read_text(), count=1)
runtime, m = re.subn(r'(?m)^__version__\s*=\s*"[^"]+"$', f'__version__ = "{version}"', f.read_text(), count=1)
if n != 1 or m != 1:
    raise SystemExit('Version fields not found; no files changed')
p.write_text(updated)
f.write_text(runtime)
PY
fi
# Snapshot tracked files, including their current edits, without unrelated untracked packages.
"$build_python" - "$root" "$work/source" <<'PY'
from pathlib import Path
import shutil, subprocess, sys, tomli as tomllib
root, target = map(Path, sys.argv[1:])
files = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0')
for rel in filter(None, files):
    p = root / rel
    if p.is_file():
        q = target / rel
        q.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, q)
project = tomllib.loads((target/'pyproject.toml').read_text())['project']
if project.get('optional-dependencies'):
    raise SystemExit('Only the complete distribution is allowed; remove dependency extras')
PY
version="$("$build_python" -c 'import tomli as tomllib; print(tomllib.load(open("pyproject.toml", "rb"))["project"]["version"])')"
package="$("$build_python" -c 'import tomli as tomllib; print(tomllib.load(open("pyproject.toml", "rb"))["project"]["name"])')"
output="${output:-$root/dist/$version}"
mkdir -p "$output"
output="$(cd "$output" && pwd)"
if compgen -G "$output/*.whl" >/dev/null || compgen -G "$output/*.tar.gz" >/dev/null; then
  echo "Output already contains distributions; choose a new --output directory: $output" >&2
  exit 2
fi
unset PYTHONPATH
"$build_python" -m build --outdir "$work/artifacts" "$work/source"
"$build_python" -m twine check "$work/artifacts/"*
"$python_bin" -m venv "$work/install"
install_python="$work/install/bin/python"
install_args=("$work/artifacts/"*.whl)
# An unpublished dependency patch can be verified locally before the maintainer
# publishes dependencies in order. The normal index-only path is unchanged.
if [[ -n "$operas_dist" ]]; then
  [[ -d "$operas_dist" ]] || { echo 'Missing Operas artifact directory' >&2; exit 2; }
  operas_wheels=("$operas_dist/"jarvis_operas-*.whl)
  [[ ${#operas_wheels[@]} -eq 1 && -f "${operas_wheels[0]}" ]] || { echo 'Expected exactly one Operas wheel in --operas-dist' >&2; exit 2; }
  install_args+=("${operas_wheels[0]}" --find-links "$operas_dist")
fi
# Deliberately resolve and install every declared dependency; never --no-deps or -e.
"$install_python" -m pip install --disable-pip-version-check "${install_args[@]}"
"$install_python" -m pip check
mkdir "$work/verify"
cd "$work/verify"
"$install_python" -I "$work/source/scripts/smoke_install.py" "$version"
if [[ "$package" == 'Jarvis-HEP' ]]; then
  "$install_python" -I "$work/source/scripts/smoke_agent_install.py" "$version"
  "$install_python" -I "$work/source/scripts/smoke_redis_install.py" "$version"
fi
if [[ -f "$work/source/tests/test_release_2015_regressions.py" ]]; then
  "$install_python" -m pip install --disable-pip-version-check pytest fakeredis colorlog
  JARVIS_OPERAS_CURVE_CACHE_ROOT="$work/verify/operas-cache" \
    "$install_python" -I -m pytest -q -o pythonpath= \
    "$work/source/tests/test_release_2015_regressions.py"
fi
"$install_python" - "$work/artifacts" "$work/source" <<'PY'
import email, hashlib, importlib.metadata as m, json, pathlib, re, sys, zipfile
out, source = map(pathlib.Path, sys.argv[1:])
text = (source/'pyproject.toml').read_text()
p = {k: re.search(r'(?m)^'+k+r'\s*=\s*"([^"]+)"', text).group(1) for k in ('name','version')}
wheels, sdists = list(out.glob('*.whl')), list(out.glob('*.tar.gz'))
assert len(wheels) == len(sdists) == 1, 'Expected one wheel and one sdist'
with zipfile.ZipFile(wheels[0]) as z:
    meta = email.message_from_bytes(z.read(next(n for n in z.namelist() if n.endswith('.dist-info/METADATA'))))
assert meta['Name'] == p['name'] and meta['Version'] == p['version']
assert not meta.get_all('Provides-Extra'), 'Extras are not allowed'
if p['name'] == 'Jarvis-HEP':
    assert any(n.startswith('jarvishep2/') for n in z.namelist())
    assert not any(n.startswith('jarvishep/') for n in z.namelist())
assert m.version(p['name']) == p['version']
report = {'distribution': p['name'], 'version': p['version'], 'installation': 'complete wheel with all dependencies',
          'checks': {'twine': 'passed', 'pip_check': 'passed', 'installed_smoke': 'passed',
                     'installed_agent_smoke': 'passed' if p['name'] == 'Jarvis-HEP' else 'not included',
                     'installed_redis_smoke': 'passed' if p['name'] == 'Jarvis-HEP' else 'not included',
                     'release_regressions': 'passed' if (source/'tests/test_release_2015_regressions.py').exists() else 'not included'},
          'artifacts': {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in wheels+sdists},
          'installed_dependencies': {d.metadata['Name']: d.version for d in m.distributions()}}
(out/'release-verification.json').write_text(json.dumps(report, sort_keys=True, indent=2)+'\n')
(out/'SHA256SUMS').write_text(''.join(f'{h}  {n}\n' for n,h in report['artifacts'].items()))
PY
cp "$work/artifacts/"* "$output/"
printf 'Verified complete artifacts: %s\n' "$output"
if [[ "$upload" -eq 1 ]]; then
  "$build_python" -m twine upload --non-interactive "$output/"*.whl "$output/"*.tar.gz
fi
