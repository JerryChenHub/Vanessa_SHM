"""Femap COM bootstrap without a machine-specific installation path.

Use installed pywin32 by default; BWB_PYWIN32_PATH selects another runtime.
The old NASA_SHM workspace runtime remains a read-only fallback. FEMAP_TLB
overrides registry/install discovery. Import loads type metadata only; connect()
is the operation that attaches to or starts Femap. BWB owns its generated cache.
"""
from pathlib import Path
import os
import re
import site
import sys

ROOT = Path(__file__).resolve().parents[2]


def _add_runtime(path):
    """Process pywin32's .pth bootstrap and prefer this explicit runtime."""
    path = Path(os.path.expandvars(str(path))).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f'pywin32 runtime directory does not exist: {path}')
    preferred = [str(p) for p in (path, path / 'win32', path / 'win32/lib', path / 'pythonwin') if p.is_dir()]
    sys.path[:] = preferred + [entry for entry in sys.path if entry not in preferred]
    before = list(sys.path)
    site.addsitedir(str(path))
    added = [entry for entry in sys.path if entry not in before]
    priority = preferred + added
    sys.path[:] = priority + [entry for entry in sys.path if entry not in priority]


if os.environ.get('BWB_PYWIN32_PATH'):
    _add_runtime(os.environ['BWB_PYWIN32_PATH'])
else:
    try:
        import pythoncom
        import win32com.client
    except ImportError as exc:
        legacy_runtime = ROOT / 'workspace/NASA_SHM/femap_preview/python_deps'
        if not legacy_runtime.is_dir():
            raise ImportError(
                'Femap scripting requires pywin32 for this Python interpreter. '
                'Install with python -m pip install pywin32, or set '
                'BWB_PYWIN32_PATH to a compatible runtime directory.'
            ) from exc
        _add_runtime(legacy_runtime)

import pythoncom
import win32com
import win32com.client
import win32com.gen_py


def _registered_tlb_paths():
    """Read COM registration without creating a Femap application object."""
    import winreg

    def value(key):
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, key) as handle:
            return winreg.QueryValueEx(handle, '')[0]

    try:
        clsid = value(r'femap.model\CLSID')
        command = value(f'CLSID\\{clsid}\\LocalServer32')
        match = re.match(r'\s*(?:"([^"]+\.exe)"|(.+?\.exe))(?:\s|$)', command, re.I)
        if match:
            yield Path(os.path.expandvars(match.group(1) or match.group(2))).parent / 'femap.tlb'
    except OSError:
        pass

    key = r'TypeLib\{08F336B3-E668-11D4-9441-001083FFF11C}'
    try:
        with winreg.OpenKey(winreg.HKEY_CLASSES_ROOT, key) as handle:
            versions = [winreg.EnumKey(handle, i) for i in range(winreg.QueryInfoKey(handle)[0])]
        versions.sort(key=lambda version: tuple(int(n, 16) for n in version.split('.')), reverse=True)
        for version in versions:
            for platform in ('win64', 'win32'):
                try:
                    yield Path(os.path.expandvars(value(f'{key}\\{version}\\0\\{platform}')))
                except OSError:
                    pass
    except OSError:
        pass


def find_femap_tlb():
    override = os.environ.get('FEMAP_TLB')
    if override:
        path = Path(os.path.expandvars(override)).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f'FEMAP_TLB does not name an existing file: {path}')
        return path

    for path in _registered_tlb_paths():
        if path.is_file():
            return path
    for variable in ('ProgramFiles', 'ProgramFiles(x86)'):
        base = os.environ.get(variable)
        if base:
            for pattern in ('Siemens/Femap*/femap.tlb', 'Femap*/femap.tlb'):
                for path in sorted(Path(base).glob(pattern), reverse=True):
                    if path.is_file():
                        return path
    raise FileNotFoundError(
        'Femap type library was not found. Install/register Femap, or set '
        'FEMAP_TLB to the full path of its femap.tlb file.'
    )


CACHE = ROOT / 'workspace/BWB/femap_demo/com_cache'
CACHE.mkdir(parents=True, exist_ok=True)
win32com.__gen_path__ = str(CACHE)
win32com.gen_py.__path__ = [str(CACHE)]
TYPELIB_PATH = find_femap_tlb()
attrs = pythoncom.LoadTypeLib(str(TYPELIB_PATH.resolve())).GetLibAttr()
module = win32com.client.gencache.EnsureModule(attrs[0], attrs[1], attrs[3], attrs[4], bForDemand=False)

def connect(start=False):
    try:
        raw = win32com.client.GetActiveObject('femap.model')
    except pythoncom.com_error:
        if not start:
            raise
        raw = win32com.client.Dispatch('femap.model')
    return module.Imodel(raw._oleobj_)

def check(rc, action):
    if rc != -1:
        raise RuntimeError(f'{action}: Femap returned {rc}')
    return rc

def ids(entity_set):
    entity_set.Reset()
    out = []
    while True:
        n = entity_set.Next()
        if n <= 0:
            break
        out.append(n)
    return out
