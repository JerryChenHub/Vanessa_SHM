"""Project locations, independent of the current working directory.

Default inputs belong to NASA_SHM; generated files belong to workspace/NASA_SHM.
Explicit relative paths passed by a caller remain relative to that caller.
"""
from pathlib import Path

CORE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CORE_DIR.parent
REPOSITORY_ROOT = PROJECT_ROOT.parent
BASE_DIR = CORE_DIR / "base"
DATA_DIR = CORE_DIR / "data"
ASSETS_DIR = CORE_DIR / "assets"
WORKSPACE_DIR = REPOSITORY_ROOT / "workspace" / "NASA_SHM"
CASES_DIR = WORKSPACE_DIR / "cases"
FIGURES_DIR = WORKSPACE_DIR / "figures"


def require_working_output(path):
    """Prevent generated cases and solver outputs from overwriting core files."""
    path = Path(path).resolve()
    if path == CORE_DIR or CORE_DIR in path.parents:
        raise ValueError(f"Generated output must be outside core: {path}")
    return path
