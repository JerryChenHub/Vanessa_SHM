"""BWB locations, independent of the current working directory."""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = PROJECT_ROOT.parent
DATA_DIR = PROJECT_ROOT / "data"
WORKSPACE_DIR = REPOSITORY_ROOT / "workspace" / "BWB"
