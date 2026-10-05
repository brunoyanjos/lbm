from __future__ import annotations

import re
from pathlib import Path

POSTPROCESS_DIR = Path(__file__).resolve().parents[1]
ROOT_DIR = POSTPROCESS_DIR.parent

DEFAULT_IDS = POSTPROCESS_DIR / "ids.txt"
DEFAULT_RUNS_DIR = ROOT_DIR / "runs"
DEFAULT_FIGURES_DIR = POSTPROCESS_DIR / "figures"


def safe_filename(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return value.strip("_") or "run"


def run_dir(runs_dir: Path, run_id: str) -> Path:
    return runs_dir / run_id


def centerline_path(runs_dir: Path, run_id: str) -> Path:
    return run_dir(runs_dir, run_id) / "outputs" / "centerline.bin"


def tke_path(runs_dir: Path, run_id: str) -> Path:
    return run_dir(runs_dir, run_id) / "outputs" / "tke.bin"


def vtk_dir(runs_dir: Path, run_id: str) -> Path:
    return run_dir(runs_dir, run_id) / "vtk"


def stdout_path(runs_dir: Path, run_id: str) -> Path:
    return run_dir(runs_dir, run_id) / "logs" / "stdout.txt"
