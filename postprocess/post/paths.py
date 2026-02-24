# post/paths.py
from __future__ import annotations

from pathlib import Path


def default_runs_dir() -> Path:
    # padrão compatível com o seu layout atual
    return Path("../runs")


def default_outcomes_dir() -> Path:
    return Path("outcomes")


def resolve_run_outputs(runs_dir: Path, run_id: str) -> Path:
    # runs/<run_id>/outputs
    return runs_dir / run_id / "outputs"


def resolve_out_dir(outcomes_dir: Path, run_id: str) -> Path:
    # outcomes/<run_id>/
    return outcomes_dir / run_id
