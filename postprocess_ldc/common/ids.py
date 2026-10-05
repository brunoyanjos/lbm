from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RunSpec:
    run_id: str
    label: str | None = None


def read_run_specs(path: Path) -> list[RunSpec]:
    if not path.exists():
        raise FileNotFoundError(f"ids file not found: {path}")

    specs: list[RunSpec] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        if "|" in line:
            run_id, label = line.split("|", maxsplit=1)
            run_id = run_id.strip()
            label = label.strip() or None
        else:
            run_id = line
            label = None

        if run_id:
            specs.append(RunSpec(run_id=run_id, label=label))

    if not specs:
        raise ValueError(f"no run ids found in: {path}")

    return specs
