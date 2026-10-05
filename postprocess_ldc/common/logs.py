from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RunConfig:
    run_id: str
    stencil: str | None = None
    nx: int | None = None
    ny: int | None = None
    re: float | None = None
    u_lid: float | None = None


def _extract(text: str, label: str) -> str | None:
    match = re.search(rf"^{re.escape(label)}\s*:\s*(.+?)\s*$", text, re.MULTILINE)
    return match.group(1).strip() if match else None


def _extract_float(text: str, label: str) -> float | None:
    value = _extract(text, label)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _extract_domain(text: str) -> tuple[int | None, int | None]:
    match = re.search(r"^Domain size\s*:\s*(\d+)\s*x\s*(\d+)\s*$", text, re.MULTILINE)
    if not match:
        return None, None
    return int(match.group(1)), int(match.group(2))


def _extract_re_from_run_id(run_id: str) -> float | None:
    match = re.search(r"(?:^|_)RE([0-9]+(?:p[0-9]+)?)(?:_|$)", run_id)
    if not match:
        return None
    try:
        return float(match.group(1).replace("p", "."))
    except ValueError:
        return None


def _extract_grid_from_run_id(run_id: str) -> tuple[int | None, int | None]:
    match = re.search(r"(?:^|_)(\d+)x(\d+)(?:_|$)", run_id)
    if not match:
        return None, None
    return int(match.group(1)), int(match.group(2))


def parse_run_config(run_id: str, stdout: Path) -> RunConfig:
    nx_from_id, ny_from_id = _extract_grid_from_run_id(run_id)

    if not stdout.exists():
        return RunConfig(
            run_id=run_id,
            nx=nx_from_id,
            ny=ny_from_id,
            re=_extract_re_from_run_id(run_id),
        )

    text = stdout.read_text(encoding="utf-8", errors="replace")
    nx, ny = _extract_domain(text)

    return RunConfig(
        run_id=run_id,
        stencil=_extract(text, "Stencil"),
        nx=nx if nx is not None else nx_from_id,
        ny=ny if ny is not None else ny_from_id,
        re=_extract_float(text, "Re") or _extract_re_from_run_id(run_id),
        u_lid=_extract_float(text, "U_lid"),
    )


def default_label(config: RunConfig) -> str:
    pieces: list[str] = []
    if config.stencil:
        pieces.append(config.stencil)
    if config.nx and config.ny:
        pieces.append(f"{config.nx}x{config.ny}")
    if config.re is not None:
        pieces.append(f"Re={config.re:g}")
    return ", ".join(pieces) if pieces else config.run_id
