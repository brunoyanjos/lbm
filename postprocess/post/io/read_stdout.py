# post/io/read_stdout.py
from __future__ import annotations

from pathlib import Path
import re


def _extract(pattern: str, text: str, cast=float):
    m = re.search(pattern, text)
    if not m:
        return None
    return cast(m.group(1))


def read_stdout_summary(path: str | Path) -> dict:
    """
    Parse stdout.txt and return dictionary with simulation metadata.
    """
    path = Path(path)

    text = path.read_text()

    meta = {}

    # --- GPU ---
    m = re.search(r"GPU\s*:\s*(.+)", text)
    if m:
        meta["GPU"] = m.group(1).strip()

    meta["Compute Capability"] = _extract(r"Compute Capability:\s*([0-9.]+)", text)

    # --- Stencil ---
    m = re.search(r"Stencil\s*:\s*(\S+)", text)
    if m:
        meta["Stencil"] = m.group(1)

    meta["Q"] = _extract(r"Q\s*:\s*(\d+)", text, int)
    meta["cs2"] = _extract(r"cs\^2\s*:\s*([0-9.eE+-]+)", text)

    # --- Domain ---
    m = re.search(r"Domain size\s*:\s*(\d+)\s*x\s*(\d+)", text)
    if m:
        meta["nx"] = int(m.group(1))
        meta["ny"] = int(m.group(2))

    meta["Re"] = _extract(r"Re\s*:\s*([0-9.eE+-]+)", text)
    meta["U_lid"] = _extract(r"U_lid\s*:\s*([0-9.eE+-]+)", text)
    meta["L_char"] = _extract(r"L_char\s*:\s*([0-9.eE+-]+)", text)

    meta["nu"] = _extract(r"nu \(visc\)\s*:\s*([0-9.eE+-]+)", text)
    meta["tau"] = _extract(r"tau\s*:\s*([0-9.eE+-]+)", text)
    meta["omega"] = _extract(r"omega\s*:\s*([0-9.eE+-]+)", text)
    meta["rho_0"] = _extract(r"rho_0\s*:\s*([0-9.eE+-]+)", text)

    # --- Performance ---
    meta["MLUPS_GPU"] = _extract(r"MLUPS_GPU=([0-9.eE+-]+)", text)
    meta["MLUPS_WALL"] = _extract(r"MLUPS_WALL=([0-9.eE+-]+)", text)

    meta["timesteps"] = _extract(r"finished after\s*(\d+)\s*timesteps", text, int)

    return meta
