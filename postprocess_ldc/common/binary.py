from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np

TKE_DTYPE = np.dtype([("tstar", "<i8"), ("ke", "<f8")])
CENTERLINE_MAGIC = b"CLN1"
CENTERLINE_HEADER_SIZE = 24


@dataclass(frozen=True)
class CenterlineData:
    nx: int
    ny: int
    xc: int
    yc: int
    step: int
    ux_xc_y: np.ndarray
    uy_yc_x: np.ndarray


def read_tke(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"TKE file not found: {path}")

    size = path.stat().st_size
    if size == 0:
        raise ValueError(f"TKE file is empty: {path}")
    if size % TKE_DTYPE.itemsize != 0:
        raise ValueError(f"invalid TKE binary size: {path}")

    data = np.fromfile(path, dtype=TKE_DTYPE)
    tstar = data["tstar"].astype(np.int64)
    ke = data["ke"].astype(np.float64)

    if np.any(~np.isfinite(ke)):
        raise ValueError(f"non-finite TKE values found in: {path}")

    return tstar, ke


def read_centerline(path: Path) -> CenterlineData:
    if not path.exists():
        raise FileNotFoundError(f"centerline file not found: {path}")

    size = path.stat().st_size
    if size < CENTERLINE_HEADER_SIZE:
        raise ValueError(f"centerline file is too small: {path}")

    with path.open("rb") as stream:
        magic = stream.read(4)
        if magic != CENTERLINE_MAGIC:
            raise ValueError(f"invalid centerline magic in {path}: {magic!r}")

        nx, ny, xc, yc, step = struct.unpack("<5i", stream.read(20))
        payload_size = size - CENTERLINE_HEADER_SIZE
        n_values = nx + ny

        if nx <= 0 or ny <= 0:
            raise ValueError(f"invalid centerline grid in {path}: {nx}x{ny}")
        if payload_size % n_values != 0:
            raise ValueError(f"invalid centerline payload size: {path}")

        real_size = payload_size // n_values
        if real_size == 4:
            dtype = np.dtype("<f4")
        elif real_size == 8:
            dtype = np.dtype("<f8")
        else:
            raise ValueError(f"unsupported centerline real size {real_size} in {path}")

        ux_xc_y = np.fromfile(stream, dtype=dtype, count=ny).astype(np.float64)
        uy_yc_x = np.fromfile(stream, dtype=dtype, count=nx).astype(np.float64)

    if ux_xc_y.size != ny or uy_yc_x.size != nx:
        raise ValueError(f"could not read full centerline arrays from: {path}")
    if np.any(~np.isfinite(ux_xc_y)) or np.any(~np.isfinite(uy_yc_x)):
        raise ValueError(f"non-finite centerline values found in: {path}")

    return CenterlineData(nx, ny, xc, yc, step, ux_xc_y, uy_yc_x)
