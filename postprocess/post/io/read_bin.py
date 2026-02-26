# post/io/read_bin.py
from __future__ import annotations

from pathlib import Path
import struct
import numpy as np

from .formats import TKE_DTYPE, CLN_MAGIC


def read_tke_bin(path: str | Path):
    """
    Returns:
        t_star (np.ndarray int64)
        ke     (np.ndarray float64)
    """
    path = Path(path)
    rec = np.fromfile(path, dtype=TKE_DTYPE)
    return rec["t_star"], rec["ke"]


def read_centerline_bin(path: str | Path):
    """
    Reduced format:

      magic[4] = "CLN1"
      int32 nx, ny, xc, yc, t
      float32 ux_xc_y[ny]
      float32 uy_yc_x[nx]

    Returns:
      meta dict, ux_xc_y (ny), uy_yc_x (nx)
    """
    path = Path(path)

    with path.open("rb") as f:
        magic = f.read(4)
        if magic != CLN_MAGIC:
            raise ValueError("Invalid centerline file")

        nx, ny, xc, yc, t = struct.unpack("<iiiii", f.read(20))

        # LER COMO FLOAT32 (real_t = float)
        ux_xc_y = np.frombuffer(f.read(4 * ny), dtype="<f4").copy()
        uy_yc_x = np.frombuffer(f.read(4 * nx), dtype="<f4").copy()

        leftover = f.read(1)
        if leftover:
            raise ValueError("Centerline file has extra bytes")

    meta = {"nx": nx, "ny": ny, "xc": xc, "yc": yc, "t": t}
    return meta, ux_xc_y, uy_yc_x
