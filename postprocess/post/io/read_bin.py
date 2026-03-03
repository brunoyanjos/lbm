# post/io/read_bin.py
from __future__ import annotations

from pathlib import Path
import struct
import numpy as np

from .formats import TKE_DTYPE, CLN_MAGIC, APR_MAGIC, APP_MAGIC


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


def read_annul_profile_bin(path: str | Path):
    """
    Format APR1:

      magic[4] = "APR1"
      int32 n, t, y_line, x_start
      float32 r[n]
      float32 ur[n]
      float32 utheta[n]

    Returns:
      meta dict, r (n), ur (n), utheta (n)
    """
    path = Path(path)

    with path.open("rb") as f:
        magic = f.read(4)
        if magic != APR_MAGIC:
            raise ValueError("Invalid annular profile file")

        n, t, y_line, x_start = struct.unpack("<iiii", f.read(16))

        r = np.frombuffer(f.read(4 * n), dtype="<f4").copy()
        ur = np.frombuffer(f.read(4 * n), dtype="<f4").copy()
        utheta = np.frombuffer(f.read(4 * n), dtype="<f4").copy()

        leftover = f.read(1)
        if leftover:
            raise ValueError("Annular profile file has extra bytes")

    meta = {"n": n, "t": t, "y_line": y_line, "x_start": x_start}
    return meta, r, ur, utheta


def read_annul_pressure_bin(path: str | Path):
    """
    Format APP1:

      magic[4] = "APP1"
      int32 n, t, y_line, x_start
      float32 r[n]
      float32 rho_prime[n]   (rho - rho0)

    Returns:
      meta dict, r (n), rho_prime (n)
    """
    path = Path(path)

    with path.open("rb") as f:
        magic = f.read(4)
        if magic != APP_MAGIC:
            raise ValueError("Invalid annular pressure file")

        n, t, y_line, x_start = struct.unpack("<iiii", f.read(16))

        r = np.frombuffer(f.read(4 * n), dtype="<f4").copy()
        rho_p = np.frombuffer(f.read(4 * n), dtype="<f4").copy()

        leftover = f.read(1)
        if leftover:
            raise ValueError("Annular pressure file has extra bytes")

    meta = {"n": n, "t": t, "y_line": y_line, "x_start": x_start}
    return meta, r, rho_p
