# post/io/read_bin.py
from __future__ import annotations

from pathlib import Path
import struct
import numpy as np

from .formats import (
    TKE_DTYPE,
    REAL_DTYPE,
    CLN_MAGIC,
    APR_MAGIC,
    APP_MAGIC,
    CHP_MAGIC,
    CUP_MAGIC,
    JCL_MAGIC,
    JSC_MAGIC,
)


def _read_exact(f, nbytes: int, what: str) -> bytes:
    b = f.read(nbytes)
    if len(b) != nbytes:
        raise ValueError(
            f"{what}: truncated file (expected {nbytes} bytes, got {len(b)})"
        )
    return b


def _expect_eof(f, what: str) -> None:
    extra = f.read(1)
    if extra:
        raise ValueError(f"{what}: file has extra bytes (format mismatch)")


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
      real ux_xc_y[ny]
      real uy_yc_x[nx]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "centerline magic")
        if magic != CLN_MAGIC:
            raise ValueError(f"Invalid centerline file (magic={magic!r})")

        nx, ny, xc, yc, t = struct.unpack(
            "<iiiii", _read_exact(f, 20, "centerline header")
        )

        ux_xc_y = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * ny, "ux_xc_y"), dtype=REAL_DTYPE
        ).copy()
        uy_yc_x = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * nx, "uy_yc_x"), dtype=REAL_DTYPE
        ).copy()

        _expect_eof(f, "centerline")

    meta = {"nx": nx, "ny": ny, "xc": xc, "yc": yc, "t": t}
    return meta, ux_xc_y, uy_yc_x


def read_annul_profile_bin(path: str | Path):
    """
    Format APR1:
      magic[4] = "APR1"
      int32 n, t, y_line, x_start
      real r[n]
      real ur[n]
      real utheta[n]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "annul_profile magic")
        if magic != APR_MAGIC:
            raise ValueError(f"Invalid annular profile file (magic={magic!r})")

        n, t, y_line, x_start = struct.unpack(
            "<iiii", _read_exact(f, 16, "annul_profile header")
        )

        r = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * n, "r"), dtype=REAL_DTYPE
        ).copy()
        ur = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * n, "ur"), dtype=REAL_DTYPE
        ).copy()
        ut = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * n, "utheta"), dtype=REAL_DTYPE
        ).copy()

        _expect_eof(f, "annul_profile")

    meta = {"n": n, "t": t, "y_line": y_line, "x_start": x_start}
    return meta, r, ur, ut


def read_annul_pressure_bin(path: str | Path):
    """
    Format APP1:
      magic[4] = "APP1"
      int32 n, t, y_line, x_start
      real r[n]
      real rho_prime[n]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "annul_pressure magic")
        if magic != APP_MAGIC:
            raise ValueError(f"Invalid annular pressure file (magic={magic!r})")

        n, t, y_line, x_start = struct.unpack(
            "<iiii", _read_exact(f, 16, "annul_pressure header")
        )

        r = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * n, "r"), dtype=REAL_DTYPE
        ).copy()
        rho_p = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * n, "rho_prime"), dtype=REAL_DTYPE
        ).copy()

        _expect_eof(f, "annul_pressure")

    meta = {"n": n, "t": t, "y_line": y_line, "x_start": x_start}
    return meta, r, rho_p


def read_channel_profile_bin(path: str | Path):
    """
    Format CHP1:
      magic[4] = "CHP1"
      int32 nx, ny, t, x_sample
      real ux_y[ny]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "channel_profile magic")
        if magic != CHP_MAGIC:
            raise ValueError(f"Invalid channel profile file (magic={magic!r})")

        nx, ny, t, x_sample = struct.unpack(
            "<iiii", _read_exact(f, 16, "channel_profile header")
        )
        ux_y = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * ny, "ux_y"), dtype=REAL_DTYPE
        ).copy()

        _expect_eof(f, "channel_profile")

    meta = {"nx": nx, "ny": ny, "t": t, "x_sample": x_sample}
    return meta, ux_y


def read_couette_profile_bin(path: str | Path):
    """
    Format CUP1:
      magic[4] = "CUP1"
      int32 nx, ny, t, x_sample
      real ux_y[ny]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "couette_profile magic")
        if magic != CUP_MAGIC:
            raise ValueError(f"Invalid couette profile file (magic={magic!r})")

        nx, ny, t, x_sample = struct.unpack(
            "<iiii", _read_exact(f, 16, "couette_profile header")
        )
        ux_y = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * ny, "ux_y"), dtype=REAL_DTYPE
        ).copy()

        _expect_eof(f, "couette_profile")

    meta = {"nx": nx, "ny": ny, "t": t, "x_sample": x_sample}
    return meta, ux_y


def read_jet_centerline_bin(path: str | Path):
    """
    Format JCL1:
      magic[4] = "JCL1"
      int32 nx, ny, t, y_line
      real ux_x[nx]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "jet_centerline magic")
        if magic != JCL_MAGIC:
            raise ValueError(f"Invalid jet centerline file (magic={magic!r})")

        nx, ny, t, y_line = struct.unpack(
            "<iiii", _read_exact(f, 16, "jet_centerline header")
        )
        ux_x = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * nx, "ux_x"), dtype=REAL_DTYPE
        ).copy()

        _expect_eof(f, "jet_centerline")

    meta = {"nx": nx, "ny": ny, "t": t, "y_line": y_line}
    return meta, ux_x


def read_jet_sections_bin(path: str | Path):
    """
    Format JSC1:
      magic[4] = "JSC1"
      int32 nx, ny, t, nsec
      int32 x_sec[nsec]
      real ux_sec[nsec][ny]
    """
    path = Path(path)
    with path.open("rb") as f:
        magic = _read_exact(f, 4, "jet_sections magic")
        if magic != JSC_MAGIC:
            raise ValueError(f"Invalid jet sections file (magic={magic!r})")

        nx, ny, t, nsec = struct.unpack(
            "<iiii", _read_exact(f, 16, "jet_sections header")
        )

        x_sec = np.frombuffer(_read_exact(f, 4 * nsec, "x_sec"), dtype="<i4").copy()

        ux_all = np.frombuffer(
            _read_exact(f, REAL_DTYPE.itemsize * nsec * ny, "ux_sec"),
            dtype=REAL_DTYPE,
        ).copy()
        ux_sec = ux_all.reshape((nsec, ny))

        _expect_eof(f, "jet_sections")

    meta = {"nx": nx, "ny": ny, "t": t, "nsec": int(nsec)}
    return meta, x_sec, ux_sec
