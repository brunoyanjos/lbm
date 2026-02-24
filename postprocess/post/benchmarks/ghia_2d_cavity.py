# post/benchmarks/ghia_2d_cavity.py
from __future__ import annotations

import numpy as np

# Dados clássicos (Ghia et al.) para cavidade 2D:
# - u ao longo da vertical centerline (x=0.5): u(y) para vários y
# - v ao longo da horizontal centerline (y=0.5): v(x) para vários x
#
# IMPORTANTE:
# Estes valores são (u/U_lid) e (v/U_lid) no benchmark.
# Seu plot atual usa velocidades divididas por (2*U_lid), então no overlay
# nós vamos dividir por 2 também para cair no mesmo eixo.

_GHIA_U_VERTICAL = {
    # Re = 7500 (Table I)
    7500: {
        "y": np.array(
            [
                1.0000,
                0.9766,
                0.9688,
                0.9609,
                0.9531,
                0.8516,
                0.7344,
                0.6172,
                0.5000,
                0.4531,
                0.2813,
                0.1719,
                0.1016,
                0.0703,
                0.0625,
                0.0547,
                0.0000,
            ],
            dtype=np.float64,
        ),
        "u": np.array(
            [
                1.00000,
                0.47244,
                0.47048,
                0.47323,
                0.47167,
                0.34228,
                0.20591,
                0.08342,
                -0.03800,
                -0.07503,
                -0.23176,
                -0.32393,
                -0.38324,
                -0.43025,
                -0.43590,
                -0.43154,
                0.00000,
            ],
            dtype=np.float64,
        ),
    }
}

_GHIA_V_HORIZONTAL = {
    # Re = 7500 (Table II)
    7500: {
        "x": np.array(
            [
                1.0000,
                0.9688,
                0.9609,
                0.9531,
                0.9453,
                0.9063,
                0.8594,
                0.8047,
                0.5000,
                0.2344,
                0.2266,
                0.1563,
                0.0938,
                0.0781,
                0.0703,
                0.0625,
                0.0000,
            ],
            dtype=np.float64,
        ),
        "v": np.array(
            [
                0.00000,
                -0.53858,
                -0.55216,
                -0.52347,
                -0.48590,
                -0.41050,
                -0.36213,
                -0.30448,
                0.00824,
                0.27348,
                0.28117,
                0.35060,
                0.41824,
                0.43564,
                0.44030,
                0.43979,
                0.00000,
            ],
            dtype=np.float64,
        ),
    }
}


def ghia_centerlines(Re: float | int):
    """
    Returns (u_centerline, v_centerline) or (None, None) if Re not available.

    u_centerline: dict with keys ["y", "u"]
    v_centerline: dict with keys ["x", "v"]
    """
    # Re pode vir float tipo 7500.0
    Re_i = int(round(float(Re)))

    u = _GHIA_U_VERTICAL.get(Re_i, None)
    v = _GHIA_V_HORIZONTAL.get(Re_i, None)

    return u, v


def to_plot_coords_u(u_centerline: dict):
    """
    Convert benchmark u(y) to your current plot coordinates:
      y_star = y - 0.5
      u_star = (u/U_lid)/2   (because you currently divide by 2*U_lid)
    """
    y_star = u_centerline["y"] - 0.5
    u_star = u_centerline["u"] / 2.0
    return u_star, y_star


def to_plot_coords_v(v_centerline: dict):
    """
    Convert benchmark v(x) to your current plot coordinates:
      x_star = x - 0.5
      v_star = (v/U_lid)/2
    """
    x_star = v_centerline["x"] - 0.5
    v_star = v_centerline["v"] / 2.0
    return x_star, v_star
