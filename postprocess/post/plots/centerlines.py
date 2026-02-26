# post/plots/centerlines.py
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from ..benchmarks.ghia_2d_cavity import (
    ghia_centerlines,
    to_plot_coords_u,
    to_plot_coords_v,
)

# Paleta oficial do projeto
PALETTE = [
    "#5EBD3E",
    "#FFB900",
    "#F78200",
    "#E23838",
    "#973999",
    "#009CDF",
]


def plot_centerlines(
    meta,
    ux_xc_y,
    uy_yc_x,
    *,
    u_lid=1.0,
    Re=None,
    show=True,
    savepath=None,
    title=None,
    show_benchmark=True,
):
    nx = meta["nx"]
    ny = meta["ny"]
    t = meta["t"]

    # coordenadas normalizadas centradas em zero
    x = np.arange(nx, dtype=np.float64)
    y = np.arange(ny, dtype=np.float64)

    x_star = x / nx - 0.5
    y_star = y / ny - 0.5

    u = float(u_lid) if float(u_lid) != 0.0 else 1.0
    ux_v = np.asarray(ux_xc_y, dtype=np.float64) / (2 * u)
    uy_h = np.asarray(uy_yc_x, dtype=np.float64) / (2 * u)

    fig, ax = plt.subplots(figsize=(7.0, 7.0))

    # ---- Solver curves ----
    ax.plot(
        ux_v,
        y_star,
        color=PALETTE[5],
        linewidth=2.0,
        label=r"$u_x(x_c,y)/U_{\mathrm{lid}}$ (solver)",
    )

    ax.plot(
        x_star,
        uy_h,
        color=PALETTE[3],
        linewidth=2.0,
        label=r"$u_y(x,y_c)/U_{\mathrm{lid}}$ (solver)",
    )

    # ---- Benchmark overlay (Ghia et al.) ----
    if show_benchmark and (Re is not None):
        u_b, v_b = ghia_centerlines(Re)

        # u(x=0.5, y) -> plot as (u*, y*)
        if u_b is not None:
            u_star_b, y_star_b = to_plot_coords_u(u_b)
            ax.scatter(
                u_star_b,
                y_star_b,
                s=28,
                marker="o",
                facecolors="none",
                edgecolors=PALETTE[0],
                linewidths=1.3,
                label="Benchmark (Ghia) u-centerline",
            )

        # v(x, y=0.5) -> plot as (x*, v*)
        if v_b is not None:
            x_star_b, v_star_b = to_plot_coords_v(v_b)
            ax.scatter(
                x_star_b,
                v_star_b,
                s=28,
                marker="s",
                facecolors="none",
                edgecolors=PALETTE[1],
                linewidths=1.3,
                label="Benchmark (Ghia) v-centerline",
            )

    ax.set_xlabel(r"Velocity (normalized)", fontsize=14)
    ax.set_ylabel(r"Position (normalized)", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect("equal", adjustable="box")

    if title is None:
        if Re is None:
            title = rf"Centerline profiles, $t={t}$"
        else:
            title = rf"Centerline profiles, $t={t}$, $Re={float(Re):g}$"
    ax.set_title(title, fontsize=13)

    ax.legend(frameon=False, fontsize=11, loc="best")
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()

    plt.close(fig)
