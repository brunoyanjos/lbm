# post/plots/jet.py
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

# Paleta oficial do projeto
PALETTE = [
    "#5EBD3E",
    "#FFB900",
    "#F78200",
    "#E23838",
    "#973999",
    "#009CDF",
]


def plot_jet_centerline(
    sim_meta: dict,
    meta: dict,
    ux_x: np.ndarray,
    *,
    show=True,
    savepath=None,
    title=None,
    normalize=True,
):
    ux_x = np.asarray(ux_x, dtype=np.float64)

    nx = int(meta.get("nx", ux_x.size))
    if nx <= 0:
        raise ValueError("invalid nx for jet centerline")
    if ux_x.size != nx:
        raise ValueError(f"ux_x size mismatch: ux_x.size={ux_x.size}, nx={nx}")

    x = np.arange(nx, dtype=np.float64)

    # normalização (se existir velocidade de referência)
    U = sim_meta.get("U_in", None)
    if U is None:
        U = sim_meta.get("U_lid", None)
    if U is None or float(U) == 0.0:
        U = 1.0
    U = float(U)

    if normalize:
        ux_plot = ux_x / U
        ylabel = r"$u_x/U_{\mathrm{ref}}$"
    else:
        ux_plot = ux_x
        ylabel = r"$u_x$"

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    ax.plot(
        x,
        ux_plot,
        color=PALETTE[5],
        linewidth=1.8,
        marker="o",
        markersize=3.0,
        markevery=max(1, len(x) // 25),
        label="Centerline (solver)",
    )

    ax.set_xlabel("x (lattice units)", fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)

    ax.tick_params(axis="both", which="major", labelsize=11)

    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title is None:
        stencil = sim_meta.get("Stencil", "UNKNOWN")
        Re = sim_meta.get("Re", None)
        t = meta.get("t", None)

        parts = [f"{stencil}", "Jet centerline"]
        if Re is not None:
            parts.append(f"Re={float(Re):g}")
        if t is not None:
            parts.append(f"t={int(t)}")

        title = " | ".join(parts)

    ax.set_title(title, fontsize=12)
    ax.legend(frameon=False, fontsize=11, loc="best")

    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    plt.close(fig)


def plot_jet_sections(
    sim_meta: dict,
    meta: dict,
    x_sec: np.ndarray,
    ux_sec: np.ndarray,
    *,
    show=True,
    savepath=None,
    title=None,
    normalize=True,
):
    ux_sec = np.asarray(ux_sec, dtype=np.float64)
    x_sec = np.asarray(x_sec, dtype=np.int32)

    if ux_sec.ndim != 2:
        raise ValueError("ux_sec must have shape (nsec, ny)")

    nsec, ny = ux_sec.shape
    y = np.arange(ny, dtype=np.float64)

    # velocidade de referência
    U = sim_meta.get("U_in", None)
    if U is None:
        U = sim_meta.get("U_lid", None)
    if U is None or float(U) == 0.0:
        U = 1.0
    U = float(U)

    if normalize:
        ux_plot = ux_sec / U
        xlabel = r"$u_x/U_{\mathrm{ref}}$"
    else:
        ux_plot = ux_sec
        xlabel = r"$u_x$"

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    for k in range(nsec):
        color = PALETTE[k % len(PALETTE)]

        ax.plot(
            ux_plot[k],
            y,
            color=color,
            linewidth=1.8,
            marker="o",
            markersize=3.0,
            markevery=max(1, len(y) // 25),
            label=f"x = {int(x_sec[k])}",
        )

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel("y (lattice units)", fontsize=13)

    ax.tick_params(axis="both", which="major", labelsize=11)

    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title is None:
        stencil = sim_meta.get("Stencil", "UNKNOWN")
        Re = sim_meta.get("Re", None)
        t = meta.get("t", None)

        parts = [f"{stencil}", "Jet sections"]
        if Re is not None:
            parts.append(f"Re={float(Re):g}")
        if t is not None:
            parts.append(f"t={int(t)}")

        title = " | ".join(parts)

    ax.set_title(title, fontsize=12)
    ax.legend(frameon=False, fontsize=11, loc="best")

    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")

    if show:
        plt.show()

    plt.close(fig)
