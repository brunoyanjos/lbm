# post/plots/couette_profile.py
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


def _couette_analytic(y: np.ndarray, H: float, Uwall: float) -> np.ndarray:
    """
    Couette plano (placa inferior parada, superior com Uwall):
      u(y) = Uwall * (y/H), y in [0,H]
    """
    y = np.asarray(y, dtype=np.float64)
    if H == 0.0:
        return y * 0.0
    return Uwall * (y / H)


def plot_couette_profile(
    sim_meta: dict,
    prof_meta: dict,
    ux_y: np.ndarray,
    *,
    show=True,
    savepath=None,
    title=None,
    show_analytic=True,
    normalize=True,
):
    ux_y = np.asarray(ux_y, dtype=np.float64)

    ny = int(prof_meta.get("ny", 0))
    if ny <= 1:
        raise ValueError("invalid ny for couette profile")

    delta = int(prof_meta.get("delta", 0))
    if delta < 0:
        raise ValueError("delta must be >= 0")

    # ux_y já vem cortado: primeiro ponto corresponde a y = delta
    y = np.arange(delta, delta + ny, dtype=np.float64)

    H = float(ny + 2 * delta - 1)
    if H <= 0.0:
        raise ValueError("invalid couette height (ny-1 must be > 0)")

    y_star = y / H

    # velocidade de referência (parede)
    U = sim_meta.get("U_wall", None)
    if U is None:
        U = sim_meta.get("U_lid", None)
    if U is None or float(U) == 0.0:
        U = 1.0
    U = float(U)

    # normalização
    if normalize:
        ux_plot = ux_y / U
        xlabel = r"$u_x/U_{\mathrm{wall}}$"
        ylabel = r"$y^{*}=y/H$"
    else:
        ux_plot = ux_y
        xlabel = r"$u_x$"
        ylabel = "y (lattice units)"

    # ordem por segurança
    order = np.argsort(y_star)
    y_star = y_star[order]
    ux_plot = ux_plot[order]

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    # Solver
    ax.plot(
        ux_plot,
        y_star,
        color=PALETTE[5],
        linewidth=1.8,
        marker="o",
        markersize=3.2,
        markevery=max(1, len(y_star) // 25),
        label="Solver (samples)",
        zorder=2,
    )

    # Analítico no domínio completo
    if show_analytic:
        y_full = np.linspace(0.0, H, 400, dtype=np.float64)
        y_full_star = y_full / H
        ua = _couette_analytic(y_full, H, U)
        if normalize:
            ua = ua / U

        ax.plot(
            ua,
            y_full_star,
            color=PALETTE[1],
            linewidth=2.0,
            linestyle="--",
            label="Analytical (Couette)",
            zorder=3,
        )

    # regiões excluídas perto das paredes: pintar no eixo Y
    y_gap_bottom = float(np.clip(y_star.min(), 0.0, 1.0))
    y_gap_top = float(np.clip(y_star.max(), 0.0, 1.0))

    if y_gap_bottom > 0.0:
        ax.axhspan(0.0, y_gap_bottom, alpha=0.10, color="gray", zorder=0)

    if y_gap_top < 1.0:
        ax.axhspan(y_gap_top, 1.0, alpha=0.10, color="gray", zorder=0)

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.tick_params(axis="both", which="major", labelsize=11)

    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_ylim(0.0, 1.0)

    if title is None:
        stencil = sim_meta.get("Stencil", "UNKNOWN")
        Re = sim_meta.get("Re", None)
        t = prof_meta.get("t", None)
        parts = [f"{stencil}", "Couette profile"]
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
