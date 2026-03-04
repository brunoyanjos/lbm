# post/plots/channel_profile.py
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


def _poiseuille_analytic(y: np.ndarray, H: float, Umax: float) -> np.ndarray:
    """
    Plano canal (Poiseuille):
      u(y) = 4 Umax (y/H) (1 - y/H)
    onde y in [0, H].
    """
    y = np.asarray(y, dtype=np.float64)
    eta = y / H if H != 0.0 else y * 0.0
    return 4.0 * Umax * eta * (1.0 - eta)


def plot_channel_profile(
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

    ny = int(prof_meta.get("ny", ux_y.size))
    if ny <= 0:
        raise ValueError("invalid ny for channel profile")
    if ux_y.size != ny:
        # evita gráfico “meio certo” caso bin truncado/meta errado
        raise ValueError(f"ux_y size mismatch: ux_y.size={ux_y.size}, ny={ny}")

    # coordenada y (nó) e y* (normalizada)
    y = np.arange(ny, dtype=np.float64)
    H = float(ny - 1)
    if H <= 0.0:
        raise ValueError("invalid channel height (ny-1 must be > 0)")
    y_star = y / H

    # velocidade de referência:
    # - tenta U_in (canal)
    # - fallback U_lid (se você reaproveitar logs)
    # - fallback 1.0
    U = sim_meta.get("U_in", None)
    if U is None:
        U = sim_meta.get("U_lid", None)
    if U is None or float(U) == 0.0:
        U = 1.0
    U = float(U)

    # Normalização
    if normalize:
        ux_plot = ux_y / U
        ylabel = r"$y^{*}=y/H$"
        xlabel = r"$u_x/U_{\mathrm{ref}}$"
    else:
        ux_plot = ux_y
        ylabel = "y (lattice units)"
        xlabel = r"$u_x$"

    # ordenação (por segurança, embora y já seja crescente)
    order = np.argsort(y_star)
    y_star = y_star[order]
    ux_plot = ux_plot[order]

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    # -------------------------
    # Solver curve
    # -------------------------
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

    # -------------------------
    # Analítico (parábola)
    # -------------------------
    if show_analytic:
        ua = _poiseuille_analytic(y, H, Umax=U)
        if normalize:
            ua = ua / U

        ax.plot(
            ua[order],
            y_star,
            color=PALETTE[1],
            linewidth=2.0,
            linestyle="--",
            label="Analytical (Poiseuille)",
            zorder=3,
        )

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.tick_params(axis="both", which="major", labelsize=11)

    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Canal: y* em [0,1]
    ax.set_ylim(0.0, 1.0)

    # (opcional) colocar y=0 embaixo como “físico”
    ax.invert_yaxis()

    if title is None:
        stencil = sim_meta.get("Stencil", "UNKNOWN")
        Re = sim_meta.get("Re", None)
        t = prof_meta.get("t", None)

        parts = [f"{stencil}", "Channel profile"]
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
