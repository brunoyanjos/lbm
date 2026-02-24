# post/plots/tke.py
from __future__ import annotations

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


def plot_tke(t_star, ke, *, title=None, show=True, savepath=None):
    fig, ax = plt.subplots()

    ax.plot(
        t_star,
        ke,
        color=PALETTE[5],  # azul principal
        linewidth=2.0,
    )

    ax.set_xlabel(r"$t^{*}$", fontsize=14)
    ax.set_ylabel(r"$E_{k}^{*}$", fontsize=14)

    ax.tick_params(axis="both", which="major", labelsize=12)

    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if title:
        ax.set_title(title, fontsize=14)

    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")
    if show:
        plt.show()

    plt.close(fig)
