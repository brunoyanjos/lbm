# post/plots/annul_profile.py
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

PALETTE = [
    "#5EBD3E",
    "#FFB900",
    "#F78200",
    "#E23838",
    "#973999",
    "#009CDF",
]


def _taylor_couette_utheta(r, Rin, Rout, Uwall):

    r = np.asarray(r, dtype=np.float64)

    denom = Rout * Rout - Rin * Rin

    A = -Uwall * Rin / denom
    B = Uwall * Rin * Rout * Rout / denom

    return A * r + B / r


def plot_annul_profile(
    sim_meta: dict,
    prof_meta: dict,
    r: np.ndarray,
    ur: np.ndarray,
    utheta: np.ndarray,
    *,
    show=True,
    savepath=None,
    title=None,
    show_ur=True,
    show_analytic=True,
):
    nx = int(sim_meta.get("nx", 0))
    if nx <= 0:
        raise ValueError("sim_meta missing nx")

    Rin = nx * 0.25
    Rout = (nx - 1) * 0.5
    gap = Rout - Rin
    if gap <= 0:
        raise ValueError("Invalid annular gap computed from nx")

    U = sim_meta.get("U_lid", None)
    if U is None or float(U) == 0.0:
        U = 1.0
    U = float(U)

    r = np.asarray(r, dtype=np.float64)
    ur = np.asarray(ur, dtype=np.float64)
    utheta = np.asarray(utheta, dtype=np.float64)

    # normalizações
    r_star = (r - Rin) / gap
    ur_star = ur / U
    ut_star = utheta / U

    # garante ordenação por r (importante se você fez push_back pulando nós)
    order = np.argsort(r_star)
    r_star = r_star[order]
    ur_star = ur_star[order]
    ut_star = ut_star[order]

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    # -------------------------
    # Solver: somente pontos disponíveis
    # -------------------------
    ax.plot(
        r_star,
        ut_star,
        color=PALETTE[5],
        linewidth=1.6,
        marker="o",
        markersize=3.2,
        label=r"$u_{\theta}/U_{\mathrm{wall}}$ (solver samples)",
        zorder=1,
    )

    # -------------------------
    # Analítico: domínio completo
    # -------------------------
    if show_analytic:
        r_full = np.linspace(Rin, Rout, 400, dtype=np.float64)
        r_full_star = (r_full - Rin) / gap
        ut_ana = _taylor_couette_utheta(r_full, Rin, Rout, U) / U

        ax.plot(
            r_full_star,
            ut_ana,
            color=PALETTE[1],
            linewidth=2.0,
            linestyle="--",
            label="Analytical (full domain)",
            zorder=3,
        )

    # -------------------------
    # Evidenciar regiões faltantes perto das paredes
    # -------------------------
    if r_star.size > 0:
        left = float(np.clip(r_star.min(), 0.0, 1.0))
        right = float(np.clip(r_star.max(), 0.0, 1.0))

        if left > 0.0:
            ax.axvspan(0.0, left, alpha=0.10, zorder=0)

        if right < 1.0:
            ax.axvspan(right, 1.0, alpha=0.10, zorder=0)

    ax.set_xlabel(r"$r^{*}=(r-R_i)/(R_o-R_i)$", fontsize=13)
    ax.set_ylabel("Velocity (normalized)", fontsize=13)
    ax.tick_params(axis="both", which="major", labelsize=11)
    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0.0, 1.0)

    if title is None:
        stencil = sim_meta.get("Stencil", "UNKNOWN")
        Re = sim_meta.get("Re", None)
        t = prof_meta.get("t", None)
        parts = [f"{stencil}", "Annular profile"]
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
