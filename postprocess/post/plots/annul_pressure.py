from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt

PALETTE = ["#5EBD3E", "#FFB900", "#F78200", "#E23838", "#973999", "#009CDF"]


def _taylor_couette_p(r, Rin, Rout, Uwall, rho0=1.0):
    """
    p(r) = rho * [ (A^2/2) r^2 + 2AB ln r - (B^2)/(2 r^2) ] + C
    Retorna sem C. Depois a gente remove uma referência para trabalhar com Δp.
    """
    r = np.asarray(r, dtype=np.float64)
    denom = Rout * Rout - Rin * Rin
    A = -Uwall * Rin / denom
    B = Uwall * Rin * Rout * Rout / denom

    return rho0 * (
        0.5 * A * A * r * r + 2.0 * A * B * np.log(r) - 0.5 * (B * B) / (r * r)
    )


def plot_annul_pressure(
    sim_meta: dict,
    prof_meta: dict,
    r: np.ndarray,
    rho_prime: np.ndarray,
    *,
    show=True,
    savepath=None,
    title=None,
    show_analytic=True,
):
    # metadados relevantes
    cs2 = sim_meta.get("cs2", None)
    if cs2 is None:
        raise ValueError("stdout meta missing cs2 (needed to compute p' = cs2*rho')")
    cs2 = float(cs2)

    rho0 = sim_meta.get("rho_0", 1.0)
    rho0 = float(rho0) if rho0 is not None else 1.0

    U = sim_meta.get("U_lid", 1.0)
    U = float(U) if U is not None and float(U) != 0.0 else 1.0

    r = np.asarray(r, dtype=np.float64)
    rho_p = np.asarray(rho_prime, dtype=np.float64)

    nx = int(sim_meta.get("nx", 0))
    if nx <= 0:
        raise ValueError("sim_meta missing nx")

    Rin = nx * 0.25
    Rout = (nx - 1) * 0.5
    gap = Rout - Rin
    if gap <= 0:
        raise ValueError("Invalid annular gap")

    # ---- ordena por r (importantíssimo agora) ----
    order = np.argsort(r)
    r = r[order]
    rho_p = rho_p[order]

    r_star = (r - Rin) / gap

    # pressão desviatória do LBM: p' = cs2 * rho'
    p_prime = cs2 * rho_p

    if p_prime.size == 0:
        raise ValueError("empty pressure profile (no fluid samples saved)")

    # trabalhar com Δp' removendo constante: referência = primeiro ponto disponível (menor r)
    dp = p_prime - p_prime[0]

    # normalização: dp* = dp / (rho0 * U^2)
    denom = rho0 * (U * U)
    dp_star = dp / denom if denom != 0.0 else dp

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    # -------------------------
    # Solver: somente amostras disponíveis
    # -------------------------
    ax.plot(
        r_star,
        dp_star,
        color=PALETTE[5],
        linewidth=1.6,
        marker="o",
        markersize=3.2,
        markevery=max(1, len(r_star) // 25),
        label=r"$\Delta p^{*}$ (solver samples)",
        zorder=1,
    )

    # -------------------------
    # Analítico: domínio completo
    # -------------------------
    if show_analytic:
        r_full = np.linspace(Rin, Rout, 500, dtype=np.float64)
        r_full_star = (r_full - Rin) / gap

        p_ana_full = _taylor_couette_p(r_full, Rin, Rout, U, rho0=rho0)

        # Referencia o analítico no MESMO r de referência do solver (r[0])
        p_ana_ref = _taylor_couette_p(np.array([r[0]]), Rin, Rout, U, rho0=rho0)[0]
        dp_ana_full = p_ana_full - p_ana_ref
        dp_ana_star = dp_ana_full / denom if denom != 0.0 else dp_ana_full

        ax.plot(
            r_full_star,
            dp_ana_star,
            color=PALETTE[1],
            linewidth=2.0,
            linestyle="--",
            label="Analytical (full domain)",
            zorder=3,
        )

    # -------------------------
    # Evidenciar regiões faltantes
    # -------------------------
    left = float(np.clip(r_star.min(), 0.0, 1.0))
    right = float(np.clip(r_star.max(), 0.0, 1.0))

    if left > 0.0:
        ax.axvspan(0.0, left, alpha=0.10, zorder=0)

    if right < 1.0:
        ax.axvspan(right, 1.0, alpha=0.10, zorder=0)

    ax.set_xlabel(r"$r^{*}=(r-R_i)/(R_o-R_i)$", fontsize=13)
    ax.set_ylabel(r"$\Delta p^{*}=\Delta p/(\rho_0 U_{\mathrm{wall}}^2)$", fontsize=13)
    ax.tick_params(axis="both", which="major", labelsize=11)

    ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0.0, 1.0)

    if title is None:
        stencil = sim_meta.get("Stencil", "UNKNOWN")
        Re = sim_meta.get("Re", None)
        t = prof_meta.get("t", None)
        parts = [f"{stencil}", "Annular pressure"]
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
