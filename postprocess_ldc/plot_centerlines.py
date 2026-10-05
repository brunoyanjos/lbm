from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from common.binary import CenterlineData, read_centerline
from common.ids import RunSpec, read_run_specs
from common.logs import default_label, parse_run_config
from common.paths import DEFAULT_FIGURES_DIR, DEFAULT_IDS, DEFAULT_RUNS_DIR, centerline_path, stdout_path
from common.plot_style import apply_style

import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Ghia et al. (1982) reference data
# -----------------------------------------------------------------------------

GHIA_U = {
    100: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.84123, 0.78871, 0.73722, 0.68717,
            0.23151, 0.00332, -0.13641, -0.20581, -0.21090,
            -0.15662, -0.10150, -0.06434, -0.04775, -0.04192,
            -0.03717, 0.00000,
        ]),
    ),

    400: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.75837, 0.68439, 0.61756, 0.55892,
            0.29093, 0.16256, 0.02135, -0.11477, -0.17119,
            -0.32726, -0.24299, -0.14612, -0.10338, -0.09266,
            -0.08186, 0.00000,
        ]),
    ),

    1000: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.65928, 0.57492, 0.51117, 0.46604,
            0.33304, 0.18719, 0.05702, -0.06080, -0.10648,
            -0.27805, -0.38289, -0.29730, -0.22220, -0.20196,
            -0.18109, 0.00000,
        ]),
    ),

    3200: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.53236, 0.48296, 0.46547, 0.46101,
            0.34682, 0.19791, 0.07156, -0.04272, -0.08636,
            -0.24427, -0.34323, -0.41933, -0.37827, -0.35344,
            -0.32407, 0.00000,
        ]),
    ),

    5000: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.48223, 0.46120, 0.45992, 0.46036,
            0.33556, 0.20087, 0.08183, -0.03039, -0.07404,
            -0.22855, -0.33050, -0.40435, -0.43643, -0.42901,
            -0.41165, 0.00000,
        ]),
    ),

    7500: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.47244, 0.47048, 0.47323, 0.47167,
            0.34228, 0.20591, 0.08342, -0.03800, -0.07503,
            -0.23176, -0.32393, -0.38324, -0.43025, -0.43590,
            -0.43154, 0.00000,
        ]),
    ),

    10000: (
        np.array([
            1.0000, 0.9766, 0.9688, 0.9609, 0.9531,
            0.8516, 0.7344, 0.6172, 0.5000, 0.4531,
            0.2813, 0.1719, 0.1016, 0.0703, 0.0625,
            0.0547, 0.0000,
        ]),
        np.array([
            1.00000, 0.47221, 0.47783, 0.48070, 0.47804,
            0.34635, 0.20673, 0.08344, 0.03111, -0.07540,
            -0.23186, -0.32709, -0.38000, -0.41657, -0.42537,
            -0.42735, 0.00000,
        ]),
    ),
}

GHIA_V = {
    100: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.05906, -0.07391, -0.08864, -0.10313,
            -0.16914, -0.22445, -0.24533, 0.05454, 0.17527,
            0.17507, 0.16077, 0.12317, 0.10890, 0.10091,
            0.09233, 0.00000,
        ]),
    ),

    400: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.12146, -0.15663, -0.19254, -0.22847,
            -0.23827, -0.44993, -0.38598, 0.05186, 0.30174,
            0.30203, 0.28124, 0.22965, 0.20920, 0.19713,
            0.18360, 0.00000,
        ]),
    ),

    1000: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.21388, -0.27669, -0.33714, -0.39188,
            -0.51550, -0.42665, -0.31966, 0.02526, 0.32235,
            0.33075, 0.37095, 0.32627, 0.30353, 0.29012,
            0.27485, 0.00000,
        ]),
    ),

    3200: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.39017, -0.47425, -0.52357, -0.54053,
            -0.44307, -0.37401, -0.31184, 0.00999, 0.28188,
            0.29030, 0.37119, 0.42768, 0.41906, 0.40917,
            0.39560, 0.00000,
        ]),
    ),

    5000: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.49774, -0.55069, -0.55408, -0.52876,
            -0.41442, -0.36214, -0.30018, 0.00945, 0.27280,
            0.28066, 0.35368, 0.42951, 0.43648, 0.43329,
            0.42447, 0.00000,
        ]),
    ),

    7500: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.53858, -0.55216, -0.52347, -0.48590,
            -0.41050, -0.36213, -0.30448, 0.00824, 0.27348,
            0.28117, 0.35060, 0.41824, 0.43564, 0.44030,
            0.43979, 0.00000,
        ]),
    ),

    10000: (
        np.array([
            1.0000, 0.9688, 0.9609, 0.9531, 0.9453,
            0.9063, 0.8594, 0.8047, 0.5000, 0.2344,
            0.2266, 0.1563, 0.0938, 0.0781, 0.0703,
            0.0625, 0.0000,
        ]),
        np.array([
            0.00000, -0.54302, -0.52987, -0.49099, -0.45863,
            -0.41496, -0.36737, -0.30719, 0.00831, 0.27224,
            0.28003, 0.35070, 0.41487, 0.43124, 0.43733,
            0.43983, 0.00000,
        ]),
    ),
}


def _coords(data: CenterlineData) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(data.nx, dtype=np.float64)
    y = np.arange(data.ny, dtype=np.float64)
    return x / max(data.nx - 1, 1), y / max(data.ny - 1, 1)


def _normalise(values: np.ndarray, u_lid: float | None) -> np.ndarray:
    if u_lid is None or abs(u_lid) < 1.0e-14:
        return values
    return values / u_lid


def plot_centerlines(ids: Path, runs_dir: Path, figures_dir: Path) -> None:
    apply_style()
    specs = read_run_specs(ids)
    out_dir = figures_dir / "centerlines"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig_ux, ax_ux = plt.subplots(figsize=(6.4, 5.0))
    fig_uy, ax_uy = plt.subplots(figsize=(6.4, 5.0))
    loaded = 0
    reference_re = None

    styles = [
        dict(linestyle="-",  linewidth=2.0),   # Traditional
        dict(linestyle="--", linewidth=2.0),   # Isothermal
    ]

    for spec, style in zip(specs, styles):
        try:
            config = parse_run_config(spec.run_id, stdout_path(runs_dir, spec.run_id))
            data = read_centerline(centerline_path(runs_dir, spec.run_id))
        except Exception as exc:
            print(f"[warning] skipping centerline for {spec.run_id}: {exc}")
            continue

        x, y = _coords(data)
        label = spec.label or default_label(config)

        ax_ux.plot(_normalise(data.ux_xc_y, config.u_lid), y, label=label, **style)
        ax_uy.plot(x, _normalise(data.uy_yc_x, config.u_lid), label=label, **style)

        if reference_re is None and config.re is not None:
            reference_re = int(round(config.re))

        loaded += 1

    if loaded == 0:
        raise RuntimeError("no valid centerline data loaded")

    # Plota a referência de Ghia por último, para ficar acima das curvas numéricas
    if reference_re is not None:

        if reference_re in GHIA_U:
            yy, ux = GHIA_U[reference_re]
            ax_ux.plot(
                ux,
                yy,
                linestyle="none",
                marker="x",
                color="black",
                markersize=6,
                label="Ghia et al. (1982)",
            )
        else:
            print(f"[warning] No Ghia vertical data for Re={reference_re}")

        if reference_re in GHIA_V:
            xx, uy = GHIA_V[reference_re]
            ax_uy.plot(
                xx,
                uy,
                linestyle="none",
                marker="x",
                color="black",
                markersize=6,
                label="Ghia et al. (1982)",
            )
        else:
            print(f"[warning] No Ghia horizontal data for Re={reference_re}")

    ax_ux.set_xlabel(r"$u_x/U_{\mathrm{lid}}$")
    ax_ux.set_ylabel(r"$y/L$")
    #ax_ux.set_title("Vertical centerline")
    ax_ux.margins(x=0.06, y=0.02)

    ax_uy.set_xlabel(r"$x/L$")
    ax_uy.set_ylabel(r"$u_y/U_{\mathrm{lid}}$")
    #ax_uy.set_title("Horizontal centerline")
    ax_uy.margins(x=0.02, y=0.08)

    ax_ux.legend()
    ax_uy.legend()

    re_tag = f"RE{reference_re}_" if reference_re is not None else ""

    for fig, name in ((fig_ux, "centerline_ux_y"), (fig_uy, "centerline_uy_x")):
        png = out_dir / f"{re_tag}{name}.png"
        pdf = out_dir / f"{re_tag}{name}.pdf"
        fig.savefig(png)
        fig.savefig(pdf)
        print(f"Saved: {png}")
        print(f"Saved: {pdf}")
        plt.close(fig)

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot combined LDC centerlines from ids.txt.")
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    args = parser.parse_args()

    plot_centerlines(args.ids, args.runs_dir, args.figures_dir)


if __name__ == "__main__":
    main()
