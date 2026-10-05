from __future__ import annotations

import argparse
from pathlib import Path

from common.binary import read_tke
from common.ids import read_run_specs
from common.logs import default_label, parse_run_config
from common.paths import DEFAULT_FIGURES_DIR, DEFAULT_IDS, DEFAULT_RUNS_DIR, stdout_path, tke_path
from common.plot_style import apply_style

import matplotlib.pyplot as plt


def plot_tke(ids: Path, runs_dir: Path, figures_dir: Path) -> None:
    apply_style()
    specs = read_run_specs(ids)
    out_dir = figures_dir / "tke"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    loaded = 0
    reference_re = None
    
    styles = [
        dict(linestyle="-", linewidth=2.0),   # Traditional
        dict(linestyle="--", linewidth=2.0),  # Isothermal
    ]

    for spec, style in zip(specs, styles):
        try:
            config = parse_run_config(spec.run_id, stdout_path(runs_dir, spec.run_id))
            tstar, ke = read_tke(tke_path(runs_dir, spec.run_id))
        except Exception as exc:
            print(f"[warning] skipping TKE for {spec.run_id}: {exc}")
            continue

        if reference_re is None and config.re is not None:
            reference_re = int(round(config.re))

        ax.plot(
            tstar,
            ke,
            label=spec.label or default_label(config),
            **style,
        )

        loaded += 1

    if loaded == 0:
        raise RuntimeError("no valid TKE data loaded")

    ax.set_xlabel(r"$t^*$")
    ax.set_ylabel(r"$\langle E_k^* \rangle$")
    ax.margins(x=0.02, y=0.08)
    if loaded > 1:
        ax.legend()

    re_tag = f"RE{reference_re}_" if reference_re is not None else ""

    png = out_dir / f"{re_tag}tke_comparison.png"
    pdf = out_dir / f"{re_tag}tke_comparison.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)

    print(f"Saved: {png}")
    print(f"Saved: {pdf}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot combined TKE from ids.txt.")
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    args = parser.parse_args()

    plot_tke(args.ids, args.runs_dir, args.figures_dir)


if __name__ == "__main__":
    main()
