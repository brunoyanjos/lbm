from __future__ import annotations

import argparse
from pathlib import Path

from common.paths import DEFAULT_FIGURES_DIR, DEFAULT_IDS, DEFAULT_RUNS_DIR
from plot_centerlines import plot_centerlines
from plot_tke import plot_tke
from plot_vti_velocity import plot_vti_velocity


def main() -> None:
    parser = argparse.ArgumentParser(description="Run all LDC post-processing plots.")
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--vti-which", choices=("latest", "first", "all"), default="latest")
    parser.add_argument("--streamlines", action="store_true")
    args = parser.parse_args()

    plot_centerlines(args.ids, args.runs_dir, args.figures_dir)
    plot_tke(args.ids, args.runs_dir, args.figures_dir)
    plot_vti_velocity(
        args.ids,
        args.runs_dir,
        args.figures_dir,
        which=args.vti_which,
        streamlines=args.streamlines,
    )


if __name__ == "__main__":
    main()
