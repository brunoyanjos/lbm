from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from common.ids import read_run_specs
from common.logs import default_label, parse_run_config
from common.paths import (
    DEFAULT_FIGURES_DIR,
    DEFAULT_IDS,
    DEFAULT_RUNS_DIR,
    safe_filename,
    stdout_path,
    vtk_dir,
)
from common.plot_style import apply_style
from common.vti import read_vti_velocity, select_vti_files

import matplotlib.pyplot as plt


def _plot_one(vti_path: Path, title: str, out_path: Path, streamlines: bool) -> None:
    data = read_vti_velocity(vti_path)

    # ------------------------------------------------------------------
    # Divergence of the velocity field
    # ------------------------------------------------------------------

    div = (
        0.5 * (data.ux[2:, 1:-1] - data.ux[:-2, 1:-1]) +
        0.5 * (data.uy[1:-1, 2:] - data.uy[1:-1, :-2])
    )

    div_rms = np.sqrt(np.mean(div**2))
    div_max = np.max(np.abs(div))

    print(
        f"{vti_path.stem}: "
        f"RMS(div u) = {div_rms:.6e}   "
        f"max|div u| = {div_max:.6e}"
    )
    
    speed = data.speed
    x = np.arange(data.nx, dtype=np.float64)
    y = np.arange(data.ny, dtype=np.float64)

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    image = ax.imshow(
        speed,
        origin="lower",
        extent=(0, data.nx - 1, 0, data.ny - 1),
        cmap="viridis",
        aspect="equal",
    )

    if streamlines and data.nx >= 2 and data.ny >= 2:
        stride = max(1, max(data.nx, data.ny) // 160)
        ax.streamplot(
            x[::stride],
            y[::stride],
            data.ux[::stride, ::stride],
            data.uy[::stride, ::stride],
            color="white",
            density=1.0,
            linewidth=0.65,
            arrowsize=0.7,
        )

    colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    colorbar.set_label(r"$|\mathbf{u}|$")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    #ax.set_title(title)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_vti_velocity(
    ids: Path,
    runs_dir: Path,
    figures_dir: Path,
    which: str,
    streamlines: bool,
) -> None:
    apply_style()
    specs = read_run_specs(ids)
    out_dir = figures_dir / "velocity_fields"
    out_dir.mkdir(parents=True, exist_ok=True)
    loaded = 0

    for spec in specs:
        try:
            config = parse_run_config(spec.run_id, stdout_path(runs_dir, spec.run_id))
            files = select_vti_files(vtk_dir(runs_dir, spec.run_id), which)
        except Exception as exc:
            print(f"[warning] skipping VTI files for {spec.run_id}: {exc}")
            continue

        label = spec.label or default_label(config)
        for vti_path in files:
            data_step = read_vti_velocity(vti_path).step
            step_label = f"step {data_step}" if data_step is not None else vti_path.stem
            stem = safe_filename(f"{spec.run_id}_{step_label}")
            out_path = out_dir / f"velocity_{stem}.png"
            _plot_one(vti_path, f"{label} - {step_label}", out_path, streamlines)
            loaded += 1

    if loaded == 0:
        raise RuntimeError("no valid VTI velocity fields loaded")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot velocity fields from ASCII .vti files.")
    parser.add_argument("--ids", type=Path, default=DEFAULT_IDS)
    parser.add_argument("--runs-dir", type=Path, default=DEFAULT_RUNS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--which", choices=("latest", "first", "all"), default="latest")
    parser.add_argument("--streamlines", action="store_true")
    args = parser.parse_args()

    plot_vti_velocity(args.ids, args.runs_dir, args.figures_dir, args.which, args.streamlines)


if __name__ == "__main__":
    main()
