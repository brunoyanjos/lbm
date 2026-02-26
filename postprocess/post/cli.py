# post/cli.py
from __future__ import annotations

import argparse
from pathlib import Path

from .io import read_tke_bin, read_centerline_bin
from .io.read_stdout import read_stdout_summary
from .plots import plot_tke, plot_centerlines
from .paths import (
    default_runs_dir,
    default_outcomes_dir,
    resolve_run_outputs,
    resolve_out_dir,
)


def main():
    parser = argparse.ArgumentParser(description="LBM postprocessing tools")

    parser.add_argument("--runs-dir", default=str(default_runs_dir()))
    parser.add_argument("--out-dir", default=str(default_outcomes_dir()))

    sub = parser.add_subparsers(dest="command", required=True)

    # ---- TKE command ----
    p_tke = sub.add_parser("tke", help="Plot kinetic energy")
    p_tke.add_argument("run_id")
    p_tke.add_argument("--show", action="store_true")

    # ---- Centerline command ----
    p_cl = sub.add_parser("centerline", help="Plot centerline velocity profiles")
    p_cl.add_argument("run_id")
    p_cl.add_argument("--show", action="store_true")

    args = parser.parse_args()

    runs_dir = Path(args.runs_dir)
    outcomes_dir = Path(args.out_dir)

    run_root = runs_dir / args.run_id
    run_path = resolve_run_outputs(runs_dir, args.run_id)

    save_dir = resolve_out_dir(outcomes_dir, args.run_id)
    save_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load stdout summary (common metadata) ----
    stdout_path = run_root / "logs" / "stdout.txt"
    if not stdout_path.exists():
        raise FileNotFoundError(f"stdout.txt not found: {stdout_path}")

    sim = read_stdout_summary(stdout_path)

    stencil = sim.get("Stencil", "UNKNOWN")
    Re = sim.get("Re", None)
    u_lid = sim.get("U_lid", None)
    mlups = sim.get("MLUPS_GPU", None)

    # fallback seguro caso stdout não tenha (não deveria, mas evita crash)
    if u_lid is None:
        u_lid = 1.0

    # ---- Command dispatch ----
    if args.command == "tke":
        tke_path = run_path / "tke.bin"
        if not tke_path.exists():
            raise FileNotFoundError(f"TKE file not found: {tke_path}")

        t_star, ke = read_tke_bin(tke_path)
        save_path = save_dir / "tke.png"

        # Título "perfeito" no mesmo estilo
        tke_title_parts = [f"{stencil}"]
        if Re is not None:
            tke_title_parts.append(f"Re={Re:g}")
        if mlups is not None:
            tke_title_parts.append(f"MLUPS={mlups:.1f}")
        tke_title = " | ".join(tke_title_parts)

        plot_tke(t_star, ke, title=tke_title, show=args.show, savepath=save_path)
        print(f"Saved figure to: {save_path}")

    elif args.command == "centerline":
        cl_path = run_path / "centerline.bin"
        if not cl_path.exists():
            raise FileNotFoundError(f"Centerline file not found: {cl_path}")

        meta, ux_xc_y, uy_yc_x = read_centerline_bin(cl_path)
        save_path = save_dir / "centerline.png"

        # Título "perfeito"
        title_parts = [f"{stencil}"]
        if Re is not None:
            title_parts.append(f"Re={Re:g}")
        if mlups is not None:
            title_parts.append(f"MLUPS={mlups:.1f}")
        title = " | ".join(title_parts)

        plot_centerlines(
            meta,
            ux_xc_y,
            uy_yc_x,
            u_lid=u_lid,
            Re=Re,
            title=title,
            show=args.show,
            savepath=save_path,
        )

        print(f"Saved figure to: {save_path}")
