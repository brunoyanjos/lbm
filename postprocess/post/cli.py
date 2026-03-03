# post/cli.py
from __future__ import annotations

import argparse
from pathlib import Path

from .io import read_tke_bin, read_centerline_bin
from .io.read_stdout import read_stdout_summary
from .plots import plot_tke, plot_centerlines
from .plots.annul_profile import plot_annul_profile
from .paths import (
    default_runs_dir,
    default_outcomes_dir,
    resolve_run_outputs,
    resolve_out_dir,
)

from .io.read_bin import read_annul_profile_bin, read_annul_pressure_bin
from .plots.annul_pressure import plot_annul_pressure


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

    # ---- Annular profile command ----
    p_ap = sub.add_parser(
        "annul-profile",
        help="Plot annular velocity profile (right side line at y=NY/2)",
    )
    p_ap.add_argument("run_id")
    p_ap.add_argument("--show", action="store_true")
    p_ap.add_argument(
        "--no-ur", action="store_true", help="Do not plot ur (sanity check)"
    )

    # ---- Annular pressure command ----
    p_apr = sub.add_parser(
        "annul-pressure",
        help="Plot annular pressure profile (Δp* vs r*)",
    )
    p_apr.add_argument("run_id")
    p_apr.add_argument("--show", action="store_true")

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
    u_lid = sim.get("U_lid", None)  # no Taylor–Couette usamos como U_wall (fallback)
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

        # Título no mesmo estilo
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

    elif args.command == "annul-profile":
        ap_path = run_path / "annul_profile.bin"
        if not ap_path.exists():
            raise FileNotFoundError(f"Annular profile file not found: {ap_path}")

        prof_meta, r, ur, ut = read_annul_profile_bin(ap_path)
        save_path = save_dir / "annul_profile.png"

        title_parts = [f"{stencil}", "Annular profile"]
        if Re is not None:
            title_parts.append(f"Re={Re:g}")
        if mlups is not None:
            title_parts.append(f"MLUPS={mlups:.1f}")
        title = " | ".join(title_parts)

        plot_annul_profile(
            sim,
            prof_meta,
            r,
            ur,
            ut,
            title=title,
            show=args.show,
            savepath=save_path,
            show_ur=(not args.no_ur),
        )

        print(f"Saved figure to: {save_path}")

    elif args.command == "annul-pressure":
        ap_path = run_path / "annul_pressure.bin"
        if not ap_path.exists():
            raise FileNotFoundError(f"Annular pressure file not found: {ap_path}")

        prof_meta, r, rho_p = read_annul_pressure_bin(ap_path)
        save_path = save_dir / "annul_pressure.png"

        title_parts = [f"{stencil}", "Annular pressure"]
        if Re is not None:
            title_parts.append(f"Re={Re:g}")
        if mlups is not None:
            title_parts.append(f"MLUPS={mlups:.1f}")
        title = " | ".join(title_parts)

        plot_annul_pressure(
            sim,
            prof_meta,
            r,
            rho_p,
            title=title,
            show=args.show,
            savepath=save_path,
            show_analytic=True,
        )
        print(f"Saved figure to: {save_path}")
