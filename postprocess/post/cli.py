# post/cli.py
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .paths import (
    default_runs_dir,
    default_outcomes_dir,
    resolve_run_outputs,
    resolve_out_dir,
)

from .io import read_tke_bin, read_centerline_bin
from .io.read_stdout import read_stdout_summary
from .io.read_bin import (
    read_annul_profile_bin,
    read_annul_pressure_bin,
    read_channel_profile_bin,
    read_couette_profile_bin,
    read_jet_centerline_bin,
    read_jet_sections_bin,
)

from .plots import (
    plot_tke,
    plot_centerlines,
    plot_channel_profile,
    plot_couette_profile,
    plot_jet_centerline,
    plot_jet_sections,
)

from .plots.annul_profile import plot_annul_profile
from .plots.annul_pressure import plot_annul_pressure


def _require_file(p: Path, what: str) -> Path:
    if not p.exists():
        raise FileNotFoundError(f"{what} not found: {p}")
    return p


def _make_title(stencil: str, label: str | None = None, sim: dict | None = None) -> str:
    parts: list[str] = [stencil]
    if label:
        parts.append(label)

    if sim is not None:
        nx = sim.get("nx", None)
        ny = sim.get("ny", None)
        if nx is not None and ny is not None:
            parts.append(f"{int(nx)}×{int(ny)}")

        Re = sim.get("Re", None)
        if Re is not None:
            parts.append(f"Re={float(Re):g}")

        mlups = sim.get("MLUPS_GPU", None)
        if mlups is not None:
            parts.append(f"MLUPS={float(mlups):.1f}")

    return " | ".join(parts)


def _warn(msg: str) -> None:
    print(f"[lbm-post] {msg}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description="LBM postprocessing tools")

    parser.add_argument("--runs-dir", type=Path, default=default_runs_dir())
    parser.add_argument("--out-dir", type=Path, default=default_outcomes_dir())

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
    p_ap.add_argument("--no-ur", action="store_true", help="Do not plot ur")

    # ---- Annular pressure command ----
    p_apr = sub.add_parser(
        "annul-pressure",
        help="Plot annular pressure profile (Δp* vs r*)",
    )
    p_apr.add_argument("run_id")
    p_apr.add_argument("--show", action="store_true")

    # ---- Channel profile ----
    p_ch = sub.add_parser("channel-profile", help="Plot channel ux(y) profile")
    p_ch.add_argument("run_id")
    p_ch.add_argument("--show", action="store_true")

    # ---- Couette profile ----
    p_cu = sub.add_parser("couette-profile", help="Plot couette ux(y) profile")
    p_cu.add_argument("run_id")
    p_cu.add_argument("--show", action="store_true")

    # ---- Jet ----
    p_jcl = sub.add_parser("jet-centerline", help="Plot jet centerline ux(x)")
    p_jcl.add_argument("run_id")
    p_jcl.add_argument("--show", action="store_true")

    p_jsc = sub.add_parser(
        "jet-sections", help="Plot jet sections ux(y) at x positions"
    )
    p_jsc.add_argument("run_id")
    p_jsc.add_argument("--show", action="store_true")

    args = parser.parse_args()

    runs_dir: Path = args.runs_dir
    outcomes_dir: Path = args.out_dir

    run_path = resolve_run_outputs(runs_dir, args.run_id)  # runs/<run_id>/outputs
    run_root = run_path.parent  # runs/<run_id>

    save_dir = resolve_out_dir(outcomes_dir, args.run_id)
    save_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load stdout summary (common metadata) ----
    stdout_path = run_root / "logs" / "stdout.txt"
    _require_file(stdout_path, "stdout.txt")

    sim = read_stdout_summary(stdout_path)

    stencil = str(sim.get("Stencil", "UNKNOWN"))

    # em alguns casos (cavidade) pode existir como U_lid; em outros casos pode não existir
    u_lid = sim.get("U_lid", None)
    if u_lid is None or float(u_lid) == 0.0:
        _warn("stdout missing/invalid U_lid; using fallback u_lid=1.0")
        u_lid = 1.0
    else:
        u_lid = float(u_lid)

    # ---- Command dispatch ----
    if args.command == "tke":
        tke_path = _require_file(run_path / "tke.bin", "TKE file")
        t_star, ke = read_tke_bin(tke_path)

        save_path = save_dir / "tke.png"
        title = _make_title(stencil, None, sim)

        plot_tke(t_star, ke, title=title, show=args.show, savepath=save_path)
        print(f"Saved figure to: {save_path}")

    elif args.command == "centerline":
        cl_path = _require_file(run_path / "centerline.bin", "Centerline file")
        meta, ux_xc_y, uy_yc_x = read_centerline_bin(cl_path)

        save_path = save_dir / "centerline.png"
        title = _make_title(stencil, None, sim)

        plot_centerlines(
            meta,
            ux_xc_y,
            uy_yc_x,
            u_lid=u_lid,
            Re=sim.get("Re", None),
            title=title,
            show=args.show,
            savepath=save_path,
        )
        print(f"Saved figure to: {save_path}")

    elif args.command == "annul-profile":
        ap_path = _require_file(run_path / "annul_profile.bin", "Annular profile file")
        prof_meta, r, ur, ut = read_annul_profile_bin(ap_path)

        save_path = save_dir / "annul_profile.png"
        title = _make_title(stencil, None, sim)

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
        ap_path = _require_file(
            run_path / "annul_pressure.bin", "Annular pressure file"
        )
        prof_meta, r, rho_p = read_annul_pressure_bin(ap_path)

        save_path = save_dir / "annul_pressure.png"
        title = _make_title(stencil, None, sim)

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

    elif args.command == "channel-profile":
        p = _require_file(run_path / "channel_profile.bin", "Channel profile")
        meta, ux_y = read_channel_profile_bin(p)

        save_path = save_dir / "channel_profile.png"
        title = _make_title(stencil, None, sim)

        plot_channel_profile(
            sim, meta, ux_y, title=title, show=args.show, savepath=save_path
        )
        print(f"Saved figure to: {save_path}")

    elif args.command == "couette-profile":
        p = _require_file(run_path / "couette_profile.bin", "Couette profile")
        meta, ux_y = read_couette_profile_bin(p)

        save_path = save_dir / "couette_profile.png"
        title = _make_title(stencil, None, sim)

        plot_couette_profile(
            sim, meta, ux_y, title=title, show=args.show, savepath=save_path
        )
        print(f"Saved figure to: {save_path}")

    elif args.command == "jet-centerline":
        p = _require_file(run_path / "jet_centerline.bin", "Jet centerline")
        meta, ux_x = read_jet_centerline_bin(p)

        save_path = save_dir / "jet_centerline.png"
        title = _make_title(stencil, None, sim)

        plot_jet_centerline(
            sim, meta, ux_x, title=title, show=args.show, savepath=save_path
        )
        print(f"Saved figure to: {save_path}")

    elif args.command == "jet-sections":
        p = _require_file(run_path / "jet_sections.bin", "Jet sections")
        meta, x_sec, ux_sec = read_jet_sections_bin(p)

        save_path = save_dir / "jet_sections.png"
        title = _make_title(stencil, None, sim)

        plot_jet_sections(
            sim, meta, x_sec, ux_sec, title=title, show=args.show, savepath=save_path
        )
        print(f"Saved figure to: {save_path}")

    else:
        raise ValueError(f"Unknown command: {args.command}")
