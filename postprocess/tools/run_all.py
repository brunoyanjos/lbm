#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


COMMANDS = [
    ("tke", "tke.bin"),
    ("centerline", "centerline.bin"),
    ("annul-profile", "annul_profile.bin"),
    ("annul-pressure", "annul_pressure.bin"),
    ("channel-profile", "channel_profile.bin"),
    ("couette-profile", "couette_profile.bin"),
    ("jet-centerline", "jet_centerline.bin"),
    ("jet-sections", "jet_sections.bin"),
]


def run(cmd: list[str]) -> int:
    p = subprocess.run(cmd)
    return p.returncode


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Run all available lbm-post plots for one or more RUN_IDs"
    )
    ap.add_argument("run_ids", nargs="*", help="RUN_IDs (or use --file)")
    ap.add_argument("--file", type=Path, help="Text file with RUN_IDs (one per line)")
    ap.add_argument("--runs-dir", type=Path, default=Path("../runs"))
    ap.add_argument("--out-dir", type=Path, default=Path("outcomes"))
    ap.add_argument(
        "--no-show", action="store_true", help="Do not show plots (default)"
    )
    ap.add_argument("--stop-on-error", action="store_true")
    args = ap.parse_args()

    run_ids: list[str] = []

    if args.file:
        txt = args.file.read_text(encoding="utf-8", errors="replace").splitlines()
        run_ids.extend(
            [l.strip() for l in txt if l.strip() and not l.strip().startswith("#")]
        )

    run_ids.extend(args.run_ids)

    if not run_ids:
        print("No RUN_IDs provided.", file=sys.stderr)
        sys.exit(2)

    for rid in run_ids:
        out_path = args.runs_dir / rid / "outputs"
        print(f"\n=== {rid} ===")
        if not out_path.exists():
            print(f"[skip] outputs dir not found: {out_path}", file=sys.stderr)
            continue

        for cmd_name, bin_name in COMMANDS:
            p = out_path / bin_name
            if not p.exists():
                # não é erro: só não tem esse output
                continue

            cmd = [
                "lbm-post",
                "--runs-dir",
                str(args.runs_dir),
                "--out-dir",
                str(args.out_dir),
                cmd_name,
                rid,
            ]
            if not args.no_show:
                cmd.append("--show")

            print(f"[run] {cmd_name} ({bin_name})")
            rc = run(cmd)
            if rc != 0:
                print(f"[error] {cmd_name} failed for {rid} (rc={rc})", file=sys.stderr)
                if args.stop_on_error:
                    sys.exit(rc)

    print("\nDone.")


if __name__ == "__main__":
    main()
