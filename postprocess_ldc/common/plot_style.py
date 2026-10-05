from __future__ import annotations

import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lbm-matplotlib")

import matplotlib.pyplot as plt


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "font.size": 14,
            "axes.labelsize": 16,
            "axes.titlesize": 15,
            "legend.fontsize": 11,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 1.1,
            "axes.grid": True,
            "grid.alpha": 0.18,
            "grid.linewidth": 0.7,
            "lines.linewidth": 2.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "legend.frameon": False,
        }
    )
