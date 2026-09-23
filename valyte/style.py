"""Shared Matplotlib plotting style and helper utilities for Valyte."""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt

FONT_MAP = {
    "arial": "Arial",
    "helvetica": "Helvetica",
    "times": "Times New Roman",
    "times new roman": "Times New Roman",
}

DEFAULT_PALETTE = [
    "#e63946", "#457b9d", "#2a9d8f", "#f4a261", "#6a4c93",
    "#8ac926", "#1982c4", "#ca6702", "#ff595e", "#6a994e",
    "#b5179e", "#219ebc", "#9b2226", "#606c38", "#0077b6",
    "#bb3e03", "#005f73", "#ee9b00", "#7209b7", "#94d2bd",
]


def resolve_font_family(font: str = "Arial") -> str:
    """Resolve font family name from supported map."""
    if not font:
        return "Arial"
    return FONT_MAP.get(font.lower(), "Arial")


def get_font_weight(bold: bool = True) -> str:
    """Return 'bold' or 'normal' weight string."""
    return "bold" if bold else "normal"


def apply_style(
    font: str = "Arial",
    bold: bool = True,
    fontsize: float = 12,
    linewidth: float = None,
):
    """Apply consistent Matplotlib rcParams across all Valyte plotting modules."""
    plt.style.use("default")
    mpl.use("agg")
    mpl.rcParams["axes.unicode_minus"] = False

    weight = get_font_weight(bold)
    font_family = resolve_font_family(font)
    axes_lw = linewidth if linewidth is not None else (0.8 if not bold else 0.8)

    mpl.rcParams["font.family"] = font_family
    mpl.rcParams["axes.linewidth"] = axes_lw
    mpl.rcParams["font.weight"] = weight
    mpl.rcParams["font.size"] = fontsize
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"
    mpl.rcParams["xtick.major.width"] = 1.2 if bold else 0.8
    mpl.rcParams["ytick.major.width"] = 1.2 if bold else 0.8
    mpl.rcParams["xtick.minor.width"] = 0.8 if bold else 0.6
    mpl.rcParams["ytick.minor.width"] = 0.8 if bold else 0.6
    mpl.rcParams["xtick.major.size"] = 4 if bold else 5
    mpl.rcParams["ytick.major.size"] = 4 if bold else 5
    mpl.rcParams["xtick.minor.size"] = 2 if bold else 3
    mpl.rcParams["ytick.minor.size"] = 2 if bold else 3


def save_plot(fig, output: str, dpi: int = 400):
    """Save matplotlib figure with specified DPI and tight layout, then close."""
    plt.tight_layout()
    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
