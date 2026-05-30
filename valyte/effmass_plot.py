"""Effective mass parabolic fit plotting."""

import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt


def plot_effective_mass(results, output="valyte_effmass.png",
                        figsize=(8, 4), dpi=400, font="Arial"):
    """Plot parabolic fits overlaid on band data near VBM and CBM.

    Parameters
    ----------
    results : dict
        Output from compute_effective_masses().
    output : str
        Output image filename.
    figsize : tuple
        Figure size in inches.
    dpi : int
        Output resolution.
    font : str
        Font family name.
    """
    if results.get("is_metal", False):
        print("No plot generated — system is metallic.")
        return

    hole_masses = results.get("hole_masses", [])
    electron_masses = results.get("electron_masses", [])

    if not hole_masses and not electron_masses:
        print("No effective mass fits to plot.")
        return

    # ── Matplotlib style (matching Valyte conventions) ──────────────
    font_map = {
        "arial": "Arial",
        "helvetica": "Helvetica",
        "times": "Times New Roman",
        "times new roman": "Times New Roman",
    }
    font = font_map.get(font.lower(), "Arial")
    mpl.rcParams["font.family"] = font
    mpl.rcParams["axes.linewidth"] = 1.4
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["font.size"] = 14
    mpl.rcParams["xtick.major.width"] = 1.2
    mpl.rcParams["ytick.major.width"] = 1.2

    has_holes = len(hole_masses) > 0
    has_electrons = len(electron_masses) > 0

    if has_holes and has_electrons:
        fig, (ax_vbm, ax_cbm) = plt.subplots(1, 2, figsize=figsize)
    elif has_holes:
        fig, ax_vbm = plt.subplots(1, 1, figsize=(figsize[0] / 2, figsize[1]))
        ax_cbm = None
    else:
        fig, ax_cbm = plt.subplots(1, 1, figsize=(figsize[0] / 2, figsize[1]))
        ax_vbm = None

    # ── Colors & line styles ─────────────────────────────────────────
    color_vbm = "#8e44ad"  # purple (valence band color from band_plot.py)
    color_cbm = "#2a9d8f"  # teal (conduction band color from band_plot.py)

    # Use distinct line styles for different directions
    linestyles = ["-", "--", "-.", ":"]
    markers = ["o", "s", "D", "^", "v", "<", ">"]

    # ── VBM panel ────────────────────────────────────────────────────
    if ax_vbm is not None and hole_masses:
        _plot_panel(ax_vbm, hole_masses, color_vbm, linestyles, markers,
                    "Hole (VBM)")

    # ── CBM panel ────────────────────────────────────────────────────
    if ax_cbm is not None and electron_masses:
        _plot_panel(ax_cbm, electron_masses, color_cbm, linestyles, markers,
                    "Electron (CBM)")

    plt.tight_layout(pad=1.5)
    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")


def _plot_panel(ax, masses, base_color, linestyles, markers, title):
    """Plot a single panel (VBM or CBM) with fit curves."""
    # Slightly vary the color for different directions
    from matplotlib.colors import to_rgb
    base_rgb = np.array(to_rgb(base_color))

    for i, m in enumerate(masses):
        k_data = m["k_data"]
        e_data = m["e_data"]
        k_fit = m["k_fit"]
        e_fit = m["e_fit"]
        direction = m["direction"]
        m_star = m["m_star"]
        r_sq = m["r_squared"]

        # Slightly shift color for each direction
        factor = 1.0 - 0.15 * (i % 4)
        color = np.clip(base_rgb * factor, 0, 1)
        color_hex = "#{:02x}{:02x}{:02x}".format(
            int(color[0] * 255), int(color[1] * 255), int(color[2] * 255)
        )

        ls = linestyles[i % len(linestyles)]
        mk = markers[i % len(markers)]

        # Data points
        ax.plot(k_data, e_data, mk, color=color_hex, markersize=7,
                markeredgewidth=1.2, markeredgecolor="white", zorder=5)

        # Fit curve
        ax.plot(k_fit, e_fit, ls=ls, color=color_hex, lw=2.0, alpha=0.8,
                zorder=4, label=direction)

        # Annotate with m* value
        # Place annotation at the curve edge
        x_ann = k_fit[int(len(k_fit) * 0.75)]
        y_ann = e_fit[int(len(e_fit) * 0.75)]

        if abs(m_star) < 1e6:
            ann_text = f"m* = {m_star:.3f} $m_0$"
        else:
            ann_text = r"m* = $\infty$"

        ax.annotate(
            ann_text,
            xy=(x_ann, y_ann),
            fontsize=9,
            fontweight="bold",
            color=color_hex,
            xytext=(5, 8 + 12 * i),
            textcoords="offset points",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=color_hex, alpha=0.85, linewidth=0.8),
        )

    ax.set_xlabel("k (1/Å)", fontsize=14, fontweight="bold", labelpad=6)
    ax.set_ylabel("Energy (eV)", fontsize=14, fontweight="bold", labelpad=6)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=8)

    ax.axhline(0, color="k", lw=0.6, ls="--", alpha=0.4)
    ax.axvline(0, color="k", lw=0.6, ls="--", alpha=0.4)

    if len(masses) > 1:
        ax.legend(fontsize=9, frameon=True, loc="best", framealpha=0.9,
                  edgecolor="0.8")
