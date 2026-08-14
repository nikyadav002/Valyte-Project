#!/usr/bin/env python3
"""Combined Band Structure & DOS side-by-side plotting."""

import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator

from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

from valyte.dos_plot import load_dos


def gradient_fill_rotated(y, x, ax=None, color=None, **kwargs):
    """Fill area under a curve with a horizontal gradient (useful for rotated DOS).

    y: DOS values (plotted on X-axis)
    x: Energies (plotted on Y-axis)
    """
    if ax is None:
        ax = plt.gca()

    if len(x) == 0 or len(y) == 0:
        return None

    line, = ax.plot(y, x, color=color, lw=1.5, **kwargs)

    fill_color = line.get_color() if color is None else color
    alpha = line.get_alpha() or 1.0
    zorder = line.get_zorder()

    # Horizontal gradient of shape (1, 100, 4)
    z = np.empty((1, 100, 4))
    rgb = mcolors.to_rgb(fill_color)
    z[:, :, :3] = rgb

    min_alpha = 0.05
    max_alpha = 0.95
    gradient_vector = np.linspace(min_alpha, max_alpha, 100)
    gradient_vector *= alpha

    if np.mean(y) < 0:
        gradient_vector = gradient_vector[::-1]

    z[0, :, -1] = gradient_vector

    ymin, ymax = x.min(), x.max()
    local_xmax = max(y.max(), abs(y.min()))
    if local_xmax == 0:
        local_xmax = 1.0

    if np.mean(y) < 0:
        extent_xmin = -local_xmax
        extent_xmax = 0
    else:
        extent_xmin = 0
        extent_xmax = local_xmax

    im = ax.imshow(
        z,
        aspect="auto",
        extent=[extent_xmin, extent_xmax, ymin, ymax],
        origin="lower",
        zorder=zorder,
    )

    xy = np.column_stack([y, x])
    verts = np.vstack([[0, x[0]], xy, [0, x[-1]], [0, x[0]]])

    clip = Polygon(verts, lw=0, facecolor="none", closed=True)
    ax.add_patch(clip)
    im.set_clip_path(clip)

    return line


def plot_combined(
    vasprun_path,
    kpoints_path=None,
    dos_path=None,
    elements=None,
    plotting_config=None,
    output="valyte_combined.png",
    ylim=(-4, 4),
    figsize=(3.2, 3.2),
    dpi=400,
    font="Arial",
    show_fermi=False,
    scale_factor=1.0,
    legend_cutoff=0.10,
    save_data=False,
    spin_resolved=False,
    bold=True,
    colors=None,
):
    """Plot combined Band Structure and DOS side-by-side."""
    if os.path.isdir(vasprun_path):
        vasprun_path = os.path.join(vasprun_path, "vasprun.xml")

    if not os.path.exists(vasprun_path):
        raise FileNotFoundError(f"{vasprun_path} not found")

    # Font and styling
    plt.style.use("default")
    font_map = {
        "arial": "Arial",
        "helvetica": "Helvetica",
        "times": "Times New Roman",
        "times new roman": "Times New Roman",
    }
    _weight = "bold" if bold else "normal"
    font = font_map.get(font.lower(), "Arial")
    mpl.rcParams["font.family"] = font
    mpl.rcParams["axes.linewidth"] = 0.8
    mpl.rcParams["font.weight"] = _weight
    mpl.rcParams["font.size"] = 12
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

    # Load Band Structure
    try:
        vr = BSVasprun(vasprun_path, parse_projected_eigen=False)
        bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)
    except Exception as e:
        raise ValueError(f"Failed to load band structure: {e}")

    bs_plotter = BSPlotter(bs)
    band_data = bs_plotter.bs_plot_data(zero_to_efermi=True)
    distances = band_data["distances"]
    band_energies = band_data["energy"]
    ticks = band_data["ticks"]

    # Load DOS (it automatically zeroes relative to Fermi level)
    dos_load_path = dos_path if dos_path else os.path.dirname(vasprun_path)
    dos_data, pdos_data = load_dos(dos_load_path, elements)

    # Resolve items to plot in DOS (orbital resolved by default)
    if not plotting_config:
        plotting_config = []
        for el in sorted(pdos_data.keys()):
            for orb in sorted(pdos_data[el].keys()):
                plotting_config.append((el, orb))

    # Setup figure
    fig, (ax_band, ax_dos) = plt.subplots(
        1, 2, figsize=figsize, sharey=True,
        gridspec_kw={"width_ratios": [2.2, 1], "wspace": 0.06}
    )

    # 1. Plot Band Structure
    color_vb = "#8e44ad"
    color_cb = "#2a9d8f"

    if isinstance(band_energies, dict):
        spin_list = list(band_energies.keys())
    else:
        spin_list = list(band_energies[0].keys())

    if spin_resolved and len(spin_list) < 2:
        print("Warning: Only one spin channel found in this calculation; --spin-resolved has no effect.")
        spin_resolved = False

    spin_styles = {}
    if spin_resolved:
        spin_styles = {
            spin_list[0]: {"color": "#3498db", "ls": "-"},
            spin_list[1]: {"color": "#e74c3c", "ls": "--"},
        }

    for i in range(len(distances)):
        d = distances[i]
        if isinstance(band_energies, dict):
            for spin in band_energies:
                for band in band_energies[spin][i]:
                    if spin_resolved:
                        style = spin_styles[spin]
                        ax_band.plot(d, band, color=style["color"], ls=style["ls"], lw=1.5, alpha=1.0)
                    else:
                        c = color_vb if np.mean(band) <= 0 else color_cb
                        ax_band.plot(d, band, color=c, lw=1.5, alpha=1.0)
        else:
            for spin in band_energies[i]:
                for band in band_energies[i][spin]:
                    if spin_resolved:
                        style = spin_styles[spin]
                        ax_band.plot(d, band, color=style["color"], ls=style["ls"], lw=1.5, alpha=1.0)
                    else:
                        c = color_vb if np.mean(band) <= 0 else color_cb
                        ax_band.plot(d, band, color=c, lw=1.5, alpha=1.0)

    if spin_resolved:
        up_line = mlines.Line2D([], [], color="#3498db", lw=1.5, ls="-", label="Spin up")
        dn_line = mlines.Line2D([], [], color="#e74c3c", lw=1.5, ls="--", label="Spin down")
        ax_band.legend(handles=[up_line, dn_line], fontsize=10, frameon=False, loc="upper left")

    # Set band structure axes
    ax_band.set_xticks(ticks["distance"])
    clean_labels = [(l or "").replace("$\\mid$", "|") for l in ticks["label"]]
    if bold:
        clean_labels = [
            l.replace("\\Gamma", "\\mathbf{\\Gamma}")
             .replace("\\Sigma", "\\mathbf{\\Sigma}")
             .replace("\\Delta", "\\mathbf{\\Delta}")
             .replace("\\Lambda", "\\mathbf{\\Lambda}")
            for l in clean_labels
        ]
    ax_band.set_xticklabels(clean_labels, fontsize=12, fontweight=_weight)

    for d in ticks["distance"]:
        ax_band.axvline(d, color="k", lw=0.8, ls="-", alpha=0.3)

    if show_fermi:
        ax_band.axhline(0, color="k", lw=0.8, ls="--", alpha=0.7)

    ax_band.set_ylabel("Energy (eV)", fontsize=14, weight=_weight, labelpad=6)
    ax_band.set_ylim(ylim)
    ax_band.set_xlim(distances[0][0], distances[-1][-1])

    # Enforce bottom ticks only on the band plot
    ax_band.tick_params(axis="x", direction="in", which="both", bottom=True, top=False)
    ax_band.tick_params(axis="y", direction="in", which="both", left=True, right=False)

    # 2. Plot DOS
    is_spin_polarized = Spin.down in dos_data.densities

    palette = [
        "#e63946", "#457b9d", "#2a9d8f", "#f4a261", "#6a4c93",
        "#8ac926", "#1982c4", "#ca6702", "#ff595e", "#6a994e",
        "#b5179e", "#219ebc", "#9b2226", "#606c38", "#0077b6",
        "#bb3e03", "#005f73", "#ee9b00", "#7209b7", "#94d2bd",
    ]

    dos_x_mask = (dos_data.energies >= ylim[0]) & (dos_data.energies <= ylim[1])

    max_visible_dos = 0
    min_visible_dos = 0
    dos_lines = []

    for i, (el, orb) in enumerate(plotting_config):
        if el not in pdos_data:
            continue

        label = el if orb == "total" else f"{el}({orb})"
        c = resolve_color(label, el, orb, i, colors, palette)


        if orb == "total":
            y_up = np.zeros_like(dos_data.energies)
            y_down = np.zeros_like(dos_data.energies)
            for o_data in pdos_data[el].values():
                y_up += o_data.get(Spin.up, np.zeros_like(dos_data.energies))
                y_down += o_data.get(Spin.down, np.zeros_like(dos_data.energies))
            label = el
        else:
            if orb not in pdos_data[el]:
                continue
            y_up = pdos_data[el][orb].get(Spin.up, np.zeros_like(dos_data.energies))
            y_down = pdos_data[el][orb].get(Spin.down, np.zeros_like(dos_data.energies))
            label = f"{el}({orb})"

        y_down = -y_down

        visible_y_up = y_up[dos_x_mask]
        visible_y_down = y_down[dos_x_mask]

        current_max_dos = 0
        if len(visible_y_up) > 0:
            max_y = np.max(visible_y_up)
            max_visible_dos = max(max_visible_dos, max_y)
            current_max_dos = max(current_max_dos, max_y)

        if is_spin_polarized and len(visible_y_down) > 0:
            min_y = np.min(visible_y_down)
            min_visible_dos = min(min_visible_dos, min_y)
            current_max_dos = max(current_max_dos, abs(min_y))

        # We plot with DOS values on X-axis, Energies on Y-axis
        line = ax_dos.plot(y_up, dos_data.energies, lw=1.5, color=c, alpha=0)[0]
        gradient_fill_rotated(y_up, dos_data.energies, ax=ax_dos, color=c, alpha=0.9)
        if is_spin_polarized:
            gradient_fill_rotated(y_down, dos_data.energies, ax=ax_dos, color=c, alpha=0.9)

        dos_lines.append(
            {
                "line": line,
                "label": label,
                "max_dos": current_max_dos,
            }
        )

    # Set limits and ticks for DOS plot
    global_max_dos = max(max_visible_dos, abs(min_visible_dos))
    if global_max_dos == 0:
        global_max_dos = 1.0

    dos_xlim_max = (global_max_dos * 1.1) / scale_factor
    dos_xlim_min = (min_visible_dos * 1.1) / scale_factor if is_spin_polarized else 0
    ax_dos.set_xlim(dos_xlim_min, dos_xlim_max)

    if show_fermi:
        ax_dos.axhline(0, color="k", lw=0.8, ls="--", alpha=0.7)

    # Ticks only at the bottom, nowhere else, and no numbers (labels)!
    ax_dos.tick_params(axis="x", direction="in", which="both", bottom=True, top=False, labelbottom=False)
    ax_dos.tick_params(axis="y", which="both", left=False, right=False)

    # Set x-label (removed DOS text as requested)
    ax_dos.set_xlabel("")

    # Render Legend
    threshold = legend_cutoff * global_max_dos
    legend_handles = []
    legend_labels = []

    for item in dos_lines:
        if item["max_dos"] >= threshold:
            item["line"].set_alpha(1.0)
            legend_handles.append(item["line"])
            legend_labels.append(item["label"])

    if legend_handles:
        legend = ax_dos.legend(
            legend_handles,
            legend_labels,
            frameon=False,
            fontsize=8.0,
            loc="upper right",
            ncol=1,
            handlelength=1.0,
            columnspacing=0.5,
            handletextpad=0.4,
            borderaxespad=0.2,
        )
        for text in legend.get_texts():
            text.set_fontweight(_weight)

    # Minor x-ticks on both panels
    ax_band.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax_dos.xaxis.set_minor_locator(AutoMinorLocator(2))

    # Force formatting y-tick numbers to bold
    plt.setp(ax_band.get_yticklabels(), fontweight=_weight)

    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")

    if save_data:
        # Save both band and DOS data
        # Band data
        band_out = output.replace(".png", "_band.dat")
        nbranches = len(distances)
        all_d = np.concatenate(distances)
        if isinstance(band_energies, dict):
            spins = list(band_energies.keys())
            nbands = len(band_energies[spins[0]][0])
            all_e = {
                sp: np.array([
                    np.concatenate([np.array(band_energies[sp][i][ib]) for i in range(nbranches)])
                    for ib in range(nbands)
                ]) for sp in spins
            }
        else:
            spins = list(band_energies[0].keys())
            nbands = len(band_energies[0][spins[0]])
            all_e = {
                sp: np.array([
                    np.concatenate([np.array(band_energies[i][sp][ib]) for i in range(nbranches)])
                    for ib in range(nbands)
                ]) for sp in spins
            }

        label_info = ", ".join(
            f"{l}={d:.4f}" for l, d in zip(ticks["label"], ticks["distance"]) if l
        )
        spin_tag = {spins[0]: ""} if len(spins) == 1 else {spins[0]: "_up", spins[1]: "_dn"}
        col_labels = ["k-dist"] + [f"band_{ib+1}{spin_tag[sp]}" for sp in spins for ib in range(nbands)]
        header = f"K-point labels: {label_info}\n" + "  ".join(col_labels)
        cols = [all_d] + [all_e[sp][ib] for sp in spins for ib in range(nbands)]
        np.savetxt(band_out, np.column_stack(cols), header=header, fmt="%.6f")

        # DOS data (only PDOS as specified)
        dos_out = output.replace(".png", "_dos.dat")
        dos_cols = [dos_data.energies]
        dos_labels = ["Energy(eV)"]
        for el, orb in plotting_config:
            if el not in pdos_data:
                continue
            if orb == "total":
                y_up = sum(o.get(Spin.up, np.zeros_like(dos_data.energies)) for o in pdos_data[el].values())
                label = el
            else:
                if orb not in pdos_data[el]:
                    continue
                y_up = pdos_data[el][orb].get(Spin.up, np.zeros_like(dos_data.energies))
                label = f"{el}({orb})"

            dos_cols.append(y_up)
            dos_labels.append(f"{label}_up")

            if is_spin_polarized:
                if orb == "total":
                    y_dn = sum(o.get(Spin.down, np.zeros_like(dos_data.energies)) for o in pdos_data[el].values())
                else:
                    y_dn = pdos_data[el][orb].get(Spin.down, np.zeros_like(dos_data.energies))
                dos_cols.append(-y_dn)
                dos_labels.append(f"{label}_dn")

        header = "  ".join(dos_labels)
        np.savetxt(dos_out, np.column_stack(dos_cols), header=header, fmt="%.6f")
        print(f"Saved band data: {band_out}")
        print(f"Saved DOS data: {dos_out}")
