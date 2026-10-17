#!/usr/bin/env python3
"""Combined Band Structure & DOS side-by-side plotting."""

import os
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator

from pymatgen.io.vasp import BSVasprun, Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

from valyte.dos_plot import load_dos, _collect_series, warn_missing_selection
from valyte.colors import assign_styles
from valyte.style import apply_style, get_font_weight, save_plot, font_scale


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
    fontsize=None,
    palette=None,
):
    """Plot combined Band Structure and DOS side-by-side.

    DOS colours follow the same rules as `valyte.dos_plot.plot_dos`: `palette`
    picks the palette ('wong', 'bright', 'muted' or a list of colours) and
    `colors` maps individual series.
    """
    if os.path.isdir(vasprun_path):
        vasprun_path = os.path.join(vasprun_path, "vasprun.xml")

    if not os.path.exists(vasprun_path):
        raise FileNotFoundError(f"{vasprun_path} not found")

    # Font and styling
    _weight = get_font_weight(bold)
    _fsbase = fontsize if fontsize is not None else 12
    _fscale = font_scale(fontsize, 12)
    apply_style(font=font, bold=bold, fontsize=_fsbase)

    # When the DOS comes from the same vasprun.xml (the default), parse the
    # file once and serve both panels from it.  A separate --dos source still
    # needs its own parse, so the lighter BSVasprun is enough for the bands.
    dos_same_source = dos_path is None
    try:
        if dos_same_source:
            vr = Vasprun(
                vasprun_path,
                parse_dos=True,
                parse_eigen=True,
                parse_projected_eigen=False,
                parse_potcar_file=False,
            )
        else:
            vr = BSVasprun(vasprun_path, parse_projected_eigen=False)
        bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)
    except Exception as e:
        raise ValueError(f"Failed to load band structure: {e}")

    bs_plotter = BSPlotter(bs)
    band_data = bs_plotter.bs_plot_data(zero_to_efermi=True)
    distances = band_data["distances"]
    band_energies = band_data["energy"]
    ticks = band_data["ticks"]

    # Both panels share a y-axis, so they must share one energy zero.
    # BSPlotter(zero_to_efermi=True) uses the VBM for a gapped system and the
    # Fermi level for a metal; hand that same value to the DOS loader instead
    # of letting it resolve a reference of its own.
    try:
        zero_energy = bs.efermi if bs.is_metal() else bs.get_vbm()["energy"]
    except Exception:
        zero_energy = None

    if dos_same_source:
        dos_data, pdos_data = load_dos(vasprun_path, elements,
                                       efermi=zero_energy, vr=vr)
    else:
        dos_data, pdos_data = load_dos(dos_path, elements, efermi=zero_energy)

    # Resolve items to plot in DOS (orbital resolved by default)
    user_selected = bool(plotting_config)
    if user_selected:
        warn_missing_selection(pdos_data, plotting_config)
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
        ax_band.legend(handles=[up_line, dn_line], fontsize=10 * _fscale, frameon=False, loc="upper left")

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
    ax_band.set_xticklabels(clean_labels, fontsize=12 * _fscale, fontweight=_weight)

    for d in ticks["distance"]:
        ax_band.axvline(d, color="k", lw=0.8, ls="-", alpha=0.3)

    if show_fermi:
        ax_band.axhline(0, color="k", lw=0.8, ls="--", alpha=0.7)

    ax_band.set_ylabel("Energy (eV)", fontsize=14 * _fscale, weight=_weight, labelpad=6)
    ax_band.set_ylim(ylim)
    ax_band.set_xlim(distances[0][0], distances[-1][-1])

    # Enforce bottom ticks only on the band plot
    ax_band.tick_params(axis="x", direction="in", which="both", bottom=True, top=False)
    ax_band.tick_params(axis="y", direction="in", which="both", left=True, right=False)

    # 2. Plot DOS
    is_spin_polarized = Spin.down in dos_data.densities

    dos_x_mask = (dos_data.energies >= ylim[0]) & (dos_data.energies <= ylim[1])

    series = _collect_series(dos_data, pdos_data, plotting_config, dos_x_mask,
                             is_spin_polarized)
    styles = assign_styles(series, palette=palette, colors=colors,
                           cutoff=legend_cutoff, explicit=user_selected)

    max_visible_dos = max((s["y_up"][dos_x_mask].max(initial=0.0) for s in series), default=0.0)
    min_visible_dos = 0
    if is_spin_polarized:
        min_visible_dos = min((s["y_down"][dos_x_mask].min(initial=0.0) for s in series),
                              default=0.0)

    # Grey background series first, then legend series largest-first.
    order = sorted(range(len(series)),
                   key=lambda i: (styles[i]["rank"] is not None, styles[i]["rank"] or 0))
    legend_items = []
    for i in order:
        s, st = series[i], styles[i]
        # DOS values on the x-axis, energies on the y-axis
        line = gradient_fill_rotated(s["y_up"], dos_data.energies, ax=ax_dos,
                                     color=st["color"], alpha=0.9, ls=st["ls"])
        if is_spin_polarized:
            gradient_fill_rotated(s["y_down"], dos_data.energies, ax=ax_dos,
                                  color=st["color"], alpha=0.9, ls=st["ls"])
        if st["legend"] and line is not None:
            legend_items.append((st["rank"], line, s["label"]))
    legend_items.sort(key=lambda t: t[0])

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
    legend_handles = [t[1] for t in legend_items]
    legend_labels = [t[2] for t in legend_items]

    if legend_handles:
        legend = ax_dos.legend(
            legend_handles,
            legend_labels,
            frameon=False,
            fontsize=8.0 * _fscale,
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
        out_stem = os.path.splitext(output)[0]
        band_out = out_stem + "_band.dat"
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
        dos_out = out_stem + "_dos.dat"
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
                # positive, as in valyte_dos.dat; only the plot draws spin-down below zero
                dos_cols.append(y_dn)
                dos_labels.append(f"{label}_dn")

        header = "  ".join(dos_labels)
        np.savetxt(dos_out, np.column_stack(dos_cols), header=header, fmt="%.6f")
        print(f"Saved band data: {band_out}")
        print(f"Saved DOS data: {dos_out}")
