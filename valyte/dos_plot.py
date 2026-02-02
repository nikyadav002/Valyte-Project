#!/usr/bin/env python3
"""DOS plotting utilities."""

import os
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin


def gradient_fill(x, y, ax=None, color=None, xlim=None, **kwargs):
    """Fill area under a curve with a vertical gradient."""
    if ax is None:
        ax = plt.gca()

    if len(x) == 0 or len(y) == 0:
        return None

    line, = ax.plot(x, y, color=color, lw=2, **kwargs)

    fill_color = line.get_color() if color is None else color
    alpha = line.get_alpha() or 1.0
    zorder = line.get_zorder()

    z = np.empty((100, 1, 4))
    rgb = mcolors.to_rgb(fill_color)
    z[:, :, :3] = rgb

    min_alpha = 0.05
    max_alpha = 0.95
    gradient_vector = np.linspace(min_alpha, max_alpha, 100)
    gradient_vector *= alpha

    if np.mean(y) < 0:
        gradient_vector = gradient_vector[::-1]

    z[:, :, -1] = gradient_vector[:, None]

    xmin, xmax = x.min(), x.max()

    local_ymax = max(y.max(), abs(y.min()))
    if local_ymax == 0:
        local_ymax = 1.0

    if np.mean(y) < 0:
        extent_ymin = -local_ymax
        extent_ymax = 0
    else:
        extent_ymin = 0
        extent_ymax = local_ymax

    im = ax.imshow(
        z,
        aspect="auto",
        extent=[xmin, xmax, extent_ymin, extent_ymax],
        origin="lower",
        zorder=zorder,
    )

    xy = np.column_stack([x, y])
    verts = np.vstack([[x[0], 0], xy, [x[-1], 0], [x[0], 0]])

    clip = Polygon(verts, lw=0, facecolor="none", closed=True)
    ax.add_patch(clip)
    im.set_clip_path(clip)

    return line


class ValyteDos:
    """Container for total DOS data."""

    def __init__(self, energies, densities, efermi):
        self.energies = np.array(energies)
        self.densities = densities
        self.efermi = float(efermi)

    @property
    def total(self):
        tot = np.zeros_like(self.energies)
        for spin in self.densities:
            tot += self.densities[spin]
        return tot

    @property
    def spin_up(self):
        return self.densities.get(Spin.up, np.zeros_like(self.energies))

    @property
    def spin_down(self):
        return self.densities.get(Spin.down, np.zeros_like(self.energies))


def load_dos(vasprun, elements=None, **_):
    """Load total and projected DOS from a vasprun.xml file."""
    if os.path.isdir(vasprun):
        vasprun = os.path.join(vasprun, "vasprun.xml")

    if not os.path.exists(vasprun):
        raise FileNotFoundError(f"{vasprun} not found")

    vr = Vasprun(vasprun)
    dos = vr.complete_dos
    efermi = dos.efermi

    try:
        bs = vr.get_band_structure()
        if not bs.is_metal():
            efermi = bs.get_vbm()["energy"]
    except Exception:
        try:
            cbm, vbm = dos.get_cbm_vbm()
            if cbm - vbm > 0.01:
                efermi = vbm
        except Exception:
            pass

    energies = dos.energies - efermi
    pdos = get_pdos(dos, elements)

    return ValyteDos(energies, dos.densities, efermi), pdos


def get_pdos(dos, elements=None):
    """Extract projected DOS for specified elements."""
    structure = dos.structure
    symbols = [str(site.specie) for site in structure]

    if not elements:
        unique = sorted(set(symbols))
        elements = {el: () for el in unique}
    else:
        if isinstance(elements, list):
            elements = {el: () for el in elements}

    pdos = {}
    for el in elements:
        el_sites = [s for s in structure if str(s.specie) == el]
        el_pdos = {}

        for site in el_sites:
            try:
                site_dos = dos.get_site_spd_dos(site)
            except Exception:
                continue

            for orb, orb_dos in site_dos.items():
                label = orb.name[0]
                if label not in el_pdos:
                    el_pdos[label] = {}

                for spin in orb_dos.densities:
                    if spin not in el_pdos[label]:
                        el_pdos[label][spin] = np.zeros_like(dos.energies)
                    el_pdos[label][spin] += orb_dos.densities[spin]

        pdos[el] = el_pdos
    return pdos


def plot_dos(
    dos,
    pdos,
    out="valyte_dos.png",
    xlim=(-6, 6),
    ylim=None,
    figsize=(5, 4),
    dpi=400,
    legend_loc="auto",
    font="Arial",
    show_fermi=False,
    show_total=True,
    plotting_config=None,
    legend_cutoff=0.10,
    scale_factor=1.0,
):
    """Plot total and projected DOS with the Valyte style."""

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
    mpl.rcParams["font.size"] = 12

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=figsize)

    is_spin_polarized = Spin.down in dos.densities

    if show_fermi:
        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.7)

    if is_spin_polarized:
        ax.axhline(0, color="k", lw=0.5, alpha=1.0)

    palette = [
        "#4b0082",
        "#0096c7",
        "#e63946",
        "#023e8a",
        "#ffb703",
        "#2a9d8f",
        "#8e44ad",
        "#118ab2",
        "#d62828",
        "#00b4d8",
        "#f4a261",
        "#003049",
        "#6a994e",
        "#48cae4",
        "#0077b6",
        "#90e0ef",
        "#ade8f4",
        "#caf0f8",
    ]
    lines, labels = [], []

    x_mask = (dos.energies >= xlim[0]) & (dos.energies <= xlim[1])

    if plotting_config:
        items_to_plot = plotting_config
    else:
        items_to_plot = []
        for el, el_pdos in pdos.items():
            for orb in el_pdos.keys():
                items_to_plot.append((el, orb))

    max_visible_y = 0
    min_visible_y = 0

    for i, (el, orb) in enumerate(items_to_plot):
        if el not in pdos:
            continue

        c = palette[i % len(palette)]

        if orb == "total":
            y_up = np.zeros_like(dos.energies)
            y_down = np.zeros_like(dos.energies)
            for o_data in pdos[el].values():
                y_up += o_data.get(Spin.up, np.zeros_like(dos.energies))
                y_down += o_data.get(Spin.down, np.zeros_like(dos.energies))
            label = el
        else:
            if orb not in pdos[el]:
                continue
            y_up = pdos[el][orb].get(Spin.up, np.zeros_like(dos.energies))
            y_down = pdos[el][orb].get(Spin.down, np.zeros_like(dos.energies))
            label = f"{el}({orb})"

        y_down = -y_down

        visible_y_up = y_up[x_mask]
        visible_y_down = y_down[x_mask]

        has_visible_data = False
        current_max_y = 0

        if len(visible_y_up) > 0:
            max_y = np.max(visible_y_up)
            max_visible_y = max(max_visible_y, max_y)
            current_max_y = max(current_max_y, max_y)
            if max_y > 1e-6:
                has_visible_data = True

        if is_spin_polarized and len(visible_y_down) > 0:
            min_y = np.min(visible_y_down)
            min_visible_y = min(min_visible_y, min_y)
            current_max_y = max(current_max_y, abs(min_y))
            if abs(min_y) > 1e-6:
                has_visible_data = True

        line, = ax.plot(dos.energies, y_up, lw=1.5, color=c, label=label, alpha=0)
        lines.append(
            {
                "line": line,
                "y_up": y_up,
                "y_down": y_down,
                "max_y": current_max_y,
                "color": c,
                "label": label,
                "has_visible": has_visible_data,
            }
        )

    global_max = max(max_visible_y, abs(min_visible_y))
    threshold = legend_cutoff * global_max

    final_lines = []
    final_labels = []

    for item in lines:
        line = item["line"]
        y_up = item["y_up"]
        y_down = item["y_down"]
        c = item["color"]
        label = item["label"]
        max_y = item["max_y"]
        has_visible = item["has_visible"]

        line.set_alpha(1.0)

        gradient_fill(dos.energies, y_up, ax=ax, color=c, alpha=0.9)
        if is_spin_polarized:
            gradient_fill(dos.energies, y_down, ax=ax, color=c, alpha=0.9)

        if has_visible and max_y >= threshold:
            final_lines.append(line)
            final_labels.append(label)

    lines = final_lines
    labels = final_labels

    if show_total:
        y_total_up = dos.spin_up
        y_total_down = -dos.spin_down

        ax.plot(dos.energies, y_total_up, color="k", lw=1.2, label="Total DOS")
        gradient_fill(dos.energies, y_total_up, ax=ax, color="k", alpha=0.15)

        if is_spin_polarized:
            ax.plot(dos.energies, y_total_down, color="k", lw=1.2)
            gradient_fill(dos.energies, y_total_down, ax=ax, color="k", alpha=0.15)

            visible_total_up = y_total_up[x_mask]
            visible_total_down = y_total_down[x_mask]
            if len(visible_total_up) > 0:
                max_visible_y = max(max_visible_y, np.max(visible_total_up))
            if len(visible_total_down) > 0:
                min_visible_y = min(min_visible_y, np.min(visible_total_down))

    if not ylim:
        if max_visible_y > 0 or min_visible_y < 0:
            upper_limit = (max_visible_y * 1.1) / scale_factor
            lower_limit = (min_visible_y * 1.1) / scale_factor if is_spin_polarized else 0
            ax.set_ylim(lower_limit, upper_limit)
    else:
        ax.set_ylim(*ylim)

    ax.set_xlim(*xlim)
    ax.set_xlabel("Energy (eV)", fontsize=14, weight="bold", labelpad=6)
    ax.set_ylabel("Density of States", fontsize=14, weight="bold", labelpad=6)

    xticks = np.arange(np.ceil(xlim[0]), np.floor(xlim[1]) + 1, 1)
    ax.set_xticks(xticks)
    tick_labels = [f"{int(x)}" if x == int(x) else f"{x}" for x in xticks]
    ax.set_xticklabels(tick_labels, fontweight="bold")
    ax.set_yticks([])

    if len(lines) > 0:
        legend = ax.legend(
            lines,
            labels,
            frameon=False,
            fontsize=13,
            loc="upper right" if legend_loc == "auto" else legend_loc,
            ncol=1,
            handlelength=1.5,
            columnspacing=0.8,
            handletextpad=0.6,
        )
        for text in legend.get_texts():
            text.set_fontweight("bold")

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tight_layout(pad=0.4)
    plt.savefig(out, dpi=dpi)
    plt.close(fig)
