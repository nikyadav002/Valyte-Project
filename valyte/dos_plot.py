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

from valyte.style import apply_style, get_font_weight, font_scale, energy_ticks, finish_figure
from valyte.colors import assign_styles, series_stats, load_color_file  # noqa: F401 (re-exported)


def _collect_series(dos, pdos, items_to_plot, x_mask, is_spin_polarized):
    """Build the DOS curves to draw for `items_to_plot` [(element, orbital)].

    orbital "total" sums all orbitals of the element. Missing elements or
    orbitals are skipped. Spin-down densities are returned negated, ready to
    draw below the axis.
    """
    zeros = np.zeros_like(dos.energies)
    series = []
    for el, orb in items_to_plot:
        if el not in pdos:
            continue
        if orb == "total":
            y_up, y_down = zeros.copy(), zeros.copy()
            for o_data in pdos[el].values():
                y_up = y_up + o_data.get(Spin.up, zeros)
                y_down = y_down + o_data.get(Spin.down, zeros)
            label = el
        else:
            if orb not in pdos[el]:
                continue
            y_up = pdos[el][orb].get(Spin.up, zeros)
            y_down = pdos[el][orb].get(Spin.down, zeros)
            label = f"{el}({orb})"
        peak, weight = series_stats(dos.energies, y_up,
                                    y_down if is_spin_polarized else None, x_mask)
        series.append({
            "label": label, "el": el, "orb": None if orb == "total" else orb,
            "peak": peak, "weight": weight, "y_up": y_up, "y_down": -y_down,
        })
    return series


def warn_missing_selection(pdos, plotting_config):
    """Warn about -e entries that the calculation has no data for."""
    for el, orb in plotting_config or []:
        if not pdos.get(el):
            print(f"Warning: no projected DOS for {el}. Check that it is in the "
                  "structure and that the run used LORBIT >= 10.")
        elif orb != "total" and orb not in pdos[el]:
            print(f"Warning: {el} has no {orb} projection "
                  f"(available: {', '.join(sorted(pdos[el]))}).")


_ORB_RANK = {"s": 0, "p": 1, "d": 2, "f": 3}


def _panel_styles(panels, group_by, x_mask, palette, colors, minor_color=None):
    """Colour every series of a panel figure in one pass, keyed by entity.

    Returns {(panel_index, series_index): style}. The entity is the orbital in
    element panels and the element in orbital panels, so it keeps one colour
    in every panel it appears in.
    """
    flat, where, el_rank = [], [], {}
    for idx, (panel_label, series) in enumerate(panels):
        for j, (s_label, y_up, y_down) in enumerate(series):
            if group_by == "orbital":
                el, orb, key = s_label, None, s_label
            elif "(" in s_label and s_label.endswith(")"):
                el, orb = s_label[:s_label.index("(")], s_label[s_label.index("(") + 1:-1]
                key = orb
            else:
                el, orb, key = s_label, None, "total"
            el_rank.setdefault(el, len(el_rank))
            peak, weight = series_stats(None, y_up, y_down, x_mask)
            flat.append({"label": s_label, "el": el, "orb": orb, "key": key,
                         "peak": peak, "weight": weight})
            where.append((idx, j))

    # Palette slots go to entities in a fixed order: s, p, d, f, element
    # totals, or elements in order of appearance.
    def entity_rank(item):
        k = item["key"]
        return _ORB_RANK.get(k, 4) if group_by != "orbital" else el_rank[k]

    order = sorted(range(len(flat)), key=lambda i: (entity_rank(flat[i]), i))
    styles = assign_styles([flat[i] for i in order], palette=palette,
                           colors=colors, explicit=True, minor_color=minor_color)
    return {where[i]: st for i, st in zip(order, styles)}


def gradient_fill(x, y, ax=None, color=None, xlim=None, fill=True, **kwargs):
    """Fill area under a curve with a vertical gradient.

    With ``fill=False`` only the curve is drawn.
    """
    if ax is None:
        ax = plt.gca()

    if len(x) == 0 or len(y) == 0:
        return None

    line, = ax.plot(x, y, color=color, lw=2, **kwargs)
    line.set_gid("dos-y")  # lets a glow skip the zero line
    if not fill:
        return line

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


def load_dos(vasprun, elements=None, efermi=None, vr=None, **_):
    """Load total and projected DOS from a vasprun.xml file.

    The energy zero is the VBM for a gapped system and the Fermi level for a
    metal.  Pass `efermi` to force a reference, which `combined` uses so that
    its band and DOS panels share one zero.
    """
    if vr is None:
        if os.path.isdir(vasprun):
            vasprun = os.path.join(vasprun, "vasprun.xml")

        if not os.path.exists(vasprun):
            raise FileNotFoundError(f"{vasprun} not found")

        vr = Vasprun(
            vasprun,
            parse_dos=True,
            parse_eigen=True,
            parse_projected_eigen=False,
            parse_potcar_file=False,
        )
    dos = vr.complete_dos

    if efermi is None:
        efermi = dos.efermi
        resolved = False

        try:
            bs = vr.get_band_structure()
            if bs.is_metal():
                resolved = True          # metal: the Fermi level is the zero
            else:
                efermi = bs.get_vbm()["energy"]
                resolved = True
        except Exception:
            try:
                cbm, vbm = dos.get_cbm_vbm()
                if cbm - vbm > 0.01:
                    efermi = vbm
                    resolved = True
                else:
                    resolved = True      # no gap on the DOS grid: treat as metal
            except Exception:
                resolved = False

        if not resolved:
            print("Warning: could not determine the VBM; energies are "
                  "referenced to the Fermi level from vasprun.xml.")

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


def _save_dos_dat(dos, pdos, items_to_plot, filepath):
    """Save DOS data (energies, total, and PDOS) to a text file."""
    is_spin_polarized = Spin.down in dos.densities
    energies = dos.energies

    cols = [energies, dos.spin_up]
    col_labels = ["Energy(eV)", "Total_up"]

    if is_spin_polarized:
        cols.append(dos.spin_down)
        col_labels.append("Total_dn")

    for el, orb in items_to_plot:
        if el not in pdos:
            continue
        if orb == "total":
            y_up = sum(o.get(Spin.up, np.zeros_like(energies)) for o in pdos[el].values())
            label = el
        else:
            if orb not in pdos[el]:
                continue
            y_up = pdos[el][orb].get(Spin.up, np.zeros_like(energies))
            label = f"{el}({orb})"

        cols.append(y_up)
        col_labels.append(f"{label}_up")

        if is_spin_polarized:
            if orb == "total":
                y_dn = sum(o.get(Spin.down, np.zeros_like(energies)) for o in pdos[el].values())
            else:
                y_dn = pdos[el][orb].get(Spin.down, np.zeros_like(energies))
            cols.append(y_dn)
            col_labels.append(f"{label}_dn")

    header = "  ".join(col_labels)
    np.savetxt(filepath, np.column_stack(cols), header=header, fmt="%.6f")
    print(f"Saved data: {filepath}")


def plot_dos(
    dos,
    pdos,
    out="valyte_dos.png",
    xlim=(-6, 6),
    ylim=None,
    figsize=(4.45, 4.2),
    dpi=400,
    legend_loc="auto",
    font=None,
    show_fermi=False,
    show_total=True,
    plotting_config=None,
    legend_cutoff=0.10,
    scale_factor=1.0,
    save_data=False,
    bold=True,
    colors=None,
    fontsize=None,
    palette=None,
    style=None,
):
    """Plot total and projected DOS with the Valyte style.

    `palette` names a colour palette ('wong', the default, 'bright' or
    'muted') or gives a list of colours; `colors` maps individual series
    (e.g. {"Fe(d)": "red"}) and always takes precedence.
    """

    _fsbase = fontsize if fontsize is not None else 12
    _fscale = font_scale(fontsize, 12)
    st = apply_style(font=font, bold=bold, fontsize=_fsbase, style=style)
    bold = bold and st.bold
    _weight = get_font_weight(bold)
    fig, ax = plt.subplots(figsize=figsize)

    is_spin_polarized = Spin.down in dos.densities

    if show_fermi:
        ax.axvline(0, color=st.fg, lw=0.8, ls="--", alpha=0.7)

    if is_spin_polarized:
        ax.axhline(0, color=st.fg, lw=0.5, alpha=1.0)

    palette = palette if palette is not None else st.palette

    x_mask = (dos.energies >= xlim[0]) & (dos.energies <= xlim[1])

    if plotting_config:
        items_to_plot = plotting_config
        warn_missing_selection(pdos, plotting_config)
    else:
        items_to_plot = []
        for el, el_pdos in pdos.items():
            for orb in el_pdos.keys():
                items_to_plot.append((el, orb))

    series = _collect_series(dos, pdos, items_to_plot, x_mask, is_spin_polarized)
    styles = assign_styles(series, palette=palette, colors=colors,
                           cutoff=legend_cutoff, explicit=bool(plotting_config),
                           minor_color=st.minor_color)

    max_visible_y = max((s["y_up"][x_mask].max(initial=0.0) for s in series), default=0.0)
    min_visible_y = 0
    if is_spin_polarized:
        min_visible_y = min((s["y_down"][x_mask].min(initial=0.0) for s in series), default=0.0)

    # Grey background series first, then the legend series largest-first, so
    # smaller channels are drawn on top of the bigger fills.
    order = sorted(range(len(series)),
                   key=lambda i: (styles[i]["rank"] is not None, styles[i]["rank"] or 0))
    legend_items = []
    for i in order:
        s, sty = series[i], styles[i]
        line = gradient_fill(dos.energies, s["y_up"], ax=ax, color=sty["color"],
                             alpha=0.9, ls=sty["ls"], fill=st.fill)
        if is_spin_polarized:
            gradient_fill(dos.energies, s["y_down"], ax=ax, color=sty["color"],
                          alpha=0.9, ls=sty["ls"], fill=st.fill)
        if sty["legend"] and line is not None:
            legend_items.append((sty["rank"], line, s["label"]))

    legend_items.sort(key=lambda t: t[0])
    lines = [t[1] for t in legend_items]
    labels = [t[2] for t in legend_items]

    if show_total:
        y_total_up = dos.spin_up
        y_total_down = -dos.spin_down

        ax.plot(dos.energies, y_total_up, color=st.total_color, lw=1.2, label="Total DOS",
                gid="dos-y")
        if st.fill:
            # No path effects on the fill's translucent outline: a style's halo
            # (xkcd) would otherwise paint over the solid total line.
            gradient_fill(dos.energies, y_total_up, ax=ax, color=st.total_color, alpha=0.15,
                          path_effects=[])

        if is_spin_polarized:
            ax.plot(dos.energies, y_total_down, color=st.total_color, lw=1.2, gid="dos-y")
            if st.fill:
                gradient_fill(dos.energies, y_total_down, ax=ax, color=st.total_color, alpha=0.15,
                              path_effects=[])

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
    ax.set_xlabel("Energy (eV)", fontsize=14 * _fscale, weight=_weight, labelpad=6)
    ax.set_ylabel("Density of States", fontsize=14 * _fscale, weight=_weight, labelpad=6)

    xticks = energy_ticks(xlim[0], xlim[1])
    ax.set_xticks(xticks)
    tick_labels = [f"{int(x)}" if x == int(x) else f"{x}" for x in xticks]
    ax.set_xticklabels(tick_labels, fontweight=_weight)
    ax.tick_params(axis="x", direction="in", which="both", bottom=True, top=False)
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_yticks([])

    if len(lines) > 0:
        legend = ax.legend(
            lines,
            labels,
            frameon=False,
            fontsize=13 * _fscale,
            loc="upper right" if legend_loc == "auto" else legend_loc,
            ncol=1,
            handlelength=1.5,
            columnspacing=0.8,
            handletextpad=0.6,
        )
        for text in legend.get_texts():
            text.set_fontweight(_weight)

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.tight_layout(pad=0.4)
    finish_figure(fig, dpi)
    plt.savefig(out, dpi=dpi)
    plt.close(fig)

    if save_data:
        _save_dos_dat(dos, pdos, items_to_plot, os.path.splitext(out)[0] + ".dat")


def plot_dos_panels(
    dos,
    pdos,
    out="valyte_dos_panels.png",
    xlim=(-6, 6),
    ylim=None,
    figsize=None,
    dpi=400,
    font=None,
    show_fermi=False,
    show_total=True,
    plotting_config=None,
    scale_factor=1.0,
    save_data=False,
    group_by="element",
    bold=True,
    colors=None,
    fontsize=None,
    palette=None,
    style=None,
):
    """Plot DOS as vertically stacked panels — one per element (or per orbital).

    Parameters
    ----------
    group_by : str
        ``"element"`` (default) — one panel per element, orbitals shown within.
        ``"orbital"`` — one panel per angular-momentum channel (s, p, d, f),
        elements shown within each.
    palette, colors
        As in `plot_dos`. Colours follow the entity across panels: with
        ``group_by="element"`` each orbital keeps one colour in every panel,
        with ``group_by="orbital"`` each element does.
    """

    _fsbase = fontsize if fontsize is not None else 11
    _fscale = font_scale(fontsize, 11)
    st = apply_style(font=font, bold=bold, fontsize=_fsbase, style=style)
    bold = bold and st.bold
    _weight = get_font_weight(bold)
    palette = palette if palette is not None else st.palette

    is_spin_polarized = Spin.down in dos.densities
    x_mask = (dos.energies >= xlim[0]) & (dos.energies <= xlim[1])
    if plotting_config:
        warn_missing_selection(pdos, plotting_config)

    # ------------------------------------------------------------------
    # Build panel groups: list of (panel_label, [(series_label, y_up, y_down), ...])
    # ------------------------------------------------------------------
    panels = []

    if group_by == "orbital":
        # Collect all orbital channels across all elements
        all_orbitals = []
        for el, el_pdos in pdos.items():
            for orb in el_pdos:
                if orb not in all_orbitals:
                    all_orbitals.append(orb)

        # Sort by angular momentum
        orb_order = {"s": 0, "p": 1, "d": 2, "f": 3}
        all_orbitals.sort(key=lambda o: orb_order.get(o, 99))

        for orb in all_orbitals:
            series = []
            for el, el_pdos in pdos.items():
                if orb not in el_pdos:
                    continue
                y_up = el_pdos[orb].get(Spin.up, np.zeros_like(dos.energies))
                y_down = el_pdos[orb].get(Spin.down, np.zeros_like(dos.energies))
                series.append((el, y_up, -y_down))
            if series:
                panels.append((orb, series))
    else:
        # Default: one panel per element
        if plotting_config:
            # Group plotting_config entries by element
            from collections import OrderedDict
            el_groups = OrderedDict()
            for el, orb in plotting_config:
                el_groups.setdefault(el, []).append(orb)

            for el, orbs in el_groups.items():
                if el not in pdos:
                    continue
                series = []
                for orb in orbs:
                    if orb == "total":
                        y_up = np.zeros_like(dos.energies)
                        y_down = np.zeros_like(dos.energies)
                        for o_data in pdos[el].values():
                            y_up += o_data.get(Spin.up, np.zeros_like(dos.energies))
                            y_down += o_data.get(Spin.down, np.zeros_like(dos.energies))
                        series.append((el, y_up, -y_down))
                    else:
                        if orb not in pdos[el]:
                            continue
                        y_up = pdos[el][orb].get(Spin.up, np.zeros_like(dos.energies))
                        y_down = pdos[el][orb].get(Spin.down, np.zeros_like(dos.energies))
                        series.append((f"{el}({orb})", y_up, -y_down))
                if series:
                    panels.append((el, series))
        else:
            for el, el_pdos in pdos.items():
                series = []
                for orb in el_pdos:
                    y_up = el_pdos[orb].get(Spin.up, np.zeros_like(dos.energies))
                    y_down = el_pdos[orb].get(Spin.down, np.zeros_like(dos.energies))
                    series.append((f"{el}({orb})", y_up, -y_down))
                if series:
                    panels.append((el, series))

    n_panels = len(panels)
    if n_panels == 0:
        print("Warning: no data to plot in panel mode.")
        return

    # Compute figure size if not given
    if figsize is None:
        panel_h = 2.4
        figsize = (5, panel_h * n_panels + 0.6)

    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=figsize,
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )
    if n_panels == 1:
        axes = [axes]

    # Pre-compute total DOS for backdrop
    y_total_up = dos.spin_up
    y_total_down = -dos.spin_down

    styles = _panel_styles(panels, group_by, x_mask, palette, colors, st.minor_color)

    for idx, (panel_label, series) in enumerate(panels):
        ax = axes[idx]

        if show_fermi:
            ax.axvline(0, color=st.fg, lw=0.8, ls="--", alpha=0.5, zorder=0)

        if is_spin_polarized:
            ax.axhline(0, color=st.fg, lw=0.4, alpha=0.8)

        # Draw semi-transparent total DOS as backdrop
        if show_total and st.fill:
            ax.fill_between(
                dos.energies, 0, y_total_up,
                color=st.backdrop, alpha=0.35, lw=0, zorder=0,
            )
            if is_spin_polarized:
                ax.fill_between(
                    dos.energies, 0, y_total_down,
                    color=st.backdrop, alpha=0.35, lw=0, zorder=0,
                )
        elif show_total:
            ax.plot(dos.energies, y_total_up, color=st.backdrop_line, lw=1.0, zorder=0)
            if is_spin_polarized:
                ax.plot(dos.energies, y_total_down, color=st.backdrop_line, lw=1.0, zorder=0)

        panel_max_y = 0
        panel_min_y = 0

        legend_lines = []
        legend_labels = []

        for j, (s_label, y_up, y_down) in enumerate(series):
            sty = styles[(idx, j)]
            c = sty["color"]

            vis_up = y_up[x_mask]
            vis_dn = y_down[x_mask]

            if len(vis_up) > 0:
                panel_max_y = max(panel_max_y, np.max(vis_up))
            if is_spin_polarized and len(vis_dn) > 0:
                panel_min_y = min(panel_min_y, np.min(vis_dn))

            l = gradient_fill(dos.energies, y_up, ax=ax, color=c, alpha=0.85, ls=sty["ls"], fill=st.fill)
            if is_spin_polarized:
                gradient_fill(dos.energies, y_down, ax=ax, color=c, alpha=0.85, ls=sty["ls"], fill=st.fill)

            if l is not None and sty["legend"]:
                legend_lines.append(l)
                legend_labels.append(s_label)

        # Y-limits
        if ylim:
            ax.set_ylim(*ylim)
        else:
            if show_total:
                vis_total_up = y_total_up[x_mask]
                vis_total_dn = y_total_down[x_mask]
                if len(vis_total_up) > 0:
                    panel_max_y = max(panel_max_y, np.max(vis_total_up))
                if is_spin_polarized and len(vis_total_dn) > 0:
                    panel_min_y = min(panel_min_y, np.min(vis_total_dn))

            upper = (panel_max_y * 1.15) / scale_factor if panel_max_y > 0 else 1
            lower = (panel_min_y * 1.15) / scale_factor if is_spin_polarized else 0
            ax.set_ylim(lower, upper)

        ax.set_xlim(*xlim)
        ax.set_yticks([])

        # Panel label — placed as a text annotation inside the panel
        ax.text(
            0.02, 0.88, panel_label,
            transform=ax.transAxes,
            fontsize=13 * _fscale, fontweight=_weight,
            va="top", ha="left",
            bbox=dict(
                facecolor=st.bg, edgecolor="none",
                alpha=0.75, boxstyle="round,pad=0.2",
            ),
        )

        # Per-panel legend for orbital breakdown
        if len(legend_lines) > 1:
            leg = ax.legend(
                legend_lines, legend_labels,
                frameon=False, fontsize=10 * _fscale,
                loc="upper right",
                ncol=1,
                handlelength=1.2,
                handletextpad=0.4,
                borderpad=0.3,
            )
            for t in leg.get_texts():
                t.set_fontweight(_weight)

        ax.tick_params(axis="x", direction="in", which="both", bottom=True, top=False)
        ax.tick_params(axis="y", which="both", left=False, right=False)
        if idx < n_panels - 1:
            ax.tick_params(axis="x", labelbottom=False)

        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

        # Light border between panels
        for spine in ("top", "bottom"):
            ax.spines[spine].set_linewidth(0.6)
            ax.spines[spine].set_color(st.separator)

    # Bottom panel: energy axis label and ticks
    bottom_ax = axes[-1]
    xticks = energy_ticks(xlim[0], xlim[1])
    bottom_ax.set_xticks(xticks)
    tick_labels = [f"{int(x)}" if x == int(x) else f"{x}" for x in xticks]
    bottom_ax.set_xticklabels(tick_labels, fontweight=_weight)
    bottom_ax.set_xlabel("Energy (eV)", fontsize=14 * _fscale, weight=_weight, labelpad=6)

    # Shared Y-label in the centre
    fig.text(
        0.01, 0.5, "Density of States",
        va="center", ha="left", rotation="vertical",
        fontsize=14 * _fscale, fontweight=_weight,
    )

    fig.subplots_adjust(left=0.10, right=0.97, top=0.97, bottom=0.08)
    finish_figure(fig, dpi)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)

    if save_data:
        items_to_plot = []
        for _, series in panels:
            for s_label, _, _ in series:
                # reconstruct (el, orb) from label
                if "(" in s_label and s_label.endswith(")"):
                    el_part = s_label[:s_label.index("(")]
                    orb_part = s_label[s_label.index("(") + 1:-1]
                    items_to_plot.append((el_part, orb_part))
                else:
                    items_to_plot.append((s_label, "total"))
        _save_dos_dat(dos, pdos, items_to_plot, os.path.splitext(out)[0] + ".dat")
