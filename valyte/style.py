"""Shared Matplotlib plotting style and helper utilities for Valyte."""

import os
from dataclasses import dataclass, field

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patheffects

FONT_MAP = {
    "arial": "Arial",
    "helvetica": "Helvetica",
    "times": "Times New Roman",
    "times new roman": "Times New Roman",
}

# Palettes live in valyte.colors; this name is kept for existing imports.
from valyte.colors import WONG as _WONG, MINOR_COLOR as _MINOR_COLOR
DEFAULT_PALETTE = list(_WONG)

# Fonts tried in order for the xkcd style; the first one installed wins.
XKCD_FONTS = ["xkcd Script", "xkcd", "Humor Sans", "Comic Neue", "Comic Sans MS"]


@dataclass(frozen=True)
class PlotStyle:
    """A named look shared by every plotting module.

    The defaults reproduce the original Valyte figures; each preset only
    overrides what it changes.
    """

    name: str
    palette: tuple = tuple(DEFAULT_PALETTE)
    vb_color: str = "#8e44ad"
    cb_color: str = "#2a9d8f"
    spin_colors: tuple = ("#3498db", "#e74c3c")
    total_color: str = "k"
    bz_colors: dict = field(default_factory=lambda: {
        "edge": "#2b2d42", "face": "#8d99ae", "path": "#e63946",
        "point": "#1d3557", "vector": "#2a9d8f",
    })
    fill: bool = True        # gradient fill under DOS curves
    bold: bool = True        # False forces normal weight regardless of --no-bold
    fonts: tuple = ()        # preferred families when --font is not given
    sketch: tuple = None     # (scale, length, randomness) at 100 dpi
    legacy_axes: bool = False
    glow: bool = False       # soft glow under data lines
    # Colours the plot code would otherwise hard-code for a white background
    fg: str = "k"                    # reference lines (E = 0, k-path breaks)
    bg: str = "white"                # label boxes, marker edges, text outlines
    muted: str = "0.35"              # thin outlines (tricolor triangle)
    minor_color: str = _MINOR_COLOR  # DOS channels below the legend cutoff
    backdrop: str = "#d4d4d4"        # total-DOS backdrop in panel plots
    backdrop_line: str = "#9a9a9a"   # the same backdrop when fills are off
    separator: str = "#888888"       # borders between DOS panels
    rc: dict = field(default_factory=dict)


def _dark_rc(bg, fg, frame, grid=None):
    """rcParams for a dark background: light text and frame, matching legends."""
    rc = {
        "figure.facecolor": bg, "axes.facecolor": bg, "savefig.facecolor": bg,
        "text.color": fg, "axes.labelcolor": fg, "axes.edgecolor": frame,
        "xtick.color": frame, "ytick.color": frame,
        "xtick.labelcolor": fg, "ytick.labelcolor": fg,
        "legend.facecolor": bg, "legend.edgecolor": frame,
    }
    if grid:
        rc.update({"axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6,
                   "axes.axisbelow": True})
    return rc


STYLES = {
    "default": PlotStyle("default"),
    "xkcd": PlotStyle(
        "xkcd",
        palette=(
            "xkcd:tomato", "xkcd:cerulean", "xkcd:kelly green", "xkcd:tangerine",
            "xkcd:purple", "xkcd:teal", "xkcd:hot pink", "xkcd:mustard",
            "xkcd:navy blue", "xkcd:brick red", "xkcd:leaf green", "xkcd:sky blue",
        ),
        vb_color="xkcd:purple",
        cb_color="xkcd:teal",
        spin_colors=("xkcd:cerulean", "xkcd:tomato"),
        bz_colors={
            "edge": "k", "face": "xkcd:light grey", "path": "xkcd:tomato",
            "point": "xkcd:navy blue", "vector": "xkcd:teal",
        },
        bold=False,
        fonts=tuple(XKCD_FONTS),
        sketch=(1.0, 100.0, 2.0),
        rc={
            "path.effects": [patheffects.withStroke(linewidth=4, foreground="w")],
            "axes.linewidth": 1.5,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
            "xtick.minor.width": 1.2,
            "ytick.minor.width": 1.2,
        },
    ),
    "legacy": PlotStyle(
        "legacy",
        # Classic saturated line colours; black is kept for the total
        palette=(
            "#FF0000", "#0000FF", "#008000", "#FF00FF", "#008080", "#808000",
            "#000080", "#800080", "#800000", "#FF8000", "#8000FF", "#0000A0",
        ),
        vb_color="#000000",
        cb_color="#FF0000",
        spin_colors=("#000000", "#FF0000"),
        bz_colors={
            "edge": "#000000", "face": "#C0C0C0", "path": "#FF0000",
            "point": "#000000", "vector": "#0000FF",
        },
        fill=False,
        bold=False,
        fonts=("Arial", "Helvetica"),
        legacy_axes=True,
        rc={
            "axes.linewidth": 1.5,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
        },
    ),
    "lines": PlotStyle("lines", fill=False),
    "dark": PlotStyle(
        "dark",
        # The default palette keeps >= 3:1 contrast on this background.
        vb_color="#B784E0",
        cb_color="#3FC1B0",
        total_color="#F2F2F2",
        bz_colors={
            "edge": "#BBBBBB", "face": "#5A5A5A", "path": "#E69F00",
            "point": "#EDEDED", "vector": "#56B4E9",
        },
        fg="#EDEDED", bg="#141414", muted="#9A9A9A", minor_color="#5C5C5C",
        backdrop="#8A8A8A", backdrop_line="#7A7A7A", separator="#5A5A5A",
        rc=_dark_rc("#141414", "#EDEDED", "#BBBBBB"),
    ),
    "neon": PlotStyle(
        "neon",
        # Five colours: neon green is left out, it matches the yellow for
        # red-green colour-blind readers.
        palette=("#08F7FE", "#FE53BB", "#F5D300", "#FF6C11", "#9467FF"),
        vb_color="#FE53BB",
        cb_color="#08F7FE",
        spin_colors=("#08F7FE", "#FE53BB"),
        total_color="#E6E9FF",
        bz_colors={
            "edge": "#8B93C9", "face": "#2A3459", "path": "#FE53BB",
            "point": "#E6E9FF", "vector": "#08F7FE",
        },
        glow=True,
        fg="#E6E9FF", bg="#212946", muted="#8B93C9", minor_color="#4A5380",
        backdrop="#6B74A8", backdrop_line="#5A6394", separator="#3A4470",
        rc=_dark_rc("#212946", "#E6E9FF", "#8B93C9", grid="#2A3459"),
    ),
    "grid": PlotStyle(
        "grid",
        rc={"axes.grid": True, "grid.color": "#DADADA", "grid.linewidth": 0.6,
            "axes.axisbelow": True},
    ),
}

STYLE_NAMES = list(STYLES)

_active = STYLES["default"]


def active_style() -> PlotStyle:
    """The style applied by the most recent apply_style() call."""
    return _active


def get_style(style=None) -> PlotStyle:
    """Look up a style by name; a PlotStyle instance is returned as-is."""
    if isinstance(style, PlotStyle):
        return style
    name = (style or "default").lower()
    if name not in STYLES:
        raise ValueError(f"Unknown style '{style}'. Choose from: {', '.join(STYLE_NAMES)}")
    return STYLES[name]


def _first_installed(families):
    from matplotlib import font_manager
    installed = {f.name for f in font_manager.fontManager.ttflist}
    return next((f for f in families if f in installed), None)


def resolve_font_family(font: str = "Arial") -> str:
    """Resolve font family name from supported map."""
    if not font:
        return "Arial"
    return FONT_MAP.get(font.lower(), "Arial")


def get_font_weight(bold: bool = True) -> str:
    """Return 'bold' or 'normal' weight string."""
    return "bold" if bold else "normal"


def energy_ticks(lo, hi, max_ticks=9):
    """Pick readable tick positions for an energy axis spanning [lo, hi].

    Wide spans keep the 1 eV integer steps the plots were designed around, so
    the default ranges are unchanged.  Narrow spans fall back to a locator that
    chooses a finer step, instead of collapsing to a single tick at zero.
    """
    import numpy as np
    from matplotlib.ticker import MaxNLocator

    ticks = np.arange(np.ceil(lo), np.floor(hi) + 1, 1)
    if len(ticks) >= 3:
        return ticks

    ticks = MaxNLocator(nbins=max_ticks - 1, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    return np.asarray([t for t in ticks if lo - 1e-9 <= t <= hi + 1e-9])


def font_scale(fontsize, reference):
    """Factor mapping a module's reference design size onto `fontsize`.

    Each plotting module was laid out against a base font size, with the
    individual label, tick and legend sizes chosen relative to it.  Scaling
    every one of those by this factor keeps the proportions intact when the
    user asks for a different base.  Returns 1.0 when no size is requested,
    so the default output is unchanged.
    """
    if fontsize is None or not reference:
        return 1.0
    return float(fontsize) / float(reference)


def apply_style(
    font: str = None,
    bold: bool = True,
    fontsize: float = 12,
    linewidth: float = None,
    style=None,
) -> PlotStyle:
    """Apply consistent Matplotlib rcParams across all Valyte plotting modules.

    Returns the active PlotStyle so callers can pick up its colours.  A font
    passed explicitly always wins over the style's preferred fonts.
    """
    global _active
    st = get_style(style)
    _active = st
    bold = bold and st.bold

    plt.style.use("default")
    mpl.use("agg")
    mpl.rcParams["axes.unicode_minus"] = False

    weight = get_font_weight(bold)
    if font or not st.fonts:
        font_family = resolve_font_family(font)
    else:
        font_family = _first_installed(st.fonts) or "DejaVu Sans"
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

    mpl.rcParams.update(st.rc)
    if st.sketch:
        mpl.rcParams["path.sketch"] = st.sketch
        # Hand-lettered maths to match the text font
        mpl.rcParams["mathtext.fontset"] = "custom"
        mpl.rcParams["mathtext.rm"] = font_family
        mpl.rcParams["mathtext.it"] = f"{font_family}:italic"
        mpl.rcParams["mathtext.bf"] = f"{font_family}:bold"
    return st


def finish_figure(fig, dpi: int = 400):
    """Apply the active style's last touches that rcParams cannot express.

    Called right before saving, after the module has set its own tick and
    legend options.
    """
    st = _active

    if st.sketch:
        # Sketch amplitude is in device pixels; scale it so the wobble looks
        # the same at any --dpi.
        f = dpi / 100.0
        for artist in fig.findobj():
            params = artist.get_sketch_params() if hasattr(artist, "get_sketch_params") else None
            if params:
                artist.set_sketch_params(params[0] * f, params[1] * f, params[2])

        # The white halo is there to separate crossing data lines. On the
        # frame and ticks, which are drawn above the data, it would erase the
        # curves and fills that run along the axes, so take it off those.
        for ax in fig.axes:
            for spine in ax.spines.values():
                spine.set_path_effects([])
            for axis in (ax.xaxis, ax.yaxis):
                for tick in axis.get_ticklines() + axis.get_ticklines(minor=True):
                    tick.set_path_effects([])

    if st.glow:
        _add_glow(fig)

    if st.legacy_axes:
        for ax in fig.axes:
            if not ax.axison or ax.name == "3d":
                continue
            # Ticks outward on the bottom and left axes only, closed black frame
            ax.tick_params(which="both", direction="out", top=False, right=False)
            for spine in ax.spines.values():
                spine.set_color("k")
                spine.set_linewidth(mpl.rcParams["axes.linewidth"])
            legend = ax.get_legend()
            if legend is not None:
                legend.set_frame_on(True)
                frame = legend.get_frame()
                frame.set_boxstyle("square", pad=0.0)
                frame.set_edgecolor("k")
                frame.set_facecolor("white")
                frame.set_linewidth(1.0)
                frame.set_alpha(1.0)
        # Outward ticks need more room than the layout was fitted for; widen
        # the margins just enough that no label is cut off at the edges.
        _fit_inside(fig)


def _add_glow(fig, layers=10, spread=1.05, strength=0.3):
    """Give every data line a soft glow: wider, fainter copies drawn beneath it.

    Skipped for reference lines (axhline/axvline), translucent helper lines,
    and DOS channels below the legend cutoff. DOS curves (tagged by the DOS
    plotters) get no glow where they lie along their zero line, which would
    otherwise smear colour along the frame.
    """
    import numpy as np
    from matplotlib.lines import Line2D

    minor = mpl.colors.to_hex(_active.minor_color)
    for ax in fig.axes:
        if ax.name == "3d":
            continue
        lines = [line for line in ax.get_lines()
                 if line.get_visible() and line.get_transform() == ax.transData
                 and (line.get_alpha() if line.get_alpha() is not None else 1) >= 0.5
                 and line.get_linewidth() > 0
                 and mpl.colors.to_hex(line.get_color()) != minor]
        # Many overlapping glows (a band structure has dozens of lines) add up
        # to solid blobs, so the glow fades as the line count grows.
        fade = min(1.0, (6 / len(lines)) ** 0.5) if lines else 1.0
        for line in lines:
            x = np.asarray(line.get_xdata(), dtype=float)
            y = np.asarray(line.get_ydata(), dtype=float)
            gid = line.get_gid() or ""
            if gid in ("dos-y", "dos-x"):
                values, (lo, hi) = (y, ax.get_ylim()) if gid == "dos-y" else (x, ax.get_xlim())
                quiet = np.abs(values) < 0.01 * abs(hi - lo)
                if gid == "dos-y":
                    y = np.where(quiet, np.nan, y)
                else:
                    x = np.where(quiet, np.nan, x)
            for k in range(1, layers + 1):
                ax.add_line(Line2D(x, y, color=line.get_color(), alpha=strength * fade / layers,
                                   linewidth=line.get_linewidth() + spread * k,
                                   linestyle=line.get_linestyle(), solid_capstyle="round",
                                   zorder=line.get_zorder() - 0.01))


def _fit_inside(fig):
    """Shrink the subplot area so every artist lies within the figure."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    box = fig.get_tightbbox(renderer)              # inches
    w, h = fig.get_size_inches()
    pad = 0.02
    over_left = max(0.0, -box.x0 + pad) / w
    over_bottom = max(0.0, -box.y0 + pad) / h
    over_right = max(0.0, box.x1 - w + pad) / w
    over_top = max(0.0, box.y1 - h + pad) / h
    if not any((over_left, over_bottom, over_right, over_top)):
        return
    sp = fig.subplotpars
    fig.subplots_adjust(left=sp.left + over_left, bottom=sp.bottom + over_bottom,
                        right=sp.right - over_right, top=sp.top - over_top)


def save_plot(fig, output: str, dpi: int = 400):
    """Save matplotlib figure with specified DPI and tight layout, then close."""
    plt.tight_layout()
    finish_figure(fig, dpi)
    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")
