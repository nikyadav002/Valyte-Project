"""Colour palettes and colour assignment for DOS curves.

Every DOS figure (single plot, stacked panels, and the DOS side of the
combined plot) gets its colours from `assign_styles`, so the rules below hold
everywhere:

* Custom colours (a mapping such as ``{"Fe(d)": "red"}``) always win.
* Series shown in the legend take palette colours in order. In automatic mode
  that order is by weight (area in the energy window), so the strongest colours
  go to the channels that matter most; series the user named explicitly keep
  the order they were given.
* Series below the legend cutoff are drawn in a neutral grey and left out of
  the legend, so no colour ever appears without a label.
* With more series than colours, colours repeat with a different line style,
  so no two legend entries look the same.
"""

import json
import os
import re

from matplotlib.colors import is_color_like, to_hex

# Bang Wong's colour-blind-safe palette (2011). Black is reserved for the
# total DOS, and yellow is left out as too faint for a line on white.
WONG = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]

# Paul Tol's qualitative schemes (https://personal.sron.nl/~pault/), without
# their trailing grey, which here means "below the legend cutoff".
BRIGHT = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
MUTED = ["#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255",
         "#44AA99", "#999933", "#AA4499"]

PALETTES = {"wong": WONG, "bright": BRIGHT, "muted": MUTED}
PALETTE_NAMES = list(PALETTES)
DEFAULT_PALETTE = "wong"

MINOR_COLOR = "#BDBDBD"
LINESTYLES = ["-", "--", ":", "-."]

_ORBITALS = ("s", "p", "d", "f")
_BARE_HEX = re.compile(r"^(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})$")
_SERIES_KEY = re.compile(r"^([A-Za-z]{1,3})\s*(?:[\(:]\s*([A-Za-z]+)\s*\)?)?$")


def normalize_color(value, what="colour"):
    """Return `value` as a Matplotlib colour, accepting hex without the '#'.

    Raises ValueError with the offending value instead of letting Matplotlib
    fail later in the middle of drawing.
    """
    if not isinstance(value, str):
        if is_color_like(value):
            return value
        raise ValueError(f"Invalid {what}: {value!r}")
    text = value.strip()
    if _BARE_HEX.match(text):
        text = "#" + text
    if not text or not is_color_like(text):
        raise ValueError(
            f"Invalid {what}: {value!r}. Use a colour name (e.g. 'red'), "
            "a hex code (e.g. '#0072B2' or '0072B2'), or 'xkcd:...' names.")
    return text


def get_palette(palette=None):
    """Resolve a palette name or an explicit list of colours into a list."""
    if palette is None:
        palette = DEFAULT_PALETTE
    if isinstance(palette, str):
        name = palette.strip().lower()
        if name not in PALETTES:
            raise ValueError(
                f"Unknown palette {palette!r}. Choose from: {', '.join(PALETTE_NAMES)}")
        return list(PALETTES[name])
    colours = [normalize_color(c, "palette colour") for c in palette]
    if not colours:
        raise ValueError("A custom palette needs at least one colour.")
    return colours


def normalize_series_key(key):
    """Canonical form of a colour-mapping key.

    'pb(P)', 'Pb (p)' and 'Pb:p' all become 'Pb(p)'; 'fe' becomes 'Fe'.
    A single lowercase s/p/d/f is an orbital, so 's' and 'S' (sulfur) stay
    distinct. Keys that match neither form are returned stripped, unchanged.
    """
    text = str(key).strip()
    if text in _ORBITALS:
        return text
    m = _SERIES_KEY.match(text)
    if not m:
        return text
    el = m.group(1).capitalize()
    orb = m.group(2)
    if orb is None:
        return el
    orb = orb.lower()
    return el if orb == "total" else f"{el}({orb})"


def normalize_mapping(mapping):
    """Validate a {series: colour} mapping and canonicalise its keys."""
    out = {}
    for key, value in (mapping or {}).items():
        out[normalize_series_key(key)] = normalize_color(value, f"colour for {key!r}")
    return out


def split_colors(colors):
    """Split the `colors` argument into (mapping, palette_override).

    `colors` may be None, a dict (per-series colours), or a list (a custom
    palette, kept for backward compatibility with earlier Valyte versions).
    """
    if colors is None:
        return {}, None
    if isinstance(colors, dict):
        return normalize_mapping(colors), None
    if isinstance(colors, (list, tuple)):
        return {}, get_palette(list(colors)) if colors else None
    raise TypeError("colors must be a dict of series colours, a list of colours, or None")


def load_color_file(filepath):
    """Load a {series: colour} mapping from a JSON object or a key/value text file.

    Text format: one ``key value`` pair per line, separated by whitespace, ':'
    or '='; blank lines and lines starting with '#' are ignored.
    """
    if not filepath:
        return {}
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Colours file not found: {filepath}")

    if filepath.lower().endswith(".json"):
        with open(filepath, "r") as f:
            try:
                raw = json.load(f)
            except json.JSONDecodeError as e:
                raise ValueError(f"{filepath}: invalid JSON ({e})") from None
        if not isinstance(raw, dict):
            raise ValueError(f"{filepath}: expected a JSON object mapping series to colours")
        return normalize_mapping({k: str(v) for k, v in raw.items()})

    mapping = {}
    with open(filepath, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Split on the last separator so keys like "Fe(d)" and "Fe:d" work.
            m = re.match(r"^(.+?)\s*(?:=|\s|:(?![spdfSPDF]\b))\s*(\S+)$", line)
            if not m:
                raise ValueError(f"{filepath}, line {lineno}: expected 'series colour', got {line!r}")
            mapping[m.group(1)] = m.group(2)
    return normalize_mapping(mapping)


def assign_styles(series, palette=None, colors=None, cutoff=0.10, explicit=False,
                  warn=print):
    """Decide colour, line style and legend membership for each DOS series.

    Parameters
    ----------
    series : list of dict
        One entry per curve, in drawing order, with keys
        ``label`` (legend text, e.g. "Pb(p)"), ``el``, ``orb`` (or None),
        ``peak`` (largest |DOS| in the window) and ``weight`` (area in the
        window). Optional ``key``: series sharing a key share a palette colour
        (used by the panel plot so each orbital keeps one colour throughout).
    palette : str or list, optional
        Palette name ('wong', 'bright', 'muted') or a list of colours.
    colors : dict or list, optional
        Per-series colour mapping, or a list used as the palette.
    cutoff : float
        Fraction of the largest peak a series needs to enter the legend in
        automatic mode.
    explicit : bool
        True when the user chose the series (e.g. with -e): every visible
        series is shown, in the given order, whatever its size.

    Returns
    -------
    styles : list of dict
        Per series: ``color``, ``ls``, ``legend`` (bool), ``rank`` (legend
        position, or None).
    """
    if not 0 <= cutoff <= 1:
        raise ValueError(f"legend cutoff must be between 0 and 1, got {cutoff}")

    mapping, palette_override = split_colors(colors)
    colours = palette_override or get_palette(palette)

    peaks = [s["peak"] for s in series]
    top = max(peaks, default=0.0)
    visible = [p > 1e-6 for p in peaks]

    def custom(s):
        for k in (s["label"], s["el"], s.get("orb")):
            if k is not None and k in mapping:
                return k
        return None

    if explicit or cutoff == 0:
        shown = [i for i, v in enumerate(visible) if v]
    else:
        shown = [i for i, s in enumerate(series)
                 if visible[i] and (s["peak"] >= cutoff * top or s["label"] in mapping)]
        shown.sort(key=lambda i: -series[i]["weight"])

    # Palette colours a custom mapping already uses are skipped, so an
    # automatic colour never duplicates a hand-picked one.
    used_keys = {custom(series[i]) for i in shown} - {None}
    claimed = {to_hex(mapping[k]).lower() for k in used_keys}
    free = [c for c in colours if to_hex(c).lower() not in claimed] or colours

    styles = [{"color": MINOR_COLOR, "ls": "-", "legend": False, "rank": None}
              for _ in series]
    slot_of_key = {}
    for rank, i in enumerate(shown):
        s = series[i]
        k = custom(s)
        if k is not None:
            styles[i].update(color=mapping[k], ls="-", legend=True, rank=rank)
            continue
        group = s.get("key", s["label"])
        if group not in slot_of_key:
            slot_of_key[group] = len(slot_of_key)
        slot = slot_of_key[group]
        styles[i].update(color=free[slot % len(free)],
                         ls=LINESTYLES[(slot // len(free)) % len(LINESTYLES)],
                         legend=True, rank=rank)

    if len(slot_of_key) > len(free):
        warn(f"Note: {len(slot_of_key)} series share {len(free)} palette colours; "
             "repeats are told apart by line style. Use -e to pick fewer, or --panels.")

    known = ({s["label"] for s in series} | {s["el"] for s in series}
             | {s["orb"] for s in series if s.get("orb")})
    typos = sorted(set(mapping) - known)
    hidden = sorted((set(mapping) & known) - used_keys)
    if typos:
        warn(f"Warning: no plotted series matches colour key(s) {', '.join(typos)}. "
             f"Available: {', '.join(sorted(known))}")
    if hidden:
        warn(f"Note: colour key(s) {', '.join(hidden)} only match series below the "
             "legend cutoff, which are drawn grey. Name them with -e to show them.")
    return styles


def series_stats(energies, y_up, y_down, mask):
    """Peak and weight of one series inside the energy window `mask`."""
    import numpy as np
    up = np.abs(np.asarray(y_up)[mask])
    dn = np.abs(np.asarray(y_down)[mask]) if y_down is not None else np.zeros(0)
    peak = float(max(up.max(initial=0.0), dn.max(initial=0.0)))
    return peak, float(up.sum() + dn.sum())
