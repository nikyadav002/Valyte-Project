"""First Brillouin zone with the suggested high-symmetry k-path."""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from valyte.style import apply_style, get_font_weight, save_plot, font_scale

BZ_EDGE = "#2b2d42"
BZ_FACE = "#8d99ae"
PATH_COLOR = "#e63946"
POINT_COLOR = "#1d3557"
VECTOR_COLOR = "#2a9d8f"

_GREEK = ("\\Gamma", "\\Sigma", "\\Delta", "\\Lambda")


def _format_label(label, bold=True):
    """Turn a raw path label such as '\\\\Gamma' into mathtext."""
    text = (label or "").strip()
    if not text:
        return ""
    if bold:
        for g in _GREEK:
            text = text.replace(g, f"\\mathbf{{{g}}}")
        if "\\" not in text:
            text = f"\\mathbf{{{text}}}"
    return f"${text}$"


def _cartesian(lattice, frac):
    """Fractional reciprocal coordinates to Cartesian (2*pi convention)."""
    return np.asarray(lattice.reciprocal_lattice.get_cartesian_coords(frac))


def plot_brillouin_zone(prim_std, path, kpoints, output="valyte_bz.png",
                        figsize=(5.5, 5.5), dpi=400, font="Arial", bold=True,
                        fontsize=None, elev=22, azim=30):
    """Draw the first Brillouin zone with the k-path overlaid.

    Parameters
    ----------
    prim_std : pymatgen Structure
        The cell the path is defined in, as returned by band.resolve_kpath.
    path : list of list of str
        Branches of point labels; a break between branches is a discontinuity.
    kpoints : dict
        {label: [kx, ky, kz]} in fractional reciprocal coordinates.
    elev, azim : float
        Viewing angles in degrees. The figure is static, so these are the only
        way to change the viewpoint.
    """
    _weight = get_font_weight(bold)
    _fsbase = fontsize if fontsize is not None else 12
    _fscale = font_scale(fontsize, 12)
    apply_style(font=font, bold=bold, fontsize=_fsbase,
                linewidth=1.4 if bold else 0.8)

    lattice = prim_std.lattice

    # get_brillouin_zone() takes the reciprocal lattice internally, so it must
    # be called on the real-space lattice.  Calling it on lattice.reciprocal_
    # lattice silently returns the Wigner-Seitz cell of the wrong lattice.
    facets = lattice.get_brillouin_zone()

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    # ── Brillouin zone ────────────────────────────────────────────────────
    poly = Poly3DCollection(
        [np.asarray(f) for f in facets],
        facecolor=BZ_FACE, edgecolor=BZ_EDGE,
        linewidths=0.7 if bold else 0.5, alpha=0.10,
    )
    poly.set_zorder(1)
    ax.add_collection3d(poly)

    # ── k-path ────────────────────────────────────────────────────────────
    drawn = set()
    for subpath in path:
        for i in range(len(subpath) - 1):
            a, b = subpath[i], subpath[i + 1]
            if a not in kpoints or b not in kpoints:
                continue
            p0, p1 = _cartesian(lattice, kpoints[a]), _cartesian(lattice, kpoints[b])
            ax.plot(*zip(p0, p1), color=PATH_COLOR,
                    lw=3.0 if bold else 2.2, solid_capstyle="round", zorder=5)
            drawn.update((a, b))

    # ── high-symmetry points and labels ───────────────────────────────────
    # Some conventions place points outside the first zone for low-symmetry
    # lattices, so the framing has to account for the k-points too or they
    # get clipped.  They are drawn where the convention puts them, unfolded.
    bz_span = max(np.linalg.norm(v) for f in facets for v in f)
    pts = [_cartesian(lattice, kpoints[l]) for l in drawn] or [np.zeros(3)]
    span = max(bz_span, max(np.linalg.norm(p) for p in pts))
    offset = 0.11 * span
    for label in drawn:
        p = _cartesian(lattice, kpoints[label])
        ax.scatter(*p, color=POINT_COLOR, s=45 * _fscale ** 2,
                   edgecolor="white", linewidth=0.8,
                   depthshade=False, zorder=6)
        # push the label radially outward so it clears the path and the zone
        direction = p / np.linalg.norm(p) if np.linalg.norm(p) > 1e-9 else np.array([0, 0, 1.0])
        lp = p + direction * offset
        ax.text(lp[0], lp[1], lp[2], _format_label(label, bold),
                color=POINT_COLOR, fontsize=13 * _fscale,
                ha="center", va="center", zorder=7)

    # ── reciprocal lattice vectors ────────────────────────────────────────
    bvecs = np.asarray(lattice.reciprocal_lattice.matrix)
    vscale = 0.62 * bz_span / max(np.linalg.norm(v) for v in bvecs)
    for vec, name in zip(bvecs, ("b_1", "b_2", "b_3")):
        v = np.asarray(vec) * vscale
        ax.quiver(0, 0, 0, *v, color=VECTOR_COLOR, arrow_length_ratio=0.15,
                  lw=1.2 if bold else 0.9, alpha=0.9, zorder=3)
        tip = v * 1.22
        ax.text(*tip, f"${name}$", color=VECTOR_COLOR,
                fontsize=10 * _fscale, ha="center", va="center", zorder=7)

    # ── framing ───────────────────────────────────────────────────────────
    # Equal box aspect with symmetric limits, or the polyhedron renders skewed.
    lim = 1.05 * span
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect([1, 1, 1])
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)

    save_plot(fig, output, dpi=dpi)
