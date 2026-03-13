import os
import re
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

# Orbital index mapping in pymatgen PROCAR order:
# 0:s  1:py  2:pz  3:px  4:dxy  5:dyz  6:dz2  7:dxz  8:x2-y2  9-15:f
_ORBITAL_INDICES = {
    "s": [0],
    "p": [1, 2, 3],
    "d": [4, 5, 6, 7, 8],
    "f": [9, 10, 11, 12, 13, 14, 15],
}


def _parse_orb_spec(spec):
    """Parse orbital spec string → (element_or_None, orbital_or_None).

    Supported formats:
      's', 'p', 'd', 'f'        → (None, 'orbital')
      'Fe'                       → ('Fe', None)
      'Fe:d' or 'Fe(d)'         → ('Fe', 'd')
    """
    spec = spec.strip()
    m = re.match(r"^([A-Za-z][a-z]?)[\:\(]([spdf])\)?$", spec)
    if m:
        return m.group(1).capitalize(), m.group(2).lower()
    if spec.lower() in _ORBITAL_INDICES:
        return None, spec.lower()
    return spec.capitalize(), None


def _get_orbital_weights(bs, spec, structure):
    """Return summed orbital projection weights, shape (nkpts, nbands).

    Uses Spin.up projections (first available spin).
    """
    spin = Spin.up if Spin.up in bs.projections else list(bs.projections.keys())[0]
    proj = bs.projections[spin]          # (nkpts, nbands, nions, norbitals)
    nkpts, nbands, nions, norbitals = proj.shape

    element, orbital = _parse_orb_spec(spec)

    if element is not None:
        atom_indices = [i for i, s in enumerate(structure) if s.specie.symbol == element]
        if not atom_indices:
            return np.zeros((nkpts, nbands))
    else:
        atom_indices = list(range(nions))

    if orbital is not None:
        orb_indices = [idx for idx in _ORBITAL_INDICES[orbital] if idx < norbitals]
    else:
        orb_indices = list(range(norbitals))

    return proj[:, :, :, :][:, :, atom_indices, :][:, :, :, orb_indices].sum(axis=(2, 3))


def _draw_triangle_legend(ax, tricolors, tri_labels):
    """Draw a ternary color triangle as an inset in the lower-right corner."""
    color_arr = np.array([to_rgb(c) for c in tricolors])

    # Equilateral triangle vertices: top, bottom-left, bottom-right
    V = np.array([[0.5, np.sqrt(3) / 2], [0.0, 0.0], [1.0, 0.0]])

    # Sample points in barycentric coordinates
    n = 60
    xs, ys, cs = [], [], []
    for w0 in np.linspace(0, 1, n):
        for w1 in np.linspace(0, 1 - w0, n):
            w2 = 1.0 - w0 - w1
            if w2 < -1e-9:
                continue
            w2 = max(w2, 0.0)
            xs.append(w0 * V[0, 0] + w1 * V[1, 0] + w2 * V[2, 0])
            ys.append(w0 * V[0, 1] + w1 * V[1, 1] + w2 * V[2, 1])
            cs.append(np.clip(w0 * color_arr[0] + w1 * color_arr[1] + w2 * color_arr[2], 0, 1))

    ax_tri = ax.inset_axes([0.70, 0.01, 0.28, 0.28])
    ax_tri.scatter(xs, ys, c=cs, s=4, linewidths=0, rasterized=True)

    tri = mpatches.Polygon(V, fill=False, edgecolor="k", linewidth=0.8)
    ax_tri.add_patch(tri)

    ax_tri.set_xlim(-0.18, 1.18)
    ax_tri.set_ylim(-0.18, np.sqrt(3) / 2 + 0.12)
    ax_tri.set_aspect("equal")
    ax_tri.axis("off")

    pad = 0.10
    ax_tri.text(V[0, 0], V[0, 1] + pad, tri_labels[0],
                ha="center", va="bottom", fontsize=7, fontweight="bold", color=tricolors[0])
    ax_tri.text(V[1, 0] - pad, V[1, 1] - pad * 0.6, tri_labels[1],
                ha="right", va="top", fontsize=7, fontweight="bold", color=tricolors[1])
    ax_tri.text(V[2, 0] + pad, V[2, 1] - pad * 0.6, tri_labels[2],
                ha="left", va="top", fontsize=7, fontweight="bold", color=tricolors[2])


def _save_band_dat(distances, energies, ticks, filepath):
    """Save band structure data (k-distances and energies) to a text file."""
    nbranches = len(distances)
    all_d = np.concatenate(distances)

    if isinstance(energies, dict):
        spins = list(energies.keys())
        nbands = len(energies[spins[0]][0])
        all_e = {}
        for spin in spins:
            all_e[spin] = np.array([
                np.concatenate([np.array(energies[spin][i][ib]) for i in range(nbranches)])
                for ib in range(nbands)
            ])
    else:
        spins = list(energies[0].keys())
        nbands = len(energies[0][spins[0]])
        all_e = {}
        for spin in spins:
            all_e[spin] = np.array([
                np.concatenate([np.array(energies[i][spin][ib]) for i in range(nbranches)])
                for ib in range(nbands)
            ])

    label_info = ", ".join(
        f"{l}={d:.4f}" for l, d in zip(ticks["label"], ticks["distance"]) if l
    )
    spin_tag = {spins[0]: ""} if len(spins) == 1 else {spins[0]: "_up", spins[1]: "_dn"}
    col_labels = ["k-dist"] + [
        f"band_{ib+1}{spin_tag[sp]}"
        for sp in spins for ib in range(nbands)
    ]
    header = f"K-point labels: {label_info}\n" + "  ".join(col_labels)

    cols = [all_d] + [all_e[sp][ib] for sp in spins for ib in range(nbands)]
    np.savetxt(filepath, np.column_stack(cols), header=header, fmt="%.6f")
    print(f"Saved data: {filepath}")


def plot_band_structure(vasprun_path, kpoints_path=None, output="valyte_band.png",
                        ylim=None, figsize=(4, 4), dpi=400, font="Arial",
                        save_data=False, spin_resolved=False):
    """Plot the electronic band structure from a VASP vasprun.xml."""

    if os.path.isdir(vasprun_path):
        vasprun_path = os.path.join(vasprun_path, "vasprun.xml")

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

    try:
        vr = BSVasprun(vasprun_path, parse_projected_eigen=False)
        bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)
    except Exception as e:
        raise ValueError(f"Failed to load band structure: {e}")

    bs_plotter = BSPlotter(bs)
    data = bs_plotter.bs_plot_data(zero_to_efermi=True)

    distances = data["distances"]
    energies = data["energy"]
    ticks = data["ticks"]

    fig, ax = plt.subplots(figsize=figsize)

    color_vb = "#8e44ad"
    color_cb = "#2a9d8f"

    # Determine spin keys in consistent order
    if isinstance(energies, dict):
        spin_list = list(energies.keys())
    else:
        spin_list = list(energies[0].keys())

    # Validate --spin-resolved: need at least two spin channels
    if spin_resolved and len(spin_list) < 2:
        print("Warning: Only one spin channel found in this calculation; "
              "--spin-resolved has no effect.")
        spin_resolved = False

    spin_styles = {}
    if spin_resolved:
        spin_styles = {
            spin_list[0]: {"color": "#3498db", "ls": "-"},   # spin up
            spin_list[1]: {"color": "#e74c3c", "ls": "--"},  # spin down
        }

    for i in range(len(distances)):
        d = distances[i]

        if isinstance(energies, dict):
            for spin in energies:
                for band in energies[spin][i]:
                    if spin_resolved:
                        style = spin_styles[spin]
                        ax.plot(d, band, color=style["color"], ls=style["ls"],
                                lw=1.5, alpha=1.0)
                    else:
                        c = color_vb if np.mean(band) <= 0 else color_cb
                        ax.plot(d, band, color=c, lw=1.5, alpha=1.0)
        else:
            for spin in energies[i]:
                for band in energies[i][spin]:
                    if spin_resolved:
                        style = spin_styles[spin]
                        ax.plot(d, band, color=style["color"], ls=style["ls"],
                                lw=1.5, alpha=1.0)
                    else:
                        c = color_vb if np.mean(band) <= 0 else color_cb
                        ax.plot(d, band, color=c, lw=1.5, alpha=1.0)

    if spin_resolved:
        up_line = mlines.Line2D([], [], color="#3498db", lw=1.5, ls="-", label="Spin up")
        dn_line = mlines.Line2D([], [], color="#e74c3c", lw=1.5, ls="--", label="Spin down")
        ax.legend(handles=[up_line, dn_line], fontsize=10, frameon=True,
                  loc="upper right")

    ax.set_xticks(ticks["distance"])
    clean_labels = [(l or "").replace("$\\mid$", "|") for l in ticks["label"]]
    ax.set_xticklabels(clean_labels, fontsize=14, fontweight="bold")

    for d in ticks["distance"]:
        ax.axvline(d, color="k", lw=0.8, ls="-", alpha=0.3)

    ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)

    ax.set_ylabel("Energy (eV)", fontsize=16, fontweight="bold", labelpad=8)
    if ylim:
        ax.set_ylim(ylim)
        yticks = np.arange(np.ceil(ylim[0]), np.floor(ylim[1]) + 1, 1)
        ax.set_yticks(yticks)
    else:
        ax.set_ylim(-4, 4)
        ax.set_yticks(np.arange(-4, 5, 1))

    ax.set_xlim(distances[0][0], distances[-1][-1])

    plt.tight_layout()
    plt.savefig(output, dpi=dpi)
    plt.close(fig)

    if save_data:
        _save_band_dat(distances, energies, ticks, "valyte_band.dat")


def plot_orbital_band_structure(
    vasprun_path,
    kpoints_path=None,
    output="valyte_band_orbital.png",
    ylim=None,
    tricolor=None,
    tricolors=None,
    tri_labels=None,
    figsize=(4, 4),
    dpi=400,
    lw=2.0,
    font="Arial",
    save_data=False,
):
    """Plot orbital-resolved band structure with a tricolor RGB map.

    Each band segment is colored by blending three base colors weighted by the
    relative orbital/element projections at that k-point.

    Parameters
    ----------
    tricolor : list of 3 str
        Orbital/element specs, e.g. ['s', 'p', 'd'] or ['Fe:d', 'O:p', 's'].
        Formats accepted: 's'|'p'|'d'|'f', 'Element', 'Element:orb', 'Element(orb)'.
    tricolors : list of 3 str
        Matplotlib color strings for each spec. Default: red, green, blue.
    tri_labels : list of 3 str
        Labels shown in the triangle legend. Defaults to tricolor specs.
    """
    if os.path.isdir(vasprun_path):
        vasprun_path = os.path.join(vasprun_path, "vasprun.xml")

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

    if tricolor is None:
        tricolor = ["s", "p", "d"]
    if tricolors is None:
        tricolors = ["#e74c3c", "#2ecc71", "#3498db"]
    if tri_labels is None:
        tri_labels = list(tricolor)

    if len(tricolor) != 3 or len(tricolors) != 3:
        raise ValueError("Exactly 3 orbital specs and 3 colors are required for tricolor mode.")

    color_arr = np.array([to_rgb(c) for c in tricolors])  # (3, 3)

    try:
        vr = BSVasprun(vasprun_path, parse_projected_eigen=True)
        bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)
    except Exception as e:
        raise ValueError(f"Failed to load band structure: {e}")

    structure = vr.final_structure

    if not bs.projections:
        raise ValueError(
            "No orbital projections found. Run VASP with LORBIT=11 (or ≥10) "
            "to generate projected eigenvalues."
        )

    bs_plotter = BSPlotter(bs)
    data = bs_plotter.bs_plot_data(zero_to_efermi=True)

    distances = data["distances"]
    energies = data["energy"]
    ticks = data["ticks"]

    # Determine how pymatgen keys the energy dict (Spin object or string)
    if isinstance(energies, dict):
        spin_keys_e = list(energies.keys())
        def get_energy(branch_i, band_i, spin_i):
            return energies[spin_keys_e[spin_i]][branch_i][band_i]
    else:
        spin_keys_e = list(energies[0].keys())
        def get_energy(branch_i, band_i, spin_i):
            return energies[branch_i][spin_keys_e[spin_i]][band_i]

    # Compute raw projection weights for each spec → (nkpts, nbands)
    raw_weights = [_get_orbital_weights(bs, spec, structure) for spec in tricolor]
    total = sum(raw_weights)
    total = np.where(total < 1e-9, 1e-9, total)
    norm_weights = [w / total for w in raw_weights]  # each (nkpts, nbands)

    fig, ax = plt.subplots(figsize=figsize)

    spins_proj = list(bs.projections.keys())

    for branch_i, branch in enumerate(bs.branches):
        kpt_start = branch["start_index"]
        kpt_end = branch["end_index"] + 1      # exclusive
        d = distances[branch_i]

        for spin_i in range(len(spins_proj)):
            for ib in range(bs.nb_bands):
                e = np.array(get_energy(branch_i, ib, spin_i))

                # RGB color at each k-point in this branch
                nk = kpt_end - kpt_start
                rgb = np.zeros((nk, 3))
                for iw, nw in enumerate(norm_weights):
                    w_branch = nw[kpt_start:kpt_end, ib]
                    rgb += np.outer(w_branch, color_arr[iw])
                rgb = np.clip(rgb, 0, 1)

                # Build line segments; color = average of the two endpoint colors
                pts = np.array([d, e]).T.reshape(-1, 1, 2)
                segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
                seg_colors = (rgb[:-1] + rgb[1:]) / 2.0

                linestyle = "-" if spin_i == 0 else "--"
                lc = LineCollection(segs, colors=seg_colors, linewidths=lw,
                                    linestyle=linestyle, alpha=0.9)
                ax.add_collection(lc)

    ax.set_xticks(ticks["distance"])
    clean_labels = [(l or "").replace("$\\mid$", "|") for l in ticks["label"]]
    ax.set_xticklabels(clean_labels, fontsize=14, fontweight="bold")

    for d in ticks["distance"]:
        ax.axvline(d, color="k", lw=0.8, ls="-", alpha=0.3)

    ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)

    ax.set_ylabel("Energy (eV)", fontsize=16, fontweight="bold", labelpad=8)
    if ylim:
        ax.set_ylim(ylim)
        ax.set_yticks(np.arange(np.ceil(ylim[0]), np.floor(ylim[1]) + 1, 1))
    else:
        ax.set_ylim(-4, 4)
        ax.set_yticks(np.arange(-4, 5, 1))

    ax.set_xlim(distances[0][0], distances[-1][-1])

    _draw_triangle_legend(ax, tricolors, tri_labels)

    plt.tight_layout()
    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")

    if save_data:
        _save_band_dat(distances, energies, ticks, "valyte_band.dat")


def plot_spin_texture_band_structure(
    vasprun_path,
    kpoints_path=None,
    output="valyte_band_spin_texture.png",
    ylim=None,
    figsize=(4, 4),
    dpi=400,
    font="Arial",
    save_data=False,
    spin_component="sz",
    cmap="seismic",
    lw=1.5,
):
    """Plot non-collinear band structure colored by a spin texture component (Sx, Sy, or Sz).

    Each band segment is colored by the expectation value of the chosen spin
    component summed over all atoms and orbitals, using a diverging colormap
    centered at zero.

    Requires a VASP non-collinear calculation (LSORBIT = .TRUE. or
    LNONCOLLINEAR = .TRUE.) with LORBIT >= 11.

    Parameters
    ----------
    spin_component : {'sx', 'sy', 'sz'}
        Spin magnetization component to visualize.
    cmap : str
        Diverging colormap name (default: 'seismic').
    lw : float
        Band line width.
    """
    if os.path.isdir(vasprun_path):
        vasprun_path = os.path.join(vasprun_path, "vasprun.xml")

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

    try:
        vr = BSVasprun(vasprun_path, parse_projected_eigen=True)
        bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)
    except Exception as e:
        raise ValueError(f"Failed to load band structure: {e}")

    # projected_magnetisation shape: (nkpoints, nbands, natoms, norbitals, 3)
    # last axis → [0]=Sx, [1]=Sy, [2]=Sz  (British spelling in pymatgen)
    mag = getattr(vr, "projected_magnetisation", None)
    if mag is None:
        raise ValueError(
            "No spin magnetization data found in vasprun.xml.\n"
            "Spin texture requires a non-collinear VASP calculation with "
            "LSORBIT = .TRUE. (or LNONCOLLINEAR = .TRUE.) and LORBIT >= 11."
        )

    comp_idx = {"sx": 0, "sy": 1, "sz": 2}[spin_component.lower()]
    mag = np.asarray(mag)  # (nkpts, nbands, natoms, norbitals, 3)
    if mag.ndim != 5 or mag.shape[-1] != 3:
        raise ValueError(
            f"Unexpected projected_magnetisation shape {mag.shape}. "
            "Expected (nkpts, nbands, natoms, norbitals, 3)."
        )
    # Sum over atoms and orbitals → (nkpts, nbands)
    spin_texture = mag[:, :, :, :, comp_idx].sum(axis=(2, 3))

    bs_plotter = BSPlotter(bs)
    data = bs_plotter.bs_plot_data(zero_to_efermi=True)

    distances = data["distances"]
    energies = data["energy"]
    ticks = data["ticks"]

    # Energy accessor — non-collinear has a single spin key
    if isinstance(energies, dict):
        spin_key = list(energies.keys())[0]
        def get_energy(branch_i, band_i):
            return energies[spin_key][branch_i][band_i]
    else:
        spin_key = list(energies[0].keys())[0]
        def get_energy(branch_i, band_i):
            return energies[branch_i][spin_key][band_i]

    # Symmetric color limits centered at zero
    clim = max(abs(float(spin_texture.max())), abs(float(spin_texture.min())))
    if clim < 1e-9:
        clim = 1.0
    norm = mpl.colors.Normalize(vmin=-clim, vmax=clim)

    fig, ax = plt.subplots(figsize=figsize)

    for branch_i, branch in enumerate(bs.branches):
        kpt_start = branch["start_index"]
        kpt_end = branch["end_index"] + 1
        d = distances[branch_i]

        for ib in range(bs.nb_bands):
            e = np.array(get_energy(branch_i, ib))
            st = spin_texture[kpt_start:kpt_end, ib]

            pts = np.array([d, e]).T.reshape(-1, 1, 2)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            seg_values = (st[:-1] + st[1:]) / 2.0

            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=lw, alpha=0.9)
            lc.set_array(seg_values)
            ax.add_collection(lc)

    comp_label = {"sx": r"$S_x$", "sy": r"$S_y$", "sz": r"$S_z$"}[spin_component.lower()]
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label(comp_label, fontsize=13, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    ax.set_xticks(ticks["distance"])
    clean_labels = [(l or "").replace("$\\mid$", "|") for l in ticks["label"]]
    ax.set_xticklabels(clean_labels, fontsize=14, fontweight="bold")

    for d in ticks["distance"]:
        ax.axvline(d, color="k", lw=0.8, ls="-", alpha=0.3)

    ax.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)

    ax.set_ylabel("Energy (eV)", fontsize=16, fontweight="bold", labelpad=8)
    if ylim:
        ax.set_ylim(ylim)
        ax.set_yticks(np.arange(np.ceil(ylim[0]), np.floor(ylim[1]) + 1, 1))
    else:
        ax.set_ylim(-4, 4)
        ax.set_yticks(np.arange(-4, 5, 1))

    ax.set_xlim(distances[0][0], distances[-1][-1])

    plt.tight_layout()
    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")

    if save_data:
        _save_band_dat(distances, energies, ticks, "valyte_band.dat")
