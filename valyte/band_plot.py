import os
import numpy as np
import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotter


def plot_band_structure(vasprun_path, kpoints_path=None, output="valyte_band.png",
                        ylim=None, figsize=(4, 4), dpi=400, font="Arial"):
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

    for i in range(len(distances)):
        d = distances[i]

        if isinstance(energies, dict):
            for spin in energies:
                for band in energies[spin][i]:
                    c = color_vb if np.mean(band) <= 0 else color_cb
                    ax.plot(d, band, color=c, lw=1.5, alpha=1.0)
        else:
            for spin in energies[i]:
                for band in energies[i][spin]:
                    c = color_vb if np.mean(band) <= 0 else color_cb
                    ax.plot(d, band, color=c, lw=1.5, alpha=1.0)

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
