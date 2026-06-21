"""Extract the electronic bandgap from a VASP vasprun.xml."""

import os
from pymatgen.io.vasp import Vasprun


def get_bandgap(vasprun_path="."):
    """Return the bandgap in eV from a vasprun.xml file.

    Parameters
    ----------
    vasprun_path : str
        Path to vasprun.xml or a directory containing it.

    Returns
    -------
    float
        Band gap in eV.
    """
    if os.path.isdir(vasprun_path):
        vasprun_path = os.path.join(vasprun_path, "vasprun.xml")

    vr = Vasprun(vasprun_path, parse_dos=False, parse_eigen=True)
    bs = vr.get_band_structure()
    return bs.get_band_gap()["energy"]
