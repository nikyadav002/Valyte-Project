"""Supercell generation."""

import os
from pymatgen.core import Structure


def create_supercell(poscar_path="POSCAR", nx=1, ny=1, nz=1, output="POSCAR_supercell"):
    """Create a supercell from a POSCAR file."""
    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"{poscar_path} not found")

    structure = Structure.from_file(poscar_path)
    supercell = structure.copy()
    supercell.make_supercell([nx, ny, nz])
    supercell.to(filename=output, fmt="poscar")

    supercell_atoms = len(supercell)
    print(f"Supercell created: {output} ({supercell_atoms} atoms)")
