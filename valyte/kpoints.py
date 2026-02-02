"""Interactive KPOINTS generation."""

import os
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Kpoints
from valyte.potcar import generate_potcar


def generate_kpoints_interactive():
    """Interactively generate a KPOINTS file based on POSCAR and user input."""
    print("\nValyte K-Point Generator")

    poscar_path = "POSCAR"
    if not os.path.exists(poscar_path):
        files = [f for f in os.listdir(".") if f.startswith("POSCAR")]
        if files:
            poscar_path = files[0]
            print(f"Found structure: {poscar_path}")
        else:
            print("Error: POSCAR file not found in current directory.")
            return

    try:
        structure = Structure.from_file(poscar_path)
    except Exception as e:
        print(f"Error reading structure: {e}")
        return

    print("\nSelect K-Mesh Scheme:")
    print("  1. Monkhorst-Pack")
    print("  2. Gamma (Default)")

    choice = input("   > ").strip()
    scheme = "MP" if choice == "1" else "Gamma"

    print("\nEnter K-Spacing (units of 2π/Å):")
    print("   Typical values: 0.03 - 0.04")

    try:
        kspacing_str = input("   > ").strip()
        kspacing = float(kspacing_str) if kspacing_str else 0.04
    except ValueError:
        print("Error: invalid number.")
        return

    if kspacing <= 0:
        print("Using Gamma-only grid (1 1 1)")
        kpts = Kpoints.gamma_automatic((1, 1, 1))
        grid = (1, 1, 1)
    else:
        recip_lattice = structure.lattice.reciprocal_lattice
        b_lengths = recip_lattice.abc

        # KSPACING-style estimate: N = |b| / kspacing
        grid = [max(1, int(l / kspacing + 0.5)) for l in b_lengths]

        if scheme == "MP":
            kpts = Kpoints.monkhorst_automatic(grid)
        else:
            kpts = Kpoints.gamma_automatic(grid)

    print("\nSummary")
    print(f"   Structure: {structure.formula}")
    print(f"   Lattice:   a={structure.lattice.a:.2f}, b={structure.lattice.b:.2f}, c={structure.lattice.c:.2f} Å")
    print(f"   K-Mesh:    {grid[0]} x {grid[1]} x {grid[2]} ({scheme})")

    output_file = "KPOINTS"
    kpts.write_file(output_file)
    print(f"\nGenerated {output_file}.")

    potcar_file = "POTCAR"
    if not os.path.exists(potcar_file):
        try:
            print("\nPOTCAR not found. Generating default POTCAR (PBE)...")
            generate_potcar(poscar_path=poscar_path, functional="PBE", output=potcar_file)
        except Exception as e:
            print(f"Warning: could not generate POTCAR: {e}")
    else:
        print("POTCAR already exists; skipping generation.")
