"""KPOINTS generation."""

import os
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Kpoints
from valyte.potcar import generate_potcar


def _resolve_poscar(poscar_path):
    """Find the requested POSCAR or the first POSCAR-like file nearby."""
    if poscar_path and os.path.exists(poscar_path):
        return poscar_path

    if poscar_path in (None, "POSCAR"):
        files = sorted(f for f in os.listdir(".") if f.startswith("POSCAR"))
        if files:
            found = files[0]
            print(f"Found structure: {found}")
            return found

    raise FileNotFoundError(f"{poscar_path or 'POSCAR'} not found")


def _normalize_scheme(scheme):
    """Normalize common KPOINTS scheme names."""
    value = (scheme or "gamma").strip().lower()
    if value in ("mp", "monkhorst", "monkhorst-pack", "monkhorst_pack"):
        return "MP"
    if value in ("gamma", "g"):
        return "Gamma"
    raise ValueError("scheme must be 'gamma' or 'mp'")


def _build_kpoints(structure, kspacing, scheme):
    """Create a Kpoints object and grid from a reciprocal-space spacing."""
    if kspacing <= 0:
        print("Using Gamma-only grid (1 1 1)")
        return Kpoints.gamma_automatic((1, 1, 1)), (1, 1, 1)

    recip_lattice = structure.lattice.reciprocal_lattice
    b_lengths = recip_lattice.abc

    # KSPACING-style estimate: N = |b| / kspacing
    grid = tuple(max(1, int(length / kspacing + 0.5)) for length in b_lengths)

    if scheme == "MP":
        return Kpoints.monkhorst_automatic(grid), grid
    return Kpoints.gamma_automatic(grid), grid


def generate_kpoints(
    poscar_path="POSCAR",
    kspacing=0.04,
    scheme="gamma",
    output="KPOINTS",
    generate_potcar_if_missing=True,
):
    """Generate a KPOINTS file from a POSCAR without prompting."""
    poscar_path = _resolve_poscar(poscar_path)
    structure = Structure.from_file(poscar_path)
    scheme = _normalize_scheme(scheme)
    kpts, grid = _build_kpoints(structure, kspacing, scheme)

    print("\nSummary")
    print(f"   Structure: {structure.formula}")
    print(f"   Lattice:   a={structure.lattice.a:.2f}, b={structure.lattice.b:.2f}, c={structure.lattice.c:.2f} A")
    print(f"   K-Mesh:    {grid[0]} x {grid[1]} x {grid[2]} ({scheme})")

    kpts.write_file(output)
    print(f"\nGenerated {output}.")

    potcar_file = "POTCAR"
    if not generate_potcar_if_missing:
        return

    if not os.path.exists(potcar_file):
        try:
            print("\nPOTCAR not found. Generating default POTCAR (PBE)...")
            generate_potcar(poscar_path=poscar_path, functional="PBE", output=potcar_file)
        except Exception as e:
            print(f"Warning: could not generate POTCAR: {e}")
    else:
        print("POTCAR already exists; skipping generation.")


def generate_kpoints_interactive():
    """Interactively generate a KPOINTS file based on POSCAR and user input."""
    print("\nValyte K-Point Generator")

    try:
        poscar_path = _resolve_poscar("POSCAR")
    except Exception as e:
        print(f"Error: {e}")
        return

    print("\nSelect K-Mesh Scheme:")
    print("  1. Monkhorst-Pack")
    print("  2. Gamma (Default)")

    choice = input("   > ").strip()
    scheme = "mp" if choice == "1" else "gamma"

    print("\nEnter K-Spacing (units of 2π/Å):")
    print("   Typical values: 0.03 - 0.04")

    try:
        kspacing_str = input("   > ").strip()
        kspacing = float(kspacing_str) if kspacing_str else 0.04
    except ValueError:
        print("Error: invalid number.")
        return

    try:
        generate_kpoints(poscar_path=poscar_path, kspacing=kspacing, scheme=scheme)
    except Exception as e:
        print(f"Error: {e}")
