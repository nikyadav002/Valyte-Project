import os
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Potcar


def _ordered_symbols(structure):
    symbols = []
    seen = set()
    for site in structure:
        sym = site.specie.symbol if hasattr(site.specie, "symbol") else str(site.specie)
        if sym not in seen:
            symbols.append(sym)
            seen.add(sym)
    return symbols


def generate_potcar(poscar_path="POSCAR", functional="PBE", output="POTCAR"):
    """Generate a POTCAR based on species in a POSCAR."""
    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"Input file '{poscar_path}' not found.")

    structure = Structure.from_file(poscar_path)
    species = _ordered_symbols(structure)

    print(f"Reading structure from {poscar_path}...")
    print(f"Detected species: {species}")

    try:
        potcar = Potcar(species, functional=functional)
    except OSError as e:
        print(f"Could not find POTCARs for functional '{functional}'.")
        print("Ensure PMG_VASP_PSP_DIR is set in ~/.pmgrc.yaml")
        raise e

    potcar.write_file(output)
    print(f"Generated {output} using functional '{functional}'")
