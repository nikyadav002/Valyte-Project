import os
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Potcar

DEFAULT_FUNCTIONAL = "PBE_64"

# VASP-recommended PAW potentials (bold entries of the potpaw.64 PBE standard
# table, https://www.vasp.at/wiki/index.php/Available_pseudopotentials).
# Elements not listed here use the plain element symbol.
RECOMMENDED_POTCARS = {
    "Li": "Li_sv", "Na": "Na_pv", "K": "K_sv", "Ca": "Ca_sv",
    "Rb": "Rb_sv", "Sr": "Sr_sv", "Cs": "Cs_sv", "Ba": "Ba_sv",
    "Fr": "Fr_sv", "Ra": "Ra_sv",
    "Sc": "Sc_sv", "Ti": "Ti_sv", "V": "V_sv", "Cr": "Cr_pv", "Mn": "Mn_pv",
    "Y": "Y_sv", "Zr": "Zr_sv", "Nb": "Nb_sv", "Mo": "Mo_sv",
    "Tc": "Tc_pv", "Ru": "Ru_pv", "Rh": "Rh_pv",
    "Hf": "Hf_pv", "Ta": "Ta_pv", "W": "W_sv",
    "Ga": "Ga_d", "Ge": "Ge_d", "In": "In_d", "Sn": "Sn_d",
    "Tl": "Tl_d", "Pb": "Pb_d", "Bi": "Bi_d", "Po": "Po_d",
    "Pr": "Pr_3", "Nd": "Nd_3", "Pm": "Pm_3", "Sm": "Sm_3", "Eu": "Eu_2",
    "Gd": "Gd_3", "Tb": "Tb_3", "Dy": "Dy_3", "Ho": "Ho_3", "Er": "Er_3",
    "Tm": "Tm_3", "Yb": "Yb_2", "Lu": "Lu_3",
}


def _ordered_symbols(structure):
    symbols = []
    seen = set()
    for site in structure:
        sym = site.specie.symbol if hasattr(site.specie, "symbol") else str(site.specie)
        if sym not in seen:
            symbols.append(sym)
            seen.add(sym)
    return symbols


def parse_potcar_overrides(items):
    """Parse ["Fe=Fe_pv", "O=O_h"] into {"Fe": "Fe_pv", "O": "O_h"}."""
    overrides = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Invalid POTCAR override '{item}'; expected Element=POTCAR, e.g. Fe=Fe_pv")
        el, name = item.split("=", 1)
        overrides[el.strip()] = name.strip()
    return overrides


def recommended_potcar_symbols(species, overrides=None):
    """Map element symbols to POTCAR names: overrides, then VASP recommendation."""
    overrides = overrides or {}
    return [overrides.get(el) or RECOMMENDED_POTCARS.get(el, el) for el in species]


def generate_potcar(poscar_path="POSCAR", functional=DEFAULT_FUNCTIONAL, output="POTCAR",
                    overrides=None):
    """Generate a POTCAR from the species in a POSCAR using VASP-recommended potentials."""
    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"Input file '{poscar_path}' not found.")

    structure = Structure.from_file(poscar_path)
    species = _ordered_symbols(structure)
    symbols = recommended_potcar_symbols(species, overrides)

    print(f"Reading structure from {poscar_path}...")
    print(f"Detected species: {species}")
    print(f"POTCARs: {', '.join(f'{el} -> {sym}' for el, sym in zip(species, symbols))}")

    try:
        potcar = Potcar(symbols, functional=functional)
    except OSError as e:
        print(f"Could not find POTCARs for functional '{functional}'.")
        print("Ensure PMG_VASP_PSP_DIR is set in ~/.pmgrc.yaml")
        raise e

    potcar.write_file(output)
    print(f"Generated {output} using functional '{functional}'")
