"""Band structure KPOINTS generation."""

import os
import json
import numpy as np
import seekpath
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
try:
    from importlib.resources import files as ilr_files
except ImportError:
    from importlib_resources import files as ilr_files

from valyte.potcar import generate_potcar


def generate_band_kpoints(poscar_path="POSCAR", npoints=40, output="KPOINTS", symprec=0.01, mode="bradcrack"):
    """Generate a line-mode KPOINTS file for band structure calculations."""

    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"{poscar_path} not found")

    mode = (mode or "bradcrack").lower()

    structure = Structure.from_file(poscar_path)

    if mode == "bradcrack":
        try:
            kpath = BradCrackKpath(structure, symprec=symprec)
            prim_std = kpath.prim
            path = kpath.path
            kpoints = kpath.kpoints

            standard_filename = "POSCAR_standard"
            prim_std.to(filename=standard_filename)
        except Exception as e:
            raise RuntimeError(f"Error generating Bradley-Cracknell path: {e}")
    else:
        try:
            if mode == "seekpath":
                mode = "hinuma"

            sga = SpacegroupAnalyzer(structure, symprec=symprec)
            prim_std = sga.get_primitive_standard_structure()
        except Exception as e:
            raise RuntimeError(f"Error during standardization: {e}")

        try:
            kpath = HighSymmKpath(prim_std, path_type=mode, symprec=symprec)

            standard_filename = "POSCAR_standard"
            prim_std.to(filename=standard_filename)

            path = kpath.kpath["path"]
            kpoints = kpath.kpath["kpoints"]
        except Exception as e:
            raise RuntimeError(f"Error generating K-path: {e}")

    try:
        with open(output, "w") as f:
            f.write("KPOINTS for Band Structure\n")
            f.write(f"{npoints}\n")
            f.write("Line-mode\n")
            f.write("Reciprocal\n")

            for subpath in path:
                for i in range(len(subpath) - 1):
                    start_label = subpath[i]
                    end_label = subpath[i + 1]

                    start_coords = kpoints[start_label]
                    end_coords = kpoints[end_label]

                    f.write(
                        f"{start_coords[0]:10.6f} {start_coords[1]:10.6f} {start_coords[2]:10.6f} ! {start_label}\n"
                    )
                    f.write(
                        f"{end_coords[0]:10.6f} {end_coords[1]:10.6f} {end_coords[2]:10.6f} ! {end_label}\n"
                    )
                    f.write("\n")

        print(f"Generated {output} ({' - '.join([' - '.join(seg) for seg in path])})")
        print(f"Generated {standard_filename} (Standardized Primitive Cell)")
        print("IMPORTANT: Use POSCAR_standard for the band calculation.")
    except Exception as e:
        raise RuntimeError(f"Error writing KPOINTS file: {e}")

    try:
        print("Generating default POTCAR (PBE)...")
        generate_potcar(poscar_path=poscar_path, functional="PBE", output="POTCAR")
    except Exception as e:
        print(f"Warning: could not generate POTCAR: {e}")
        print("Proceeding without POTCAR generation.")


class BradCrackKpath:
    """Bradley-Cracknell K-path generation using SeeK-path output."""

    def __init__(self, structure, symprec=0.01):
        self.structure = structure
        self.symprec = symprec

        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        self._spg_data = sga.get_symmetry_dataset()

        cell = (
            structure.lattice.matrix,
            structure.frac_coords,
            [s.specie.number for s in structure],
        )

        self._seek_data = seekpath.get_path(cell, symprec=symprec)

        prim_lattice = self._seek_data["primitive_lattice"]
        prim_pos = self._seek_data["primitive_positions"]
        prim_types = self._seek_data["primitive_types"]

        z_to_specie = {s.specie.number: s.specie for s in structure}
        prim_species = [z_to_specie[z] for z in prim_types]

        self.prim = Structure(prim_lattice, prim_species, prim_pos)

        conv_lattice = self._seek_data["conv_lattice"]
        conv_pos = self._seek_data["conv_positions"]
        conv_types = self._seek_data["conv_types"]
        conv_species = [z_to_specie[z] for z in conv_types]
        self.conv = Structure(conv_lattice, conv_species, conv_pos)

        self._get_bradcrack_path()

    def _get_bradcrack_path(self):
        a, b, c = self.conv.lattice.abc
        angles = self.conv.lattice.angles

        angles_r = [round(x, 3) for x in angles]
        unique_val = min(angles_r, key=angles_r.count)
        unique = angles_r.index(unique_val)

        spg_symbol = self._spg_data["international"]
        spg_number = self._spg_data["number"]

        lattice_type = self.get_lattice_type(spg_number)
        bravais = self._get_bravais_lattice(spg_symbol, lattice_type, a, b, c, unique)

        json_file = ilr_files("valyte.data").joinpath("bradcrack.json")
        with json_file.open("r") as f:
            data = json.load(f)

        if bravais not in data:
            raise ValueError(f"Bravais lattice code '{bravais}' not found in BradCrack data.")

        self.bradcrack_data = data[bravais]
        self.kpoints = self.bradcrack_data["kpoints"]
        self.path = self.bradcrack_data["path"]

    def get_lattice_type(self, number):
        if 1 <= number <= 2:
            return "triclinic"
        if 3 <= number <= 15:
            return "monoclinic"
        if 16 <= number <= 74:
            return "orthorhombic"
        if 75 <= number <= 142:
            return "tetragonal"
        if 143 <= number <= 167:
            if number in [146, 148, 155, 160, 161, 166, 167]:
                return "rhombohedral"
            return "trigonal"
        if 168 <= number <= 194:
            return "hexagonal"
        if 195 <= number <= 230:
            return "cubic"
        return "unknown"

    def _get_bravais_lattice(self, spg_symbol, lattice_type, a, b, c, unique):
        if lattice_type == "triclinic":
            return "triclinic"

        if lattice_type == "monoclinic":
            if "P" in spg_symbol:
                if unique == 0:
                    return "mon_p_a"
                if unique == 1:
                    return "mon_p_b"
                if unique == 2:
                    return "mon_p_c"
            if "C" in spg_symbol:
                if unique == 0:
                    return "mon_c_a"
                if unique == 1:
                    return "mon_c_b"
                if unique == 2:
                    return "mon_c_c"

        if lattice_type == "orthorhombic":
            if "P" in spg_symbol:
                return "orth_p"
            if "A" in spg_symbol or "C" in spg_symbol:
                if a > b:
                    return "orth_c_a"
                if b > a:
                    return "orth_c_b"
            if "F" in spg_symbol:
                inv_a2 = 1 / a**2
                inv_b2 = 1 / b**2
                inv_c2 = 1 / c**2
                if (inv_a2 < inv_b2 + inv_c2) and (inv_b2 < inv_c2 + inv_a2) and (inv_c2 < inv_a2 + inv_b2):
                    return "orth_f_1"
                if inv_c2 > inv_a2 + inv_b2:
                    return "orth_f_2"
                if inv_b2 > inv_a2 + inv_c2:
                    return "orth_f_3"
                if inv_a2 > inv_c2 + inv_b2:
                    return "orth_f_4"
            if "I" in spg_symbol:
                if a > b and a > c:
                    return "orth_i_a"
                if b > a and b > c:
                    return "orth_i_b"
                if c > a and c > b:
                    return "orth_i_c"

        if lattice_type == "tetragonal":
            if "P" in spg_symbol:
                return "tet_p"
            if "I" in spg_symbol:
                return "tet_i_a" if a > c else "tet_i_c"

        if lattice_type in ["trigonal", "hexagonal", "rhombohedral"]:
            if "R" in spg_symbol:
                return "trig_r_a" if a > np.sqrt(2) * c else "trig_r_c"
            if "P" in spg_symbol:
                if unique == 0:
                    return "trig_p_a"
                if unique == 2:
                    return "trig_p_c"

        if lattice_type == "cubic":
            if "P" in spg_symbol:
                return "cubic_p"
            if "I" in spg_symbol:
                return "cubic_i"
            if "F" in spg_symbol:
                return "cubic_f"

        return "unknown"
