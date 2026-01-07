"""
Band structure KPOINTS generation module for Valyte.
"""

import os
import json
import numpy as np
import seekpath
import spglib
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
try:
    from importlib.resources import files as ilr_files
except ImportError:
    import importlib_resources as ilr_files


def generate_band_kpoints(poscar_path="POSCAR", npoints=40, output="KPOINTS", symprec=0.01, mode="bradcrack"):
    """
    Generates KPOINTS file in line-mode for band structure calculations.
    Uses SeeK-path method for high-symmetry path determination.
    
    IMPORTANT: Writes a standardized POSCAR (POSCAR_standard) that MUST be used
    for the band structure calculation to ensure K-points are valid.
    
    Args:
        poscar_path (str): Path to input POSCAR file.
        npoints (int): Number of points per segment (default: 40).
        output (str): Output filename for KPOINTS.
        symprec (float): Symmetry precision for standardization (default: 0.01).
        mode (str): Standardization convention (default: "bradcrack").
    """
    
    if not os.path.exists(poscar_path):
        raise FileNotFoundError(f"{poscar_path} not found")
    
    # Read structure
    structure = Structure.from_file(poscar_path)
    
    # --- K-Point Generation Logic ---
    if mode == "bradcrack":
        try:
            kpath = BradCrackKpath(structure, symprec=symprec)
            prim_std = kpath.prim
            path = kpath.path
            kpoints = kpath.kpoints
            
            # Write standardized POSCAR from BradCrack logic
            standard_filename = "POSCAR_standard"
            prim_std.to(filename=standard_filename)
        except Exception as e:
            print(f"❌ Error generating BradCrack path: {e}")
            return

    else: 
        # Fallback to Pymatgen logic for other modes
        try:
            # Map 'seekpath' alias to 'hinuma' which pymatgen uses (wrapper around seekpath)
            if mode == "seekpath":
                mode = "hinuma"

            # Standardize structure first using SpacegroupAnalyzer
            sga = SpacegroupAnalyzer(structure, symprec=symprec)
            prim_std = sga.get_primitive_standard_structure()
        except Exception as e:
            print(f"❌ Error during standardization: {e}")
            return

        # Get high-symmetry path for the STANDARDIZED structure
        try:
            kpath = HighSymmKpath(prim_std, path_type=mode, symprec=symprec)
            
            # Write the standardized primitive structure
            standard_filename = "POSCAR_standard"
            prim_std.to(filename=standard_filename)
            
            # Get the path
            path = kpath.kpath["path"]
            kpoints = kpath.kpath["kpoints"]
        except Exception as e:
            print(f"❌ Error generating K-path: {e}")
            return

    # Write KPOINTS file
    try:
        with open(output, "w") as f:
            f.write("KPOINTS for Band Structure\n")
            f.write(f"{npoints}\n")
            f.write("Line-mode\n")
            f.write("Reciprocal\n")
            
            for subpath in path:
                for i in range(len(subpath) - 1):
                    start_label = subpath[i]
                    end_label = subpath[i+1]
                    
                    start_coords = kpoints[start_label]
                    end_coords = kpoints[end_label]
                    
                    f.write(f"{start_coords[0]:10.6f} {start_coords[1]:10.6f} {start_coords[2]:10.6f} ! {start_label}\n")
                    f.write(f"{end_coords[0]:10.6f} {end_coords[1]:10.6f} {end_coords[2]:10.6f} ! {end_label}\n")
                    f.write("\n") # Optional newline between segments

        print(f"✅ Generated {output} ({' - '.join([' - '.join(seg) for seg in path])})")
        print(f"✅ Generated {standard_filename} (Standardized Primitive Cell)")
        print(f"\n⚠️  IMPORTANT: You MUST use '{standard_filename}' for your band calculation!")
        print(f"   The K-points are generated for this standardized orientation.")
        print(f"   Using your original POSCAR may result in incorrect paths or 'Reciprocal lattice' errors.")

    except Exception as e:
        print(f"❌ Error writing KPOINTS file: {e}")


class BradCrackKpath:
    """
    Native implementation of Bradley-Cracknell K-path generation.
    Replicates logic from Sumo/SeeK-path to determine standard paths.
    """
    def __init__(self, structure, symprec=0.01):
        self.structure = structure
        self.symprec = symprec
        
        # Use SpacegroupAnalyzer for basic data
        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        self._spg_data = sga.get_symmetry_dataset()
        
        # Use SeeK-path to get primitive/conventional structures matches Sumo Kpath.__init__
        
        # refine_cell logic from Sumo base class
        # atom_numbers = [site.specie.number for site in structure] 
        # But pymatgen structure to spglib cell tuple:
        # cell = (lattice, positions, numbers)
        cell = (structure.lattice.matrix, structure.frac_coords, [s.specie.number for s in structure])
        
        # Sumo uses spglib.refine_cell on the cell first? 
        # "std = spglib.refine_cell(sym._cell, symprec=symprec)" 
        # pymatgen sga._cell is (lattice, positions, numbers)
        
        # seekpath.get_path takes the cell structure
        # output is dictionary
        self._seek_data = seekpath.get_path(cell, symprec=symprec)
        
        # Reconstruct primitive structure from seekpath output
        prim_lattice = self._seek_data["primitive_lattice"]
        prim_pos = self._seek_data["primitive_positions"]
        prim_types = self._seek_data["primitive_types"]
        # Map types back to species? 
        # We need a map from number to Element. 
        # unique_species from sga?
        # Let's just use explicit element list from input structure, assuming types are consistent?
        # Or better, use sga to map Z to elements.
        
        # Setup element mapping
        # Create a map from atomic number to Element object from input structure
        z_to_specie = {s.specie.number: s.specie for s in structure}
        prim_species = [z_to_specie[z] for z in prim_types]
        
        self.prim = Structure(prim_lattice, prim_species, prim_pos)
        
        conv_lattice = self._seek_data["conv_lattice"]
        conv_pos = self._seek_data["conv_positions"]
        conv_types = self._seek_data["conv_types"]
        conv_species = [z_to_specie[z] for z in conv_types]
        self.conv = Structure(conv_lattice, conv_species, conv_pos)
        
        # Now determine Bravais lattice for BradCrack
        self._get_bradcrack_path()

    def _get_bradcrack_path(self):
        
        # Determine lattice parameters from CONVENTIONAL cell
        a, b, c = self.conv.lattice.abc
        angles = self.conv.lattice.angles
        # finding unique axis for monoclinic
        # logic from BradCrackKpath.__init__
        # "unique = angles.index(min(angles, key=angles.count))"
        # usually 90, 90, beta. So unique is beta (non-90) index? No.
        # Monoclinic: alpha=gamma=90, beta!=90. 90 appears twice. non-90 appears once.
        # min count of angle values?
        # if angles are [90, 90, 105], counts are {90:2, 105:1}. min count is 1. value is 105. index is 2.
        # so unique is index of non-90 degree angle.
        
        # Round angles to avoid float issues
        angles_r = [round(x, 3) for x in angles]
        unique_val = min(angles_r, key=angles_r.count)
        unique = angles_r.index(unique_val)

        # Get Space Group Symbol and Number
        # From seekpath? or sga? 
        # Sumo uses: "spg_symbol = self.spg_symbol" which is "self._spg_data['international']"
        # spglib dataset returns 'international'
        spg_symbol = self._spg_data["international"]
        spg_number = self._spg_data["number"]
        
        lattice_type = self.get_lattice_type(spg_number)
        
        bravais = self._get_bravais_lattice(spg_symbol, lattice_type, a, b, c, unique)
        
        # Load JSON
        
        json_file = ilr_files("valyte.data").joinpath("bradcrack.json")
        with open(json_file, 'r') as f:
            data = json.load(f)
            
        if bravais not in data:
            raise ValueError(f"Bravais lattice code '{bravais}' not found in BradCrack data.")
            
        self.bradcrack_data = data[bravais]
        self.kpoints = self.bradcrack_data["kpoints"]
        self.path = self.bradcrack_data["path"]

    def get_lattice_type(self, number):
        # Logic from Sumo
        if 1 <= number <= 2: return "triclinic"
        if 3 <= number <= 15: return "monoclinic"
        if 16 <= number <= 74: return "orthorhombic"
        if 75 <= number <= 142: return "tetragonal"
        if 143 <= number <= 167: 
            if number in [146, 148, 155, 160, 161, 166, 167]: return "rhombohedral"
            return "trigonal"
        if 168 <= number <= 194: return "hexagonal"
        if 195 <= number <= 230: return "cubic"
        return "unknown"

    def _get_bravais_lattice(self, spg_symbol, lattice_type, a, b, c, unique):
        # Logic from Sumo BradCrackKpath._get_bravais_lattice
        if lattice_type == "triclinic": return "triclinic"
        
        elif lattice_type == "monoclinic":
            if "P" in spg_symbol:
                if unique == 0: return "mon_p_a"
                elif unique == 1: return "mon_p_b"
                elif unique == 2: return "mon_p_c"
            elif "C" in spg_symbol:
                if unique == 0: return "mon_c_a"
                elif unique == 1: return "mon_c_b"
                elif unique == 2: return "mon_c_c"
        
        elif lattice_type == "orthorhombic":
            if "P" in spg_symbol: return "orth_p"
            elif "A" in spg_symbol or "C" in spg_symbol:
                if a > b: return "orth_c_a"
                elif b > a: return "orth_c_b"
            elif "F" in spg_symbol:
                # 1/a^2 etc conditions... need to replicate exact math
                # Copied from Sumo source view
                inv_a2 = 1/a**2; inv_b2 = 1/b**2; inv_c2 = 1/c**2
                if (inv_a2 < inv_b2 + inv_c2) and (inv_b2 < inv_c2 + inv_a2) and (inv_c2 < inv_a2 + inv_b2):
                    return "orth_f_1"
                elif inv_c2 > inv_a2 + inv_b2: return "orth_f_2"
                elif inv_b2 > inv_a2 + inv_c2: return "orth_f_3"
                elif inv_a2 > inv_c2 + inv_b2: return "orth_f_4"
            elif "I" in spg_symbol:
                if a > b and a > c: return "orth_i_a"
                elif b > a and b > c: return "orth_i_b"
                elif c > a and c > b: return "orth_i_c"

        elif lattice_type == "tetragonal":
            if "P" in spg_symbol: return "tet_p"
            elif "I" in spg_symbol:
                if a > c: return "tet_i_a"
                else: return "tet_i_c"

        elif lattice_type in ["trigonal", "hexagonal", "rhombohedral"]:
            if "R" in spg_symbol:
                if a > np.sqrt(2)*c: return "trig_r_a"
                else: return "trig_r_c"
            elif "P" in spg_symbol:
                if unique == 0: return "trig_p_a"
                elif unique == 2: return "trig_p_c"

        elif lattice_type == "cubic":
            if "P" in spg_symbol: return "cubic_p"
            elif "I" in spg_symbol: return "cubic_i"
            elif "F" in spg_symbol: return "cubic_f"
            
        return "unknown"
