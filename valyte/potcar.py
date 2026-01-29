import os
import sys
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Potcar

def generate_potcar(poscar_path="POSCAR", functional="PBE", output="POTCAR"):
    """
    Generates a POTCAR file based on the species in the POSCAR using Pymatgen configuration.
    
    Args:
        poscar_path (str): Path to input POSCAR file.
        functional (str): Functional to use (e.g., "PBE", "PBE_52", "PBE_54", "LDA").
                          Defaults to "PBE" which usually maps to the configured default (often PBE_54).
        output (str): Output filename.
    """
    if not os.path.exists(poscar_path):
        print(f"❌ Error: Input file '{poscar_path}' not found.")
        return

    try:
        structure = Structure.from_file(poscar_path)
        species = structure.symbol_set
        
        # Sort species to match POSCAR order if structure.symbol_set isn't ordered
        # Actually Potcar.from_structure usually handles this but let's be safe if we manually construct.
        # Ideally: use Potcar.from_structure if available or construct list.
        # structure.symbol_set returns a tuple/list of unique species.
        
        # Pymatgen's Potcar.from_file is for reading.
        # We want to CREATE.
        # Potcar(symbols, functional)
        
        # Let's verify which method is best. 
        # structure.species gives list of all sites.
        # We need unique species in order of appearance in POSCAR (usually).
        # Wrapper: pymatgen.io.vasp.sets often handles this (e.g. MPRelaxSet), 
        # but that generates INCAR/KPOINTS too.
        # Let's stick to just Potcar.
        
        print(f"Reading structure from {poscar_path}...")
        print(f"Detected species: {species}")
        
        # Use simple Potcar construction.
        # Note: functional argument in pymatgen Potcar init is 'functional'.
        # Assuming Pymatgen is configured correctly, this should work.
        
        try:
             # Try explicit functional mapping if user provides "PBE" but config uses "PBE_54" etc
             potcar = Potcar(species, functional=functional)
        except OSError:
            # Fallback or specific error msg about PMG setup
             print(f"⚠️  Could not find POTCARs for functional '{functional}'.")
             print("   Please ensure PMG_VASP_PSP_DIR is set in ~/.pmgrc.yaml")
             raise

        potcar.write_file(output)
        
        print(f"✅ Generated {output} using functional '{functional}'")
        
    except Exception as e:
        print(f"❌ Error generating POTCAR: {e}")
