#!/usr/bin/env python3
"""Valyte CLI entry point."""

import os
import sys
import argparse
import re
import warnings

# Suppress pymatgen warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

from valyte.supercell import create_supercell
from valyte.band import generate_band_kpoints
from valyte.band_plot import plot_band_structure
from valyte.dos_plot import load_dos, plot_dos
from valyte.kpoints import generate_kpoints_interactive
from valyte.potcar import generate_potcar


def _normalize_element(symbol: str) -> str:
    if not symbol:
        return symbol
    return symbol[0].upper() + symbol[1:].lower()


def parse_element_selection(inputs):
    """Parse element/orbital selections like "Fe", "Fe(d)", "O(p)"."""
    if not inputs:
        return None, None

    elements_to_load = []
    elements_seen = set()
    plotting_config = []

    pattern = re.compile(r"^([A-Za-z]+)(?:\(([spdf])\))?$", re.IGNORECASE)

    for item in inputs:
        match = pattern.match(item.strip())
        if not match:
            print(f"Warning: could not parse '{item}'. Ignoring.")
            continue

        el = _normalize_element(match.group(1))
        orb = match.group(2).lower() if match.group(2) else None

        if el not in elements_seen:
            elements_to_load.append(el)
            elements_seen.add(el)

        plotting_config.append((el, orb or "total"))

    return elements_to_load, plotting_config


def main():
    parser = argparse.ArgumentParser(description="Valyte: VASP Post-Processing Tool")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # DOS
    dos_parser = subparsers.add_parser("dos", help="Plot Density of States (DOS)")
    dos_parser.add_argument("filepath", nargs="?", help="Path to vasprun.xml or directory containing it")
    dos_parser.add_argument("--vasprun", help="Explicit path to vasprun.xml (alternative to positional argument)")
    dos_parser.add_argument("-e", "--elements", nargs="+", help="Elements/Orbitals to plot (e.g., 'Fe O' or 'Fe(d) O(p)')")
    dos_parser.add_argument("-o", "--output", default="valyte_dos.png", help="Output filename")
    dos_parser.add_argument("--xlim", nargs=2, type=float, default=[-6, 6], help="Energy range (min max)")
    dos_parser.add_argument("--ylim", nargs=2, type=float, help="DOS range (min max)")
    dos_parser.add_argument("--scale", type=float, default=1.0, help="Scaling factor for Y-axis (zoom in)")
    dos_parser.add_argument("--fermi", action="store_true", help="Draw dashed line at Fermi level (E=0)")
    dos_parser.add_argument("--pdos", action="store_true", help="Plot only Projected DOS (hide Total DOS)")
    dos_parser.add_argument("--legend-cutoff", type=float, default=0.10, help="Threshold for legend visibility (0.0-1.0)")
    dos_parser.add_argument("--font", default="Arial", help="Font family")

    # Supercell
    supercell_parser = subparsers.add_parser("supercell", help="Create a supercell")
    supercell_parser.add_argument("nx", type=int, help="Supercell size x")
    supercell_parser.add_argument("ny", type=int, help="Supercell size y")
    supercell_parser.add_argument("nz", type=int, help="Supercell size z")
    supercell_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    supercell_parser.add_argument("-o", "--output", default="POSCAR_supercell", help="Output filename")

    # Band
    band_parser = subparsers.add_parser("band", help="Band structure utilities")
    band_subparsers = band_parser.add_subparsers(dest="band_command", help="Band commands")

    # Band plotting (default)
    band_parser.add_argument("--vasprun", default=".", help="Path to vasprun.xml or directory")
    band_parser.add_argument("--kpoints", help="Path to KPOINTS file (for labels)")
    band_parser.add_argument("-o", "--output", default="valyte_band.png", help="Output filename")
    band_parser.add_argument("--ylim", nargs=2, type=float, help="Energy range (min max)")
    band_parser.add_argument("--font", default="Arial", help="Font family")

    # Band KPOINTS generation
    kpt_gen_parser = band_subparsers.add_parser("kpt-gen", help="Generate KPOINTS for band structure")
    kpt_gen_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    kpt_gen_parser.add_argument("-n", "--npoints", type=int, default=40, help="Points per segment")
    kpt_gen_parser.add_argument("-o", "--output", default="KPOINTS", help="Output filename")
    kpt_gen_parser.add_argument("--symprec", type=float, default=0.01, help="Symmetry precision (default: 0.01)")
    kpt_gen_parser.add_argument("--mode", default="bradcrack", help="Standardization mode (default: bradcrack)")

    # KPOINTS (interactive)
    subparsers.add_parser("kpt", help="Interactive K-Point Generation (SCF)")

    # POTCAR
    potcar_parser = subparsers.add_parser("potcar", help="Generate POTCAR")
    potcar_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    potcar_parser.add_argument("-o", "--output", default="POTCAR", help="Output filename")
    potcar_parser.add_argument("--functional", default="PBE", help="Functional (default: PBE)")

    args = parser.parse_args()

    if args.command == "dos":
        target_path = args.filepath if args.filepath else args.vasprun
        if not target_path:
            target_path = "."

        elements, plotting_config = parse_element_selection(args.elements)

        try:
            dos_data, pdos_data = load_dos(target_path, elements)
            plot_dos(
                dos_data,
                pdos_data,
                out=args.output,
                xlim=tuple(args.xlim),
                ylim=tuple(args.ylim) if args.ylim else None,
                font=args.font,
                show_fermi=args.fermi,
                show_total=not args.pdos,
                plotting_config=plotting_config,
                legend_cutoff=args.legend_cutoff,
                scale_factor=args.scale,
            )
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "supercell":
        try:
            create_supercell(
                poscar_path=args.input,
                nx=args.nx,
                ny=args.ny,
                nz=args.nz,
                output=args.output,
            )
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "kpt":
        try:
            generate_kpoints_interactive()
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "potcar":
        try:
            generate_potcar(
                poscar_path=args.input,
                functional=args.functional,
                output=args.output,
            )
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "band":
        if args.band_command == "kpt-gen":
            try:
                generate_band_kpoints(
                    poscar_path=args.input,
                    npoints=args.npoints,
                    output=args.output,
                    symprec=args.symprec,
                    mode=args.mode,
                )
            except Exception as e:
                print(f"Error: {e}")
                sys.exit(1)
        elif args.band_command == "plot" or args.band_command is None:
            try:
                target_path = args.vasprun or "."
                if os.path.isdir(target_path):
                    target_path = os.path.join(target_path, "vasprun.xml")

                kpoints_path = args.kpoints
                if not kpoints_path:
                    base_dir = os.path.dirname(target_path)
                    potential_kpoints = os.path.join(base_dir, "KPOINTS")
                    if os.path.exists(potential_kpoints):
                        kpoints_path = potential_kpoints

                plot_band_structure(
                    vasprun_path=target_path,
                    kpoints_path=kpoints_path,
                    output=args.output,
                    ylim=tuple(args.ylim) if args.ylim else None,
                    font=args.font,
                )
            except Exception:
                import traceback

                traceback.print_exc()
                sys.exit(1)
        else:
            band_parser.print_help()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
