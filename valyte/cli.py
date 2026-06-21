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
from valyte.band_plot import (
    DEFAULT_TRICOLORS,
    plot_band_structure,
    plot_orbital_band_structure,
    plot_spin_texture_band_structure,
)
from valyte.dos_plot import load_dos, plot_dos
from valyte.kpoints import generate_kpoints, generate_kpoints_interactive
from valyte.potcar import generate_potcar
from valyte.ipr import run_ipr, run_ipr_interactive
from valyte.geoopt import check_convergence
from valyte.effmass import compute_effective_masses, print_results, save_results_dat
from valyte.effmass_plot import plot_effective_mass
from valyte.converge import run_converge
from valyte.bandgap import get_bandgap


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
    dos_parser.add_argument("--width", type=float, default=5.0, help="Plot width in inches (default: 5)")
    dos_parser.add_argument("--height", type=float, default=4.0, help="Plot height in inches (default: 4)")
    dos_parser.add_argument("--save-data", action="store_true", help="Save DOS data to valyte_dos.dat")

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
    band_parser.add_argument(
        "--tricolor", nargs=3, metavar=("SPEC1", "SPEC2", "SPEC3"),
        help="Enable tricolor orbital-resolved mode with 3 orbital/element specs. "
             "Formats: 's'|'p'|'d'|'f', 'Fe', 'Fe:d', 'O(p)'. Example: --tricolor s p d",
    )
    band_parser.add_argument(
        "--tricolors", "--colors", nargs=3, metavar=("COLOR1", "COLOR2", "COLOR3"),
        default=DEFAULT_TRICOLORS,
        help="Colors for the 3 tricolor specs (default: red blue green)",
    )
    band_parser.add_argument(
        "--tri-labels", nargs=3, metavar=("LBL1", "LBL2", "LBL3"),
        help="Labels for the triangle legend (defaults to the spec strings)",
    )
    band_parser.add_argument("--lw", type=float, default=2.0, help="Line width for tricolor bands (default: 2.0)")
    band_parser.add_argument("--width", type=float, default=4.0, help="Plot width in inches (default: 4)")
    band_parser.add_argument("--height", type=float, default=4.0, help="Plot height in inches (default: 4)")
    band_parser.add_argument("--save-data", action="store_true", help="Save band data to valyte_band.dat")
    band_parser.add_argument(
        "--spin-resolved", action="store_true",
        help="Plot spin-up and spin-down channels in distinct colors (blue/red). "
             "Only effective for spin-polarized calculations.",
    )
    band_parser.add_argument(
        "--spin-texture", choices=["sx", "sy", "sz"], metavar="COMPONENT",
        help="Plot non-collinear spin texture colored by the chosen component "
             "(sx, sy, or sz). Requires LSORBIT=.TRUE. and LORBIT>=11.",
    )
    band_parser.add_argument(
        "--spin-cmap", default="seismic",
        help="Colormap for --spin-texture (default: seismic).",
    )

    # Band KPOINTS generation
    kpt_gen_parser = band_subparsers.add_parser("kpt-gen", help="Generate KPOINTS for band structure")
    kpt_gen_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    kpt_gen_parser.add_argument("-n", "--npoints", type=int, default=40, help="Points per segment")
    kpt_gen_parser.add_argument("-o", "--output", default="KPOINTS", help="Output filename")
    kpt_gen_parser.add_argument("--symprec", type=float, default=0.01, help="Symmetry precision (default: 0.01)")
    kpt_gen_parser.add_argument("--mode", default="bradcrack", help="Standardization mode (default: bradcrack)")

    # KPOINTS
    kpt_parser = subparsers.add_parser("kpt", help="Generate KPOINTS for SCF/relaxation")
    kpt_parser.add_argument("-i", "--input", help="Input POSCAR file (default: POSCAR)")
    kpt_parser.add_argument("-o", "--output", help="Output KPOINTS file (default: KPOINTS)")
    kpt_parser.add_argument("--spacing", type=float, help="K-spacing in 2*pi/A (default: 0.04)")
    kpt_parser.add_argument("--scheme", help="K-mesh scheme: gamma or mp (default: gamma)")
    kpt_parser.add_argument("--no-potcar", action="store_true", help="Do not auto-generate POTCAR if missing")

    # POTCAR
    potcar_parser = subparsers.add_parser("potcar", help="Generate POTCAR")
    potcar_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    potcar_parser.add_argument("-o", "--output", default="POTCAR", help="Output filename")
    potcar_parser.add_argument("--functional", default="PBE", help="Functional (default: PBE)")

    # IPR
    ipr_parser = subparsers.add_parser("ipr", help="Compute IPR from PROCAR")
    ipr_parser.add_argument("-i", "--input", default="PROCAR", help="Input PROCAR file")
    ipr_parser.add_argument(
        "-b", "--bands", nargs="+",
        help="Band indices/ranges, e.g. '5', '5-8', or '5 8-10'",
    )
    ipr_parser.add_argument("-o", "--output", default="ipr_procar.dat", help="Output data filename")
    ipr_parser.add_argument("--details", action="store_true", help="Print per-k-point IPR values")

    # Force check
    force_check_parser = subparsers.add_parser("force-check", help="Geometry optimization force/energy convergence check")
    force_check_parser.add_argument("outcar", nargs="?", default="OUTCAR", help="Path to OUTCAR file (default: OUTCAR)")
    force_check_parser.add_argument("--ediffg", type=float, default=None, help="Force convergence threshold in eV/Å (e.g. 0.02)")

    # Convergence monitor
    conv_parser = subparsers.add_parser("converge", help="Monitor VASP relaxation/SCF convergence")
    conv_parser.add_argument("path", nargs="?", default=".", help="Directory or OSZICAR path (default: .)")
    conv_parser.add_argument("--electronic", action="store_true", help="Show SCF convergence instead of ionic")
    conv_parser.add_argument("--forces", action="store_true", help="Include max-force panel (requires OUTCAR)")
    conv_parser.add_argument("--stress", action="store_true", help="Include pressure panel (requires OUTCAR)")
    conv_parser.add_argument("--ethresh", type=float, default=1e-4, help="Energy convergence threshold (eV)")
    conv_parser.add_argument("--fthresh", type=float, default=0.02, help="Force convergence threshold (eV/Å)")
    conv_parser.add_argument("--start", type=int, default=1, help="First ionic step to show")
    conv_parser.add_argument("--end", type=int, default=None, help="Last ionic step to show")
    conv_parser.add_argument("-o", "--output", default="valyte_converge.png", help="Output plot filename")
    conv_parser.add_argument("--save-data", action="store_true", help="Save parsed data to valyte_converge.dat")
    conv_parser.add_argument("--no-plot", action="store_true", help="Print terminal summary only")
    conv_parser.add_argument("--mag", action="store_true", help="Include magnetization subplot")

    # Effective mass
    effmass_parser = subparsers.add_parser("effmass", help="Compute carrier effective masses at VBM/CBM")
    effmass_parser.add_argument("--vasprun", default=".", help="Path to vasprun.xml or directory containing it")
    effmass_parser.add_argument("--kpoints", default=None, help="Path to KPOINTS file for high-symmetry labels")
    effmass_parser.add_argument("--npoints", type=int, default=3, help="Fitting points on each side of extremum (default: 3)")
    effmass_parser.add_argument("--band-index", type=int, nargs="+", default=None, help="Manual 1-indexed band indices to fit")
    effmass_parser.add_argument("--plot", action="store_true", help="Save parabolic fit plot")
    effmass_parser.add_argument("-o", "--output", default="valyte_effmass.png", help="Output plot filename (with --plot)")
    effmass_parser.add_argument("--save-data", action="store_true", help="Save results to valyte_effmass.dat")
    effmass_parser.add_argument("--tol", type=float, default=1e-3, help="Degeneracy tolerance in eV (default: 1e-3)")

    # Bandgap
    bandgap_parser = subparsers.add_parser("bandgap", help="Print electronic bandgap")
    bandgap_parser.add_argument("filepath", nargs="?", default=".", help="Path to vasprun.xml or directory containing it")

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
                figsize=(args.width, args.height),
                font=args.font,
                show_fermi=args.fermi,
                show_total=not args.pdos,
                plotting_config=plotting_config,
                legend_cutoff=args.legend_cutoff,
                scale_factor=args.scale,
                save_data=args.save_data,
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
            noninteractive = any([
                args.input is not None,
                args.output is not None,
                args.spacing is not None,
                args.scheme is not None,
                args.no_potcar,
            ])
            if noninteractive:
                generate_kpoints(
                    poscar_path=args.input or "POSCAR",
                    kspacing=args.spacing if args.spacing is not None else 0.04,
                    scheme=args.scheme or "gamma",
                    output=args.output or "KPOINTS",
                    generate_potcar_if_missing=not args.no_potcar,
                )
            else:
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

    elif args.command == "ipr":
        try:
            if args.bands:
                run_ipr(
                    procar_file=args.input,
                    band_text=" ".join(args.bands),
                    output=args.output,
                    show_details=args.details,
                )
            else:
                run_ipr_interactive(procar_file=args.input, output=args.output)
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "force-check":
        try:
            check_convergence(
                outcar_path=args.outcar,
                ediffg=args.ediffg,
            )
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "effmass":
        try:
            results = compute_effective_masses(
                vasprun_path=args.vasprun,
                kpoints_path=args.kpoints,
                npoints=args.npoints,
                band_index=args.band_index,
                tol=args.tol,
            )
            print_results(results)

            if args.plot:
                plot_effective_mass(results, output=args.output)

            if args.save_data:
                save_results_dat(results)
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "converge":
        try:
            run_converge(
                path=args.path,
                electronic=args.electronic,
                forces=args.forces,
                stress=args.stress,
                ethresh=args.ethresh,
                fthresh=args.fthresh,
                start=args.start,
                end=args.end,
                output=args.output,
                save_data=args.save_data,
                no_plot=args.no_plot,
                mag=args.mag,
            )
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "bandgap":
        try:
            bg = get_bandgap(args.filepath)
            print(f"Bandgap = {bg:.4f} eV")
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

                if args.spin_texture:
                    output = args.output if args.output != "valyte_band.png" \
                        else f"valyte_band_{args.spin_texture}.png"
                    plot_spin_texture_band_structure(
                        vasprun_path=target_path,
                        kpoints_path=kpoints_path,
                        output=output,
                        ylim=tuple(args.ylim) if args.ylim else None,
                        figsize=(args.width, args.height),
                        font=args.font,
                        save_data=args.save_data,
                        spin_component=args.spin_texture,
                        cmap=args.spin_cmap,
                    )
                elif args.tricolor:
                    tri_labels = args.tri_labels if args.tri_labels else list(args.tricolor)
                    output = args.output if args.output != "valyte_band.png" else "valyte_band_orbital.png"
                    plot_orbital_band_structure(
                        vasprun_path=target_path,
                        kpoints_path=kpoints_path,
                        output=output,
                        ylim=tuple(args.ylim) if args.ylim else None,
                        tricolor=args.tricolor,
                        tricolors=args.tricolors,
                        tri_labels=tri_labels,
                        lw=args.lw,
                        figsize=(args.width, args.height),
                        font=args.font,
                        save_data=args.save_data,
                    )
                else:
                    plot_band_structure(
                        vasprun_path=target_path,
                        kpoints_path=kpoints_path,
                        output=args.output,
                        ylim=tuple(args.ylim) if args.ylim else None,
                        figsize=(args.width, args.height),
                        font=args.font,
                        save_data=args.save_data,
                        spin_resolved=args.spin_resolved,
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
