#!/usr/bin/env python3
"""Valyte CLI entry point."""

import os
import sys
import argparse
import re
import warnings

# Suppress pymatgen warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

from valyte import __version__
from valyte.supercell import create_supercell
from valyte.band import generate_band_kpoints, resolve_kpath
from valyte.bz_plot import plot_brillouin_zone
from valyte.band_plot import (
    DEFAULT_TRICOLORS,
    plot_band_structure,
    plot_orbital_band_structure,
    plot_spin_texture_band_structure,
)
from valyte.dos_plot import load_dos, plot_dos, plot_dos_panels
from valyte.colors import (PALETTE_NAMES, load_color_file,
                           normalize_color, normalize_mapping)
from matplotlib.colors import is_color_like
from pymatgen.core.periodic_table import Element
from valyte.combined_plot import plot_combined
from valyte.kpoints import generate_kpoints, generate_kpoints_interactive
from valyte.potcar import generate_potcar, parse_potcar_overrides
from valyte.ipr import run_ipr, run_ipr_interactive
from valyte.effmass import compute_effective_masses, print_results, save_results_dat
from valyte.effmass_plot import plot_effective_mass
from valyte.converge import run_converge
from valyte.bandgap import get_bandgap
from valyte.style import STYLE_NAMES

# seekpath still reads spglib datasets as dicts, which spglib >= 2.5 deprecates.
# spglib registers its own filter for this on import, so ours must come after
# the imports above (which pull in spglib) to take precedence.
warnings.filterwarnings("ignore", category=DeprecationWarning, message=r"dict interface")


def _apply_format(output, fmt):
    if not fmt:
        return output
    base, _ = os.path.splitext(output)
    return f"{base}.{fmt}"


def _normalize_element(symbol: str) -> str:
    if not symbol:
        return symbol
    return symbol[0].upper() + symbol[1:].lower()


_SELECTION = re.compile(
    r"^([A-Za-z]{1,3})"                                  # element
    r"(?:\(\s*([spdfSPDF])\s*\)|:([spdfSPDF])(?=$|[:=]))?"  # (orb) or :orb
    r"(?:[:=](.+))?$"                                    # optional :colour or =colour
)


def parse_element_selection(inputs):
    """Parse -e selections: "Fe", "Fe(d)", "Fe:d", optionally with a colour
    ("Fe(d):red", "Fe:d=#0072B2", "O=0072B2").

    Returns (elements_to_load, plotting_config, inline_colors). Raises
    ValueError for anything that is not an element symbol, so a typo is
    reported instead of silently dropped.
    """
    if not inputs:
        return None, None, {}

    elements_to_load = []
    plotting_config = []
    inline_colors = {}

    for item in inputs:
        match = _SELECTION.match(item.strip())
        el = _normalize_element(match.group(1)) if match else None
        if not match or not Element.is_valid_symbol(el):
            raise ValueError(
                f"Could not read -e entry {item!r}. Use an element (Fe), an element "
                "and orbital (Fe(d) or Fe:d), optionally with a colour (Fe(d):red).")

        orb = (match.group(2) or match.group(3) or "").lower() or None
        color = match.group(4)

        if el not in elements_to_load:
            elements_to_load.append(el)

        orb_str = orb or "total"
        if (el, orb_str) not in plotting_config:
            plotting_config.append((el, orb_str))

        if color:
            label = el if orb_str == "total" else f"{el}({orb_str})"
            inline_colors[label] = normalize_color(color, f"colour for {item!r}")

    return elements_to_load, plotting_config, inline_colors


def _build_colors(colors_file, colors_args, inline_colors):
    """Combine --colors-file, -c/--colors and inline -e colours.

    Returns (mapping, palette): mapping is {series: colour} or None; palette is
    the list of plain colours given to -c (which replaces the palette) or None.
    Later sources win: file, then -c mappings, then inline -e colours.
    """
    mapping = dict(load_color_file(colors_file)) if colors_file else {}
    positional = []
    for item in colors_args or []:
        item = item.strip()
        if "=" in item:
            key, value = item.split("=", 1)
        elif is_color_like(item) or is_color_like("#" + item):
            positional.append(normalize_color(item))
            continue
        elif ":" in item:
            key, value = item.rsplit(":", 1)
        else:
            raise ValueError(
                f"Could not read -c entry {item!r}: not a colour, and not a "
                "mapping like Fe(d)=red.")
        mapping[key] = value
    mapping.update(inline_colors or {})
    return (normalize_mapping(mapping) or None), (positional or None)


def _fraction(text):
    """argparse type for values between 0 and 1."""
    value = float(text)
    if not 0 <= value <= 1:
        raise argparse.ArgumentTypeError(f"must be between 0 and 1, got {text}")
    return value


def _dos_colour_args(args):
    """Parse -e, -c, --colors-file and --palette for dos/combined.

    Returns (elements, plotting_config, colors, palette). Errors in any of
    them are reported as a single "Error:" line and exit, never a traceback.
    """
    try:
        elements, plotting_config, inline_colors = parse_element_selection(args.elements)
        colors, custom_palette = _build_colors(args.colors_file, args.colors, inline_colors)
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")
        sys.exit(1)
    if custom_palette and args.palette:
        print(f"Note: the colours given to -c replace --palette {args.palette}.")
    # None lets the plot use the style's own palette (wong for the default style)
    return elements, plotting_config, colors, custom_palette or args.palette


def main():
    parser = argparse.ArgumentParser(description="Valyte: VASP Post-Processing Tool")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
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
    dos_parser.add_argument("--legend-cutoff", type=_fraction, default=0.10, help="Smallest peak, as a fraction of the largest, that gets a colour and a legend entry; smaller channels are drawn grey (0.0-1.0, default 0.10). Ignored for channels named with -e")
    dos_parser.add_argument("--font", default=None, help="Font family (default: Arial, or the style's own font)")
    dos_parser.add_argument("--width", type=float, default=4.45, help="Plot width in inches (default: 4.45)")
    dos_parser.add_argument("--height", type=float, default=4.2, help="Plot height in inches (default: 4.2)")
    dos_parser.add_argument("--save-data", action="store_true", help="Also save DOS data as a .dat named after -o")
    dos_parser.add_argument("--format", choices=["png", "pdf", "svg"], help="Output figure format")
    dos_parser.add_argument("--panels", action="store_true",
                            help="Split DOS into stacked panels (one per element)")
    dos_parser.add_argument("--panel-by", choices=["element", "orbital"], default="element",
                            help="Grouping mode for panels: 'element' (default) or 'orbital'")
    dos_parser.add_argument("--no-bold", action="store_true", help="Use normal font weight and thinner lines/ticks")
    dos_parser.add_argument("--dpi", type=int, default=400, help="Output figure resolution DPI (default: 400)")
    dos_parser.add_argument("--fontsize", type=float, default=None, help="Base font size in points; all labels scale with it")
    dos_parser.add_argument("--style", choices=STYLE_NAMES, default="default", help="Plot style: default, xkcd (hand-drawn), legacy (classic boxed axes, outward ticks), lines (no DOS fills), dark (for dark slides), neon (glowing lines on dark), grid (light grid lines)")
    dos_parser.add_argument("--colors-file", help="Path to JSON or text file mapping elements/orbitals to colors")
    dos_parser.add_argument("-c", "--colors", nargs="+", help="Custom colours: per series (Fe(d)=red O(p)=0072B2) and/or a plain list that replaces the palette (red blue)")
    dos_parser.add_argument("--palette", choices=PALETTE_NAMES, default=None, help="Colour palette for the projected DOS: wong (colour-blind safe), bright, or muted. Default: the style's own palette (wong for --style default and lines). -c/--colors and --colors-file override single series")


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
    band_parser.add_argument("--font", default=None, help="Font family (default: Arial, or the style's own font)")
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
    band_parser.add_argument("--width", type=float, default=4.2, help="Plot width in inches (default: 4.2)")
    band_parser.add_argument("--height", type=float, default=4.2, help="Plot height in inches (default: 4.2)")
    band_parser.add_argument("--save-data", action="store_true", help="Also save band data as a .dat named after -o")
    band_parser.add_argument("--format", choices=["png", "pdf", "svg"], help="Output figure format")
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
    band_parser.add_argument("--no-bold", action="store_true", help="Use normal font weight and thinner lines/ticks")
    band_parser.add_argument("--dpi", type=int, default=400, help="Output figure resolution DPI (default: 400)")
    band_parser.add_argument("--fontsize", type=float, default=None, help="Base font size in points; all labels scale with it")
    band_parser.add_argument("--style", choices=STYLE_NAMES, default="default", help="Plot style: default, xkcd (hand-drawn), legacy (classic boxed axes, outward ticks), dark (for dark slides), neon (glowing lines on dark), grid (light grid lines). lines only changes DOS plots")

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
    potcar_parser.add_argument("--functional", default="PBE_64", help="POTCAR set (default: PBE_64)")
    potcar_parser.add_argument("--potcar", nargs="+", metavar="EL=NAME",
                               help="Override the VASP-recommended POTCAR per element, e.g. Fe=Fe_pv O=O_h")

    # IPR
    ipr_parser = subparsers.add_parser("ipr", help="Compute IPR from PROCAR")
    ipr_parser.add_argument("-i", "--input", default="PROCAR", help="Input PROCAR file")
    ipr_parser.add_argument(
        "-b", "--bands", nargs="+",
        help="Band indices/ranges, e.g. '5', '5-8', or '5 8-10'",
    )
    ipr_parser.add_argument("-o", "--output", default="ipr_procar.dat", help="Output data filename")
    ipr_parser.add_argument("--details", action="store_true", help="Print per-k-point IPR values")

    # Brillouin zone
    bz_parser = subparsers.add_parser("bz", help="Plot the suggested k-path in the first Brillouin zone")
    bz_parser.add_argument("-i", "--input", default="POSCAR", help="Input POSCAR file")
    bz_parser.add_argument("-o", "--output", default="valyte_bz.png", help="Output plot filename")
    bz_parser.add_argument("--symprec", type=float, default=0.01, help="Symmetry precision (default: 0.01)")
    bz_parser.add_argument("--mode", default="bradcrack",
                           choices=["bradcrack", "seekpath", "hinuma", "setyawan_curtarolo", "latimer_munro"],
                           help="K-path convention (default: bradcrack)")
    bz_parser.add_argument("--elev", type=float, default=22, help="Viewing elevation in degrees (default: 22)")
    bz_parser.add_argument("--azim", type=float, default=30, help="Viewing azimuth in degrees (default: 30)")
    bz_parser.add_argument("--width", type=float, default=5.5, help="Plot width in inches (default: 5.5)")
    bz_parser.add_argument("--height", type=float, default=5.5, help="Plot height in inches (default: 5.5)")
    bz_parser.add_argument("--format", choices=["png", "pdf", "svg"], help="Output figure format")
    bz_parser.add_argument("--font", default=None, help="Font family: Arial, Helvetica, Times New Roman (default: matplotlib math font)")
    bz_parser.add_argument("--no-bold", action="store_true", help="Use normal font weight and thinner lines/ticks")
    bz_parser.add_argument("--dpi", type=int, default=400, help="Output figure resolution DPI (default: 400)")
    bz_parser.add_argument("--fontsize", type=float, default=None, help="Base font size in points; all labels scale with it")
    bz_parser.add_argument("--style", choices=STYLE_NAMES, default="default", help="Plot style: default, xkcd (hand-drawn), legacy (classic boxed axes, outward ticks), dark (for dark slides), neon (glowing lines on dark), grid (light grid lines). lines only changes DOS plots")

    # Convergence monitor
    conv_parser = subparsers.add_parser("converge", help="Monitor VASP structural relaxation convergence")
    conv_parser.add_argument("path", nargs="?", default=".", help="Directory or OSZICAR path (default: .)")
    conv_parser.add_argument("--save-data", action="store_true", help="Save parsed data to valyte_converge.dat")

    # Effective mass
    effmass_parser = subparsers.add_parser("effmass", help="Compute carrier effective masses at VBM/CBM")
    effmass_parser.add_argument("--vasprun", default=".", help="Path to vasprun.xml or directory containing it")
    effmass_parser.add_argument("--kpoints", default=None, help="Path to KPOINTS file for high-symmetry labels")
    effmass_parser.add_argument("--npoints", type=int, default=3, help="Fitting points on each side of extremum (default: 3)")
    effmass_parser.add_argument("--band-index", type=int, nargs="+", default=None, help="Manual 1-indexed band indices to fit")
    effmass_parser.add_argument("--plot", action="store_true", help="Save parabolic fit plot")
    effmass_parser.add_argument("-o", "--output", default="valyte_effmass.png", help="Output plot filename (with --plot)")
    effmass_parser.add_argument("--save-data", action="store_true", help="Save results as a .dat named after -o (default: valyte_effmass.dat)")
    effmass_parser.add_argument("--tol", type=float, default=1e-3, help="Degeneracy tolerance in eV (default: 1e-3)")
    effmass_parser.add_argument("--format", choices=["png", "pdf", "svg"], help="Output figure format")
    effmass_parser.add_argument("--font", default=None, help="Font family (default: Arial, or the style's own font)")
    effmass_parser.add_argument("--no-bold", action="store_true", help="Use normal font weight and thinner lines/ticks")
    effmass_parser.add_argument("--dpi", type=int, default=400, help="Output figure resolution DPI (default: 400)")
    effmass_parser.add_argument("--fontsize", type=float, default=None, help="Base font size in points; all labels scale with it")
    effmass_parser.add_argument("--style", choices=STYLE_NAMES, default="default", help="Plot style: default, xkcd (hand-drawn), legacy (classic boxed axes, outward ticks), dark (for dark slides), neon (glowing lines on dark), grid (light grid lines). lines only changes DOS plots")

    # Bandgap
    bandgap_parser = subparsers.add_parser("bandgap", help="Print electronic bandgap")
    bandgap_parser.add_argument("filepath", nargs="?", default=".", help="Path to vasprun.xml or directory containing it")

    # Combined Band + DOS Plot
    combined_parser = subparsers.add_parser("combined", help="Plot Band Structure and DOS side-by-side")
    combined_parser.add_argument("filepath", nargs="?", help="Path to vasprun.xml or directory containing it")
    combined_parser.add_argument("--vasprun", help="Explicit path to vasprun.xml (alternative to positional argument)")
    combined_parser.add_argument("--kpoints", help="Path to KPOINTS file (for labels)")
    combined_parser.add_argument("--dos", help="Path to vasprun.xml or directory for DOS data (if different from band structure)")
    combined_parser.add_argument("-e", "--elements", nargs="+", help="Elements/Orbitals to plot in DOS (e.g., 'Fe' or 'Fe(d)')")
    combined_parser.add_argument("-o", "--output", default="valyte_combined.png", help="Output filename")
    combined_parser.add_argument("--ylim", nargs=2, type=float, default=[-4, 4], help="Energy range (min max)")
    combined_parser.add_argument("--scale", type=float, default=1.0, help="Scaling factor for DOS axis")
    combined_parser.add_argument("--fermi", action="store_true", help="Draw dashed horizontal line at Fermi level (E=0)")
    combined_parser.add_argument("--legend-cutoff", type=_fraction, default=0.10, help="Smallest peak, as a fraction of the largest, that gets a colour and a legend entry; smaller channels are drawn grey (0.0-1.0, default 0.10). Ignored for channels named with -e")
    combined_parser.add_argument("--font", default=None, help="Font family (default: Arial, or the style's own font)")
    combined_parser.add_argument("--width", type=float, default=3.2, help="Plot width in inches (default: 3.2)")
    combined_parser.add_argument("--height", type=float, default=3.2, help="Plot height in inches (default: 3.2)")
    combined_parser.add_argument("--save-data", action="store_true", help="Save band and DOS data to text files")
    combined_parser.add_argument("--format", choices=["png", "pdf", "svg"], help="Output figure format")
    combined_parser.add_argument("--spin-resolved", action="store_true", help="Plot spin-up/spin-down channels in distinct colors")
    combined_parser.add_argument("--no-bold", action="store_true", help="Use normal font weight and thinner lines/ticks")
    combined_parser.add_argument("--dpi", type=int, default=400, help="Output figure resolution DPI (default: 400)")
    combined_parser.add_argument("--fontsize", type=float, default=None, help="Base font size in points; all labels scale with it")
    combined_parser.add_argument("--style", choices=STYLE_NAMES, default="default", help="Plot style: default, xkcd (hand-drawn), legacy (classic boxed axes, outward ticks), lines (no DOS fills), dark (for dark slides), neon (glowing lines on dark), grid (light grid lines)")
    combined_parser.add_argument("--colors-file", help="Path to JSON or text file mapping elements/orbitals to colors")
    combined_parser.add_argument("-c", "--colors", nargs="+", help="Custom colours: per series (Fe(d)=red O(p)=0072B2) and/or a plain list that replaces the palette (red blue)")
    combined_parser.add_argument("--palette", choices=PALETTE_NAMES, default=None, help="Colour palette for the projected DOS: wong (colour-blind safe), bright, or muted. Default: the style's own palette (wong for --style default and lines). -c/--colors and --colors-file override single series")


    args = parser.parse_args()

    if args.command == "dos":
        target_path = args.filepath if args.filepath else args.vasprun
        if not target_path:
            target_path = "."

        elements, plotting_config, colors_map, palette = _dos_colour_args(args)

        try:
            dos_data, pdos_data = load_dos(target_path, elements)

            if args.panels:
                default_out = "valyte_dos_panels.png"
                out_file = args.output if args.output != "valyte_dos.png" else default_out
                plot_dos_panels(
                    dos_data,
                    pdos_data,
                    out=_apply_format(out_file, args.format),
                    xlim=tuple(args.xlim),
                    ylim=tuple(args.ylim) if args.ylim else None,
                    figsize=(args.width, args.height) if (args.width != 4.45 or args.height != 4.2) else None,
                    font=args.font,
                    show_fermi=args.fermi,
                    show_total=not args.pdos,
                    plotting_config=plotting_config,
                    scale_factor=args.scale,
                    save_data=args.save_data,
                    group_by=args.panel_by,
                    bold=not args.no_bold,
                    dpi=args.dpi,
                    colors=colors_map,
                    palette=palette,
                    fontsize=args.fontsize,
                    style=args.style,
                )
            else:
                plot_dos(
                    dos_data,
                    pdos_data,
                    out=_apply_format(args.output, args.format),
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
                    bold=not args.no_bold,
                    dpi=args.dpi,
                    colors=colors_map,
                    palette=palette,
                    fontsize=args.fontsize,
                    style=args.style,
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
                overrides=parse_potcar_overrides(args.potcar),
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
                plot_effective_mass(results, output=_apply_format(args.output, args.format), dpi=args.dpi,
                                    font=args.font,
                                    bold=not args.no_bold, fontsize=args.fontsize,
                                    style=args.style)

            if args.save_data:
                save_results_dat(results, os.path.splitext(args.output)[0] + ".dat")
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "bz":
        try:
            prim_std, path, kpoints = resolve_kpath(
                args.input, symprec=args.symprec, mode=args.mode)
            plot_brillouin_zone(
                prim_std, path, kpoints,
                output=_apply_format(args.output, args.format),
                figsize=(args.width, args.height),
                dpi=args.dpi,
                bold=not args.no_bold,
                fontsize=args.fontsize,
                font=args.font,
                style=args.style,
                elev=args.elev,
                azim=args.azim,
            )
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    elif args.command == "converge":
        try:
            run_converge(
                path=args.path,
                save_data=args.save_data,
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

    elif args.command == "combined":
        target_path = args.filepath if args.filepath else args.vasprun
        if not target_path:
            target_path = "."

        elements, plotting_config, colors_map, palette = _dos_colour_args(args)

        if os.path.isdir(target_path):
            vasprun_file = os.path.join(target_path, "vasprun.xml")
        else:
            vasprun_file = target_path

        kpoints_file = args.kpoints
        if not kpoints_file:
            base_dir = os.path.dirname(vasprun_file) or "."
            potential_kpoints = os.path.join(base_dir, "KPOINTS")
            if os.path.exists(potential_kpoints):
                kpoints_file = potential_kpoints

        try:
            plot_combined(
                vasprun_path=vasprun_file,
                kpoints_path=kpoints_file,
                dos_path=args.dos,
                elements=elements,
                plotting_config=plotting_config,
                output=_apply_format(args.output, args.format),
                ylim=tuple(args.ylim),
                figsize=(args.width, args.height),
                font=args.font,
                show_fermi=args.fermi,
                scale_factor=args.scale,
                legend_cutoff=args.legend_cutoff,
                save_data=args.save_data,
                spin_resolved=args.spin_resolved,
                bold=not args.no_bold,
                dpi=args.dpi,
                colors=colors_map,
                palette=palette,
                fontsize=args.fontsize,
                style=args.style,
            )

        except Exception:
            import traceback
            traceback.print_exc()
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
        else:
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
                        output=_apply_format(output, args.format),
                        ylim=tuple(args.ylim) if args.ylim else None,
                        figsize=(args.width, args.height),
                        font=args.font,
                        save_data=args.save_data,
                        spin_component=args.spin_texture,
                        cmap=args.spin_cmap,
                        bold=not args.no_bold,
                        dpi=args.dpi,
                        fontsize=args.fontsize,
                        style=args.style,
                    )
                elif args.tricolor:
                    tri_labels = args.tri_labels if args.tri_labels else list(args.tricolor)
                    output = args.output if args.output != "valyte_band.png" else "valyte_band_orbital.png"
                    plot_orbital_band_structure(
                        vasprun_path=target_path,
                        kpoints_path=kpoints_path,
                        output=_apply_format(output, args.format),
                        ylim=tuple(args.ylim) if args.ylim else None,
                        tricolor=args.tricolor,
                        tricolors=args.tricolors,
                        tri_labels=tri_labels,
                        lw=args.lw,
                        figsize=(args.width, args.height),
                        font=args.font,
                        save_data=args.save_data,
                        bold=not args.no_bold,
                        dpi=args.dpi,
                        fontsize=args.fontsize,
                        style=args.style,
                    )
                else:
                    plot_band_structure(
                        vasprun_path=target_path,
                        kpoints_path=kpoints_path,
                        output=_apply_format(args.output, args.format),
                        ylim=tuple(args.ylim) if args.ylim else None,
                        figsize=(args.width, args.height),
                        font=args.font,
                        save_data=args.save_data,
                        spin_resolved=args.spin_resolved,
                        bold=not args.no_bold,
                        dpi=args.dpi,
                        fontsize=args.fontsize,
                        style=args.style,
                    )
            except Exception:
                import traceback

                traceback.print_exc()
                sys.exit(1)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
