<p align="center">
  <img src="https://raw.githubusercontent.com/nikyadav002/Valyte-Project/main/valyte/Logo.png" alt="Valyte Logo" width="100%"/>
</p>

<p align="center">
  <a href="https://pypi.org/project/valyte/"><img src="https://img.shields.io/pypi/v/valyte?color=7c3aed&label=PyPI" alt="PyPI Version"></a>
  <a href="https://pypi.org/project/valyte/"><img src="https://img.shields.io/pypi/pyversions/valyte?color=7c3aed" alt="Python Versions"></a>
  <a href="https://valyte-project.readthedocs.io/en/latest/"><img src="https://readthedocs.org/projects/valyte-project/badge/?version=latest" alt="Docs"></a>
  <a href="https://github.com/nikyadav002/Valyte-Project/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-2a9d8f" alt="License"></a>
</p>

<p align="center">
  <strong>VASP pre- and post-processing from one command line.</strong>
</p>

---

If you run VASP, you know how this goes. The calculation finishes, and then you
go looking for the plotting script you wrote six months ago and can't quite
remember how to call. Valyte started as a way to stop doing that.

Point it at a finished run and it reads what VASP left behind (`vasprun.xml`,
`PROCAR`, `OSZICAR`, `OUTCAR`, `POSCAR`) and hands you a figure or a data file.
One command per job. The defaults try to give you something you could drop into
a paper as-is, and when they don't suit you, almost all of them can be changed.

## Installation

```bash
pip install valyte
```

To upgrade:

```bash
pip install --upgrade valyte
```

Or from source, if you'd like to poke at the code:

```bash
git clone https://github.com/nikyadav002/Valyte-Project
cd Valyte-Project
pip install -e .
```

### Requirements

Python 3.9 or newer. Everything else (`numpy`, `matplotlib`, `pymatgen`,
`scipy`, `seekpath`) comes along with the install, so there's nothing else to
set up.

The one exception is `valyte potcar`, which needs pymatgen to know where your
pseudopotentials live. The [pymatgen POTCAR setup notes](https://pymatgen.org/installation.html#potcar-setup)
cover that, and it's a one-time thing.

## Quick start

Change into a directory with your VASP output and try any of these:

```bash
valyte dos                      # density of states
valyte band                     # band structure
valyte converge                 # relaxation convergence table
```

If you'd rather see the whole thing end to end, the
[Getting Started guide](https://valyte-project.readthedocs.io/en/latest/getting-started/)
walks through a calculation from setup to finished figure.

## Gallery

<p align="center">
  <img src="https://raw.githubusercontent.com/nikyadav002/Valyte-Project/main/valyte/valyte_dos.png" alt="DOS Plot Example" width="47%"/>
  <img src="https://raw.githubusercontent.com/nikyadav002/Valyte-Project/main/valyte/valyte_band.png" alt="Band Structure Example" width="38%"/>
</p>

<p align="center">
  <em>Left: Orbital-resolved density of states with gradient fills. Right: Color-coded band structure with VBM at 0 eV.</em>
</p>

## Commands

### Setting up a calculation

| Command | What it does |
|---|---|
| `valyte supercell nx ny nz` | Build a supercell from a POSCAR |
| `valyte kpt` | Write a KPOINTS grid (Monkhorst-Pack or Gamma). Runs interactively if you pass no flags |
| `valyte band kpt-gen` | Write a line-mode KPOINTS along a high-symmetry path, Bradley-Cracknell by default |
| `valyte potcar` | Concatenate a POTCAR for the species in a POSCAR |

### Analysing the output

| Command | What it does |
|---|---|
| `valyte dos` | Total and projected DOS, orbital resolved, with gradient fills |
| `valyte dos --panels` | The same DOS split into stacked panels, one per element (`--panel-by orbital` for orbitals instead) |
| `valyte band` | Band structure with the VBM placed at 0 eV |
| `valyte band --tricolor s p d` | Orbital-projected bands colored by three specs. Accepts `s`, `Fe`, `Fe:d`, `O(p)` |
| `valyte band --spin-resolved` | Spin-up and spin-down channels in separate colors |
| `valyte band --spin-texture sz` | Non-collinear spin texture, bands colored by `sx`, `sy` or `sz` |
| `valyte combined` | Band structure and DOS side by side on a shared energy axis |
| `valyte effmass` | Carrier effective masses at the VBM and CBM by parabolic fitting |
| `valyte ipr` | Inverse participation ratio from a PROCAR, for judging localization |
| `valyte bandgap` | Print the band gap and nothing else |
| `valyte converge` | Per-step energy, force and pressure table for a structural relaxation |

A couple of these need the right flags set in VASP: `--tricolor` and
`--spin-texture` read projections out of `vasprun.xml`, so the run needs
`LORBIT >= 11`, and spin texture needs a non-collinear calculation on top of
that. If a plot comes out empty, that's usually why.

### Shared options

The four commands that draw figures (`dos`, `band`, `combined`, `effmass`) all
take the same output flags, so once you know them they work everywhere:

| Flag | Effect |
|---|---|
| `--save-data` | Also write the plotted numbers to a `.dat` file |
| `--format {png,pdf,svg}` | Figure format |
| `--dpi` | Resolution for raster output (default 400) |
| `--no-bold` | Lighter type and thinner lines, closer to a journal house style |

A few small exceptions: `valyte effmass` only draws its fit if you ask for it
with `--plot`, `valyte converge` takes `--save-data` even though it prints a
table rather than a figure, and `valyte ipr` always writes `ipr_procar.dat`
(pass `-o` to name it something else).

Every command has `--help` if you just want to see what's there, and the
[CLI reference](https://valyte-project.readthedocs.io/en/latest/cli-reference/)
has all of it on one page.

## Documentation

**[valyte-project.readthedocs.io](https://valyte-project.readthedocs.io/en/latest/)**

| Page | What you'll find |
|---|---|
| [Getting Started](https://valyte-project.readthedocs.io/en/latest/getting-started/) | Installation, prerequisites, and your first plot |
| [Band Structure](https://valyte-project.readthedocs.io/en/latest/band/) | Standard, tricolor, spin-resolved, and spin-texture modes |
| [Density of States](https://valyte-project.readthedocs.io/en/latest/dos/) | Total and projected DOS with orbital resolution |
| [Effective Mass](https://valyte-project.readthedocs.io/en/latest/effmass/) | Carrier effective masses from parabolic fitting |
| [Convergence](https://valyte-project.readthedocs.io/en/latest/converge/) | Relaxation convergence: energy, force, and pressure |
| [IPR](https://valyte-project.readthedocs.io/en/latest/ipr/) | Wavefunction localization analysis |
| [Pre-processing](https://valyte-project.readthedocs.io/en/latest/preprocessing/) | Supercells, k-points, and POTCAR generation |
| [CLI Reference](https://valyte-project.readthedocs.io/en/latest/cli-reference/) | Every command and flag in one searchable page |
| [FAQ](https://valyte-project.readthedocs.io/en/latest/faq/) | Common issues and troubleshooting |

## Contributing

Bug reports, ideas and pull requests are all genuinely welcome, and you don't
need to be sure it's a real bug before saying something. If Valyte fell over on
your files, the VASP output that caused it is the most useful thing you can
send, because it's almost always an edge case in someone's calculation that I
haven't seen yet.

- [Open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new) for a bug or an idea
- [Open a pull request](https://github.com/nikyadav002/Valyte-Project/pulls) if you've already got a fix

[CONTRIBUTING.md](CONTRIBUTING.md) covers the development setup and how the code
is laid out, if you want to dig in.

## Acknowledgements

Valyte stands on pymatgen, seekpath, numpy, scipy and matplotlib, and wouldn't
be much without them. Thanks to everyone who maintains those, and to
[sumo](https://github.com/SMTG-Bham/sumo) and the other open-source VASP tools
that worked this out before I did.

## License

Released under the [MIT License](LICENSE).

---

<p align="center">
  Built by <a href="https://github.com/nikyadav002">Nikhil Singh</a>
</p>
