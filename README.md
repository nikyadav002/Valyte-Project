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
  <strong>Publication-quality VASP pre- and post-processing from a single CLI.</strong>
</p>

---

Valyte turns raw VASP output into clean, publication-ready plots and analysis — band structures, density of states, effective masses, convergence diagnostics, and more — all from one command-line tool. No boilerplate scripts, no Jupyter notebooks, no manual formatting.

## ✨ Highlights

- **One command, one plot** — `valyte band`, `valyte dos`, `valyte converge` — each produces a publication-quality figure with zero configuration.
- **Orbital & spin resolution** — Tricolor orbital-projected bands, spin-resolved channels, and non-collinear spin textures out of the box.
- **Smart pre-processing** — Automatic high-symmetry k-paths (Bradley–Cracknell), supercell generation, and POTCAR handling.
- **Data export** — Every plot command supports `--save-data` to export raw data for custom post-processing.
- **Beautiful defaults** — Gradient fills, adaptive legends, and clean typography that look great in papers without tweaking.

## 📸 Gallery

<p align="center">
  <img src="https://raw.githubusercontent.com/nikyadav002/Valyte-Project/main/valyte/valyte_dos.png" alt="DOS Plot Example" width="47%"/>
  <img src="https://raw.githubusercontent.com/nikyadav002/Valyte-Project/main/valyte/valyte_band.png" alt="Band Structure Example" width="38%"/>
</p>

<p align="center">
  <em>Left: Orbital-resolved density of states with gradient fills. Right: Color-coded band structure with VBM at 0 eV.</em>
</p>

## 🚀 Quick Start

```bash
pip install valyte
```

Then, from a directory containing your VASP output files:

```bash
valyte dos                      # Plot density of states
valyte band                     # Plot band structure
valyte converge --forces        # Check relaxation convergence
```

That's it. See the [Getting Started guide](https://valyte-project.readthedocs.io/en/latest/getting-started/) for a complete walkthrough.

---

## Features

### Pre-processing

| Command | Description |
|---|---|
| `valyte supercell` | Generate supercells from POSCAR files |
| `valyte kpt` | Interactive or batch KPOINTS generation (Monkhorst-Pack / Gamma) |
| `valyte band kpt-gen` | Automatic high-symmetry k-path (Bradley–Cracknell by default) |
| `valyte potcar` | Generate POTCAR from POSCAR species |

### Post-processing

| Command | Description |
|---|---|
| `valyte dos` | Total and projected DOS with orbital resolution and gradient fills |
| `valyte band` | Color-coded band structure with VBM aligned to 0 eV |
| `valyte band --tricolor` | Orbital-resolved tricolor band structure |
| `valyte band --spin-resolved` | Spin-polarized band structure — spin-up and spin-down channels |
| `valyte band --spin-texture` | Non-collinear spin texture — bands colored by Sₓ, Sᵧ, or S_z |
| `valyte ipr` | Inverse Participation Ratio from PROCAR |
| `valyte effmass` | Carrier effective masses at VBM/CBM from parabolic fitting |
| `valyte converge` | Ionic and SCF convergence monitor with multi-panel plots |

---

## Installation

### From PyPI (recommended)

```bash
pip install valyte
```

### Update to the latest version

```bash
pip install --upgrade valyte
```

### From source (for development)

```bash
git clone https://github.com/nikyadav002/Valyte-Project
cd Valyte-Project
pip install -e .
```

### Requirements

- Python ≥ 3.8
- Dependencies (`numpy`, `matplotlib`, `pymatgen`, `scipy`, `click`, `seekpath`) are installed automatically.
- For `valyte potcar`: requires [pymatgen pseudopotential setup](https://pymatgen.org/installation.html#potcar-setup).

---

## Documentation

**[📖 Full Documentation → valyte-project.readthedocs.io](https://valyte-project.readthedocs.io/en/latest/)**

The documentation site includes:

| Page | What you'll find |
|---|---|
| [Getting Started](https://valyte-project.readthedocs.io/en/latest/getting-started/) | Installation, prerequisites, and your first plot |
| [Band Structure](https://valyte-project.readthedocs.io/en/latest/band/) | Standard, tricolor, spin-resolved, and spin-texture modes |
| [Density of States](https://valyte-project.readthedocs.io/en/latest/dos/) | Total and projected DOS with orbital resolution |
| [Effective Mass](https://valyte-project.readthedocs.io/en/latest/effmass/) | Carrier effective masses from parabolic fitting |
| [Convergence](https://valyte-project.readthedocs.io/en/latest/converge/) | Ionic, electronic, force, and pressure convergence |
| [IPR](https://valyte-project.readthedocs.io/en/latest/ipr/) | Wavefunction localization analysis |
| [Pre-processing](https://valyte-project.readthedocs.io/en/latest/preprocessing/) | Supercells, k-points, and POTCAR generation |
| [CLI Reference](https://valyte-project.readthedocs.io/en/latest/cli-reference/) | Every command and flag in one searchable page |
| [FAQ](https://valyte-project.readthedocs.io/en/latest/faq/) | Common issues and troubleshooting |

---

## Contributing

Contributions are welcome! Whether it's a bug report, feature request, or pull request — all help is appreciated.

- 🐛 **Report a bug** — [Open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new)
- 💡 **Request a feature** — [Open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new)
- 🔧 **Submit a fix** — Fork, branch, and [open a pull request](https://github.com/nikyadav002/Valyte-Project/pulls)

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

---

## License

This project is licensed under the [MIT License](LICENSE).

---

<p align="center">
  Built by <a href="https://github.com/nikyadav002">Nikhil Singh</a>
</p>
