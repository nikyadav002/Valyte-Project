# Getting Started

This guide will take you from installation to your first publication-quality plot in under two minutes.

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

---

## Prerequisites

### Python

Valyte requires **Python 3.8 or later**. All Python dependencies are installed automatically:

| Package | Purpose |
|---|---|
| `numpy` | Numerical operations |
| `matplotlib` | Plotting |
| `pymatgen` | Crystal structure parsing |
| `scipy` | Curve fitting |
| `click` | CLI framework |
| `seekpath` | K-path generation |

### VASP output files

Valyte reads standard VASP output files. Here's what each command needs:

| Command | Required files |
|---|---|
| `valyte dos` | `vasprun.xml` |
| `valyte band` | `vasprun.xml`, `KPOINTS` |
| `valyte converge` | `OSZICAR` (optionally `OUTCAR`) |
| `valyte ipr` | `PROCAR` |
| `valyte effmass` | `vasprun.xml`, `KPOINTS` |

### Pymatgen setup (for `valyte potcar` only)

The `valyte potcar` command requires access to VASP pseudopotential files through pymatgen. Set this up by adding the following to `~/.pmgrc.yaml`:

```yaml
PMG_VASP_PSP_DIR: /path/to/your/pseudopotentials
```

See the [pymatgen documentation](https://pymatgen.org/installation.html#potcar-setup) for full setup instructions.

!!! note
    This setup is only required for `valyte potcar`. All other commands work without pymatgen configuration.

---

## Verify your installation

```bash
valyte --help
```

You should see a list of available commands:

```
Usage: valyte [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  band       Band structure plotting and KPOINTS generation
  converge   Convergence monitoring
  dos        Density of states plotting
  effmass    Effective mass calculation
  ipr        Inverse Participation Ratio
  kpt        Interactive K-points generation
  potcar     POTCAR generation
  supercell  Supercell generation
```

---

## Your first plot

### Density of states

Navigate to a directory containing a `vasprun.xml` from a completed VASP calculation:

```bash
cd /path/to/your/vasp/run
valyte dos
```

This produces `valyte_dos.png` — a publication-quality DOS plot with orbital resolution and gradient fills. No flags needed.

### Band structure

For a band structure calculation (line-mode `KPOINTS`):

```bash
valyte band
```

This produces `valyte_band.png` with the VBM set to 0 eV, valence bands in purple, and conduction bands in teal.

### Convergence check

For a relaxation or single-point calculation:

```bash
valyte converge
```

This reads `OSZICAR` and produces a multi-panel convergence plot showing energy and ΔE across ionic steps.

---

## Common workflows

### "I just finished a relaxation"

```bash
valyte converge --forces              # Check energy + force convergence
valyte converge --forces --stress     # Also check pressure
valyte converge --no-plot             # Quick terminal summary (no plot)
```

### "I have a band structure calculation"

```bash
valyte band                           # Standard plot
valyte band --tricolor s p d          # Orbital-resolved (needs LORBIT ≥ 10)
valyte band --spin-resolved           # Spin-polarized
valyte effmass                        # Effective masses at VBM/CBM
```

### "I need to set up a new calculation"

```bash
valyte supercell 2 2 2                # Generate a 2×2×2 supercell
valyte kpt                            # Interactive k-point grid
valyte band kpt-gen                   # High-symmetry k-path for bands
valyte potcar                         # Generate POTCAR
```

### "I want to export the data for custom plotting"

Every plot command supports `--save-data`:

```bash
valyte dos --save-data                # → valyte_dos.dat
valyte band --save-data               # → valyte_band.dat
valyte effmass --save-data            # → valyte_effmass.dat
valyte converge --save-data           # → valyte_converge.dat
```

---

## Next steps

- **[Band Structure](band.md)** — All five band structure modes explained
- **[DOS](dos.md)** — Element and orbital selection, smart features
- **[Effective Mass](effmass.md)** — Parabolic fitting for carrier masses
- **[Convergence](converge.md)** — Multi-panel convergence diagnostics
- **[IPR](ipr.md)** — Wavefunction localization analysis
- **[Pre-processing](preprocessing.md)** — Supercells, k-points, and POTCAR
- **[CLI Reference](cli-reference.md)** — Every command and flag in one page
- **[FAQ](faq.md)** — Common issues and troubleshooting
