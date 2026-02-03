<p align="center">
  <img src="valyte/Logo.png" alt="Valyte Logo" width="100%"/>
</p>

# Valyte

**Valyte** is a comprehensive CLI tool for VASP workflows, providing both pre-processing and post-processing capabilities with a focus on clean, publication-quality outputs and modern aesthetics.

## Features

### Pre-processing
- **Supercell Creation**: Generate supercells from POSCAR files.
- **Interactive K-Point Generation**: Create KPOINTS files with automatic grid calculation based on K-spacing.
- **Band KPOINTS Generation**: Automatic high-symmetry path detection for band structure calculations.

### Post-processing
- **DOS Plotting**: 
  - Smart Plotting: Automatically handles total DOS and Projected DOS (PDOS).
  - Orbital-Resolved: Plots individual orbitals (s, p, d, f) by default.
  - Adaptive Legend: Intelligently hides the legend if PDOS contributions are low.
  - Gradient Fill: Aesthetically pleasing gradient fills for DOS peaks.
- **Band Structure Plotting**:
  - VBM alignment to 0 eV.
  - Color-coded bands (Purple for VB, Teal for CB).
  - High-symmetry path labels from KPOINTS.
- **IPR Analysis**: Compute Inverse Participation Ratio from PROCAR to analyze wavefunction localization.
- **Publication Quality**: Clean aesthetics, custom fonts (Arial, Helvetica, Times New Roman), high DPI output.

## Installation

Install Valyte directly from PyPI:

```bash
pip install valyte
```

Or install from source:

```bash
git clone https://github.com/nikyadav002/Valyte-Project
cd Valyte-Project
pip install -e .
```

## Examples

<p align="center">
  <img src="valyte/valyte_dos.png" alt="DOS Plot Example" width="47%"/>
  <img src="valyte/valyte_band.png" alt="Band Structure Example" width="38%"/>
</p>

## Updating Valyte

To update to the latest version:

```bash
pip install --upgrade valyte
```

## Usage

The main command is `valyte`.

<details>
<summary><strong>Click to view detailed usage instructions</strong></summary>

<br>

### 🧊 Create Supercell

Generate a supercell from a POSCAR file:

```bash
valyte supercell nx ny nz [options]
```

**Example:**
```bash
# Create a 2×2×2 supercell
valyte supercell 2 2 2

# Specify input and output files
valyte supercell 3 3 1 -i POSCAR_primitive -o POSCAR_3x3x1
```

**Options:**
- `-i`, `--input`: Input POSCAR file (default: `POSCAR`).
- `-o`, `--output`: Output filename (default: `POSCAR_supercell`).

---

### 📉 Band Structure

#### 1. Generate KPOINTS

Automatically generate a KPOINTS file with high-symmetry paths for band structure calculations.

> [!TIP]
> **Smart K-Path Generation (New in v0.1.7+)**: Valyte now automatically determines the standard path (e.g., `\Gamma - Y - V` for Monoclinic cells) using the **Bradley-Cracknell** convention by default. This ensures clean, publication-ready labels without external dependencies.

```bash
valyte band kpt-gen [options]
```

**Options:**
- `-i`, `--input`: Input POSCAR file (default: `POSCAR`).
- `-n`, `--npoints`: Points per segment (default: `40`).
- `-o`, `--output`: Output filename (default: `KPOINTS`).
- `--mode`: Path convention. Options: `bradcrack` (Default), `seekpath`, `latimer_munro`, `setyawan_curtarolo`.

**Example:**
```bash
# Default (Smart/BradCrack)
valyte band kpt-gen -n 60

# Explicitly use Seekpath convention
valyte band kpt-gen --mode seekpath
```

> [!IMPORTANT]
> The command will generate a **`POSCAR_standard`** file. You **MUST** use this structure for your band structure calculation (i.e., `cp POSCAR_standard POSCAR`) because the K-path corresponds to this specific orientation. Using your original POSCAR may result in incorrect paths.

### 🕸️ Generate K-Points (Interactive)

Generate a `KPOINTS` file for SCF calculations interactively.

```bash
valyte kpt
```

This command will prompt you for:
1. **K-Mesh Scheme**: Monkhorst-Pack or Gamma.
2. **K-Spacing**: Value in $2\pi/\AA$ (e.g., 0.04).

It automatically calculates the optimal grid based on your `POSCAR` structure.

---

### 🧪 Generate POTCAR

Generate a `POTCAR` file based on the species in your `POSCAR`.

```bash
valyte potcar [options]
```

**Options:**
- `-i`, `--input`: Input POSCAR file (default: `POSCAR`).
- `-o`, `--output`: Output filename (default: `POTCAR`).
- `--functional`: Functional to use (default: `PBE`). Options include `PBE`, `PBE_52`, `PBE_54`, `LDA`, etc.

**Example:**
```bash
# Generate POTCAR using default PBE functional
valyte potcar

# Use a specific functional
valyte potcar --functional PBE_54

# Specify input and output files
valyte potcar -i POSCAR_relaxed -o POTCAR_new
```

> [!IMPORTANT]
> **Pymatgen Configuration Required**: This command requires Pymatgen to be configured with your VASP pseudopotential directory. Set `PMG_VASP_PSP_DIR` in your `~/.pmgrc.yaml` file. See the [Pymatgen documentation](https://pymatgen.org/installation.html#potcar-setup) for setup instructions.

---

#### 2. Plot Band Structure

Plot the electronic band structure from `vasprun.xml`.

```bash
valyte band [options]
```

**Options:**
- `--vasprun`: Path to `vasprun.xml` (default: current directory).
- `--kpoints`: Path to `KPOINTS` file for path labels (default: looks for `KPOINTS` in same dir).
- `--ylim`: Energy range, e.g., `--ylim -4 4`.
- `-o, --output`: Output filename (default: `valyte_band.png`).

**Example:**
```bash
valyte band --ylim -3 3 -o my_bands.png
```

---

### 📊 Plot DOS

```bash
valyte dos [path/to/vasprun.xml] [options]
```

You can provide the path as a positional argument, use the `--vasprun` flag, or omit it to use the current directory.

**Examples:**
```bash
# Plot all orbitals for all elements (Default)
valyte dos

# Plot specific elements (Total PDOS)
valyte dos -e Fe O

# Plot specific orbitals
valyte dos -e "Fe(d)" "O(p)"

# Plot mixed (Fe Total and Fe d-orbital)
valyte dos -e Fe "Fe(d)"
```

**Options:**
- `-e`, `--elements`: Specific elements or orbitals to plot.
    - Example: `-e Fe O` (Plots Total PDOS for Fe and O).
    - Example: `-e Fe(d) O(p)` (Plots Fe d-orbital and O p-orbital).
    - Example: `-e Fe Fe(d)` (Plots Fe Total and Fe d-orbital).
- `--xlim`: Energy range (default: `-6 6`).
- `--ylim`: DOS range (e.g., `--ylim 0 10`).
- `--scale`: Scaling factor for Y-axis (e.g., `--scale 3` divides DOS by 3).
- `--fermi`: Draw a dashed line at the Fermi level (E=0). Default is OFF.
- `--pdos`: Plot only Projected DOS (hide Total DOS).
- `--legend-cutoff`: Threshold for legend visibility (default: `0.10` = 10%).
- `-o`, `--output`: Output filename (default: `valyte_dos.png`).
- `--font`: Font family (default: `Arial`).

**Example:**

```bash
valyte dos ./vasp_data --xlim -5 5 -o my_dos.png
```

---

### 📐 Compute IPR (Inverse Participation Ratio)

Compute the Inverse Participation Ratio (IPR) from VASP `PROCAR` to analyze wavefunction localization.

```bash
valyte ipr
```

This interactive command will:
1. Read the `PROCAR` file from the current directory.
2. Display the number of k-points, bands, and atoms.
3. Prompt you for **band indices** to analyze (e.g., `5 6 7` or `5-7`).
4. Optionally show per-k-point IPR values.
5. Save results to `ipr_procar.dat`.

**Output Columns:**
- **Band**: Band index.
- **Energy**: Average energy (eV) across k-points.
- **IPR**: Inverse Participation Ratio (higher = more localized).
- **N_eff**: Effective number of atoms (1/IPR).

> [!TIP]
> Use IPR to identify localized defect states. A state localized on a single atom has IPR ≈ 1.0 and N_eff ≈ 1. Delocalized band states have small IPR and large N_eff.

</details>
