<p align="center">
  <img src="valyte/Logo.png" alt="Valyte Logo" width="100%"/>
</p>

<p align="center">
  <a href="https://pypi.org/project/valyte/"><img src="https://img.shields.io/pypi/v/valyte?color=7c3aed&label=PyPI" alt="PyPI Version"></a>
  <a href="https://pypi.org/project/valyte/"><img src="https://img.shields.io/pypi/pyversions/valyte?color=7c3aed" alt="Python Versions"></a>
  <a href="https://valyte.readthedocs.io/"><img src="https://readthedocs.org/projects/valyte/badge/?version=latest" alt="Docs"></a>
  <a href="https://github.com/nikyadav002/Valyte-Project/blob/main/LICENSE"><img src="https://img.shields.io/github/license/nikyadav002/Valyte-Project?color=2a9d8f" alt="License"></a>
</p>

# Valyte

**Valyte** is a CLI tool for VASP workflows — pre-processing and post-processing — built for clean, publication-quality output.

**[Full Documentation → valyte.readthedocs.io](https://valyte.readthedocs.io/)**

---

## Features

### Pre-processing

| Command | Description |
|---|---|
| `valyte supercell` | Generate supercells from POSCAR files |
| `valyte kpt` | Interactive KPOINTS generation (Monkhorst-Pack / Gamma) |
| `valyte band kpt-gen` | Automatic high-symmetry k-path (Bradley-Cracknell by default) |
| `valyte potcar` | Generate POTCAR from POSCAR species |

### Post-processing

| Command | Description |
|---|---|
| `valyte dos` | Total and projected DOS with orbital resolution and gradient fills |
| `valyte band` | Color-coded band structure with VBM aligned to 0 eV |
| `valyte band --tricolor` | **Orbital-resolved tricolor band structure** |
| `valyte ipr` | Inverse Participation Ratio from PROCAR |

---

## Installation

```bash
pip install valyte
```

To update:

```bash
pip install --upgrade valyte
```

Or from source:

```bash
git clone https://github.com/nikyadav002/Valyte-Project
cd Valyte-Project
pip install -e .
```

---

## Gallery

<p align="center">
  <img src="valyte/valyte_dos.png" alt="DOS Plot Example" width="47%"/>
  <img src="valyte/valyte_band.png" alt="Band Structure Example" width="38%"/>
</p>

---

## Usage

Run `valyte --help` or `valyte <command> --help` at any time.

<details>
<summary><strong>Supercell</strong></summary>

<br>

```bash
valyte supercell nx ny nz [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `POSCAR_supercell` | Output filename |

```bash
valyte supercell 2 2 2
valyte supercell 3 3 1 -i POSCAR_primitive -o POSCAR_3x3x1
```

</details>

<details>
<summary><strong>Band Structure — KPOINTS generation</strong></summary>

<br>

Automatically generate a line-mode KPOINTS file with high-symmetry paths.

```bash
valyte band kpt-gen [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-n`, `--npoints` | `40` | Points per segment |
| `-o`, `--output` | `KPOINTS` | Output filename |
| `--mode` | `bradcrack` | Path convention: `bradcrack`, `seekpath`, `latimer_munro`, `setyawan_curtarolo` |

```bash
valyte band kpt-gen -n 60
valyte band kpt-gen --mode seekpath
```

> **Important:** A `POSCAR_standard` file is also written. You **must** use this standardized structure for the band calculation (`cp POSCAR_standard POSCAR`) — the k-path corresponds to this specific cell orientation.

</details>

<details>
<summary><strong>Band Structure — Standard plot</strong></summary>

<br>

Bands are colored purple (valence) and teal (conduction), with the VBM set to 0 eV.

```bash
valyte band [options]
```

| Option | Default | Description |
|---|---|---|
| `--vasprun` | `.` | Path to `vasprun.xml` or directory |
| `--kpoints` | auto-detected | Path to `KPOINTS` file for labels |
| `--ylim` | `-4 4` | Energy window |
| `-o`, `--output` | `valyte_band.png` | Output filename |
| `--font` | `Arial` | Font: `Arial`, `Helvetica`, `Times New Roman` |

```bash
valyte band --ylim -3 3 -o my_bands.png
```

</details>

<details>
<summary><strong>Band Structure — Tricolor orbital-resolved plot</strong></summary>

<br>

Each band segment is colored by blending three base colors weighted by the relative orbital/element projection at each k-point. A ternary triangle legend is drawn in the corner.

```bash
valyte band --tricolor SPEC1 SPEC2 SPEC3 [options]
```

**Spec formats:**

| Format | Example | Selects |
|---|---|---|
| Orbital | `s`, `p`, `d`, `f` | That orbital across all atoms |
| Element | `Fe`, `O` | All orbitals for that element |
| Element + orbital | `Fe:d`, `O(p)` | Specific orbital for that element |

| Option | Default | Description |
|---|---|---|
| `--tricolor` | — | 3 specs (required) |
| `--tricolors` | `#e74c3c #2ecc71 #3498db` | 3 colors (red, green, blue) |
| `--tri-labels` | spec strings | 3 labels for the triangle legend |
| `--lw` | `2.0` | Line width |

> **Requirement:** VASP must be run with `LORBIT = 11` (or ≥ 10) to write projected eigenvalues into `vasprun.xml`.

```bash
# s / p / d — default red, green, blue
valyte band --tricolor s p d --ylim -4 4

# Element-resolved (e.g. MoSSe)
valyte band --tricolor Mo S Se \
  --tricolors "#e74c3c" "#2ecc71" "#3498db" \
  --tri-labels Mo S Se --ylim -3 3

# Element + orbital resolved
valyte band --tricolor Fe:d O:p s --ylim -5 5 -o orbital_band.png
```

</details>

<details>
<summary><strong>DOS — Density of States</strong></summary>

<br>

```bash
valyte dos [path/to/vasprun.xml] [options]
```

| Option | Default | Description |
|---|---|---|
| `-e`, `--elements` | all | Elements/orbitals to plot |
| `--xlim` | `-6 6` | Energy range |
| `--ylim` | auto | DOS range |
| `--scale` | `1.0` | Divide DOS by this factor |
| `--fermi` | off | Draw dashed line at E = 0 |
| `--pdos` | off | Show only projected DOS |
| `--legend-cutoff` | `0.10` | Hide legend if PDOS fraction < threshold |
| `-o`, `--output` | `valyte_dos.png` | Output filename |
| `--font` | `Arial` | Font family |

```bash
valyte dos                            # All orbitals, all elements
valyte dos -e Fe O                    # Total PDOS for Fe and O
valyte dos -e "Fe(d)" "O(p)"          # Specific orbitals
valyte dos -e Fe "Fe(d)"              # Fe total + Fe d-orbital
valyte dos ./run --xlim -5 5 --fermi -o my_dos.png
```

</details>

<details>
<summary><strong>K-Points — Interactive SCF grid</strong></summary>

<br>

```bash
valyte kpt
```

Prompts for K-mesh scheme (Monkhorst-Pack or Gamma) and K-spacing in 2π/Å, then calculates the optimal grid from your `POSCAR`.

</details>

<details>
<summary><strong>POTCAR</strong></summary>

<br>

```bash
valyte potcar [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `POTCAR` | Output filename |
| `--functional` | `PBE` | `PBE`, `PBE_52`, `PBE_54`, `LDA`, etc. |

```bash
valyte potcar
valyte potcar --functional PBE_54
valyte potcar -i POSCAR_relaxed -o POTCAR_new
```

> **Pymatgen configuration required:** Set `PMG_VASP_PSP_DIR` in `~/.pmgrc.yaml`. See the [Pymatgen docs](https://pymatgen.org/installation.html#potcar-setup).

</details>

<details>
<summary><strong>IPR — Inverse Participation Ratio</strong></summary>

<br>

Compute the Inverse Participation Ratio from `PROCAR` to quantify wavefunction localization.

```bash
valyte ipr
```

Interactive — reads `PROCAR`, shows system info, prompts for band indices, saves results to `ipr_procar.dat`.

**Output columns:** Band | Energy (eV) | IPR | N_eff (= 1/IPR)

A state localized on a single atom has IPR ≈ 1 and N_eff ≈ 1. Delocalized band states have small IPR and large N_eff. Use IPR to identify defect states in supercell calculations.

</details>
