<p align="center">
  <img src="Logo.png" alt="Valyte Logo" width="100%"/>
</p>

# Valyte

**Valyte** is a CLI tool for VASP workflows — pre-processing and post-processing — built for clean, publication-quality output.

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
  <img src="valyte_dos.png" alt="DOS Plot Example" width="47%"/>
  <img src="valyte_band.png" alt="Band Structure Example" width="38%"/>
</p>

---

## Usage

### Supercell

Generate a supercell from a POSCAR file:

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

---

### Band Structure

#### 1. Generate KPOINTS

Automatically generate a line-mode KPOINTS file with high-symmetry paths.

```bash
valyte band kpt-gen [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-n`, `--npoints` | `40` | Points per segment |
| `-o`, `--output` | `KPOINTS` | Output filename |
| `--mode` | `bradcrack` | Path convention |

!!! tip "Smart K-Path Generation"
    Valyte uses the **Bradley-Cracknell** convention by default, automatically determining the correct high-symmetry path (e.g., Γ–Y–V for monoclinic cells) without external dependencies.

```bash
valyte band kpt-gen -n 60
valyte band kpt-gen --mode seekpath
```

!!! important
    The command also writes `POSCAR_standard`. You **must** use this standardized structure for the band calculation:

    ```bash
    cp POSCAR_standard POSCAR
    ```

    The k-path corresponds to this specific cell orientation — using your original POSCAR will give incorrect paths.

---

#### 2. Standard Band Structure Plot

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
| `--font` | `Arial` | `Arial`, `Helvetica`, `Times New Roman` |

```bash
valyte band --ylim -3 3 -o my_bands.png
```

---

#### 3. Tricolor Orbital-Resolved Band Structure

Each band segment is colored by blending three base colors weighted by the relative orbital/element projection at that k-point. A ternary triangle legend is drawn in the corner, showing which vertex corresponds to which orbital.

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
| `--tricolor` | — | 3 specs (required to activate this mode) |
| `--tricolors` | `#e74c3c #2ecc71 #3498db` | 3 colors (red, green, blue) |
| `--tri-labels` | spec strings | 3 labels for the triangle legend |
| `--lw` | `2.0` | Line width |

!!! important "VASP requirement"
    Run VASP with `LORBIT = 11` (or ≥ 10) in your `INCAR` so that `vasprun.xml` contains projected eigenvalues.

```bash
# s / p / d orbital contributions (default red, green, blue)
valyte band --tricolor s p d --ylim -4 4

# Element-resolved (e.g. MoSSe heterostructure)
valyte band --tricolor Mo S Se \
  --tricolors "#e74c3c" "#2ecc71" "#3498db" \
  --tri-labels Mo S Se --ylim -3 3

# Element + orbital resolved
valyte band --tricolor Fe:d O:p s --ylim -5 5 -o orbital_band.png
```

---

### DOS — Density of States

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

---

### K-Points — Interactive SCF grid

```bash
valyte kpt
```

Prompts for K-mesh scheme (Monkhorst-Pack or Gamma) and K-spacing in 2π/Å, then calculates the optimal grid from your `POSCAR`.

---

### POTCAR

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

!!! important "Pymatgen configuration required"
    Set `PMG_VASP_PSP_DIR` in `~/.pmgrc.yaml` to point to your VASP pseudopotential directory.
    See the [Pymatgen documentation](https://pymatgen.org/installation.html#potcar-setup) for setup instructions.

---

### IPR — Inverse Participation Ratio

Compute the Inverse Participation Ratio from `PROCAR` to quantify wavefunction localization.

```bash
valyte ipr
```

Interactive workflow:

1. Reads `PROCAR` from the current directory
2. Displays k-points, bands, and atom count
3. Prompts for band indices (e.g. `5 6 7` or `5-7`)
4. Optionally shows per-k-point IPR values
5. Saves results to `ipr_procar.dat`

**Output columns:** Band | Energy (eV) | IPR | N_eff (= 1/IPR)

!!! tip
    Use IPR to identify localized defect states in supercell calculations. A state localized on a single atom has IPR ≈ 1 and N_eff ≈ 1. Delocalized band states have small IPR and large N_eff.
