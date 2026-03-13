# Band Structure

Valyte provides several band structure modes:

| Sub-command / flag | Description |
|---|---|
| `valyte band kpt-gen` | Generate KPOINTS with high-symmetry paths |
| `valyte band` | Standard color-coded band structure plot |
| `valyte band --tricolor` | Orbital-resolved tricolor band structure |
| `valyte band --spin-resolved` | Spin-polarized plot — spin-up (blue) and spin-down (red) |
| `valyte band --spin-texture` | Non-collinear spin texture colored by Sₓ, S_y, or S_z |

---

## 1. Generate KPOINTS

Automatically generate a line-mode KPOINTS file with high-symmetry paths for your structure.

```bash
valyte band kpt-gen [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-n`, `--npoints` | `40` | K-points per segment |
| `-o`, `--output` | `KPOINTS` | Output filename |
| `--mode` | `bradcrack` | Path convention: `bradcrack`, `seekpath`, `latimer_munro`, `setyawan_curtarolo` |
| `--symprec` | `0.01` | Symmetry precision |

!!! tip "Smart K-Path Generation"
    Valyte uses the **Bradley-Cracknell** convention by default, automatically determining the correct high-symmetry path for your Bravais lattice (e.g., Γ–Y–V for monoclinic, Γ–X–M–Γ–R for cubic) without any external dependencies.

```bash
# Default — Bradley-Cracknell, 40 points per segment
valyte band kpt-gen

# Denser path
valyte band kpt-gen -n 60

# Use Seekpath convention instead
valyte band kpt-gen --mode seekpath
```

!!! important "Use POSCAR_standard for the calculation"
    The command also writes a `POSCAR_standard` file. You **must** use this standardized primitive cell for your VASP calculation — the k-path is defined in its reciprocal space.

    ```bash
    cp POSCAR_standard POSCAR
    ```

    Running the band calculation with your original POSCAR will produce incorrect k-point labels and potentially wrong paths.

!!! note
    A `POTCAR` is also automatically generated (PBE) after KPOINTS creation.

---

## 2. Standard Band Structure Plot

Plot the electronic band structure from `vasprun.xml`. Bands below the Fermi level (VBM set to 0 eV) are colored purple; conduction bands are teal.

```bash
valyte band [options]
```

| Option | Default | Description |
|---|---|---|
| `--vasprun` | `.` | Path to `vasprun.xml` or its parent directory |
| `--kpoints` | auto-detected | Path to `KPOINTS` file (for high-symmetry labels) |
| `--ylim` | `-4 4` | Energy window in eV, e.g. `--ylim -3 3` |
| `-o`, `--output` | `valyte_band.png` | Output filename |
| `--font` | `Arial` | Font family: `Arial`, `Helvetica`, `Times New Roman` |
| `--save-data` | off | Save band data to `valyte_band.dat` |

```bash
# Basic plot from current directory
valyte band

# Custom energy window and output
valyte band --ylim -3 3 -o my_bands.png

# Explicit paths
valyte band --vasprun ./band_run --kpoints ./band_run/KPOINTS --ylim -4 4

# Save raw data alongside the plot
valyte band --ylim -3 3 --save-data
```

### Exported data format (`valyte_band.dat`)

A plain-text, whitespace-delimited file with one row per k-point:

```
# K-point labels: G=0.0000, X=0.3124, M=0.5412, G=0.8660
# k-dist  band_1  band_2  band_3  ...
0.000000  -5.123  -3.456  -1.234  ...
...
```

- **`k-dist`** — cumulative k-path distance (same x-axis as the plot)
- **`band_N`** — energy in eV relative to the VBM (E = 0)
- For spin-polarized calculations columns are labeled `band_N_up` / `band_N_dn`
- The header comment lists every high-symmetry label and its k-distance position

---

## 3. Tricolor Orbital-Resolved Plot

Plot an orbital- or element-resolved band structure where each band segment is continuously colored by blending three base colors according to the relative projection weights at each k-point. A ternary triangle legend in the corner maps vertex colors to orbital contributions.

```bash
valyte band --tricolor SPEC1 SPEC2 SPEC3 [options]
```

!!! important "VASP requirement"
    Run VASP with `LORBIT = 11` (or ≥ 10) in your `INCAR`. This writes orbital projection data into `vasprun.xml`, which is required for this mode.

### Spec formats

| Format | Example | Selects |
|---|---|---|
| Orbital type | `s`, `p`, `d`, `f` | That orbital across all atoms |
| Element | `Fe`, `O` | All orbitals summed over that element |
| Element + orbital | `Fe:d`, `O(p)` | Specific orbital for a specific element |

### Options

| Option | Default | Description |
|---|---|---|
| `--tricolor` | — | 3 specs — **required** to activate this mode |
| `--tricolors` | `#e74c3c #2ecc71 #3498db` | 3 matplotlib colors (red, green, blue) |
| `--tri-labels` | spec strings | Labels shown at the triangle legend vertices |
| `--lw` | `2.0` | Band line width |
| `--ylim`, `-o`, `--font` | — | Same as standard band options |
| `--save-data` | off | Save band data to `valyte_band.dat` (same format as standard plot) |

### How the colors work

For each k-point and each band, the three projection weights are normalized to sum to 1:

$$w_1 + w_2 + w_3 = 1$$

The RGB color is then computed as a linear blend:

$$\text{color} = w_1 \cdot c_1 + w_2 \cdot c_2 + w_3 \cdot c_3$$

A segment with pure $w_1 = 1$ shows color $c_1$; mixed contributions produce intermediate colors. The triangle legend visualizes this mapping.

### Examples

```bash
# s / p / d contributions — default red, green, blue
valyte band --tricolor s p d --ylim -4 4

# Element-resolved for a MoSSe heterostructure
valyte band --tricolor Mo S Se \
  --tricolors "#e74c3c" "#2ecc71" "#3498db" \
  --tri-labels Mo S Se \
  --ylim -3 3 -o band_tricolor.png

# Element + orbital resolved
valyte band --tricolor Fe:d O:p s --ylim -5 5 -o orbital_band.png

# Custom line width and font
valyte band --tricolor s p d --lw 2.5 --font Helvetica --ylim -4 4
```

---

## 4. Spin-Resolved Band Structure (collinear)

For spin-polarized calculations, plot both spin channels on the same axes with distinct colors:

- **Spin-up** — solid blue line
- **Spin-down** — dashed red line

A legend is added automatically. If the calculation is non-magnetic (single spin channel), Valyte prints a warning and falls back to the standard plot.

```bash
valyte band --spin-resolved [options]
```

### Options

| Option | Default | Description |
|---|---|---|
| `--spin-resolved` | off | Enable spin-resolved mode |
| `--vasprun` | `.` | Path to `vasprun.xml` or parent directory |
| `--kpoints` | auto-detected | Path to `KPOINTS` file (for labels) |
| `--ylim` | `-4 4` | Energy window in eV |
| `-o`, `--output` | `valyte_band.png` | Output filename |
| `--font` | `Arial` | Font family |
| `--width` / `--height` | `4.0` | Figure dimensions in inches |
| `--save-data` | off | Save band data to `valyte_band.dat` |

### Examples

```bash
# Basic spin-resolved plot
valyte band --spin-resolved

# Custom energy window
valyte band --spin-resolved --ylim -3 3

# With explicit path and output
valyte band --spin-resolved --vasprun ./band_run --ylim -4 4 -o spin_bands.png
```

---

## 5. Non-Collinear Spin Texture

For non-collinear calculations (spin-orbit coupling or non-collinear magnetism), color each band segment by the expectation value of a spin component. The color encodes the spin polarization direction along the k-path using a diverging colormap centered at zero.

```bash
valyte band --spin-texture {sx,sy,sz} [options]
```

!!! important "VASP requirements"
    - `LSORBIT = .TRUE.` or `LNONCOLLINEAR = .TRUE.` in `INCAR`
    - `LORBIT >= 11` to write projected eigenvalues and magnetization into `vasprun.xml`

### Options

| Option | Default | Description |
|---|---|---|
| `--spin-texture` | — | Component: `sx`, `sy`, or `sz` — **required** to activate this mode |
| `--spin-cmap` | `seismic` | Diverging colormap (any matplotlib colormap name) |
| `--vasprun` | `.` | Path to `vasprun.xml` or parent directory |
| `--kpoints` | auto-detected | Path to `KPOINTS` file (for labels) |
| `--ylim` | `-4 4` | Energy window in eV |
| `-o`, `--output` | `valyte_band_sz.png` | Output filename (auto-named by component if not set) |
| `--font` | `Arial` | Font family |
| `--width` / `--height` | `4.0` | Figure dimensions in inches |
| `--lw` | `2.0` | Band line width |
| `--save-data` | off | Save band energies to `valyte_band.dat` |

### How it works

The `projected_magnetisation` array from `vasprun.xml` has shape `(nkpts, nbands, natoms, norbitals, 3)`. The chosen component (index 0, 1, or 2 for x, y, z) is summed over all atoms and orbitals to give a scalar spin texture value per k-point per band:

$$\langle S_\alpha \rangle (k, n) = \sum_{i,m} M_{\alpha,km,in}$$

This is mapped to color via a diverging colormap (positive = one color, negative = the opposite, zero = white).

### Examples

```bash
# Out-of-plane spin texture (most common for Rashba / topological systems)
valyte band --spin-texture sz

# In-plane components
valyte band --spin-texture sx
valyte band --spin-texture sy

# Custom colormap and energy window
valyte band --spin-texture sz --spin-cmap RdBu_r --ylim -2 2

# All three components in one go
for comp in sx sy sz; do
    valyte band --spin-texture $comp --ylim -2 2
done
```
