# CLI Reference

A complete reference for every Valyte command and flag. Use `valyte --help` or `valyte <command> --help` at any time for quick help.

---

## Global options

```bash
valyte [OPTIONS] COMMAND [ARGS]...
```

| Option | Description |
|---|---|
| `--version` | Show the Valyte version and exit |
| `--help` | Show available commands and exit |

---

## `valyte supercell`

Generate a supercell from a POSCAR file.

```bash
valyte supercell nx ny nz [options]
```

| Argument / Option | Default | Description |
|---|---|---|
| `nx ny nz` | *(required)* | Supercell dimensions along a, b, c |
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `POSCAR_supercell` | Output filename |

**Examples:**

```bash
valyte supercell 2 2 2
valyte supercell 3 3 1 -i POSCAR_primitive -o POSCAR_3x3x1
```

---

## `valyte kpt`

Generate KPOINTS for SCF or relaxation calculations.

```bash
valyte kpt [options]
```

With no options, `valyte kpt` runs interactively and prompts for:

1. K-mesh scheme — Monkhorst-Pack or Gamma
2. K-spacing in 2π/Å (e.g., `0.04`)

The optimal grid is calculated from the lattice vectors in your `POSCAR`.

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `KPOINTS` | Output KPOINTS file |
| `--spacing` | `0.04` | K-spacing in 2π/Å |
| `--scheme` | `gamma` | K-mesh scheme: `gamma` or `mp` |
| `--no-potcar` | off | Do not auto-generate POTCAR if missing |

**Examples:**

```bash
valyte kpt
valyte kpt --spacing 0.04 --scheme gamma
valyte kpt -i POSCAR_relaxed -o KPOINTS_scf --spacing 0.03 --scheme mp
valyte kpt --spacing 0.05 --no-potcar
```

---

## `valyte potcar`

Generate a POTCAR file from the species listed in your POSCAR.

```bash
valyte potcar [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `POTCAR` | Output filename |
| `--functional` | `PBE` | Pseudopotential functional: `PBE`, `PBE_52`, `PBE_54`, `LDA`, `LDA_52`, `LDA_54`, `PW91`, etc. |

!!! important "Pymatgen configuration required"
    Set `PMG_VASP_PSP_DIR` in `~/.pmgrc.yaml`. See [pymatgen docs](https://pymatgen.org/installation.html#potcar-setup).

**Examples:**

```bash
valyte potcar
valyte potcar --functional PBE_54
valyte potcar -i POSCAR_relaxed -o POTCAR_new
```

---

## `valyte band kpt-gen`

Generate a line-mode KPOINTS file with high-symmetry paths.

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

!!! important "Use POSCAR_standard for the calculation"
    A `POSCAR_standard` file is also written. You **must** use this standardized cell for your band calculation — the k-path corresponds to this specific cell orientation.

**Examples:**

```bash
valyte band kpt-gen
valyte band kpt-gen -n 60
valyte band kpt-gen --mode seekpath
```

---

## `valyte band`

Plot the electronic band structure from `vasprun.xml`.

```bash
valyte band [options]
```

### Common options (all band modes)

| Option | Default | Description |
|---|---|---|
| `--vasprun` | `.` | Path to `vasprun.xml` or its parent directory |
| `--kpoints` | auto-detected | Path to `KPOINTS` file for high-symmetry labels |
| `--ylim` | `-4 4` | Energy window in eV |
| `-o`, `--output` | `valyte_band.png` | Output filename |
| `--font` | `Arial` | Font family: `Arial`, `Helvetica`, `Times New Roman` |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--save-data` | off | Save band data to `valyte_band.dat` |

### Tricolor mode

```bash
valyte band --tricolor SPEC1 SPEC2 SPEC3 [options]
```

| Option | Default | Description |
|---|---|---|
| `--tricolor` | — | 3 projection specs (required) |
| `--tricolors` | `#e74c3c #2ecc71 #3498db` | 3 colors for the triangle legend |
| `--tri-labels` | spec strings | Labels at triangle legend vertices |
| `--lw` | `2.0` | Band line width |

Spec formats: `s`, `p`, `d`, `f` (orbital), `Fe`, `O` (element), `Fe:d`, `O(p)` (element + orbital).

!!! note "Requires `LORBIT = 11` (or ≥ 10) in INCAR."

### Spin-resolved mode

```bash
valyte band --spin-resolved [options]
```

| Option | Default | Description |
|---|---|---|
| `--spin-resolved` | off | Enable spin-resolved mode |

Spin-up: solid blue. Spin-down: dashed red. Falls back to standard plot for non-magnetic calculations.

### Spin-texture mode

```bash
valyte band --spin-texture {sx,sy,sz} [options]
```

| Option | Default | Description |
|---|---|---|
| `--spin-texture` | — | Component to visualize: `sx`, `sy`, or `sz` |
| `--spin-cmap` | `seismic` | Diverging colormap (any matplotlib name) |
| `--lw` | `2.0` | Band line width |

!!! note "Requires `LSORBIT = .TRUE.` and `LORBIT ≥ 11` in INCAR."

**Examples:**

```bash
# Standard plot
valyte band --ylim -3 3 -o my_bands.png

# Tricolor orbital-resolved
valyte band --tricolor s p d --ylim -4 4

# Element-resolved
valyte band --tricolor Mo S Se --tricolors "#e74c3c" "#2ecc71" "#3498db" --ylim -3 3

# Spin-resolved
valyte band --spin-resolved --ylim -3 3

# Spin texture
valyte band --spin-texture sz --spin-cmap RdBu_r --ylim -2 2

# Vector output for publication
valyte band --format pdf

# Export data
valyte band --save-data
```

---

## `valyte dos`

Plot total and projected density of states from `vasprun.xml`.

```bash
valyte dos [path/to/vasprun.xml] [options]
```

| Option | Default | Description |
|---|---|---|
| `-e`, `--elements` | all | Elements/orbitals to plot |
| `--xlim` | `-6 6` | Energy range in eV |
| `--ylim` | auto | DOS range |
| `--scale` | `1.0` | Divide DOS values by this factor |
| `--fermi` | off | Draw dashed line at E = 0 |
| `--pdos` | off | Show only projected DOS (hide total) |
| `--legend-cutoff` | `0.10` | Hide legend if max PDOS fraction < threshold |
| `-o`, `--output` | `valyte_dos.png` | Output filename |
| `--font` | `Arial` | Font family |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--save-data` | off | Save DOS data to `valyte_dos.dat` |

Element/orbital formats: `Fe` (element total), `Fe(d)` or `Fe:d` (specific orbital), `-e Fe O` (multiple elements).

**Examples:**

```bash
valyte dos
valyte dos -e Fe O
valyte dos -e "Fe(d)" "O(p)"
valyte dos ./run --xlim -5 5 --fermi -o my_dos.png
valyte dos -e Fe O --format svg
valyte dos -e Fe O --save-data
```

---

## `valyte ipr`

Compute the Inverse Participation Ratio from a `PROCAR` file.

```bash
valyte ipr [options]
```

Interactive mode — reads `PROCAR`, shows system info, prompts for band indices, and saves results to `ipr_procar.dat`.

Use `--bands` for non-interactive/batch runs.

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `PROCAR` | Input PROCAR file |
| `-b`, `--bands` | interactive prompt | Band indices/ranges, e.g. `5`, `5-8`, or `5 8-10` |
| `-o`, `--output` | `ipr_procar.dat` | Output data filename |
| `--details` | off | Print per-k-point IPR values |

Band index input formats:

```
5           → single band
5 6 7       → multiple bands
5-7         → range (inclusive)
5 8-10 13   → mixed
```

Output columns: `Band`, `Energy (eV)`, `IPR`, `N_eff (= 1/IPR)`.

!!! note "Requires `LORBIT = 11` (or ≥ 10) to write atom-projected orbital weights into PROCAR."

**Examples:**

```bash
valyte ipr
valyte ipr --bands 5-8
valyte ipr -i PROCAR -b 5 8-10 13 -o defect_ipr.dat
valyte ipr --bands 6 --details
```

---

## `valyte effmass`

Compute carrier effective masses at VBM and CBM by parabolic fitting.

```bash
valyte effmass [options]
```

| Option | Default | Description |
|---|---|---|
| `--vasprun` | `.` | Path to `vasprun.xml` or directory |
| `--kpoints` | auto-detected | Path to `KPOINTS` file |
| `--npoints` | `3` | K-points on each side of extremum for fitting |
| `--band-index` | auto | Manual 1-indexed band indices to fit |
| `--tol` | `1e-3` | Degeneracy tolerance in eV |
| `--plot` | off | Save parabolic fit plot |
| `-o`, `--output` | `valyte_effmass.png` | Output plot filename (with `--plot`) |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--save-data` | off | Save results to `valyte_effmass.dat` |

!!! warning "Requires a line-mode band structure calculation. SCF calculations will raise an error."

**Examples:**

```bash
valyte effmass
valyte effmass --npoints 5
valyte effmass --plot -o effmass.png
valyte effmass --plot --format pdf
valyte effmass --save-data
```

---

## `valyte converge`

Monitor and visualize VASP convergence from `OSZICAR` and optionally `OUTCAR`.

```bash
valyte converge [path] [options]
```

| Option | Default | Description |
|---|---|---|
| `path` | `.` | Directory or path to `OSZICAR` |
| `--electronic` | off | Show SCF convergence instead of ionic |
| `--forces` | off | Add max-force panel (requires `OUTCAR`) |
| `--stress` | off | Add pressure panel (requires `OUTCAR`) |
| `--mag` | off | Add magnetization subplot |
| `--ethresh` | from `OUTCAR` / `1e-4` | Energy convergence reference line (eV) |
| `--fthresh` | `0.02` | Force convergence reference line (eV/Å) |
| `--start` | `1` | First ionic step to include |
| `--end` | last | Last ionic step to include |
| `--no-plot` | off | Terminal summary only — no plot |
| `-o`, `--output` | `valyte_converge.png` | Output plot filename |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--save-data` | off | Save data to `valyte_converge.dat` |

**Examples:**

```bash
valyte converge
valyte converge --forces --stress
valyte converge --electronic
valyte converge --no-plot
valyte converge /path/to/run --start 5 --end 30
valyte converge --forces --format pdf
valyte converge --forces -o converge.png --save-data
```
