# Combined Band Structure & DOS

Plot the electronic band structure and the projected density of states (PDOS) side-by-side with aligned energy axes ($E - E_F$).

```bash
valyte combined [path] [options]
```

The path can be a positional argument, passed via `--vasprun`, or omitted to use the current directory.

!!! important "No Total DOS"
    In accordance with publication best practices, the total density of states is never plotted in the combined side-by-side view. Only the element- or orbital-projected contributions (PDOS) are shown on the right-hand panel.

---

## Options

| Option | Default | Description |
|---|---|---|
| `--dos` | same as band | Path to `vasprun.xml` or directory for DOS data (if different from band structure) |
| `-e`, `--elements` | all | Elements/orbitals to plot in the DOS panel |
| `--ylim` | `-4 4` | Shared energy range in eV (shared Y-axis limits) |
| `--scale` | `1.0` | Divide all DOS values by this factor (zoom in on small features) |
| `--fermi` | off | Draw a dashed horizontal line at the Fermi level (E = 0) |
| `--legend-cutoff` | `0.10` | Hide the legend if the max PDOS contribution is below this fraction of total |
| `-o`, `--output` | `valyte_combined.png` | Output filename |
| `--font` | `Arial` | Font family: `Arial`, `Helvetica`, `Times New Roman` |
| `--width` | `3.2` | Plot width in inches |
| `--height` | `3.2` | Plot height in inches |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--save-data` | off | Save both band and DOS data to text files |
| `--spin-resolved` | off | Plot spin-up/spin-down channels in distinct colors |
| `--no-bold` | off | Use normal font weight and thinner lines/ticks (scientific style) |

---

## Element and orbital formats

The `-e` / `--elements` flag accepts a mix of formats in a single command:

| Format | Example | Plots |
|---|---|---|
| Element | `Fe` | Total PDOS for Fe (all orbitals summed) |
| Orbital | `Fe(d)` or `Fe:d` | Only the d-orbital contribution of Fe |
| Multiple | `-e Fe O` | Total PDOS for Fe and O |
| Mixed | `-e Fe "Fe(d)"` | Fe total and Fe d-orbital on the same plot |

---

## Examples

```bash
# Side-by-side plot of band structure and element-projected DOS (default)
valyte combined

# Custom energy range, Fermi line, and output filename
valyte combined --ylim -3 3 --fermi -o band_dos_combined.png

# Show specific elements/orbitals in the DOS panel
valyte combined -e "Bi(p)" "O(p)"

# Plot spin-up/spin-down channels in distinct colors
valyte combined --spin-resolved

# Save raw band structure and DOS data alongside the plot
valyte combined --save-data
```

---

## Exported data format

When using `--save-data`, two separate files are written:

1. **`valyte_combined_band.dat`**:
   - Contains the band structure energies relative to Fermi level (or VBM) along the cumulative $k$-path.
   - Headers contain $k$-point labels and distances.

2. **`valyte_combined_dos.dat`**:
   - Contains the energy coordinates and the selected projected density of states (PDOS) columns for the right panel.
