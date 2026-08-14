# Density of States

Plot total and projected density of states from `vasprun.xml` with orbital resolution, adaptive legends, and gradient fills.

```bash
valyte dos [path/to/vasprun.xml] [options]
```

The path can be a positional argument, passed via `--vasprun`, or omitted to use the current directory.

!!! tip "Combine with band structure"
    For a complete electronic structure analysis, pair your DOS plot with a [band structure](band.md). Both read from the same `vasprun.xml`.

---

## Options

| Option | Default | Description |
|---|---|---|
| `-e`, `--elements` | all | Elements/orbitals to plot, optionally with inline colors (`Fe(d):red`, `O(p):#2a9d8f`) |
| `--colors-file` | none | Path to JSON or text file mapping elements/orbitals/species to colors |
| `-c`, `--colors` | none | Custom colors mapping or list (`Fe(d)=red O(p)=blue` or `red blue`) |
| `--xlim` | `-6 6` | Energy range in eV |
| `--ylim` | auto | DOS range |
| `--scale` | `1.0` | Divide all DOS values by this factor (zoom in on small features) |
| `--fermi` | off | Draw a dashed vertical line at the Fermi level (E = 0) |
| `--pdos` | off | Show only projected DOS (hide total DOS) |
| `--legend-cutoff` | `0.10` | Hide the legend if the max PDOS contribution is below this fraction of total |
| `-o`, `--output` | `valyte_dos.png` | Output filename |
| `--font` | `Arial` | Font family: `Arial`, `Helvetica`, `Times New Roman` |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--save-data` | off | Save DOS data to `valyte_dos.dat` |
| `--panels` | off | Split DOS into vertically stacked panels (one per element) |
| `--panel-by` | `element` | Grouping mode for panels: `element` or `orbital` |
| `--no-bold` | off | Use normal font weight and thinner lines/ticks (scientific style) |

---

## Element, orbital, and color formats

The `-e` / `--elements` flag accepts element names, orbital projections, and inline color assignments (using names or hex codes):

| Format | Example | Description |
|---|---|---|
| Element | `Fe` | Total PDOS for Fe (all orbitals summed) |
| Orbital | `Fe(d)` or `Fe:d` | Only the d-orbital contribution of Fe |
| Inline Color | `Fe(d):red` or `O(p):#2a9d8f` | Assign specific colors (names or hex codes) |
| Hex Code | `Fe(d):e63946` | Hex code without `#` (convenient in CLI) |

---

## Color Customization

You can customize colors via a config file, CLI `--colors`, or inline in `-e`:

### Using a colors file (`--colors-file`)
Create a JSON or simple key-value file mapping elements, orbital channels, or specific species (`Fe(d)`) to colors:

**`colors.json`**:
```json
{
  "Fe(d)": "#e63946",
  "Fe(s)": "orange",
  "O(p)": "#00b4d8"
}
```

```bash
valyte dos -e "Fe(d)" "O(p)" --colors-file colors.json
```

### Using CLI mapping (`--colors`)
```bash
# Mapping mode
valyte dos -e "Fe(d)" "O(p)" --colors Fe(d)=red O(p)=blue

# Ordered list mode
valyte dos -e "Fe(d)" "O(p)" --colors red blue
```

---

## Examples

```bash
# All orbitals for all elements (default)
valyte dos

# Total PDOS for specific elements
valyte dos -e Fe O

# Specific orbitals with custom colors file
valyte dos -e "Fe(d)" "O(p)" --colors-file colors.json

# Inline custom colors with hex codes
valyte dos -e Fe(d):e63946 O(p):2a9d8f

# Custom path, energy range, Fermi line, output
valyte dos ./vasp_run --xlim -5 5 --fermi -o my_dos.png


# Show only projected DOS, no total DOS
valyte dos --pdos -e Fe O

# Scale down a very tall DOS
valyte dos --scale 3 --ylim 0 10

# Vector output for publication
valyte dos -e Fe O --format pdf

# Save raw data alongside the plot
valyte dos -e Fe O --save-data
```

---

## Panels

Use `--panels` to split the DOS into a vertically stacked figure with one panel per element. Each panel shows the orbital breakdown for that element, with a semi-transparent total DOS backdrop for context.

```bash
# One panel per element (default grouping)
valyte dos --panels

# Group by orbital instead — one panel per s, p, d channel
valyte dos --panels --panel-by orbital

# Panels with Fermi line and specific elements
valyte dos --panels --fermi -e Fe O

# Panels without total DOS backdrop
valyte dos --panels --pdos

# Custom energy range and PDF output
valyte dos --panels --xlim -4 4 --format pdf
```

!!! tip "When to use panels"
    Panels are most useful for multi-element systems (e.g., ternary or quaternary compounds) where a single-axis DOS plot becomes crowded. Each element gets its own y-axis scaling, making small contributions visible.

---

## Exported data format (`valyte_dos.dat`)

A plain-text, whitespace-delimited file with one row per energy point:

```
# Energy(eV)  Total_up  Fe_up  O_up  ...
-6.000000  0.000000  0.000000  0.000000  ...
...
```

- **`Energy(eV)`** — energy relative to the VBM (E = 0)
- **`Total_up`** — total DOS spin-up channel
- **`Total_dn`** — total DOS spin-down channel (spin-polarized calculations only)
- Remaining columns are one per plotted element/orbital (e.g. `Fe_up`, `Fe(d)_up`) matching exactly what is shown in the plot
- Spin-down columns are labeled `_dn` and appear immediately after their `_up` counterpart

---

## Smart features

**Adaptive legend**
:   If the maximum PDOS contribution across all plotted species is below the `--legend-cutoff` threshold (default 10% of the total DOS peak), the legend is automatically hidden to avoid cluttering a flat plot. See [FAQ → Legend is not showing](faq.md#legend-is-not-showing) if you need to adjust this.

**Gradient fill**
:   Peaks are filled with a smooth gradient tied to the line color, giving plots a clean publication aesthetic.

**Orbital-resolved by default**
:   When no `-e` flag is given, Valyte plots all available elements broken down by s, p, d (and f if present) orbitals automatically.

**Panel mode**
:   The `--panels` flag tiles each element (or orbital) into its own subplot with independent y-scaling, making it easy to compare contributions that differ by orders of magnitude.
