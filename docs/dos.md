# DOS — Density of States

Plot total and projected density of states from `vasprun.xml` with orbital resolution, adaptive legends, and gradient fills.

```bash
valyte dos [path/to/vasprun.xml] [options]
```

The path can be a positional argument, passed via `--vasprun`, or omitted to use the current directory.

---

## Options

| Option | Default | Description |
|---|---|---|
| `-e`, `--elements` | all | Elements/orbitals to plot (see formats below) |
| `--xlim` | `-6 6` | Energy range in eV |
| `--ylim` | auto | DOS range |
| `--scale` | `1.0` | Divide all DOS values by this factor (zoom in on small features) |
| `--fermi` | off | Draw a dashed vertical line at the Fermi level (E = 0) |
| `--pdos` | off | Show only projected DOS (hide total DOS) |
| `--legend-cutoff` | `0.10` | Hide the legend if the max PDOS contribution is below this fraction of total |
| `-o`, `--output` | `valyte_dos.png` | Output filename |
| `--font` | `Arial` | Font family: `Arial`, `Helvetica`, `Times New Roman` |

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
# All orbitals for all elements (default)
valyte dos

# Total PDOS for specific elements
valyte dos -e Fe O

# Specific orbitals only
valyte dos -e "Fe(d)" "O(p)"

# Fe total PDOS alongside Fe d-orbital
valyte dos -e Fe "Fe(d)"

# Custom path, energy range, Fermi line, output
valyte dos ./vasp_run --xlim -5 5 --fermi -o my_dos.png

# Show only projected DOS, no total DOS
valyte dos --pdos -e Fe O

# Scale down a very tall DOS
valyte dos --scale 3 --ylim 0 10
```

---

## Smart features

**Adaptive legend:** If the maximum PDOS contribution across all plotted species is below the `--legend-cutoff` threshold (default 10% of the total DOS peak), the legend is automatically hidden to avoid cluttering a flat plot.

**Gradient fill:** Peaks are filled with a smooth gradient tied to the line color, giving plots a clean publication aesthetic.

**Orbital-resolved by default:** When no `-e` flag is given, Valyte plots all available elements broken down by s, p, d (and f if present) orbitals automatically.
