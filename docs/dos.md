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
| `-c`, `--colors` | none | Custom colours: per series (`Fe(d)=red O(p)=0072B2`) and/or a plain list that replaces the palette (`red blue`) |
| `--palette` | `wong` | Colour palette for the projected DOS: `wong`, `bright` or `muted` (see [Colours](#colours)) |
| `--xlim` | `-6 6` | Energy range in eV |
| `--ylim` | auto | DOS range |
| `--scale` | `1.0` | Divide all DOS values by this factor (zoom in on small features) |
| `--fermi` | off | Draw a dashed vertical line at the Fermi level (E = 0) |
| `--pdos` | off | Show only projected DOS (hide total DOS) |
| `--legend-cutoff` | `0.10` | Channels whose peak is below this fraction of the largest are drawn grey and left out of the legend. Channels named with `-e` are always shown |
| `-o`, `--output` | `valyte_dos.png` | Output filename |
| `--font` | `Arial` | Font family: `Arial`, `Helvetica`, `Times New Roman` |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--dpi` | `400` | Output figure resolution in DPI |
| `--save-data` | off | Save DOS data to a `.dat` named after `-o` (default `valyte_dos.dat`, or `valyte_dos_panels.dat` with `--panels`) |
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
| Hex Code | `Fe(d):0072B2` or `Fe:d=0072B2` | Hex code without `#` (convenient in CLI) |

---

## Colours

### Palettes

Projected DOS curves take their colours from a palette, chosen with `--palette`:

| Palette | Colours | Character |
|---|---|---|
| `wong` (default) | <span style="color:#0072B2">■</span> <span style="color:#D55E00">■</span> <span style="color:#009E73">■</span> <span style="color:#CC79A7">■</span> <span style="color:#E69F00">■</span> <span style="color:#56B4E9">■</span> | Bang Wong's colour-blind-safe set: blue, vermillion, bluish green, reddish purple, orange, sky blue |
| `bright` | <span style="color:#4477AA">■</span> <span style="color:#EE6677">■</span> <span style="color:#228833">■</span> <span style="color:#CCBB44">■</span> <span style="color:#66CCEE">■</span> <span style="color:#AA3377">■</span> | Paul Tol's bright scheme, designed to stay distinct with colour blindness |
| `muted` | <span style="color:#CC6677">■</span> <span style="color:#332288">■</span> <span style="color:#DDCC77">■</span> <span style="color:#117733">■</span> <span style="color:#88CCEE">■</span> <span style="color:#882255">■</span> <span style="color:#44AA99">■</span> <span style="color:#999933">■</span> <span style="color:#AA4499">■</span> | Paul Tol's muted scheme: softer tones, and the most colours (9) |

```bash
valyte dos                    # wong
valyte dos --palette bright
valyte dos --palette muted
```

The total DOS is always black, so no palette uses black.

### How colours are assigned

- **Biggest channels first.** Channels that make the legend take palette colours in order of how many states they have in the energy window, and the legend lists them in that order. The most important channels get the strongest colours.
- **Small channels are grey.** Channels whose peak is below `--legend-cutoff` (10% of the largest by default) are drawn in light grey and left out of the legend, so every coloured curve on the plot has a label.
- **`-e` is explicit.** Channels you name with `-e` are always coloured and listed, in the order you typed them, however small they are.
- **Never two identical entries.** With more channels than colours, colours repeat with dashed, then dotted, lines. Valyte prints a note when that happens; `-e` or `--panels` usually make a clearer figure.
- **Panels keep colours per entity.** In `--panels`, each orbital keeps one colour in every element panel (s, p, d, f in palette order); with `--panel-by orbital`, each element does.

### Custom colours

Custom colours always override the palette. A key can be a series (`Fe(d)`), an element (`Fe`, all its shown channels), or an orbital (`d`). Keys are not case-sensitive for elements (`fe(D)` works), but a single lowercase `s`, `p`, `d` or `f` means the orbital, while `S` and `P` mean sulfur and phosphorus.

**Inline in `-e`:**
```bash
valyte dos -e "Fe(d):red" "O(p):0072B2"
```

**With `-c` / `--colors`.** Mappings set single series; a plain list of colours replaces the palette. Both can be used together:
```bash
valyte dos -c "Fe(d)=red" "O(p)=0072B2"     # per series
valyte dos -c red blue green               # replaces the palette
valyte dos -c black orange "Fe(d)=red"      # both
```

**From a file (`--colors-file`)**, as a JSON object or as `key colour` lines (`#` starts a comment):

**`colors.json`**:
```json
{
  "Fe(d)": "#D55E00",
  "Fe(s)": "orange",
  "O(p)": "0072B2"
}
```

**`colors.txt`**:
```
# series   colour
Fe(d)      #D55E00
O:p        0072B2
Pb = black
```

```bash
valyte dos -e "Fe(d)" "O(p)" --colors-file colors.json
```

When the same series is set in several places, the later one wins: `--colors-file`, then `-c`, then inline `-e` colours. Palette colours already used by a custom mapping are skipped, so an automatic colour never duplicates one you picked.

Colours can be any Matplotlib colour: names (`red`), hex with or without `#` (`#0072B2`, `0072B2`), or prefixed names (`tab:blue`, `xkcd:sky blue`). An invalid colour, an unknown palette or a missing colours file stops with a short error message. Keys that match no plotted series are reported with the list of valid ones.

---

## Examples

```bash
# All orbitals for all elements (default)
valyte dos

# Total PDOS for specific elements
valyte dos -e Fe O

# Specific orbitals with custom colors file
valyte dos -e "Fe(d)" "O(p)" --colors-file colors.json

# Inline custom colours with hex codes
valyte dos -e Fe(d):D55E00 O(p):0072B2

# Another palette
valyte dos --palette muted

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
- Spin-down values are written as positive numbers; only the plot draws them below zero

---

## Smart features

**Adaptive legend**
:   Only channels that matter get a colour and a legend entry; the rest are drawn in grey, and the legend lists channels from largest to smallest. See [Colours](#colours) for the exact rules and [FAQ → A channel is grey or missing from the legend](faq.md#a-channel-is-grey-or-missing-from-the-legend).

**Gradient fill**
:   Peaks are filled with a smooth gradient tied to the line color, giving plots a clean publication aesthetic.

**Orbital-resolved by default**
:   When no `-e` flag is given, Valyte plots all available elements broken down by s, p, d (and f if present) orbitals automatically.

**Panel mode**
:   The `--panels` flag tiles each element (or orbital) into its own subplot with independent y-scaling, making it easy to compare contributions that differ by orders of magnitude.
