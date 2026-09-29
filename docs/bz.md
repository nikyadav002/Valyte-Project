# Brillouin Zone

Draw the first Brillouin zone of a structure with the suggested high-symmetry
k-path on it. Useful for checking the path before you commit to a band
structure run, since it shows exactly which points `valyte band kpt-gen` would
write and where they sit in reciprocal space.

```bash
valyte bz [options]
```

---

## Options

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input structure |
| `-o`, `--output` | `valyte_bz.png` | Output plot filename |
| `--mode` | `bradcrack` | K-path convention: `bradcrack`, `seekpath`, `hinuma`, `setyawan_curtarolo`, `latimer_munro` |
| `--symprec` | `0.01` | Symmetry precision |
| `--elev` | `22` | Viewing elevation in degrees |
| `--azim` | `30` | Viewing azimuth in degrees |
| `--width` | `5.5` | Plot width in inches |
| `--height` | `5.5` | Plot height in inches |
| `--format` | from `-o` extension | Output figure format: `png`, `pdf`, or `svg` |
| `--dpi` | `400` | Output figure resolution |
| `--fontsize` | per-plot default | Base font size in points; all labels scale with it |
| `--no-bold` | off | Normal font weight and thinner lines |

---

## What's drawn

- The first Brillouin zone, the Wigner-Seitz cell of the reciprocal lattice
- The k-path, as connected red segments. A break between path branches is a
  discontinuity and is left unconnected
- Each high-symmetry point, labelled (Γ, X, L, W and so on)
- The reciprocal basis vectors b₁, b₂, b₃ from the origin, scaled down to keep
  them clear of the zone

The zone is drawn for the cell the path is defined in, which is the same
standardized primitive cell `valyte band kpt-gen` writes as `POSCAR_standard`.
That cell depends on `--mode`, so the zone shape can differ slightly between
conventions for the same input.

---

## Usage examples

```bash
# Default path and view
valyte bz

# A different convention
valyte bz --mode seekpath

# Look down on the zone
valyte bz --elev 60 --azim 120

# Larger type for a slide
valyte bz --fontsize 18 --width 7 --height 7

# Vector output
valyte bz --format pdf
```

---

## Viewing angle

The figure is static, so `--elev` and `--azim` are the only way to change the
viewpoint. If a label sits behind an edge or a path segment is hidden, rotate
the view rather than reaching for another tool.

---

## Points outside the zone

For low-symmetry lattices some conventions place high-symmetry points outside
the first Brillouin zone, on the reciprocal unit cell rather than the
Wigner-Seitz cell. Triclinic cells do this in every convention, and hexagonal
cells do it for `bradcrack`. Valyte draws those points where the convention puts
them rather than folding them back, so the figure matches the KPOINTS file. The
axes are sized to include them, so nothing is clipped.
