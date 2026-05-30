# Effective Mass

Compute carrier effective masses (m\*/m₀) at the VBM and CBM by fitting a parabola to the band dispersion along each high-symmetry direction. Works from `vasprun.xml` alone.

```bash
valyte effmass [options]
```

---

## How it works

For each band extremum (VBM for holes, CBM for electrons), Valyte fits:

$$E(k) = \frac{\hbar^2 k^2}{2 m^*} + \text{const}$$

using a least-squares parabolic fit to a small window of k-points around the extremum. The effective mass is extracted from the quadratic coefficient:

$$\frac{m^*}{m_0} = \frac{\hbar^2}{2a \cdot m_0}$$

where $a$ is the fitted coefficient in eV·Å².

A fit is produced for each high-symmetry direction in the band structure — so if the VBM sits at Γ and the path passes through Γ twice (e.g. Γ→X and Γ→Y), you get separate masses for each direction.

!!! warning "Non-parabolic bands"
    If the band is not well described by a parabola in the fitting window, a warning is printed. Narrow-gap semiconductors, d-bands, and bands with strong spin-orbit coupling are common cases where the parabolic approximation breaks down.

---

## Options

| Option | Default | Description |
|---|---|---|
| `--vasprun` | `.` | Path to `vasprun.xml` or directory containing it |
| `--kpoints` | auto-detected | Path to `KPOINTS` file for high-symmetry labels |
| `--npoints` | `3` | K-points on each side of extremum used for fitting |
| `--band-index` | auto | Manual 1-indexed band indices to fit |
| `--tol` | `1e-3` | Degeneracy tolerance in eV (for degenerate bands at extremum) |
| `--plot` | off | Save parabolic fit plot |
| `-o`, `--output` | `valyte_effmass.png` | Output plot filename (used with `--plot`) |
| `--save-data` | off | Save results to `valyte_effmass.dat` |

---

## Usage

```bash
# Auto-detect VBM and CBM, print masses
valyte effmass

# Wider fitting window (more k-points → smoother fit)
valyte effmass --npoints 5

# Save parabolic fit figure
valyte effmass --plot -o effmass.png

# Export numerical results
valyte effmass --save-data

# Fit specific bands manually (1-indexed)
valyte effmass --band-index 42 43
```

---

## Output

### Terminal

```
Effective mass
══════════════

  Band gap:     1.124 eV (indirect)
  VBM:          k = (0.000, 0.000, 0.000)  [Γ]    Band 42
  CBM:          k = (0.500, 0.000, 0.500)  [X]    Band 43

  Hole effective masses (VBM):
    Γ → X    :     m* =    0.284 m₀
    Γ → L    :     m* =    0.192 m₀
    Γ → K    :     m* =    0.231 m₀

  Electron effective masses (CBM):
    X → Γ    :     m* =    0.916 m₀
    X → U    :     m* =    0.191 m₀
```

### Data file (`--save-data`)

`valyte_effmass.dat` — plain text with columns: Type, Band, Spin, Direction, m\*/m₀.

---

## Requirements

- Line-mode band structure calculation (`IBRION = -1`, line-mode `KPOINTS`)
- `vasprun.xml` written with `LWAVE = .FALSE.` or `LCHARG = .FALSE.` is fine — projections are not required
- Self-consistent (SCF) calculations will raise an error — a k-path is needed

!!! tip
    Run `valyte band kpt-gen` to generate the line-mode KPOINTS, do a non-self-consistent band calculation, then point `valyte effmass` at the resulting `vasprun.xml`.
