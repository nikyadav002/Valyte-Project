# FAQ & Troubleshooting

Common issues and solutions when using Valyte.

---

## Installation

### `valyte: command not found`

**Cause:** The `valyte` script is not on your `PATH`. This usually happens when pip installs to a user-local directory.

**Fix:**

```bash
# Check if valyte is installed
pip show valyte

# Option 1: Install with --user and add to PATH
pip install --user valyte
export PATH="$HOME/.local/bin:$PATH"  # Add to your .bashrc or .zshrc

# Option 2: Install in a virtual environment (recommended)
python -m venv venv
source venv/bin/activate
pip install valyte
```

### `ModuleNotFoundError: No module named 'valyte'`

**Cause:** You may have multiple Python installations, and `valyte` is installed in a different one.

**Fix:** Ensure you're using the same Python that pip installed to:

```bash
python -m valyte.cli --help
# or
python -m pip install valyte
```

---

## Band structure

### My k-point labels are wrong or missing

**Cause:** You ran the band calculation with the original `POSCAR` instead of the standardized cell.

**Fix:** After generating k-points, always use the standardized cell:

```bash
valyte band kpt-gen
cp POSCAR_standard POSCAR    # Use this for the band calculation
```

The k-path is defined for the standardized primitive cell. Using a different cell orientation will produce incorrect or mismatched labels.

### My tricolor plot is all one color

**Cause:** VASP did not write orbital projections into `vasprun.xml`.

**Fix:** Add `LORBIT = 11` to your `INCAR` and re-run the calculation. This writes the projected eigenvalues that tricolor mode needs.

```
# In INCAR
LORBIT = 11
```

### `ValueError: No spin channels found for spin-resolved mode`

**Cause:** You used `--spin-resolved` on a non-spin-polarized calculation.

**Fix:** Either:
- Run VASP with `ISPIN = 2` in your `INCAR` for a spin-polarized calculation, or
- Use `valyte band` without `--spin-resolved` for non-magnetic systems

Valyte will automatically fall back to a standard plot with a warning.

### `ValueError: No projected_magnetisation data found`

**Cause:** You used `--spin-texture` on a collinear or non-SOC calculation.

**Fix:** Spin texture requires a non-collinear calculation:

```
# In INCAR
LSORBIT = .TRUE.
LORBIT = 11
```

---

## Density of states

### The DOS looks noisy or has sharp spikes

**Cause:** Insufficient k-point sampling in the SCF calculation.

**Fix:** Increase the k-point density. For metals, a very dense grid (e.g., k-spacing < 0.03 in 2π/Å) is often needed for smooth DOS.

### Legend is not showing

**Cause:** The `--legend-cutoff` threshold is hiding the legend because the projected contributions are small relative to the total DOS.

**Fix:** Lower the cutoff:

```bash
valyte dos -e Fe O --legend-cutoff 0.01
```

---

## Effective mass

### `Error: This does not appear to be a line-mode band calculation`

**Cause:** You're pointing `valyte effmass` at a self-consistent (SCF) calculation rather than a non-self-consistent band structure calculation.

**Fix:** Run a band structure calculation first:

1. Generate k-points: `valyte band kpt-gen`
2. Copy `POSCAR_standard` to `POSCAR`
3. Run VASP with `ICHARG = 11` and the generated line-mode `KPOINTS`
4. Then run `valyte effmass`

### Effective masses are unreasonably large or negative

**Cause:** The parabolic approximation breaks down for non-parabolic bands (narrow-gap semiconductors, d-bands, strong SOC).

**Fix:**

- Try reducing the fitting window: `valyte effmass --npoints 2`
- Use `--plot` to visually inspect the parabolic fit
- Consider whether a parabolic effective mass is physically meaningful for your system

---

## POTCAR

### `ValueError: No POTCAR for element ...` or pymatgen errors

**Cause:** Pymatgen cannot find your VASP pseudopotential files.

**Fix:** Set up `~/.pmgrc.yaml`:

```yaml
PMG_VASP_PSP_DIR: /path/to/your/pseudopotentials
```

The directory should contain subdirectories like `POT_GGA_PAW_PBE/`, `POT_GGA_PAW_PBE_52/`, etc. See the [pymatgen documentation](https://pymatgen.org/installation.html#potcar-setup).

---

## Convergence

### `FileNotFoundError: OSZICAR not found`

**Cause:** Valyte cannot find `OSZICAR` in the specified directory.

**Fix:**

```bash
# Point to the correct directory
valyte converge /path/to/calculation/directory

# Or run from the calculation directory
cd /path/to/calculation
valyte converge
```

### Force panel shows no data

**Cause:** `OUTCAR` is missing or incomplete.

**Fix:** The `--forces` and `--stress` flags require `OUTCAR`. Without it, only energy and ΔE panels (from `OSZICAR`) are available. Make sure the `OUTCAR` file exists in the same directory as `OSZICAR`.

---

## General

### Can I use Valyte on a remote cluster over SSH?

Yes. Valyte uses matplotlib's `Agg` backend by default, so it does not require a display. Plots are saved to files.

For a quick convergence check without generating a plot:

```bash
valyte converge --no-plot
```

### How do I export data for custom plotting?

Every plot command supports `--save-data`:

```bash
valyte dos --save-data           # → valyte_dos.dat
valyte band --save-data          # → valyte_band.dat
valyte effmass --save-data       # → valyte_effmass.dat
valyte converge --save-data      # → valyte_converge.dat
```

These are plain-text, whitespace-delimited files that can be loaded into any plotting tool (Python, MATLAB, gnuplot, Origin, etc.).

### Something else is broken

[Open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new) with:

- Valyte version (`pip show valyte`)
- Python version (`python --version`)
- The full error traceback
- The command you ran
- Relevant VASP input/output files if possible
