# Convergence Monitor

Monitor and visualize the convergence of VASP relaxation, single-point, or MD calculations. Parses `OSZICAR` for the basic ionic/electronic summary and optionally `OUTCAR` for forces and pressure.

```bash
valyte converge [path] [options]
```

The default mode (no extra flags) needs only `OSZICAR` and produces a 2-panel plot — useful over SSH where OUTCAR may be large or incomplete.

---

## Options

| Option | Default | Description |
|---|---|---|
| `path` | `.` | Directory containing `OSZICAR`/`OUTCAR`, or direct path to `OSZICAR` |
| `--electronic` | off | Show SCF convergence instead of ionic |
| `--forces` | off | Add max-force panel (requires `OUTCAR`) |
| `--stress` | off | Add external pressure panel (requires `OUTCAR`) |
| `--mag` | off | Add magnetization subplot (spin-polarized runs) |
| `--ethresh` | from `OUTCAR` / `1e-4` | Energy convergence reference line (eV) |
| `--fthresh` | `0.02` | Force convergence reference line (eV/Å) |
| `--start` | `1` | First ionic step to include in the plot |
| `--end` | last | Last ionic step to include |
| `--no-plot` | off | Print terminal summary only — no plot generated |
| `-o`, `--output` | `valyte_converge.png` | Output plot filename |
| `--save-data` | off | Save parsed data to `valyte_converge.dat` |

---

## Modes

### Ionic convergence (default)

Produces a vertically stacked, shared-x-axis multi-panel plot:

- **Panel 1** — Total energy (E₀, σ→0) vs ionic step
- **Panel 2** — |ΔE| vs ionic step (log scale) with reference line at `--ethresh`
- **Panel 3** — Max atomic force vs ionic step (log scale), only with `--forces`
- **Panel 4** — External pressure vs ionic step, only with `--stress`
- **Panel 5** — Magnetization, only with `--mag` (spin-polarized runs)

```bash
# Basic — OSZICAR only, 2 panels
valyte converge

# With forces
valyte converge --forces

# Full output
valyte converge --forces --stress --mag
```

### Electronic (SCF) convergence

`--electronic` produces a single-panel plot showing |ΔE| per SCF iteration on a continuous x-axis across all ionic steps. Ionic step boundaries are marked with vertical dashed lines.

```bash
valyte converge --electronic
```

This is useful for diagnosing slow or oscillating SCF convergence — for example, identifying ionic steps where the SCF needed many iterations.

---

## Usage examples

```bash
# Run from the calculation directory
valyte converge

# Point at a specific directory
valyte converge /path/to/run

# Quick SSH check — no plot, just the summary
valyte converge --no-plot

# Include forces (reads OUTCAR)
valyte converge --forces --fthresh 0.01

# Zoom into steps 10–50
valyte converge --start 10 --end 50

# SCF convergence view
valyte converge --electronic

# Save the plot and data
valyte converge --forces -o converge.png --save-data
```

---

## Terminal output

A convergence summary is always printed before the plot:

```
Convergence
═══════════

  Calculation type:   Relaxation (IBRION = 2, NSW = 60)
  Ionic steps:        23 / 60
  Status:             ✓ Converged

  Energy:
    Final E0        = -148.29374182 eV
    ΔE (last step)  =   2.41e-06 eV

  Forces:
    Max |F| (final) =   0.0087 eV/Å
    Threshold       =   0.0200 eV/Å    ✓

  Timing:
    Total walltime  =  4h 12m 33s
    Avg per step    =  11m 00s
```

Calculation type, total steps, and convergence thresholds are read automatically from `OUTCAR` or `INCAR` if present.

---

## Data file (`--save-data`)

`valyte_converge.dat` — plain text with columns: Step, E0 (eV), ΔE (eV), F_max (eV/Å), P (kB), mag.

---

## File handling

- Reads `OSZICAR` or `OSZICAR.gz` transparently
- Reads `OUTCAR` or `OUTCAR.gz` transparently (only when `--forces` or `--stress` is set)
- Falls back to `INCAR` for convergence parameters if `OUTCAR` is absent
- Handles incomplete files — works on running calculations
- Handles restarted calculations where the ionic step counter resets
