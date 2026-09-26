# Convergence Monitor

Monitor the convergence of a VASP **structural relaxation**. Parses `OSZICAR`
for the per-ionic-step energies and `OUTCAR` for the max force and pressure,
then prints a table to the terminal. No plot is produced.

```bash
valyte converge [path]
```

---

## Options

| Option | Default | Description |
|---|---|---|
| `path` | `.` | Directory containing `OSZICAR`/`OUTCAR`, or direct path to `OSZICAR` |
| `--save-data` | off | Save parsed data to `valyte_converge.dat` |

The convergence criterion is not a flag — it is read from `EDIFFG`, so the
summary always reflects the criterion the calculation actually used.

---

## Convergence criterion

`EDIFFG` is read from `INCAR`, falling back to the `INCAR` echo in `OUTCAR`
when no `INCAR` is present. VASP's sign convention decides which column is
checked:

| `EDIFFG` | Criterion | Checked against |
|---|---|---|
| negative (e.g. `-0.02`) | Force, `|EDIFFG|` eV/Å | `Max \|F\|` column |
| positive (e.g. `1E-4`) | Energy, `EDIFFG` eV | `ΔE` column |
| absent | none | `Status` falls back to step count |

A ✓ next to a value marks a step meeting the criterion.

---

## Terminal output

```
Convergence
═══════════

  Calculation type:   Relaxation (IBRION = 2, NSW = 60)
  Ionic steps:        4 / 60
  Status:             ✓ Converged

   Step            E0 (eV)          ΔE (eV)     Max |F| (eV/Å)       P (kB)
  ─────────────────────────────────────────────────────────────────────────
      1      -303.72457000        -3.04e+02             0.5412       -12.34
      2      -303.81234000        -8.78e-02             0.0873        -3.21
      3      -303.81501000        -2.67e-03           0.0154 ✓        -0.87
      4      -303.81509000        -8.10e-05           0.0087 ✓        -0.12

  EDIFFG           =  -0.02  →  force criterion
  Force threshold  =  0.0200 eV/Å
  Final Max |F|    =  0.0087 eV/Å    ✓
```

With a positive `EDIFFG` the ✓ marks move to the `ΔE` column and the footer
reports the energy criterion instead:

```
  EDIFFG           =  0.0001  →  energy criterion
  Energy threshold =  1.00e-04 eV
  Final ΔE         =  -8.10e-05 eV       ✓
```

`Step`, `E0`, and `ΔE` come from `OSZICAR`; `Max |F|` and `P` need `OUTCAR` and
show `—` when it is absent or does not carry the value. Calculation type and
total step count are read from `OUTCAR`, falling back to `INCAR`.

---

## Usage examples

```bash
# Run from the calculation directory
valyte converge

# Point at a specific directory
valyte converge /path/to/run

# Point directly at an OSZICAR
valyte converge /path/to/run/OSZICAR

# Save the parsed data alongside the summary
valyte converge --save-data
```

---

## Data file (`--save-data`)

`valyte_converge.dat` — plain text with columns: Step, E0 (eV), ΔE (eV),
F_max (eV/Å), P (kB), mag. Columns sourced from `OUTCAR` are written as `—`
when it is not available.

---

## File handling

Valyte handles common edge cases automatically:

- Reads `OSZICAR` or `OSZICAR.gz` and `OUTCAR` or `OUTCAR.gz` transparently
- Works without `OUTCAR` — the force and pressure columns show `—`
- Falls back to `INCAR` for calculation parameters if `OUTCAR` is absent
- Handles incomplete files — works on running calculations. A trailing SCF cycle
  with no ionic summary yet is not counted as a step.
- Handles restarted calculations where the ionic step counter resets

See [FAQ → OSZICAR not found](faq.md#filenotfounderror-oszicar-not-found) or
[FAQ → Force column shows no data](faq.md#force-column-shows-no-data) for
troubleshooting.
