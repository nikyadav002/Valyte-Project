# Pre-processing

Valyte provides three pre-processing utilities for setting up VASP calculations.

---

## Supercell

Generate a supercell from a POSCAR file.

```bash
valyte supercell nx ny nz [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `POSCAR_supercell` | Output filename |

```bash
# 2×2×2 supercell
valyte supercell 2 2 2

# 3×3×1 slab supercell with custom filenames
valyte supercell 3 3 1 -i POSCAR_primitive -o POSCAR_3x3x1
```

---

## K-Points — Interactive SCF grid

Generate a `KPOINTS` file for SCF or relaxation calculations interactively.

```bash
valyte kpt
```

The command will prompt for:

1. **K-mesh scheme** — Monkhorst-Pack or Gamma
2. **K-spacing** — in units of 2π/Å (e.g. `0.04`)

The optimal k-grid dimensions are automatically calculated from the lattice vectors in your `POSCAR`.

!!! tip
    A k-spacing of 0.03–0.05 (2π/Å) is typical for DFT calculations. Use smaller values (denser grids) for metals or accurate total energy convergence.

---

## POTCAR

Generate a `POTCAR` file from the species listed in your POSCAR.

```bash
valyte potcar [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `POTCAR` | Output filename |
| `--functional` | `PBE` | Pseudopotential functional |

**Available functionals:** `PBE`, `PBE_52`, `PBE_54`, `LDA`, `LDA_52`, `LDA_54`, `PW91`, and others supported by Pymatgen.

```bash
# Default PBE
valyte potcar

# Use PBE_54
valyte potcar --functional PBE_54

# Custom input and output
valyte potcar -i POSCAR_relaxed -o POTCAR_new
```

!!! important "Pymatgen configuration required"
    This command uses Pymatgen to locate your VASP pseudopotential files. You must set `PMG_VASP_PSP_DIR` in `~/.pmgrc.yaml` to point to the directory containing your PAW potentials.

    ```yaml
    PMG_VASP_PSP_DIR: /path/to/your/pseudopotentials
    ```

    See the [Pymatgen documentation](https://pymatgen.org/installation.html#potcar-setup) for full setup instructions.
