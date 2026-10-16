# Pre-processing

Valyte provides three pre-processing utilities for setting up VASP calculations — supercell generation, k-point grids, and automatic POTCAR generation.

!!! tip "Band structure k-paths"
    For line-mode KPOINTS (band structure calculations), see [`valyte band kpt-gen`](band.md#1-generate-kpoints) instead.

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

## K-Points — SCF grid

Generate a `KPOINTS` file for SCF or relaxation calculations.

```bash
valyte kpt
```

With no flags, the command runs interactively and prompts for:

1. **K-mesh scheme** — Monkhorst-Pack or Gamma
2. **K-spacing** — in units of 2π/Å (e.g. `0.04`)

The optimal k-grid dimensions are automatically calculated from the lattice vectors in your `POSCAR`.

For batch scripts or cluster workflows, provide the settings directly:

```bash
valyte kpt --spacing 0.04 --scheme gamma
valyte kpt -i POSCAR_relaxed -o KPOINTS_scf --spacing 0.03 --scheme mp
valyte kpt --spacing 0.05 --no-potcar
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | `POSCAR` | Input POSCAR file |
| `-o`, `--output` | `KPOINTS` | Output KPOINTS file |
| `--spacing` | `0.04` | K-spacing in 2π/Å |
| `--scheme` | `gamma` | K-mesh scheme: `gamma` or `mp` |
| `--no-potcar` | off | Skip automatic POTCAR generation |

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
| `--functional` | `PBE_64` | POTCAR set |
| `--potcar` | — | Per-element overrides, e.g. `Fe=Fe_pv O=O_h` |

Each element gets the potential VASP recommends in the [potpaw.64 PBE table](https://www.vasp.at/wiki/index.php/Available_pseudopotentials) (for example `Ti_sv`, `Ga_d`, `Ba_sv`, `Gd_3`); elements without a recommended variant use the plain symbol. The choice for each element is printed when the POTCAR is written.

**Available functionals:** `PBE_64`, `PBE_54`, `PBE_52`, `PBE`, `LDA_64`, `LDA_54`, `LDA`, `PW91`, and others supported by Pymatgen.

```bash
# VASP-recommended potentials, PBE_64
valyte potcar

# Use Fe_pv instead of the recommended Fe
valyte potcar --potcar Fe=Fe_pv

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

    See the [pymatgen documentation](https://pymatgen.org/installation.html#potcar-setup) for full setup instructions, or the [FAQ → POTCAR errors](faq.md#valueerror-no-potcar-for-element--or-pymatgen-errors) for troubleshooting.
