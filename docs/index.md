<p align="center">
  <img src="Logo.png" alt="Valyte Logo" width="100%"/>
</p>

# Valyte

**Publication-quality VASP pre- and post-processing from a single CLI.**

Valyte turns raw VASP output into clean, publication-ready plots and analysis — band structures, density of states, effective masses, convergence diagnostics, and more — all from one command-line tool.

---

## What can Valyte do?

### Pre-processing

Set up VASP calculations with minimal effort:

| Command | Description |
|---|---|
| [`valyte supercell`](preprocessing.md#supercell) | Generate supercells from POSCAR files |
| [`valyte kpt`](preprocessing.md#k-points-interactive-scf-grid) | Interactive KPOINTS generation (Monkhorst-Pack / Gamma) |
| [`valyte band kpt-gen`](band.md#1-generate-kpoints) | Automatic high-symmetry k-path (Bradley–Cracknell by default) |
| [`valyte potcar`](preprocessing.md#potcar) | Generate POTCAR from POSCAR species |

### Post-processing

Generate publication-quality analysis and figures:

| Command | Description |
|---|---|
| [`valyte dos`](dos.md) | Total and projected DOS with orbital resolution and gradient fills |
| [`valyte band`](band.md#2-standard-band-structure-plot) | Color-coded band structure with VBM aligned to 0 eV |
| [`valyte band --tricolor`](band.md#3-tricolor-orbital-resolved-plot) | Orbital-resolved tricolor band structure |
| [`valyte band --spin-resolved`](band.md#4-spin-resolved-band-structure-collinear) | Spin-polarized plot — spin-up and spin-down channels |
| [`valyte band --spin-texture`](band.md#5-non-collinear-spin-texture) | Non-collinear spin texture — bands colored by Sₓ, Sᵧ, or S_z |
| [`valyte ipr`](ipr.md) | Inverse Participation Ratio from PROCAR |
| [`valyte effmass`](effmass.md) | Carrier effective masses at VBM/CBM from parabolic fitting |
| [`valyte converge`](converge.md) | Ionic and SCF convergence monitor with multi-panel plots |

---

## Gallery

<p align="center">
  <img src="valyte_dos.png" alt="DOS Plot Example" width="47%"/>
  <img src="valyte_band.png" alt="Band Structure Example" width="38%"/>
</p>

<p align="center">
  <em>Left: Orbital-resolved density of states with gradient fills. Right: Color-coded band structure with VBM at 0 eV.</em>
</p>

---

## Quick start

```bash
pip install valyte
```

Then, from a directory containing your VASP output:

```bash
valyte dos                      # Plot density of states
valyte band                     # Plot band structure
valyte converge --forces        # Check relaxation convergence
```

→ **[Getting Started guide](getting-started.md)** for a complete walkthrough.

---

## Explore the documentation

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } **Getting Started**

    ---

    Installation, prerequisites, and your first plot in under two minutes.

    [:octicons-arrow-right-24: Get started](getting-started.md)

-   :material-chart-line:{ .lg .middle } **Band Structure**

    ---

    Standard, tricolor, spin-resolved, and spin-texture band plots.

    [:octicons-arrow-right-24: Band modes](band.md)

-   :material-chart-bell-curve-cumulative:{ .lg .middle } **Density of States**

    ---

    Total and projected DOS with orbital resolution and gradient fills.

    [:octicons-arrow-right-24: DOS plotting](dos.md)

-   :material-scale-balance:{ .lg .middle } **Effective Mass**

    ---

    Carrier effective masses from parabolic fitting at VBM/CBM.

    [:octicons-arrow-right-24: Effective mass](effmass.md)

-   :material-check-circle:{ .lg .middle } **Convergence**

    ---

    Multi-panel convergence monitoring for relaxations and SCF.

    [:octicons-arrow-right-24: Convergence](converge.md)

-   :material-console:{ .lg .middle } **CLI Reference**

    ---

    Every command and flag in one searchable page.

    [:octicons-arrow-right-24: CLI reference](cli-reference.md)

</div>
