"""Inverse participation ratio (IPR) from VASP PROCAR."""

import os
import re
import numpy as np


def read_procar_header(filename="PROCAR"):
    """Read just the PROCAR header: (nkpts, nbands, natoms, nspin).

    Stops as soon as the counts are known, so it is cheap on a multi-GB file.
    `nspin` counts how many times the header repeats, which is 2 for ISPIN=2.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found")

    first = None
    nspin = 0
    with _open_procar(filename) as fh:
        for line in fh:
            if "k-points" in line and "bands" in line:
                nspin += 1
                if first is None:
                    nums = [int(x) for x in re.findall(r"\d+", line)]
                    if len(nums) < 3:
                        raise ValueError(
                            "PROCAR header does not contain k-points, bands, and ions")
                    first = nums[:3]
            elif first is not None and nspin >= 2:
                break

    if first is None:
        raise ValueError("Could not find PROCAR header with k-points/bands/ions")

    nkpts, nbands, natoms = first
    return nkpts, nbands, natoms, min(max(nspin, 1), 2)


def _open_procar(filename):
    if filename.endswith(".gz"):
        import gzip
        return gzip.open(filename, "rt", errors="replace")
    return open(filename, "r", errors="replace")


def read_procar(filename="PROCAR", wanted_bands=None):
    """Stream a PROCAR and extract per-atom charge projections.

    Returns (proj, energies, weights, nkpts, nbands, natoms, nspin) where
    proj[spin][ik][ib] is a list of natoms per-atom weights (empty for bands
    that were not requested) and energies[spin][ik][ib] is the eigenvalue.

    `wanted_bands` is a set of 0-indexed band numbers; anything else is skipped
    without being stored, which keeps memory flat on large files.

    Handles ISPIN=2, where the whole set of k-point blocks repeats, and
    non-collinear runs, where each band carries four ion blocks (charge, then
    mx/my/mz) of which only the first is the charge density.
    """
    nkpts, nbands, natoms, nspin = read_procar_header(filename)

    proj = [[[[] for _ in range(nbands)] for _ in range(nkpts)]
            for _ in range(nspin)]
    energies = [[[0.0 for _ in range(nbands)] for _ in range(nkpts)]
                for _ in range(nspin)]
    weights = [None] * nkpts

    re_weight = re.compile(r"weight\s*=\s*([-+0-9.eEdD]+)")

    kcount = -1
    isp = ik = ib = -1
    in_charge_block = False      # only the first ion block per band is charge
    overrun = bad_rows = bad_energies = 0

    with _open_procar(filename) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("k-point"):
                kcount += 1
                if kcount >= nspin * nkpts:
                    overrun += 1
                    isp = ik = -1          # refuse to wrap onto earlier data
                    continue
                isp, ik = divmod(kcount, nkpts)
                ib = -1
                in_charge_block = False
                m = re_weight.search(line)
                if m and isp == 0:
                    try:
                        weights[ik] = float(m.group(1).replace("D", "e").replace("d", "e"))
                    except ValueError:
                        pass
                continue

            if ik < 0:
                continue                    # inside an overrun block

            if line.startswith("band"):
                ib += 1
                in_charge_block = True
                parts = line.split()
                if len(parts) > 4 and 0 <= ib < nbands:
                    try:
                        energies[isp][ik][ib] = float(parts[4])
                    except ValueError:
                        bad_energies += 1
                continue

            # A "tot" line closes the charge block; later blocks are mx/my/mz
            # or phase factors and must not be mixed into the charge weights.
            if line.startswith("tot"):
                in_charge_block = False
                continue

            if (in_charge_block and line[0].isdigit()
                    and 0 <= ib < nbands
                    and (wanted_bands is None or ib in wanted_bands)):
                parts = line.split()
                try:
                    proj[isp][ik][ib].append(float(parts[-1]))
                except (ValueError, IndexError):
                    bad_rows += 1

    if overrun:
        print(f"Warning: PROCAR holds {overrun} more k-point block(s) than the "
              f"header declares ({nspin} x {nkpts}); the extra blocks were ignored.")
    if bad_rows:
        print(f"Warning: {bad_rows} ion row(s) could not be parsed and were skipped.")
    if bad_energies:
        print(f"Warning: {bad_energies} band energ(ies) could not be parsed "
              f"and are reported as 0.0 eV.")

    parsed = [w for w in weights if w is not None]
    if not parsed or any(w == 0 for w in parsed):
        # No weights, or a hybrid run where the band-structure k-points carry
        # weight 0.  Weighting there would silently drop exactly the k-points
        # of interest, so fall back to a plain mean.
        if parsed and any(w == 0 for w in parsed):
            print("Note: PROCAR contains zero-weight k-points; using an "
                  "unweighted average over k-points.")
        weights = [1.0 / nkpts] * nkpts
    else:
        mean_w = sum(parsed) / len(parsed)
        weights = [mean_w if w is None else w for w in weights]
        total = sum(weights)
        weights = [w / total for w in weights]

    return proj, energies, weights, nkpts, nbands, natoms, nspin


def compute_ipr_atomic(projections):
    """Compute atomic IPR from atomic projections."""
    projections = np.array(projections)
    tot = projections.sum()
    if tot < 1e-12:
        return 0.0
    weights = projections / tot
    return np.sum(weights ** 2)


def analyze_bands(proj, energies, weights, nkpts, band_indices,
                  nspin=1, verbose=True):
    """Compute k-averaged atomic IPR for selected bands.

    The average over k-points is weighted by the k-point weights, so an
    irreducible-wedge mesh is handled correctly.
    """
    results = []
    spin_labels = [None] if nspin == 1 else ["up", "dn"]

    for isp, spin_label in enumerate(spin_labels):
        for iband in band_indices:
            ipr_k = []
            e_k = []

            if verbose:
                tag = f"Band {iband}" if spin_label is None \
                    else f"Band {iband} (spin {spin_label})"
                print(f"\n{tag}")

            for ik in range(nkpts):
                ipr = compute_ipr_atomic(proj[isp][ik][iband - 1])
                ipr_k.append(ipr)
                e_k.append(energies[isp][ik][iband - 1])

                if verbose:
                    neff = (1 / ipr) if ipr > 0 else 0.0
                    print(
                        f"  k-point {ik + 1:3d}  "
                        f"E = {energies[isp][ik][iband - 1]:8.4f} eV  "
                        f"IPR = {ipr:8.4f}  "
                        f"N_eff = {neff:8.2f}"
                    )

            if ipr_k:
                w = np.asarray(weights[:len(ipr_k)], dtype=float)
                avg_ipr = float(np.average(np.asarray(ipr_k), weights=w))
                avg_e = float(np.average(np.asarray(e_k), weights=w))
            else:
                avg_ipr = avg_e = 0.0
            neff = (1 / avg_ipr) if avg_ipr > 0 else 0.0

            if verbose:
                print("  ------------------------------")
                print(f"  Avg IPR = {avg_ipr:.4f}")
                print(f"  N_eff  = {neff:.2f}")

            results.append((iband, spin_label, avg_e, avg_ipr, neff))

    return results


def save_results(results, filename="ipr_procar.dat"):
    """Save IPR results to a file."""
    with open(filename, "w") as f:
        has_spin = any(sp is not None for _, sp, _, _, _ in results)
        if has_spin:
            f.write("# Band  Spin  Energy(eV)   IPR   N_eff\n")
        else:
            f.write("# Band  Energy(eV)   IPR   N_eff\n")
        for band, spin, e, ipr, neff in results:
            if has_spin:
                f.write(f"{band:5d}  {spin or '-':>4s}  "
                        f"{e:10.6f}  {ipr:8.4f}  {neff:8.2f}\n")
            else:
                f.write(f"{band:5d}  {e:10.6f}  {ipr:8.4f}  {neff:8.2f}\n")

    print(f"\nResults written to {filename}")


def print_summary(results):
    """Print a compact IPR summary table."""
    has_spin = any(sp is not None for _, sp, _, _, _ in results)
    print("\nIPR summary")
    if has_spin:
        print("  Band  Spin    Energy(eV)      IPR      N_eff")
        for band, spin, e, ipr, neff in results:
            print(f"  {band:4d}  {spin or '-':>4s}    "
                  f"{e:10.6f}   {ipr:7.4f}   {neff:8.2f}")
    else:
        print("  Band    Energy(eV)      IPR      N_eff")
        for band, _spin, e, ipr, neff in results:
            print(f"  {band:4d}    {e:10.6f}   {ipr:7.4f}   {neff:8.2f}")


def _parse_band_indices(text):
    tokens = re.split(r"[\s,]+", text.strip())
    indices = []
    seen = set()

    for token in tokens:
        if not token:
            continue
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2:
                continue
            try:
                start = int(parts[0])
                end = int(parts[1])
            except ValueError:
                continue
            if start > end:
                start, end = end, start
            for i in range(start, end + 1):
                if i not in seen:
                    indices.append(i)
                    seen.add(i)
        else:
            try:
                i = int(token)
            except ValueError:
                continue
            if i not in seen:
                indices.append(i)
                seen.add(i)

    return indices


def _print_procar_info(nkpts, nbands, natoms, nspin):
    print("PROCAR info")
    print(f"  k-points : {nkpts}")
    print(f"  bands    : {nbands}")
    print(f"  atoms    : {natoms}")
    if nspin == 2:
        print("  spin     : collinear (2 channels)")


def _filter_band_indices(band_indices, nbands):
    """Keep valid 1-indexed band numbers and warn about skipped values."""
    filtered = [b for b in band_indices if 1 <= b <= nbands]
    if not filtered:
        raise ValueError(f"No bands in range 1..{nbands}.")

    if len(filtered) != len(band_indices):
        print("Warning: some bands were out of range and were skipped.")

    return filtered


def run_ipr(
    procar_file="PROCAR",
    band_text=None,
    output="ipr_procar.dat",
    show_details=False,
):
    """Run IPR analysis without prompting."""
    nkpts, nbands, natoms, nspin = read_procar_header(procar_file)
    _print_procar_info(nkpts, nbands, natoms, nspin)

    if not band_text:
        raise ValueError("No band indices provided.")

    band_indices = _parse_band_indices(band_text)
    if not band_indices:
        raise ValueError("No valid band indices found.")

    filtered = _filter_band_indices(band_indices, nbands)
    proj, energies, weights, nkpts, nbands, natoms, nspin = read_procar(
        procar_file, wanted_bands={b - 1 for b in filtered})
    results = analyze_bands(proj, energies, weights, nkpts, filtered,
                            nspin=nspin, verbose=show_details)
    if not show_details:
        print_summary(results)
    save_results(results, output)
    return results


def run_ipr_interactive(procar_file="PROCAR", output="ipr_procar.dat"):
    """Interactive IPR workflow."""
    try:
        nkpts, nbands, natoms, nspin = read_procar_header(procar_file)
    except Exception as e:
        print(f"Error: {e}")
        return

    _print_procar_info(nkpts, nbands, natoms, nspin)

    band_text = input("Band indices (e.g., 5 6 7 or 5-7): ").strip()
    if not band_text:
        print("No band indices provided.")
        return

    show_details = input("Show per-k-point values? [y/N]: ").strip().lower() == "y"

    try:
        band_indices = _parse_band_indices(band_text)
        if not band_indices:
            raise ValueError("No valid band indices found.")

        filtered = _filter_band_indices(band_indices, nbands)
        proj, energies, weights, nkpts, nbands, natoms, nspin = read_procar(
            procar_file, wanted_bands={b - 1 for b in filtered})
        results = analyze_bands(proj, energies, weights, nkpts, filtered,
                                nspin=nspin, verbose=show_details)
        if not show_details:
            print_summary(results)
        save_results(results, output)
    except Exception as e:
        print(f"Error: {e}")
