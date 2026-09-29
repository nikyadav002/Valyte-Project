"""Inverse participation ratio (IPR) from VASP PROCAR."""

import os
import re
import numpy as np


def read_procar(filename="PROCAR"):
    """Read PROCAR and extract per-atom total projections.

    Returns (proj, energies, weights, nkpts, nbands, natoms, nspin) where
    proj[spin][ik][ib] is a list of natoms per-atom weights and
    energies[spin][ik][ib] is the eigenvalue in eV.

    Handles ISPIN=2 (the k-point blocks repeat for the second channel) and
    non-collinear runs, where each band carries four ion blocks (charge,
    mx, my, mz) of which only the first is the charge density.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found")

    with open(filename, "r") as f:
        lines = f.readlines()

    headers = [l for l in lines if "k-points" in l and "bands" in l]
    if not headers:
        raise ValueError("Could not find PROCAR header with k-points/bands/ions")

    numbers = [int(x) for x in re.findall(r"\d+", headers[0])]
    if len(numbers) < 3:
        raise ValueError("PROCAR header does not contain k-points, bands, and ions")

    nkpts, nbands, natoms = numbers[0], numbers[1], numbers[2]
    nspin = max(1, len(headers))

    proj = [[[[] for _ in range(nbands)] for _ in range(nkpts)]
            for _ in range(nspin)]
    energies = [[[0.0 for _ in range(nbands)] for _ in range(nkpts)]
                for _ in range(nspin)]
    weights = [1.0] * nkpts

    re_weight = re.compile(r"weight\s*=\s*([-+0-9.eEdD]+)")

    kcount = -1
    isp = 0
    ik = -1
    ib = -1

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("k-point"):
            kcount += 1
            isp = min(kcount // nkpts, nspin - 1)
            ik = kcount % nkpts
            ib = -1
            m = re_weight.search(line)
            if m and isp == 0:
                try:
                    weights[ik] = float(m.group(1).replace("D", "e"))
                except ValueError:
                    pass
            continue

        if line.startswith("band"):
            ib += 1
            parts = line.split()
            if len(parts) > 4 and 0 <= ib < nbands and ik >= 0:
                try:
                    energies[isp][ik][ib] = float(parts[4])
                except ValueError:
                    pass
            continue

        if line[0].isdigit() and ik >= 0 and 0 <= ib < nbands:
            bucket = proj[isp][ik][ib]
            # Keep only the first natoms ion rows: later blocks are the
            # magnetisation components (non-collinear) or phase factors.
            if len(bucket) < natoms:
                parts = line.split()
                try:
                    bucket.append(float(parts[-1]))
                except (ValueError, IndexError):
                    pass

    total_w = sum(weights)
    if total_w > 0:
        weights = [w / total_w for w in weights]
    else:
        weights = [1.0 / nkpts] * nkpts

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
    proj, energies, weights, nkpts, nbands, natoms, nspin = read_procar(procar_file)

    print("PROCAR info")
    print(f"  k-points : {nkpts}")
    print(f"  bands    : {nbands}")
    print(f"  atoms    : {natoms}")
    if nspin == 2:
        print(f"  spin     : collinear (2 channels)")

    if not band_text:
        raise ValueError("No band indices provided.")

    band_indices = _parse_band_indices(band_text)
    if not band_indices:
        raise ValueError("No valid band indices found.")

    filtered = _filter_band_indices(band_indices, nbands)
    results = analyze_bands(proj, energies, weights, nkpts, filtered,
                            nspin=nspin, verbose=show_details)
    if not show_details:
        print_summary(results)
    save_results(results, output)
    return results


def run_ipr_interactive(procar_file="PROCAR", output="ipr_procar.dat"):
    """Interactive IPR workflow."""
    try:
        proj, energies, weights, nkpts, nbands, natoms, nspin = read_procar(procar_file)
    except Exception as e:
        print(f"Error: {e}")
        return

    print("PROCAR info")
    print(f"  k-points : {nkpts}")
    print(f"  bands    : {nbands}")
    print(f"  atoms    : {natoms}")
    if nspin == 2:
        print(f"  spin     : collinear (2 channels)")

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
        results = analyze_bands(proj, energies, weights, nkpts, filtered,
                            nspin=nspin, verbose=show_details)
        if not show_details:
            print_summary(results)
        save_results(results, output)
    except Exception as e:
        print(f"Error: {e}")
