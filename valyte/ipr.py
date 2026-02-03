"""Inverse participation ratio (IPR) from VASP PROCAR."""

import os
import re
import numpy as np


def read_procar(filename="PROCAR"):
    """Read PROCAR and extract atomic total projections."""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found")

    with open(filename, "r") as f:
        lines = f.readlines()

    header_line = next((l for l in lines if "k-points" in l), None)
    if header_line is None:
        raise ValueError("Could not find PROCAR header with k-points/bands/ions")

    numbers = [int(x) for x in re.findall(r"\d+", header_line)]
    if len(numbers) < 3:
        raise ValueError("PROCAR header does not contain k-points, bands, and ions")

    nkpts, nbands, natoms = numbers[0], numbers[1], numbers[2]

    proj = [[[] for _ in range(nbands)] for _ in range(nkpts)]
    energies = [[0.0 for _ in range(nbands)] for _ in range(nkpts)]

    ik = -1
    ib = -1

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("k-point"):
            ik += 1
            ib = -1
            continue

        if line.startswith("band"):
            ib += 1
            parts = line.split()
            if len(parts) > 4:
                energies[ik][ib] = float(parts[4])
            continue

        if line[0].isdigit() and ik >= 0 and ib >= 0:
            parts = line.split()
            proj[ik][ib].append(float(parts[-1]))

    return proj, energies, nkpts, nbands, natoms


def compute_ipr_atomic(projections):
    """Compute atomic IPR from atomic projections."""
    projections = np.array(projections)
    tot = projections.sum()
    if tot < 1e-12:
        return 0.0
    weights = projections / tot
    return np.sum(weights ** 2)


def analyze_bands(proj, energies, nkpts, band_indices, verbose=True):
    """Compute k-averaged atomic IPR for selected bands."""
    results = []

    for iband in band_indices:
        ipr_k = []
        e_k = []

        if verbose:
            print(f"\nBand {iband}")

        for ik in range(nkpts):
            ipr = compute_ipr_atomic(proj[ik][iband - 1])
            ipr_k.append(ipr)
            e_k.append(energies[ik][iband - 1])

            if verbose:
                neff = (1 / ipr) if ipr > 0 else 0.0
                print(
                    f"  k-point {ik + 1:3d}  "
                    f"E = {energies[ik][iband - 1]:8.4f} eV  "
                    f"IPR = {ipr:8.4f}  "
                    f"N_eff = {neff:8.2f}"
                )

        avg_ipr = float(np.mean(ipr_k)) if ipr_k else 0.0
        avg_e = float(np.mean(e_k)) if e_k else 0.0
        neff = (1 / avg_ipr) if avg_ipr > 0 else 0.0

        if verbose:
            print("  ------------------------------")
            print(f"  Avg IPR = {avg_ipr:.4f}")
            print(f"  N_eff  = {neff:.2f}")

        results.append((iband, avg_e, avg_ipr, neff))

    return results


def save_results(results, filename="ipr_procar.dat"):
    """Save IPR results to a file."""
    with open(filename, "w") as f:
        f.write("# Band  Energy(eV)   IPR   N_eff\n")
        for band, e, ipr, neff in results:
            f.write(f"{band:5d}  {e:10.6f}  {ipr:8.4f}  {neff:8.2f}\n")

    print(f"\nResults written to {filename}")


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


def run_ipr_interactive():
    """Interactive IPR workflow."""
    procar_file = "PROCAR"

    try:
        proj, energies, nkpts, nbands, natoms = read_procar(procar_file)
    except Exception as e:
        print(f"Error: {e}")
        return

    print("PROCAR info")
    print(f"  k-points : {nkpts}")
    print(f"  bands    : {nbands}")
    print(f"  atoms    : {natoms}")

    band_text = input("Band indices (e.g., 5 6 7 or 5-7): ").strip()
    if not band_text:
        print("No band indices provided.")
        return

    band_indices = _parse_band_indices(band_text)
    if not band_indices:
        print("No valid band indices found.")
        return

    filtered = [b for b in band_indices if 1 <= b <= nbands]
    if not filtered:
        print(f"No bands in range 1..{nbands}.")
        return

    if len(filtered) != len(band_indices):
        print("Warning: some bands were out of range and were skipped.")

    show_details = input("Show per-k-point values? [y/N]: ").strip().lower() == "y"

    results = analyze_bands(proj, energies, nkpts, filtered, verbose=show_details)
    save_results(results, "ipr_procar.dat")
