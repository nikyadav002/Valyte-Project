"""Geometry optimization force convergence checker."""

import os
import re
import numpy as np


def parse_outcar(outcar_path="OUTCAR"):
    """
    Parse OUTCAR and return a list of ionic steps, each with:
        {'energy': float, 'max_force': float}
    """
    with open(outcar_path, "r") as f:
        content = f.read()

    nions_match = re.search(r"NIONS\s*=\s*(\d+)", content)
    if not nions_match:
        raise ValueError("Could not find NIONS in OUTCAR.")
    nions = int(nions_match.group(1))

    lines = content.splitlines()

    steps = []
    i = 0
    while i < len(lines):
        if "TOTAL-FORCE (eV/Angst)" in lines[i]:
            i += 1  # skip dashes line
            if i >= len(lines):
                break
            i += 1  # now on first force line

            forces = []
            for _ in range(nions):
                if i >= len(lines):
                    break
                parts = lines[i].split()
                if len(parts) >= 6:
                    fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                    forces.append(np.sqrt(fx**2 + fy**2 + fz**2))
                i += 1

            if not forces:
                continue

            max_force = max(forces)

            # Find the next "energy without entropy=" after this block
            energy = None
            j = i
            while j < len(lines) and j < i + 200:
                m = re.search(r"energy\s+without\s+entropy\s*=\s*([-\d.]+)", lines[j])
                if m:
                    energy = float(m.group(1))
                    break
                j += 1

            if energy is not None:
                steps.append({"energy": energy, "max_force": max_force})
        else:
            i += 1

    return steps


def check_convergence(outcar_path="OUTCAR", ediffg=None):
    """
    Print per-ionic-step force and energy convergence from OUTCAR
    and save data to valyte_forces.dat.

    Parameters
    ----------
    outcar_path : str
        Path to OUTCAR file.
    ediffg : float or None
        Force convergence threshold in eV/Å. If provided, marks converged steps
        and prints a convergence summary.
    """
    if not os.path.isfile(outcar_path):
        raise FileNotFoundError(f"No such file: {outcar_path}")

    steps = parse_outcar(outcar_path)
    if not steps:
        print("No ionic steps found in OUTCAR.")
        return

    e0 = steps[0]["energy"]

    header = f"{'Step':>5}  {'Max Force (eV/Å)':>18}  {'Energy (eV)':>18}  {'ΔE (eV)':>14}"
    print(header)
    print("-" * len(header))

    rows = []
    for idx, step in enumerate(steps):
        f_max = step["max_force"]
        e = step["energy"]
        de = e - e0
        converged = ediffg is not None and f_max < abs(ediffg)
        marker = " *" if converged else ""
        print(f"{idx:>5}  {f_max:>18.6f}  {e:>18.8f}  {de:>14.6g}{marker}")
        rows.append((idx, f_max, e, de))

    # Summary
    last = steps[-1]
    print()
    if ediffg is not None:
        threshold = abs(ediffg)
        if last["max_force"] < threshold:
            print(f"Converged: max force {last['max_force']:.6f} eV/Å < EDIFFG {threshold} eV/Å")
        else:
            print(f"NOT converged: max force {last['max_force']:.6f} eV/Å > EDIFFG {threshold} eV/Å")
    else:
        print(f"Final max force: {last['max_force']:.6f} eV/Å")
        print(f"Final energy:    {last['energy']:.8f} eV")

    dat_file = "valyte_forces.dat"
    with open(dat_file, "w") as f:
        f.write("# Step  MaxForce(eV/A)  Energy(eV)  DeltaE(eV)\n")
        for idx, f_max, e, de in rows:
            f.write(f"{idx:5d}  {f_max:18.6f}  {e:18.8f}  {de:14.6g}\n")
    print(f"Data saved to {dat_file}")
