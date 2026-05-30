"""Effective mass calculation from VASP band structure."""

import os
import numpy as np
from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.core import Spin


# ── Physical constants & precomputed conversion factor ──────────────────
# m* / m₀  =  ℏ² / (2·a · eV_to_J · Å_to_m²) / m₀
#   ℏ   = 1.0545718e-34  J·s
#   m₀  = 9.1093837e-31  kg
#   eV  = 1.6021766e-19  J
#   Å   = 1e-10 m
# ℏ² / (2·m₀·eV·Å²) = (1.0545718e-34)² / (2 * 9.1093837e-31 * 1.6021766e-19 * 1e-20)
#                     = 3.8099821
# Equivalently: m*/m₀ = 7.6199682 / (2a), where a is in eV·Å²
_EFFMASS_CONST = 3.8099821  # ℏ² / (2·m₀) in eV·Å²


def _resolve_vasprun(path):
    """Resolve path to a vasprun.xml file."""
    if os.path.isdir(path):
        path = os.path.join(path, "vasprun.xml")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"vasprun.xml not found: {path}")
    return path


def _resolve_kpoints(kpoints_path, vasprun_path):
    """Auto-detect KPOINTS file if not explicitly given."""
    if kpoints_path:
        return kpoints_path
    base_dir = os.path.dirname(os.path.abspath(vasprun_path))
    candidate = os.path.join(base_dir, "KPOINTS")
    if os.path.isfile(candidate):
        return candidate
    return None


def _get_label_for_kpoint(bs, kpt_index):
    """Get the high-symmetry label for a k-point index, or empty string."""
    # bs.kpoints is a list of Kpoint objects; each has a .label attribute
    kpt = bs.kpoints[kpt_index]
    if kpt.label:
        label = kpt.label
        # Replace GAMMA with the unicode symbol
        if label.upper() in ("GAMMA", "\\GAMMA", "\\\\GAMMA"):
            return "Γ"
        return label
    return ""


def _get_segment_label(bs, branch, extremum_index):
    """Build direction label like 'Γ → X' for a branch relative to an extremum.

    Returns (label_string, sign) where sign indicates direction:
      +1 means the segment goes forward from the extremum,
      -1 means the segment ends at the extremum (backward direction).
    """
    start_idx = branch["start_index"]
    end_idx = branch["end_index"]

    start_label = _get_label_for_kpoint(bs, start_idx) or "?"
    end_label = _get_label_for_kpoint(bs, end_idx) or "?"

    if extremum_index == start_idx:
        # Extremum is at the start → segment goes forward
        return f"{start_label} → {end_label}", +1
    elif extremum_index == end_idx:
        # Extremum is at the end → we fit backwards from the end
        return f"{end_label} → {start_label}", -1
    else:
        # Extremum is interior to the segment
        # Determine which end is closer for labeling
        return f"{start_label} → {end_label}", 0


def _cartesian_k_distances(bs, kpt_indices):
    """Convert k-point indices to Cartesian coordinates and return them.

    Returns an array of Cartesian coordinates in 1/Å, shape (N, 3).
    """
    rec_lattice = bs.lattice_rec
    coords = []
    for idx in kpt_indices:
        frac = bs.kpoints[idx].frac_coords
        cart = rec_lattice.get_cartesian_coords(frac)
        coords.append(cart)
    return np.array(coords)


def _signed_k_distance(cart_coords, extremum_pos_in_array):
    """Compute signed distance along a 1D k-path from the extremum.

    Parameters
    ----------
    cart_coords : ndarray, shape (N, 3)
        Cartesian coordinates of k-points in the fitting window.
    extremum_pos_in_array : int
        Index of the extremum within cart_coords.

    Returns
    -------
    k_dist : ndarray, shape (N,)
        Signed distances in 1/Å.  Negative before the extremum, positive after.
    """
    origin = cart_coords[extremum_pos_in_array]
    diffs = cart_coords - origin

    # Cumulative distance along the path for sign assignment
    k_dist = np.zeros(len(cart_coords))
    for i in range(len(cart_coords)):
        dist = np.linalg.norm(diffs[i])
        if i < extremum_pos_in_array:
            k_dist[i] = -dist
        elif i > extremum_pos_in_array:
            k_dist[i] = dist
        # i == extremum_pos_in_array → 0.0

    return k_dist


def _parabolic_fit(k_dist, energies):
    """Fit E(k) = a·k² + b·k + c and return (a, r_squared, fit_curve).

    Parameters
    ----------
    k_dist : ndarray
        Signed k-distances in 1/Å.
    energies : ndarray
        Band energies in eV at those k-points.

    Returns
    -------
    a : float
        Quadratic coefficient in eV·Å².
    r_squared : float
        Coefficient of determination.
    k_fine : ndarray
        Dense k-grid for plotting the fit curve.
    e_fine : ndarray
        Fitted energies on the dense grid.
    """
    coeffs = np.polyfit(k_dist, energies, 2)
    a, b, c = coeffs

    # R²
    e_pred = np.polyval(coeffs, k_dist)
    ss_res = np.sum((energies - e_pred) ** 2)
    ss_tot = np.sum((energies - np.mean(energies)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0

    # Dense curve for plotting
    k_fine = np.linspace(k_dist.min(), k_dist.max(), 200)
    e_fine = np.polyval(coeffs, k_fine)

    return a, r_squared, k_fine, e_fine


def _find_branches_for_kpt(bs, kpt_index):
    """Find all branches that contain a given k-point index.

    Returns list of (branch_dict, position) where position is:
      'start', 'end', or 'interior'.
    """
    results = []
    for branch in bs.branches:
        s = branch["start_index"]
        e = branch["end_index"]
        if kpt_index == s:
            results.append((branch, "start"))
        elif kpt_index == e:
            results.append((branch, "end"))
        elif s < kpt_index < e:
            results.append((branch, "interior"))
    return results


def _do_one_fit(bs, kpt_indices, extremum_pos, band_idx, spin, efermi,
                direction_label):
    """Low-level: fit a parabola to one set of k-points around an extremum.

    Returns a dict with fit results, or None if insufficient points.
    """
    if len(kpt_indices) < 3:
        return None  # need at least 3 points for a quadratic fit

    # Get energies for this band at these k-points (shifted by Fermi)
    bands = bs.bands[spin]  # shape: (nbands, nkpts)
    energies = np.array([bands[band_idx][ki] - efermi for ki in kpt_indices])

    # Get Cartesian k-distances
    cart_coords = _cartesian_k_distances(bs, kpt_indices)
    k_dist = _signed_k_distance(cart_coords, extremum_pos)

    # Perform the fit
    a, r_squared, k_fine, e_fine = _parabolic_fit(k_dist, energies)

    # Compute m*/m₀
    if abs(a) < 1e-12:
        m_star = float("inf")
    else:
        m_star = _EFFMASS_CONST / a

    # Warning for non-parabolic bands
    warning = None
    if r_squared < 0.95:
        warning = "Non-parabolic band — effective mass is approximate"

    spin_label = None
    if Spin.down in bs.bands:
        spin_label = "up" if spin == Spin.up else "down"

    return {
        "band_index": band_idx,
        "spin": spin_label,
        "direction": direction_label,
        "m_star": m_star,
        "r_squared": r_squared,
        "k_fit": k_fine,
        "e_fit": e_fine,
        "k_data": k_dist,
        "e_data": energies,
        "warning": warning,
    }


def _fits_for_branch(bs, branch, extremum_index, band_idx, spin, npoints,
                     efermi):
    """Generate directional parabolic fits along a branch.

    For boundary extrema (start / end of the branch) a single fit toward the
    opposite end is produced.  For an interior extremum two separate fits are
    produced — one toward each end of the branch — because the curvature can
    differ in the two k-directions.

    Returns a *list* of fit-result dicts (0, 1, or 2 items).
    """
    start_idx = branch["start_index"]
    end_idx = branch["end_index"]

    start_label = _get_label_for_kpoint(bs, start_idx) or "?"
    end_label = _get_label_for_kpoint(bs, end_idx) or "?"
    ext_label = _get_label_for_kpoint(bs, extremum_index)
    if not ext_label:
        frac = bs.kpoints[extremum_index].frac_coords
        ext_label = f"({frac[0]:.2f},{frac[1]:.2f},{frac[2]:.2f})"

    fit_results = []

    if extremum_index == start_idx:
        # ── Extremum at start → fit forward toward end ──────────────
        fit_end = min(start_idx + npoints, end_idx)
        kpt_indices = list(range(start_idx, fit_end + 1))
        direction = f"{ext_label} → {end_label}"
        r = _do_one_fit(bs, kpt_indices, 0, band_idx, spin, efermi,
                        direction)
        if r is not None:
            fit_results.append(r)

    elif extremum_index == end_idx:
        # ── Extremum at end → fit backward toward start ─────────────
        fit_start = max(end_idx - npoints, start_idx)
        kpt_indices = list(range(fit_start, end_idx + 1))
        direction = f"{ext_label} → {start_label}"
        r = _do_one_fit(bs, kpt_indices, len(kpt_indices) - 1,
                        band_idx, spin, efermi, direction)
        if r is not None:
            fit_results.append(r)

    else:
        # ── Interior extremum → two directional fits ────────────────
        # Forward (toward end)
        fwd_end = min(extremum_index + npoints, end_idx)
        fwd_indices = list(range(extremum_index, fwd_end + 1))
        if len(fwd_indices) >= 3:
            direction = f"→ {end_label}"
            r = _do_one_fit(bs, fwd_indices, 0, band_idx, spin, efermi,
                            direction)
            if r is not None:
                fit_results.append(r)

        # Backward (toward start)
        bwd_start = max(extremum_index - npoints, start_idx)
        bwd_indices = list(range(bwd_start, extremum_index + 1))
        if len(bwd_indices) >= 3:
            direction = f"→ {start_label}"
            r = _do_one_fit(bs, bwd_indices, len(bwd_indices) - 1,
                            band_idx, spin, efermi, direction)
            if r is not None:
                fit_results.append(r)

    return fit_results


def compute_effective_masses(vasprun_path=".", kpoints_path=None,
                             npoints=3, band_index=None, tol=1e-3):
    """Compute carrier effective masses from a VASP band structure.

    Parameters
    ----------
    vasprun_path : str
        Path to vasprun.xml or directory containing it.
    kpoints_path : str or None
        Path to KPOINTS file for high-symmetry labels.
    npoints : int
        Number of k-points on each side of extremum for fitting.
    band_index : list of int or None
        Manual 1-indexed band indices to fit. If None, auto-detect VBM/CBM.
    tol : float
        Degeneracy tolerance in eV.

    Returns
    -------
    dict
        Results containing band gap, VBM/CBM info, and effective masses.
    """
    vasprun_path = _resolve_vasprun(vasprun_path)
    kpoints_path = _resolve_kpoints(kpoints_path, vasprun_path)

    vr = BSVasprun(vasprun_path, parse_projected_eigen=False)
    bs = vr.get_band_structure(kpoints_filename=kpoints_path, line_mode=True)

    # Check for line-mode calculation
    if not bs.branches:
        raise ValueError(
            "No k-path branches found. This appears to be a self-consistent "
            "(non-line-mode) calculation.\n"
            "Effective mass calculation requires a line-mode band structure "
            "(line-mode KPOINTS)."
        )

    # Check for metals
    if bs.is_metal():
        return {"is_metal": True}

    efermi = bs.efermi
    vbm_data = bs.get_vbm()
    cbm_data = bs.get_cbm()

    band_gap = bs.get_band_gap()
    gap_value = band_gap["energy"]
    gap_type = "direct" if band_gap["direct"] else "indirect"

    # Spin channels to process
    spin_channels = list(bs.bands.keys())
    is_spin_polarized = len(spin_channels) > 1

    results = {
        "band_gap": gap_value,
        "gap_type": gap_type,
        "vbm": {},
        "cbm": {},
        "hole_masses": [],
        "electron_masses": [],
        "is_metal": False,
    }

    # ── Process VBM and CBM ────────────────────────────────────────────
    for extremum_type, ext_data, mass_key in [
        ("vbm", vbm_data, "hole_masses"),
        ("cbm", cbm_data, "electron_masses"),
    ]:
        ext_energy = ext_data["energy"]
        ext_kpoint_index = ext_data["kpoint_index"]

        # ext_data["kpoint_index"] can be a list (one per spin or degenerate)
        # ext_data["band_index"] is a dict: {Spin: [band_indices]}

        # Pick the first k-point index for reporting
        if isinstance(ext_kpoint_index, (list, tuple)):
            primary_kpt_idx = ext_kpoint_index[0]
        else:
            primary_kpt_idx = ext_kpoint_index

        kpt_frac = tuple(np.round(bs.kpoints[primary_kpt_idx].frac_coords, 4))
        kpt_label = _get_label_for_kpoint(bs, primary_kpt_idx)

        # Collect all (spin, band_index, kpt_index) combinations to fit.
        # Crucially, we iterate over ALL k-point indices where the extremum
        # occurs, not just the first.  High-symmetry points like Γ appear
        # at the boundary of multiple k-path segments with distinct indices.
        fit_targets = []

        all_kpt_indices = (list(ext_kpoint_index)
                          if isinstance(ext_kpoint_index, (list, tuple))
                          else [ext_kpoint_index])

        if band_index is not None:
            # Manual band indices (1-indexed from user → 0-indexed)
            for bi_user in band_index:
                bi = bi_user - 1
                for spin in spin_channels:
                    for ki in all_kpt_indices:
                        fit_targets.append((spin, bi, ki))
        else:
            # Auto-detect: use pymatgen's band_index dict
            band_idx_dict = ext_data["band_index"]
            for spin in spin_channels:
                if spin not in band_idx_dict:
                    continue

                band_indices_for_spin = set(band_idx_dict[spin])

                for bi in band_idx_dict[spin]:
                    # Add fit targets for ALL k-point indices at this energy
                    for ki in all_kpt_indices:
                        e = bs.bands[spin][bi][ki]
                        if abs(e - ext_energy) < tol:
                            fit_targets.append((spin, bi, ki))

                    # Also check for degenerate bands not in pymatgen's list
                    for other_bi in range(bs.nb_bands):
                        if other_bi in band_indices_for_spin:
                            continue
                        for ki in all_kpt_indices:
                            e_other = bs.bands[spin][other_bi][ki]
                            if abs(e_other - ext_energy) < tol:
                                fit_targets.append((spin, other_bi, ki))

        # Deduplicate targets
        seen = set()
        unique_targets = []
        for target in fit_targets:
            if target not in seen:
                seen.add(target)
                unique_targets.append(target)
        fit_targets = unique_targets

        # Report the primary band index (first one, 1-indexed)
        primary_band_indices = []
        for spin in spin_channels:
            if spin in ext_data["band_index"]:
                primary_band_indices.extend(ext_data["band_index"][spin])
        if primary_band_indices:
            primary_band_display = primary_band_indices[0] + 1  # 1-indexed
        else:
            primary_band_display = "?"

        results[extremum_type] = {
            "kpoint_frac": kpt_frac,
            "kpoint_label": kpt_label,
            "band_index_display": primary_band_display,
            "energy": ext_energy,
        }

        # ── Fit along each branch for each target ────────────────────
        seen_directions = set()
        for spin, bi, kpt_idx in fit_targets:
            branches = _find_branches_for_kpt(bs, kpt_idx)
            if not branches:
                continue

            for branch, _position in branches:
                for fit_result in _fits_for_branch(
                    bs, branch, kpt_idx, bi, spin, npoints, efermi
                ):
                    # Deduplicate by (band, spin, direction) to avoid
                    # reporting the same physical direction twice when the
                    # same k-point appears at multiple indices.
                    dedup_key = (bi, fit_result["spin"],
                                 fit_result["direction"])
                    if dedup_key in seen_directions:
                        continue
                    seen_directions.add(dedup_key)

                    fit_result["band_index_display"] = bi + 1
                    results[mass_key].append(fit_result)

    return results


def print_results(results):
    """Print effective mass results to stdout."""
    if results.get("is_metal", False):
        print("Effective mass")
        print("═" * 14)
        print()
        print("  This system is metallic (no band gap).")
        print("  Effective mass calculation requires a semiconductor/insulator.")
        return

    print("Effective mass")
    print("═" * 14)
    print()

    gap = results["band_gap"]
    gap_type = results["gap_type"]
    vbm = results["vbm"]
    cbm = results["cbm"]

    print(f"  Band gap:     {gap:.3f} eV ({gap_type})")

    kv = vbm["kpoint_frac"]
    lv = f"  [{vbm['kpoint_label']}]" if vbm["kpoint_label"] else ""
    print(f"  VBM:          k = ({kv[0]:.3f}, {kv[1]:.3f}, {kv[2]:.3f}){lv}"
          f"    Band {vbm['band_index_display']}")

    kc = cbm["kpoint_frac"]
    lc = f"  [{cbm['kpoint_label']}]" if cbm["kpoint_label"] else ""
    print(f"  CBM:          k = ({kc[0]:.3f}, {kc[1]:.3f}, {kc[2]:.3f}){lc}"
          f"    Band {cbm['band_index_display']}")

    # ── Hole masses ─────────────────────────────────────────────────
    hole_masses = results["hole_masses"]
    if hole_masses:
        # Group by spin if spin-polarized
        spins_present = sorted(set(m["spin"] for m in hole_masses if m["spin"]),
                               key=lambda s: (s != "up", s))

        if spins_present:
            for spin_label in spins_present:
                spin_display = "Spin Up" if spin_label == "up" else "Spin Down"
                print(f"\n  Hole effective masses (VBM) — {spin_display}:")
                for m in hole_masses:
                    if m["spin"] != spin_label:
                        continue
                    _print_mass_line(m)
        else:
            print("\n  Hole effective masses (VBM):")
            for m in hole_masses:
                _print_mass_line(m)

    # ── Electron masses ─────────────────────────────────────────────
    electron_masses = results["electron_masses"]
    if electron_masses:
        spins_present = sorted(set(m["spin"] for m in electron_masses if m["spin"]),
                               key=lambda s: (s != "up", s))

        if spins_present:
            for spin_label in spins_present:
                spin_display = "Spin Up" if spin_label == "up" else "Spin Down"
                print(f"\n  Electron effective masses (CBM) — {spin_display}:")
                for m in electron_masses:
                    if m["spin"] != spin_label:
                        continue
                    _print_mass_line(m)
        else:
            print("\n  Electron effective masses (CBM):")
            for m in electron_masses:
                _print_mass_line(m)

    print()


def _print_mass_line(m):
    """Print a single effective mass result line."""
    band_suffix = ""
    # Show band index if multiple bands are being reported
    direction = m["direction"]
    m_star = m["m_star"]
    r_sq = m["r_squared"]

    if abs(m_star) > 1e6:
        print(f"    {direction:12s}:     m* = {'∞':>8s} m₀")
    else:
        print(f"    {direction:12s}:     m* = {m_star:>8.3f} m₀")

    if m.get("warning"):
        print(f"      ⚠ {m['warning']}")


def save_results_dat(results, filepath="valyte_effmass.dat"):
    """Save effective mass results to a text file."""
    if results.get("is_metal", False):
        with open(filepath, "w") as f:
            f.write("# Valyte Effective Mass Results\n")
            f.write("# System is metallic — no band gap found.\n")
        print(f"Saved data: {filepath}")
        return

    gap = results["band_gap"]
    gap_type = results["gap_type"]
    vbm = results["vbm"]
    cbm = results["cbm"]

    kv = vbm["kpoint_frac"]
    kc = cbm["kpoint_frac"]
    lv = f" [{vbm['kpoint_label']}]" if vbm["kpoint_label"] else ""
    lc = f" [{cbm['kpoint_label']}]" if cbm["kpoint_label"] else ""

    with open(filepath, "w") as f:
        f.write("# Valyte Effective Mass Results\n")
        f.write(f"# Band gap: {gap:.3f} eV ({gap_type})\n")
        f.write(f"# VBM: k = ({kv[0]:.3f}, {kv[1]:.3f}, {kv[2]:.3f}){lv}"
                f"  Band {vbm['band_index_display']}\n")
        f.write(f"# CBM: k = ({kc[0]:.3f}, {kc[1]:.3f}, {kc[2]:.3f}){lc}"
                f"  Band {cbm['band_index_display']}\n")
        f.write("#\n")
        f.write(f"# {'Type':<6s}  {'Band':>4s}  {'Spin':<4s}  "
                f"{'Direction':<14s}  {'m*/m0':>10s}  {'R²':>8s}\n")

        for m in results["hole_masses"]:
            spin_str = m["spin"] if m["spin"] else "-"
            m_str = f"{m['m_star']:.4f}" if abs(m["m_star"]) < 1e6 else "inf"
            f.write(f"  {'hole':<6s}  {m['band_index_display']:>4d}  "
                    f"{spin_str:<4s}  {m['direction']:<14s}  "
                    f"{m_str:>10s}  {m['r_squared']:>8.4f}\n")

        for m in results["electron_masses"]:
            spin_str = m["spin"] if m["spin"] else "-"
            m_str = f"{m['m_star']:.4f}" if abs(m["m_star"]) < 1e6 else "inf"
            f.write(f"  {'elec':<6s}  {m['band_index_display']:>4d}  "
                    f"{spin_str:<4s}  {m['direction']:<14s}  "
                    f"{m_str:>10s}  {m['r_squared']:>8.4f}\n")

    print(f"Saved data: {filepath}")
