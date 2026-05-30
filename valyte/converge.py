"""VASP convergence monitoring — OSZICAR/OUTCAR parsing and plotting."""

import gzip
import os
import re
import sys

import matplotlib as mpl
mpl.use("agg")
mpl.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np


# ── File helpers ──────────────────────────────────────────────────────────────

def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "r", errors="replace")


def _find_file(directory, name):
    for candidate in (name, name + ".gz"):
        p = os.path.join(directory, candidate)
        if os.path.isfile(p):
            return p
    return None


# ── OSZICAR parser ────────────────────────────────────────────────────────────

_RE_ELEC = re.compile(
    r"^\s*(?:DAV|RMM|CG|DIAG)\s*:\s*(\d+)\s+"
    r"([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)\s+"
    r"([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)",
    re.IGNORECASE,
)

_RE_IONIC = re.compile(
    r"^\s*(\d+)\s+F=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r"\s+E0=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r"\s+d\s*E\s*=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
)

_RE_IONIC_MD = re.compile(
    r"^\s*(\d+)\s+T=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r".*?E=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r"\s+F=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r"\s+E0=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r".*?EK=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)",
    re.DOTALL,
)

_RE_MAG = re.compile(r"mag=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)")


def _float(s):
    return float(s.replace("D", "e").replace("d", "e"))


def parse_oszicar(path):
    """Parse an OSZICAR file.

    Returns a list of ionic step dicts, each containing:
        number      int    VASP ionic step counter (1-based)
        F           float  free energy (eV)
        E0          float  energy sigma→0 (eV)
        dE          float  energy change (eV)
        mag         float or None
        is_md       bool
        T           float or None  (MD temperature K)
        EK          float or None  (kinetic energy eV)
        elec_steps  list of {step, E, dE}

    Trailing electronic steps with no ionic summary are attached to a
    synthetic entry with _incomplete=True and ionic fields set to None.
    """
    steps = []
    pending_elec = []
    prev_ionic_number = 0

    with _open(path) as fh:
        for raw in fh:
            line = raw.rstrip()

            # MD ionic line (more fields — try first)
            m = _RE_IONIC_MD.match(line)
            if m:
                step = {
                    "number": int(m.group(1)),
                    "F":  _float(m.group(4)),
                    "E0": _float(m.group(5)),
                    "dE": _float(m.group(4)) - _float(m.group(3)),
                    "mag": None,
                    "is_md": True,
                    "T":  _float(m.group(2)),
                    "EK": _float(m.group(6)),
                    "elec_steps": list(pending_elec),
                }
                mm = _RE_MAG.search(line)
                if mm:
                    step["mag"] = _float(mm.group(1))
                pending_elec.clear()
                steps.append(step)
                prev_ionic_number = step["number"]
                continue

            # Regular ionic summary
            m = _RE_IONIC.match(line)
            if m:
                step = {
                    "number": int(m.group(1)),
                    "F":  _float(m.group(2)),
                    "E0": _float(m.group(3)),
                    "dE": _float(m.group(4)),
                    "mag": None,
                    "is_md": False,
                    "T":  None,
                    "EK": None,
                    "elec_steps": list(pending_elec),
                }
                mm = _RE_MAG.search(line)
                if mm:
                    step["mag"] = _float(mm.group(1))
                pending_elec.clear()
                steps.append(step)
                prev_ionic_number = step["number"]
                continue

            # Electronic step
            m = _RE_ELEC.match(line)
            if m:
                pending_elec.append({
                    "step": int(m.group(1)),
                    "E":   _float(m.group(2)),
                    "dE":  abs(_float(m.group(3))),
                })

    if pending_elec:
        steps.append({
            "number": prev_ionic_number + 1,
            "F": None, "E0": None, "dE": None,
            "mag": None, "is_md": False, "T": None, "EK": None,
            "elec_steps": list(pending_elec),
            "_incomplete": True,
        })

    return steps


# ── OUTCAR parsers ────────────────────────────────────────────────────────────

_RE_IBRION = re.compile(r"IBRION\s*=\s*(-?\d+)")
_RE_NSW    = re.compile(r"NSW\s*=\s*(\d+)")
_RE_EDIFF  = re.compile(r"EDIFF\s*=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)")
_RE_EDIFFG = re.compile(r"EDIFFG\s*=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)")
_RE_NIONS  = re.compile(r"NIONS\s*=\s*(\d+)")
_RE_WTIME  = re.compile(r"Total CPU time used \(sec\):\s*([\d.]+)")
_RE_WTIME2 = re.compile(r"Elapsed time \(sec\):\s*([\d.]+)")


def parse_outcar_header(path):
    """Read only the INCAR-echo section of OUTCAR to get job parameters.

    Returns dict: ibrion, nsw, ediff, ediffg, nions, walltime_s.
    """
    info = {"ibrion": None, "nsw": None, "ediff": None,
            "ediffg": None, "nions": None, "walltime_s": None}
    found_header = False

    with _open(path) as fh:
        for line in fh:
            if "INCAR" in line and "input" in line.lower():
                found_header = True

            if info["ibrion"] is None:
                m = _RE_IBRION.search(line)
                if m:
                    info["ibrion"] = int(m.group(1))

            if info["nsw"] is None:
                m = _RE_NSW.search(line)
                if m:
                    info["nsw"] = int(m.group(1))

            if info["ediff"] is None:
                m = _RE_EDIFF.search(line)
                if m:
                    info["ediff"] = _float(m.group(1))

            if info["ediffg"] is None:
                m = _RE_EDIFFG.search(line)
                if m:
                    info["ediffg"] = _float(m.group(1))

            if info["nions"] is None:
                m = _RE_NIONS.search(line)
                if m:
                    info["nions"] = int(m.group(1))

            m = _RE_WTIME.search(line)
            if m:
                info["walltime_s"] = float(m.group(1))

            if info["walltime_s"] is None:
                m = _RE_WTIME2.search(line)
                if m:
                    info["walltime_s"] = float(m.group(1))

            if (found_header and all(info[k] is not None
                    for k in ("ibrion", "nsw", "ediff", "nions"))):
                break

    return info


def parse_outcar_forces(path, nions):
    """Single-pass OUTCAR parse for per-ionic-step max force and pressure.

    Returns list of {max_force, pressure}, one per ionic step.
    """
    re_force_hdr = re.compile(r"TOTAL-FORCE \(eV/Angst\)")
    re_pressure  = re.compile(
        r"external pressure\s*=\s*([-+]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*kB"
    )
    re_sep = re.compile(r"^-{10,}")

    steps = []
    in_block = False
    skip_sep = False
    buf = []
    pending_p = None

    with _open(path) as fh:
        for line in fh:
            if re_force_hdr.search(line):
                in_block = True
                skip_sep = True
                buf = []
                continue

            if in_block:
                if skip_sep:
                    skip_sep = False
                    continue
                if re_sep.match(line.strip()):
                    in_block = False
                    max_f = max(buf) if buf else None
                    steps.append({"max_force": max_f, "pressure": pending_p})
                    pending_p = None
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                        buf.append((fx*fx + fy*fy + fz*fz) ** 0.5)
                    except ValueError:
                        pass
                continue

            m = re_pressure.search(line)
            if m:
                pending_p = float(m.group(1))

    return steps


# ── INCAR fallback ────────────────────────────────────────────────────────────

def _parse_incar(path):
    info = {}
    re_tag = re.compile(r"^\s*([A-Z]+\w*)\s*=\s*(.+?)(?:!.*)?$")
    with _open(path) as fh:
        for line in fh:
            m = re_tag.match(line)
            if not m:
                continue
            tag, val = m.group(1).strip(), m.group(2).strip()
            try:
                if tag == "IBRION":
                    info["ibrion"] = int(val)
                elif tag == "NSW":
                    info["nsw"] = int(val)
                elif tag == "EDIFF":
                    info["ediff"] = _float(val.split()[0])
                elif tag == "EDIFFG":
                    info["ediffg"] = _float(val.split()[0])
            except (ValueError, IndexError):
                pass
    return info


# ── Utilities ─────────────────────────────────────────────────────────────────

def _fmt_time(seconds):
    s = int(seconds)
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    if h:
        return f"{h}h {m:02d}m {s:02d}s"
    if m:
        return f"{m}m {s:02d}s"
    return f"{s}s"


def _calc_type(ibrion, nsw):
    if ibrion is None:
        return "Unknown"
    if ibrion == 0:
        return "MD"
    if ibrion == -1 or nsw == 0:
        return "Single-point"
    return f"Relaxation (IBRION = {ibrion}, NSW = {nsw})"


# ── Terminal summary ──────────────────────────────────────────────────────────

def print_summary(steps, outcar_info, force_steps=None, fthresh=0.02):
    ionic = [s for s in steps if not s.get("_incomplete")]
    ibrion = outcar_info.get("ibrion")
    nsw    = outcar_info.get("nsw")
    ediffg = outcar_info.get("ediffg")
    wtime  = outcar_info.get("walltime_s")
    is_md  = any(s["is_md"] for s in ionic) if ionic else False

    print("Convergence")
    print("═" * 11)
    print()
    print(f"  Calculation type:   {_calc_type(ibrion, nsw)}")

    n_done = len(ionic)
    if nsw and not is_md:
        print(f"  Ionic steps:        {n_done} / {nsw}")
    else:
        print(f"  Ionic steps:        {n_done}")

    if ionic:
        last = ionic[-1]
        if not is_md:
            if ediffg is not None and ediffg < 0:
                fok = (force_steps and force_steps[-1]["max_force"] is not None
                       and force_steps[-1]["max_force"] < abs(ediffg))
                mark = "✓ Converged" if fok else "✗ Not converged"
            elif ediffg is not None and ediffg > 0:
                eok = last["dE"] is not None and abs(last["dE"]) < ediffg
                mark = "✓ Converged" if eok else "✗ Not converged"
            elif nsw and n_done < nsw:
                mark = "Running / incomplete"
            else:
                mark = "—"
            print(f"  Status:             {mark}")

        print()
        print("  Energy:")
        if last["E0"] is not None:
            print(f"    Final E0        = {last['E0']:.8f} eV")
        if last["dE"] is not None:
            print(f"    ΔE (last step)  =  {last['dE']:.2e} eV")

    if force_steps:
        last_f = force_steps[-1]
        print()
        print("  Forces:")
        if last_f["max_force"] is not None:
            tick = "✓" if last_f["max_force"] < fthresh else "✗"
            print(f"    Max |F| (final) =  {last_f['max_force']:.4f} eV/Å")
            print(f"    Threshold       =  {fthresh:.4f} eV/Å    {tick}")
        if last_f["pressure"] is not None:
            print()
            print("  Pressure:")
            print(f"    Final P         =  {last_f['pressure']:.2f} kB")

    if wtime is not None:
        print()
        print("  Timing:")
        print(f"    Total walltime  =  {_fmt_time(wtime)}")
        if n_done > 0:
            print(f"    Avg per step    =  {_fmt_time(wtime / n_done)}")

    print()


# ── Matplotlib style ──────────────────────────────────────────────────────────

def _apply_style(font="Arial"):
    font_map = {
        "arial": "Arial", "helvetica": "Helvetica",
        "times": "Times New Roman", "times new roman": "Times New Roman",
    }
    mpl.rcParams["font.family"]          = font_map.get(font.lower(), "Arial")
    mpl.rcParams["axes.linewidth"]       = 1.4
    mpl.rcParams["font.weight"]          = "bold"
    mpl.rcParams["font.size"]            = 12
    mpl.rcParams["xtick.major.width"]    = 1.2
    mpl.rcParams["ytick.major.width"]    = 1.2
    mpl.rcParams["xtick.direction"]      = "in"
    mpl.rcParams["ytick.direction"]      = "in"
    mpl.rcParams["xtick.minor.visible"]  = True
    mpl.rcParams["ytick.minor.visible"]  = True
    mpl.rcParams["xtick.minor.width"]    = 0.8
    mpl.rcParams["ytick.minor.width"]    = 0.8


_PRIMARY = "#4b0082"
_REFLINE = "#888888"


def _style_ax(ax):
    for sp in ax.spines.values():
        sp.set_linewidth(1.4)
    ax.tick_params(which="both", direction="in", top=True, right=True)


def _thresh_label(ax, y, label, log=False):
    """Draw a dashed reference line and annotate at the right edge."""
    ax.axhline(y, color=_REFLINE, lw=1.0, ls="--", zorder=1)
    xlim = ax.get_xlim()
    ypos = y * 1.5 if (log and y > 0) else y
    ax.text(xlim[1], ypos, label, ha="right", va="bottom",
            fontsize=9, color=_REFLINE, fontstyle="italic")


# ── Ionic convergence plot ────────────────────────────────────────────────────

def plot_ionic(steps, force_steps=None, ethresh=1e-4, fthresh=0.02,
               show_mag=False, show_stress=False, start=1, end=None,
               output="valyte_converge.png", dpi=400, font="Arial"):

    plt.style.use("default")
    _apply_style(font)

    # Complete ionic steps only
    ionic = [s for s in steps if not s.get("_incomplete") and s["E0"] is not None]
    if not ionic:
        print("No complete ionic steps to plot.")
        return

    # Apply start/end window
    if end is not None:
        ionic = [s for s in ionic if start <= s["number"] <= end]
    else:
        ionic = [s for s in ionic if s["number"] >= start]

    if not ionic:
        print("No ionic steps in the requested range.")
        return

    xs  = [s["number"] for s in ionic]
    e0s = [s["E0"]     for s in ionic]
    des = [abs(s["dE"]) if s["dE"] is not None else None for s in ionic]

    has_mag = show_mag and any(s["mag"] is not None for s in ionic)

    # Force/pressure arrays aligned to the filtered ionic list.
    # force_steps[j] corresponds to the j-th COMPLETE ionic step overall,
    # so we need the original position in the unfiltered complete list.
    all_ionic = [s for s in steps if not s.get("_incomplete") and s["E0"] is not None]
    ionic_to_orig = {id(s): i for i, s in enumerate(all_ionic)}

    fmax_list, pres_list = [], []
    if force_steps is not None:
        for s in ionic:
            orig_i = ionic_to_orig.get(id(s))
            if orig_i is not None and orig_i < len(force_steps):
                fmax_list.append(force_steps[orig_i]["max_force"])
                pres_list.append(force_steps[orig_i]["pressure"])
            else:
                fmax_list.append(None)
                pres_list.append(None)

    has_forces = force_steps is not None and any(f is not None for f in fmax_list)
    has_pressure = show_stress and any(p is not None for p in pres_list)

    nrows = 2 + (1 if has_forces else 0) + (1 if has_pressure else 0) + (1 if has_mag else 0)

    fig, axes = plt.subplots(
        nrows, 1, sharex=True,
        figsize=(6, 2.5 * nrows),
        gridspec_kw={"hspace": 0.08},
    )
    axes = [axes] if nrows == 1 else list(axes)
    row = 0

    # Panel 1 — Energy
    ax = axes[row]; row += 1
    ax.plot(xs, e0s, color=_PRIMARY, lw=1.6, marker="o", ms=3.5, zorder=3)
    ax.set_ylabel("Energy (eV)", fontweight="bold")
    _style_ax(ax)

    # Panel 2 — |dE|
    ax = axes[row]; row += 1
    vx = [x for x, d in zip(xs, des) if d is not None and d > 0]
    vd = [d for d in des if d is not None and d > 0]
    if vx:
        ax.plot(vx, vd, color=_PRIMARY, lw=1.6, marker="o", ms=3.5, zorder=3)
    ax.set_yscale("log")
    ax.set_ylabel("|ΔE| (eV)", fontweight="bold")
    _style_ax(ax)
    if ethresh > 0:
        ax.axhline(ethresh, color=_REFLINE, lw=1.0, ls="--", zorder=1)
        xlim = ax.get_xlim()
        ax.text(xlim[1], ethresh * 1.5, f"{ethresh:.0e} eV",
                ha="right", va="bottom", fontsize=9,
                color=_REFLINE, fontstyle="italic")

    # Panel 3 (optional) — Max force
    if has_forces:
        ax = axes[row]; row += 1
        fx = [x for x, f in zip(xs, fmax_list) if f is not None]
        fy = [f for f in fmax_list if f is not None]
        if fy:
            ax.plot(fx, fy, color=_PRIMARY, lw=1.6, marker="o", ms=3.5, zorder=3)
        ax.set_yscale("log")
        ax.set_ylabel("Max |F| (eV/Å)", fontweight="bold")
        _style_ax(ax)
        if fthresh > 0:
            ax.axhline(fthresh, color=_REFLINE, lw=1.0, ls="--", zorder=1)
            xlim = ax.get_xlim()
            ax.text(xlim[1], fthresh * 1.5, f"{fthresh:.2f} eV/Å",
                    ha="right", va="bottom", fontsize=9,
                    color=_REFLINE, fontstyle="italic")

    # Panel 4 (optional) — Pressure
    if has_pressure:
        ax = axes[row]; row += 1
        px = [x for x, p in zip(xs, pres_list) if p is not None]
        py = [p for p in pres_list if p is not None]
        if py:
            ax.plot(px, py, color=_PRIMARY, lw=1.6, marker="o", ms=3.5, zorder=3)
        ax.axhline(0, color=_REFLINE, lw=1.0, ls="--", zorder=1)
        ax.set_ylabel("Pressure (kB)", fontweight="bold")
        _style_ax(ax)

    # Panel (optional) — Magnetization
    if has_mag:
        ax = axes[row]; row += 1
        mx = [s["number"] for s in ionic if s["mag"] is not None]
        my = [s["mag"]    for s in ionic if s["mag"] is not None]
        ax.plot(mx, my, color="#e63946", lw=1.6, marker="o", ms=3.5, zorder=3)
        ax.axhline(0, color=_REFLINE, lw=1.0, ls="--", zorder=1)
        ax.set_ylabel("Mag. (μ_B)", fontweight="bold")
        _style_ax(ax)

    axes[-1].set_xlabel("Ionic Step", fontweight="bold")
    axes[-1].set_xlim(min(xs) - 0.5, max(xs) + 0.5)

    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"Saved plot: {output}")


# ── Electronic convergence plot ───────────────────────────────────────────────

def plot_electronic(steps, ethresh=None, output="valyte_converge.png",
                    dpi=400, font="Arial"):

    plt.style.use("default")
    _apply_style(font)

    cumulative_x, cumulative_y = [], []
    boundaries = []
    cursor = 0

    for s in steps:
        elec = s.get("elec_steps", [])
        if not elec:
            continue
        start_x = cursor
        for e in elec:
            if e["dE"] > 0:
                cumulative_x.append(cursor)
                cumulative_y.append(e["dE"])
            cursor += 1
        if start_x > 0:
            boundaries.append((start_x - 0.5, str(s["number"])))

    if not cumulative_x:
        print("No electronic step data found.")
        return

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(cumulative_x, cumulative_y, color=_PRIMARY, lw=1.2,
            marker="o", ms=2.5, zorder=3)
    ax.set_yscale("log")

    ylim = ax.get_ylim()
    for xb, label in boundaries:
        ax.axvline(xb, color="#cccccc", lw=0.8, ls="--", zorder=1)
        ax.text(xb + 0.5, ylim[1] * 0.7, label, fontsize=7,
                color="#888888", va="top")

    if ethresh is not None and ethresh > 0:
        ax.axhline(ethresh, color=_REFLINE, lw=1.0, ls="--", zorder=1)
        xlim = ax.get_xlim()
        ax.text(xlim[1], ethresh * 1.5, f"{ethresh:.0e} eV",
                ha="right", va="bottom", fontsize=9,
                color=_REFLINE, fontstyle="italic")

    ax.set_xlabel("Electronic Step", fontweight="bold")
    ax.set_ylabel("|ΔE| (eV)", fontweight="bold")
    _style_ax(ax)

    plt.tight_layout()
    plt.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"Saved plot: {output}")


# ── Data export ───────────────────────────────────────────────────────────────

def save_converge_dat(steps, force_steps=None, filepath="valyte_converge.dat"):
    ionic = [s for s in steps if not s.get("_incomplete") and s["E0"] is not None]
    with open(filepath, "w") as f:
        f.write("# Step  E0(eV)  dE(eV)  F_max(eV/A)  P(kB)  mag\n")
        for i, s in enumerate(ionic):
            fmax = pres = "—"
            if force_steps and i < len(force_steps):
                fs = force_steps[i]
                if fs["max_force"] is not None:
                    fmax = f"{fs['max_force']:.6f}"
                if fs["pressure"] is not None:
                    pres = f"{fs['pressure']:.3f}"
            mag = f"{s['mag']:.4f}" if s["mag"] is not None else "—"
            dE  = f"{s['dE']:.6e}"  if s["dE"]  is not None else "—"
            e0  = f"{s['E0']:.8f}"  if s["E0"]  is not None else "—"
            f.write(f"  {s['number']:5d}  {e0}  {dE}  {fmax}  {pres}  {mag}\n")
    print(f"Saved data: {filepath}")


# ── Main entry point ──────────────────────────────────────────────────────────

def run_converge(path=".", electronic=False, forces=False, stress=False,
                 ethresh=1e-4, fthresh=0.02, start=1, end=None,
                 output="valyte_converge.png", save_data=False,
                 no_plot=False, mag=False):

    # Resolve paths
    if os.path.isfile(path):
        directory = os.path.dirname(os.path.abspath(path))
        oszicar_path = path
    else:
        directory = os.path.abspath(path)
        oszicar_path = _find_file(directory, "OSZICAR")

    if oszicar_path is None:
        print(f"Error: no OSZICAR found in {directory}")
        sys.exit(1)

    outcar_path = _find_file(directory, "OUTCAR")
    incar_path  = _find_file(directory, "INCAR")

    steps = parse_oszicar(oszicar_path)
    if not steps:
        print("No data found in OSZICAR.")
        sys.exit(1)

    outcar_info = {"ibrion": None, "nsw": None, "ediff": None,
                   "ediffg": None, "nions": None, "walltime_s": None}

    if outcar_path:
        try:
            outcar_info = parse_outcar_header(outcar_path)
        except Exception:
            pass

    if incar_path:
        try:
            incar_info = _parse_incar(incar_path)
            for k in ("ibrion", "nsw", "ediff", "ediffg"):
                if outcar_info.get(k) is None:
                    outcar_info[k] = incar_info.get(k)
        except Exception:
            pass

    # Use EDIFF from job files as ethresh default
    if outcar_info.get("ediff") is not None and ethresh == 1e-4:
        ethresh = outcar_info["ediff"]

    # Parse forces/stress only when requested
    force_steps = None
    if (forces or stress) and outcar_path:
        nions = outcar_info.get("nions")
        if nions:
            try:
                force_steps = parse_outcar_forces(outcar_path, nions)
            except Exception as e:
                print(f"Warning: could not parse forces from OUTCAR: {e}")
        else:
            print("Warning: NIONS not found — skipping force parsing.")

    print_summary(steps, outcar_info,
                  force_steps=force_steps if forces else None,
                  fthresh=fthresh)

    if no_plot:
        if save_data:
            save_converge_dat(steps, force_steps, "valyte_converge.dat")
        return

    if electronic:
        ediff = outcar_info.get("ediff") or ethresh
        plot_electronic(steps, ethresh=ediff, output=output)
    else:
        plot_ionic(
            steps,
            force_steps=force_steps if (forces or stress) else None,
            ethresh=ethresh,
            fthresh=fthresh,
            show_mag=mag,
            show_stress=stress,
            start=start,
            end=end,
            output=output,
        )

    if save_data:
        save_converge_dat(steps, force_steps, "valyte_converge.dat")
