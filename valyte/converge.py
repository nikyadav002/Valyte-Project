"""VASP convergence monitoring — OSZICAR/OUTCAR parsing."""

import gzip
import os
import re
import sys


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

_RE_IONIC = re.compile(
    r"^\s*(\d+)\s+F=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r"\s+E0=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
    r"\s+d\s*E\s*=\s*([-+]?\d*\.?\d+(?:[eEdD][+-]?\d+)?)"
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

    Only completed ionic steps are returned; a trailing SCF cycle with no
    ionic summary yet is ignored.
    """
    steps = []

    with _open(path) as fh:
        for raw in fh:
            line = raw.rstrip()

            m = _RE_IONIC.match(line)
            if m:
                step = {
                    "number": int(m.group(1)),
                    "F":  _float(m.group(2)),
                    "E0": _float(m.group(3)),
                    "dE": _float(m.group(4)),
                    "mag": None,
                }
                mm = _RE_MAG.search(line)
                if mm:
                    step["mag"] = _float(mm.group(1))
                steps.append(step)
                continue

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


def parse_outcar_forces(path):
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
        return "Relaxation"
    if nsw:
        return f"Relaxation (IBRION = {ibrion}, NSW = {nsw})"
    return f"Relaxation (IBRION = {ibrion})"


# ── Terminal summary ──────────────────────────────────────────────────────────

def print_summary(steps, outcar_info, force_steps=None):
    ionic = steps
    ibrion = outcar_info.get("ibrion")
    nsw    = outcar_info.get("nsw")
    ediffg = outcar_info.get("ediffg")
    wtime  = outcar_info.get("walltime_s")

    # EDIFFG sets the criterion: negative is a force threshold in eV/Å,
    # positive is an energy threshold in eV.
    fthresh = abs(ediffg) if ediffg is not None and ediffg < 0 else None
    ethresh = ediffg if ediffg is not None and ediffg > 0 else None

    print("Convergence")
    print("═" * 11)
    print()
    print(f"  Calculation type:   {_calc_type(ibrion, nsw)}")

    n_done = len(ionic)
    if nsw:
        print(f"  Ionic steps:        {n_done} / {nsw}")
    else:
        print(f"  Ionic steps:        {n_done}")

    if ionic:
        last = ionic[-1]

        if fthresh is not None:
            fok = (force_steps and force_steps[-1]["max_force"] is not None
                   and force_steps[-1]["max_force"] < fthresh)
            mark = "✓ Converged" if fok else "✗ Not converged"
        elif ethresh is not None:
            eok = last["dE"] is not None and abs(last["dE"]) < ethresh
            mark = "✓ Converged" if eok else "✗ Not converged"
        elif nsw and n_done < nsw:
            mark = "Running / incomplete"
        else:
            mark = "—"
        print(f"  Status:             {mark}")

        # ── Per-step table ────────────────────────────────────────────────
        # Step/E0/ΔE come from OSZICAR; Max |F| and P need OUTCAR and show
        # "—" when it is absent or does not carry the value.
        print()
        hdr = (f"  {'Step':>5s}   {'E0 (eV)':>16s}   {'ΔE (eV)':>14s}"
               f"   {'Max |F| (eV/Å)':>16s}   {'P (kB)':>10s}")
        print(hdr)
        print("  " + "─" * (len(hdr) - 2))

        for i, s in enumerate(ionic):
            e0_str = f"{s['E0']:.8f}" if s["E0"] is not None else "—"

            if s["dE"] is not None:
                de_str = f"{s['dE']:.2e}"
                if ethresh is not None and abs(s["dE"]) < ethresh:
                    de_str += " ✓"
            else:
                de_str = "—"

            fs = force_steps[i] if force_steps and i < len(force_steps) else None
            if fs and fs["max_force"] is not None:
                f_val = fs["max_force"]
                f_str = f"{f_val:.4f}"
                if fthresh is not None and f_val < fthresh:
                    f_str += " ✓"
            else:
                f_str = "—"

            p_str = f"{fs['pressure']:.2f}" if fs and fs["pressure"] is not None else "—"

            print(f"  {s['number']:5d}   {e0_str:>16s}   {de_str:>14s}"
                  f"   {f_str:>16s}   {p_str:>10s}")

        # ── Criterion footer ──────────────────────────────────────────────
        if fthresh is not None:
            print()
            print(f"  EDIFFG           =  {ediffg:g}  →  force criterion")
            print(f"  Force threshold  =  {fthresh:.4f} eV/Å")
            if force_steps and len(force_steps) >= len(ionic):
                last_f = force_steps[len(ionic) - 1]
                if last_f["max_force"] is not None:
                    tick = "✓" if last_f["max_force"] < fthresh else "✗"
                    print(f"  Final Max |F|    =  {last_f['max_force']:.4f} eV/Å    {tick}")
        elif ethresh is not None:
            print()
            print(f"  EDIFFG           =  {ediffg:g}  →  energy criterion")
            print(f"  Energy threshold =  {ethresh:.2e} eV")
            if last["dE"] is not None:
                tick = "✓" if abs(last["dE"]) < ethresh else "✗"
                print(f"  Final ΔE         =  {last['dE']:.2e} eV       {tick}")
        else:
            print()
            print("  EDIFFG not found in INCAR or OUTCAR — no convergence criterion.")

    if wtime is not None:
        print()
        print("  Timing:")
        print(f"    Total walltime  =  {_fmt_time(wtime)}")
        if n_done > 0:
            print(f"    Avg per step    =  {_fmt_time(wtime / n_done)}")

    print()


# ── Data export ───────────────────────────────────────────────────────────────

def save_converge_dat(steps, force_steps=None, filepath="valyte_converge.dat"):
    ionic = [s for s in steps if s["E0"] is not None]
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

def run_converge(path=".", save_data=False):

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
            # EDIFFG is read from INCAR when present; the OUTCAR echo is the
            # fallback for directories where only OUTCAR was kept.
            if incar_info.get("ediffg") is not None:
                outcar_info["ediffg"] = incar_info["ediffg"]
            for k in ("ibrion", "nsw", "ediff"):
                if outcar_info.get(k) is None:
                    outcar_info[k] = incar_info.get(k)
        except Exception:
            pass

    # Forces and pressure, when an OUTCAR is available
    force_steps = None
    if outcar_path:
        try:
            force_steps = parse_outcar_forces(outcar_path)
        except Exception as e:
            print(f"Warning: could not parse forces from OUTCAR: {e}")

    print_summary(steps, outcar_info, force_steps=force_steps)

    if save_data:
        save_converge_dat(steps, force_steps, "valyte_converge.dat")
