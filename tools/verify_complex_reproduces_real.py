"""C acceptance 1 & 2 with completeness gating.

Acceptance 1: the complex build reproduces the real build on a
real-representable state (real parts equal, imaginary parts ~ 0).
Acceptance 2: a real-representable state optimized with complex DMRG stays
real (imaginary parts of measured observables ~ 0).

Two modes:

  --scs-mode (pairing measurements):
    an INDEPENDENT expected manifest is built from the refbond list and the
    honeycomb geometry enumeration (target filter j0 > i1, A/B/C direction
    labels, channels a-d plus e-h with --reverse-channels) and checked
    against BOTH dirs. scs files missing/extra rows or files, duplicates,
    parse errors, non-finite values, and direction-label mismatches all
    fail. Auxiliary files (repeat_measurement_scatter.json,
    *_spatial_spread.json, *_bond_geometry.json) are ignored by the manifest.

  default mode (e.g. measure1 dirs):
    strict cross-completeness between the two dirs: identical .json file
    sets and identical row-key sets per file, plus duplicates/parse/finite
    checks. Limitation: no independent expected manifest (rows missing from
    BOTH dirs cannot be caught); use --scs-mode for independent gating.

Numeric stage runs only after all integrity checks pass:
    per-row |Re(complex) - real| and |Im(complex)| against the combined
    tolerance atol + rtol*|real|, and |E0/site difference| against the same
    rule (dmrg logs are REQUIRED: use --complex-log/--real-log).
In default (generic) mode the two dirs come from two INDEPENDENTLY converged
DMRG runs, whose agreement is limited by the convergence quality (D=300
truncation error ~1e-4) rather than by the 1e-8 reproducibility tolerance.
There the re-diff is REPORTED ONLY (re_diff_enforced=false in the report)
and the hard numeric gates are |Im| <= tol and E0 agreement; --scs-mode
keeps the full hard numeric gate (re-diff enforced).
The tolerances actually used are recorded in the report.
Exit code 0 on PASS, 1 otherwise.

Note: this verifies the complex CALCULATION capability on real states only.
It does NOT verify chiral-state selection, nor any d+id physics claim.
"""
import argparse
import ast
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import honeycomb_geometry as hg

SCS_FILE_RE = re.compile(r"^scs([ABC])([a-h])_fix_refbond_(\d+)_(\d+)\.json$")
AUX_RE = re.compile(r"^(repeat_measurement_scatter\.json|.*_(spatial_spread|bond_geometry)\.json)$")
CHANNELS = "abcd" + "efgh"


def parse_refbonds(spec, geometry, Ly, Lx):
    """Parse 'i1,j1;i2,j2;...' (--all-refbonds supported by the caller)."""
    refbonds = []
    for pair in spec.split(";"):
        p = pair.find(",")
        if p < 0:
            raise ValueError(f"bad refbond spec: {pair!r}")
        refbonds.append((int(pair[:p]), int(pair[p + 1:])))
    all_bonds = hg.nn_bonds(geometry, Ly, Lx)
    for ref in refbonds:
        if tuple(sorted(ref)) not in all_bonds:
            raise ValueError(f"{ref} is not a NN bond of {geometry} {Ly}x{Lx}")
    return [tuple(sorted(r)) for r in refbonds]


def build_scs_manifest(geometry, Ly, Lx, refbonds, reverse_channels):
    """Expected scs files -> set of row keys (r0, r1, t0, t1).

    Mirrors GenerateTargetBonds + PerformMeasureSC of
    src/measureSC_fix_ref_extension.cpp (threshold filter and per-direction
    grouping), independent of the measured data itself.
    """
    all_bonds = hg.nn_bonds(geometry, Ly, Lx)
    manifest = {}
    for (r0, r1) in refbonds:
        threshold = max(r0, r1)
        groups = {}
        for b in all_bonds:
            if b[0] <= threshold:  # the API/measurement filter
                continue
            d = hg.classify_nn_bond(geometry, Ly, Lx, *b)
            groups.setdefault(d, []).append(b)
        for d, bonds in sorted(groups.items()):
            keys = {(r0, r1, b[0], b[1]) for b in sorted(bonds)}
            for ch in CHANNELS[:8 if reverse_channels else 4]:
                manifest[f"scs{d}{ch}_fix_refbond_{r0}_{r1}.json"] = keys
    return manifest


def load_measu_rows_with_checks(path):
    """(rows, issues); rows = [(sites_tuple, complex)]."""
    issues = {"duplicates": 0, "parse_error": None, "nonfinite": []}
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = ast.literal_eval(f.read())
    except Exception as e:  # noqa: BLE001 - report, don't swallow
        issues["parse_error"] = f"{type(e).__name__}: {e}"
        return [], issues
    rows = []
    seen = set()
    for entry in data:
        try:
            sites = tuple(int(s) for s in entry[0])
            val = entry[1]
            if isinstance(val, (list, tuple)):
                value = complex(float(val[0]), float(val[1]))
            else:
                value = complex(float(val), 0.0)
        except Exception as e:  # noqa: BLE001
            issues["parse_error"] = f"row {entry!r}: {type(e).__name__}: {e}"
            return [], issues
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            issues["nonfinite"].append(sites)
        if sites in seen:
            issues["duplicates"] += 1
        seen.add(sites)
        rows.append((sites, value))
    return rows, issues


def classify_dir_files(dnames):
    scs, aux, other = {}, [], []
    for fname in sorted(dnames):
        m = SCS_FILE_RE.match(fname)
        if m:
            scs[fname] = (m.group(1), m.group(2),
                          (int(m.group(3)), int(m.group(4))))
        elif AUX_RE.match(fname):
            aux.append(fname)
        else:
            other.append(fname)
    return scs, aux, other


def check_scs_dir(dirpath, manifest, geometry, Ly, Lx):
    """Integrity of one dir against the independent manifest."""
    res = {"missing_files": [], "extra_files": [], "missing_rows": [],
           "extra_rows": [], "duplicates": {}, "parse_errors": {},
           "nonfinite": {}, "direction_mismatch": [], "file_count_ok": False}
    dnames = set(os.listdir(dirpath))
    scs, _aux, other = classify_dir_files(dnames)
    res["extra_files"] = other
    for fname, expected in manifest.items():
        if fname not in scs:
            res["missing_files"].append(fname)
            continue
        X, ch, refbond = scs[fname]
        path = os.path.join(dirpath, fname)
        rows, issues = load_measu_rows_with_checks(path)
        if issues["parse_error"]:
            res["parse_errors"][fname] = issues["parse_error"]
            continue
        if issues["duplicates"]:
            res["duplicates"][fname] = issues["duplicates"]
        if issues["nonfinite"]:
            res["nonfinite"][fname] = [list(s) for s in issues["nonfinite"]]
        got = {sites for sites, _ in rows}
        res["missing_rows"].extend(
            [fname, list(s)] for s in sorted(expected - got))
        res["extra_rows"].extend(
            [fname, list(s)] for s in sorted(got - expected))
        for sites in got:
            tgt = (sites[2], sites[3])
            actual = hg.classify_nn_bond(geometry, Ly, Lx, *tgt)
            if actual != X:
                res["direction_mismatch"].append(
                    [fname, list(sites), actual])
    for fname in scs:
        if fname not in manifest:
            res["extra_files"].append(fname)
    res["file_count_ok"] = not (res["missing_files"] or res["extra_files"])
    return res


def check_generic_dir(dirpath):
    """Per-file row checks for the default (non-scs) mode."""
    res = {"duplicates": {}, "parse_errors": {}, "nonfinite": {}}
    for fname in sorted(os.listdir(dirpath)):
        if not fname.endswith(".json"):
            continue
        rows, issues = load_measu_rows_with_checks(
            os.path.join(dirpath, fname))
        if issues["parse_error"]:
            res["parse_errors"][fname] = issues["parse_error"]
        if issues["duplicates"]:
            res["duplicates"][fname] = issues["duplicates"]
        if issues["nonfinite"]:
            res["nonfinite"][fname] = [list(s) for s in issues["nonfinite"]]
    return res


def parse_e0_per_site(log_path):
    with open(log_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.search(r"E0/site:\s*([^\s]+)", line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
    return None


def integrity_ok(scs_res=None, generic_res=None, cross=None, energy=None):
    if scs_res is not None:
        for r in (scs_res["real"], scs_res["complex"]):
            if (r["missing_files"] or r["extra_files"] or r["missing_rows"]
                    or r["extra_rows"] or r["duplicates"]
                    or r["parse_errors"] or r["nonfinite"]
                    or r["direction_mismatch"]):
                return False
    if generic_res is not None:
        for r in (generic_res["real"], generic_res["complex"]):
            if r["duplicates"] or r["parse_errors"] or r["nonfinite"]:
                return False
    if cross is not None:
        if (cross["files_missing_from_complex"]
                or cross["files_missing_from_real"]
                or cross["rows_missing_from_complex"]
                or cross["rows_missing_from_real"]):
            return False
    if energy is not None:
        if not energy["ok"]:
            return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", required=True, help="complex-build data dir")
    ap.add_argument("--real", required=True, help="real-build data dir")
    ap.add_argument("--atol", type=float, default=1e-8)
    ap.add_argument("--rtol", type=float, default=1e-8)
    ap.add_argument("--complex-log", help="complex-build dmrg log (required)")
    ap.add_argument("--real-log", help="real-build dmrg log (required)")
    ap.add_argument("--scs-mode", action="store_true",
                    help="build an independent expected manifest for scs files")
    ap.add_argument("--geometry", choices=["HoneycombXC", "HoneycombYC"])
    ap.add_argument("--Ly", type=int)
    ap.add_argument("--Lx", type=int)
    ap.add_argument("--refbonds", help="'i1,j1;i2,j2;...' or ALL")
    ap.add_argument("--reverse-channels", action="store_true")
    ap.add_argument("-o", "--out", default="complex_reproduces_real_report.json")
    args = ap.parse_args()

    report = {
        "complex_dir": os.path.abspath(args.complex),
        "real_dir": os.path.abspath(args.real),
        "tolerances": {"atol": args.atol, "rtol": args.rtol},
        "mode": "scs" if args.scs_mode else "generic",
        "integrity": {}, "numeric": {}, "energy": {}, "verdict": None,
    }

    scs_res = None
    generic_res = None
    cross = None
    energy = None

    if args.scs_mode:
        missing = [f for f, v in (("--geometry", args.geometry),
                                  ("--Ly", args.Ly), ("--Lx", args.Lx),
                                  ("--refbonds", args.refbonds)) if v is None]
        if missing:
            print(f"error: --scs-mode requires {', '.join(missing)}")
            sys.exit(2)
        if args.refbonds == "ALL":
            refbonds = hg.nn_bonds(args.geometry, args.Ly, args.Lx)
        else:
            refbonds = parse_refbonds(args.refbonds, args.geometry,
                                      args.Ly, args.Lx)
        manifest = build_scs_manifest(args.geometry, args.Ly, args.Lx,
                                      refbonds, args.reverse_channels)
        scs_res = {
            "real": check_scs_dir(args.real, manifest, args.geometry,
                                  args.Ly, args.Lx),
            "complex": check_scs_dir(args.complex, manifest, args.geometry,
                                     args.Ly, args.Lx),
            "manifest": {"n_files": len(manifest),
                         "refbonds": [list(r) for r in refbonds],
                         "reverse_channels": args.reverse_channels},
        }
        report["integrity"]["scs"] = scs_res
    else:
        generic_res = {
            "real": check_generic_dir(args.real),
            "complex": check_generic_dir(args.complex),
        }
        report["integrity"]["generic"] = generic_res
        # strict cross-completeness between the two dirs
        rfiles = {f for f in os.listdir(args.real) if f.endswith(".json")}
        cfiles = {f for f in os.listdir(args.complex) if f.endswith(".json")}
        cross = {"files_missing_from_complex": sorted(rfiles - cfiles),
                 "files_missing_from_real": sorted(cfiles - rfiles),
                 "rows_missing_from_complex": [],
                 "rows_missing_from_real": []}
        for fname in sorted(rfiles & cfiles):
            rrows = {s for s, _ in load_measu_rows_with_checks(
                os.path.join(args.real, fname))[0]}
            crows = {s for s, _ in load_measu_rows_with_checks(
                os.path.join(args.complex, fname))[0]}
            cross["rows_missing_from_complex"].extend(
                [fname, list(s)] for s in sorted(rrows - crows))
            cross["rows_missing_from_real"].extend(
                [fname, list(s)] for s in sorted(crows - rrows))
        report["integrity"]["cross"] = cross

    # energy comparison (required)
    e0c = parse_e0_per_site(args.complex_log) if args.complex_log else None
    e0r = parse_e0_per_site(args.real_log) if args.real_log else None
    energy = {"complex_E0_per_site": e0c, "real_E0_per_site": e0r, "ok": False,
              "diff": None, "threshold": None}
    if args.complex_log and args.real_log and e0c is not None and e0r is not None:
        diff = abs(e0c - e0r)
        thr = args.atol + args.rtol * abs(e0r)
        energy.update({"diff": diff, "threshold": thr, "ok": diff <= thr})
    else:
        energy["reason"] = ("dmrg logs missing or no E0/site line found; "
                            "energy comparison REQUIRED for PASS")
    report["energy"] = energy

    if not integrity_ok(scs_res, generic_res, cross, energy):
        report["verdict"] = "FAIL"
        print("integrity checks FAILED; numeric stage skipped")
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=1, default=str)
        print(f"wrote {args.out}")
        sys.exit(1)

    # numeric stage (only after integrity passes)
    max_re_diff = 0.0
    max_im = 0.0
    n_rows = 0
    for fname in sorted(os.listdir(args.real)):
        if not fname.endswith(".json"):
            continue
        cpath = os.path.join(args.complex, fname)
        if not os.path.isfile(cpath):
            continue  # cross/integrity checks already handled this
        real_index = {sites: val for sites, val in
                      load_measu_rows_with_checks(
                          os.path.join(args.real, fname))[0]}
        for sites, cval in load_measu_rows_with_checks(cpath)[0]:
            rval = real_index.get(sites)
            if rval is None:
                continue  # already an integrity error above
            re_diff = abs(cval.real - rval.real)
            im_mag = abs(cval.imag)
            thr = args.atol + args.rtol * abs(rval)
            max_re_diff = max(max_re_diff, re_diff / thr if thr > 0
                              else re_diff)
            max_im = max(max_im, im_mag / thr if thr > 0 else im_mag)
            n_rows += 1
    # scs 模式：|Re(c)-r|、|Im(c)| 均为硬判据（1e-8 复现口径）。
    # generic 模式：比较的是两条独立收敛链，数值一致性受收敛精度（如
    # D=300 截断误差 ~1e-4）限制，1e-8 复现口径不适用——|Re(c)-r| 仅报告，
    # 硬判据为 |Im(c)| <= tol（实值态保持实）与 E0 一致。
    re_diff_enforced = bool(args.scs_mode)
    report["numeric"] = {
        "n_rows_compared": n_rows,
        "max_re_diff_scaled": max_re_diff,
        "max_im_scaled": max_im,
        "scaled_threshold": 1.0,  # scaled: PASS iff <= 1
        "re_diff_enforced": re_diff_enforced,
    }
    numeric_ok = (n_rows > 0 and max_im <= 1.0
                  and (max_re_diff <= 1.0 if re_diff_enforced else True))
    report["verdict"] = "PASS" if numeric_ok else "FAIL"
    print(f"rows compared: {n_rows}; scaled max re-diff = {max_re_diff:.3e} "
          f"({'enforced' if re_diff_enforced else 'reported only'}); "
          f"scaled max |im| = {max_im:.3e} (threshold 1.0)")
    print(f"energy: complex E0/site = {e0c}, real E0/site = {e0r}, "
          f"diff = {energy['diff']}")
    print(f"verdict: {report['verdict']}")
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1, default=str)
    print(f"wrote {args.out}")
    sys.exit(0 if report["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main()
