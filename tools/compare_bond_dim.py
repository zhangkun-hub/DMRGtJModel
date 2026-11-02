"""A4: per-bond-dimension convergence record of the same observable.

Protocol: run dmrg with a single --D=<d> per stage, measure after each stage
into a separate directory, then compare across stages.

    python compare_bond_dim.py --ref dir_D1 --others dir_D2 [dir_D3 ...] \
        [-o bond_dim_convergence_record.json]

For each measurement JSON present in the reference dir, rows are aligned by
their site index list and the per-row value differences against the other
dirs are recorded (plus the max |difference| per file).

This record measures the sensitivity of observables to the bond dimension.
It is NOT an independent statistical noise estimate.
"""
import argparse
import json
import os

from measure_json_util import load_measu_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True,
                    help="reference (usually smallest D) measurement dir")
    ap.add_argument("--others", nargs="+", required=True,
                    help="other D measurement dirs")
    ap.add_argument("-o", "--out", default="bond_dim_convergence_record.json")
    args = ap.parse_args()

    record = {
        "reference_dir": os.path.abspath(args.ref),
        "other_dirs": [os.path.abspath(d) for d in args.others],
        "files": [],
    }
    for fname in sorted(os.listdir(args.ref)):
        if not fname.endswith(".json"):
            continue
        ref_rows = load_measu_rows(os.path.join(args.ref, fname))
        file_rec = {
            "file": fname,
            "n_rows": len(ref_rows),
            "max_abs_diff": 0.0,
            "max_abs_diff_dir": None,
            "comparisons": [],
        }
        for d in args.others:
            other_path = os.path.join(d, fname)
            if not os.path.isfile(other_path):
                file_rec["comparisons"].append(
                    {"dir": os.path.abspath(d), "status": "missing"})
                continue
            other_index = {tuple(sites): val
                           for sites, val in load_measu_rows(other_path)}
            rows = []
            for sites, ref_val in ref_rows:
                key = tuple(sites)
                if key not in other_index:
                    rows.append({"sites": list(sites), "status": "missing"})
                    continue
                delta = other_index[key] - ref_val
                rows.append({"sites": list(sites),
                             "diff": [delta.real, delta.imag]})
                mag = abs(delta)
                if mag > file_rec["max_abs_diff"]:
                    file_rec["max_abs_diff"] = mag
                    file_rec["max_abs_diff_dir"] = os.path.abspath(d)
            file_rec["comparisons"].append(
                {"dir": os.path.abspath(d), "rows": rows})
        record["files"].append(file_rec)
        print(f"{fname}: max_abs_diff = {file_rec['max_abs_diff']:.3e}")
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
