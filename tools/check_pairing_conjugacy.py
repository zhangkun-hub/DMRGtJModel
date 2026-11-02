"""B: numerically verify the pairing-matrix conjugation relation.

Channels a-d measure P(a->b) = <Delta^dag_a Delta_b> for reference bond a
and target bond b (creation at the reference pair, annihilation at the
target pair). Channels e-h (--reverse-channels) measure the reversed roles
Q(a->b) = <Delta_a Delta^dag_b>. The conjugation relation

    P_{ba}(b, a) = conj(P_{ab}(a, b))

is then verified numerically through

    P(a->b) == conj(Q(a->b))

for every measured pair (both sides come from the same run of reference
bond a, so no extra measurement is needed; the API only measures target
pairs beyond the reference pair, and the "lower triangle" of the pairing
matrix follows analytically from this relation).

The script also cross-checks that every file label A/B/C matches the
honeycomb bond direction of its target bonds (ClassifyNearestNeighbor
mirrored in honeycomb_geometry.py).

Usage:
    python check_pairing_conjugacy.py --dir <measure_dir> \
        --geometry HoneycombXC --Ly 4 --Lx 8 \
        [--rtol 1e-6] [-o pairing_conjugacy_check.json]
"""
import argparse
import json
import os
import re

import honeycomb_geometry as hg
from measure_json_util import load_measu_rows

FILE_RE = re.compile(r"^scs([ABC])([a-h])_fix_refbond_(\d+)_(\d+)\.json$")
CHANNEL_PAIR = {"a": "e", "b": "f", "c": "g", "d": "h"}


def collect_files(measure_dir):
    """Return {(refbond, X, ch): {target_bond: value}}."""
    data = {}
    for fname in sorted(os.listdir(measure_dir)):
        m = FILE_RE.match(fname)
        if not m:
            continue
        X, ch = m.group(1), m.group(2)
        refbond = (int(m.group(3)), int(m.group(4)))
        rows = {}
        for sites, value in load_measu_rows(os.path.join(measure_dir, fname)):
            rows[(sites[2], sites[3])] = value  # target bond
        data[(refbond, X, ch)] = rows
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="measurement output dir")
    ap.add_argument("--geometry", required=True,
                    choices=["HoneycombXC", "HoneycombYC"])
    ap.add_argument("--Ly", type=int, required=True)
    ap.add_argument("--Lx", type=int, required=True)
    ap.add_argument("--rtol", type=float, default=1e-6,
                    help="relative-residual threshold for the verdict")
    ap.add_argument("-o", "--out", default="pairing_conjugacy_check.json")
    args = ap.parse_args()

    data = collect_files(args.dir)
    if not data:
        print(f"no scs*_fix_refbond_*.json files found in {args.dir}")
        return

    record = {
        "measure_dir": os.path.abspath(args.dir),
        "geometry": args.geometry, "Ly": args.Ly, "Lx": args.Lx,
        "rtol": args.rtol,
        "direction_label_check": {"errors": []},
        "conjugation_check": {"n_pairs": 0, "n_missing": 0,
                              "max_abs_residual": 0.0,
                              "max_rel_residual": 0.0,
                              "worst": []},
    }

    # 1) file label vs. actual bond direction
    for (refbond, X, ch), rows in data.items():
        hg.classify_nn_bond(args.geometry, args.Ly, args.Lx, *refbond)
        for target in rows:
            tgt_dir = hg.classify_nn_bond(args.geometry, args.Ly, args.Lx,
                                          *target)
            if tgt_dir != X:
                record["direction_label_check"]["errors"].append(
                    {"refbond": list(refbond), "file_label": X,
                     "target": list(target), "actual_dir": tgt_dir})

    # 2) P(a->b) vs conj(Q(a->b))
    worst = []
    for (refbond, X, ch), p_rows in data.items():
        if ch not in CHANNEL_PAIR:
            continue  # e-h channels are the Q side
        q_rows = data.get((refbond, X, CHANNEL_PAIR[ch]))
        if q_rows is None:
            record["conjugation_check"]["n_missing"] += len(p_rows)
            continue
        for target, p_val in p_rows.items():
            record["conjugation_check"]["n_pairs"] += 1
            q_val = q_rows.get(target)
            if q_val is None:
                record["conjugation_check"]["n_missing"] += 1
                continue
            residual = abs(p_val - q_val.conjugate())
            scale = max(abs(p_val), abs(q_val), 1e-300)
            rel = residual / scale
            if residual > record["conjugation_check"]["max_abs_residual"]:
                record["conjugation_check"]["max_abs_residual"] = residual
            if rel > record["conjugation_check"]["max_rel_residual"]:
                record["conjugation_check"]["max_rel_residual"] = rel
            worst.append({"refbond": list(refbond), "dir": X, "ch": ch,
                          "target": list(target), "P": [p_val.real, p_val.imag],
                          "conjQ": [q_val.conjugate().real,
                                    q_val.conjugate().imag],
                          "abs_residual": residual, "rel_residual": rel})
    worst.sort(key=lambda w: -w["rel_residual"])
    record["conjugation_check"]["worst"] = worst[:20]

    n_pairs = record["conjugation_check"]["n_pairs"]
    n_missing = record["conjugation_check"]["n_missing"]
    max_rel = record["conjugation_check"]["max_rel_residual"]
    n_err = len(record["direction_label_check"]["errors"])
    verdict = (
        "PASS" if n_err == 0 and n_pairs > 0 and max_rel <= args.rtol
        else "FAIL" if n_err > 0 or n_pairs == 0
        else "NEEDS_ATTENTION"
    )
    record["verdict"] = verdict
    print(f"direction-label errors : {n_err}")
    print(f"conjugation pairs      : {n_pairs} checked, {n_missing} missing")
    print(f"max abs/rel residual   : "
          f"{record['conjugation_check']['max_abs_residual']:.3e} / "
          f"{max_rel:.3e}")
    print(f"verdict                : {verdict}")
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
