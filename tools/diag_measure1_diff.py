"""Compare two measurement dirs row by row and print the largest deviations.

Usage:
  python3 diag_measure1_diff.py <real_dir> <complex_dir> [top_n]

Diagnostic for the N=64 smoke verify v2 (measure1 generic mode): prints the
rows with the largest scaled |Re(complex)-real| deviation and per-file maxima,
so a FAIL can be classified as (a) uniform convergence noise across files or
(b) systematic deviation concentrated in one observable / a few rows.

Scaled metric identical to tools/verify_complex_reproduces_real.py:
  scaled = |Re(complex) - real| / (atol + rtol*|real|), atol = rtol = 1e-8.

Python 3.6 compatible; read-only.
"""
import ast
import os
import sys
from collections import defaultdict


def load_rows(path):
    rows = {}
    for e in ast.literal_eval(open(path, encoding="utf-8").read()):
        v = e[1]
        val = complex(v[0], v[1]) if isinstance(v, (list, tuple)) \
            else complex(float(v), 0.0)
        rows[tuple(int(x) for x in e[0])] = val
    return rows


def main():
    if len(sys.argv) < 3:
        print("usage: diag_measure1_diff.py <real_dir> <complex_dir> [top_n]")
        sys.exit(2)
    real_dir, cplx_dir = sys.argv[1], sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    atol = rtol = 1e-8
    diffs = []
    for fn in sorted(os.listdir(real_dir)):
        if not fn.endswith(".json"):
            continue
        cpath = os.path.join(cplx_dir, fn)
        if not os.path.isfile(cpath):
            continue
        rr = load_rows(os.path.join(real_dir, fn))
        cc = load_rows(cpath)
        for sites in rr:
            if sites not in cc:
                continue
            rv, cv = rr[sites], cc[sites]
            d = abs(cv.real - rv.real)
            thr = atol + rtol * abs(rv)
            diffs.append((d / thr if thr > 0 else d, fn, list(sites), rv, cv))
    diffs.sort(key=lambda t: -t[0])
    print("top %d scaled diffs (scaled = |re-diff|/(atol+rtol*|real|)):"
          % top_n)
    for sc, fn, sites, rv, cv in diffs[:top_n]:
        print("%12.1f %-34s %-14s real=%.10g cplx=%.10g"
              % (sc, fn, sites, rv.real, cv.real))
    agg = defaultdict(lambda: [0.0, 0])
    for sc, fn, _s, _rv, _cv in diffs:
        a = agg[fn]
        a[0] = max(a[0], sc)
        a[1] += 1
    print("\nper-file max scaled diff / n rows:")
    for fn, (m, n) in sorted(agg.items(), key=lambda kv: -kv[1][0]):
        print("%12.1f %5d  %s" % (m, n, fn))


if __name__ == "__main__":
    main()
