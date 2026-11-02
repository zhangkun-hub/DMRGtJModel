"""A3: extract the per-site truncation-error archive from DMRG stdout logs.

The dmrg binary prints one line per two-site update (master rank):

    Site (   l,   r) E0 = <16位> TruncErr = <sci> D = <5位> Iter = <3位>
        LanczT = <8位> TotT = <8位> S = <10位>

plus sweep markers `sweep <n>` and the final `E0/site: <v>` line.
This script turns one or several log files into a per-run JSON archive.

The truncation-error record is an MPS-quality archive. It is NOT a direct
error estimate of the correlation functions, and must not be merged into
any noise-floor analysis.

Usage:
    python extract_trunc_err.py LOG [LOG...] [-o trunc_err_record.json]
"""
import argparse
import json
import re

SITE_LINE_RE = re.compile(r"^\s*Site\s*\(?\s*(\d+)\s*,?\s*(\d*)\s*\)?(.*)$")
KV_RE = re.compile(r"([A-Za-z]+)\s*=\s*([^\s]+)")


def _to_number(s):
    try:
        return float(s)
    except ValueError:
        return s  # keep as string, e.g. "-nan"


def parse_site_line(line):
    m = SITE_LINE_RE.match(line)
    if not m:
        return None
    site_l = int(m.group(1))
    site_r = int(m.group(2)) if m.group(2) else site_l
    kv = {k: _to_number(v) for k, v in KV_RE.findall(m.group(3))}
    return {"site_l": site_l, "site_r": site_r, **kv}


def parse_log(path):
    sweeps = []
    cur = None
    e0_per_site = None
    d_max_sequence = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("sweep"):
                m = re.match(r"sweep\s+(\d+)", stripped)
                if m:
                    cur = {"sweep": int(m.group(1)), "sites": []}
                    sweeps.append(cur)
                continue
            if stripped.startswith("E0/site"):
                m = re.search(r"E0/site:\s*([^\s]+)", line)
                if m:
                    e0_per_site = _to_number(m.group(1))
                continue
            if stripped.startswith("D_max"):
                m = re.search(r"D_max\s*=\s*([^\s]+)", line)
                if m:
                    d_max_sequence.append(int(float(m.group(1))))
                continue
            site = parse_site_line(line)
            if site is not None:
                if cur is None:
                    cur = {"sweep": None, "sites": []}
                    sweeps.append(cur)
                cur["sites"].append(site)
    return {
        "E0_per_site": e0_per_site,
        "D_max_sequence": d_max_sequence,
        "sweeps": sweeps,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+", help="DMRG stdout/slurm log files")
    ap.add_argument("-o", "--out", default="trunc_err_record.json")
    args = ap.parse_args()

    record = {"records": []}
    for path in args.logs:
        rec = {"log_file": path}
        rec.update(parse_log(path))
        record["records"].append(rec)
        n_sites = sum(len(s["sites"]) for s in rec["sweeps"])
        print(f"{path}: sweeps={len(rec['sweeps'])} site_updates={n_sites} "
              f"E0/site={rec['E0_per_site']}")
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=1)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
