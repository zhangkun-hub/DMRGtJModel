"""Negative/positive test suite for tools/verify_complex_reproduces_real.py.

Builds a clean scs-mode synthetic setup (geometry XC Ly=4 Lx=8, refbond
(0,1), real + complex dirs + dmrg logs) and asserts:

  * the clean case PASSes (exit 0);
  * every integrity violation FAILs (exit 1): deleted file, deleted row,
    row deleted from BOTH dirs (manifest still expects it), extra row,
    duplicated row, corrupted JSON, NaN value, stray extra file;
  * generic (non-scs) mode: clean PASS and cross-missing row FAIL.

Run: python tools/tests/test_verify_negatives.py
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(HERE, os.pardir)
SCRIPT = os.path.join(TOOLS, "verify_complex_reproduces_real.py")
sys.path.insert(0, TOOLS)
import honeycomb_geometry as hg  # noqa: E402

GEOM, LY, LX = "HoneycombXC", 4, 8
REFBONDS = "0,1"


def build_expected():
    bonds = hg.nn_bonds(GEOM, LY, LX)
    threshold = 1
    groups = {}
    for b in bonds:
        if b[0] <= threshold:
            continue
        d = hg.classify_nn_bond(GEOM, LY, LX, *b)
        groups.setdefault(d, []).append(b)
    return {d: sorted(bs) for d, bs in sorted(groups.items())}


def write_dirs(base, real_rows, complex_rows):
    """Write the scs files for refbond (0,1), channels a-d."""
    rdir = os.path.join(base, "real")
    cdir = os.path.join(base, "complex")
    os.makedirs(rdir, exist_ok=True)
    os.makedirs(cdir, exist_ok=True)
    for d, bonds in build_expected().items():
        for ch in "abcd":
            fname = f"scs{d}{ch}_fix_refbond_0_1.json"
            with open(os.path.join(rdir, fname), "w") as f:
                json.dump([[list(k), v.real] for k, v in
                           sorted(real_rows[(d, ch)].items())], f)
            with open(os.path.join(cdir, fname), "w") as f:
                json.dump([[list(k), [v.real, v.imag]] for k, v in
                           sorted(complex_rows[(d, ch)].items())], f)


def write_logs(base):
    with open(os.path.join(base, "real.log"), "w") as f:
        f.write("E0/site: -0.45\n")
    with open(os.path.join(base, "complex.log"), "w") as f:
        f.write("E0/site: -0.4500000001\n")


def make_clean(base):
    rows = {}
    val = 0.1
    for d, bonds in build_expected().items():
        for ch in "abcd":
            entries = {}
            for b in bonds:
                entries[(0, 1, b[0], b[1])] = complex(val, 0.0)
                val += 0.01
            rows[(d, ch)] = entries
    complex_rows = {k: {kk: v + complex(0.0, 1e-12) for kk, v in vv.items()}
                    for k, vv in rows.items()}
    write_dirs(base, rows, complex_rows)
    write_logs(base)


def run_verify(base, extra_args=None):
    cmd = [sys.executable, SCRIPT,
           "--complex", os.path.join(base, "complex"),
           "--real", os.path.join(base, "real"),
           "--complex-log", os.path.join(base, "complex.log"),
           "--real-log", os.path.join(base, "real.log"),
           "--scs-mode", "--geometry", GEOM, "--Ly", str(LY),
           "--Lx", str(LX), "--refbonds", REFBONDS,
           "-o", os.path.join(base, "report.json")]
    if extra_args:
        cmd += extra_args
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc


def read_verdict(base):
    with open(os.path.join(base, "report.json")) as f:
        return json.load(f)["verdict"]


def expect(base, want_pass, label, extra_args=None):
    proc = run_verify(base, extra_args)
    got_pass = (proc.returncode == 0
                and read_verdict(base) == "PASS")
    status = "PASS" if got_pass else "FAIL"
    ok = (status == "PASS") == want_pass
    print(f"  [{status}] {label}  (rc={proc.returncode})")
    if not ok:
        print(f"    expected {'PASS' if want_pass else 'FAIL'}, stdout: "
              f"{proc.stdout.strip()[:400]}")
    assert ok, label


def main():
    tmp = tempfile.mkdtemp(prefix="verify_negatives_")
    try:
        # clean positive case
        base = os.path.join(tmp, "clean")
        make_clean(base)
        expect(base, True, "clean scs-mode case")

        # 1) delete a complete JSON file (complex side)
        base = os.path.join(tmp, "del_file")
        make_clean(base)
        os.remove(os.path.join(base, "complex",
                               "scsAa_fix_refbond_0_1.json"))
        expect(base, False, "deleted file (complex)")

        # 2) delete one row (complex side)
        base = os.path.join(tmp, "del_row")
        make_clean(base)
        path = os.path.join(base, "complex", "scsAa_fix_refbond_0_1.json")
        with open(path) as f:
            rows = json.load(f)
        with open(path, "w") as f:
            json.dump(rows[:-1], f)
        expect(base, False, "deleted row (complex)")

        # 3) delete the same row from BOTH dirs (manifest still expects it)
        base = os.path.join(tmp, "del_row_both")
        make_clean(base)
        for side in ("real", "complex"):
            path = os.path.join(base, side, "scsAa_fix_refbond_0_1.json")
            with open(path) as f:
                rows = json.load(f)
            with open(path, "w") as f:
                json.dump(rows[:-1], f)
        expect(base, False, "row deleted from both dirs")

        # 4) extra row (complex side)
        base = os.path.join(tmp, "extra_row")
        make_clean(base)
        path = os.path.join(base, "complex", "scsAa_fix_refbond_0_1.json")
        with open(path) as f:
            rows = json.load(f)
        rows.append([[0, 1, 999, 998], [0.5, 0.0]])
        with open(path, "w") as f:
            json.dump(rows, f)
        expect(base, False, "extra row (complex)")

        # 5) duplicated row (complex side)
        base = os.path.join(tmp, "dup_row")
        make_clean(base)
        path = os.path.join(base, "complex", "scsAa_fix_refbond_0_1.json")
        with open(path) as f:
            rows = json.load(f)
        rows.append(rows[0])
        with open(path, "w") as f:
            json.dump(rows, f)
        expect(base, False, "duplicated row (complex)")

        # 6) corrupted JSON (complex side)
        base = os.path.join(tmp, "corrupt")
        make_clean(base)
        path = os.path.join(base, "complex", "scsAa_fix_refbond_0_1.json")
        with open(path, "w") as f:
            f.write("[[[0, 1, 2, 3], [0.5,")
        expect(base, False, "corrupted JSON (complex)")

        # 7) non-finite value (complex side): 1e999 -> inf, valid literal
        base = os.path.join(tmp, "nan")
        make_clean(base)
        path = os.path.join(base, "complex", "scsAa_fix_refbond_0_1.json")
        with open(path) as f:
            rows = json.load(f)
        rows[0][1][0] = 1e999
        with open(path, "w") as f:
            f.write(json.dumps(rows).replace("Infinity", "1e999"))
        expect(base, False, "non-finite value (complex)")

        # 8) stray extra file (complex side)
        base = os.path.join(tmp, "extra_file")
        make_clean(base)
        with open(os.path.join(base, "complex", "stray.json"), "w") as f:
            f.write("[[[0], 1.0]]")
        expect(base, False, "stray extra file (complex)")

        # 9) missing E0 line in the complex log
        base = os.path.join(tmp, "no_e0")
        make_clean(base)
        with open(os.path.join(base, "complex.log"), "w") as f:
            f.write("sweep 1\n")
        expect(base, False, "E0/site missing from complex log")

        # generic mode: clean PASS
        base = os.path.join(tmp, "generic_clean")
        os.makedirs(os.path.join(base, "real"))
        os.makedirs(os.path.join(base, "complex"))
        with open(os.path.join(base, "real", "nf.json"), "w") as f:
            json.dump([[[0], 1.25], [[1], -0.75]], f)
        with open(os.path.join(base, "complex", "nf.json"), "w") as f:
            json.dump([[[0], [1.25, 1e-12]], [[1], [-0.75, -2e-13]]], f)
        write_logs(base)
        proc = subprocess.run(
            [sys.executable, SCRIPT, "--complex",
             os.path.join(base, "complex"), "--real",
             os.path.join(base, "real"), "--complex-log",
             os.path.join(base, "complex.log"), "--real-log",
             os.path.join(base, "real.log"), "-o",
             os.path.join(base, "report.json")],
            capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout
        print("  [PASS] generic-mode clean case")

        # generic mode: row missing from complex only -> cross check FAIL
        base = os.path.join(tmp, "generic_missing")
        shutil.copytree(os.path.join(tmp, "generic_clean"),
                        os.path.join(tmp, "generic_missing"))
        with open(os.path.join(base, "complex", "nf.json"), "w") as f:
            json.dump([[[0], [1.25, 1e-12]]], f)
        proc = subprocess.run(
            [sys.executable, SCRIPT, "--complex",
             os.path.join(base, "complex"), "--real",
             os.path.join(base, "real"), "--complex-log",
             os.path.join(base, "complex.log"), "--real-log",
             os.path.join(base, "real.log"), "-o",
             os.path.join(base, "report.json")],
            capture_output=True, text=True)
        assert proc.returncode != 0, "generic missing row must FAIL"
        print("  [FAIL] generic-mode missing row (cross check)")

        # generic mode: large re-diff (independent-run convergence noise)
        # must PASS now (re-diff is reported, not enforced)
        base = os.path.join(tmp, "generic_re_diff_only")
        shutil.copytree(os.path.join(tmp, "generic_clean"),
                        os.path.join(tmp, "generic_re_diff_only"))
        with open(os.path.join(base, "complex", "nf.json"), "w") as f:
            json.dump([[[0], [1.2501, 1e-12]], [[1], [-0.7501, -2e-13]]], f)
        proc = subprocess.run(
            [sys.executable, SCRIPT, "--complex",
             os.path.join(base, "complex"), "--real",
             os.path.join(base, "real"), "--complex-log",
             os.path.join(base, "complex.log"), "--real-log",
             os.path.join(base, "real.log"), "-o",
             os.path.join(base, "report.json")],
            capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout
        with open(os.path.join(base, "report.json")) as f:
            rep = json.load(f)
        assert rep["numeric"]["re_diff_enforced"] is False
        print("  [PASS] generic-mode re-diff reported only "
              "(independent-run noise)")

        # generic mode: large imaginary part must still FAIL (hard gate)
        base = os.path.join(tmp, "generic_im_fail")
        shutil.copytree(os.path.join(tmp, "generic_clean"),
                        os.path.join(tmp, "generic_im_fail"))
        with open(os.path.join(base, "complex", "nf.json"), "w") as f:
            json.dump([[[0], [1.25, 0.5]], [[1], [-0.75, -0.3]]], f)
        proc = subprocess.run(
            [sys.executable, SCRIPT, "--complex",
             os.path.join(base, "complex"), "--real",
             os.path.join(base, "real"), "--complex-log",
             os.path.join(base, "complex.log"), "--real-log",
             os.path.join(base, "real.log"), "-o",
             os.path.join(base, "report.json")],
            capture_output=True, text=True)
        assert proc.returncode != 0, "generic large |Im| must FAIL"
        print("  [FAIL] generic-mode large |Im| (hard gate)")

        # scs mode: re-diff is still enforced -> large re-diff must FAIL
        base = os.path.join(tmp, "scs_re_diff_fail")
        make_clean(base)
        path = os.path.join(base, "complex", "scsAa_fix_refbond_0_1.json")
        with open(path) as f:
            rows = json.load(f)
        rows[0][1][0] += 0.01
        with open(path, "w") as f:
            json.dump(rows, f)
        expect(base, False, "scs-mode large re-diff (hard gate)")

        print("test_verify_negatives: ALL PASS")
    finally:
        shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
