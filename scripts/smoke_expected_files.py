"""H2 helper: expected measurement output file set for a scs stage.

Mirrors GenerateTargetBonds + PerformMeasureSC of
src/measureSC_fix_ref_extension.cpp (threshold filter t0 <= r1, A/B/C
direction grouping, channels a-d plus e-h under --reverse-channels),
independently of the measured data, using tools/honeycomb_geometry.py
(which is verified against the C++ link construction).

Usage:
    python smoke_expected_files.py list <geometry> <Ly> <Lx> \
        --refbond=<i,j>[;<k,l>...] [--reverse-channels]
      prints the expected file names (one per line)

    python smoke_expected_files.py check <dir> <geometry> <Ly> <Lx> \
        --refbond=<i,j>[;<k,l>...] [--reverse-channels]
      asserts <dir> contains exactly the expected scs files plus
      repeat_measurement_scatter.json (no missing, no extra); prints the
      expected file names and exits non-zero on any violation.
"""
import os
import sys

# 显式 UTF-8 输出：避免 Windows GBK 本地环境下 stderr 编码与接收侧不一致
sys.stdout.reconfigure(encoding="utf-8", newline="\n")
sys.stderr.reconfigure(encoding="utf-8")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "tools"))
import honeycomb_geometry as hg  # noqa: E402

CHANNELS = "abcd" + "efgh"


def parse_args(argv):
    """list <geometry> <Ly> <Lx> [args...] | check <dir> <geometry> <Ly> <Lx> [args...]"""
    mode = argv[1]
    if mode == "check":
        check_dir = argv[2]
        geometry, ly, lx = argv[3], int(argv[4]), int(argv[5])
        rest = argv[6:]
    else:
        check_dir = None
        geometry, ly, lx = argv[2], int(argv[3]), int(argv[4])
        rest = argv[5:]
    refbonds = []
    reverse = False
    for arg in rest:
        if arg == "--reverse-channels":
            reverse = True
        elif arg.startswith("--refbond="):
            for pair in arg[len("--refbond="):].split(";"):
                i, j = pair.split(",")
                refbonds.append((int(i), int(j)))
        else:
            raise SystemExit(f"unexpected argument: {arg}")
    return mode, check_dir, geometry, ly, lx, refbonds, reverse


def expected_files(geometry, ly, lx, refbonds, reverse):
    all_bonds = hg.nn_bonds(geometry, ly, lx)
    files = []
    for (r0, r1) in refbonds:
        if tuple(sorted((r0, r1))) not in all_bonds:
            raise SystemExit(f"{(r0, r1)} 不是 {geometry} {ly}x{lx} 的 NN 键")
        threshold = max(r0, r1)
        groups = {}
        for b in all_bonds:
            if b[0] <= threshold:  # C++ 测量过滤（API 约束）
                continue
            d = hg.classify_nn_bond(geometry, ly, lx, *b)
            groups.setdefault(d, []).append(b)
        for d in sorted(groups):
            for ch in CHANNELS[:8 if reverse else 4]:
                files.append(f"scs{d}{ch}_fix_refbond_{r0}_{r1}.json")
    return files


def main():
    mode, check_dir, geometry, ly, lx, refbonds, reverse = parse_args(sys.argv)
    files = expected_files(geometry, ly, lx, refbonds, reverse)
    if mode == "list":
        print("\n".join(files))
        return
    if mode == "check":
        present = sorted(f for f in os.listdir(check_dir) if f.endswith(".json"))
        expected = sorted(files + ["repeat_measurement_scatter.json"])
        missing = [f for f in expected if f not in present]
        extra = [f for f in present if f not in expected]
        if missing:
            raise SystemExit(f"缺失输出文件: {missing}")
        if extra:
            raise SystemExit(f"多余输出文件: {extra}")
        print("\n".join(files))
        return
    raise SystemExit(f"unknown mode: {mode}")


if __name__ == "__main__":
    main()
