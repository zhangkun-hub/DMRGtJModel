"""Static CMake target/source integrity check (C1 regression).

Asserts every add_executable/add_library source listed in CMakeLists.txt
(same-line and continuation-line) exists on disk. Also asserts that the
removed legacy measureSC target is gone from CMakeLists.txt.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir, os.pardir))


def main():
    lines = open(os.path.join(ROOT, "CMakeLists.txt"),
                 encoding="utf-8").read().splitlines()
    missing = []
    targets = []
    i = 0
    while i < len(lines):
        m = re.match(r"\s*add_(executable|library)\((\w+)\s*(.*)", lines[i])
        if not m:
            i += 1
            continue
        tgt, rest = m.group(2), m.group(3)
        targets.append(tgt)
        sources = re.findall(r"src/[\w./-]+\.cpp", rest)
        j = i + 1
        while j < len(lines) and ")" not in rest:
            nxt = lines[j]
            sources += re.findall(r"src/[\w./-]+\.cpp", nxt)
            rest += nxt
            j += 1
        for s in sources:
            if not os.path.isfile(os.path.join(ROOT, s)):
                missing.append((tgt, s))
        i = j
    assert missing == [], f"missing sources: {missing}"
    assert "measureSC" not in targets, \
        "legacy measureSC target (removed source) must not reappear"
    assert "measureSC_fix_ref_ext" in targets
    assert "dmrg" in targets and "measure1" in targets and "measureS" in targets
    print(f"test_cmake_sources: ALL PASS ({len(targets)} targets, "
          f"0 missing sources)")


if __name__ == "__main__":
    main()
