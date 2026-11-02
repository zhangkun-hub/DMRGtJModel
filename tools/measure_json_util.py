"""Shared helpers to load DMRGtJModel measurement JSON files.

The measurement JSONs are plain nested-list text files, one row per measured
event:

    complex build:  [[[r0, r1, t0, t1], [re, im]], ...]
    real build:     [[[r0, r1, t0, t1], 0.123], ...]

This mirrors the consumer-side contract (ast.literal_eval based parsing,
scalar == real value, two-element list == [re, im]).
"""
import ast


def load_measu_rows(path):
    """Parse a measurement JSON file into a list of (sites, value).

    sites: list of ints (site indices of the operators).
    value: complex; real-build scalars are treated as reals.
    """
    with open(path, "r", encoding="utf-8") as f:
        data = ast.literal_eval(f.read())
    rows = []
    for entry in data:
        sites = [int(s) for s in entry[0]]
        val = entry[1]
        if isinstance(val, (list, tuple)):
            value = complex(float(val[0]), float(val[1]))
        else:
            value = complex(float(val), 0.0)
        rows.append((sites, value))
    return rows
