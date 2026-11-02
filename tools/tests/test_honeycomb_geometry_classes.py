"""Tests for the A1 spatial-grouping geometry (F2).

Validates the exact lattice displacement classes (dx_mid, dy_mid) that
src/measureSC_fix_ref_extension.cpp now uses for the spatial-spread record:

  dcol    = t0 / Ly - r0 / Ly                     (site-level column displacement)
  dx_mid  = (x_t0 + x_t1) - (x_r0 + x_r1)         (x = site / Ly, no wrap: x is open)
  dy_mid  = wrap((y_t0 + y_t1) - (y_r0 + y_r1), 2*Ly) into [-Ly, Ly)   (y = site % Ly)
  dx      = dx_mid * dx_unit, dy = dy_mid * dy_unit
  dx_unit = (px(2*Ly) - px(0)) / 4               (half column step)
  dy_unit = (py(2) - py(0)) / 4                  (half row step)

This Python mirror is the local test of the C++ math (the C++ binary cannot
be compiled on this machine; cluster compilation is the acceptance gate).
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 os.pardir))
import honeycomb_geometry as hg

S3 = math.sqrt(3.0)


def site_coord(geometry, Ly, site):
    y = site % Ly
    x = site // Ly
    is_A = ((x + y) % 2) == 0
    if geometry == "HoneycombXC":
        return x * S3 / 2.0 - (0.0 if is_A else 1.0 / (2.0 * S3)), y / 2.0
    return x / 2.0, y * S3 / 2.0 + (1.0 / (2.0 * S3) if is_A else 0.0)


def cpp_nn_loop(geometry, Ly, Lx):
    """Literal mirror of the C++ nearest_neighbor_links construction."""
    links = []
    for i in range(Ly * Lx * 2):
        y, x = i % Ly, i // Ly
        Ty = (y + 1) % Ly + Ly * x
        Tx = y + Ly * ((x + 1) % (2 * Lx))
        if geometry == "HoneycombYC":
            if x < 2 * Lx - 1:
                if (x + y) % 2 == 0:
                    links.append((min(i, Ty), max(i, Ty)))
                    links.append((min(i, Tx), max(i, Tx)))
                else:
                    links.append((min(i, Tx), max(i, Tx)))
            elif y % 2 == 1:
                links.append((min(i, Ty), max(i, Ty)))
        else:  # XC
            if x == 0:
                links.append((min(i, Ty), max(i, Ty)))
                if y % 2 == 0:
                    links.append((min(i, Tx), max(i, Tx)))
            elif x == 2 * Lx - 1:
                links.append((min(i, Ty), max(i, Ty)))
            elif (x + y) % 2 == 0:
                links.append((min(i, Ty), max(i, Ty)))
                links.append((min(i, Tx), max(i, Tx)))
            else:
                links.append((min(i, Ty), max(i, Ty)))
    return links


def wrap_class(raw, period, half):
    """C++ wrap: canonical [0, period) then shift into [-half, half)."""
    v = ((raw % period) + period) % period
    if v >= half:
        v -= period
    return v


def check_geometry(geometry, Ly, Lx):
    bonds = hg.nn_bonds(geometry, Ly, Lx)
    # (1) mirror vs C++ loop
    cpp = cpp_nn_loop(geometry, Ly, Lx)
    assert len(cpp) == len(set(cpp)), "C++ loop produces duplicate links"
    assert sorted(set(cpp)) == bonds, "Python nn_bonds differs from C++ loop"

    dx_unit = (site_coord(geometry, Ly, 2 * Ly)[0]
               - site_coord(geometry, Ly, 0)[0]) / 4.0
    dy_unit = (site_coord(geometry, Ly, 2)[1]
               - site_coord(geometry, Ly, 0)[1]) / 4.0

    n_pairs = 0
    classes = {}
    dcol_classes = {}
    for ref in bonds:
        r0, r1 = ref
        Xr = r0 // Ly + r1 // Ly
        Yr = r0 % Ly + r1 % Ly
        for tgt in bonds:
            t0, t1 = tgt
            if t0 <= r1:  # measurement filter (API constraint)
                continue
            n_pairs += 1
            dcol = t0 // Ly - r0 // Ly
            dX = t0 // Ly + t1 // Ly - Xr
            dY = wrap_class((t0 % Ly + t1 % Ly) - Yr, 2 * Ly, Ly)
            # (2) exact physical mapping, min-image range
            ty = (site_coord(geometry, Ly, t0)[1]
                  + site_coord(geometry, Ly, t1)[1]) / 2.0
            ry = (site_coord(geometry, Ly, r0)[1]
                  + site_coord(geometry, Ly, r1)[1]) / 2.0
            tx = (site_coord(geometry, Ly, t0)[0]
                  + site_coord(geometry, Ly, t1)[0]) / 2.0
            rx = (site_coord(geometry, Ly, r0)[0]
                  + site_coord(geometry, Ly, r1)[0]) / 2.0
            # dy derived from the exact class, no floating-point wrap
            dy = dY * dy_unit
            dx = dX * dx_unit
            # physical displacement is within the min-image range
            assert abs(dy) <= Ly * 2.0 * dy_unit / 2.0 + 1e-9, \
                f"min-image violation {ref}->{tgt}"
            # raw physical difference equals the class-derived value:
            # raw_dy = DeltaY_raw * dy_unit exactly, so rounding raw_dy/dy_unit
            # recovers DeltaY_raw and the canonical wrap must equal dY.
            raw_dy = ty - ry
            diff_wrapped = wrap_class(round(raw_dy / dy_unit), 2 * Ly, Ly)
            assert diff_wrapped == dY, \
                f"wrap mismatch {ref}->{tgt}: raw={raw_dy} dY={dY}"
            assert abs(dx - (tx - rx)) < 1e-9, f"dx mismatch {ref}->{tgt}"
            classes[(tuple(ref), dX, dY)] = classes.get(
                (tuple(ref), dX, dY), 0) + 1
            dcol_classes.setdefault((tuple(ref), dcol), set()).add(dY)
    # (3) same dcol can contain several dy_mid classes -> the new grouping
    #     separates them (counterexample must exist)
    merged = [k for k, v in dcol_classes.items() if len(v) > 1]
    assert merged, f"{geometry}: no same-dcol/different-dY counterexample"
    # (4) some classes contain several bonds (translations along the
    #     periodic direction) -> spread is meaningful
    multi = sum(1 for v in classes.values() if v > 1)
    assert multi > 0, f"{geometry}: no multi-bond displacement class"
    # (5) canonical wrap symmetry: raw +Ly and -Ly give the same class
    assert wrap_class(Ly, 2 * Ly, Ly) == wrap_class(-Ly, 2 * Ly, Ly)
    print(f"{geometry} Ly={Ly} Lx={Lx}: pairs={n_pairs} classes="
          f"{len(classes)} multi_bond_classes={multi} "
          f"same_dcol_counterexamples={len(merged)}")
    return bonds


def check_directions(geometry, Ly, Lx, bonds):
    # (6) every NN bond direction A/B/C occurs, classified consistently
    counts = {"A": 0, "B": 0, "C": 0}
    for i, j in bonds:
        counts[hg.classify_nn_bond(geometry, Ly, Lx, i, j)] += 1
    assert all(v > 0 for v in counts.values()), \
        f"{geometry}: missing direction, counts={counts}"


def main():
    for geometry in ("HoneycombXC", "HoneycombYC"):
        for Ly, Lx in ((4, 8), (4, 12)):
            bonds = check_geometry(geometry, Ly, Lx)
            check_directions(geometry, Ly, Lx, bonds)
    print("test_honeycomb_geometry_classes: ALL PASS")


if __name__ == "__main__":
    main()
