"""Honeycomb lattice geometry helpers mirroring src/honeycomb_lattice.h.

Mirrors the C++ HoneycombXCCylinder / HoneycombYCCylinder nearest-neighbor
link construction and ClassifyNearestNeighbor. Used by the pairing-matrix
tools to map bonds to the A/B/C directions without re-running C++.

Site index convention (column-major): site = x * Ly + y, x in [0, 2*Lx).
Both geometries are periodic along y (Ty links wrap) and open along x
(Tx links never wrap; the last column has no Tx link).
"""


def nn_bonds(geometry, Ly, Lx):
    """Mirror of the C++ nearest_neighbor_links construction.

    Returns a sorted list of sorted pairs (i, j) with i < j.
    """
    if geometry not in ("HoneycombXC", "HoneycombYC"):
        raise ValueError(f"unknown geometry: {geometry}")
    links = set()
    for i in range(Ly * Lx * 2):
        y, x = i % Ly, i // Ly
        Ty = (y + 1) % Ly + Ly * x
        Tx = y + Ly * ((x + 1) % (2 * Lx))
        if geometry == "HoneycombYC":
            if x < 2 * Lx - 1:
                if (x + y) % 2 == 0:
                    links.add((min(i, Ty), max(i, Ty)))
                    links.add((min(i, Tx), max(i, Tx)))
                else:
                    links.add((min(i, Tx), max(i, Tx)))
            elif y % 2 == 1:
                links.add((min(i, Ty), max(i, Ty)))
        else:  # HoneycombXC
            if x == 0:
                links.add((min(i, Ty), max(i, Ty)))
                if y % 2 == 0:
                    links.add((min(i, Tx), max(i, Tx)))
            elif x == 2 * Lx - 1:
                links.add((min(i, Ty), max(i, Ty)))
            elif (x + y) % 2 == 0:
                links.add((min(i, Ty), max(i, Ty)))
                links.add((min(i, Tx), max(i, Tx)))
            else:
                links.add((min(i, Ty), max(i, Ty)))
    return sorted(links)


def classify_nn_bond(geometry, Ly, Lx, i, j):
    """Mirror of ClassifyNearestNeighbor. Returns 'A', 'B' or 'C'.

    The input pair is sorted internally (the C++ call sites pass i < j).
    """
    if geometry not in ("HoneycombXC", "HoneycombYC"):
        raise ValueError(f"unknown geometry: {geometry}")
    if i == j:
        raise ValueError(f"degenerate bond: {i}")
    if i > j:
        i, j = j, i
    y0, x0 = i % Ly, i // Ly
    y1, x1 = j % Ly, j // Ly
    if geometry == "HoneycombYC":
        if y0 == y1:
            x_min = min(x0, x1)
            return "A" if (x_min + y0) % 2 == 0 else "B"
        if x0 != x1:
            raise ValueError(f"invalid YC bond: ({i}, {j})")
        return "C"  # vertical bond
    # HoneycombXC
    if x0 == x1:
        y_min = y0 if (y0 + Ly - y1) % Ly == 1 else y1
        return "A" if (y_min + x0) % 2 == 0 else "B"
    if y0 != y1:
        raise ValueError(f"invalid XC bond: ({i}, {j})")
    return "C"  # horizontal bond
