// SPDX-License-Identifier: LGPL-3.0-only
//
// Honeycomb lattice helper classes adapted from
// /Users/wanghaoxin/GitHub/HoneycombHeisenbergJ1J2J3/src/honeycomb_lattice.h
// with style adjustments for the SquareLatticetJmodel project.

#ifndef TJMODEL_SRC_HONEYCOMB_LATTICE_H
#define TJMODEL_SRC_HONEYCOMB_LATTICE_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "lattice_base.h"

class HoneycombLattice : public LatticeBase {
  public:
    HoneycombLattice() = default;
    HoneycombLattice(size_t Ly, size_t Lx) : LatticeBase(Ly, Lx, 2) {
    }

    static constexpr double kSqrt3 = 1.73205080756887729353;
    static constexpr double kSqrt3Over2 = kSqrt3 / 2.0;

    enum class BondDirection : size_t { kA = 0, kB = 1, kC = 2 };

    virtual BondDirection ClassifyNearestNeighbor(size_t site0, size_t site1) const = 0;
    std::string DirectionColor(BondDirection direction) const;
};

class HoneycombYCCylinder : public HoneycombLattice {
  public:
    HoneycombYCCylinder() = default;
    HoneycombYCCylinder(size_t Ly, size_t Lx);

    std::pair<double, double> SiteCoordinate(size_t site) const override {
      const size_t y = site % Ly;
      const size_t x = site / Ly;
      const bool a_sublattice = ((x + y) % 2) == 0;

      const double px = static_cast<double>(x) / 2.0;
      const double py = static_cast<double>(y) * kSqrt3Over2 +
          (a_sublattice ? 1.0 / (2.0 * kSqrt3) : 0.0);
      return {px, py};
    }

    BondDirection ClassifyNearestNeighbor(size_t site0, size_t site1) const override;
};

class HoneycombXCCylinder : public HoneycombLattice {
  public:
    HoneycombXCCylinder() = default;
    HoneycombXCCylinder(size_t Ly, size_t Lx);

    std::pair<double, double> SiteCoordinate(size_t site) const override {
      const size_t y = site % Ly;
      const size_t x = site / Ly;
      const bool a_sublattice = ((x + y) % 2) == 0;

      const double px = static_cast<double>(x) * kSqrt3Over2 -
          (a_sublattice ? 0.0 : 1.0 / (2.0 * kSqrt3));
      const double py = static_cast<double>(y) / 2.0;
      return {px, py};
    }

    BondDirection ClassifyNearestNeighbor(size_t site0, size_t site1) const override;
};

inline HoneycombYCCylinder::HoneycombYCCylinder(size_t Ly, size_t Lx)
  : HoneycombLattice(Ly, Lx) {
  if (Ly % 2 != 0) {
    std::cout << "Ly = " << Ly << std::endl;
    std::cout << "Unexpected odd Ly!" << std::endl;
    exit(1);
  }
  nearest_neighbor_links.reserve(3 * N / 2);
  next_nearest_neighbor_links.reserve(3 * N);
  third_nearest_neighbor_links.reserve(3 * N / 2);
  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly;
    const size_t x = i / Ly;
    const size_t Tx = y + Ly * ((x + 1) % (2 * Lx));
    const size_t Txx = y + Ly * ((x + 2) % (2 * Lx));
    const size_t Ty = (y + 1) % Ly + Ly * x;
    const size_t Txy = (y + 1) % Ly + Ly * ((x + 1) % (2 * Lx));
    const size_t Txmy = (y + Ly - 1) % Ly + Ly * ((x + 1) % (2 * Lx));
    const size_t Txxy = (y + 1) % Ly + Ly * ((x + 2) % (2 * Lx));
    const size_t Txxmy = (y + Ly - 1) % Ly + Ly * ((x + 2) % (2 * Lx));

    if (x < 2 * Lx - 1) {
      if ((x + y) % 2 == 0) {
        nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
        nearest_neighbor_links.emplace_back(std::min(i, Tx), std::max(i, Tx));
      } else {
        nearest_neighbor_links.emplace_back(std::min(i, Tx), std::max(i, Tx));
      }
    } else if (y % 2 == 1) {
      nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
    }

    if (x < 2 * Lx - 2) {
      next_nearest_neighbor_links.emplace_back(std::min(i, Txx), std::max(i, Txx));
      next_nearest_neighbor_links.emplace_back(std::min(i, Txy), std::max(i, Txy));
      next_nearest_neighbor_links.emplace_back(std::min(i, Txmy), std::max(i, Txmy));
    } else if (x == 2 * Lx - 2) {
      next_nearest_neighbor_links.emplace_back(std::min(i, Txy), std::max(i, Txy));
      next_nearest_neighbor_links.emplace_back(std::min(i, Txmy), std::max(i, Txmy));
    }

    if (x < 2 * Lx - 2) {
      if ((x + y) % 2 == 1) {
        third_nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
        third_nearest_neighbor_links.emplace_back(std::min(i, Txxmy), std::max(i, Txxmy));
      } else {
        third_nearest_neighbor_links.emplace_back(std::min(i, Txxy), std::max(i, Txxy));
      }
    } else if ((x + y) % 2 == 1) {
      third_nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
    }
  }
}

inline HoneycombXCCylinder::HoneycombXCCylinder(size_t Ly, size_t Lx)
  : HoneycombLattice(Ly, Lx) {
  if (Ly % 2 != 0) {
    std::cout << "Ly = " << Ly << std::endl;
    std::cout << "Unexpected odd Ly!" << std::endl;
    exit(1);
  }
  nearest_neighbor_links.reserve(3 * N / 2);
  next_nearest_neighbor_links.reserve(3 * N);
  third_nearest_neighbor_links.reserve(3 * N / 2);
  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly;
    const size_t x = i / Ly;
    const size_t Tx = y + Ly * ((x + 1) % (2 * Lx));
    const size_t Ty = (y + 1) % Ly + Ly * x;
    const size_t Tyy = (y + 2) % Ly + Ly * x;
    const size_t Txy = (y + 1) % Ly + Ly * ((x + 1) % (2 * Lx));
    const size_t Txyy = (y + 2) % Ly + Ly * ((x + 1) % (2 * Lx));
    const size_t Txmyy = (y + Ly - 2) % Ly + Ly * ((x + 1) % (2 * Lx));
    const size_t Txmy = (y + Ly - 1) % Ly + Ly * ((x + 1) % (2 * Lx));

    if (x == 0) {
      nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
      if (y % 2 == 0) {
        nearest_neighbor_links.emplace_back(std::min(i, Tx), std::max(i, Tx));
      }
    } else if (x == 2 * Lx - 1) {
      nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
    } else if ((x + y) % 2 == 0) {
      nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
      nearest_neighbor_links.emplace_back(std::min(i, Tx), std::max(i, Tx));
    } else {
      nearest_neighbor_links.emplace_back(std::min(i, Ty), std::max(i, Ty));
    }

    if (x < 2 * Lx - 1) {
      next_nearest_neighbor_links.emplace_back(std::min(i, Tyy), std::max(i, Tyy));
      next_nearest_neighbor_links.emplace_back(std::min(i, Txy), std::max(i, Txy));
      next_nearest_neighbor_links.emplace_back(std::min(i, Txmy), std::max(i, Txmy));
    } else {
      next_nearest_neighbor_links.emplace_back(std::min(i, Tyy), std::max(i, Tyy));
    }

    if (x < 2 * Lx - 1) {
      if ((x + y) % 2 == 0) {
        third_nearest_neighbor_links.emplace_back(std::min(i, Txyy), std::max(i, Txyy));
        third_nearest_neighbor_links.emplace_back(std::min(i, Txmyy), std::max(i, Txmyy));
      } else {
        third_nearest_neighbor_links.emplace_back(std::min(i, Tx), std::max(i, Tx));
      }
    }
  }
}

inline HoneycombLattice::BondDirection HoneycombYCCylinder::ClassifyNearestNeighbor(
  size_t site0,
  size_t site1) const {
  if (site0 == site1) {
    std::cerr << "[HoneycombYCCylinder] degenerate bond: " << site0 << '\n';
    std::abort();
  }

  if (site0 > site1) {
    std::swap(site0, site1);
  }

  const size_t y0 = site0 % Ly;
  const size_t x0 = site0 / Ly;
  const size_t y1 = site1 % Ly;
  const size_t x1 = site1 / Ly;

  if (y0 == y1) {
    const size_t x_min = std::min(x0, x1);
    const size_t parity = (x_min + y0) % 2;
    return parity == 0 ? BondDirection::kA : BondDirection::kB;
  }

  if (x0 != x1) {
    std::cerr << "[HoneycombYCCylinder] invalid bond: (" << site0 << ", "
        << site1 << ")" << '\n';
    std::abort();
  }

  return BondDirection::kC;//vertical bond
}

inline HoneycombLattice::BondDirection HoneycombXCCylinder::ClassifyNearestNeighbor(
  size_t site0,
  size_t site1) const {
  if (site0 == site1) {
    std::cerr << "[HoneycombXCCylinder] degenerate bond: " << site0 << '\n';
    std::abort();
  }

  if (site0 > site1) {
    std::swap(site0, site1);
  }

  const size_t y0 = site0 % Ly;
  const size_t x0 = site0 / Ly;
  const size_t y1 = site1 % Ly;
  const size_t x1 = site1 / Ly;

  if (x0 == x1) {
    const size_t y_min = (y0 + Ly - y1) % Ly == 1 ? y0 : y1;
    const size_t parity = (y_min + x0) % 2;
    return parity == 0 ? BondDirection::kA : BondDirection::kB;
  }

  if (y0 != y1) {
    std::cerr << "[HoneycombXCCylinder] invalid bond: (" << site0 << ", "
        << site1 << ")" << '\n';
    std::abort();
  }

  return BondDirection::kC; // horizontal bond
}

inline std::string HoneycombLattice::DirectionColor(BondDirection direction) const {
  switch (direction) {
    case BondDirection::kA:
      return "#1f77b4";
    case BondDirection::kB:
      return "#2ca02c";
    case BondDirection::kC:
      return "#ff7f0e";
  }
  return "#444444";
}

#endif // TJMODEL_SRC_HONEYCOMB_LATTICE_H
