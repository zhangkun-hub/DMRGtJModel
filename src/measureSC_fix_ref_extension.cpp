/*
 * measureSC_fix_ref.cpp
 * Extended version:
 *   - Support --refbond=i,j for arbitrary nearest-neighbor reference bond
 *   - Backward compatible with --start
 * usage:
        mpirun -n 8 ./measureSC --refbond=i,j
 * Extended flags:
 *   --refbond=i1,j1;i2,j2;...  measure multiple reference bonds in one run
 *   --all-refbonds             use every NN bond as a reference bond
 *   --reverse-channels         also measure channels e-h (conjugates of a-d)
 */

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <functional>
#include <memory>
#include <unordered_set>
#include <utility>
#include <array>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <iostream>

#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"
#include "qlten/utility/timer.h"
#include "honeycomb_lattice.h"
#include "my_measure.h"

using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
using qlmps::SiteVec;
using qlmps::MeasureFourSiteOpGroup;
using qlmps::kMpsPath;
using qlmps::DumpMeasuRes;
using qlten::Timer;
using qlmps::MeasureBonds;

namespace {

using Bond = std::array<size_t, 2>;
using BondList = std::vector<Bond>;

struct MeasureTask {
  std::string label_prefix;
  BondList bonds;
};

/* --------------------------------------------------
 * Honeycomb lattice construction
 * -------------------------------------------------- */
std::unique_ptr<HoneycombLattice>
BuildHoneycombLattice(const CaseParams &params) {
  if (params.Geometry == "HoneycombXC") {
    return std::make_unique<HoneycombXCCylinder>(params.Ly, params.Lx);
  }
  if (params.Geometry == "HoneycombYC") {
    return std::make_unique<HoneycombYCCylinder>(params.Ly, params.Lx);
  }
  return nullptr;
}

bool ContainsNearestNeighbor(const HoneycombLattice &lattice,
                             size_t i, size_t j) {
  if (i > j) std::swap(i, j);
  for (const auto &link : lattice.nearest_neighbor_links) {
    size_t a = std::get<0>(link);
    size_t b = std::get<1>(link);
    if (a > b) std::swap(a, b);
    if (a == i && b == j) return true;
  }
  return false;
}

/* --------------------------------------------------
 * [NEW] Validate user-provided reference bond
 * -------------------------------------------------- */
Bond ValidateUserRefBond(const Bond &input,
                         size_t total_sites,
                         const CaseParams &params,
                         const HoneycombLattice *honeycomb) {
  size_t i = input[0];
  size_t j = input[1];

  if (i == j || i >= total_sites || j >= total_sites) {
    std::cerr << "[measureSC] Invalid --refbond indices\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (i > j) std::swap(i, j);

  if (honeycomb) {
    if (!ContainsNearestNeighbor(*honeycomb, i, j)) {
      std::cerr << "[measureSC] --refbond is not a honeycomb NN bond\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    // Square lattice NN check
    size_t xi = i / params.Ly, yi = i % params.Ly;
    size_t xj = j / params.Ly, yj = j % params.Ly;
    size_t dx = (xi > xj) ? xi - xj : xj - xi;
    size_t dy = (yi > yj) ? yi - yj : yj - yi;
    if (dx + dy != 1) {
      std::cerr << "[measureSC] --refbond is not a square NN bond\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
  return Bond{i, j};
}

/* --------------------------------------------------
 * Original reference-bond logic (unchanged)
 * -------------------------------------------------- */
Bond DetermineSquareReferenceBond(size_t start,
                                  size_t total_sites,
                                  size_t Ly,
                                  const CaseParams &) {
  if (start + 1 >= total_sites) {
    std::cerr << "Invalid start site\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return Bond{start, start + 1};
}

Bond DetermineHoneycombReferenceBond(size_t start,
                                     size_t total_sites,
                                     size_t Ly,
                                     const CaseParams &,
                                     const HoneycombLattice &lattice) {
  if (start + 1 < total_sites &&
      ContainsNearestNeighbor(lattice, start, start + 1)) {
    return Bond{start, start + 1};
  }
  if (start > 0 &&
      ContainsNearestNeighbor(lattice, start - 1, start)) {
    return Bond{start - 1, start};
  }
  std::cerr << "Failed to find honeycomb ref bond near start\n";
  MPI_Abort(MPI_COMM_WORLD, 1);
  return Bond{0, 1};
}

Bond DetermineReferenceBond(size_t start,
                            size_t total_sites,
                            size_t Ly,
                            const CaseParams &params,
                            const HoneycombLattice *honeycomb) {
  if (honeycomb) {
    return DetermineHoneycombReferenceBond(start, total_sites, Ly, params, *honeycomb);
  }
  return DetermineSquareReferenceBond(start, total_sites, Ly, params);
}

/* --------------------------------------------------
 * Target-bond generation (unchanged)
 * -------------------------------------------------- */
std::vector<MeasureTask>
GenerateTargetBonds(const CaseParams &params,
                    size_t ref_threshold,
                    const HoneycombLattice *honeycomb,
                    bool ignore_threshold = false) {
  std::vector<MeasureTask> tasks;

  if (honeycomb) {
    std::array<BondList, 3> groups;
    for (const auto &link : honeycomb->nearest_neighbor_links) {
      size_t i = std::get<0>(link);
      size_t j = std::get<1>(link);
      if (i > j) std::swap(i, j);
      if (!ignore_threshold && (i <= ref_threshold || j <= ref_threshold))
        continue;
      auto dir = honeycomb->ClassifyNearestNeighbor(i, j);
      groups[static_cast<int>(dir)].push_back({i, j});
    }
    const char *labels[3] = {"scsA", "scsB", "scsC"};
    for (int d = 0; d < 3; ++d) {
      if (!groups[d].empty())
        tasks.push_back({labels[d], std::move(groups[d])});
    }
  }
  return tasks;
}

/* --------------------------------------------------
 * A1: geometric descriptors (no MPS contraction needed).
 * These are data-quality records, NOT statistical errors.
 * -------------------------------------------------- */
struct BondGeometry {
  Bond bond;           // target bond {t0, t1}
  long dcol;           // site-level column displacement: t0/Ly - r0/Ly
  long dx_mid, dy_mid; // exact lattice displacement class of the bond midpoints
  double dx, dy;       // physical midpoint displacement (dy = dy_mid * dy_unit)
  double mx, my;       // target bond midpoint position (physical coords)
};

std::vector<BondGeometry> ComputeBondGeometry(
    const HoneycombLattice &lattice,
    const Bond &ref_bond,
    const BondList &target_bonds) {
  const size_t Ly = lattice.Ly;
  // Half column step and half row step, taken from same-parity site pairs
  // (sites 0 and 2 differ by two rows, sites 0 and 2*Ly by two columns, so
  // the A/B sublattice offsets cancel; valid for Ly>=4 and Lx>=2).
  const double dx_unit = (lattice.SiteCoordinate(2 * Ly).first -
                          lattice.SiteCoordinate(0).first) / 4.0;
  const double dy_unit = (lattice.SiteCoordinate(2).second -
                          lattice.SiteCoordinate(0).second) / 4.0;

  const long ref_xsum = ref_bond[0] / Ly + ref_bond[1] / Ly;
  const long ref_ysum = ref_bond[0] % Ly + ref_bond[1] % Ly;

  std::vector<BondGeometry> rows;
  rows.reserve(target_bonds.size());
  for (const auto &bond : target_bonds) {
    const long xsum = bond[0] / Ly + bond[1] / Ly;
    const long ysum = bond[0] % Ly + bond[1] % Ly;
    const long dx_mid = xsum - ref_xsum;
    // Wrap the y-sum difference into [-Ly, Ly): pure integer arithmetic,
    // no floating-point minimum-image wrap (which is ill-defined exactly at
    // half the period).
    long dy_mid = ysum - ref_ysum;
    const long period = 2 * static_cast<long>(Ly);
    dy_mid = ((dy_mid % period) + period) % period;
    if (dy_mid >= static_cast<long>(Ly)) dy_mid -= period;

    const auto tc0 = lattice.SiteCoordinate(bond[0]);
    const auto tc1 = lattice.SiteCoordinate(bond[1]);
    const double mx = (tc0.first + tc1.first) / 2.0;
    const double my = (tc0.second + tc1.second) / 2.0;
    rows.push_back({bond,
                    static_cast<long>(bond[0] / Ly) -
                        static_cast<long>(ref_bond[0] / Ly),
                    dx_mid, dy_mid,
                    dx_mid * dx_unit, dy_mid * dy_unit,
                    mx, my});
  }
  return rows;
}

// Value writer mirroring DumpAvgVal: scalar in the real build, [re, im] in
// the complex build, so the new files stay structurally identical to the
// main measurement data.
inline void DumpRealLike(std::ofstream &ofs, double v) {
#ifdef USE_REAL
  ofs << std::setw(14) << std::setprecision(12) << v;
#else
  ofs << "[" << std::setw(14) << std::setprecision(12) << v << ", "
      << std::setw(14) << std::setprecision(12) << 0.0 << "]";
#endif
}

// One row per midpoint lattice displacement class:
// [[dx_mid, dy_mid], [spread, 0]].
// spread = RMS deviation of the target-bond midpoints within the class.
void DumpSpatialGroupSpread(const std::string &basename,
                            const std::vector<BondGeometry> &rows) {
  std::map<std::pair<long, long>, std::vector<BondGeometry>> groups;
  for (const auto &row : rows) {
    groups[{row.dx_mid, row.dy_mid}].push_back(row);
  }
  std::ofstream ofs(basename + ".json");
  if (!ofs.is_open()) {
    std::cerr << "[measureSC] Failed to open " << basename << ".json\n";
    return;
  }
  ofs << "[\n";
  size_t gi = 0;
  for (const auto &[key, group] : groups) {
    double mx_mean = 0.0, my_mean = 0.0;
    for (const auto &row : group) {
      mx_mean += row.mx;
      my_mean += row.my;
    }
    mx_mean /= static_cast<double>(group.size());
    my_mean /= static_cast<double>(group.size());
    double var = 0.0;
    for (const auto &row : group) {
      const double dx = row.mx - mx_mean;
      const double dy = row.my - my_mean;
      var += dx * dx + dy * dy;
    }
    const double spread = std::sqrt(var / static_cast<double>(group.size()));
    ofs << "  [[" << key.first << ", " << key.second << "], ";
    DumpRealLike(ofs, spread);
    ofs << "]" << (gi + 1 < groups.size() ? "," : "") << "\n";
    ++gi;
  }
  ofs << "]\n";
}

// Per-bond geometry rows, sorted the same way as the main measurement files
// (by the first target site):
// [[r0, r1, t0, t1], [dcol, dx_mid, dy_mid, dx, dy, mx, my]].
void DumpBondGeometryFile(const std::string &basename,
                          const std::vector<BondGeometry> &rows) {
  auto sorted = rows;
  std::sort(sorted.begin(), sorted.end(),
            [](const BondGeometry &a, const BondGeometry &b) {
              return a.bond < b.bond;
            });
  std::ofstream ofs(basename + ".json");
  if (!ofs.is_open()) {
    std::cerr << "[measureSC] Failed to open " << basename << ".json\n";
    return;
  }
  ofs << "[\n";
  for (size_t i = 0; i < sorted.size(); ++i) {
    const auto &row = sorted[i];
    ofs << "  [[" << row.bond[0] << ", " << row.bond[1] << "], ["
        << std::setw(14) << std::setprecision(12) << static_cast<double>(row.dcol)
        << ", " << std::setw(14) << std::setprecision(12) << static_cast<double>(row.dx_mid)
        << ", " << std::setw(14) << std::setprecision(12) << static_cast<double>(row.dy_mid)
        << ", " << std::setw(14) << std::setprecision(12) << row.dx
        << ", " << std::setw(14) << std::setprecision(12) << row.dy
        << ", " << std::setw(14) << std::setprecision(12) << row.mx
        << ", " << std::setw(14) << std::setprecision(12) << row.my
        << "]]" << (i + 1 < sorted.size() ? "," : "") << "\n";
  }
  ofs << "]\n";
}

/* --------------------------------------------------
 * A2: marker that no repeat-sampling mechanism exists.
 * Explicitly NOT_APPLICABLE instead of fabricating a scatter record.
 * -------------------------------------------------- */
void DumpRepeatScatterManifest(const std::vector<Bond> &ref_bonds) {
  std::ofstream ofs("repeat_measurement_scatter.json");
  if (!ofs.is_open()) {
    std::cerr << "[measureSC] Failed to open repeat_measurement_scatter.json\n";
    return;
  }
  ofs << "{\n"
      << "  \"status\": \"NOT_APPLICABLE\",\n"
      << "  \"reason\": \"no repeat sampling mechanism exists: each target bond is measured by a single deterministic contraction (MeasureFourSiteOpGroup); no stochastic sampling or independent repeated measurements are performed\",\n"
      << "  \"refbonds\": [";
  for (size_t i = 0; i < ref_bonds.size(); ++i) {
    ofs << "[" << ref_bonds[i][0] << ", " << ref_bonds[i][1] << "]"
        << (i + 1 < ref_bonds.size() ? ", " : "");
  }
  ofs << "]\n}\n";
}

} // namespace

/* ==================================================
 * Argument parser
 * ================================================== */
struct InputArgs {
  bool has_start = false;
  size_t start = 1;
  bool has_refbond = false;
  std::vector<Bond> refbonds;     // [NEW] one or more reference bonds
  bool all_refbonds = false;      // [NEW] use every NN bond as a reference bond
  bool reverse_channels = false;  // [NEW] also measure channels e-h (conjugates of a-d)
};

InputArgs ParseArgs(int argc, char **argv) {
  InputArgs args;
  for (int i = 1; i < argc; ++i) {
    std::string s(argv[i]);
    if (s.rfind("--start=", 0) == 0) {
      args.start = std::stoul(s.substr(8));
      args.has_start = true;
    }
    if (s == "--all-refbonds") {
      args.all_refbonds = true;
    }
    if (s == "--reverse-channels") {
      args.reverse_channels = true;
    }
    if (s.rfind("--refbond=", 0) == 0) {
      auto body = s.substr(10);
      size_t pos = 0;
      while (true) {
        auto next = body.find(';', pos);
        std::string pair_str =
            body.substr(pos, next == std::string::npos ? std::string::npos
                                                       : next - pos);
        auto p = pair_str.find(',');
        if (p == std::string::npos || pair_str.empty()) {
          std::cerr << "--refbond format: i1,j1;i2,j2;...\n";
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        args.refbonds.push_back({std::stoul(pair_str.substr(0, p)),
                                 std::stoul(pair_str.substr(p + 1))});
        args.has_refbond = true;
        if (next == std::string::npos) break;
        pos = next + 1;
      }
    }
  }
  return args;
}

/* ==================================================
 * Main measurement routine
 * ================================================== */
void PerformMeasureSC(int rank,
                      int mpi_size,
                      const InputArgs &args,
                      const CaseParams &params,
                      const std::array<Tensor, 4> &opsA,
                      const std::array<Tensor, 4> &opsB,
                      const std::array<Tensor, 4> &opsC,
                      const std::array<Tensor, 4> &opsD,
                      const std::array<Tensor, 4> &opsE,
                      const std::array<Tensor, 4> &opsF,
                      const std::array<Tensor, 4> &opsG,
                      const std::array<Tensor, 4> &opsH,
                      const Tensor &f) {
  size_t N = ComputeTotalSites(params);
  auto honeycomb = BuildHoneycombLattice(params);

  // [NEW] Build the reference-bond list. A single --refbond keeps the
  // original one-reference-bond behavior.
  std::vector<Bond> ref_bonds;
  if (args.all_refbonds) {
    if (!honeycomb) {
      std::cerr << "[measureSC] --all-refbonds requires a honeycomb lattice\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (const auto &link : honeycomb->nearest_neighbor_links) {
      size_t i = std::get<0>(link);
      size_t j = std::get<1>(link);
      if (i > j) std::swap(i, j);
      ref_bonds.push_back({i, j});
    }
    std::sort(ref_bonds.begin(), ref_bonds.end());
  } else if (args.has_refbond) {
    for (const auto &bond : args.refbonds) {
      ref_bonds.push_back(ValidateUserRefBond(bond, N, params, honeycomb.get()));
    }
  } else {
    ref_bonds.push_back(DetermineReferenceBond(args.start, N, params.Ly,
                                               params, honeycomb.get()));
  }

  SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  FiniteMPST mps(sites);

  const std::array<std::reference_wrapper<const std::array<Tensor, 4>>, 8>
      op_groups{opsA, opsB, opsC, opsD, opsE, opsF, opsG, opsH};
  // Channels a-d: creation at the reference pair, annihilation at the target
  // pair. Channels e-h: the reversed roles, enabling the conjugation check
  // P_{ab}(i,j) = P*_{ba}(j,i).
  const char *suffix[8] = {"a", "b", "c", "d", "e", "f", "g", "h"};
  const size_t num_channels = args.reverse_channels ? 8 : 4;

  // A2: record that no repeat-sampling mechanism exists (do not fabricate).
  if (rank == 0) {
    DumpRepeatScatterManifest(ref_bonds);
  }

  size_t idx = 0;
  for (const auto &ref_bond : ref_bonds) {
    size_t ref_threshold = std::max(ref_bond[0], ref_bond[1]);
    auto tasks = GenerateTargetBonds(params, ref_threshold, honeycomb.get());

    std::string postfix =
        "_fix_refbond_" + std::to_string(ref_bond[0]) + "_" +
        std::to_string(ref_bond[1]);

    // A1: geometric descriptors (spatial group spread + per-bond geometry).
    // The discrete lattice classes need same-parity site pairs (Lx>=2,
    // Ly>=4); skip on degenerate geometries instead of emitting wrong numbers.
    const bool geometry_ok = (honeycomb->Lx >= 2 && honeycomb->Ly >= 4);
    if (rank == 0) {
      if (!geometry_ok && !tasks.empty()) {
        std::cerr << "[measureSC] A1 geometry records require Lx>=2 and "
                     "Ly>=4; skipped\n";
      }
      if (geometry_ok) {
        for (const auto &task : tasks) {
          auto geometry = ComputeBondGeometry(*honeycomb, ref_bond, task.bonds);
          DumpSpatialGroupSpread(task.label_prefix + postfix + "_spatial_spread",
                                 geometry);
          DumpBondGeometryFile(task.label_prefix + postfix + "_bond_geometry",
                               geometry);
        }
      }
    }

    for (const auto &task : tasks) {
      for (int o = 0; o < static_cast<int>(num_channels); ++o) {
        if ((idx * num_channels + o) % mpi_size ==
            static_cast<size_t>(rank)) {
          auto res = MeasureFourSiteOpGroup(
              mps, kMpsPath, op_groups[o].get(),
              ref_bond, task.bonds, f);
          DumpMeasuRes(res, task.label_prefix + suffix[o] + postfix);
        }
      }
      ++idx;
    }
  }
}

/* ==================================================
 * main
 * ================================================== */
int main(int argc, char **argv) {
  MPI_Init(nullptr, nullptr);

  int rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (argc < 2) {
    if (rank == 0)
      std::cerr << "Usage: measureSC params.json [--start=i | "
                   "--refbond=i,j[;k,l;...] | --all-refbonds] "
                   "[--reverse-channels]\n";
    MPI_Finalize();
    return 1;
  }

  CaseParams params(argv[1]);
  InputArgs args = ParseArgs(argc, argv);

  tJOperators ops;
  std::array<Tensor, 4> opsA{ops.bupc, ops.bdnc, ops.bupa, ops.bdna};
  std::array<Tensor, 4> opsB{ops.bdnc, ops.bupc, ops.bupa, ops.bdna};
  std::array<Tensor, 4> opsC{ops.bupc, ops.bdnc, ops.bdna, ops.bupa};
  std::array<Tensor, 4> opsD{ops.bdnc, ops.bupc, ops.bdna, ops.bupa};
  // Reverse channels: annihilation at the reference pair, creation at the
  // target pair (conjugates of a-d).
  std::array<Tensor, 4> opsE{ops.bupa, ops.bdna, ops.bupc, ops.bdnc};
  std::array<Tensor, 4> opsF{ops.bdna, ops.bupa, ops.bupc, ops.bdnc};
  std::array<Tensor, 4> opsG{ops.bupa, ops.bdna, ops.bdnc, ops.bupc};
  std::array<Tensor, 4> opsH{ops.bdna, ops.bupa, ops.bdnc, ops.bupc};

  PerformMeasureSC(rank, mpi_size, args,
                   params, opsA, opsB, opsC, opsD,
                   opsE, opsF, opsG, opsH, ops.f);

  MPI_Finalize();
  return 0;
}