/*
 * measureSC_fix_ref.cpp
 * Extended version:
 *   - Support --refbond=i,j for arbitrary nearest-neighbor reference bond
 *   - Backward compatible with --start
 * usage:
        mpirun -n 8 ./measureSC --refbond=i,j
 */

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <mpi.h>
#include <algorithm>
#include <ctime>
#include <functional>
#include <memory>
#include <unordered_set>
#include <utility>
#include <array>
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

} // namespace

/* ==================================================
 * Argument parser
 * ================================================== */
struct InputArgs {
  bool has_start = false;
  size_t start = 1;
  bool has_refbond = false;
  Bond refbond;
};

InputArgs ParseArgs(int argc, char **argv) {
  InputArgs args;
  for (int i = 1; i < argc; ++i) {
    std::string s(argv[i]);
    if (s.rfind("--start=", 0) == 0) {
      args.start = std::stoul(s.substr(8));
      args.has_start = true;
    }
    if (s.rfind("--refbond=", 0) == 0) {
      auto body = s.substr(10);
      auto p = body.find(',');
      if (p == std::string::npos) {
        std::cerr << "--refbond format: i,j\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      args.refbond = {std::stoul(body.substr(0, p)),
                      std::stoul(body.substr(p + 1))};
      args.has_refbond = true;
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
                      const Tensor &f) {
  size_t N = ComputeTotalSites(params);
  auto honeycomb = BuildHoneycombLattice(params);

  Bond ref_bond;
  if (args.has_refbond) {
    ref_bond = ValidateUserRefBond(args.refbond, N, params, honeycomb.get());
  } else {
    ref_bond = DetermineReferenceBond(args.start, N, params.Ly,
                                      params, honeycomb.get());
  }

  size_t ref_threshold = std::max(ref_bond[0], ref_bond[1]);
  auto tasks = GenerateTargetBonds(params, ref_threshold, honeycomb.get());

  SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  FiniteMPST mps(sites);

  const std::array<std::reference_wrapper<const std::array<Tensor, 4>>, 4>
      op_groups{opsA, opsB, opsC, opsD};
  const char *suffix[4] = {"a", "b", "c", "d"};

  std::string postfix =
      "_fix_refbond_" + std::to_string(ref_bond[0]) + "_" +
      std::to_string(ref_bond[1]);

  size_t idx = 0;
  for (const auto &task : tasks) {
    for (int o = 0; o < 4; ++o) {
      if ((idx * 4 + o) % mpi_size == static_cast<size_t>(rank)) {
        auto res = MeasureFourSiteOpGroup(
            mps, kMpsPath, op_groups[o].get(),
            ref_bond, task.bonds, f);
        DumpMeasuRes(res, task.label_prefix + suffix[o] + postfix);
      }
    }
    ++idx;
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
      std::cerr << "Usage: measureSC params.json [--start=i | --refbond=i,j]\n";
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

  PerformMeasureSC(rank, mpi_size, args,
                   params, opsA, opsB, opsC, opsD, ops.f);

  MPI_Finalize();
  return 0;
}