/*
    measureSC.cpp
    for measure pair correlation function. memory optimized and parallel version.
    usage:
        mpirun -n 8 ./measureSC --start=start_site
*/
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <algorithm>
#include <ctime>
#include <functional>
#include <memory>
#include <unordered_set>
#include <utility>
#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"
#include "qlten/utility/timer.h"
#include "honeycomb_lattice.h"

using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
using qlmps::SiteVec;
using qlmps::MeasureFourSiteOpGroup;
using qlmps::kMpsPath;
using qlmps::DumpMeasuRes;
using qlten::Timer;

namespace {

struct PairHash {
  size_t operator()(const std::pair<size_t, size_t> &p) const noexcept {
    size_t lhs = std::hash<size_t>{}(p.first);
    size_t rhs = std::hash<size_t>{}(p.second);
    // Fowler–Noll–Vo inspired mix to reduce collisions.
    return lhs ^ (rhs + 0x9e3779b97f4a7c15ULL + (lhs << 6) + (lhs >> 2));
  }
};

using Bond = std::array<size_t, 2>;
using BondList = std::vector<Bond>;

struct MeasureTask {
  std::string label_prefix;
  BondList bonds;
};

std::unique_ptr<HoneycombLattice> BuildHoneycombLattice(const CaseParams &params) {
  if (params.Geometry == "HoneycombXC") {
    return std::make_unique<HoneycombXCCylinder>(params.Ly, params.Lx);
  }
  if (params.Geometry == "HoneycombYC") {
    return std::make_unique<HoneycombYCCylinder>(params.Ly, params.Lx);
  }
  return nullptr;
}

bool ContainsNearestNeighbor(const HoneycombLattice &lattice,
                             size_t site0,
                             size_t site1) {
  if (site0 > site1) {
    std::swap(site0, site1);
  }
  for (const auto &link : lattice.nearest_neighbor_links) {
    size_t lhs = std::get<0>(link);
    size_t rhs = std::get<1>(link);
    if (lhs > rhs) {
      std::swap(lhs, rhs);
    }
    if (lhs == site0 && rhs == site1) {
      return true;
    }
  }
  return false;
}

Bond DetermineSquareReferenceBond(size_t requested_start,
                                  size_t total_sites,
                                  size_t Ly,
                                  const CaseParams &params) {
  if (requested_start + 1 >= total_sites) {
    std::cerr << "[measureSC_fix_ref] Invalid reference start " << requested_start
              << " (exceeds system size)." << std::endl;
    exit(1);
  }
  size_t site0 = requested_start;
  const size_t y0 = site0 % Ly;
  if (y0 == Ly - 1) {
    if (site0 == 0) {
      std::cerr << "[measureSC_fix_ref] Cannot shift reference bond for start=0." << std::endl;
      exit(1);
    }
    std::cerr << "[measureSC_fix_ref] Warning: reference bond {" << requested_start
              << ", " << requested_start + 1
              << "} is not a valid vertical bond on geometry " << params.Geometry
              << ". Using {" << requested_start - 1 << ", " << requested_start
              << "} instead." << std::endl;
    --site0;
  }
  const size_t site1 = site0 + 1;
  if (site1 >= total_sites || (site0 / Ly) != (site1 / Ly)) {
    std::cerr << "[measureSC_fix_ref] Failed to build a valid square-lattice reference bond." << std::endl;
    exit(1);
  }
  return Bond{site0, site1};
}

Bond DetermineHoneycombReferenceBond(size_t requested_start,
                                     size_t total_sites,
                                     size_t Ly,
                                     const CaseParams &params,
                                     const HoneycombLattice &lattice) {
  auto is_valid = [&](size_t site) -> bool {
    if (site + 1 >= total_sites) {
      return false;
    }
    if ((site / Ly) != ((site + 1) / Ly)) {
      return false;
    }
    if (!ContainsNearestNeighbor(lattice, site, site + 1)) {
      return false;
    }
    if (params.Geometry == "HoneycombYC") {
      const size_t x = site / Ly;
      const size_t y = site % Ly;
      if ((x % 2 != 0) || (y % 2 != 0)) {
        return false;
      }
      const auto direction = lattice.ClassifyNearestNeighbor(site, site + 1);
      if (direction != HoneycombLattice::BondDirection::kC) {
        return false;
      }
    }
    return true;
  };

  if (requested_start >= total_sites) {
    std::cerr << "[measureSC_fix_ref] Invalid reference start " << requested_start
              << " (exceeds system size)." << std::endl;
    exit(1);
  }

  size_t candidate = requested_start;
  if (params.Geometry == "HoneycombXC") {
    if (!is_valid(candidate)) {
      if (candidate == 0 || !is_valid(candidate - 1)) {
        std::cerr << "[measureSC_fix_ref] HoneycombXC requires reference bond to be a real link; starting site="
                  << requested_start << " fails both {s,s+1} and {s-1,s}." << std::endl;
        exit(1);
      }
      candidate -= 1;
      std::cerr << "[measureSC_fix_ref] Warning: reference bond {" << requested_start
                << ", " << requested_start + 1 << "} is invalid on HoneycombXC."
                << " Using {" << candidate << ", " << candidate + 1 << "} instead." << std::endl;
    }
  } else {
    if (!is_valid(candidate)) {
      bool found = false;
      const size_t max_delta = total_sites - 1;
      for (size_t delta = 1; delta <= max_delta; ++delta) {
        if (candidate >= delta) {
          const size_t below = candidate - delta;
          if (is_valid(below)) {
            candidate = below;
            found = true;
            break;
          }
        }
        const size_t above = candidate + delta;
        if (above < total_sites && is_valid(above)) {
          candidate = above;
          found = true;
          break;
        }
      }
      if (!found) {
        std::cerr << "[measureSC_fix_ref] Failed to locate a valid honeycomb reference bond near site "
                  << requested_start << " for geometry " << params.Geometry << std::endl;
        exit(1);
      }
      std::cerr << "[measureSC_fix_ref] Warning: reference bond {" << requested_start
                << ", " << requested_start + 1 << "} adjusted to {" << candidate
                << ", " << candidate + 1 << "} for geometry " << params.Geometry << std::endl;
    }
  }

  return Bond{candidate, candidate + 1};
}

Bond DetermineReferenceBond(size_t requested_start,
                            size_t total_sites,
                            size_t Ly,
                            const CaseParams &params,
                            const HoneycombLattice *honeycomb) {
  if (honeycomb != nullptr) {
    return DetermineHoneycombReferenceBond(requested_start, total_sites, Ly, params, *honeycomb);
  }
  return DetermineSquareReferenceBond(requested_start, total_sites, Ly, params);
}

std::array<BondList, 3> BuildHoneycombTargetBonds(const HoneycombLattice &lattice,
                                                  size_t ref_threshold) {
  std::array<BondList, 3> targets;
  std::unordered_set<std::pair<size_t, size_t>, PairHash> seen;
  seen.reserve(lattice.nearest_neighbor_links.size());
  for (const auto &link : lattice.nearest_neighbor_links) {
    size_t site0 = std::get<0>(link);
    size_t site1 = std::get<1>(link);
    if (site0 > site1) {
      std::swap(site0, site1);
    }
    if (site0 <= ref_threshold || site1 <= ref_threshold) {
      continue;
    }
    if (!seen.insert({site0, site1}).second) {
      continue;
    }
    const auto direction = lattice.ClassifyNearestNeighbor(site0, site1);
    size_t index = 0;
    switch (direction) {
      case HoneycombLattice::BondDirection::kA:
        index = 0;
        break;
      case HoneycombLattice::BondDirection::kB:
        index = 1;
        break;
      case HoneycombLattice::BondDirection::kC:
        index = 2;
        break;
    }
    targets[index].push_back({site0, site1});
  }
  return targets;
}

} // namespace

void PerformMeasureSC(int rank,
                      int mpi_size,
                      size_t start,
                      const CaseParams &params,
                      const std::array<Tensor, 4> &sc_phys_ops_a,
                      const std::array<Tensor, 4> &sc_phys_ops_b,
                      const std::array<Tensor, 4> &sc_phys_ops_c,
                      const std::array<Tensor, 4> &sc_phys_ops_d,
                      const Tensor &f);

// When used to measure, note if should not set start too small to exceed canonical center.
bool Parser(const int argc, char *argv[],
            size_t &start,
            size_t &end) {
  int nOptionIndex = 1;

  std::string arguement1 = "--start=";
  std::string arguement2 = "--end=";
  bool start_argument_has(false), end_argument_has(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement1.size()];
      start = atoi(para_string.c_str());
      start_argument_has = true;
    } else if (strncmp(argv[nOptionIndex], arguement2.c_str(), arguement2.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement2.size()];
      end = atoi(para_string.c_str());
      end_argument_has = true;
    }
    nOptionIndex++;
  }

  if (start_argument_has != end_argument_has) {
    std::cout << "Only setting one start/end argument, exit(1)." << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    std::cout << "Note: no start/end argument, set it by default (L/4, 3*L/4+2)." << std::endl;
  }

  return start_argument_has;
}

int main(int argc, char *argv[]) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  if (mpi_size > 8) {
    std::cout << "too many MPI processes. Set -np 8" << std::endl;
  }
  clock_t startTime, endTime;
  startTime = clock();

  size_t start;
  size_t endx;
  bool start_argument_has = Parser(argc, argv, start, endx);

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = ComputeTotalSites(params);
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    start = 1; // default start site for pairing measurement
    endx = start + Lx / 2 + 2;
  }

  tJOperators ops;

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif
  Timer foursite_timer("measure four site operators");
  std::array<Tensor, 4> sc_phys_ops_a = {ops.bupc, ops.bdnc, ops.bupa, ops.bdna};
  std::array<Tensor, 4> sc_phys_ops_b = {ops.bdnc, ops.bupc, ops.bupa, ops.bdna};
  std::array<Tensor, 4> sc_phys_ops_c = {ops.bupc, ops.bdnc, ops.bdna, ops.bupa}; //<B>=<C>, by spin inversion symmetry
  std::array<Tensor, 4> sc_phys_ops_d = {ops.bdnc, ops.bupc, ops.bdna, ops.bupa}; //<A>=<D>, by spin inversion symmetry

  PerformMeasureSC(rank, mpi_size, start, params, sc_phys_ops_a, sc_phys_ops_b, sc_phys_ops_c, sc_phys_ops_d, ops.f);

  std::cout << "measured SC correlation function.<====" << std::endl;
  foursite_timer.PrintElapsed();
#if SYM_LEVEL == 0
  /**
   * Effectively measure local pairing:
   * (A-B)/sqrt(2) gives the singlet pair
   * (A+B)/sqrt(2) gives the triplet pair (except the boundary, boundary should count for an additional sign)
   * The \pm 1 components of the local triplet pairs are forbidden by spin U1 symmetry
   * which should degenerate with C and D if SU(2) symmetry is restored.
   */
  sc_phys_ops_a = {ops.id, ops.id, ops.bupa, ops.bdna};
  sc_phys_ops_b = {ops.id, ops.id, ops.bdna, ops.bupa};
  sc_phys_ops_c = {ops.id, ops.id, ops.bupc, ops.bupc};
  sc_phys_ops_d = {ops.id, ops.id, ops.bdnc, ops.bdnc};
  PerformMeasureSC(rank, mpi_size, 0, params, sc_phys_ops_a, sc_phys_ops_b, sc_phys_ops_c, sc_phys_ops_d, ops.f);
#endif
  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}
void PerformMeasureSC(int rank,
                      int mpi_size,
                      size_t start,
                      const CaseParams &params,
                      const std::array<Tensor, 4> &sc_phys_ops_a,
                      const std::array<Tensor, 4> &sc_phys_ops_b,
                      const std::array<Tensor, 4> &sc_phys_ops_c,
                      const std::array<Tensor, 4> &sc_phys_ops_d,
                      const Tensor &f) {
  const size_t Lx = params.Lx, Ly = params.Ly;
  const size_t total_sites = ComputeTotalSites(params);
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(total_sites);
  FiniteMPST mps(sites);
  auto honeycomb_lattice = BuildHoneycombLattice(params);
  const Bond ref_bond = DetermineReferenceBond(start, total_sites, Ly, params, honeycomb_lattice.get());
  const size_t ref_threshold = std::max(ref_bond[0], ref_bond[1]);

  std::vector<MeasureTask> tasks;
  tasks.reserve(4);

  if (honeycomb_lattice) {
    auto targets = BuildHoneycombTargetBonds(*honeycomb_lattice, ref_threshold);
    const std::array<std::string, 3> labels{"scsA", "scsB", "scsC"};
    for (size_t idx = 0; idx < targets.size(); ++idx) {
      if (!targets[idx].empty()) {
        tasks.push_back(MeasureTask{labels[idx], std::move(targets[idx])});
      }
    }
  } else {
    BondList target_y_bond;
    BondList target_x_bond;
    target_y_bond.reserve(Lx * Ly);
    target_x_bond.reserve(Lx * Ly);

    for (size_t x = 0; x < Lx; ++x) {
      for (size_t y = 0; y < Ly; ++y) {
        const size_t site0 = x * Ly + y;
        const size_t site1 = (y + 1) % Ly + Ly * x;
        if (site0 <= ref_threshold || site1 <= ref_threshold) {
          continue;
        }
        if (params.Geometry == "SquareOBC") {
          if (y < Ly - 1) {
            target_y_bond.push_back({site0, site1});
          }
        } else if (params.Geometry == "SquareCylinder") {
          if (y == 0) {
            target_y_bond.push_back({site0, site1});
            const size_t wrap_site = site0 + Ly - 1;
            if (wrap_site <= ref_threshold) {
              continue;
            }
            target_y_bond.push_back({site0, wrap_site});
          } else if (y < Ly - 1) {
            target_y_bond.push_back({site0, site1});
          }
        }
      }
    }

    for (size_t x = 0; x + 1 < Lx; ++x) {
      for (size_t y = 0; y < Ly; ++y) {
        const size_t site0 = x * Ly + y;
        const size_t site1 = y + Ly * (x + 1);
        if (site0 <= ref_threshold || site1 <= ref_threshold) {
          continue;
        }
        target_x_bond.push_back({site0, site1});
      }
    }

    if (!target_y_bond.empty()) {
      tasks.push_back(MeasureTask{"scsyy", std::move(target_y_bond)});
    }
    if (!target_x_bond.empty()) {
      tasks.push_back(MeasureTask{"scsyx", std::move(target_x_bond)});
    }
  }

  if (tasks.empty()) {
    std::cerr << "[measureSC_fix_ref] No target bonds generated for geometry "
              << params.Geometry << "." << std::endl;
    exit(1);
  }

  const std::array<std::reference_wrapper<const std::array<Tensor, 4>>, 4> op_groups{
      std::cref(sc_phys_ops_a), std::cref(sc_phys_ops_b), std::cref(sc_phys_ops_c), std::cref(sc_phys_ops_d)};
  const std::array<std::string, 4> op_suffix{"a", "b", "c", "d"};

  const std::string file_name_postfix = "_fix_refer" + std::to_string(ref_bond[0]);
  size_t task_index = 0;
  for (const auto &task : tasks) {
    if (task.bonds.empty()) {
      ++task_index;
      continue;
    }
    for (size_t op_index = 0; op_index < op_groups.size(); ++op_index) {
      const size_t slot = task_index * op_groups.size() + op_index;
      if (mpi_size == 0) {
        continue;
      }
      if (rank == static_cast<int>(slot % static_cast<size_t>(mpi_size))) {
        auto measure_res = MeasureFourSiteOpGroup(
            mps, kMpsPath, op_groups[op_index].get(), ref_bond, task.bonds, f);
        DumpMeasuRes(measure_res, task.label_prefix + op_suffix[op_index] + file_name_postfix);
      }
    }
    ++task_index;
  }
}
