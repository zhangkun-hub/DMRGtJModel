/*
    measureSC.cpp
    for measure pair correlation function. memory optimized and parallel version.
    usage:
        mpirun -n 8 ./measureSC --start=start_site
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

/**
 * @brief Factory function to create a HoneycombLattice object based on the geometry parameter.
 * 
 * @param params Case parameters containing geometry information.
 * @return std::unique_ptr<HoneycombLattice> Pointer to the created lattice object, or nullptr if not a honeycomb geometry.
 */
std::unique_ptr<HoneycombLattice> BuildHoneycombLattice(const CaseParams &params) {
  if (params.Geometry == "HoneycombXC") {
    return std::make_unique<HoneycombXCCylinder>(params.Ly, params.Lx);
  }
  if (params.Geometry == "HoneycombYC") {
    return std::make_unique<HoneycombYCCylinder>(params.Ly, params.Lx);
  }
  return nullptr;
}

/**
 * @brief Checks if two sites are nearest neighbors in the given lattice.
 * 
 * @param lattice The honeycomb lattice structure.
 * @param site0 First site index.
 * @param site1 Second site index.
 * @return true If site0 and site1 are nearest neighbors.
 * @return false Otherwise.
 */
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

/**
 * @brief Determines the reference bond for square lattice geometry based on a requested start site.
 * 
 * @param requested_start The user-requested starting site for the reference bond.
 * @param total_sites Total number of sites in the system.
 * @param Ly Width of the system (cylinder circumference).
 * @param params Case parameters.
 * @return Bond The determined reference bond {site0, site1}.
 */
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

/**
 * @brief Determines the reference bond for honeycomb lattice geometry.
 * 
 * Validates and potentially adjusts the requested start site to ensure the bond is a valid
 * nearest-neighbor link compatible with the specific honeycomb geometry (XC or YC).
 * 
 * @param requested_start The user-requested starting site.
 * @param total_sites Total number of sites.
 * @param Ly Lattice width.
 * @param params Case parameters.
 * @param lattice The honeycomb lattice structure.
 * @return Bond The valid reference bond.
 */
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

/**
 * @brief Main dispatcher to determine reference bond based on geometry.
 */
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

/**
 * @brief Builds target bond lists for honeycomb lattice.
 * 
 * Groups bonds by direction (A, B, C) and filters them based on the reference threshold.
 * 
 * @param lattice Honeycomb lattice structure.
 * @param ref_threshold Site index threshold; bonds involving sites <= this are ignored (unless ignore_threshold is true).
 * @param ignore_threshold If true, includes all bonds regardless of position relative to reference.
 * @return std::array<BondList, 3> Array of bond lists for directions A, B, and C.
 */
std::array<BondList, 3> BuildHoneycombTargetBonds(const HoneycombLattice &lattice,
                                                  size_t ref_threshold,
                                                  bool ignore_threshold = false) {
  std::array<BondList, 3> targets;
  std::unordered_set<std::pair<size_t, size_t>, PairHash> seen;
  seen.reserve(lattice.nearest_neighbor_links.size());
  for (const auto &link : lattice.nearest_neighbor_links) {
    size_t site0 = std::get<0>(link);
    size_t site1 = std::get<1>(link);
    if (site0 > site1) {
      std::swap(site0, site1);
    }
    if (!ignore_threshold && (site0 <= ref_threshold || site1 <= ref_threshold)) {
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

/**
 * @brief Generates all target bonds to be measured.
 * 
 * Supports both Square and Honeycomb geometries.
 * 
 * @param params Case parameters.
 * @param ref_threshold Threshold for filtering bonds "behind" the reference bond.
 * @param honeycomb_lattice Pointer to HoneycombLattice if applicable.
 * @param ignore_threshold Whether to ignore the reference threshold (e.g. for local pairing measurement).
 * @return std::vector<MeasureTask> List of measurement tasks grouped by bond type/direction.
 */
std::vector<MeasureTask> GenerateTargetBonds(const CaseParams &params,
                                             size_t ref_threshold,
                                             const HoneycombLattice* honeycomb_lattice,
                                             bool ignore_threshold = false) {
  const size_t Lx = params.Lx;
  const size_t Ly = params.Ly;
  std::vector<MeasureTask> tasks;
  tasks.reserve(4);

  if (honeycomb_lattice) {
    auto targets = BuildHoneycombTargetBonds(*honeycomb_lattice, ref_threshold, ignore_threshold);
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
        if (!ignore_threshold && (site0 <= ref_threshold || site1 <= ref_threshold)) {
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
            if (!ignore_threshold && wrap_site <= ref_threshold) {
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
        if (!ignore_threshold && (site0 <= ref_threshold || site1 <= ref_threshold)) {
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
  return tasks;
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
            size_t &start) {
  int nOptionIndex = 1;

  std::string arguement1 = "--start=";
  bool start_argument_has(false);
  while (nOptionIndex < argc) {
    if (strncmp(argv[nOptionIndex], arguement1.c_str(), arguement1.size()) == 0) {
      std::string para_string = &argv[nOptionIndex][arguement1.size()];
      start = atoi(para_string.c_str());
      start_argument_has = true;
    }
    nOptionIndex++;
  }

  if (!start_argument_has) {
    std::cout << "Note: no start argument, set it by default (1)." << std::endl;
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
  bool start_argument_has = Parser(argc, argv, start);

  CaseParams params(argv[1]);

  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = ComputeTotalSites(params);
  if (GetNumofMps() != N) {
    std::cout << "The number of mps files are inconsistent with mps size!" << std::endl;
    exit(1);
  }

  if (!start_argument_has) {
    start = 1; // default start site for pairing measurement
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

  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  FiniteMPST mps(sites);
  auto honeycomb_lattice = BuildHoneycombLattice(params);

  // Generate all bonds, ignoring reference threshold
  auto local_tasks = GenerateTargetBonds(params, 0, honeycomb_lattice.get(), true);

  struct MeasureOp {
    std::string suffix;
    Tensor op1;
    Tensor op2;
  };

  // Corresponding to old logic:
  // a: bupa, bdna
  // b: bdna, bupa
  // c: bupc, bupc (forbidden by spin U1 symmetry, should be 0)
  // d: bdnc, bdnc (forbidden by spin U1 symmetry, should be 0)
  // We only measure a and b.
  std::vector<MeasureOp> meas_ops = {
    {"a", ops.bupa, ops.bdna},
    {"b", ops.bdna, ops.bupa}
  };

  int task_counter = 0;
  for (const auto& task : local_tasks) {
    if (task.bonds.empty()) continue;
    
    // Convert bonds to pair format
    std::vector<std::pair<size_t, size_t>> bond_pairs;
    bond_pairs.reserve(task.bonds.size());
    for(const auto& b : task.bonds) {
      bond_pairs.push_back({b[0], b[1]});
    }
    
    for (const auto& op_pair : meas_ops) {
      // Use MPI to distribute these independent measurements
      // Each (Task + OpPair) is a unit of work
      if (task_counter++ % mpi_size == rank) {
        auto res = MeasureBonds(mps, kMpsPath, bond_pairs, op_pair.op1, op_pair.op2, ops.f);
        // Assuming "local" suffix to distinguish from 4-point SC correlations
        // The old code reused PerformMeasureSC with ref_bond=0, producing files with "_fix_refer0"
        // We can match that or use something clearer.
        // Let's stick to what the user might expect or something descriptive.
        // Old filename format: task.label_prefix + op_suffix + "_fix_refer" + ref_bond[0]
        // Here ref_bond[0] would be 0 effectively.
        DumpMeasuRes(res, task.label_prefix + op_pair.suffix + "_local_pairing");
      }
    }
  }
#endif

  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  MPI_Finalize();
  return 0;
}

/**
 * @brief Executes the superconducting correlation measurement.
 * 
 * Calculates the reference bond and threshold, generates target bonds, 
 * and then distributes the 4-point correlation measurements across MPI ranks.
 * 
 * @param rank Current MPI rank.
 * @param mpi_size Total number of MPI processes.
 * @param start Requested start site for reference bond.
 * @param params Case parameters.
 * @param sc_phys_ops_a Operator group A.
 * @param sc_phys_ops_b Operator group B.
 * @param sc_phys_ops_c Operator group C.
 * @param sc_phys_ops_d Operator group D.
 * @param f Fermionic parity operator.
 */
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

  auto tasks = GenerateTargetBonds(params, ref_threshold, honeycomb_lattice.get(), false);

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
