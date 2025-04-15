/*
    measureSC.cpp
    for measure pair correlation function. memory optimized and parallel version.
    usage:
        mpirun -n 8 ./measureSC --start=start_site
*/
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"
#include "qlten/utility/timer.h"

using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
using qlmps::SiteVec;
using qlmps::MeasureFourSiteOpGroup;
using qlmps::kMpsPath;
using qlmps::DumpMeasuRes;
using qlten::Timer;

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
  size_t N = Lx * Ly;
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
  sc_phys_ops_c = {ops.id, ops.id, ops.bupc, ops.bdnc};
  sc_phys_ops_d = {ops.id, ops.id, ops.bdnc, ops.bupc};
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
  const size_t N = Lx * Ly;
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  FiniteMPST mps(sites);
  std::array<size_t, 2> ref_site{start, start + 1}; // y bond
  std::vector<std::array<size_t, 2>> target_x_bond;
  std::vector<std::array<size_t, 2>> target_y_bond;
  std::vector<size_t> Tx(Lx * Ly), Ty(Lx * Ly);
  for (size_t i = 0; i < Lx * Ly; ++i) {
    size_t y = i % Ly, x = i / Ly;
    Tx[i] = y + Ly * ((x + 1) % Lx);
    Ty[i] = (y + 1) % Ly + Ly * x;
  }

  target_y_bond.reserve(Ly * Lx);
  for (size_t x = 0; x < Lx; x++) {
    for (size_t y = 0; y < Ly; y++) {
      size_t target_site0 = x * Ly + y;
      size_t target_site1 = Ty[target_site0];
      if (target_site0 > start + 1 && target_site1 > start + 1) {
        if (params.Geometry == "OBC" && y < Ly - 1) {
          target_y_bond.push_back({target_site0, target_site1});
        } else if (params.Geometry == "Cylinder") {
          if (y == 0) {
            target_y_bond.push_back({target_site0, target_site1});
            target_y_bond.push_back({target_site0, target_site0 + Ly - 1});
          } else if (y < Ly - 1) {
            target_y_bond.push_back({target_site0, target_site1});
          }
        }
      }
    }
  }
  target_x_bond.reserve(Ly * Lx);
  for (size_t x = 0; x < Lx - 1; x++) {
    for (size_t y = 0; y < Ly; y++) {
      size_t target_site0 = x * Ly + y;
      size_t target_site1 = Tx[target_site0];
      if (target_site0 > start + 1 && target_site1 > start + 1) {
        target_x_bond.push_back({target_site0, target_site1});
      }
    }
  }

  std::string file_name_postfix = "_fix_refer" + std::to_string(ref_site[0]);
  if (rank == 0) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_a, ref_site, target_y_bond, f);
    DumpMeasuRes(measure_res, "scsyya" + file_name_postfix);
  }
  if (rank == 1 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_b, ref_site, target_y_bond, f);
    DumpMeasuRes(measure_res, "scsyyb" + file_name_postfix);
  }
  if (rank == 2 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_c, ref_site, target_y_bond, f);
    DumpMeasuRes(measure_res, "scsyyc" + file_name_postfix);
  }
  if (rank == 3 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_d, ref_site, target_y_bond, f);
    DumpMeasuRes(measure_res, "scsyyd" + file_name_postfix);
  }
  if (rank == 4 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_a, ref_site, target_x_bond, f);
    DumpMeasuRes(measure_res, "scsyxa" + file_name_postfix);
  }
  if (rank == 5 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_b, ref_site, target_x_bond, f);
    DumpMeasuRes(measure_res, "scsyxb" + file_name_postfix);
  }
  if (rank == 6 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_c, ref_site, target_x_bond, f);
    DumpMeasuRes(measure_res, "scsyxc" + file_name_postfix);
  }
  if (rank == 7 % mpi_size) {
    auto measure_res = MeasureFourSiteOpGroup(mps, kMpsPath, sc_phys_ops_d, ref_site, target_x_bond, f);
    DumpMeasuRes(measure_res, "scsyxd" + file_name_postfix);
  }
}
