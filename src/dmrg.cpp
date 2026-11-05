// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2023-04-1
*
* Description: DMRG for t-J model with t_2 and J_2. Renormalized operator realization.
*/


#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include <algorithm>
#include <random>
#include <cassert>
#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"
#include "squarelattice.h"
#include "honeycomb_lattice.h"
#include "tJmodel.h"
#ifndef USE_REAL
#include "complex_perturb_init.h"
#endif // USE_REAL

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);
  CaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = ComputeTotalSites(params);
  tJModelParamters model_params(params);
  if (rank == 0) {
    model_params.Print();
  }
  clock_t startTime, endTime;
  startTime = clock();
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
  std::vector<size_t> input_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      input_D_set);

  double complex_init_eps = 0.0;
  unsigned long complex_init_seed = 12345;
  if (!ParserComplexInitArgs(argc, argv, complex_init_eps, complex_init_seed)) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  const bool use_complex_perturb = (complex_init_eps > 0.0);
#ifdef USE_REAL
  if (use_complex_perturb) {
    if (rank == kMPIMasterRank) {
      std::cerr << "[error] --complex-init-eps requires complex build (RealCode=OFF)." << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }
#endif

  size_t DMRG_time = input_D_set.size();
  std::vector<size_t> MaxLanczIterSet(DMRG_time);
  if (has_bond_dimension_parameter) {
    MaxLanczIterSet.back() = params.MaxLanczIter;
    if (DMRG_time > 1) {
      size_t MaxLanczIterSetSpace;
      MaxLanczIterSet[0] = 3;
      MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
      if (rank == kMPIMasterRank)
        std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
      for (size_t i = 1; i < DMRG_time - 1; i++) {
        MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
        if (rank == kMPIMasterRank)
          std::cout << MaxLanczIterSet[i] << ", ";
      }
      if (rank == kMPIMasterRank)
        std::cout << MaxLanczIterSet.back() << "]" << std::endl;
    } else {
      if (rank == kMPIMasterRank)
        std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
    }
  }

  qlmps::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);

  const bool is_honeycomb = IsHoneycombGeometry(params);
  std::vector<std::pair<size_t, size_t>> honeycomb_nn_links;

  if (!is_honeycomb && Ly == 2) {
    SquareLadder lattice(Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "SquareCylinder") {
    SquareCylinder lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "SquareOBC") {
    SquareOBC lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "SquareRotated") {
    SquareRotatedCylinder lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "SquareTorus") {
    SquareTorus lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "HoneycombXC") {
    if (std::abs(params.phi) > 1e-15 && rank == kMPIMasterRank) {
      std::cerr << "[warning] Honeycomb geometry ignores phi; set phi=0." << std::endl;
    }
    if (std::abs(params.phi) > 1e-15) {
      model_params.phi = 0.0;
    }
    HoneycombXCCylinder lattice(Ly, Lx);
    cout << "honeycomb lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
    for (const auto &link : lattice.nearest_neighbor_links) {
      honeycomb_nn_links.emplace_back(std::get<0>(link), std::get<1>(link));
    }
  } else if (params.Geometry == "HoneycombYC") {
    if (std::abs(params.phi) > 1e-15 && rank == kMPIMasterRank) {
      std::cerr << "[warning] Honeycomb geometry ignores phi; set phi=0." << std::endl;
    }
    if (std::abs(params.phi) > 1e-15) {
      model_params.phi = 0.0;
    }
    HoneycombYCCylinder lattice(Ly, Lx);
    cout << "honeycomb lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
    for (const auto &link : lattice.nearest_neighbor_links) {
      honeycomb_nn_links.emplace_back(std::get<0>(link), std::get<1>(link));
    }
  } else {
    if (rank == kMPIMasterRank) {
      std::cerr << "[error] Unsupported geometry '" << params.Geometry << "'" << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }
#if SYM_LEVEL == 0
  qlmps::tJOperators<TenElemT, QNT> ops;
  for (size_t i = 0; i < N; i++) {
    mpo_gen.AddTerm(-params.mu, ops.nf, i);
  }
  mpo_gen.AddTerm(params.BondSingletPairPerturbation, ops.bupa, 0, ops.bdna, 1);
  mpo_gen.AddTerm(-params.BondSingletPairPerturbation, ops.bdna, 0, ops.bupa, 1);
  if (rank == kMPIMasterRank && std::abs(params.BondSingletPairPerturbation) > 0) {
    std::cout << "Add singlet pair perturbation on bond (0,1)" << std::endl;
  }
#endif
  auto mat_repr_mpo = mpo_gen.GenMatReprMPO();
  using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
  FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(N);
  std::vector<size_t> A_sublattice, B_sublattice;

  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly;
    const size_t x = i / Ly;
    const bool a_sublattice = ((x + y) % 2 == 0); 

    if (a_sublattice) {
        A_sublattice.push_back(i);
        stat_labs[i] = 0;
    } else {
        B_sublattice.push_back(i);
        stat_labs[i] = 1;
    }   
  }

  if (params.Numhole > 1) {
    size_t nhalf = params.Numhole / 2;
    size_t stepA = A_sublattice.size() / nhalf;
    size_t stepB = B_sublattice.size() / nhalf;

    size_t offsetA = stepA / 2;
    size_t offsetB = 0;

    for (size_t k = 0; k < nhalf; ++k) {
        size_t idxA = (offsetA + k * stepA) % A_sublattice.size();
        size_t idxB = (offsetB + k * stepB) % B_sublattice.size();

        stat_labs[A_sublattice[idxA]] = 2;
        stat_labs[B_sublattice[idxB]] = 2;
      }

    if (2 * nhalf < params.Numhole) {
      size_t idxB = B_sublattice.size() - 1;
      stat_labs[B_sublattice[idxB]] = 2;
    }
  }
  else if (params.Numhole == 1){
    size_t idxB = B_sublattice.size() / 2;
    stat_labs[B_sublattice[idxB]] = 2;
  }

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif

  // --- Complex perturbation: pre-master safety checks (all ranks) ---
#ifndef USE_REAL
  std::vector<size_t> perturb_labs;
  size_t swap_s1 = 0, swap_s2 = 0;
  // Amplitudes computed in the pre-master block below and consumed by the
  // master-rank MPS construction later — declared at the outer guard scope.
  double amp_a_val = 0.0;
  double amp_b_re = 0.0;
  double amp_b_im = 0.0;
  if (use_complex_perturb) {
    if (!is_honeycomb) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] --complex-init-eps only supports honeycomb geometries." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    int mps_path_exists = 0;
    if (rank == kMPIMasterRank) {
      mps_path_exists = IsPathExist(kMpsPath) ? 1 : 0;
    }
    MPI_Bcast(&mps_path_exists, 1, MPI_INT, kMPIMasterRank, comm);
    if (mps_path_exists) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] --complex-init-eps requires no existing ./mps directory. "
                  << "Remove or rename it first." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    std::sort(honeycomb_nn_links.begin(), honeycomb_nn_links.end());
    honeycomb_nn_links.erase(
        std::unique(honeycomb_nn_links.begin(), honeycomb_nn_links.end()),
        honeycomb_nn_links.end());

    if (honeycomb_nn_links.empty()) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] Honeycomb NN link table is empty (unexpected for a valid lattice)." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    std::vector<std::pair<size_t, size_t>> valid_swap_links;
    for (const auto &lnk : honeycomb_nn_links) {
      assert(lnk.first < N && lnk.second < N && lnk.first != lnk.second);
      if (stat_labs[lnk.first] != 2 && stat_labs[lnk.second] != 2 &&
          stat_labs[lnk.first] != stat_labs[lnk.second]) {
        valid_swap_links.push_back(lnk);
      }
    }

    if (valid_swap_links.empty()) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] No valid NN swap link found for complex perturbation." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    std::mt19937 rng(complex_init_seed);
    std::uniform_int_distribution<size_t> dist(0, valid_swap_links.size() - 1);
    size_t chosen_idx = dist(rng);
    swap_s1 = valid_swap_links[chosen_idx].first;
    swap_s2 = valid_swap_links[chosen_idx].second;

    perturb_labs.assign(stat_labs.begin(), stat_labs.end());
    size_t orig_lab_s1 = perturb_labs[swap_s1];
    size_t orig_lab_s2 = perturb_labs[swap_s2];
    std::swap(perturb_labs[swap_s1], perturb_labs[swap_s2]);

    size_t count_before[3] = {0, 0, 0};
    size_t count_after[3] = {0, 0, 0};
    for (size_t idx = 0; idx < N; ++idx) {
      assert(stat_labs[idx] <= 2);
      assert(perturb_labs[idx] <= 2);
      count_before[stat_labs[idx]]++;
      count_after[perturb_labs[idx]]++;
    }
    if (count_before[0] != count_after[0] ||
        count_before[1] != count_after[1] ||
        count_before[2] != count_after[2]) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] Swap changed total label counts — QN conservation violated." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    double norm_factor = std::hypot(1.0, complex_init_eps);
    amp_a_val = 1.0 / norm_factor;
    amp_b_re = complex_init_eps / (std::sqrt(2.0) * norm_factor);
    amp_b_im = amp_b_re;
    if (!std::isfinite(amp_a_val) || !std::isfinite(amp_b_re)) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] eps=" << complex_init_eps
                  << " causes overflow in amplitude computation." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }

    if (rank == kMPIMasterRank) {
      std::cout << "[complex-init] swap: site " << swap_s1 << " (label " << orig_lab_s1
                << ") <-> site " << swap_s2 << " (label " << orig_lab_s2 << ")" << std::endl;
      std::cout << "[complex-init] after swap: site " << swap_s1 << " (label " << perturb_labs[swap_s1]
                << "), site " << swap_s2 << " (label " << perturb_labs[swap_s2] << ")" << std::endl;
      std::cout << "[complex-init] eps=" << complex_init_eps
                << ", seed=" << complex_init_seed
                << ", norm_factor=" << norm_factor
                << ", amp_a=" << amp_a_val
                << ", amp_b=(" << amp_b_re << ", " << amp_b_im << ")" << std::endl;
    }
  }
#endif // USE_REAL

  int init_error = 0;
  if (rank == kMPIMasterRank) {
#ifndef USE_REAL
    if (use_complex_perturb) {
      try {
        TenElemT amp_a(amp_a_val);
        TenElemT amp_b(amp_b_re, amp_b_im);
        std::vector<std::vector<size_t>> stat_labs_set = {
            std::vector<size_t>(stat_labs.begin(), stat_labs.end()),
            perturb_labs
        };
        qlmps::SuperpositionTwoProductInitMps(mps, stat_labs_set, amp_a, amp_b);
        cout << "Initial mps as complex superposition state." << endl;
        mps.Dump(kMpsPath, true);
        if (!IsPathExist(kMpsPath) || GetNumofMps() != N) {
          std::cerr << "[error] MPS dump verification failed: expected " << N
                    << " tensors in " << kMpsPath << ", found " << GetNumofMps()
                    << ". Manual cleanup may be needed." << std::endl;
          init_error = 1;
        }
      } catch (const std::exception &e) {
        std::cerr << "[error] Complex init construction failed: " << e.what() << std::endl;
        init_error = 1;
      } catch (...) {
        std::cerr << "[error] Complex init construction failed with unknown exception." << std::endl;
        init_error = 1;
      }
    } else
#endif // USE_REAL
    {
      if (IsPathExist(kMpsPath)) {
        if (N == GetNumofMps()) {
          cout << "The number of mps files is consistent with mps size." << endl;
          cout << "Directly use mps from files." << endl;
        } else {
          qlmps::DirectStateInitMps(mps, stat_labs);
          cout << "Initial mps as direct product state." << endl;
          mps.Dump(kMpsPath, true);
        }
      } else {
        qlmps::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." << endl;
        mps.Dump(kMpsPath, true);
      }
    }
  }

  if (use_complex_perturb) {
    MPI_Bcast(&init_error, 1, MPI_INT, kMPIMasterRank, comm);
    if (init_error != 0) {
      if (rank == kMPIMasterRank) {
        std::cerr << "[error] MPS initialization failed. All ranks exiting." << std::endl;
      }
      MPI_Finalize();
      return EXIT_FAILURE;
    }
  }
  MPI_Barrier(comm);

  double e0;
  if (!has_bond_dimension_parameter) {
    qlmps::FiniteVMPSSweepParams sweep_params(
        params.Sweeps,
        params.Dmin, params.Dmax, params.CutOff,
        qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
        params.noise
    );
    e0 = qlmps::FiniteDMRG(mps, mat_repr_mpo, sweep_params, comm);
  } else {
    for (size_t i = 0; i < DMRG_time; i++) {
      size_t D = input_D_set[i];
      if (rank == 1) {
        std::cout << "D_max = " << D << std::endl;
      }
      qlmps::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          D, D, params.CutOff,
          qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i]),
          params.noise
      );
      e0 = qlmps::FiniteDMRG(mps, mat_repr_mpo, sweep_params, comm);
    }
  }
  std::cout << "E0/site: " << e0 / N << std::endl;
  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  MPI_Finalize();
  return 0;
}
