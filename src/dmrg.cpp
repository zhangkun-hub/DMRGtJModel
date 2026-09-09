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
#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"
#include "squarelattice.h"
#include "honeycomb_lattice.h"
#include "tJmodel.h"

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

  std::vector<long unsigned int> stat_labs(N); // stat_lab[i]=0,1,2(0->up, 1->dn, 2->hole)
  std::vector<size_t> A_sublattice, B_sublattice; // 创建存储honeycomb子晶格A/B格点的容器

  // 给子晶格A/B容器填充相应的格点，并按照 Néel 序填充up/dn电子
  for (size_t i = 0; i < N; ++i) {
    const size_t y = i % Ly;
    const size_t x = i / Ly;
    const bool a_sublattice = ((x + y) % 2 == 0); 

    // 满足(x + y) % 2 == 0的格点属于A子晶格
    if (a_sublattice) {
        A_sublattice.push_back(i);
        stat_labs[i] = 0; // A子晶格填充 ↑
    } else {
        B_sublattice.push_back(i);
        stat_labs[i] = 1; // B子晶格填充 ↓ 
    }   
  }

  // 在A/B子晶格内等间距填充hole
  if (params.Numhole > 1) {
    // 子晶格内等间距步长
    size_t nhalf = params.Numhole / 2; // 子晶格填充hole数目
    size_t stepA = A_sublattice.size() / nhalf;
    size_t stepB = B_sublattice.size() / nhalf;

    // 引入错位 offset，使得子晶格之间的空穴也等间距分布
    size_t offsetA = stepA / 2;  // 半步错位
    size_t offsetB = 0;

    for (size_t k = 0; k < nhalf; ++k) {
        size_t idxA = (offsetA + k * stepA) % A_sublattice.size();
        size_t idxB = (offsetB + k * stepB) % B_sublattice.size();

        stat_labs[A_sublattice[idxA]] = 2;
        stat_labs[B_sublattice[idxB]] = 2;
      }

    // 处理奇数个hole情形
    if (2 * nhalf < params.Numhole) {
      size_t idxB = B_sublattice.size() - 1;
      stat_labs[B_sublattice[idxB]] = 2;
    }
  }
  // 处理1个hole情形
  else if (params.Numhole == 1){
    size_t idxB = B_sublattice.size() / 2;
    stat_labs[B_sublattice[idxB]] = 2;
  }

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif

  if (rank == kMPIMasterRank) {
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
  MPI_Barrier(comm);

  double e0;
  if (!has_bond_dimension_parameter) {
    qlmps::FiniteVMPSSweepParams sweep_params(
        params.Sweeps,
        params.Dmin, params.Dmax, params.CutOff,
        qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter),
        params.noise   // 👈 直接用 JSON 里的 noise
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
          params.noise   // 👈 直接用 JSON 里的 noise
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


