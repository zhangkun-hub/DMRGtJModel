// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2023-04-1
*
* Description: Measurement of the one-point function for the MPS of ground state of t-J model
* Usage: mpirun -n 3 ./measureS params.json
*/

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "tJ_hilbert_space.h"
#include "tJ_operators.h"
#include "params_case.h"
#include "myutil.h"
#include "my_measure.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

//forward declaration
template<typename TenElemT, typename QNT>
MPO<QLTensor<TenElemT, QNT>> MpoSquare(const MPO<QLTensor<TenElemT, QNT>> &);

template<typename TenElemT, typename QNT>
MPO<QLTensor<TenElemT, QNT>> MpoProduct(const MPO<QLTensor<TenElemT, QNT>> &,
                                        const MPO<QLTensor<TenElemT, QNT>> &);

template<typename TenElemT, typename QNT>
TenElemT ExpectationValue(FiniteMPS<TenElemT, QNT>,
                          const MPO<QLTensor<TenElemT, QNT>> &,
                          const std::string);
int main(int argc, char *argv[]) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, mpi_size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpi_size);

  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t Ly = params.Ly;
  size_t N = Lx * Ly;

  clock_t startTime, endTime;
  startTime = clock();

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif

  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);
  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);
  qlmps::MPOGenerator<TenElemT, U1U1QN> Sp(sites, qn0);
  qlmps::MPOGenerator<TenElemT, U1U1QN> Sm(sites, qn0);
  qlmps::MPOGenerator<TenElemT, U1U1QN> Sz(sites, qn0);
  tJOperators ops;
  for (size_t i = 0; i < N; i++) {
    Sz.AddTerm(1.0, ops.sz, i);
    Sp.AddTerm(1.0, ops.sp, i);
    Sm.AddTerm(1.0, ops.sm, i);
  }
  auto Sz_mpo = Sz.Gen();
  auto Sp_mpo = Sp.Gen();
  auto Sm_mpo = Sm.Gen();

  auto Sz_square = MpoSquare(Sz_mpo);
  auto SpSm = MpoProduct(Sp_mpo, Sm_mpo);
  auto SmSp = MpoProduct(Sm_mpo, Sp_mpo);

  TenElemT res_sz2, res_spsm, res_smsp;
  if (rank == 0) {
    res_sz2 = ExpectationValue(mps, Sz_square, kMpsPath);
    std::cout << "Sz square = " << res_sz2 << std::endl;
    hp_numeric::MPI_Recv(&res_spsm, 1, 1, 15, comm);
    std::cout << "Sp*Sm = " << res_spsm << std::endl;
    hp_numeric::MPI_Recv(&res_smsp, 1, 2, 16, comm);
    std::cout << "Sm*Sp = " << res_smsp << std::endl;
    auto S_square = res_sz2 + (res_spsm + res_smsp) / 2.0;
    std::cout << "S^2 = " << S_square << std::endl;
  } else if (rank == 1) {
    res_spsm = ExpectationValue(mps, SpSm, kMpsPath);
    hp_numeric::MPI_Send(&res_spsm, 1, 0, 15, comm);
  } else if (rank == 2) {
    res_smsp = ExpectationValue(mps, SmSp, kMpsPath);
    hp_numeric::MPI_Send(&res_smsp, 1, 0, 16, comm);
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  MPI_Finalize();
  return 0;
}

template<typename TenElemT, typename QNT>
MPO<QLTensor<TenElemT, QNT>> MpoSquare(const MPO<QLTensor<TenElemT, QNT>> &mpo) {

/*
 *     2          2         2            2
 *     |          |         |            |
 *0----|----3 0-------3 0--------3  0---------3
 *     |          |        |             |
 *     1          1        1             1
 *     2          2         2            2
 *     |          |         |            |
 *0----|----3 0-------3 0--------3  0---------3
 *     |          |        |             |
 *     1          1        1             1
 *
 */
  using Tensor = QLTensor<TenElemT, QNT>;
  MPO<QLTensor<TenElemT, QNT>> res = mpo;
  size_t N = mpo.size();
  for (size_t i = 0; i < N; i++) {
    res[i] = Tensor();
    Contract(&mpo[i], &mpo[i], {{2}, {1}}, &res[i]);
    res[i].FuseIndex(0, 3);
    res[i].FuseIndex(2, 4);
    res[i].Transpose({1, 2, 3, 0});
    assert(res[i].GetIndexes()[1] == mpo[i].GetIndexes()[1]);
    assert(res[i].GetIndexes()[2] == mpo[i].GetIndexes()[2]);
  }
  return res;
}

template<typename TenElemT, typename QNT>
TenElemT ExpectationValue(FiniteMPS<TenElemT, QNT> mps,
                          const MPO<QLTensor<TenElemT, QNT>> &mpo,
                          const std::string mps_path) {
  using TenT = QLTensor<TenElemT, QNT>;
  auto N = mps.size();

  TenT renv;
  //Write a trivial right environment tensor to disk
  mps.LoadTen(N - 1, GenMPSTenName(mps_path, N - 1));
  auto mps_trivial_index = mps.back().GetIndexes()[2];
  auto mpo_trivial_index_inv = InverseIndex(mpo.back().GetIndexes()[3]);
  auto mps_trivial_index_inv = InverseIndex(mps_trivial_index);
  renv = TenT({mps_trivial_index, mpo_trivial_index_inv, mps_trivial_index_inv});
  renv({0, 0, 0}) = 1;

  //bulk right environment tensors
  for (size_t i = 1; i <= N; ++i) {
    std::cout << "right block length = " << i << std::endl;
    if (i > 1) { mps.LoadTen(N - i, GenMPSTenName(mps_path, N - i)); }
    renv = std::move(UpdateSiteRenvs(renv, mps[N - i], mpo[N - i]));
    mps.dealloc(N - i);
  }

  //Write a trivial left environment tensor to disk
  TenT lenv;
  mps.LoadTen(0, GenMPSTenName(mps_path, 0));
  mps_trivial_index = mps.front().GetIndexes()[0];
  mpo_trivial_index_inv = InverseIndex(mpo.front().GetIndexes()[0]);
  mps_trivial_index_inv = InverseIndex(mps_trivial_index);
  lenv = TenT({mps_trivial_index, mpo_trivial_index_inv, mps_trivial_index_inv});
  lenv({0, 0, 0}) = 1;
  mps.dealloc(0);

  TenT res_ten;
  Contract(&lenv, &renv, {{0, 1, 2}, {0, 1, 2}}, &res_ten);

  assert(mps.empty());
  return res_ten();
}

template<typename TenElemT, typename QNT>
MPO<QLTensor<TenElemT, QNT>> MpoProduct(const MPO<QLTensor<TenElemT, QNT>> &mpo1,
                                        const MPO<QLTensor<TenElemT, QNT>> &mpo2) {

/*
 *     2          2         2            2
 *     |          |         |            |
 *0----|----3 0-------3 0--------3  0---------3 (mpo2)
 *     |          |        |             |
 *     1          1        1             1
 *     2          2         2            2
 *     |          |         |            |
 *0----|----3 0-------3 0--------3  0---------3 (mpo1)
 *     |          |        |             |
 *     1          1        1             1
 *
 */
  using Tensor = QLTensor<TenElemT, QNT>;
  MPO<QLTensor<TenElemT, QNT>> res = mpo1;
  size_t N = mpo1.size();
  for (size_t i = 0; i < N; i++) {
    res[i] = Tensor();
    Contract(&mpo1[i], &mpo2[i], {{2}, {1}}, &res[i]);
    res[i].FuseIndex(0, 3);
    res[i].FuseIndex(2, 4);
    res[i].Transpose({1, 2, 3, 0});
    assert(res[i].GetIndexes()[1] == mpo1[i].GetIndexes()[1]);
    assert(res[i].GetIndexes()[2] == mpo2[i].GetIndexes()[2]);
  }
  MPI_Finalize();
  return res;
}
