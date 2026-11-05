#ifndef TJMODEL_SRC_COMPLEX_PERTURB_INIT_H
#define TJMODEL_SRC_COMPLEX_PERTURB_INIT_H

#ifndef USE_REAL

#include "qlten/qlten.h"
#include "qlmps/qlmps.h"
#include <cassert>
#include <iostream>

namespace qlmps {

template<typename TenElemT, typename QNT>
void SuperpositionTwoProductInitMps(
    FiniteMPS<TenElemT, QNT> &mps,
    const std::vector<std::vector<size_t>> &stat_labs_set,
    const TenElemT &amp_a,
    const TenElemT &amp_b
) {
  using TenT = qlten::QLTensor<TenElemT, QNT>;
  using IndexT = qlten::Index<QNT>;
  using QNSctT = qlten::QNSector<QNT>;

  assert(stat_labs_set.size() == 2);
  auto N = mps.size();
  assert(N == stat_labs_set[0].size());
  assert(N == stat_labs_set[1].size());
  for (size_t i = 0; i < N; ++i) { mps.dealloc(i); }
  auto sites_info = mps.GetSitesInfo();
  auto pb_out_set = sites_info.sites;

  auto total_qn_a = pb_out_set[0].GetQNSctFromActualCoor(stat_labs_set[0][0]).GetQn();
  for (size_t i = 1; i < N; ++i) {
    total_qn_a += pb_out_set[i].GetQNSctFromActualCoor(stat_labs_set[0][i]).GetQn();
  }

  auto total_qn_b = pb_out_set[0].GetQNSctFromActualCoor(stat_labs_set[1][0]).GetQn();
  for (size_t i = 1; i < N; ++i) {
    total_qn_b += pb_out_set[i].GetQNSctFromActualCoor(stat_labs_set[1][i]).GetQn();
  }
  if (!(total_qn_a == total_qn_b)) {
    std::cerr << "[fatal] Two product states have different total QN "
              << "(swap did not preserve quantum numbers)." << std::endl;
    assert(false);
  }

  const auto &div = total_qn_a;
  auto zero_qn = div - div;

  IndexT lvb, rvb;
  std::vector<QNSctT> rvb_qnscts;

  // --- Head tensor (site 0) ---
  std::vector<QNT> head_expected_qns(2);
  for (size_t j = 0; j < 2; ++j) {
    auto stat_lab = stat_labs_set[j][0];
    auto rvb_qn = div - pb_out_set[0].GetQNSctFromActualCoor(stat_lab).GetQn();
    head_expected_qns[j] = rvb_qn;
    rvb_qnscts.push_back(QNSctT(rvb_qn, 1));
  }
  lvb = IndexT({QNSctT(zero_qn, 1)}, qlten::TenIndexDirType::IN);
  rvb = IndexT(rvb_qnscts, qlten::TenIndexDirType::OUT);
  rvb_qnscts.clear();

  for (size_t j = 0; j < 2; ++j) {
    auto actual_qn = rvb.GetQNSctFromActualCoor(j).GetQn();
    if (!(actual_qn == head_expected_qns[j])) {
      std::cerr << "[fatal] Head rvb: actual coordinate " << j
                << " maps to unexpected QN. Index construction may reorder sectors."
                << std::endl;
      assert(false);
    }
  }

  mps(0) = new TenT({lvb, pb_out_set[0], rvb});
  (mps[0])({0, stat_labs_set[0][0], 0}) = amp_a;
  (mps[0])({0, stat_labs_set[1][0], 1}) = amp_b;

  // --- Middle tensors (sites 1 .. N-2) ---
  for (size_t i = 1; i < N - 1; ++i) {
    lvb = qlten::InverseIndex(rvb);

    std::vector<QNT> mid_expected_qns(2);
    for (size_t j = 0; j < 2; ++j) {
      auto stat_lab = stat_labs_set[j][i];
      auto lvb_qn = lvb.GetQNSctFromActualCoor(j).GetQn();
      auto rvb_qn = lvb_qn -
          pb_out_set[i].GetQNSctFromActualCoor(stat_lab).GetQn();
      mid_expected_qns[j] = rvb_qn;
      rvb_qnscts.push_back(QNSctT(rvb_qn, 1));
    }
    rvb = IndexT(rvb_qnscts, qlten::TenIndexDirType::OUT);
    rvb_qnscts.clear();

    for (size_t j = 0; j < 2; ++j) {
      auto actual_qn = rvb.GetQNSctFromActualCoor(j).GetQn();
      if (!(actual_qn == mid_expected_qns[j])) {
        std::cerr << "[fatal] Site " << i << " rvb: actual coordinate " << j
                  << " maps to unexpected QN." << std::endl;
        assert(false);
      }
    }

    mps(i) = new TenT({lvb, pb_out_set[i], rvb});
    for (size_t j = 0; j < 2; ++j) {
      (mps[i])({j, stat_labs_set[j][i], j}) = TenElemT(1.0);
    }
  }

  // --- Tail tensor (site N-1) ---
  lvb = qlten::InverseIndex(rvb);
  rvb = IndexT({QNSctT(zero_qn, 1)}, qlten::TenIndexDirType::OUT);
  mps(N - 1) = new TenT({lvb, pb_out_set[N - 1], rvb});
  for (size_t j = 0; j < 2; ++j) {
    (mps[N - 1])({j, stat_labs_set[j][N - 1], 0}) = TenElemT(1.0);
  }

  // Centralize to site 0 and normalize
  mps.Centralize(0);
  mps[0].Normalize();
}

} // namespace qlmps

#endif // USE_REAL
#endif // TJMODEL_SRC_COMPLEX_PERTURB_INIT_H
