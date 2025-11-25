// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2023-04-1
*
* Description: Measurement of the one-point function for the MPS of ground state of t-J model
* Usage: ./measure1 params.json
*/

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include "tJ_hilbert_space.h"
#include "params_case.h"
#include "myutil.h"
#include "squarelattice.h"
#include "my_measure.h"

using namespace qlmps;
using namespace qlten;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  size_t N = ComputeTotalSites(params);
  size_t ref_site = 0;
  if (argc > 2) {
    ref_site = std::atoi(argv[2]);
  }
  clock_t startTime, endTime;
  startTime = clock();

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif
  qlmps::tJOperators<TenElemT, QNT> ops;
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, QNT>;
  FiniteMPST mps(sites);

  Timer one_site_timer("measure one site operators");
  MeasureOneSiteOp(mps, kMpsPath, {ops.nf, ops.sz}, {"nf", "sz"});
  one_site_timer.PrintElapsed();

  auto f = ops.f;
  auto bdna = ops.bdna, bupa = ops.bupa, bdnc = ops.bdnc, bupc = ops.bupc;
  Timer two_site_timer("measure two site operators");
  auto szsz_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, ops.sz, ops.sz, ref_site);
  DumpMeasuRes(szsz_corr, "sz" + std::to_string(ref_site) + "sz");
  auto spsm_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, ops.sp, ops.sm, ref_site);
  DumpMeasuRes(spsm_corr, "sp" + std::to_string(ref_site) + "sm");
  auto smsp_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, ops.sm, ops.sp, ref_site);
  DumpMeasuRes(smsp_corr, "sm" + std::to_string(ref_site) + "sp");
  auto nn_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, ops.nf, ops.nf, ref_site);
  DumpMeasuRes(nn_corr, "nf" + std::to_string(ref_site) + "nf");
  auto single_particle_corr1 = MeasureTwoSiteOpGroup(mps, kMpsPath, bdna, bdnc, ref_site, f);
  auto single_particle_corr2 = MeasureTwoSiteOpGroup(mps, kMpsPath, bdnc, bdna, ref_site, f);
  auto single_particle_corr3 = MeasureTwoSiteOpGroup(mps, kMpsPath, bupa, bupc, ref_site, f);
  auto single_particle_corr4 = MeasureTwoSiteOpGroup(mps, kMpsPath, bupc, bupa, ref_site, f);
  DumpMeasuRes(single_particle_corr1, "single_particle_corr1");
  DumpMeasuRes(single_particle_corr2, "single_particle_corr2");
  DumpMeasuRes(single_particle_corr3, "single_particle_corr3");
  DumpMeasuRes(single_particle_corr4, "single_particle_corr4");

  if (params.Geometry == "SquareOBC") {
    SquareOBC lattice(params.Ly, params.Lx);
    auto convert_links = [](const std::vector<std::tuple<size_t, size_t>> &links) {
      std::vector<std::pair<size_t, size_t>> res;
      res.reserve(links.size());
      for (const auto &l : links) res.emplace_back(std::get<0>(l), std::get<1>(l));
      return res;
    };

    auto nn_links = convert_links(lattice.nearest_neighbor_links);
    auto nnn_links = convert_links(lattice.next_nearest_neighbor_links);

    // Measure NN
    auto nn_down = MeasureBonds(mps, kMpsPath, nn_links, bdnc, bdna, f);
    DumpMeasuRes(nn_down, "hopping_nn_down");

    auto nn_up = MeasureBonds(mps, kMpsPath, nn_links, bupc, bupa, f);
    DumpMeasuRes(nn_up, "hopping_nn_up");

    // Measure NNN
    auto nnn_down = MeasureBonds(mps, kMpsPath, nnn_links, bdnc, bdna, f);
    DumpMeasuRes(nnn_down, "hopping_nnn_down");

    auto nnn_up = MeasureBonds(mps, kMpsPath, nnn_links, bupc, bupa, f);
    DumpMeasuRes(nnn_up, "hopping_nnn_up");
  }

  two_site_timer.PrintElapsed();
  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  return 0;
}
