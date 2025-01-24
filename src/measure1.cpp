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
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"

using namespace qlmps;
using namespace qlten;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = Lx * params.Ly;

  clock_t startTime, endTime;
  startTime = clock();

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif
  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);

  Timer one_site_timer("measure one site operators");
  MeasureOneSiteOp(mps, kMpsPath, {nf, sz}, {"nf", "sz"});
  one_site_timer.PrintElapsed();

  Timer two_site_timer("measure two site operators");
  size_t ref_site = 0;
  auto szsz_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, sz, sz, ref_site);
  DumpMeasuRes(szsz_corr, "sz0sz");
  auto spsm_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, sp, sm, ref_site);
  DumpMeasuRes(spsm_corr, "sp0sm");
  auto smsp_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, sm, sp, ref_site);
  DumpMeasuRes(smsp_corr, "sm0sp");
  auto nn_corr = MeasureTwoSiteOpGroup(mps, kMpsPath, nf, nf, ref_site);
  DumpMeasuRes(nn_corr, "nf0nf");
  auto single_particle_corr1 = MeasureTwoSiteOpGroup(mps, kMpsPath, bdna, bdnc, ref_site, f);
  auto single_particle_corr2 = MeasureTwoSiteOpGroup(mps, kMpsPath, bdnc, bdna, ref_site, f);
  auto single_particle_corr3 = MeasureTwoSiteOpGroup(mps, kMpsPath, bupa, bupc, ref_site, f);
  auto single_particle_corr4 = MeasureTwoSiteOpGroup(mps, kMpsPath, bupc, bupa, ref_site, f);
  DumpMeasuRes(single_particle_corr1, "single_particle_corr1");
  DumpMeasuRes(single_particle_corr2, "single_particle_corr2");
  DumpMeasuRes(single_particle_corr3, "single_particle_corr3");
  DumpMeasuRes(single_particle_corr4, "single_particle_corr4");

  two_site_timer.PrintElapsed();
  endTime = clock();
  std::cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
  return 0;
}
