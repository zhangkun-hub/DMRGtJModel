// SPDX-License-Identifier: LGPL-3.0-only
/*
* Author: Hao-Xin Wang <wanghx18@mails.tsinghua.edu.cn>
* Creation Date: 2023-04-1
*
* Description: Single processor DMRG for t-J model with t_2 and J_2
*/


#include "gqmps2/gqmps2.h"
#include "gqten/gqten.h"
#include <time.h>
#include <stdlib.h>
#include "gqdouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"
#include "squarelattice.h"
#include "tJmodel.h"

using namespace gqmps2;
using namespace gqten;
using namespace std;


int main(int argc, char *argv[]) {
CaseParams params(argv[1]);
size_t Lx = params.Lx, Ly = params.Ly;
size_t N= Lx*Ly;
tJModelParamters model_params(params);
model_params.Print();
clock_t startTime,endTime;
startTime = clock();
OperatorInitial();
const SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(N, pb_out);
gqmps2::MPO<Tensor> mpo(N);
const std::string kMpoPath = "mpo";
const std::string kMpoTenBaseName = "mpo_ten";
if (IsPathExist(kMpoPath)) {
  for(size_t i=0; i<mpo.size();i++){
      std::string filename = kMpoPath + "/" +
       kMpoTenBaseName + std::to_string(i) + "." + kGQTenFileSuffix;
       mpo.LoadTen(i,filename);
  }

  cout << "MPO loaded." << endl;
}else{
  cout << "No mpo directory. exiting" << std::endl;
  exit(0);
}


using FiniteMPST = gqmps2::FiniteMPS<TenElemT, U1U1QN>;
FiniteMPST mps(sites);

  std::vector<long unsigned int> stat_labs(N);
  size_t site_number_per_hole;
  if(params.Numhole>0){
    site_number_per_hole = N/params.Numhole;
  }
  else{
    site_number_per_hole = N+99;
  }

  int sz_label = 0;
  srand((unsigned)time(NULL));
  for (size_t i = 0; i < N; ++i) {
    if ( i%site_number_per_hole == site_number_per_hole-1){
      stat_labs[i] = 2;}
    else {
      stat_labs[i] = sz_label%2;
      sz_label++;
    }
  }

gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);
gqmps2::SweepParams sweep_params(
    params.Sweeps,
    params.Dmin, params.Dmax, params.CutOff,
    gqmps2::LanczosParams(params.LanczErr, params.MaxLanczIter)
);


if(IsPathExist(kMpsPath)){
  if(N==GetNumofMps()){
    cout << "The number of mps files is consistent with mps size." <<endl;
    cout << "Directly use mps from files." <<endl;
  }else{
    gqmps2::DirectStateInitMps(mps, stat_labs);
    cout << "Initial mps as direct product state." <<endl;
    mps.Dump(sweep_params.mps_path, true);
  }
}else{
  gqmps2::DirectStateInitMps(mps, stat_labs);
  cout << "Initial mps as direct product state." <<endl;
  mps.Dump(sweep_params.mps_path, true);
}

auto e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
std::cout << "E0/site: " << e0 / N << std::endl;
endTime = clock();
cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
return 0;
}


