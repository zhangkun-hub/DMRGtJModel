#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include <stdlib.h>
#include "gqdouble.h"
#include "tJ_operators.h"
#include "params_case.h"
#include "myutil.h"
#include "squarelattice.h"
#include "tJmodel.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = Lx * Ly;
  tJModelParamters model_params(params);
  model_params.Print();
  clock_t startTime, endTime;
  startTime = clock();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_out);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  if (Ly == 2) {
    SquareLadder lattice(Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "Cylinder") {
    SquareCylinder lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "OBC") {
    SquareOBC lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "Rotated") {
    SquareRotatedCylinder lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "Torus") {
    SquareTorus lattice(Ly, Lx);
    cout << "lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  }

#ifndef USE_GPU
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);
#endif
  qlmps::MPO<Tensor> mpo = mpo_gen.Gen();
  cout << "MPO generated." << endl;

  if (!IsPathExist(kMpoPath)) {
    CreatPath(kMpoPath);
  }
  for (size_t i = 0; i < mpo.size(); i++) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i)
        + "." + kQLTenFileSuffix;
    mpo.DumpTen(i, filename);
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}


