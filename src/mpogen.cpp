#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include <cstdlib>
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
  CaseParams params(argv[1]);
  size_t Lx = params.Lx, Ly = params.Ly;
  size_t N = ComputeTotalSites(params);
  tJModelParamters model_params(params);
  model_params.Print();
  clock_t startTime, endTime;
  startTime = clock();
  const SiteVec<TenElemT, QNT> sites = tj_site.GenUniformSites<TenElemT>(N);
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
    if (std::abs(params.phi) > 1e-15) {
      std::cerr << "[warning] Honeycomb geometry ignores phi; set phi=0." << std::endl;
      model_params.phi = 0.0;
    }
    HoneycombXCCylinder lattice(Ly, Lx);
    cout << "honeycomb lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else if (params.Geometry == "HoneycombYC") {
    if (std::abs(params.phi) > 1e-15) {
      std::cerr << "[warning] Honeycomb geometry ignores phi; set phi=0." << std::endl;
      model_params.phi = 0.0;
    }
    HoneycombYCCylinder lattice(Ly, Lx);
    cout << "honeycomb lattice construct" << std::endl;
    ConstructtJMPOGenerator(mpo_gen, lattice, model_params);
  } else {
    std::cerr << "[error] Unsupported geometry '" << params.Geometry << "'" << std::endl;
    return EXIT_FAILURE;
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


