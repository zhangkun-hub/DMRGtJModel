/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2023/04/03.
 *
 */


#ifndef TJMODEL_SRC_PARAMS_CASE_H
#define TJMODEL_SRC_PARAMS_CASE_H

#include "qlmps/case_params_parser.h"
using qlmps::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Geometry = ParseStr("Geometry");
    Ly = ParseInt("Ly");
    Lx = ParseInt("Lx");
    Numhole = ParseInt("NumHole");
    t = ParseDouble("t");
    t2 = ParseDouble("t2");
    J = ParseDouble("J");
    J2 = ParseDouble("J2");
    phi = ParseDouble("phi");
#if SYM_LEVEL == 0
    mu = ParseDouble("mu");
#endif
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    Threads = ParseInt("Threads");
    Perturbation = ParseDouble("Perturbation");
#if SYM_LEVEL == 0
    BondSingletPairPerturbation = ParseDouble("BondSingletPairPerturbation");
#endif
    wavelength = ParseDouble("wavelength");
    noise = ParseDoubleVec("noise");
  }

  std::string Geometry; // Cylinder, Torus, OBC, Rotated, Ladder (useless actually)
  size_t Ly;
  size_t Lx;
  size_t Numhole;
  double t;
  double t2;
  double J;
  double J2;
  double phi;     // twist boundary condition, Phys. Rev. B 107, 075127 Eq.(5)
#if SYM_LEVEL == 0
  double mu;
#endif
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  size_t Threads;
  double Perturbation;
  double BondSingletPairPerturbation;
  double wavelength;
  std::vector<double> noise;

};

struct tJModelParamters {
 public :
  tJModelParamters(const double t,
                   const double t2,
                   const double J,
                   const double J2,
                   const double phi
  ) : t(t), t2(t2), J(J), J2(J2), phi(phi) {}

  tJModelParamters(CaseParams &params) :
      tJModelParamters(params.t,
                       params.t2,
                       params.J,
                       params.J2,
                       params.phi
      ) {}

  inline void Print(void) {
    using std::cout;
    cout << " ****** tJV Model Parameter List ****** " << "\n";
    cout << "NN hopping t  = " << t << "\n";
    cout << "NNN hopping t2 = " << t2 << "\n";
    cout << "NN super-exchange J  = " << J << "\n";
    cout << "NNN super-exchangeJ2 = " << J2 << "\n";
    cout << "Twist angle phi = " << phi << "\n";
  }

  double t;
  double t2;
  double J;
  double J2;
  double phi;
};
#endif //TJMODEL_SRC_PARAMS_CASE_H
