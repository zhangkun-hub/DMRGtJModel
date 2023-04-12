/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2023/04/03.
 *
 */


#ifndef TJMODEL_SRC_PARAMS_CASE_H
#define TJMODEL_SRC_PARAMS_CASE_H

#include "gqmps2/case_params_parser.h"
using gqmps2::CaseParamsParserBasic;

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
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    Threads = ParseInt("Threads");
    Perturbation=ParseDouble("Perturbation");
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
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  size_t Threads;
  double Perturbation;
  double wavelength;
  std::vector<double> noise;

};



struct tJModelParamters {
 public :
  tJModelParamters(const double t,
                   const double t2,
                   const double J,
                   const double J2
  ) : t(t), t2(t2), J(J), J2(J2) {}

  tJModelParamters(CaseParams& params) :
      tJModelParamters(params.t,
                        params.t2,
                        params.J,
                        params.J2
      ) {}

  inline void Print(void){
    using std::cout;
    cout << " ****** tJV Model Parameter List ****** " << "\n";
    cout << "t  = " << t << "\n";
    cout << "t2 = " << t2 << "\n";
    cout << "J  = " << J  << "\n";
    cout << "J2 = " << J2 << "\n";
  }

  double t;
  double t2;
  double J;
  double J2;
};
#endif //TJMODEL_SRC_PARAMS_CASE_H
