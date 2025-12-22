/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2023/04/03.
 *
 */


#ifndef TJMODEL_SRC_PARAMS_CASE_H
#define TJMODEL_SRC_PARAMS_CASE_H

#include "qlmps/case_params_parser.h"
#include <iostream>

using qlmps::CaseParamsParserBasic;

namespace {

inline std::string NormalizeGeometry(const std::string &raw_geometry) {
  auto warn = [](const std::string &message) {
    std::cerr << "[geometry-warning] " << message << std::endl;
  };

  if (raw_geometry == "SquareCylinder" || raw_geometry == "SquareOBC" ||
      raw_geometry == "SquareRotated" || raw_geometry == "SquareTorus" ||
      raw_geometry == "HoneycombXC" || raw_geometry == "HoneycombYC") {
    return raw_geometry;
  }

  if (raw_geometry == "Cylinder") {
    warn("Legacy geometry 'Cylinder' mapped to 'SquareCylinder'");
    return "SquareCylinder";
  }

  if (raw_geometry == "OBC") {
    warn("Legacy geometry 'OBC' mapped to 'SquareOBC'");
    return "SquareOBC";
  }

  if (raw_geometry == "Rotated") {
    warn("Legacy geometry 'Rotated' mapped to 'SquareRotated'");
    return "SquareRotated";
  }

  if (raw_geometry == "Torus") {
    warn("Legacy geometry 'Torus' mapped to 'SquareTorus'");
    return "SquareTorus";
  }

  return raw_geometry;
}

} // namespace

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    std::string raw_geometry = ParseStr("Geometry");
    Geometry = NormalizeGeometry(raw_geometry);
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

  std::string Geometry; // SquareCylinder, SquareOBC, SquareRotated, SquareTorus, HoneycombXC, HoneycombYC
  size_t Ly;
  size_t Lx;
  size_t Numhole;
  double t;
  double t2;
  double J;
  double J2;
  /**
   * @brief Twist angle for flux insertion / twisted boundary condition.
   *
   * Currently implemented for square geometries that populate winding link lists
   * (`direct_nn_links`, `wind_nn_links`, `direct_nnn_links`, `wind_nnn_links_*`).
   *
   * For honeycomb geometries (`HoneycombXC`, `HoneycombYC`), `phi` is currently
   * ignored at runtime (a warning is emitted and `phi` is treated as zero).
   *
   * Reference: Phys. Rev. B 107, 075127 Eq.(5).
   */
  double phi;
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
#if SYM_LEVEL == 0
  double BondSingletPairPerturbation;
#endif
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
