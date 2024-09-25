#ifndef TJ_MODEL_H
#define TJ_MODEL_H

#include "qlmps/one_dim_tn/mpo/mpo.h"
#include "qlmps/one_dim_tn/mpo/mpogen/mpogen.h"
#include "gqdouble.h"
#include "squarelattice.h"
#include "params_case.h"
#include "operators.h"

inline void ConstructtJMPOGenerator(
    qlmps::MPOGenerator<TenElemT, U1U1QN> &mpo_gen,
    const SquareLattice &lattice,
    const tJModelParamters &model_params
) {
  const double t1 = model_params.t,
      t2 = model_params.t2,
      J = model_params.J,
      J2 = model_params.J2,
      phi = model_params.phi;
  OperatorInitial();
  if (std::fabs(phi) < 1e-15) {
    for (const Link &link: lattice.nearest_neighbor_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1, {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1, {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1, {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1, {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {nf, nf}, {site1, site2});
    }

    for (const Link &link: lattice.next_nearest_neighbor_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t2, {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2, {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2, {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2, {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J2, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2, {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2, {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J2 / 4, {nf, nf}, {site1, site2});
    }
  } else {
    using Complex = std::complex<double>;
    for (const Link &link: lattice.direct_nn_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1, {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1, {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1, {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1, {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {nf, nf}, {site1, site2});
    }

    for (const Link &link: lattice.wind_nn_links) {
      // y of site1 == 0, y of site2 == Ly - 1
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, phi)), {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, -phi)), {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, -phi)), {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, phi)), {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2 * std::exp(Complex(0, 2 * phi)), {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2 * std::exp(Complex(0, -2 * phi)), {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {nf, nf}, {site1, site2});
    }

    for (const Link &link: lattice.direct_nnn_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t2, {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2, {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2, {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2, {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J2, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2, {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2, {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J2 / 4, {nf, nf}, {site1, site2});
    }

    for (const Link &link: lattice.wind_nnn_links_1) {
      // y of site1 == Ly - 1, y of site2 == 0
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, -phi)), {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, phi)), {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, phi)), {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, -phi)), {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J2, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2 * std::exp(Complex(0, -2 * phi)), {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2 * std::exp(Complex(0, 2 * phi)), {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J2 / 4, {nf, nf}, {site1, site2});
    }

    for (const Link &link: lattice.wind_nnn_links_2) {
      // y of site1 == 0, y of site2 == Ly - 1
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, phi)), {bupc, bupa}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, -phi)), {bdnc, bdna}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, -phi)), {bupa, bupc}, {site1, site2}, {f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, phi)), {bdna, bdnc}, {site1, site2}, {f});
      mpo_gen.AddTerm(J2, {sz, sz}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2 * std::exp(Complex(0, 2 * phi)), {sp, sm}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2 * std::exp(Complex(0, -2 * phi)), {sm, sp}, {site1, site2});
      mpo_gen.AddTerm(-J2 / 4, {nf, nf}, {site1, site2});
    }
  }

}

#endif //TJ_MODEL_H