#ifndef TJ_MODEL_H
#define TJ_MODEL_H

#include <algorithm>
#include <type_traits>

#include "qlmps/one_dim_tn/mpo/mpo.h"
#include "qlmps/one_dim_tn/mpo/mpogen/mpogen.h"
#include "tJ_hilbert_space.h"
#include "lattice_base.h"
#include "squarelattice.h"
#include "params_case.h"

template <typename LatticeT>
inline void ConstructtJMPOGenerator(
    qlmps::MPOGenerator<TenElemT, QNT> &mpo_gen,
    const LatticeT &lattice,
    const tJModelParamters &model_params
) {
  static_assert(std::is_base_of<LatticeBase, LatticeT>::value,
                "ConstructtJMPOGenerator expects a lattice derived from LatticeBase");
  using Link = LatticeBase::Link;
  const double t1 = model_params.t,
      t2 = model_params.t2,
      J = model_params.J,
      J2 = model_params.J2,
      phi = model_params.phi;
  qlmps::tJOperators<TenElemT, QNT> ops;
  if (std::fabs(phi) < 1e-15) {
    for (const Link &link : lattice.nearest_neighbor_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      if (site2 < site1) {
        std::swap(site1, site2);
      }
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1, {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {ops.nf, ops.nf}, {site1, site2});
    }

    for (const Link &link : lattice.next_nearest_neighbor_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      if (site2 < site1) {
        std::swap(site1, site2);
      }
      assert(site1 < site2);
      mpo_gen.AddTerm(-t2, {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t2, {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t2, {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t2, {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J2, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2, {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2, {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J2 / 4, {ops.nf, ops.nf}, {site1, site2});
    }
  } else {
#ifdef USE_REAL
    std::cout << "phi : " << phi << " is non-zero. Please use complex code! (-DRealCode=OFF when compile)" << std::endl;
    exit(1);
#else
    using Complex = std::complex<double>;
    for (const Link &link : lattice.direct_nn_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1, {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {ops.nf, ops.nf}, {site1, site2});
    }

    for (const Link &link : lattice.wind_nn_links) {
      // y of site1 == 0, y of site2 == Ly - 1
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, phi)), {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, -phi)), {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, -phi)), {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, phi)), {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2 * std::exp(Complex(0, 2 * phi)), {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2 * std::exp(Complex(0, -2 * phi)), {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {ops.nf, ops.nf}, {site1, site2});
    }

    for (const Link &link : lattice.direct_nnn_links) {
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1, {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1, {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2, {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {ops.nf, ops.nf}, {site1, site2});
    }

    for (const Link &link : lattice.wind_nnn_links_1) {
      // y of site1 == Ly - 1, y of site2 == 0
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, -phi)), {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, phi)), {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, phi)), {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t2 * std::exp(Complex(0, -phi)), {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J2, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2 * std::exp(Complex(0, -2 * phi)), {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J2 / 2 * std::exp(Complex(0, 2 * phi)), {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J2 / 4, {ops.nf, ops.nf}, {site1, site2});
    }

    for (const Link &link : lattice.wind_nnn_links_2) {
      // y of site1 == 0, y of site2 == Ly - 1
      size_t site1 = std::get<0>(link);
      size_t site2 = std::get<1>(link);
      assert(site1 < site2);
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, phi)), {ops.bupc, ops.bupa}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, -phi)), {ops.bdnc, ops.bdna}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, -phi)), {ops.bupa, ops.bupc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(-t1 * std::exp(Complex(0, phi)), {ops.bdna, ops.bdnc}, {site1, site2}, {ops.f});
      mpo_gen.AddTerm(J, {ops.sz, ops.sz}, {site1, site2});
      mpo_gen.AddTerm(J / 2 * std::exp(Complex(0, 2 * phi)), {ops.sp, ops.sm}, {site1, site2});
      mpo_gen.AddTerm(J / 2 * std::exp(Complex(0, -2 * phi)), {ops.sm, ops.sp}, {site1, site2});
      mpo_gen.AddTerm(-J / 4, {ops.nf, ops.nf}, {site1, site2});
    }
#endif
  }
}

#endif //TJ_MODEL_H