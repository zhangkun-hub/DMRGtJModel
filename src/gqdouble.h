#ifndef TJMODEL_SRC_GQDOUBLE_H
#define TJMODEL_SRC_GQDOUBLE_H
#include "boost/serialization/complex.hpp"
#include "qlten/qlten.h"

using qlten::QNCard;
using qlten::U1QNVal;
using qlten::TenIndexDirType;

using TenElemT = qlten::QLTEN_Complex;
using U1U1QN = qlten::special_qn::U1U1QN;
using Tensor = qlten::QLTensor<TenElemT, U1U1QN>;

using QNSctT = qlten::QNSector<U1U1QN>;
using IndexT = qlten::Index<U1U1QN>;

const auto qn0 = U1U1QN(
    {QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
);

const IndexT pb_out = IndexT({
                                 QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(1))}), 1),
                                 QNSctT(U1U1QN({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(-1))}), 1),
                                 QNSctT(U1U1QN({QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}), 1)},
                             TenIndexDirType::OUT
);
const auto pb_in = qlten::InverseIndex(pb_out);

#endif // TJMODEL_SRC_GQDOUBLE_H