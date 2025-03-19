#ifndef TJMODEL_SRC_GQDOUBLE_H
#define TJMODEL_SRC_GQDOUBLE_H
#include "qlten/qlten.h"

using qlten::QNCard;
using qlten::U1QNVal;
using qlten::TenIndexDirType;

#ifdef USE_REAL
using TenElemT = qlten::QLTEN_Double;
#else
using TenElemT = qlten::QLTEN_Complex;
#endif
#if SYM_LEVEL == 1
using QNT = qlten::special_qn::U1U1QN;
#else
using QNT = qlten::special_qn::U1QN; // only Sz
#endif
using Tensor = qlten::QLTensor<TenElemT, QNT>;

using QNSctT = qlten::QNSector<QNT>;
using IndexT = qlten::Index<QNT>;

#if SYM_LEVEL == 1
const auto qn0 = QNT(
    {QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
);

const IndexT pb_out = IndexT({
                                 QNSctT(QNT({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(1))}), 1),
                                 QNSctT(QNT({QNCard("N", U1QNVal(1)), QNCard("Sz", U1QNVal(-1))}), 1),
                                 QNSctT(QNT({QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}), 1)},
                             TenIndexDirType::OUT
);
#else
const auto qn0 = QNT(0);
const IndexT pb_out = IndexT({
                                 QNSctT(QNT({QNCard("Sz", U1QNVal(1))}), 1),
                                 QNSctT(QNT({QNCard("Sz", U1QNVal(-1))}), 1),
                                 QNSctT(QNT({QNCard("Sz", U1QNVal(0))}), 1)},
                             TenIndexDirType::OUT
);

#endif
const auto pb_in = qlten::InverseIndex(pb_out);

#endif // TJMODEL_SRC_GQDOUBLE_H