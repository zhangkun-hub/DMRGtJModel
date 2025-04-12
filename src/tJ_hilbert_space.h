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
const auto qn0 = QNT(
    {QNCard("N", U1QNVal(0)), QNCard("Sz", U1QNVal(0))}
);
#else
using QNT = qlten::special_qn::U1QN; // only Sz
const auto qn0 = QNT(0);
#endif
using Tensor = qlten::QLTensor<TenElemT, QNT>;

using QNSctT = qlten::QNSector<QNT>;
using IndexT = qlten::Index<QNT>;

qlmps::sites::tJSite<QNT> tj_site;
using tJOperators = qlmps::tJOperators<TenElemT, QNT>;

#endif // TJMODEL_SRC_GQDOUBLE_H