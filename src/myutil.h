#ifndef TJMODEL_SRC_MYUTIL_H
#define TJMODEL_SRC_MYUTIL_H
#include <stdlib.h>
#include <vector>

#include "params_case.h"

size_t GetNumofMps();
void Show(std::vector<size_t> v);
bool ParserBondDimension(int, char *[], std::vector<size_t>& );
size_t ComputeTotalSites(const CaseParams &params);
bool IsHoneycombGeometry(const CaseParams &params);


//below is useless for t-J model
std::vector<size_t> GenerateDirectStateLabel(const size_t ly,
                                             const size_t lx,
                                             const size_t num_hole,
                                             const size_t ky_int);
std::vector<std::vector<size_t>> GenAllOrderedMomentumPairs(const size_t ly);
std::vector<std::vector<size_t>> GenAllMomentumPairs(const size_t ly);
std::vector<std::vector<size_t>> GenAllMomentumPairs3(const size_t ly);


#endif //TJMODEL_SRC_MYUTIL_H