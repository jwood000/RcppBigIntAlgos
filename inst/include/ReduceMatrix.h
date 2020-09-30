#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include "SolutionUtils.h"
#include <vector>
#include <algorithm>
#include <cstdint>

void ReduceMatrix(std::vector<std::bitset<wordSize>> &nullMat,
                  std::vector<std::size_t> &myCols,
                  std::size_t nCols, std::size_t nRows);

#endif
