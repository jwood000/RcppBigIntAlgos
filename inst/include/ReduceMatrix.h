#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include "GmpxxCopy.h"
#include <vector>
#include <numeric>
#include <cstdint>
#include <cmath>
#include <bitset>

constexpr unsigned long int oneThousand = 1000u;
constexpr size_t wordSize = sizeof(std::bitset<1>) * 8;

void ReduceMatrix(std::vector<std::bitset<wordSize>> &nullMat,
                  std::vector<std::size_t> &myCols,
                  std::size_t nCols, std::size_t nRows);

#endif
