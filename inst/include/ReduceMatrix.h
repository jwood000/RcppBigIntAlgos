#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include <vector>
#include <algorithm>
#include <cstdint>

constexpr int unrollSize = 8;
constexpr std::uint8_t u8one = static_cast<std::uint8_t>(1u) << 0;

void ReduceMatrix(std::vector<std::uint8_t> &nullMat,
                  std::vector<std::size_t> &myCols,
                  int nCols, int nRows);

#endif
