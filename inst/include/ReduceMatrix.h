#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include <vector>
#include <algorithm>
#include <cstdint>

void reduceMatrix(std::vector<std::uint8_t> &nullMat,
                  std::vector<std::size_t> &myCols,
                  int nCols, int nRows);

#endif
