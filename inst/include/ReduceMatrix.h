#ifndef REDUCE_MATRIX_H
#define REDUCE_MATRIX_H

#include <vector>
#include <algorithm>
#include <cstdint>

void reduceMatrix(std::size_t nCols, std::size_t nRows, 
                  std::vector<uint8_t> &nullMat,
                  std::vector<std::size_t> &myCols);

#endif
