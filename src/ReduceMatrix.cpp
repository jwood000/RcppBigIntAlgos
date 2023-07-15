#include "ReduceMatrix.h"

void ReduceMatrix(std::vector<std::bitset<wordSize>> &nullMat,
                  std::vector<std::size_t> &myCols,
                  std::size_t nCols, std::size_t nRows) {

    const std::size_t matSize = nullMat.size();
    const std::size_t adjustedCols = (nCols + wordSize - 1) / wordSize;
    std::size_t rowInd = 0;

    for (std::size_t j = 0; j < nCols; ++j) {
        std::vector<std::size_t> rows;

        for (std::size_t d = j / wordSize, i = rowInd + d,
             m = j % wordSize; i < matSize; i += adjustedCols) {

            if (nullMat[i].test(m)) {
                rows.push_back(i - d);
            }
        }

        if (!rows.empty()) {
            std::vector<std::size_t> cols;
            auto&& firstRow = rows.front();

            if (firstRow != rowInd) {
                for (std::size_t k = j / wordSize; k < adjustedCols; ++k) {
                    std::swap(nullMat[firstRow + k], nullMat[rowInd + k]);
                }
            }

            if (rows.size() > 1) {
                for (std::size_t k = j / wordSize; k < adjustedCols; ++k) {
                    if (nullMat[rowInd + k].any()) {
                        cols.push_back(k);
                    }
                }

                for (std::size_t i = 1; i < rows.size(); ++i) {
                    auto&& r = rows[i];

                    for (auto k: cols) {
                        nullMat[r + k] ^= nullMat[rowInd + k];
                    }
                }
            }

            rowInd += adjustedCols;
        }
    }

    if (rowInd < matSize && rowInd != 0) {
        nullMat.resize(rowInd);
    }

    if (rowInd > 0) {
        for (std::size_t i = 0, k = 0; i < rowInd; ) {
            bool allZero = true;

            for (std::size_t j = 0; j < adjustedCols; ++j) {
                if (nullMat[i + j].any()) {
                    allZero = false;
                    break;
                }
            }

            if (allZero) {
                rowInd -= nCols;
            } else {
                if (!nullMat[i + k / wordSize].test(k % wordSize)) {
                    for (std::size_t d2 = (k + 1) / wordSize;
                         d2 < adjustedCols; ++d2) {
                        if (nullMat[i + d2].any()) {
                            const std::size_t d1 = k / wordSize;
                            const std::size_t col1 = k % wordSize;
                            std::size_t col2 = 0;

                            while (!nullMat[i + d2].test(col2))
                                ++col2;

                            for (std::size_t m = 0;
                                 m < rowInd; m += adjustedCols) {

                                if (nullMat[m + d1].test(col1) !=
                                    nullMat[m + d2].test(col2)) {
                                    nullMat[m + d1].flip(col1);
                                    nullMat[m + d2].flip(col2);
                                }
                            }

                            std::swap(myCols[d2 * wordSize + col2], myCols[k]);
                            break;
                        }
                    }
                }

                i += adjustedCols;
                ++k;
            }
        }

        nullMat.resize(rowInd);
    }
}
