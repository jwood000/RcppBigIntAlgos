#include "ReduceMatrix.h"

void reduceMatrix(std::size_t nCols, std::size_t nRows, 
                  std::vector<uint8_t> &nullMat,
                  std::vector<std::size_t> &myCols) {

    std::size_t matSize = nullMat.size();
    std::size_t rowInd = 0;

    for (std::size_t j = 0; j < nCols; ++j) {
        std::vector<int> myOnes;

        for (std::size_t i = rowInd; i < matSize; i += nCols)
            if (nullMat[i + j])
                myOnes.push_back(i);

        if (!myOnes.empty()) {
            std::size_t firstRow = myOnes.front();

            if (firstRow != rowInd)
                for (std::size_t k = 0; k < nCols; ++k)
                    std::swap(nullMat[firstRow + k], nullMat[rowInd + k]);

            for (auto it = myOnes.begin() + 1; it != myOnes.end(); ++it)
                for (std::size_t k = 0; k < nCols; ++k)
                    nullMat[*it + k] = nullMat[*it + k] ^ nullMat[rowInd + k];

            rowInd += nCols;
        }
    }
    
    if (rowInd < matSize && rowInd != 0)
        nullMat.resize(rowInd);

    if (rowInd > 0) {
        std::size_t i = 0;
        std::size_t k = 0;

        while (i < rowInd) {
            bool allZero = true;

            for (std::size_t j = 0; j < nCols; j++) {
                if (nullMat[i + j]) {
                    allZero = false;
                    break;
                }
            }

            if (allZero) {
                rowInd -= nCols;
                nullMat.resize(rowInd);
                continue;
            }

            if (nullMat[i + k] != 1) {
                for (std::size_t j = 0; j < nCols; ++j) {
                    if (nullMat[i + j]) {
                        for (std::size_t m = 0; m < rowInd; m += nCols)
                            std::swap(nullMat[m + k], nullMat[m + j]);

                        std::swap(myCols[j], myCols[k]);
                        break;
                    }
                }
            }

            i += nCols;
            ++k;
        }
    }
}
