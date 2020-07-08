#include "ReduceMatrix.h"

constexpr int unrollSize = 8;

void reduceMatrix(std::vector<std::uint8_t> &nullMat,
                  std::vector<std::size_t> &myCols,
                  int nCols, int nRows) {

    int matSize = nullMat.size();
    int rowInd = 0;
    
    for (int j = 0; j < nCols; ++j) {
        std::vector<int> rows;

        for (int i = rowInd; i < matSize; i += nCols)
            if (nullMat[i + j])
                rows.push_back(i);
        
        if (!rows.empty()) {
            std::vector<int> cols;
            const int firstRow = rows.front();

            if (firstRow != rowInd)
                for (int k = 0; k < nCols; ++k)
                    if (nullMat[firstRow + k] != nullMat[rowInd + k])
                        std::swap(nullMat[firstRow + k], nullMat[rowInd + k]);
            
            for (int k = 0; k < nCols; ++k)
                if (nullMat[rowInd + k])
                    cols.push_back(k);
            
            const int lastUnroll = cols.size() - (cols.size() % unrollSize);
            
            for (int i = 1, rowSize = rows.size(), colSize = cols.size(); i < rowSize; ++i) {
                for (int k = 0; k < lastUnroll; k += unrollSize) {
                    nullMat[rows[i] + cols[k]] = nullMat[rows[i] + cols[k]] ^ 1u;
                    nullMat[rows[i] + cols[k + 1]] = nullMat[rows[i] + cols[k + 1]] ^ 1u;
                    nullMat[rows[i] + cols[k + 2]] = nullMat[rows[i] + cols[k + 2]] ^ 1u;
                    nullMat[rows[i] + cols[k + 3]] = nullMat[rows[i] + cols[k + 3]] ^ 1u;
                    nullMat[rows[i] + cols[k + 4]] = nullMat[rows[i] + cols[k + 4]] ^ 1u;
                    nullMat[rows[i] + cols[k + 5]] = nullMat[rows[i] + cols[k + 5]] ^ 1u;
                    nullMat[rows[i] + cols[k + 6]] = nullMat[rows[i] + cols[k + 6]] ^ 1u;
                    nullMat[rows[i] + cols[k + 7]] = nullMat[rows[i] + cols[k + 7]] ^ 1u;
                }
                
                for (int k = lastUnroll; k < colSize; ++k)
                    nullMat[rows[i] + cols[k]] = nullMat[rows[i] + cols[k]] ^ 1u;
            }
            
            rowInd += nCols;
        }
    }
    
    if (rowInd < matSize && rowInd != 0)
        nullMat.resize(rowInd);

    if (rowInd > 0) {
        int i = 0;
        int k = 0;

        while (i < rowInd) {
            bool allZero = true;
            
            for (int j = 0; j < nCols; ++j) {
                if (nullMat[i + j]) {
                    allZero = false;
                    break;
                }
            }

            if (allZero) {
                rowInd -= nCols;
                continue;
            }

            if (nullMat[i + k] != 1) {
                for (int j = 0; j < nCols; ++j) {
                    if (nullMat[i + j]) {
                        for (int m = 0; m < rowInd; m += nCols)
                            if (nullMat[m + k] != nullMat[m + j])
                                std::swap(nullMat[m + k], nullMat[m + j]);

                        std::swap(myCols[j], myCols[k]);
                        break;
                    }
                }
            }

            i += nCols;
            ++k;
        }
        
        nullMat.resize(rowInd);
    }
}
