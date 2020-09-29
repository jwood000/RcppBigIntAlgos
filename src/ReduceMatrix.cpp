#include "ReduceMatrix.h"

void ReduceMatrix(std::vector<std::uint8_t> &nullMat,
                  std::vector<std::size_t> &myCols,
                  int nCols, int nRows) {

    const int matSize = nullMat.size();
    int rowInd = 0;
    
    for (int j = 0; j < nCols; ++j) {
        std::vector<int> rows;
        
        for (int i = rowInd + j; i < matSize; i += nCols)
            if (nullMat[i])
                rows.push_back(i - j);
        
        if (!rows.empty()) {
            std::vector<int> cols;
            const int firstRow = rows.front();

            if (firstRow != rowInd)
                for (int k = j; k < nCols; ++k)
                    if (nullMat[firstRow + k] != nullMat[rowInd + k])
                        std::swap(nullMat[firstRow + k], nullMat[rowInd + k]);
            
            if (rows.size() > 1) {
                for (int k = j; k < nCols; ++k)
                    if (nullMat[rowInd + k])
                        cols.push_back(k);
                
                for (std::size_t i = 1, r = rows[1]; i < rows.size(); ++i, r = rows[i])
                    for (auto k: cols)
                        nullMat[r + k] ^= u8one;
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
            const auto it = std::find(nullMat.cbegin(),
                                      nullMat.cend(), u8one);
            
            if (it == nullMat.end()) {
                rowInd -= nCols;
            } else {
                if (!nullMat[i + k]) {
                    for (int j = k + 1; j < nCols; ++j) {
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
        }
        
        nullMat.resize(rowInd);
    }
}
