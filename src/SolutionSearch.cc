#include "SolutionSearch.h"

void solutionSearch(std::vector<uint8_t> mat, std::size_t matNRows,
                    std::size_t matNCols, mpz_t n, mpz_t *const mpzFacBase,
                    mpz_t *const testInterval, mpz_t *const factors) {
    
    std::vector<uint8_t> nullMat;
    
    for (std::size_t j = 0; j < matNCols; ++j) {
        std::size_t i = 0;
        
        while ((i < matNRows) && ((mat[i + j] % 2) == 0))
            i += matNCols;
        
        if (i < mat.size())
            for (std::size_t k = 0; k < mat.size(); k += matNCols)
                nullMat.push_back(mat[k + j] % 2);
    }
    
    const std::size_t nCols = matNRows;
    const std::size_t nRows = nullMat.size() / nCols;
    std::vector<std::size_t> myCols(nCols, 0);
    
    for (std::size_t i = 0; i < nCols; ++i)
        myCols[i] = i;
    
    reduceMatrix(nCols, nRows, nullMat, myCols);
    const std::size_t newNrow = nullMat.size() / nCols;
    
    std::vector<std::vector<std::size_t>> myList(nCols, std::vector<std::size_t>());
    std::vector<std::size_t> freeVariables;
    freeVariables.reserve(nCols);
    
    if (nCols > newNrow) {
        for (std::size_t i = newNrow; i < nCols; ++i)
            freeVariables.push_back(myCols[i]);
        
        for (const auto fV: freeVariables)
            myList[fV].push_back(fV);
    }
    
    if (newNrow > 0) {
        for (std::size_t i = newNrow; i > 0; --i) {
            std::vector<std::size_t> temp;
            
            for (std::size_t j = 0; j < nCols; ++j)
                if (nullMat[(i - 1) * nCols + j])
                    temp.push_back(j);
                
            if (temp.size() == 1) {
                for (std::size_t j = 0; j < newNrow; ++j)
                    nullMat[(j * nCols) + temp.front()] = 0;
                
                myList[myCols[i - 1]].clear();
                myList[myCols[i - 1]].push_back(0);
            } else {
                temp.clear();
                bool allBiggerNewNrow = true;
                
                for (std::size_t j = i; j < nCols; ++j) {
                    if (nullMat[(i - 1) * nCols + j]) {
                        temp.push_back(j);
                        
                        if (allBiggerNewNrow && j < newNrow)
                            allBiggerNewNrow = false;
                    }
                }
                
                if (allBiggerNewNrow) {
                    myList[myCols[i - 1]].clear();
                    
                    for (const auto t: temp)
                        for (std::size_t j = 0; j < myList[myCols[t]].size(); j++)
                            myList[myCols[i - 1]].push_back(myList[myCols[t]][j]);
                } else {
                    for (const auto t: temp)
                        myList[myCols[i - 1]] = outersect(myList[myCols[i - 1]], myList[myCols[t]]);
                }
            }
        }
    }
    
    const std::size_t lenFree = freeVariables.size();
    
    if (lenFree > 0) {
        mpz_t mpzTemp1, mpzTemp2, mpzMin, xMpz, yMpz;
        mpz_init(mpzTemp1); mpz_init(mpzTemp2);
        mpz_init(mpzMin); mpz_init(xMpz); mpz_init(yMpz);
        
        mpz_ui_pow_ui (mpzTemp1, 2, lenFree);
        mpz_sub_ui (mpzTemp1, mpzTemp1, 1);
        
        bool bExit = false;
        const std::size_t myLim = (mpz_cmp_ui(mpzTemp1, hundredThousand) > 0)
                                    ? hundredThousand : mpz_get_ui(mpzTemp1);
        
        for (std::size_t i = 1; i <= myLim; ++i) {
            if (bExit) {break;}
            
            std::vector<std::size_t> ansVec;
            std::vector<uint8_t> posVec(nCols, 0u);
            std::vector<uint8_t> posAns = myIntToBit(i, lenFree);
            
            for (std::size_t j = 0; j < myList.size(); ++j) {
                for (const auto it: myList[j]) {
                    for (std::size_t k = 0; k < lenFree; ++k) {
                        if (freeVariables[k] == it) {
                            posVec[j] ^= posAns[k];
                            break;
                        }
                    }
                }
                
                if (posVec[j])
                    ansVec.push_back(j);
            }
            
            if (!ansVec.empty()) {
                std::size_t myCheck = 0;
                std::vector<std::size_t> yExponents(matNCols, 0);
                
                for (std::size_t j = 0; j < matNCols; ++j) {
                    for (const auto it: ansVec)
                        yExponents[j] += mat[it * matNCols + j];
                    
                    myCheck += (yExponents[j] % 2);
                    yExponents[j] >>= 1;
                }
                
                if (myCheck == 0) {
                    yExponents.erase(yExponents.begin());
                    mpz_set_ui(xMpz, 1);
                    mpz_set_ui(yMpz, 1);
                    
                    for (const auto aV: ansVec) {
                        mpz_mul(xMpz, xMpz, testInterval[aV]);
                        mpz_mod(xMpz, xMpz, n);
                    }
                    
                    for (std::size_t j = 0; j < yExponents.size(); ++j) {
                        mpz_pow_ui(mpzTemp1, mpzFacBase[j], yExponents[j]);
                        mpz_mul(yMpz, yMpz, mpzTemp1);
                        mpz_mod(yMpz, yMpz, n);
                    }
                    
                    mpz_sub(mpzTemp1, xMpz, yMpz);
                    mpz_gcd(mpzTemp1, mpzTemp1, n);
                    mpz_add(mpzTemp2, xMpz, yMpz);
                    mpz_gcd(mpzTemp2, mpzTemp2, n);
                    
                    if (mpz_cmp(mpzTemp1, mpzTemp2) < 0) {
                        mpz_set(mpzMin, mpzTemp1);
                    } else {
                        mpz_set(mpzMin, mpzTemp2);
                    }
                    
                    if (mpz_cmp_ui(mpzMin, 1) > 0) {
                        if (mpz_cmp(mpzTemp1, mpzTemp2) < 0) {
                            mpz_set(factors[0], mpzTemp1);
                            mpz_set(factors[1], mpzTemp2);
                        } else {
                            mpz_set(factors[1], mpzTemp1);
                            mpz_set(factors[0], mpzTemp2);
                        }
                        
                        bExit = true;
                    }
                }
            }
        }
        
        mpz_clear(mpzTemp1); mpz_clear(mpzTemp2);
        mpz_clear(mpzMin); mpz_clear(xMpz); mpz_clear(yMpz);
    }
}
