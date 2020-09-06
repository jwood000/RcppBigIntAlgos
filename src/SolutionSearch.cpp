#include "SolutionSearch.h"
#include "Cpp14MakeUnique.h"
#include <RcppThread.h>

void ProcessFreeMat(const std::vector<std::uint8_t> &nullMat,
                    const std::vector<std::size_t> &myCols,
                    std::vector<std::uint8_t> &freeMat,
                    std::size_t newNrow, std::size_t nCols) {
    
    const std::size_t freeMatSize = freeMat.size();

    for (int i = newNrow - 1; i >= 0; --i) {
        std::vector<int> temp;
        
        for (int j = i + 1; j < static_cast<int>(nCols); ++j)
            if (nullMat[i * nCols + j])
                temp.push_back(j);

        if (!temp.empty()) {
            if (temp.front() >= newNrow) {
                for (const auto t: temp)
                    for (std::size_t j = 0; j < freeMatSize; j += nCols)
                        if (freeMat[myCols[t] + j])
                            freeMat[myCols[i] + j] = 1u;
            } else {
                for (const auto t: temp)
                    for (std::size_t j = 0; j < freeMatSize; j += nCols)
                        freeMat[myCols[i] + j] = freeMat[myCols[i] + j] ^ freeMat[myCols[t] + j];
            }
        }
    }
}

bool GetSolution(const std::vector<std::uint8_t> &freeMat,
                 const std::vector<std::uint8_t> &mat,
                 const std::vector<std::size_t> &freeVariables,
                 const std::vector<mpz_class> &mpzFacBase,
                 const std::vector<mpz_class> &testInterval,
                 mpz_t *const factors, mpz_t myNum, std::size_t nCols, 
                 std::size_t matNCols, std::size_t ind, 
                 std::size_t lenFree, std::size_t threadInd) {
    
    mpz_t mpzTemp1, mpzTemp2, mpzMin, xMpz, yMpz;
    mpz_init(mpzTemp1); mpz_init(mpzTemp2);
    mpz_init(yMpz); mpz_init(mpzMin); mpz_init(xMpz);
    
    std::vector<std::size_t> ansVec;
    std::vector<std::uint8_t> posVec(nCols, 0u);
    std::vector<std::uint8_t> posAns = MyIntToBit(ind, lenFree);
    
    bool bSuccess = false;
    const int lastUnroll = nCols - (nCols % unrollSize);
    
    for (int i = 0; i < static_cast<int>(freeVariables.size()); ++i) {
        for (int k = 0, j = i * nCols; k < lastUnroll; k += unrollSize, j += unrollSize) {
            if (freeMat[j]) {posVec[k] = posVec[k] ^ posAns[i];}
            if (freeMat[j + 1]) {posVec[k + 1] = posVec[k + 1] ^ posAns[i];}
            if (freeMat[j + 2]) {posVec[k + 2] = posVec[k + 2] ^ posAns[i];}
            if (freeMat[j + 3]) {posVec[k + 3] = posVec[k + 3] ^ posAns[i];}
            if (freeMat[j + 4]) {posVec[k + 4] = posVec[k + 4] ^ posAns[i];}
            if (freeMat[j + 5]) {posVec[k + 5] = posVec[k + 5] ^ posAns[i];}
            if (freeMat[j + 6]) {posVec[k + 6] = posVec[k + 6] ^ posAns[i];}
            if (freeMat[j + 7]) {posVec[k + 7] = posVec[k + 7] ^ posAns[i];}
        }
        
        for (int k = lastUnroll, j = i * nCols + lastUnroll; k < static_cast<int>(nCols); ++k, ++j)
            if (freeMat[j]) {posVec[k] = posVec[k] ^ posAns[i];}
    }
    
    for (std::size_t k = 0; k < nCols; ++k)
        if (posVec[k])
            ansVec.push_back(k);

    if (!ansVec.empty()) {
        std::size_t myCheck = 0;
        std::vector<std::size_t> yExponents(matNCols, 0);
        
        for (std::size_t j = 0; j < matNCols; ++j) {
            for (const auto aV: ansVec)
                yExponents[j] += mat[aV * matNCols + j];
            
            myCheck += (yExponents[j] % 2u);
            yExponents[j] >>= 1;
        }
        
        if (myCheck == 0) {
            yExponents.erase(yExponents.begin());
            mpz_set_ui(xMpz, 1);
            mpz_set_ui(yMpz, 1);
            
            for (const auto aV: ansVec) {
                mpz_mul(xMpz, xMpz, testInterval[aV].get_mpz_t());
                mpz_mod(xMpz, xMpz, myNum);
            }
            
            for (std::size_t j = 0; j < yExponents.size(); ++j) {
                mpz_pow_ui(mpzTemp1, mpzFacBase[j].get_mpz_t(), yExponents[j]);
                mpz_mul(yMpz, yMpz, mpzTemp1);
                mpz_mod(yMpz, yMpz, myNum);
            }
            
            mpz_sub(mpzTemp1, xMpz, yMpz);
            mpz_gcd(mpzTemp1, mpzTemp1, myNum);
            mpz_add(mpzTemp2, xMpz, yMpz);
            mpz_gcd(mpzTemp2, mpzTemp2, myNum);
            
            if (mpz_cmp(mpzTemp1, mpzTemp2) < 0) {
                mpz_set(mpzMin, mpzTemp1);
            } else {
                mpz_set(mpzMin, mpzTemp2);
            }
            
            if (mpz_cmp_ui(mpzMin, 1) > 0) {
                if (mpz_cmp(mpzTemp1, mpzTemp2) < 0) {
                    mpz_set(factors[threadInd * 2], mpzTemp1);
                    mpz_set(factors[threadInd * 2 + 1], mpzTemp2);
                } else {
                    mpz_set(factors[threadInd * 2 + 1], mpzTemp1);
                    mpz_set(factors[threadInd * 2], mpzTemp2);
                }
                
                bSuccess = true;
            }
        }
    }
    
    mpz_clear(mpzTemp1); mpz_clear(mpzTemp2);
    mpz_clear(mpzMin); mpz_clear(xMpz); mpz_clear(yMpz);
    return bSuccess;
}

void SolutionSearch(const std::vector<std::uint8_t> &mat, std::size_t matNRows,
                    std::size_t matNCols, mpz_t myNum, const std::vector<mpz_class> &mpzFacBase,
                    const std::vector<mpz_class> &testInterval, mpz_t *const factors,
                    std::size_t nThreads) {
    
    const std::size_t matSize = mat.size();
    const std::size_t nCols = matNRows;
    
    std::vector<std::uint8_t> nullMat;
    nullMat.reserve(matSize);
    
    for (std::size_t j = 0; j < matNCols; ++j) {
        std::size_t i = 0;
        
        while ((i < matSize) && ((mat[i + j] % 2u) == 0))
            i += matNCols;
        
        if (i < matSize)
            for (std::size_t k = 0; k < matSize; k += matNCols)
                nullMat.push_back(mat[k + j] % 2u);
    }
    
    const std::size_t nRows = nullMat.size() / nCols;
    
    std::vector<std::size_t> myCols(nCols, 0);
    std::iota(myCols.begin(), myCols.end(), 0);
    
    ReduceMatrix(nullMat, myCols, static_cast<int>(nCols), static_cast<int>(nRows));
    const std::size_t newNrow = nullMat.size() / nCols;
    std::vector<std::size_t> freeVariables;
    
    if (nCols > newNrow && newNrow > 0) {
        for (std::size_t i = newNrow; i < nCols; ++i)
            freeVariables.push_back(myCols[i]);
        
        std::sort(freeVariables.begin(), freeVariables.end());
        const std::size_t myMin = freeVariables.front();
        
        const std::size_t lenFree = freeVariables.size();
        std::vector<std::uint8_t> freeMat(lenFree * nCols, 0u);
        
        std::transform(freeVariables.begin(), freeVariables.end(),
                       freeVariables.begin(), [myMin](std::size_t f) {return f - myMin;});
        
        // freeVariables isn't guranteed to be contiguous. That is,
        // we would have freeVarabiables = {10, 14, 15, 17}. This means
        // that lenFree = 4, and since the dimensions of freeMat is
        // based off of lenFree and not the range of (fV), we must
        // take care not to access memory we don't own.
        
        for (std::size_t i = 0; i < freeVariables.size(); ++i)
            freeMat[i * nCols + freeVariables[i] + myMin] = 1u;
        
        ProcessFreeMat(nullMat, myCols, freeMat, newNrow, nCols);
        mpz_t mpzTemp1;
        mpz_init(mpzTemp1);

        mpz_ui_pow_ui(mpzTemp1, 2, lenFree);
        mpz_sub_ui(mpzTemp1, mpzTemp1, 1);
        const std::size_t myLim = (mpz_cmp_ui(mpzTemp1, oneThousand) > 0)
                                    ? oneThousand : mpz_get_ui(mpzTemp1);

        if (nThreads > 1) {
            auto vecFactors = FromCpp14::make_unique<mpz_t[]>(nThreads * 2);
            std::vector<std::future<bool>> myFutures(nThreads);
            std::vector<bool> vecSuccess(nThreads);
            bool bSuccess = false;

            for (std::size_t i = 1; i < myLim && !bSuccess;) {
                RcppThread::ThreadPool pool(nThreads);

                for (std::size_t j = 0; j < nThreads; ++j, ++i) {
                    myFutures[j] = pool.pushReturn(std::cref(GetSolution), std::cref(freeMat), std::cref(mat),
                                                   std::cref(freeVariables), mpzFacBase, testInterval,
                                                   vecFactors.get(), myNum, nCols, matNCols, i, lenFree, j);
                }

                pool.join();

                for (std::size_t j = 0; j < nThreads; ++j)
                    vecSuccess[j] = myFutures[j].get();

                bSuccess = std::any_of(vecSuccess.begin(), vecSuccess.end(), [](bool v) {return v;});
            }

            for (std::size_t j = 0; j < nThreads; ++j) {
                if (vecSuccess[j]) {
                    mpz_set(factors[0], vecFactors[j * 2]);
                    mpz_set(factors[1], vecFactors[j * 2 + 1]);
                    break;
                }
            }

            for (std::size_t j = 0; j < (nThreads * 2); ++j)
                mpz_clear(vecFactors[j]);

            vecFactors.reset();
        } else {
            bool bSuccess = false;

            for (std::size_t i = 1; i < myLim && !bSuccess; ++i) {
                bSuccess = GetSolution(freeMat, mat, freeVariables, mpzFacBase,
                                       testInterval, factors, myNum, nCols,
                                       matNCols, i, lenFree, 0);
            }
        }
    }
}
