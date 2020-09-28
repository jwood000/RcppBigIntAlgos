#include "SolutionSearch.h"
#include "StatsUtils.h"
#include <random>

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
            if (temp.front() >= static_cast<int>(newNrow)) {
                for (const auto t: temp)
                    for (std::size_t j = 0; j < freeMatSize; j += nCols)
                        if (freeMat[myCols[t] + j])
                            freeMat[myCols[i] + j] = 1u;
            } else {
                for (const auto t: temp)
                    for (std::size_t j = 0; j < freeMatSize; j += nCols)
                        freeMat[myCols[i] + j] ^= freeMat[myCols[t] + j];
            }
        }
    }
}

bool GetSolution(const std::vector<std::uint8_t> &freeMat,
                 const std::vector<std::uint8_t> &mat,
                 const std::vector<std::size_t> &freeVariables,
                 const std::vector<mpz_class> &mpzFacBase,
                 const std::vector<mpz_class> &testInterval,
                 std::vector<mpz_class> &factors, const mpz_class &myNum,
                 std::size_t nCols, std::size_t matNCols, unsigned long int ind, 
                 std::size_t lenFree, std::size_t threadInd) {
    
    std::vector<std::size_t> ansVec;
    std::vector<std::uint8_t> posVec(nCols, 0u);
    const std::vector<std::uint8_t> posAns = MyIntToBit(ind, lenFree);
    
    bool bSuccess = false;
    
    for (std::size_t i = 0; i < freeVariables.size(); ++i)
        for (std::size_t k = 0, j = i * nCols; k < nCols; ++k, ++j)
            if (freeMat[j])
                posVec[k] ^= posAns[i];
    
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
            mpz_class mpzTemp1, mpzTemp2, mpzMin, xMpz, yMpz;
            yExponents.erase(yExponents.begin());
            xMpz = 1;
            yMpz = 1;
            
            for (const auto aV: ansVec) {
                xMpz *= testInterval[aV];
                xMpz %= myNum;
            }
            
            for (std::size_t j = 0; j < yExponents.size(); ++j) {
                mpz_pow_ui(mpzTemp1.get_mpz_t(),
                           mpzFacBase[j].get_mpz_t(), yExponents[j]);
                yMpz *= mpzTemp1;
                yMpz %= myNum;
            }
            
            mpzTemp1 = gcd(xMpz - yMpz, myNum);
            mpzTemp2 = gcd(xMpz + yMpz, myNum);
            
            if (cmp(mpzTemp1, mpzTemp2) < 0) {
                mpzMin = mpzTemp1;
            } else {
                mpzMin = mpzTemp2;
            }
            
            if (cmp(mpzMin, 1) > 0) {
                if (cmp(mpzTemp1, mpzTemp2) < 0) {
                    factors[threadInd * 2] = mpzTemp1;
                    factors[threadInd * 2 + 1] = mpzTemp2;
                } else {
                    factors[threadInd * 2 + 1] = mpzTemp1;
                    factors[threadInd * 2] = mpzTemp2;
                }
                
                bSuccess = true;
            }
        }
    }
    
    return bSuccess;
}

void SolutionSearch(const std::vector<std::uint8_t> &mat, std::size_t matNRows,
                    std::size_t matNCols, const mpz_class &myNum,
                    const std::vector<mpz_class> &mpzFacBase,
                    const std::vector<mpz_class> &testInterval,
                    std::vector<mpz_class> &factors,
                    std::size_t nThreads, bool bShowStats) {
    
    const auto t0 = std::chrono::steady_clock::now();
    
    if (bShowStats) {
        RcppThread::Rcout << "|  Mat Algebra Time  |\n|--------------------|" << std::endl;
    }
    
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
    
    if (bShowStats) {
        OneColumnStats(std::chrono::steady_clock::now() - t0);
    }
    
    ReduceMatrix(nullMat, myCols, 
                 static_cast<int>(nCols),
                 static_cast<int>(nRows));
    
    if (bShowStats) {
        OneColumnStats(std::chrono::steady_clock::now() - t0);
    }
    
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
        mpz_class mpzTemp1, cppNum(myNum);;

        mpz_ui_pow_ui(mpzTemp1.get_mpz_t(), 2, lenFree);
        --mpzTemp1;
        
        const unsigned long int myLim = (cmp(mpzTemp1, std::numeric_limits<unsigned long int>::max()) > 0)
                                    ? std::numeric_limits<unsigned long int>::max() : mpzTemp1.get_ui();
        
        const std::size_t sampSize = nThreads * (((myLim > oneThousand)
                                                      ? oneThousand : myLim) / nThreads);
        
        bool bSuccess = false;
        std::mt19937 mersenne_engine(42);
        std::uniform_int_distribution<unsigned long int> dist(1, myLim);
        
        auto gen = [&dist, &mersenne_engine](){
            return dist(mersenne_engine);
        };
        
        std::vector<unsigned long int> sample(sampSize);
        std::generate(sample.begin(), sample.end(), gen);
        
        if (bShowStats) {
            OneColumnStats(std::chrono::steady_clock::now() - t0);
        }

        if (nThreads > 1) {
            std::vector<mpz_class> vecFactors(nThreads * 2);
            std::vector<std::future<bool>> myFutures(nThreads);
            std::vector<bool> vecSuccess(nThreads);

            for (std::size_t i = 0; i < sampSize && !bSuccess;) {
                RcppThread::ThreadPool pool(nThreads);

                for (std::size_t thrd = 0; thrd < nThreads; ++thrd, ++i) {
                    myFutures[thrd] = pool.pushReturn(std::cref(GetSolution), std::cref(freeMat),
                                                      std::cref(mat), std::cref(freeVariables),
                                                      std::cref(mpzFacBase), std::cref(testInterval),
                                                      std::ref(vecFactors), std::cref(cppNum), nCols,
                                                      matNCols, sample[i], lenFree, thrd);
                }
                
                pool.join();
                
                for (std::size_t j = 0; j < nThreads; ++j)
                    vecSuccess[j] = myFutures[j].get();
                
                bSuccess = std::any_of(vecSuccess.begin(), vecSuccess.end(),
                                       [](bool v) {return v;});
            }

            for (std::size_t j = 0; j < nThreads; ++j) {
                if (vecSuccess[j]) {
                    factors[0] = vecFactors[j * 2];
                    factors[1] = vecFactors[j * 2 + 1];
                    break;
                }
            }
        } else {
            for (std::size_t i = 0; i < sampSize && !bSuccess; ++i) {
                bSuccess = GetSolution(freeMat, mat, freeVariables, mpzFacBase,
                                       testInterval, factors, cppNum, nCols,
                                       matNCols, sample[i], lenFree, 0);
            }
        }
    }
    
    if (bShowStats) {
        OneColumnStats(std::chrono::steady_clock::now() - t0);
        RcppThread::Rcout << "\n" << std::endl;
    }
}
