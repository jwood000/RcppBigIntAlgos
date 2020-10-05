#include "SolutionSearch.h"
#include "StatsUtils.h"
#include <random>

std::vector<std::uint8_t> MyIntToBit(std::size_t x, std::size_t dig) {
    
    std::vector<std::uint8_t> binaryVec(dig);
    
    for (std::size_t i = 0; x > 0; ++i) {
        binaryVec[i] = x % 2;
        x >>= 1;
    }
    
    return binaryVec;
}

void ProcessFreeMat(const std::vector<std::bitset<wordSize>> &nullMat,
                    const std::vector<std::size_t> &myCols,
                    std::vector<std::uint8_t> &freeMat,
                    std::size_t newNrow, std::size_t nCols) {
    
    const std::size_t freeMatSize = freeMat.size();
    const std::size_t adjustedCols = (nCols + wordSize - 1) / wordSize;

    for (int i = newNrow - 1; i >= 0; --i) {
        std::vector<std::size_t> nonTriv;
        
        for (std::size_t j = i + 1, myRow = i * adjustedCols; j < nCols; ++j)
            if (nullMat[myRow + j / wordSize].test(j % wordSize))
                nonTriv.push_back(j);
        
        if (!nonTriv.empty()) {
            if (nonTriv.front() >= newNrow) {
                for (std::size_t t = 0, col1 = myCols[i]; t < nonTriv.size(); ++t) {
                    for (std::size_t j = 0,
                         col2 = myCols[nonTriv[t]]; j < freeMatSize; j += nCols) {
                        if (freeMat[col2 + j])
                            freeMat[col1 + j] = 1u;
                    }
                }
            } else {
                for (std::size_t t = 0, col1 = myCols[i]; t < nonTriv.size(); ++t) {
                    for (std::size_t j = 0,
                         col2 = myCols[nonTriv[t]]; j < freeMatSize; j += nCols) {
                        freeMat[col1 + j] ^= freeMat[col2 + j];
                    }
                }
            }
        }
    }
}

char GetSolution(const std::vector<std::uint8_t> &freeMat,
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
    
    char bSuccess = 0;
    
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
                
                bSuccess = 1;
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
        RcppThread::Rcout << "|  Mat Algebra Time  |    Mat Dimension   |\n"
                          << "|--------------------|--------------------|" << std::endl;
        TwoColumnStats(std::chrono::steady_clock::now() - t0, matNCols, matNRows);
    }
    
    const std::size_t matSize = mat.size();
    const std::size_t nCols = matNRows;
    const std::size_t adjustedCols = (nCols + wordSize - 1) / wordSize;
    std::size_t nRows = 0;
    
    std::vector<std::bitset<wordSize>> nullMat;
    const std::size_t maxNullSize = (matSize + wordSize - 1u) / wordSize;
    nullMat.reserve(maxNullSize);
    
    for (std::size_t j = 0; j < matNCols; ++j) {
        std::size_t i = 0;
        
        while ((i < matSize) && ((mat[i + j] % 2u) == 0))
            i += matNCols;
        
        if (i < matSize) {
            for (std::size_t k = 0; k < matSize;) {
                std::bitset<wordSize> num;
                
                for (std::size_t r = 0; r < wordSize && k < matSize; ++r, k += matNCols)
                    if (mat[k + j] % 2u)
                        num.set(r);
                
                nullMat.push_back(num);
            }
            
            ++nRows;
        }
    }
    
    std::vector<std::size_t> myCols(nCols, 0);
    std::iota(myCols.begin(), myCols.end(), 0);
    
    if (bShowStats) {
        TwoColumnStats(std::chrono::steady_clock::now() - t0, nRows, nCols);
    }
    
    ReduceMatrix(nullMat, myCols, nCols, nRows);
    
    if (bShowStats) {
        TwoColumnStats(std::chrono::steady_clock::now() - t0, nRows, nCols);
    }
    
    const std::size_t newNrow =  nullMat.size() / adjustedCols;
    std::vector<std::size_t> freeVariables;

    if (nCols > newNrow && newNrow > 0) {
        for (std::size_t i = newNrow; i < nCols; ++i)
            freeVariables.push_back(myCols[i]);

        std::sort(freeVariables.begin(), freeVariables.end());
        const std::size_t myMin = freeVariables.front();

        const std::size_t lenFree = freeVariables.size();
        std::vector<std::uint8_t> freeMat(lenFree * nCols);

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
        mpz_class mpzTemp1, cppNum(myNum);

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
            TwoColumnStats(std::chrono::steady_clock::now() - t0, nRows, nCols);
        }

        if (nThreads > 1) {
            std::vector<mpz_class> vecFactors(nThreads * 2);
            std::vector<std::future<char>> myFutures(nThreads);
            std::vector<char> vecSuccess(nThreads);

            for (std::size_t i = 0; i < sampSize && !bSuccess;) {
                RcppThread::ThreadPool pool(nThreads);

                for (std::size_t thrd = 0; thrd < nThreads; ++thrd, ++i) {
                    myFutures[thrd] = pool.pushReturn(GetSolution, std::cref(freeMat),
                                                      std::cref(mat), std::cref(freeVariables),
                                                      std::cref(mpzFacBase), std::cref(testInterval),
                                                      std::ref(vecFactors), std::cref(cppNum), nCols,
                                                      matNCols, sample[i], lenFree, thrd);
                }

                pool.join();

                for (std::size_t j = 0; j < nThreads; ++j)
                    vecSuccess[j] = myFutures[j].get();

                bSuccess = std::any_of(vecSuccess.begin(), vecSuccess.end(),
                                       [](char v) {return v;});
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
        TwoColumnStats(std::chrono::steady_clock::now() - t0, nRows, nCols);
        RcppThread::Rcout << "\n" << std::endl;
    }
}
