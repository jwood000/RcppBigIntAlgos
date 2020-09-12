#include "SolutionSearch.h"
#include "Cpp14MakeUnique.h"
#include <RcppThread.h>
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
            if (temp.front() >= newNrow) {
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
                 mpz_t *const factors, mpz_class myNum, std::size_t nCols, 
                 std::size_t matNCols, std::size_t ind, 
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
                mpz_pow_ui(mpzTemp1.get_mpz_t(), mpzFacBase[j].get_mpz_t(), yExponents[j]);
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
                    mpz_set(factors[threadInd * 2], mpzTemp1.get_mpz_t());
                    mpz_set(factors[threadInd * 2 + 1], mpzTemp2.get_mpz_t());
                } else {
                    mpz_set(factors[threadInd * 2 + 1], mpzTemp1.get_mpz_t());
                    mpz_set(factors[threadInd * 2], mpzTemp2.get_mpz_t());
                }
                
                bSuccess = true;
            }
        }
    }
    
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
        mpz_class mpzTemp1, cppNum(myNum);;

        mpz_ui_pow_ui(mpzTemp1.get_mpz_t(), 2, lenFree);
        --mpzTemp1;
        
        const std::size_t myLim = (cmp(mpzTemp1, oneThousand) > 0)
                                    ? oneThousand : mpzTemp1.get_ui();

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
                                                   vecFactors.get(), cppNum, nCols, matNCols, i, lenFree, j);
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
            
            // First create an instance of an engine.
            std::random_device rnd_device;
            // Specify the engine and distribution.
            std::mt19937 mersenne_engine {rnd_device()};  // Generates random integers
            std::uniform_int_distribution<std::size_t> dist{1, myLim};
            
            auto gen = [&dist, &mersenne_engine](){
                return dist(mersenne_engine);
            };
            
            std::vector<std::size_t> sample(myLim - 1);
            std::generate(sample.begin(), sample.end(), gen);
            
            for (std::size_t i = 0; i < sample.size() && !bSuccess; ++i) {
                bSuccess = GetSolution(freeMat, mat, freeVariables, mpzFacBase,
                                       testInterval, factors, cppNum, nCols,
                                       matNCols, sample[i], lenFree, 0);
            }
        }
    }
}
