#include "Polynomial.h"
#include "SolutionSearch.h"
#include "StatsUtils.h"
#include <fstream>
#include <set>

constexpr std::size_t InitialNumPolys = 50u;
constexpr auto minOneMin = std::chrono::seconds(60);
constexpr double dblMsOneMin = std::chrono::duration_cast<std::chrono::milliseconds>(minOneMin).count();

void Polynomial::MergeMaster(vec2dint &powsOfSmoothsBig, vec2dint &powsOfPartialsBig,
                             std::vector<std::size_t> &coFactorIndexVecBig,
                             hash64vec &partFactorsMapBig, hash64mpz &partIntvlMapBig,
                             hash64size_t &keepingTrackBig, std::vector<mpz_class> &smoothIntervalBig,
                             std::vector<mpz_class> &largeCoFactorsBig, 
                             std::vector<mpz_class> &partialIntervalBig) {
    
    RcppThread::Rcout << powsOfSmooths.size() << "\n";
    
    powsOfSmoothsBig.insert(powsOfSmoothsBig.end(),
                            std::make_move_iterator(powsOfSmooths.begin()),
                            std::make_move_iterator(powsOfSmooths.end())
    );
    
    powsOfPartialsBig.insert(powsOfPartialsBig.end(),
                             std::make_move_iterator(powsOfPartials.begin()),
                             std::make_move_iterator(powsOfPartials.end())
    );
    
    smoothIntervalBig.insert(smoothIntervalBig.end(),
                             std::make_move_iterator(smoothInterval.begin()),
                             std::make_move_iterator(smoothInterval.end())
    );
    
    largeCoFactorsBig.insert(largeCoFactorsBig.end(),
                             std::make_move_iterator(largeCoFactors.begin()),
                             std::make_move_iterator(largeCoFactors.end())
    );
    
    partialIntervalBig.insert(partialIntervalBig.end(),
                              std::make_move_iterator(partialInterval.begin()),
                              std::make_move_iterator(partialInterval.end())
    );
    
    coFactorIndexVecBig.insert(coFactorIndexVecBig.end(),
                               coFactorIndexVec.begin(), coFactorIndexVec.end());
    
    // hash64vec &partFactorsMapBig, hash64mpz &partIntvlMapBig,
    // hash64int &keepingTrackBig
    // const uint64_t myKey = static_cast<uint64_t>(IntVal.get_d());
    // const auto pFacIt = partFactorsMap.find(myKey);
    // 
    // if (pFacIt != partFactorsMap.end()) {
    //     const auto trackIt = keepingTrack.find(myKey);
    // 
    //     if (trackIt != keepingTrack.end()) {
    //         coFactorIndexVec.push_back(trackIt->second);
    //     } else {
    //         keepingTrack[myKey] = coFactorInd;
    //         largeCoFactors.push_back(IntVal);
    //         coFactorIndexVec.push_back(coFactorInd);
    //         ++coFactorInd;
    //     }
    // 
    //     primeIndexVec.insert(primeIndexVec.begin(),
    //                          pFacIt->second.cbegin(), pFacIt->second.cend());
    // 
    //     powsOfPartials.push_back(primeIndexVec);
    //     const auto intervalIt = partIntvlMap.find(myKey);
    // 
    //     Temp += VarB;
    //     partialInterval.push_back(Temp * intervalIt->second);
    // 
    //     partFactorsMap.erase(pFacIt);
    //     partIntvlMap.erase(intervalIt);
    //     ++nPartial;
    // } else {
    //     partFactorsMap[myKey] = primeIndexVec;
    //     partIntvlMap[myKey] = Temp + VarB;
    // }
    // 
    // for (const auto &pFac: partFactorsMap) {
    //     const auto pFacIt = partFactorsMapBig.find(pFac.first);
    //     
    //     if (pFacIt != partFactorsMap.end()) {
    //         const auto trackIt = keepingTrack.find(pFac.first);
    //         
    //         if (trackIt != keepingTrack.end()) {
    //             coFactorIndexVec.push_back(trackIt->second);
    //         } else {
    //             keepingTrack[pFac.first] = coFactorInd;
    //             largeCoFactors.push_back(IntVal);
    //             coFactorIndexVec.push_back(coFactorInd);
    //             ++coFactorInd;
    //         }
    //         
    //         primeIndexVec.insert(primeIndexVec.begin(),
    //                              pFacIt->second.cbegin(), pFacIt->second.cend());
    //         
    //         powsOfPartials.push_back(primeIndexVec);
    //         const auto intervalIt = partIntvlMap.find(pFac.first);
    //         
    //         Temp += VarB;
    //         partialInterval.push_back(Temp * intervalIt->second);
    //         
    //         partFactorsMap.erase(pFacIt);
    //         partIntvlMap.erase(intervalIt);
    //         ++nPartial;
    //     } else {
    //         partFactorsMap[pFac.first] = primeIndexVec;
    //         partIntvlMap[pFac.first] = Temp + VarB;
    //     }
    // }
    
}

Polynomial::Polynomial(std::size_t _mpzContainerSize,
                       std::size_t _facSize, bool _bShowStats, mpz_class myNum) : 
            mpzFacSize(_facSize), SaPThresh(_facSize),
            facSize(_facSize), bShowStats(_bShowStats) {
    
    powsOfSmooths.reserve(_facSize);
    powsOfPartials.reserve(_facSize);
    myStart.assign(_facSize * 2, 0u);

    smoothInterval.reserve(_mpzContainerSize);
    largeCoFactors.reserve(_mpzContainerSize);
    partialInterval.reserve(_mpzContainerSize);

    nPolys = 0;
    nPartial = 0;
    nSmooth = 0;
    coFactorInd = 0;
    
    if (bShowStats) {
        RcppThread::Rcout << "\nSummary Statistics for Factoring:\n" << "    "
                          << myNum.get_str() << "\n\n"
                          << "|        Time        | Complete | Polynomials |   Smooths"
                          << "  |  Partials  |\n";
        RcppThread::Rcout << "|--------------------|----------|-------------|---------"
                          << "---|------------|" << std::endl;
    }
}

Polynomial::Polynomial(std::size_t _facSize) : 
    SaPThresh(_facSize), facSize(_facSize), bShowStats(false) {

    myStart.assign(_facSize * 2, 0u);
    nPolys = 0;
    nPartial = 0;
    nSmooth = 0;
    coFactorInd = 0;
}

void GetNPrimes(std::vector<mpz_class> &mpzFacBase, mpz_class NextPrime,
                mpz_class myNum, std::size_t numPrimesNeeded) {
    
    for (std::size_t i = 0; i < numPrimesNeeded; ++i) {
        for (bool LegendreTest = true; LegendreTest; ) {
            mpz_nextprime(NextPrime.get_mpz_t(), NextPrime.get_mpz_t());
            
            if (mpz_legendre(myNum.get_mpz_t(), NextPrime.get_mpz_t()) == 1)
                LegendreTest = false;
        }
        
        mpzFacBase.push_back(NextPrime);
    }
}

void Polynomial::SievePolys(const std::vector<std::size_t> &SieveDist,
                            const std::vector<int> &facBase,
                            const std::vector<int> &LnFB, 
                            const std::vector<mpz_class> &mpzFacBase,
                            mpz_class LowBound, mpz_class myNum, int theCut,
                            int DoubleLenB, int vecMaxSize,
                            std::size_t strt, std::size_t polyLimit) {
    
    RcppThread::Rcout << GetSmoothSize() << "\n";
    
    for (std::size_t poly = 0; poly < polyLimit; ++poly) {
        ++mpzFacSize;

        SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                   coFactorIndexVec, myStart, partFactorsMap, partIntvlMap,
                   keepingTrack, smoothInterval, largeCoFactors,
                   partialInterval, mpzFacBase[mpzFacSize - 1], LowBound, myNum,
                   nPartial, nSmooth, coFactorInd, theCut, DoubleLenB,
                   mpzFacSize, vecMaxSize, strt);
    }
    
    RcppThread::Rcout << GetSmoothSize() << "\n";
}

void Polynomial::InitialParSieve(const std::vector<std::size_t> &SieveDist,
                                 const std::vector<int> &facBase, 
                                 const std::vector<int> &LnFB,
                                 std::vector<mpz_class> &mpzFacBase, mpz_class NextPrime,
                                 mpz_class LowBound, mpz_class myNum, int theCut,
                                 int DoubleLenB, int vecMaxSize, std::size_t strt,
                                 typeTimePoint checkPoint0) {
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    
    auto showStatsTime = (checkPoint1 - checkPoint0);
    GetNPrimes(mpzFacBase, NextPrime, myNum, InitialNumPolys);
    
    SievePolys(SieveDist, facBase, LnFB, mpzFacBase,
               LowBound, myNum, theCut, DoubleLenB,
               vecMaxSize, strt, InitialNumPolys);
    
    nPolys = InitialNumPolys;
    const auto checkPoint3 = std::chrono::steady_clock::now();
    
    if ((checkPoint3 - checkPoint1) > checkInterTime) {
        // Check for user interrupt and udpate timepoint
        RcppThread::checkUserInterrupt();
        checkPoint1 = std::chrono::steady_clock::now();
    }
    
    if (bShowStats && (checkPoint3 - checkPoint2) > showStatsTime) {
        MakeStats(SaPThresh, nPolys, nSmooth,
                  nPartial, checkPoint3 - checkPoint0);
        
        checkPoint2 = std::chrono::steady_clock::now();
        UpdateStatTime(nSmooth + nPartial, facSize,
                       checkPoint3 - checkPoint0, showStatsTime);
    }
}

inline std::size_t SetNumThreads(std::size_t currLim, std::size_t SaPThresh, std::size_t nThreads) {
    
    const std::size_t MaxThreads = ((2 * currLim) > SaPThresh) ? 0u :
                                         (currLim == 0u) ? nThreads : 
                                         (SaPThresh - currLim) / currLim;
    
    return std::min(nThreads, MaxThreads);
}

template <typename typeTimeDiff>
inline std::size_t SetNumPolys(std::size_t currLim, std::size_t SaPThresh, std::size_t nThreads,
                               int nPolys, int polysPerSec, int SaPOnLast,
                               typeTimeDiff totalTime, typeTimeDiff polyTime) {
    
    const int percentComplete = ((100 * currLim) / SaPThresh) + 1;
    const auto nPercentTime = nThreads * totalTime / percentComplete;
    
    const double dblPolyTime = std::chrono::duration_cast<
                    std::chrono::milliseconds>(polyTime).count();
    
    const double dblOnePerc = std::chrono::duration_cast<
                    std::chrono::milliseconds>(nPercentTime).count();
    
    const double numerator = (nPercentTime > minOneMin) ? dblMsOneMin :
                             (nPercentTime < std::chrono::seconds(1)) ?
                             (5 * dblOnePerc) : dblOnePerc;
    
    const double myRatio = numerator / dblPolyTime;
    
    const int calcNumPolys = ((SaPThresh - currLim) < (2 * SaPOnLast)) ?
        static_cast<double>(polysPerSec) * static_cast<double>(SaPThresh - currLim) /
            static_cast<double>(SaPOnLast) : myRatio * static_cast<double>(polysPerSec);
    
    return calcNumPolys;
}

void Polynomial::FactorParallel(const std::vector<std::size_t> &SieveDist,
                                const std::vector<int> &facBase, const std::vector<int> &LnFB,
                                std::vector<mpz_class> &mpzFacBase, mpz_class NextPrime,
                                mpz_class LowBound, mpz_class myNum, int theCut,
                                int DoubleLenB, int vecMaxSize, std::size_t strt,
                                typeTimePoint checkPoint0, std::size_t nThreads) {
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    auto showStatsTime = (checkPoint1 - checkPoint0);
    
    RcppThread::Rcout << "Start\n";
    
    this->InitialParSieve(SieveDist, facBase, LnFB, mpzFacBase,
                          NextPrime, LowBound, myNum, theCut,
                          DoubleLenB, vecMaxSize, strt, checkPoint0);
    
    RcppThread::Rcout << "Finished initial sieve: " << nSmooth << " " << nPartial << " " << nPolys << "\n";
    
    nThreads = SetNumThreads(nSmooth + nPartial, SaPThresh, nThreads);
    std::size_t polysPerThread = SetNumPolys(nSmooth + nPartial, SaPThresh, 1,
                                             nPolys, InitialNumPolys, 0, showStatsTime,
                                             showStatsTime);
    polysPerThread = 96215 / nThreads;
    RcppThread::Rcout << smoothInterval.size() << " " << partialInterval.size() << "\n";
    RcppThread::Rcout << polysPerThread << " " << nThreads << "\n";
    
    // This variable keeps track of how many smooths + partials
    // were found in one iteration.
    std::size_t SaPOnLast = nSmooth + nPartial;
    
    // while ((nSmooth + nPartial) <= SaPThresh) {
        std::vector<std::unique_ptr<Polynomial>> vecPoly;
        std::vector<std::thread> myThreads;
        
        const std::size_t startStep = mpzFacBase.size();
        GetNPrimes(mpzFacBase, mpzFacBase.back(), myNum, polysPerThread * nThreads);
        
        for (std::size_t i = 0, step = startStep; i < nThreads; ++i, step += polysPerThread) {
            vecPoly.push_back(FromCpp14::make_unique<Polynomial>(facSize));
            vecPoly[i]->SetMpzFacSize(step);

            myThreads.emplace_back(&Polynomial::SievePolys, vecPoly[i].get(),
                                   std::cref(SieveDist), std::cref(facBase), std::cref(LnFB),
                                   std::cref(mpzFacBase), LowBound, myNum, theCut, DoubleLenB,
                                   vecMaxSize, strt, polysPerThread);
        }

        for (auto &thr: myThreads)
            thr.join();
        
        RcppThread::Rcout << vecPoly.size() << " " << powsOfSmooths.size() << "\n";
        
        for (std::size_t i = 0; i < nThreads; ++i) {
            vecPoly[i]->MergeMaster(powsOfSmooths, powsOfPartials, coFactorIndexVec,
                                    partFactorsMap, partIntvlMap, keepingTrack,
                                    smoothInterval, largeCoFactors, partialInterval);
        }
        
        RcppThread::Rcout << smoothInterval.size() << " " << partialInterval.size() << "\n";
        
        // 
        // const auto checkPoint3 = std::chrono::steady_clock::now();
        // 
        // if ((checkPoint3 - checkPoint1) > checkInterTime) {
        //     // Check for user interrupt and udpate timepoint
        //     RcppThread::checkUserInterrupt();
        //     checkPoint1 = std::chrono::steady_clock::now();
        // }
        // 
        // if (bShowStats && (checkPoint3 - checkPoint2) > showStatsTime) {
        //     MakeStats(SaPThresh, nPolys, nSmooth,
        //               nPartial, checkPoint3 - checkPoint0);
        //     
        //     checkPoint2 = std::chrono::steady_clock::now();
        //     UpdateStatTime(nSmooth + nPartial, facSize,
        //                    checkPoint3 - checkPoint0, showStatsTime);
        // }
    // }
}

void Polynomial::FactorSerial(const std::vector<std::size_t> &SieveDist, const std::vector<int> &facBase,
                              const std::vector<int> &LnFB, std::vector<mpz_class> &mpzFacBase,
                              mpz_class NextPrime, mpz_class LowBound, mpz_class myNum,
                              int theCut, int DoubleLenB, int vecMaxSize, std::size_t strt,
                              typeTimePoint checkPoint0) {
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    auto showStatsTime = (checkPoint1 - checkPoint0);

    while ((nSmooth + nPartial) <= SaPThresh) {
        for (bool LegendreTest = true; LegendreTest; ) {
            mpz_nextprime(NextPrime.get_mpz_t(), NextPrime.get_mpz_t());

            if (mpz_legendre(myNum.get_mpz_t(), NextPrime.get_mpz_t()) == 1)
                LegendreTest = false;
        }

        mpzFacBase.push_back(NextPrime);
        ++mpzFacSize;

        SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                   coFactorIndexVec, myStart, partFactorsMap, partIntvlMap,
                   keepingTrack, smoothInterval, largeCoFactors,
                   partialInterval, NextPrime, LowBound, myNum,
                   nPartial, nSmooth, coFactorInd, theCut, DoubleLenB,
                   mpzFacSize, vecMaxSize, strt);

        ++nPolys;
        const auto checkPoint3 = std::chrono::steady_clock::now();

        if ((checkPoint3 - checkPoint1) > checkInterTime) {
            // Check for user interrupt and udpate timepoint
            RcppThread::checkUserInterrupt();
            checkPoint1 = std::chrono::steady_clock::now();
        }

        if (bShowStats && (checkPoint3 - checkPoint2) > showStatsTime) {
            MakeStats(SaPThresh, nPolys, nSmooth,
                      nPartial, checkPoint3 - checkPoint0);

            checkPoint2 = std::chrono::steady_clock::now();
            UpdateStatTime(nSmooth + nPartial, facSize,
                           checkPoint3 - checkPoint0, showStatsTime);
        }
    }
    
    SaPThresh += 10;
}

void Polynomial::GetSolution(const std::vector<mpz_class> &mpzFacBase,
                             const std::vector<int> &facBase, mpz_t *const factors,
                             mpz_t mpzNum, std::size_t nThreads,
                             typeTimePoint checkPoint0) {
    
    // Not every prime in mpzFacBase is utilized. For example, with our
    // attempt at factoring rsa99, there were over 4 million elements
    // in mpzFacBase, however there were only ~40000 smooths + partials.
    // This caused an inordinate memory allocation haulting the
    // factorization. The set, setIndex, will be used to drastically
    // reduce the size of the factor base and related structures making
    // the linear algebra portion much more efficient.
    std::set<int> setIndex;
    auto indIt = setIndex.begin();

    for (std::size_t i = 0; i < nSmooth; ++i) {
        if (powsOfSmooths[i].front() > facSize) {
            setIndex.emplace_hint(indIt, powsOfSmooths[i].front());
            indIt = setIndex.end();
        }
    }

    // Here we get the 4 largest powers as we know these relate to
    // primes not contained in facBase.
    for (std::size_t i = 0; i < nPartial; ++i) {
        std::partial_sort(powsOfPartials[i].begin(), powsOfPartials[i].begin() + 4,
                          powsOfPartials[i].end(), std::greater<int>());

        if (powsOfPartials[i].front() > facSize)
            setIndex.insert(powsOfPartials[i].front());

        if (powsOfPartials[i][2] > facSize)
            setIndex.insert(powsOfPartials[i][2]);
    }

    const std::size_t nRows = nSmooth + nPartial;
    const std::size_t nonTrivSize = setIndex.size() + facSize;

    const std::size_t nCols = nonTrivSize + coFactorInd + 1;
    std::vector<mpz_class> nonTrivialFacs(nCols);

    for (std::size_t i = 0; i < facSize; ++i)
        nonTrivialFacs[i] = facBase[i];

    std::unordered_map<int, int> mapIndex;
    mapIndex.reserve(setIndex.size());

    int setInd = facSize;

    for (auto s: setIndex) {
        nonTrivialFacs[setInd++] = mpzFacBase[s - 1];
        mapIndex[s] = setInd;
    }

    for (std::size_t i = 0, j = setInd; i < coFactorInd; ++i, ++j)
        nonTrivialFacs[j] = largeCoFactors[i];

    std::vector<mpz_class> newTestInt(nRows);
    std::vector<std::uint8_t> mat(nRows * nCols, 0u);

    for (std::size_t r = 0, row = 0; r < nSmooth; ++r, row += nCols) {
        for (std::size_t j = 0; j < 2; ++j)
            ++mat[row + mapIndex[powsOfSmooths[r][j]]];
        
        for (auto it = powsOfSmooths[r].begin() + 2; it != powsOfSmooths[r].end(); ++it)
            ++mat[row + *it];
        
        newTestInt[r] = smoothInterval[r];
    }

    for (std::size_t i = 0, r = nSmooth,
         row = nCols * nSmooth; i < nPartial; ++i, ++r, row += nCols) {

        // Remap the powers corresponding to powers not contained in facBase
        for (std::size_t j = 0; j < 4; ++j)
            ++mat[row + mapIndex[powsOfPartials[i][j]]];
        
        for (auto it = powsOfPartials[i].begin() + 4; it != powsOfPartials[i].end(); ++it)
            ++mat[row + *it];

        mat[row + nonTrivSize + coFactorIndexVec[i] + 1] = 2u;
        newTestInt[r] = partialInterval[i];
    }

    SolutionSearch(mat, nRows, nCols, mpzNum, nonTrivialFacs,
                   newTestInt, factors, nThreads);

    if (bShowStats && mpz_cmp_ui(factors[0], 0)) {
        const auto checkPoint2 = std::chrono::steady_clock::now();

        MakeStats(nSmooth + nPartial, nPolys, nSmooth,
                  nPartial, checkPoint2 - checkPoint0);

        RcppThread::Rcout << "\n" << std::endl;
    }
}
