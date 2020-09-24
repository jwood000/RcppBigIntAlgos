#include "Polynomial.h"
#include "SolutionSearch.h"
#include "StatsUtils.h"
#include <fstream>
#include <set>

constexpr std::size_t InitialNumPolys = 50u;
constexpr std::size_t bigFacsTWO = 2u;
constexpr std::size_t bigFacsFOUR = 4u;

constexpr double dbl1000ms = 1000;
constexpr auto sOneMin = std::chrono::seconds(60);
constexpr double dblMsOneMin = std::chrono::duration_cast<std::chrono::milliseconds>(sOneMin).count();

void Polynomial::MergeMaster(vec2dint &powsOfSmoothsBig, vec2dint &powsOfPartialsBig,
                             hash64vec &partFactorsMapBig, hash64mpz &partIntvlMapBig,
                             std::vector<mpz_class> &smoothIntervalBig,
                             std::vector<uint64_t> &largeCoFactorsBig, 
                             std::vector<mpz_class> &partialIntervalBig) {
    
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
    
    partialIntervalBig.insert(partialIntervalBig.end(),
                              std::make_move_iterator(partialInterval.begin()),
                              std::make_move_iterator(partialInterval.end())
    );
    
    largeCoFactorsBig.insert(largeCoFactorsBig.end(),
                             std::make_move_iterator(largeCoFactors.begin()),
                             std::make_move_iterator(largeCoFactors.end())
    );
    
    std::vector<uint64_t> deleteLater;
    
    // First identify intersection
    for (const auto &pFac: partFactorsMap) {
        const uint64_t myKey = pFac.first;
        const auto pFacBigIt = partFactorsMapBig.find(myKey);

        if (pFacBigIt != partFactorsMapBig.end()) {
            largeCoFactorsBig.push_back(myKey);
            pFacBigIt->second.insert(pFacBigIt->second.end(),
                                     pFac.second.cbegin(), pFac.second.cend());
            
            powsOfPartialsBig.push_back(pFacBigIt->second);
            partialIntervalBig.push_back(partIntvlMap[myKey] * partIntvlMapBig[myKey]);
            deleteLater.push_back(myKey);
        }
    }
    
    for (auto myKey: deleteLater) {
        partFactorsMap.erase(myKey);
        partFactorsMapBig.erase(myKey);
        partIntvlMap.erase(myKey);
        partIntvlMapBig.erase(myKey);
    }
    
    partFactorsMapBig.insert(partFactorsMap.begin(), partFactorsMap.end());
    partIntvlMapBig.insert(partIntvlMap.begin(), partIntvlMap.end());
}

Polynomial::Polynomial(std::size_t _mpzContainerSize,
                       std::size_t _facSize, bool _bShowStats, const mpz_class &myNum) : 
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
}

void GetNPrimes(std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                const mpz_class &myNum, std::size_t numPrimesNeeded) {
    
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
                            const std::vector<int> &facBase, const std::vector<int> &LnFB, 
                            const std::vector<mpz_class> &mpzFacBase,
                            const mpz_class &LowBound, const mpz_class &myNum,
                            int theCut, int DoubleLenB, int vecMaxSize,
                            std::size_t strt, std::size_t polyLimit) {
    
    for (std::size_t poly = 0; poly < polyLimit; ++poly) {
        ++mpzFacSize;

        SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                   myStart, partFactorsMap, partIntvlMap, smoothInterval,
                   largeCoFactors, partialInterval, mpzFacBase[mpzFacSize - 1],
                   LowBound, myNum, nPartial, nSmooth, theCut, DoubleLenB,
                   mpzFacSize, vecMaxSize, strt);
    }
}

void Polynomial::InitialParSieve(const std::vector<std::size_t> &SieveDist,
                                 const std::vector<int> &facBase, const std::vector<int> &LnFB,
                                 std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                                 const mpz_class &LowBound, const mpz_class &myNum, int theCut,
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

inline void SetNumThreads(std::size_t currLim, std::size_t prevLim,
                          std::size_t SaPThresh, std::size_t &nThreads,
                          double polysPerThread) {
    
    if (currLim >= SaPThresh) {
        nThreads = 0;
    } else {
        
        double test = polysPerThread * static_cast<double>(SaPThresh - currLim)
            / static_cast<double>(currLim - prevLim);
        
        while (test < InitialNumPolys && nThreads) {
            nThreads /= 2;
            test *= 2;
        }
    }
}

template <typename typeTimeDiff>
inline std::size_t SetNumPolys(std::size_t currLim, std::size_t prevLim,
                               std::size_t SaPThresh, std::size_t &nThreads,
                               std::size_t oldThreads, double nPolys,
                               typeTimeDiff polyTime, typeTimeDiff totalTime) {
    
    const double percentComplete = (100 * currLim) / SaPThresh;
    const double percentSection = std::max((100.0 * (currLim - prevLim)) / SaPThresh, 1.0);
    const double percentRemaining = 100 - percentComplete;
    
    const auto onePercentTime = totalTime / static_cast<std::size_t>(percentComplete + 1);
    
    const double dblPolyTime = std::chrono::duration_cast<
        std::chrono::milliseconds>(polyTime).count();
    
    double myRatio = 1.0;
    const double fudge = 2 - (percentComplete / 100);
    const double thrdRatio = static_cast<double>(oldThreads) / 
        static_cast<double>(nThreads);
    
    if (percentRemaining < percentSection) {
        myRatio = (percentRemaining / (fudge * percentSection));
    } else {
        if (onePercentTime > sOneMin) {
            myRatio = (dblMsOneMin / dblPolyTime);
        } else {
            
            // Try to get 25%. If it ends up greater than 15 seconds
            // or we overshoot percentRemaining, we dial it back
            if (percentRemaining > 33.0) {
                myRatio = (25.0 / percentSection);
                
                if ((thrdRatio * myRatio * polyTime) > std::chrono::seconds(15)) {
                    myRatio = dblMsOneMin / (dblPolyTime * 4);
                }
            } else if (percentRemaining < percentSection) {
                myRatio = (dbl1000ms / dblPolyTime);
                
                if (percentRemaining < (thrdRatio * myRatio * percentSection)) {
                    myRatio = (percentRemaining / (fudge * percentSection));
                }
            }
        }
    }
    
    myRatio *= thrdRatio;
    std::size_t calcNumPolys = static_cast<std::size_t>(myRatio * nPolys);
    
    if (calcNumPolys < InitialNumPolys) {
        nThreads = 0;
    }
    
    return calcNumPolys;
}

void Polynomial::FactorParallel(const std::vector<std::size_t> &SieveDist,
                                const std::vector<int> &facBase, const std::vector<int> &LnFB,
                                std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                                const mpz_class &LowBound, const mpz_class &myNum, int theCut,
                                int DoubleLenB, int vecMaxSize, std::size_t strt,
                                typeTimePoint checkPoint0, std::size_t nThreads) {
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    auto showStatsTime = checkPoint1 - checkPoint0;
    
    this->InitialParSieve(SieveDist, facBase, LnFB, mpzFacBase,
                          NextPrime, LowBound, myNum, theCut,
                          DoubleLenB, vecMaxSize, strt, checkPoint0);
    
    auto checkPoint3 = std::chrono::steady_clock::now();
    auto polyTime = checkPoint3 - checkPoint1;
    
    SetNumThreads(nSmooth + nPartial, 0u, SaPThresh, nThreads, 1u);
    std::size_t prevThread = 1;
    
    std::size_t PrevSaP = 0;
    std::size_t polysPerThread = InitialNumPolys;
    
    while ((nSmooth + nPartial) <= SaPThresh && nThreads > 1) {
        std::vector<std::unique_ptr<Polynomial>> vecPoly;
        std::vector<std::thread> myThreads;
        
        const std::size_t startStep = mpzFacBase.size();
        NextPrime = mpzFacBase.back();
        
        polysPerThread = SetNumPolys(nSmooth + nPartial, PrevSaP, SaPThresh,
                                     nThreads, prevThread, polysPerThread,
                                     polyTime, checkPoint3 - checkPoint0);
        
        if (nThreads) {
            PrevSaP = nSmooth + nPartial;
            prevThread = nThreads;
            
            GetNPrimes(mpzFacBase, NextPrime, myNum, polysPerThread * nThreads);
            mpzFacSize = mpzFacBase.size();
            
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
            
            for (std::size_t i = 0; i < nThreads; ++i) {
                vecPoly[i]->MergeMaster(powsOfSmooths, powsOfPartials, partFactorsMap,
                                        partIntvlMap, smoothInterval,
                                        largeCoFactors, partialInterval);
            }
            
            nSmooth = smoothInterval.size();
            nPartial = partialInterval.size();
            nPolys += (nThreads * polysPerThread);
            
            SetNumThreads(nSmooth + nPartial, PrevSaP,
                          SaPThresh, nThreads, polysPerThread);
            
            const auto checkPoint4 = std::chrono::steady_clock::now();
            polyTime = checkPoint4 - checkPoint3;
            checkPoint3 = checkPoint4;
    
            if ((checkPoint4 - checkPoint1) > checkInterTime) {
                // Check for user interrupt and udpate timepoint
                RcppThread::checkUserInterrupt();
                checkPoint1 = std::chrono::steady_clock::now();
            }
    
            if (bShowStats && (checkPoint4 - checkPoint2) > showStatsTime) {
                MakeStats(SaPThresh, nPolys, nSmooth,
                          nPartial, checkPoint4 - checkPoint0);
    
                checkPoint2 = std::chrono::steady_clock::now();
                UpdateStatTime(nSmooth + nPartial, facSize,
                               checkPoint4 - checkPoint0, showStatsTime);
            }
        }
    }
    
    SaPThresh += 10;
}

void Polynomial::FactorSerial(const std::vector<std::size_t> &SieveDist,
                              const std::vector<int> &facBase, const std::vector<int> &LnFB,
                              std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                              const mpz_class &LowBound, const mpz_class &myNum, int theCut,
                              int DoubleLenB, int vecMaxSize, std::size_t strt,
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
                   myStart, partFactorsMap, partIntvlMap, smoothInterval,
                   largeCoFactors, partialInterval, NextPrime, LowBound,
                   myNum, nPartial, nSmooth, theCut, DoubleLenB,
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
                             const std::vector<int> &facBase, std::vector<mpz_class> &factors,
                             const mpz_class &mpzNum, std::size_t nThreads,
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
        if (powsOfSmooths[i].front() > static_cast<int>(facSize)) {
            setIndex.emplace_hint(indIt, powsOfSmooths[i].front());
            indIt = setIndex.end();
        }
    }

    // Here we get the 4 largest powers as we know these relate to
    // primes not contained in facBase.
    for (std::size_t i = 0; i < nPartial; ++i) {
        std::partial_sort(powsOfPartials[i].begin(), powsOfPartials[i].begin() + 4,
                          powsOfPartials[i].end(), std::greater<int>());

        if (powsOfPartials[i].front() > static_cast<int>(facSize))
            setIndex.insert(powsOfPartials[i].front());

        if (powsOfPartials[i][2] > static_cast<int>(facSize))
            setIndex.insert(powsOfPartials[i][2]);
    }

    const std::size_t nRows = nSmooth + nPartial;
    const std::size_t nonTrivSize = setIndex.size() + facSize;
    
    std::size_t coFactorInd = 0;
    std::vector<std::size_t> coFactorIndexVec;
    coFactorIndexVec.reserve(largeCoFactors.size());
    
    std::unordered_map<uint64_t, std::size_t> keepingTrack;
    keepingTrack.reserve(largeCoFactors.size());
    
    std::vector<double> uniLargeCoFacs;
    uniLargeCoFacs.reserve(largeCoFactors.size());
    
    for (auto myKey: largeCoFactors) {
        const auto trackIt = keepingTrack.find(myKey);
        
        if (trackIt != keepingTrack.end()) {
            coFactorIndexVec.push_back(trackIt->second);
        } else {
            keepingTrack[myKey] = coFactorInd;
            uniLargeCoFacs.push_back(static_cast<double>(myKey));
            coFactorIndexVec.push_back(coFactorInd);
            ++coFactorInd;
        }
    }
    
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
        nonTrivialFacs[j] = uniLargeCoFacs[i];

    std::vector<mpz_class> newTestInt(nRows);
    std::vector<std::uint8_t> mat(nRows * nCols, 0u);

    for (std::size_t r = 0, row = 0; r < nSmooth; ++r, row += nCols) {
        for (std::size_t j = 0; j < bigFacsTWO; ++j)
            ++mat[row + mapIndex[powsOfSmooths[r][j]]];
        
        for (auto it = powsOfSmooths[r].begin() + bigFacsTWO; it != powsOfSmooths[r].end(); ++it)
            ++mat[row + *it];
        
        newTestInt[r] = smoothInterval[r];
    }

    for (std::size_t i = 0, r = nSmooth,
         row = nCols * nSmooth; i < nPartial; ++i, ++r, row += nCols) {

        // Remap the powers corresponding to powers not contained in facBase
        for (std::size_t j = 0; j < bigFacsFOUR; ++j)
            ++mat[row + mapIndex[powsOfPartials[i][j]]];
        
        for (auto it = powsOfPartials[i].begin() + bigFacsFOUR; it != powsOfPartials[i].end(); ++it)
            ++mat[row + *it];
        
        // We must add one to account for negative values
        mat[row + nonTrivSize + coFactorIndexVec[i] + 1] += bigFacsTWO;
        newTestInt[r] = partialInterval[i];
    }
    
    SolutionSearch(mat, nRows, nCols, mpzNum, nonTrivialFacs,
                   newTestInt, factors, nThreads);

    if (bShowStats && cmp(factors[0], 0)) {
        const auto checkPoint2 = std::chrono::steady_clock::now();

        MakeStats(nSmooth + nPartial, nPolys, nSmooth,
                  nPartial, checkPoint2 - checkPoint0);

        RcppThread::Rcout << "\n" << std::endl;
    }
}
