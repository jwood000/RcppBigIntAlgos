#include "Polynomial.h"
#include "SolutionSearch.h"
#include <set>

// Determined empirically
constexpr std::size_t MinPolysPerThrd = 40u;
constexpr std::size_t bigFacsTWO = 2u;
constexpr std::size_t bigFacsFOUR = 4u;

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
        const auto pFacBigIt = partFactorsMapBig.find(pFac.first);

        if (pFacBigIt != partFactorsMapBig.end()) {
            largeCoFactorsBig.push_back(pFac.first);
            pFacBigIt->second.insert(pFacBigIt->second.end(),
                                     pFac.second.cbegin(), pFac.second.cend());
            
            powsOfPartialsBig.push_back(pFacBigIt->second);
            partialIntervalBig.push_back(partIntvlMap[pFac.first]
                                             * partIntvlMapBig[pFac.first]);
            deleteLater.push_back(pFac.first);
        }
    }
    
    for (const auto myKey: deleteLater) {
        partFactorsMap.erase(myKey);
        partFactorsMapBig.erase(myKey);
        partIntvlMap.erase(myKey);
        partIntvlMapBig.erase(myKey);
    }
    
    partFactorsMapBig.insert(partFactorsMap.cbegin(), partFactorsMap.cend());
    partIntvlMapBig.insert(partIntvlMap.cbegin(), partIntvlMap.cend());
}

Polynomial::Polynomial(std::size_t _facSize, bool _bShowStats, const mpz_class &myNum) : 
            mpzFacSize(_facSize), SaPThresh(_facSize),
            facSize(_facSize), bShowStats(_bShowStats) {
    
    powsOfSmooths.reserve(_facSize);
    powsOfPartials.reserve(_facSize);
    myStart.assign(_facSize * 2, 0);
    
    nPolys = 0;
    nPartial = 0;
    nSmooth = 0;
    
    if (bShowStats) {
        RcppThread::Rcout << "|      MPQS Time     | Complete | Polynomials |   Smooths"
                          << "  |  Partials  |\n|--------------------|----------|------"
                          << "-------|------------|------------|" << std::endl;
    }
}

Polynomial::Polynomial(std::size_t _facSize) : 
    SaPThresh(_facSize), facSize(_facSize), bShowStats(false) {

    myStart.assign(_facSize * 2, 0);
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
                            const std::vector<int> &facBase,
                            const std::vector<logType> &LnFB,
                            const std::vector<mpz_class> &mpzFacBase,
                            const mpz_class &myNum, int LowBound,
                            logType theCut, int TwiceLenB, int vecMaxSize,
                            std::size_t strt, std::size_t vecMaxStrt,
                            std::size_t polyLimit) {
    
    for (std::size_t poly = 0; poly < polyLimit; ++poly) {
        ++mpzFacSize;

        MPQS::SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                         myStart, partFactorsMap, partIntvlMap, smoothInterval,
                         largeCoFactors, partialInterval, mpzFacBase[mpzFacSize - 1],
                         myNum, LowBound, theCut, TwiceLenB, mpzFacSize, vecMaxSize,
                         strt, vecMaxStrt);
    }
}

void Polynomial::InitialParSieve(const std::vector<std::size_t> &SieveDist,
                                 const std::vector<int> &facBase,
                                 const std::vector<logType> &LnFB,
                                 std::vector<mpz_class> &mpzFacBase,
                                 mpz_class &NextPrime, const mpz_class &myNum,
                                 int LowBound, logType theCut, int TwiceLenB,
                                 int vecMaxSize, std::size_t strt,
                                 std::size_t vecMaxStrt, typeTimePoint checkPoint0) {
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    
    auto showStatsTime = (checkPoint1 - checkPoint0);
    GetNPrimes(mpzFacBase, NextPrime, myNum, MinPolysPerThrd);
    
    SievePolys(SieveDist, facBase, LnFB, mpzFacBase, myNum, LowBound,
               theCut, TwiceLenB, vecMaxSize, strt, vecMaxStrt, MinPolysPerThrd);
    
    nPolys = MinPolysPerThrd;
    nSmooth = smoothInterval.size();
    nPartial = partialInterval.size();
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

template <typename typeTimeDiff>
void SetThreadsPolys(std::size_t currLim, std::size_t SaPThresh, 
                     std::size_t &nThreads, std::size_t maxThreads,
                     std::size_t &polysPerThread, std::size_t nPolys,
                     typeTimeDiff totalTime) {
    
    const double percentComplete = (0.01 + 100.0 * currLim) / SaPThresh;
    const double percentRemaining = 100 - percentComplete;
    
    // This was obtained experimentally. The idea is that if we have
    // only completed a small percentage, we really don't know how
    // quickly the growth can be. Thus, we penalize these more
    // heavily. With this setup, when percentRemaining = 100, fudge
    // is ~ 3.178 and when percentRemaining = 0, fudge ~= 1.386
    const double fudge = std::log(percentRemaining * 0.2 + 4);
    
    const std::size_t nPolysRemaining = percentRemaining * 
        (static_cast<double>(nPolys) / (fudge * percentComplete));
    
    nThreads = std::min(maxThreads, nPolysRemaining / MinPolysPerThrd);
    
    if (nThreads) {
        const std::size_t maxPolys = nPolysRemaining / nThreads;
        const auto timePerPoly = nThreads * totalTime / nPolys;
        const auto maxTime = maxPolys * timePerPoly;
        
        if (maxTime > (2 * fifteenSeconds)) {
            polysPerThread = 2 * fifteenSeconds / timePerPoly;
        } else if (maxTime > fifteenSeconds) {
            polysPerThread = fifteenSeconds / timePerPoly;
        } else {
            polysPerThread =  maxPolys;
        }
    }
}

void Polynomial::FactorParallel(const std::vector<std::size_t> &SieveDist,
                                const std::vector<int> &facBase, 
                                const std::vector<logType> &LnFB,
                                std::vector<mpz_class> &mpzFacBase,
                                mpz_class &NextPrime, const mpz_class &myNum,
                                int LowBound, logType theCut, int TwiceLenB,
                                int vecMaxSize, std::size_t strt,
                                std::size_t vecMaxStrt, typeTimePoint checkPoint0,
                                std::size_t nThreads) {
    
    auto checkPoint1 = std::chrono::steady_clock::now();
    auto checkPoint2 = checkPoint1;
    
    this->InitialParSieve(SieveDist, facBase, LnFB, mpzFacBase, NextPrime,
                          myNum, LowBound, theCut, TwiceLenB, vecMaxSize,
                          strt, vecMaxStrt, checkPoint0);
    
    auto showStatsTime = std::chrono::steady_clock::now() - checkPoint0;
    std::size_t polysPerThread = MinPolysPerThrd;
    const std::size_t maxThreads = nThreads;
    
    SetThreadsPolys(nSmooth + nPartial, SaPThresh, nThreads,
                    maxThreads, polysPerThread, nPolys, showStatsTime);
    
    while ((nSmooth + nPartial) <= SaPThresh && nThreads > 1) {
        std::vector<std::unique_ptr<Polynomial>> vecPoly;
        std::vector<std::thread> myThreads;
        
        const std::size_t startStep = mpzFacBase.size();
        NextPrime = mpzFacBase.back();

        GetNPrimes(mpzFacBase, NextPrime, myNum, polysPerThread * nThreads);
        mpzFacSize = mpzFacBase.size();
        
        for (std::size_t i = 0, step = startStep; i < nThreads; ++i, step += polysPerThread) {
            vecPoly.push_back(FromCpp14::make_unique<Polynomial>(facSize));
            vecPoly[i]->SetMpzFacSize(step);

            myThreads.emplace_back(&Polynomial::SievePolys, vecPoly[i].get(),
                                   std::cref(SieveDist), std::cref(facBase), std::cref(LnFB),
                                   std::cref(mpzFacBase), std::cref(myNum), LowBound, theCut,
                                   TwiceLenB, vecMaxSize, strt, vecMaxStrt, polysPerThread);
        }

        for (auto &thr: myThreads)
            thr.join();
        
        for (std::size_t i = 0; i < nThreads; ++i) {
            vecPoly[i]->MergeMaster(powsOfSmooths, powsOfPartials,
                                    partFactorsMap, partIntvlMap, smoothInterval,
                                    largeCoFactors, partialInterval);
        }
        
        nSmooth = smoothInterval.size();
        nPartial = partialInterval.size();
        nPolys += (nThreads * polysPerThread);
        const auto checkPoint3 = std::chrono::steady_clock::now();
        
        SetThreadsPolys(nSmooth + nPartial, SaPThresh, nThreads, maxThreads,
                        polysPerThread, nPolys, checkPoint3 - checkPoint0);
        
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

void Polynomial::FactorSerial(const std::vector<std::size_t> &SieveDist,
                              const std::vector<int> &facBase,
                              const std::vector<logType> &LnFB,
                              std::vector<mpz_class> &mpzFacBase,
                              mpz_class &NextPrime, const mpz_class &myNum,
                              int LowBound, logType theCut, int TwiceLenB,
                              int vecMaxSize, std::size_t strt,
                              std::size_t vecMaxStrt, typeTimePoint checkPoint0) {
    
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

        MPQS::SinglePoly(SieveDist, facBase, LnFB, powsOfSmooths, powsOfPartials,
                         myStart, partFactorsMap, partIntvlMap, smoothInterval,
                         largeCoFactors, partialInterval, NextPrime, myNum,
                         LowBound, theCut, TwiceLenB, mpzFacSize, vecMaxSize,
                         strt, vecMaxStrt);

        ++nPolys;
        nSmooth = smoothInterval.size();
        nPartial = partialInterval.size();
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
    
    if (bShowStats) {
        const auto checkPoint2 = std::chrono::steady_clock::now();
        
        MakeStats(nSmooth + nPartial, nPolys, nSmooth,
                  nPartial, checkPoint2 - checkPoint0);
        
        RcppThread::Rcout << "\n" << std::endl;
    }
    
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
                   newTestInt, factors, nThreads, bShowStats);
}
