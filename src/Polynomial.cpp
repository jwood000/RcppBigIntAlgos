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
                             std::vector<std::uint64_t> &largeCoFactorsBig, 
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
    
    std::vector<std::uint64_t> deleteLater;
    
    // First identify intersection
    for (const auto &pFac: partFactorsMap) {
        auto&& pFacBigIt = partFactorsMapBig.find(pFac.first);

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

Polynomial::Polynomial(std::size_t _facSize, std::size_t _vecMaxSize,
                       bool _bShowStats, const mpz_class &myNum) : 
            mpzFacSize(_facSize), SaPThresh(_facSize),
            facSize(_facSize), bShowStats(_bShowStats) {
    
    powsOfSmooths.reserve(_facSize);
    powsOfPartials.reserve(_facSize);
    myStart.resize(_facSize);
    myLogs.resize(_vecMaxSize);
    
    nPolys = 0;
    nPartial = 0;
    nSmooth = 0;
    
    if (bShowStats) {
        RcppThread::Rcout << "|      MPQS Time     | Complete | Polynomials |   Smooths"
                          << "  |  Partials  |\n|--------------------|----------|------"
                          << "-------|------------|------------|" << std::endl;
    }
}

Polynomial::Polynomial(std::size_t _facSize, std::size_t _vecMaxSize) : 
    SaPThresh(_facSize), facSize(_facSize), bShowStats(false) {

    myStart.resize(_facSize);
    myLogs.resize(_vecMaxSize);
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

        SinglePoly(SieveDist, facBase, LnFB, mpzFacBase[mpzFacSize - 1],
                   myNum, LowBound, theCut, TwiceLenB, vecMaxSize, strt, vecMaxStrt);
    }
}

void Polynomial::SieveListsInit(const std::vector<int> &facBase,
                                const std::vector<logType> &LnFB,
                                const std::vector<std::size_t> &SieveDist,
                                const mpz_class &firstSqrDiff, const mpz_class &VarA,
                                const mpz_class &VarB, std::size_t strt,
                                int LowBound, int vecMaxSize) {
    
    mpz_class Temp;
    
    for (std::size_t i = strt, facSize = facBase.size(); i < facSize; ++i) {
        auto&& myPrime = facBase[i];
        Temp = VarA % myPrime;
        const int AUtil = int_invert(Temp.get_si(), myPrime);
        
        mpz_ui_sub(Temp.get_mpz_t(), SieveDist[i], VarB.get_mpz_t());
        Temp *= AUtil;
        mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), myPrime);
        int myMin = Temp.get_si();
        
        mpz_ui_sub(Temp.get_mpz_t(), myPrime - SieveDist[i], VarB.get_mpz_t());
        Temp *= AUtil;
        mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), myPrime);
        int myMax = Temp.get_si();
        
        const int q = (LowBound % myPrime) + myPrime;
        mpz_mod_ui(Temp.get_mpz_t(), firstSqrDiff.get_mpz_t(), myPrime);
        
        if (myMin > myMax) {std::swap(myMin, myMax);}
        myStart[i].InitialSet(Temp.get_si(), q, myMin, myMax, myPrime);
        
        if (myPrime < vecMaxSize) {
            myStart[i].SmallSieve(myLogs, vecMaxSize, myPrime, LnFB[i]);
        } else {
            myStart[i].LargeSieve(myLogs, vecMaxSize, myPrime, LnFB[i]);
        }
    }
}

void Polynomial::SinglePoly(const std::vector<std::size_t> &SieveDist,
                            const std::vector<int> &facBase,
                            const std::vector<logType> &LnFB,
                            const mpz_class &NextPrime, const mpz_class &myNum,
                            int LowBound, logType theCut, int TwiceLenB,
                            int vecMaxSize, std::size_t strt, std::size_t vecMaxStrt) {
    
    mpz_class VarA, VarB, VarC, IntVal;
    TonelliShanksC(myNum, NextPrime, VarC);
    
    IntVal = VarC * 2u;
    mpz_invert(IntVal.get_mpz_t(), IntVal.get_mpz_t(), NextPrime.get_mpz_t());
    
    VarA = NextPrime * NextPrime;
    VarB = (IntVal * (myNum - VarC * VarC) + VarC) % VarA;
    VarC = (VarB * VarB - myNum) / VarA;
    IntVal = LowBound * (VarA * LowBound) + VarB * 2 * LowBound + VarC;
    std::fill(myLogs.begin(), myLogs.end(), 0);
    
    SieveListsInit(facBase, LnFB, SieveDist, IntVal,
                   VarA, VarB, strt, LowBound, vecMaxSize);
    
    for (int chunk = 0; chunk < TwiceLenB; chunk += vecMaxSize) {
        std::vector<int> largeLogs;
        
        for (int i = 0; i < vecMaxSize; ++i)
            if (myLogs[i] > theCut)
                largeLogs.push_back(i + chunk);
            
        for (const auto lrgLog: largeLogs) {
            std::vector<int> primeIndexVec;
            const int myIntVal = LowBound + lrgLog;
            IntVal = (VarA * myIntVal) * myIntVal + (VarB * myIntVal) * 2 + VarC;
            
            // Add the index referring to A^2.. (i.e. add it twice)
            primeIndexVec.insert(primeIndexVec.end(), 2, mpzFacSize);
            
            // If Negative, we push zero (i.e. the index referring to -1)
            if (sgn(IntVal) < 0) {
                IntVal = abs(IntVal);
                primeIndexVec.push_back(0);
            }
            
            for (std::size_t j = 0; j < strt; ++j) {
                while (mpz_divisible_ui_p(IntVal.get_mpz_t(), facBase[j])) {
                    IntVal /= facBase[j];
                    primeIndexVec.push_back(j + 1);
                }
            }
            
            for (std::uint32_t j = strt, facSize = facBase.size(),
                 ind = vecMaxSize - (lrgLog - chunk); j < facSize; ++j) {
                if (myStart[j].IsDivisible(facBase[j], ind)) {
                    do {
                        IntVal /= facBase[j];
                        primeIndexVec.push_back(j + 1);
                    } while (mpz_divisible_ui_p(IntVal.get_mpz_t(), facBase[j]));
                }
            }
            
            if (cmp(IntVal, 1u) == 0) {
                // Found a smooth number
                smoothInterval.push_back(VarA * myIntVal + VarB);
                powsOfSmooths.push_back(primeIndexVec);
            } else if (cmp(IntVal, Significand53) < 0) {
                const uint64_t myKey = static_cast<uint64_t>(IntVal.get_d());
                auto&& pFacIt = partFactorsMap.find(myKey);
                
                if (pFacIt != partFactorsMap.end()) {
                    largeCoFactors.push_back(myKey);
                    primeIndexVec.insert(primeIndexVec.begin(),
                                         pFacIt->second.cbegin(), pFacIt->second.cend());
                    
                    powsOfPartials.push_back(primeIndexVec);
                    auto&& intervalIt = partIntvlMap.find(myKey);
                    partialInterval.push_back((VarA * myIntVal + VarB) *
                        intervalIt->second);
                    
                    partFactorsMap.erase(pFacIt);
                    partIntvlMap.erase(intervalIt);
                } else {
                    partFactorsMap[myKey] = primeIndexVec;
                    partIntvlMap[myKey] = VarA * myIntVal + VarB;
                }
            }
        }
        
        if (chunk + vecMaxSize < TwiceLenB) {
            std::fill(myLogs.begin(), myLogs.end(), 0);
            
            for (std::size_t i = strt; i < vecMaxStrt; ++i)
                myStart[i].SmallSieve(myLogs, vecMaxSize, facBase[i], LnFB[i]);
            
            for (std::size_t i = vecMaxStrt, facSize = facBase.size(); i < facSize; ++i)
                myStart[i].LargeSieve(myLogs, vecMaxSize, facBase[i], LnFB[i]);
        }
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
            vecPoly.push_back(FromCpp14::make_unique<Polynomial>(facSize, vecMaxSize));
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

        SinglePoly(SieveDist, facBase, LnFB, NextPrime, myNum,
                   LowBound, theCut, TwiceLenB, vecMaxSize, strt, vecMaxStrt);

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
    
    std::unordered_map<std::uint64_t, std::size_t> keepingTrack;
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
