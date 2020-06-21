#include "SolutionSearch.h"
#include "Cpp14MakeUnique.h"
#include "StatsUtils.h"
#include "GrowMPZArray.h"
#include <unordered_map>
#include <fstream>

void QuadraticSieve(mpz_t myNum, mpz_t *const factors,
                    std::size_t nThreads, bool bShowStats) {
    
    const auto trialTest0 = std::chrono::steady_clock::now();
    const std::size_t digCount = mpz_sizeinbase(myNum, 10);
    const std::size_t bits = mpz_sizeinbase(myNum, 2);
    
    const double lognum = bits / log2(std::exp(1.0));
    const double sqrLogLog = std::sqrt(lognum * std::log(lognum));
    const double dblDigCount = digCount;
    
    // These values were obtained from "The Multiple Polynomial
    // Quadratic Sieve" by Robert D. Silverman
    // DigSize <- c(24, 30, 36, 42, 48, 54, 60, 66)
    // FBSize <- c(100, 200, 400, 900, 1200, 2000, 3000, 4500)
    // MSize <- c(5,20,35,100,125,250,350,500)
    // CounterSize <- c(10, 7, 5, 3.5, 2.25, 1.25, 1, 0.9)
    //
    // rawCoef <- round(unname(lm(FBSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept          x^1          x^2          x^3          x^4
    // 3637.0670996 -391.8275012   15.1541456   -0.2475566    0.0016806
    
    double fudge1 = -0.4;
    double LimB = std::exp((0.5 + fudge1) * sqrLogLog);
    
    const double dblMyTarget = std::ceil(-391.8275012 * dblDigCount + 15.1541456
                                             * std::pow(dblDigCount, 2.0) - 0.2475566
                                             * std::pow(dblDigCount, 3.0) + 0.0016806
                                             * std::pow(dblDigCount, 4.0) + 3637.0671);
    
    std::size_t myTarget = static_cast<std::size_t>(dblMyTarget);
    
    while (LimB < myTarget) {
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        fudge1 += 0.001;
    }
    
    const std::vector<std::size_t> facBase = getPrimesQuadRes(myNum, LimB, fudge1,
                                                              sqrLogLog, myTarget);
    
    // rawCoef <- round(unname(lm(MSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept           x^1           x^2           x^3           x^4
    // -213.1466450    23.3947361    -0.9494686     0.0165574    -0.0000767
    //
    // Note that the smallest that dblDigCount can be in order for dblLenB to be
    // valid is 22. This is no problem as we factor numbers less than 23 digits
    // with Pollard's Rho algorithm.
    
    const double dblLenB = std::ceil(23.394736 * dblDigCount - 0.9494686
                                         * std::pow(dblDigCount, 2.0) + 0.0165574
                                         * std::pow(dblDigCount, 3.0) - 0.0000767
                                         * std::pow(dblDigCount, 4.0) - 213.1466450);
    
    const std::size_t LenB = static_cast<std::size_t>(dblLenB) * 1000;
    const std::size_t facSize = facBase.size();
    const std::size_t DoubleLenB = 2 * LenB + 1;
    
    std::vector<int> myInterval(DoubleLenB);
    std::iota(myInterval.begin(), myInterval.end(), -1 * static_cast<int>(LenB));
    
    mpz_t TS[10];
    
    for (std::size_t i = 0; i < 10; ++i)
        mpz_init(TS[i]);
    
    const std::vector<std::size_t> SieveDist = setSieveDist(myNum, TS,
                                                            facBase, facSize);
    
    // CounterSize <- c(10, 7, 5, 3.5, 2.25, 1.25, 1, 0.9)
    // rawCoef <- round(unname(lm(CounterSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //  intercept        x^1        x^2        x^3        x^4 
    // 40.1734307 -2.2578629  0.0571943 -0.0007403  0.0000039
    
    double dblFacMultiplier = std::max(std::ceil(-2.2578629 * dblDigCount + 0.0571943
                                                     * std::pow(dblDigCount, 2.0) - 0.0007403
                                                     * std::pow(dblDigCount, 3.0) + 0.0000039
                                                     * std::pow(dblDigCount, 4.0) + 40.1734307), 1.0);
    
    const std::size_t mpzContainerSize = static_cast<double>(facSize) * dblFacMultiplier;
    
    // smoothInterval will be the values that are associated
    // with smooth numbers. We make it 1.5x the size of facSize
    // in order to guarantee we will have enough space
    auto smoothInterval = FromCpp14::make_unique<mpz_t[]>(mpzContainerSize);

    // Similarly, we create an array that will house the
    // interval entries corresponding to the partially smooth facs
    auto partialInterval = FromCpp14::make_unique<mpz_t[]>(mpzContainerSize);

    // This array will store the large cofactors that will
    // later be added to mpzFacBase
    auto largeCoFactors = FromCpp14::make_unique<mpz_t[]>(mpzContainerSize);

    for (std::size_t i = 0; i < mpzContainerSize; ++i) {
        mpz_init(smoothInterval[i]);
        mpz_init(partialInterval[i]);
        mpz_init(largeCoFactors[i]);
    }
    
    std::size_t mpzFSzLIMIT = facSize * 10;
    
    // This array will be passed to solutionSeach.
    auto mpzFacBase = FromCpp14::make_unique<mpz_t[]>(mpzFSzLIMIT);
    
    for (std::size_t i = 0; i < mpzFSzLIMIT; ++i)
        mpz_init(mpzFacBase[i]);
    
    mpz_t temp;
    mpz_init(temp);
    mpz_mul_ui(temp, myNum, 2);
    mpz_sqrt(temp, temp);
    mpz_mul_ui(temp, temp, LenB);
    
    const double fudge2 = (digCount < 45) ? 1.410 :
                          (digCount < 50) ? 1.440 :
                          (digCount < 60) ? 1.500 : 
                          (digCount < 65) ? 1.515 : 1.530;
    
    const double theCut = fudge2 * mpz_sizeinbase(temp, 10);
    
    std::vector<double> LnFB(facSize);
    std::vector<double> myLogs(DoubleLenB, 0.0);
    LnFB.reserve(facSize);
    
    for (std::size_t i = 0; i < facSize; ++i) {
        mpz_set_ui(mpzFacBase[i], facBase[i]);
        LnFB[i] = std::log(static_cast<double>(facBase[i]));
    }
    
    mpz_t sqrtInt, firstInt;
    mpz_init(sqrtInt); mpz_init(firstInt);
    mpz_sqrt(sqrtInt, myNum);
    mpz_mul_ui(firstInt, sqrtInt, LenB);
    
    const std::size_t minPrime = static_cast<std::size_t>(mpz_sizeinbase(firstInt, 10) * 2);
    
    mpz_pow_ui(temp, firstInt, 2u);
    mpz_sub(temp, temp, myNum);
    
    sieveLists(facSize, facBase, DoubleLenB, temp, LnFB,
               myLogs, minPrime, SieveDist, firstInt);
    
    std::vector<std::size_t> largeLogs;

    for (std::size_t i = 0; i < DoubleLenB; ++i)
        if (myLogs[i] > theCut)
            largeLogs.push_back(i);
    
    // The map below is for almost smooth numbers, i.e. numbers
    // that are almost completely factored by our factor base
    hash64vec partFactorsMap;
    hash64mpz_t partIntvlMap;
    
    vec2dsize_t powsOfPartials;
    hash64size_t keepingTrack;
    std::vector<std::size_t> coFactorIndexVec;
    
    std::size_t nPartial = 0u;
    std::size_t coFactorInd = 0u;
    std::size_t nSmooth = 0u;
    
    mpz_t largeInt, intVal;
    mpz_init(largeInt); mpz_init(intVal);
    
    vec2dsize_t powsOfSmooths;
    powsOfSmooths.reserve(facSize);
    
    for (auto lrglog: largeLogs) {
        std::vector<std::size_t> primeIndexVec;
        
        mpz_add_ui(largeInt, firstInt, lrglog);
        mpz_pow_ui(temp, largeInt, 2u);
        mpz_sub(intVal, temp, myNum);

        if (mpz_sgn(intVal) < 0) {
            mpz_abs(intVal, intVal);
            primeIndexVec.push_back(0u);
        }

        for (std::size_t i = 0; i < facSize; ++i) {
            while (mpz_divisible_ui_p(intVal, facBase[i])) {
                mpz_divexact_ui(intVal, intVal, facBase[i]);
                primeIndexVec.push_back(i + 1);
            }
        }

        if (mpz_cmp_ui(intVal, 1) == 0) {
            // Found a smooth number
            mpz_set(smoothInterval[nSmooth], largeInt);
            powsOfSmooths.push_back(primeIndexVec);
            ++nSmooth;
        } else if (mpz_cmp_d(intVal, Significand53) < 0) {
            const uint64_t myKey = static_cast<uint64_t>(mpz_get_d(intVal));
            const auto pFacIt = partFactorsMap.find(myKey);

            if (pFacIt != partFactorsMap.end()) {
                const auto trackIt = keepingTrack.find(myKey);

                if (trackIt != keepingTrack.end()) {
                    coFactorIndexVec.push_back(trackIt->second);
                } else {
                    keepingTrack[myKey] = coFactorInd;
                    mpz_set(largeCoFactors[coFactorInd], intVal);
                    coFactorIndexVec.push_back(coFactorInd);
                    ++coFactorInd;
                }

                for (const auto p: pFacIt->second)
                    primeIndexVec.push_back(p);

                powsOfPartials.push_back(primeIndexVec);
                const auto intervalIt = partIntvlMap.find(myKey);

                mpz_mul(partialInterval[nPartial],
                        largeInt, intervalIt->second);

                partFactorsMap.erase(pFacIt);
                partIntvlMap.erase(intervalIt);
                ++nPartial;
            } else {
                partFactorsMap[myKey] = primeIndexVec;
                mpz_init(partIntvlMap[myKey]);
                mpz_set(partIntvlMap[myKey], largeInt);
            }
        }
    }
    
    const auto trialTest1 = std::chrono::steady_clock::now();

    mpz_t A, B, C, Atemp, Atemp2, Btemp, lowBound;
    mpz_init(A); mpz_init(B); mpz_init(C);
    mpz_init(Atemp); mpz_init(Atemp2); mpz_init(Btemp);

    mpz_init(lowBound);
    mpz_set_si(lowBound, myInterval.front());

    mpz_mul_2exp(Atemp, myNum, 1);
    mpz_sqrt(Atemp, Atemp);
    mpz_div_ui(Atemp, Atemp, LenB);
    mpz_sqrt(Atemp, Atemp);

    if (mpz_cmp_ui(Atemp, facBase.back()) < 0)
        mpz_set_ui(Atemp, facBase.back());

    std::size_t nPolys = 0;
    std::size_t extraFacs = 0;
    std::size_t mpzFacSize = facSize;
    
    auto showStatsTime = (trialTest1 - trialTest0);
    
    if (bShowStats) {
        const std::size_t onePercentComplete = facSize / 100;
        updateStatTime(nSmooth + nPartial, onePercentComplete, showStatsTime, showStatsTime);
    }
    
    auto checkPoint0 = std::chrono::steady_clock::now();
    auto checkPoint1 = checkPoint0;
    auto checkPoint2 = checkPoint0;
    
    std::size_t lastLim = 0;
    std::size_t currLim = nSmooth + nPartial;
    std::vector<std::size_t> polySieveD(facSize * 2, 0u);
    
    if (bShowStats) {
        auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(myNum, 10) + 2);
        mpz_get_str(buffer.get(), 10, myNum);
        const std::string strMyNum = buffer.get();
        
        RcppThread::Rcout << "\nSummary Statistics for Factoring:\n" << "    "
                          << strMyNum << "\n\n"
                          << "|        Time        | Complete | Polynomials |   Smooths"
                          << "  |  Partials  |\n";
        RcppThread::Rcout << "|--------------------|----------|-------------|---------"
                          << "---|------------|" << std::endl;
    }
    
    while (mpz_cmp_ui(factors[0], 0) == 0) {
        const std::size_t loopLimit = facSize + extraFacs;
        
        // Find enough smooth numbers to guarantee a non-trivial solution
        while (currLim <= loopLimit) {
            bool LegendreTest = true;
            
            while (LegendreTest) {
                mpz_nextprime(Atemp, Atemp);
                
                if (mpz_legendre(myNum, Atemp) == 1)
                    LegendreTest = false;
            }
            
            if (mpzFacSize >= mpzFSzLIMIT) {
                mpzFSzLIMIT <<= 1;
                Grow(mpzFacBase, mpzFacSize, mpzFSzLIMIT);
            }
            
            mpz_set(mpzFacBase[mpzFacSize], Atemp);
            ++mpzFacSize;
             
            SinglePoly(myInterval, polySieveD, smoothInterval.get(), SieveDist, TS, facBase,
                       mpzFacBase.get(), LnFB, myLogs, largeCoFactors.get(), partialInterval.get(),
                       powsOfSmooths, powsOfPartials, partFactorsMap, partIntvlMap, keepingTrack,
                       coFactorIndexVec, nPartial, nSmooth, coFactorInd, intVal, myNum, Atemp,
                       Btemp, temp, A, B, C, Atemp2, lowBound, minPrime, theCut, DoubleLenB,
                       mpzFacSize, facSize);
            
            ++nPolys;
            currLim = nSmooth + nPartial;
            
            const auto checkPoint3 = std::chrono::steady_clock::now();
            
            if ((checkPoint3 - checkPoint1) > checkInterTime) {
                // Check for user interrupt and udpate timepoint
                RcppThread::checkUserInterrupt();
                checkPoint1 = std::chrono::steady_clock::now();
            }
            
            if (bShowStats && (checkPoint3 - checkPoint2) > showStatsTime) {
                makeStats(currLim, loopLimit, nPolys, nSmooth,
                          nPartial, lastLim, checkPoint3 - checkPoint0);
                
                checkPoint2 = std::chrono::steady_clock::now();
                updateStatTime(currLim, facSize, checkPoint3 - checkPoint0, showStatsTime);
            }
        }

        const std::size_t matRow = nSmooth + nPartial;
        const std::size_t matWidth = coFactorInd + mpzFacSize + 1;

        auto newTestInt = FromCpp14::make_unique<mpz_t[]>(matRow);
        std::vector<std::uint8_t> newMat(matRow * matWidth, 0);

        for (std::size_t r = 0, row = 0; r < nSmooth; row += matWidth, ++r) {
            for (const auto p: powsOfSmooths[r])
                ++newMat[row + p];

            mpz_init_set(newTestInt[r], smoothInterval[r]);
        }

        if ((coFactorInd + mpzFacSize) >= mpzFSzLIMIT) {
            mpzFSzLIMIT <<= 1;
            Grow(mpzFacBase, mpzFacSize, mpzFSzLIMIT);
        }

        for (std::size_t i = 0, j = mpzFacSize; i < coFactorInd; ++i, ++j)
            mpz_set(mpzFacBase[j], largeCoFactors[i]);

        for (std::size_t i = 0, r = nSmooth,
             row = matWidth * nSmooth; i < nPartial; ++i, ++r, row += matWidth) {

            for (const auto p: powsOfPartials[i])
                ++newMat[row + p];

            newMat[row + mpzFacSize + 1 + coFactorIndexVec[i]] = 2u;
            mpz_init_set(newTestInt[r], partialInterval[i]);
        }

        solutionSearch(newMat, matRow, matWidth, myNum,
                       mpzFacBase.get(), newTestInt.get(), factors);

        for (std::size_t i = 0; i < matRow; ++i)
            mpz_clear(newTestInt[i]);

        newTestInt.reset();
        extraFacs += 5;

        if (bShowStats && mpz_cmp_ui(factors[0], 0)) {
            const auto checkPoint3 = std::chrono::steady_clock::now();

            makeStats(loopLimit, loopLimit, nPolys, nSmooth,
                      nPartial, lastLim, checkPoint3 - checkPoint0);

            RcppThread::Rcout << "\n" << std::endl;
        }
    }
    
    for (std::size_t i = 0; i < mpzContainerSize; ++i) {
        mpz_clear(smoothInterval[i]);
        mpz_clear(partialInterval[i]);
        mpz_clear(largeCoFactors[i]);
    }
    
    for (std::size_t i = 0; i < mpzFSzLIMIT; ++i)
        mpz_clear(mpzFacBase[i]);
    
    mpzFacBase.reset();
    smoothInterval.reset();
    partialInterval.reset();
    largeCoFactors.reset();
    
    for (auto v: partIntvlMap)
        mpz_clear(v.second);
    
    partIntvlMap.clear();
    partFactorsMap.clear();
    keepingTrack.clear();
    coFactorIndexVec.clear();
    
    mpz_clear(temp); mpz_clear(sqrtInt); mpz_clear(A);
    mpz_clear(C); mpz_clear(Atemp); mpz_clear(largeInt);
    mpz_clear(lowBound); mpz_clear(B); mpz_clear(firstInt);
    mpz_clear(Atemp2); mpz_clear(intVal); mpz_clear(Btemp);

    for (std::size_t i = 0; i < 10; ++i)
        mpz_clear(TS[i]);
}
