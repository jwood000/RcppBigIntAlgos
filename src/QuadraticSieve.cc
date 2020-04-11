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
    
    mpz_t currP, nextP, resTest, CP1;
    mpz_init(currP); mpz_init(nextP);
    mpz_init(CP1); mpz_init(resTest);
    
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
    
    const std::size_t myTarget = static_cast<std::size_t>(dblMyTarget);
    
    while (LimB < myTarget) {
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        fudge1 += 0.001;
    }
    
    std::vector<std::size_t> facBase = getPrimesQuadRes(myNum, LimB);
    
    while (facBase.size() < myTarget) {
        fudge1 += 0.005;
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        mpz_set_ui(currP, facBase.back());
        mpz_nextprime(nextP, currP);
        
        while (mpz_cmp_ui(nextP, LimB) < 0) {
            mpz_set(currP, nextP);
            mpz_nextprime(nextP, currP);
            mpz_sub_ui(CP1, currP, 1);
            mpz_div_2exp(CP1, CP1, 1);
            mpz_powm(resTest,myNum,CP1,currP);
            
            if (mpz_cmp_ui(resTest, 1) == 0)
                facBase.push_back(mpz_get_ui(currP));
        }
    }
    
    mpz_clear(currP); mpz_clear(nextP);
    mpz_clear(CP1); mpz_clear(resTest);
    
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
    const int64_t Lower = -1 * static_cast<int64_t>(LenB);
    const int64_t Upper = LenB;
    const std::size_t LenB2 = 2 * LenB + 1;
    
    std::vector<int64_t> myInterval(LenB2);
    std::iota(myInterval.begin(), myInterval.end(), Lower);
    
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

    auto largeInterval = FromCpp14::make_unique<mpz_t[]>(LenB2);
    auto sqrDiff = FromCpp14::make_unique<mpz_t[]>(LenB2);

    for (std::size_t i = 0; i < LenB2; ++i) {
        mpz_init(largeInterval[i]);
        mpz_init(sqrDiff[i]);
    }
    
    mpz_t sqrtInt;
    mpz_init(sqrtInt);
    mpz_sqrt(sqrtInt, myNum);
    mpz_sub_ui(largeInterval[0], sqrtInt, Upper);

    for (std::size_t i = 1; i < LenB2; ++i)
        mpz_add_ui(largeInterval[i], largeInterval[i - 1], 1u);

    mpz_t temp;
    mpz_init(temp);

    for (std::size_t i = 0; i < LenB2; ++i) {
        mpz_pow_ui(temp, largeInterval[i], 2u);
        mpz_sub(sqrDiff[i], temp, myNum);
    }

    mpz_mul_ui(temp, myNum, 2);
    mpz_sqrt(temp, temp);
    mpz_mul_ui(temp, temp, Upper);
    
    // rawCoef <- unname(lm(c(0.7, 0.9, 1.05, 1.1, 1.2, 1.5) ~ c(25, seq(34, 49, 5), 65))$coefficients)
    // names(rawCoef) <- c("intercept", "x^1")
    // rawCoef
    //  intercept        x^1 
    // 0.23420680 0.01970609
    const double fudge2 = 0.23420680 + digCount * 0.01970609;
    const double theCut = fudge2 * mpz_sizeinbase(temp, 10);
    const double partFacLimit = static_cast<double>((fudge2 * 200) * facBase.back());
    
    Rcpp::Rcout << partFacLimit << "\n";
    
    std::vector<double> LnFB(facSize);
    std::vector<double> myLogs(LenB2, 0.0);
    LnFB.reserve(facSize);

    for (std::size_t i = 0; i < facSize; ++i) {
        mpz_set_ui(mpzFacBase[i], facBase[i]);
        LnFB[i] = std::log(static_cast<double>(facBase[i]));
    }

    mpz_mul_ui(temp, sqrtInt, Upper);
    const std::size_t minPrime = static_cast<std::size_t>(mpz_sizeinbase(temp, 10) * 2);
    
    sieveLists(facSize, facBase, LenB2, sqrDiff.get(), LnFB,
               myLogs, minPrime, SieveDist, largeInterval[0]);
    
    std::vector<std::size_t> largeLogs;

    for (std::size_t i = 0; i < LenB2; ++i)
        if (myLogs[i] > theCut)
            largeLogs.push_back(i);
    
    std::size_t largeLogsSize = largeLogs.size();
    const std::size_t colWidth = facSize + 1;
    std::vector<uint8_t> myMat(largeLogsSize * colWidth, 0u);
    
    // Used to identify the row of the matrix that
    // contains a smooth number
    std::vector<std::size_t> smoothFacsRow;
    
    // The map below is for almost smooth numbers, i.e. numbers
    // that are almost completely factored by our factor base
    std::unordered_map<uint64_t, std::vector<std::size_t>> partFactorsMap;
    std::unordered_map<uint64_t, mpz_t> partIntvlMap;
    
    std::vector<std::vector<std::size_t>> powsOfPartials;
    std::unordered_map<uint64_t, std::size_t> keepingTrack;
    std::vector<std::size_t> coFactorIndexVec;
    
    std::size_t nPartial = 0u;
    std::size_t coFactorInd = 0u;
    std::size_t nSmooth = 0u;
    
    // std::ofstream partialFile;
    // std::ofstream smoothFile;
    // partialFile.open("cutFudge2Partial.csv");
    
    mpz_t myMax, realMax;
    mpz_init(myMax); mpz_init(realMax);
    mpz_set_ui(myMax, 0u);
    mpz_set_ui(realMax, 0u);
    
    // std::size_t count = 0;
    // 
    // std::size_t numBreak = 56865667;
    // std::ofstream valFile;
    // valFile.open("partialCofactors.csv");
    
    if (largeLogsSize > 0) {
        for (std::size_t j = 0; j < largeLogsSize; ++j) {
            std::vector<std::size_t> primeIndexVec;

            if (mpz_sgn(sqrDiff[largeLogs[j]]) < 0) {
                myMat[j * colWidth] = 1;
                mpz_abs(sqrDiff[largeLogs[j]], sqrDiff[largeLogs[j]]);
                primeIndexVec.push_back(0u);
            }
            
            for (std::size_t i = 0; i < facSize; ++i) {
                while (mpz_divisible_ui_p(sqrDiff[largeLogs[j]], facBase[i])) {
                    mpz_divexact_ui(sqrDiff[largeLogs[j]],
                                    sqrDiff[largeLogs[j]], facBase[i]);
                    ++myMat[j * colWidth + i + 1];
                    primeIndexVec.push_back(i + 1);
                }
            }

            if (mpz_cmp_ui(sqrDiff[largeLogs[j]], 1) == 0) {
                // Found a smooth number
                mpz_set(smoothInterval[nSmooth], largeInterval[largeLogs[j]]);
                
                // auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(smoothInterval[nSmooth], 10) + 2);
                // mpz_get_str(buffer.get(), 10, smoothInterval[nSmooth]);
                // std::string strMyNum = buffer.get();
                // 
                // smoothFile << strMyNum << "\n";
                
                smoothFacsRow.push_back(j * colWidth);
                ++nSmooth;
            } else if (mpz_cmp_d(sqrDiff[largeLogs[j]], partFacLimit) < 0) {
                
                const uint64_t myKey = static_cast<uint64_t>(mpz_get_d(sqrDiff[largeLogs[j]]));
                const auto pFacIt = partFactorsMap.find(myKey);
                
                // if (mpz_cmp(sqrDiff[largeLogs[j]], myMax) > 0) {
                //     mpz_set(myMax, sqrDiff[largeLogs[j]]);
                // }
                // 
                // if (mpz_cmp_ui(sqrDiff[largeLogs[j]], numBreak) > 0) {
                //     ++count;
                // }

                if (pFacIt != partFactorsMap.end()) {
                    const auto trackIt = keepingTrack.find(myKey);
                    
                    // if (mpz_cmp(sqrDiff[largeLogs[j]], realMax) > 0) {
                    //     mpz_set(realMax, sqrDiff[largeLogs[j]]);
                    // }
                    // double myVal = mpz_get_d(sqrDiff[largeLogs[j]]);
                    // valFile << myVal << "\n";

                    if (trackIt != keepingTrack.end()) {
                        coFactorIndexVec.push_back(trackIt->second);
                    } else {
                        keepingTrack[myKey] = coFactorInd;
                        mpz_set(largeCoFactors[coFactorInd], sqrDiff[largeLogs[j]]);
                        coFactorIndexVec.push_back(coFactorInd);
                        ++coFactorInd;
                    }

                    for (const auto p: pFacIt->second)
                        primeIndexVec.push_back(p);

                    powsOfPartials.push_back(primeIndexVec);
                    const auto intervalIt = partIntvlMap.find(myKey);
                    
                    // auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(intervalIt->second, 10) + 2);
                    // mpz_get_str(buffer.get(), 10, intervalIt->second);
                    // std::string strMyNum = buffer.get();
                    // 
                    // partialFile << strMyNum << "\n";
                    // 
                    // buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(largeInterval[largeLogs[j]], 10) + 2);
                    // mpz_get_str(buffer.get(), 10, largeInterval[largeLogs[j]]);
                    // strMyNum = buffer.get();
                    // 
                    // partialFile << strMyNum << "\n";
                    
                    mpz_mul(partialInterval[nPartial],
                            largeInterval[largeLogs[j]],
                            intervalIt->second);

                    partFactorsMap.erase(pFacIt);
                    partIntvlMap.erase(intervalIt);
                    ++nPartial;
                } else {
                    partFactorsMap[myKey] = primeIndexVec;
                    mpz_init(partIntvlMap[myKey]);
                    mpz_set(partIntvlMap[myKey], largeInterval[largeLogs[j]]);
                }
            }
        }
    }
    
    const auto trialTest1 = std::chrono::steady_clock::now();

    for (std::size_t i = 0; i < LenB2; ++i)
        mpz_clear(largeInterval[i]);

    largeInterval.reset();

    std::size_t lenM = smoothFacsRow.size();
    std::vector<std::vector<uint8_t>> powsOfSmooths;
    std::vector<uint8_t> tempMat(lenM * colWidth, 0u);

    for (std::size_t i = 0, rowInd = 0; i < lenM; ++i, rowInd += colWidth)
        for (std::size_t j = 0; j <= facSize; ++j)
            tempMat[rowInd + j] = myMat[smoothFacsRow[i] + j];

    powsOfSmooths.push_back(tempMat);

    mpz_t A, B, C, Atemp, Atemp2, Btemp, lowBound;
    mpz_init(A); mpz_init(B); mpz_init(C);
    mpz_init(Atemp); mpz_init(Atemp2); mpz_init(Btemp);

    mpz_init(lowBound);
    mpz_set_si(lowBound, myInterval.front());

    mpz_mul_2exp(Atemp, myNum, 1);
    mpz_sqrt(Atemp, Atemp);
    mpz_div_ui(Atemp, Atemp, Upper);
    mpz_sqrt(Atemp, Atemp);

    if (mpz_cmp_ui(Atemp, facBase.back()) < 0)
        mpz_set_ui(Atemp, facBase.back());

    std::size_t nPolys = 0;
    std::size_t extraFacs = 0;
    std::size_t mpzFacSize = facSize;
    std::vector<double> myIntervalSqrd(LenB2);

    for (std::size_t i = 0; i < LenB2; ++i)
        myIntervalSqrd[i] = static_cast<double>(myInterval[i] * myInterval[i]);
    
    auto showStatsTime = (trialTest1 - trialTest0);
    const std::size_t onePercentComplete = facSize / 100;
    
    if (bShowStats) {
        updateStatTime(nSmooth + nPartial,
                       onePercentComplete, showStatsTime, showStatsTime);
    }
    
    auto checkPoint0 = std::chrono::steady_clock::now();
    auto checkPoint1 = checkPoint0;
    auto checkPoint2 = checkPoint0;
    
    std::size_t lastLim = 0;
    std::size_t currLim = nSmooth + nPartial;
    
    std::vector<std::size_t> powsFacInd;
    powsFacInd.push_back(0);
    
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
    
    // Rcpp::Rcout << theCut << std::endl;
    
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
            
            TonelliShanksC(myNum, Atemp, TS);

            if (mpz_cmp(TS[1], TS[2]) > 0) {
                mpz_set(Btemp, TS[1]);
            } else {
                mpz_set(Btemp, TS[2]);
            }

            mpz_mul_2exp(temp, Btemp, 1);
            mpz_invert(temp, temp, Atemp);
            mpz_pow_ui(B, Btemp, 2);

            mpz_sub(B, myNum, B);
            mpz_mul(B, B, temp);
            mpz_add(B, B, Btemp);

            mpz_pow_ui(A, Atemp, 2);
            mpz_mod(B, B, A);

            mpz_pow_ui(C, B, 2);
            mpz_sub(C, C, myNum);
            mpz_divexact(C, C, A);

            std::vector<std::size_t> polySieveD(facSize * 2, 0u);

            for (std::size_t i = 0, row = 0; i < facSize; ++i, row += 2) {
                mpz_invert(Atemp2, A, mpzFacBase[i]);

                mpz_ui_sub(temp, SieveDist[row], B);
                mpz_mul(temp, temp, Atemp2);
                mpz_mod_ui(temp, temp, facBase[i]);
                polySieveD[row] = mpz_get_ui(temp);

                mpz_ui_sub(temp, SieveDist[row + 1], B);
                mpz_mul(temp, temp, Atemp2);
                mpz_mod_ui(temp, temp, facBase[i]);
                polySieveD[row + 1] = mpz_get_ui(temp);
            }

            for (std::size_t i = 0; i < LenB2; ++i) {
                mpz_mul_si(temp, B, myInterval[i]);
                mpz_mul_2exp(temp, temp, 1);
                mpz_add(temp, temp, C);
                mpz_set_d(Atemp2, myIntervalSqrd[i]);
                mpz_mul(Atemp2, Atemp2, A);
                mpz_add(sqrDiff[i], Atemp2, temp);
            }

            sieveLists(facSize, facBase, LenB2, sqrDiff.get(),
                       LnFB, myLogs, minPrime, polySieveD, lowBound);

            largeLogs.clear();

            for (std::size_t i = 0; i < LenB2; ++i)
                if (myLogs[i] > theCut)
                    largeLogs.push_back(i);

            largeLogsSize = largeLogs.size();
            myMat = std::vector<uint8_t>(largeLogsSize * colWidth, 0);

            smoothFacsRow.clear();

            if (largeLogsSize > 0) {
                for (std::size_t j = 0; j < largeLogsSize; ++j) {
                    std::vector<std::size_t> primeIndexVec;

                    // Add the index referring to A^2.. (i.e. add it twice)
                    primeIndexVec.push_back(mpzFacSize);
                    primeIndexVec.push_back(mpzFacSize);

                    // If Negative, we push zero (i.e. the index referring to -1)
                    if (mpz_sgn(sqrDiff[largeLogs[j]]) < 1) {
                        myMat[j * colWidth] = 1;
                        mpz_abs(sqrDiff[largeLogs[j]], sqrDiff[largeLogs[j]]);
                        primeIndexVec.push_back(0u);
                    }

                    for (std::size_t i = 0; i < facSize; ++i) {
                        while (mpz_divisible_ui_p(sqrDiff[largeLogs[j]], facBase[i])) {
                            mpz_divexact_ui(sqrDiff[largeLogs[j]],
                                            sqrDiff[largeLogs[j]], facBase[i]);
                            ++myMat[j * colWidth + i + 1];
                            primeIndexVec.push_back(i + 1);
                        }
                    }

                    mpz_mul_si(temp, A, myInterval[largeLogs[j]]);

                    if (mpz_cmp_ui(sqrDiff[largeLogs[j]], 1) == 0) {
                        
                        // Found a smooth number
                        smoothFacsRow.push_back(j * colWidth);
                        mpz_add(smoothInterval[nSmooth], temp, B);
                        
                        // auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(smoothInterval[nSmooth], 10) + 2);
                        // mpz_get_str(buffer.get(), 10, smoothInterval[nSmooth]);
                        // std::string strMyNum = buffer.get();
                        // 
                        // smoothFile << strMyNum << "\n";
                        
                        ++nSmooth;
                    } else if (mpz_cmp_d(sqrDiff[largeLogs[j]], partFacLimit) < 0) {
                        const uint64_t myKey = static_cast<uint64_t>(mpz_get_d(sqrDiff[largeLogs[j]]));
                        const auto pFacIt = partFactorsMap.find(myKey);
                        
                        // if (mpz_cmp(sqrDiff[largeLogs[j]], myMax) > 0) {
                        //     mpz_set(myMax, sqrDiff[largeLogs[j]]);
                        // }
                        // 
                        // if (mpz_cmp_ui(sqrDiff[largeLogs[j]], numBreak) > 0) {
                        //     ++count;
                        // }

                        if (pFacIt != partFactorsMap.end()) {
                            const auto trackIt = keepingTrack.find(myKey);
                            // double myVal = mpz_get_d(sqrDiff[largeLogs[j]]);
                            // valFile << myVal << "\n";
                            // if (mpz_cmp(sqrDiff[largeLogs[j]], realMax) > 0) {
                            //     mpz_set(realMax, sqrDiff[largeLogs[j]]);
                            // }
                            
                            if (trackIt != keepingTrack.end()) {
                                coFactorIndexVec.push_back(trackIt->second);
                            } else {
                                keepingTrack[myKey] = coFactorInd;
                                mpz_set(largeCoFactors[coFactorInd], sqrDiff[largeLogs[j]]);
                                coFactorIndexVec.push_back(coFactorInd);
                                ++coFactorInd;
                            }

                            for (const auto p: pFacIt->second)
                                primeIndexVec.push_back(p);

                            powsOfPartials.push_back(primeIndexVec);
                            const auto intervalIt = partIntvlMap.find(myKey);

                            mpz_add(temp, temp, B);
                            // auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(intervalIt->second, 10) + 2);
                            // mpz_get_str(buffer.get(), 10, intervalIt->second);
                            // std::string strMyNum = buffer.get();
                            // 
                            // partialFile << strMyNum << "\n";
                            // 
                            // buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(temp, 10) + 2);
                            // mpz_get_str(buffer.get(), 10, temp);
                            // strMyNum = buffer.get();
                            // 
                            // partialFile << strMyNum << "\n";
                            // 
                            mpz_mul(partialInterval[nPartial],
                                    temp, intervalIt->second);

                            partFactorsMap.erase(pFacIt);
                            partIntvlMap.erase(intervalIt);
                            ++nPartial;
                        } else {
                            partFactorsMap[myKey] = primeIndexVec;
                            mpz_init(partIntvlMap[myKey]);
                            mpz_add(partIntvlMap[myKey], temp, B);
                        }
                    }
                }
            }

            lenM = smoothFacsRow.size();
            
            if (lenM > 0) {
                tempMat = std::vector<uint8_t>(lenM * colWidth, 0u);
    
                for (std::size_t i = 0, rowInd = 0; i < lenM; ++i, rowInd += colWidth)
                    for (std::size_t j = 0; j <= facSize; ++j)
                        tempMat[rowInd + j] = myMat[smoothFacsRow[i] + j];
    
                powsOfSmooths.push_back(tempMat);
                powsFacInd.push_back(mpzFacSize);
            }
            
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
                updateStatTime(currLim, facSize,
                               checkPoint3 - checkPoint0, showStatsTime);
            }
        }
        
        // gmp_printf("myMax : %Zd\n", myMax);
        // gmp_printf("realMax : %Zd\n", realMax);
        // Rcpp::Rcout << count << ' ' << partFactorsMap.size() << std::endl;
        
        const std::size_t matRow = nSmooth + nPartial;
        const std::size_t matWidth = coFactorInd + mpzFacSize + 1;

        auto newTestInt = FromCpp14::make_unique<mpz_t[]>(matRow);
        std::vector<uint8_t> newMat(matRow * matWidth, 0);

        std::size_t k = 0;
        std::size_t r = 0;
        std::size_t row = 0;

        for (const auto &pows: powsOfSmooths) {
            for (std::size_t i = 0; i < pows.size(); i += colWidth, row += matWidth, ++r) {
                mpz_init_set(newTestInt[r], smoothInterval[r]);

                for (std::size_t j = 0; j < colWidth; ++j)
                    newMat[row + j] = pows[i + j];

                newMat[row + powsFacInd[k]] = (k > 0) ? 2u : 0u;
            }

            ++k;
        }

        if ((coFactorInd + mpzFacSize) >= mpzFSzLIMIT) {
            mpzFSzLIMIT <<= 1;
            Grow(mpzFacBase, mpzFacSize, mpzFSzLIMIT);
        }

        for (std::size_t i = 0, j = mpzFacSize; i < coFactorInd; ++i, ++j)
            mpz_set(mpzFacBase[j], largeCoFactors[i]);

        for (std::size_t i = 0, r = nSmooth,
             row = matWidth * nSmooth; i < nPartial; ++i, ++r, row += matWidth) {

            mpz_init_set(newTestInt[r], partialInterval[i]);

            for (const auto p: powsOfPartials[i])
                ++newMat[row + p];

            newMat[row + mpzFacSize + 1 + coFactorIndexVec[i]] = 2u;
        }
        
        // partialFile.close();
        // smoothFile.close();
        // valFile.close();
        
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
    
    for (std::size_t i = 0; i < LenB2; ++i)
        mpz_clear(sqrDiff[i]);
    
    sqrDiff.reset();
    
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
    coFactorIndexVec.clear();
    
    for (auto v: partIntvlMap)
        mpz_clear(v.second);
    
    partIntvlMap.clear();
    partFactorsMap.clear();
    keepingTrack.clear();
    
    mpz_clear(temp); mpz_clear(sqrtInt); mpz_clear(A);
    mpz_clear(C); mpz_clear(Atemp); mpz_clear(Btemp);
    mpz_clear(lowBound); mpz_clear(B); mpz_clear(Atemp2);

    for (std::size_t i = 0; i < 10; ++i)
        mpz_clear(TS[i]);
}
