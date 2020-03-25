#include "SolutionSearch.h"
#include "Cpp14MakeUnique.h"
#include "TonelliShanks.h"
#include <RcppThread.h>
#include <unordered_map>

void QuadraticSieve(mpz_t myNum, double fudge1, double fudge2,
                    int64_t LenB, mpz_t *const factors) {
    
    std::size_t digCount = mpz_sizeinbase(myNum, 10);
    std::size_t bits = mpz_sizeinbase(myNum, 2);
    std::size_t myTarget;
    
    double lognum = bits / log2(std::exp(1.0));
    double sqrLogLog = std::sqrt(lognum * std::log(lognum));
    
    mpz_t currP, nextP, resTest, CP1;
    mpz_init(currP);
    mpz_init(nextP);
    mpz_init(CP1);
    mpz_init(resTest);
    v1d facBase;
    
    double dblDigCount = digCount;
    double dblMyTarget;
    double LimB;
    
    // These values were obtained from "The Multiple Polynomial
    // Quadratic Sieve" by Robert D. Silverman
    // DigSize <- c(24,30,36,42,48,54,60,66)
    // FBSize <- c(100,200,400,900,1200,2000,3000,4500)
    // MSize <- c(5,20,35,100,125,250,350,500)
    //
    // rawCoef <- round(unname(lm(FBSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept          x^1          x^2          x^3          x^4
    // 3637.0670996 -391.8275012   15.1541456   -0.2475566    0.0016806
    
    if (fudge1 == 0) {
        fudge1 = -0.4;
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        
        dblMyTarget = -391.8275012 * dblDigCount + 15.1541456 * std::pow(dblDigCount, 2.0);
        dblMyTarget += -0.2475566 * std::pow(dblDigCount, 3.0) + 0.0016806 * std::pow(dblDigCount, 4.0) + 3637.0671;
        dblMyTarget = std::ceil(dblMyTarget);
        myTarget = static_cast<std::size_t>(dblMyTarget);
        
        while (LimB < myTarget) {
            LimB = std::exp((0.5 + fudge1) * sqrLogLog);
            fudge1 += 0.001;
        }
        
        facBase = getPrimesQuadRes(myNum, LimB);
        
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
                    facBase.push_back(mpz_get_si(currP));
            }
        }
    } else {
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        facBase = getPrimesQuadRes(myNum, LimB);
    }
    
    // rawCoef <- round(unname(lm(MSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept           x^1           x^2           x^3           x^4
    // -213.1466450    23.3947361    -0.9494686     0.0165574    -0.0000767
    
    if (LenB == 0) {
        double dblLenB = 23.394736 * dblDigCount - 0.9494686 * std::pow(dblDigCount, 2.0);
        dblLenB += 0.0165574 * std::pow(dblDigCount, 3.0) - 0.0000767 * std::pow(dblDigCount, 4.0) - 213.1466450;
        dblLenB = std::ceil(dblLenB);
        LenB = static_cast<int64_t>(dblLenB);
        LenB *= 1000;
    }

    const std::size_t facSize = facBase.size();
    std::size_t facSize2 = facSize;

    // facBase2 will be used if multiple polynomials are needed.
    // With every iteration, a prime will be added to the factor
    // base for additional factorization. The original facBase
    // should not be tampered with, hence the need for facBase2.
    v1d facBase2 = facBase;

    mpz_t sqrtInt;
    mpz_init(sqrtInt);
    mpz_sqrt(sqrtInt, myNum);
    v1d myInterval;

    int64_t Lower = -1 * LenB;
    int64_t Upper = LenB;
    int64_t LenB2 = 2 * LenB + 1;

    uint64_t uLenB2 = LenB2;
    myInterval.reserve(LenB2);

    for (int64_t i = Lower; i <= Upper; ++i)
        myInterval.push_back(i);

    v2d SieveDist(facSize, v1d(2));
    SieveDist[0][0] = SieveDist[0][1] = 1;
    setSieveDist(myNum, facBase, facSize, SieveDist);

    // testInterval will be the values that are associated
    // with smooth numbers. We make it 10x the size of facSize
    // in order to guarantee we will have enough space
    auto testInterval = FromCpp14::make_unique<mpz_t[]>(facSize * 10);

    // Similarly, we create an array that will house the
    // interval entries corresponding to the partially smooth facs
    auto partialInterval = FromCpp14::make_unique<mpz_t[]>(facSize * 10);

    // This array will be passed to solutionSeach. It contains
    // facBase as well as the Atemps and large common cofactors
    auto mpzFacBase = FromCpp14::make_unique<mpz_t[]>(facSize * 10);

    // This array will store the large cofactors that will
    // later be added to mpzFacBase
    auto largeCoFactors = FromCpp14::make_unique<mpz_t[]>(facSize * 10);

    for (std::size_t i = 0; i < (facSize * 10); ++i) {
        mpz_init(testInterval[i]);
        mpz_init(partialInterval[i]);
        mpz_init(mpzFacBase[i]);
        mpz_init(largeCoFactors[i]);
    }

    auto largeInterval = FromCpp14::make_unique<mpz_t[]>(LenB2);
    auto sqrDiff = FromCpp14::make_unique<mpz_t[]>(LenB2);

    for (std::size_t i = 0; i < uLenB2; ++i) {
        mpz_init(largeInterval[i]);
        mpz_init(sqrDiff[i]);
    }

    mpz_sub_ui(largeInterval[0], sqrtInt, Upper);

    for (std::size_t i = 1; i < uLenB2; ++i)
        mpz_add_ui(largeInterval[i], largeInterval[i - 1], 1u);

    mpz_t temp;
    mpz_init(temp);

    for (std::size_t i = 0; i < uLenB2; ++i) {
        mpz_pow_ui(temp, largeInterval[i], 2u);
        mpz_sub(sqrDiff[i], temp, myNum);
    }

    v3d FBDivSieve = sieveLists(facSize, facBase, LenB2, sqrDiff.get());

    mpz_mul_ui(temp, myNum, 2);
    mpz_sqrt(temp, temp);
    mpz_mul_ui(temp, temp, Upper);

    if (fudge2 == 0) {
        if (digCount < 30)
            fudge2 = 0.7;
        else if (digCount < 35)
            fudge2 = 0.9;
        else if (digCount < 40)
            fudge2 = 1.05;
        else if (digCount < 45)
            fudge2 = 1.1;
        else if (digCount < 50)
            fudge2 = 1.2;
        else
            fudge2 = 1.5;
    }
    
    const double theCut = fudge2 * mpz_sizeinbase(temp, 10);
    std::vector<double> LnFB;
    std::vector<double> myLogs(LenB2, 0);
    LnFB.reserve(facSize);

    for (std::size_t i = 0; i < facSize; ++i) {
        mpz_set_ui(mpzFacBase[i], facBase[i]);
        LnFB.push_back(std::log(static_cast<double>(facBase[i])));
    }

    mpz_mul_ui(temp, sqrtInt, Upper);
    int64_t minPrime = static_cast<int64_t>(mpz_sizeinbase(temp, 10) * 2);

    std::vector<std::vector<std::size_t>> indexDiv(LenB2,
                                                   std::vector<std::size_t>());

    std::size_t strt = (mpz_even_p(sqrDiff[0]) != 0) ? 0 : 1;

    for (std::size_t j = strt; j < uLenB2; j += 2)
        indexDiv[j].push_back(0);

    for (std::size_t i = 1; i < facSize; ++i) {
        for (std::size_t k = 0; k <= 1; ++k)
            for (std::size_t j = 0; j < FBDivSieve[i][k].size(); ++j)
                indexDiv[FBDivSieve[i][k][j]].push_back(i);

        if (facBase[i] > minPrime)
            for (std::size_t k = 0; k <= 1; ++k)
                for (std::size_t j = 0; j < FBDivSieve[i][k].size(); ++j)
                    myLogs[FBDivSieve[i][k][j]] += LnFB[i];
    }

    std::vector<std::size_t> largeLogs;

    for (std::size_t i = 0; i < uLenB2; i++)
        if (myLogs[i] > theCut)
            largeLogs.push_back(i);

    std::size_t largeLogsSize = largeLogs.size();
    const std::size_t colWidth = facSize + 1;
    std::vector<uint8_t> myMat(largeLogsSize * colWidth, 0u);

    for (std::size_t i = 0; i < largeLogsSize; ++i) {
        if (mpz_sgn(sqrDiff[largeLogs[i]]) < 0) {
            myMat[i * colWidth] = 1;
            mpz_abs(sqrDiff[largeLogs[i]], sqrDiff[largeLogs[i]]);
        }
    }

    mpz_t rem, quot;
    mpz_init(rem);
    mpz_init(quot);
    std::vector<std::size_t> sFacs;
    std::size_t numSmooth = 0;

    // The map below is for almost smooth numbers, i.e. numbers
    // that are almost completely factored by our factor base
    std::unordered_map<uint64_t, std::vector<std::size_t>> partialFactors;
    std::unordered_map<uint64_t, mpz_t> intervalPartial;
    
    std::size_t partialCount = 0;
    std::vector<std::vector<std::size_t>> powersOfPartials;
    std::unordered_map<uint64_t, std::size_t> keepingTrack;
    std::vector<std::size_t> coFactorIndexVec;
    
    std::size_t coFactorInd = 0u;

    if (largeLogsSize > 0) {
        for (std::size_t j = 0; j < largeLogsSize; ++j) {
            std::vector<std::size_t> primeIndexVec;
            
            if (myMat[j * colWidth])
                primeIndexVec.push_back(0u);

            for (std::size_t i = 0; i < indexDiv[largeLogs[j]].size(); ++i) {
                bool divides = true;
                const int64_t primeFactor = facBase[indexDiv[largeLogs[j]][i]];
                
                while (divides) {
                    mpz_fdiv_qr_ui(quot, rem, sqrDiff[largeLogs[j]], primeFactor);
                    divides = (mpz_cmp_ui(rem, 0) == 0);
                    
                    if (divides) {
                        mpz_set(sqrDiff[largeLogs[j]], quot);
                        ++myMat[j * colWidth + indexDiv[largeLogs[j]][i] + 1];
                        primeIndexVec.push_back(indexDiv[largeLogs[j]][i] + 1);
                    }
                }
            }

            if (mpz_cmp_ui(sqrDiff[largeLogs[j]], 1) == 0) {
                // Found a smooth number
                mpz_set(testInterval[numSmooth], largeInterval[largeLogs[j]]);
                sFacs.push_back(j * colWidth);
                ++numSmooth;
            } else {
                const uint64_t myKey = makeKey(sqrDiff[largeLogs[j]]);
                const auto pFacIt = partialFactors.find(myKey);

                if (pFacIt != partialFactors.end()) {
                    keepingTrack[myKey] = 0u;
                    mpz_set(largeCoFactors[partialCount], sqrDiff[largeLogs[j]]);
                    
                    Rprintf("DIFFIDFd\n");
                    
                    for (const auto p: pFacIt->second)  
                        primeIndexVec.push_back(p);
                    
                    powersOfPartials.push_back(primeIndexVec);
                    const auto intervalIt = intervalPartial.find(myKey);
                    
                    mpz_mul(partialInterval[partialCount], 
                            largeInterval[largeLogs[j]],
                            intervalIt->second);
                    
                    partialFactors.erase(pFacIt);
                    intervalPartial.erase(intervalIt);
                    ++partialCount;
                    
                    // const auto trackIt = keepingTrack.find(myKey);
                    // 
                    // if (trackIt != keepingTrack.end()) {
                    //     coFactorIndexVec.push_back(trackIt->second);
                    // } else {
                    //     keepingTrack[myKey] = coFactorInd;
                    //     mpz_set(largeCoFactors[coFactorInd], sqrDiff[largeLogs[j]]);
                    //     coFactorIndexVec.push_back(coFactorInd);
                    //     ++coFactorInd;
                    // }
                    // 
                    // for (const auto p: pFacIt->second)  
                    //     primeIndexVec.push_back(p);
                    // 
                    // powersOfPartials.push_back(primeIndexVec);
                    // const auto intervalIt = intervalPartial.find(myKey);
                    // 
                    // mpz_add(temp, temp, B2);
                    // mpz_mul(partialInterval[partialCount], 
                    //         temp, intervalIt->second);
                    // 
                    // partialFactors.erase(pFacIt);
                    // intervalPartial.erase(intervalIt);
                    // ++partialCount;
                } else {
                    partialFactors[myKey] = primeIndexVec;
                    mpz_set(intervalPartial[myKey], largeInterval[largeLogs[j]]);
                }
            }
        }
    }

    std::size_t lenM = sFacs.size();
    std::vector<std::vector<uint8_t>> powsOfSmooths;
    std::vector<uint8_t> tempMat(lenM * colWidth, 0u);

    for (std::size_t i = 0, rowInd = 0; i < lenM; ++i, rowInd += colWidth)
        for (std::size_t j = 0; j <= facSize; ++j)
            tempMat[rowInd + j] = myMat[sFacs[i] + j];
    
    powsOfSmooths.push_back(tempMat);

    mpz_t A, B, B2, C, Atemp, Btemp, Atemp2, Aalt, AaltTemp, B2alt;
    mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(Atemp2);
    mpz_init(Atemp); mpz_init(Btemp); mpz_init(B2);
    mpz_init(Aalt); mpz_init(B2alt); mpz_init(AaltTemp);

    mpz_mul_2exp(Atemp, myNum, 1);
    mpz_sqrt(Atemp, Atemp);
    mpz_div_ui(Atemp, Atemp, Upper);
    mpz_sqrt(Atemp, Atemp);

    int64_t maxFBase = *std::max_element(facBase.begin(), facBase.end());

    if (mpz_cmp_ui(Atemp, maxFBase) < 0)
        mpz_set_ui(Atemp, maxFBase);

    mpz_t quadRes[2];
    mpz_init(quadRes[0]);
    mpz_init(quadRes[1]);
    
    std::size_t numPolys = 0;
    std::size_t extraFacs = 0;
    
    v1d myAtemps, myQuadRes;
    std::vector<double> myIntervalSqrd(LenB2);

    if (facSize2 > lenM)
        for (std::size_t i = 0; i < uLenB2; i++)
            myIntervalSqrd[i] = static_cast<double>(myInterval[i] * myInterval[i]);
    
    while (mpz_cmp_ui(factors[0], 0) == 0) {
        // Find enough smooth numbers to guarantee a non-trivial solution
        while ((numSmooth + partialCount) <= facSize + extraFacs) {
            bool LegendreTest = true;

            while (LegendreTest) {
                mpz_nextprime(Atemp, Atemp);
                
                if (mpz_legendre(myNum, Atemp) == 1)
                    LegendreTest = false;
            }

            myAtemps.push_back(static_cast<int64_t>(mpz_get_d(Atemp)));
            facBase2.push_back(static_cast<int64_t>(mpz_get_d(Atemp)));
            ++facSize2;

            mpz_pow_ui(A, Atemp, 2);
            TonelliShanksC(myNum, Atemp, quadRes);

            if (mpz_cmp(quadRes[0], quadRes[1]) > 0) {
                mpz_set(Btemp, quadRes[0]);
            } else {
                mpz_set(Btemp, quadRes[1]);
            }

            myQuadRes.push_back(static_cast<int64_t>(mpz_get_d(Btemp)));

            mpz_mul_2exp(temp, Btemp, 1);
            mpz_invert(temp, temp, Atemp);
            mpz_pow_ui(B2, Btemp, 2);
            mpz_sub(B2, myNum, B2);
            mpz_mul(B2, B2, temp);
            mpz_add(B2, B2, Btemp);
            mpz_mod(B2, B2, A);

            mpz_pow_ui(C, B2, 2);
            mpz_sub(C, C, myNum);
            mpz_divexact(C, C, A);
            v2d polySieveD = v2d(facSize, v1d(2));

            for (std::size_t i = 0; i < facSize; ++i) {
                mpz_invert(Atemp2, A, mpzFacBase[i]);

                for (std::size_t j = 0; j <= 1; ++j) {
                    mpz_ui_sub(temp, SieveDist[i][j], B2);
                    mpz_mul(temp, temp, Atemp2);
                    mpz_mod_ui(temp, temp, facBase[i]);
                    polySieveD[i][j] = static_cast<int64_t>(mpz_get_si(temp));
                }
            }

            for (std::size_t i = 0; i < uLenB2; ++i) {
                mpz_mul_si(temp, B2, myInterval[i]);
                mpz_mul_2exp(temp, temp, 1);
                mpz_add(temp, temp, C);
                mpz_set_d(Atemp2, myIntervalSqrd[i]);
                mpz_mul(Atemp2, Atemp2, A);
                mpz_add(sqrDiff[i], Atemp2, temp);
            }

            FBDivSieve = sieveLists(facSize, facBase, LenB2, sqrDiff.get());
            std::size_t strt = (mpz_even_p(sqrDiff[0]) != 0) ? 0 : 1;
            indexDiv = std::vector<std::vector<std::size_t>>(LenB2, std::vector<std::size_t>());

            for (std::size_t j = strt; j < uLenB2; j += 2)
                indexDiv[j].push_back(0);

            std::fill(myLogs.begin(), myLogs.end(), 0.0);

            for (std::size_t i = 1; i < facSize; ++i) {
                for (std::size_t k = 0; k <= 1; ++k)
                    for (std::size_t j = 0; j < FBDivSieve[i][k].size(); ++j)
                        indexDiv[FBDivSieve[i][k][j]].push_back(i);

                if (facBase[i] > minPrime)
                    for (std::size_t k = 0; k <= 1; ++k)
                        for (std::size_t j = 0; j < FBDivSieve[i][k].size(); ++j)
                            myLogs[FBDivSieve[i][k][j]] += LnFB[i];
            }

            largeLogs.clear();

            for (int64_t i = 0; i < LenB2; i++)
                if (myLogs[i] > theCut)
                    largeLogs.push_back(i);

            largeLogsSize = largeLogs.size();
            myMat = std::vector<uint8_t>(largeLogsSize * colWidth, 0);

            for (std::size_t i = 0; i < largeLogsSize; ++i) {
                if (mpz_sgn(sqrDiff[largeLogs[i]]) < 1) {
                    myMat[i * colWidth] = 1;
                    mpz_abs(sqrDiff[largeLogs[i]], sqrDiff[largeLogs[i]]);
                }
            }

            sFacs.clear();

            if (largeLogsSize > 0) {
                for (std::size_t j = 0; j < largeLogsSize; ++j) {
                    std::vector<std::size_t> primeIndexVec;
                    
                    // Add the index referring to A^2.. (i.e. add it twice)
                    primeIndexVec.push_back(facSize2);
                    primeIndexVec.push_back(facSize2);
                    
                    // If Negative, we push zero (i.e. the index referring to -1)
                    if (myMat[j * colWidth])
                        primeIndexVec.push_back(0u);

                    for (std::size_t i = 0; i < indexDiv[largeLogs[j]].size(); ++i) {
                        bool divides = true;
                        const int64_t primeFactor = facBase[indexDiv[largeLogs[j]][i]];

                        while (divides) {
                            mpz_fdiv_qr_ui(quot, rem, sqrDiff[largeLogs[j]], primeFactor);
                            divides = (mpz_cmp_ui(rem, 0) == 0);

                            if (divides) {
                                mpz_set(sqrDiff[largeLogs[j]], quot);
                                ++myMat[j * colWidth + indexDiv[largeLogs[j]][i] + 1];
                                primeIndexVec.push_back(indexDiv[largeLogs[j]][i] + 1);
                            }
                        }
                    }
                    
                    mpz_mul_si(temp, A, myInterval[largeLogs[j]]);

                    if (mpz_cmp_ui(sqrDiff[largeLogs[j]], 1) == 0) {
                        // Found a smooth number
                        sFacs.push_back(j * colWidth);
                        mpz_add(testInterval[numSmooth], temp, B2);
                        ++numSmooth;
                    } else {
                        const uint64_t myKey = makeKey(sqrDiff[largeLogs[j]]);
                        const auto pFacIt = partialFactors.find(myKey);
                        
                        if (pFacIt != partialFactors.end()) {
                            const auto trackIt = keepingTrack.find(myKey);
                            
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
                            
                            powersOfPartials.push_back(primeIndexVec);
                            const auto intervalIt = intervalPartial.find(myKey);
                            
                            mpz_add(temp, temp, B2);
                            mpz_mul(partialInterval[partialCount], 
                                    temp, intervalIt->second);
                            
                            partialFactors.erase(pFacIt);
                            intervalPartial.erase(intervalIt);
                            ++partialCount;
                        } else {
                            partialFactors[myKey] = primeIndexVec;
                            mpz_add(intervalPartial[myKey], temp, B2);
                        }
                    }
                }
            }
            
            lenM = sFacs.size();
            tempMat = std::vector<uint8_t>(lenM * colWidth, 0u);

            for (std::size_t i = 0, rowInd = 0; i < lenM; ++i, rowInd += colWidth)
                for (std::size_t j = 0; j <= facSize; ++j)
                    tempMat[rowInd + j] = myMat[sFacs[i] + j];
            
            powsOfSmooths.push_back(tempMat);
            ++numPolys;
        }
        
        const std::size_t matRow = numSmooth + partialCount;
        const std::size_t matWidth = coFactorInd + facSize2 + 1;
        
        auto newTestInt = FromCpp14::make_unique<mpz_t[]>(matRow);
        std::vector<uint8_t> newMat(matRow * matWidth, 0);

        std::size_t fSize = facSize;
        std::size_t facBaseIndex = facSize;

        // This loop is one off from the loop below. Also note
        // the limit is k < numPolys NOT "<=".
        for (std::size_t k = 0; k < numPolys; ++k, ++facBaseIndex)
            mpz_set_si(mpzFacBase[facBaseIndex], facBase2[facBaseIndex]);

        for (std::size_t k = 0, r = 0, row = 0; k <= numPolys; ++k, ++fSize) {
            for (std::size_t i = 0; i < powsOfSmooths[k].size(); i += colWidth, row += matWidth, ++r) {
                mpz_init_set(newTestInt[r], testInterval[r]);

                for (std::size_t j = 0; j < colWidth; j++)
                    newMat[row + j] = powsOfSmooths[k][i + j];
                
                newMat[row + fSize] = (k) ? 2u : 0u;
            }
        }
        
        for (std::size_t i = 0; i < coFactorInd; ++i, ++facBaseIndex)
            mpz_set(mpzFacBase[facBaseIndex], largeCoFactors[i]);
        
        for (std::size_t i = 0, r = numSmooth,
             row = matWidth * numSmooth; i < partialCount; ++i, ++r, row += matWidth) {
            
            mpz_init_set(newTestInt[r], partialInterval[i]);
            
            for (const auto p: powersOfPartials[i])
                ++newMat[row + p];
            
            newMat[row + fSize + coFactorIndexVec[i]] = 2u;
        }
        
        solutionSearch(newMat, matRow, matWidth, myNum,
                       mpzFacBase.get(), newTestInt.get(), factors);
        
        extraFacs += 5;
    }
    
    for (std::size_t i = 0; i < uLenB2; i++) {
        mpz_clear(largeInterval[i]);
        mpz_clear(sqrDiff[i]);
    }

    mpz_clear(rem); mpz_clear(quot); mpz_clear(temp); mpz_clear(sqrtInt);
    mpz_clear(A); mpz_clear(B); mpz_clear(C); mpz_clear(Atemp);
    mpz_clear(quadRes[0]); mpz_clear(quadRes[1]);
}
