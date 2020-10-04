#include "Polynomial.h"
#include <Rcpp.h>

constexpr std::size_t ParBitCutOff = 160u;

// There are 2 variables below that depend on the size of myNum in a
// non-Trivial way. They are dblMyTarget and dblLenB.
// Below, we show how to reproduce these numbers. They were all
// obtained through experimentation and may need to be readjusted in
// the future. In particular, they were targeted to the primary dev
// machine (i.e. 2017 Mac i7 16GB). There could be a more general
// method for obtaining these values that better suits any platform.
// 
// These values were modified from what is suggested in:
// "The Multiple Polynomial Quadratic Sieve" by Robert D. Silverman
//
// DigSize <- c(24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84)
// FBSize <- c(133, 266, 532, 1197, 1596, 2660, 3990, 5985, 7980, 10640, 13965)
// MSize <- c(7, 28, 49, 140, 175, 350, 490, 700, 925, 1150, 1400)
//
// dblMyTarget:
// rawCoef <- round(unname(lm(FBSize ~ poly(DigSize, 4,
//                                      raw = TRUE))$coefficients), 7)
// names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
// rawCoef
//     intercept           x^1           x^2           x^3           x^4
//  -1196.999999    25.9297462  - 4.7601496     0.0794196  -0.000179411
//
//
// dblLenB:
// rawCoef <- round(unname(lm(MSize ~ poly(DigSize, 4,
//                                      raw = TRUE))$coefficients), 7)
// names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
// rawCoef
//    intercept           x^1           x^2           x^3           x^4
// -422.2004662    44.2665760    -1.7020526     0.0281537    -0.0001304
//
// Note that the smallest that dblDigCount can be in order for dblLenB to be
// valid is 22. This is no problem as we factor numbers less than 23 digits
// with Pollard's Rho algorithm.

// lm(MSize[10:11] ~ DigSize[10:11])
// 
// Call:
//     lm(formula = MSize[10:11] ~ DigSize[10:11])
//     
//     Coefficients:
//     (Intercept)  DigSize[10:11]  
//        -2300.00              45
//
// ***********************************************************************

inline double GetIntervalSize(double dig) {
    if (dig < 85) {
        return std::ceil(34.0943408 * dig - 1.4227370
                             * std::pow(dig, 2.0) + 0.0253568
                             * std::pow(dig, 3.0) - 0.0001237
                             * std::pow(dig, 4.0) - 300.8135198);
    } else {
        return std::ceil(41.67 * dig - 2100.0);
    }
}

void QuadraticSieve(const mpz_class &myNum, std::vector<mpz_class> &factors,
                    std::size_t nThreads, bool bShowStats) {
    
    const auto checkPoint0 = std::chrono::steady_clock::now();
    const std::size_t digCount = mpz_sizeinbase(myNum.get_mpz_t(), 10);
    const std::size_t bits = mpz_sizeinbase(myNum.get_mpz_t(), 2);
    const bool IsParallel = (nThreads == 1 || bits < ParBitCutOff) ? false : true;
    
    const double lognum = bits / std::log2(std::exp(1.0));
    const double sqrLogLog = std::sqrt(lognum * std::log(lognum));
    const double dblDigCount = digCount;
    
    double fudge1 = -0.4;
    double LimB = std::exp((0.5 + fudge1) * sqrLogLog);
    const double dblMyTarget = std::ceil(125.9297462 * dblDigCount - 4.7601496
                                             * std::pow(dblDigCount, 2.0) + 0.0794196
                                             * std::pow(dblDigCount, 3.0) - 0.0001794
                                             * std::pow(dblDigCount, 4.0) - 1196.9999);
    
    std::size_t myTarget = static_cast<std::size_t>(dblMyTarget);
    
    while (LimB < myTarget) {
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        fudge1 += 0.001;
    }
    
    const std::vector<int> facBase = GetPrimesQuadRes(myNum, LimB,
                                                      fudge1, sqrLogLog, myTarget);
    
    std::vector<int> primesAndPows;
    std::vector<int> LnFB;
        
    const double dblLenB = GetIntervalSize(dblDigCount);
    const std::size_t LenB = static_cast<std::size_t>(dblLenB) * 1000;
    const int DoubleLenB = 2 * LenB + 1;
    
    mpz_class Temp = sqrt(myNum);
    Temp -= static_cast<unsigned long int>(LenB);
    const int minPrime = static_cast<int>(mpz_sizeinbase(Temp.get_mpz_t(), 10) * 2);
    
    const std::size_t facSize = facBase.size();
    const int vecMaxSize = std::min(static_cast<int>(((facBase.back()
                                     + L1Cache - 1) / L1Cache) * L1Cache), DoubleLenB);
    
    const std::vector<std::size_t> SieveDist = SetSieveDist(facBase, primesAndPows, LnFB,myNum,
                                                            static_cast<std::size_t>(DoubleLenB),
                                                            minPrime);
    
    // This array will be passed to solutionSeach.
    std::vector<mpz_class> mpzFacBase;
    
    for (const auto fb: facBase)
        mpzFacBase.push_back(fb);
    
    Temp = myNum * 2;
    Temp = sqrt(Temp);
    Temp *= static_cast<unsigned long int>(LenB);
    
    const double fudge2 = (digCount < 45) ? 1.450 :
                          (digCount < 50) ? 1.505 :
                          (digCount < 60) ? 1.560 : 
                          (digCount < 65) ? 1.570 : 1.575;
    
    const int theCut = std::ceil(100.0 * fudge2 *
                                 static_cast<double>(mpz_sizeinbase(Temp.get_mpz_t(), 10)));
    
    mpz_class LowBound;
    LowBound = -1 * static_cast<int>(LenB);
    
    mpz_class NextPrime;
    mpz_mul_2exp(NextPrime.get_mpz_t(), myNum.get_mpz_t(), 1);
    NextPrime = sqrt(NextPrime);
    NextPrime /= static_cast<unsigned long int>(LenB);
    NextPrime = sqrt(NextPrime);
    
    if (cmp(NextPrime, facBase.back()) < 0)
        NextPrime = facBase.back();
    
    Polynomial myPoly(facSize, SieveDist.size(), bShowStats, myNum);
    bool xtraTime = true;
    
    if (IsParallel) {
        myPoly.FactorParallel(SieveDist, primesAndPows, facBase, LnFB, mpzFacBase, NextPrime,
                              LowBound, myNum, theCut, DoubleLenB, vecMaxSize,
                              checkPoint0, nThreads);

        if (myPoly.ContinueToSolution()) {
            myPoly.GetSolution(mpzFacBase, facBase, factors,
                               myNum, nThreads, checkPoint0);
            myPoly.MakeStatsFalse();

            if (bShowStats && cmp(factors[0], 0) == 0) {
                RcppThread::Rcout << "|      Extra Time      |\n|--------------------|" << std::endl;
                xtraTime = false;
            }
        }

        NextPrime = mpzFacBase.back();
    }

    auto t0 = std::chrono::steady_clock::now();
    bool bUpdateXtra = false;

    while (cmp(factors[0], 0) == 0) {
        myPoly.FactorSerial(SieveDist, primesAndPows, facBase, LnFB, mpzFacBase, NextPrime,
                            LowBound, myNum, theCut, DoubleLenB, vecMaxSize,
                            checkPoint0);
        // factors[0] = 11;
        // factors[1] = 13;
        myPoly.GetSolution(mpzFacBase, facBase, factors,
                           myNum, nThreads, checkPoint0);
        myPoly.MakeStatsFalse();
        NextPrime = mpzFacBase.back();

        if (bShowStats && cmp(factors[0], 0) == 0 && xtraTime) {
            t0 = std::chrono::steady_clock::now();
            RcppThread::Rcout << " Unsuccessful initial factorization... more smooths needed \n" << std::endl;
            RcppThread::Rcout << "|     Extra Time     |\n|--------------------|" << std::endl;
            xtraTime = false;
            bUpdateXtra = true;
        }

        if (bShowStats && bUpdateXtra) {
            OneColumnStats(std::chrono::steady_clock::now() - t0);

            if (cmp(factors[0], 0) > 0) {
                RcppThread::Rcout << "\n" << std::endl;
            }
        }
    }
}
