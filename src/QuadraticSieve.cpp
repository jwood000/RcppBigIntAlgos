#include "Polynomial.h"

constexpr std::size_t ParBitCutOff = 160u;

void QuadraticSieve(const mpz_class &myNum, std::vector<mpz_class> &factors,
                    std::size_t nThreads, bool bShowStats) {
    
    const auto checkPoint0 = std::chrono::steady_clock::now();
    const std::size_t digCount = mpz_sizeinbase(myNum.get_mpz_t(), 10);
    const std::size_t bits = mpz_sizeinbase(myNum.get_mpz_t(), 2);
    const bool IsParallel = (nThreads == 1 || bits < ParBitCutOff) ? false : true;
    
    const double lognum = bits / std::log2(std::exp(1.0));
    const double sqrLogLog = std::sqrt(lognum * std::log(lognum));
    const double dblDigCount = digCount;
    
    // These values were obtained from "The Multiple Polynomial
    // Quadratic Sieve" by Robert D. Silverman
    // DigSize <- c(24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84)
    // FBSize <- c(100, 200, 400, 900, 1200, 2000, 3000, 4500)
    // MSize <- c(5,20,35,100,125,250,350,500)
    // CounterSize <- c(10, 7, 5, 3.5, 2.25, 1.25, 1, 0.9)
    //
    // rawCoef <- round(unname(lm(FBSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //     intercept           x^1           x^2           x^3           x^4
    // -1197.0000000   125.9297462    -4.7601496     0.0794196    -0.0001794
    
    double fudge1 = -0.4;
    double LimB = std::exp((0.5 + fudge1) * sqrLogLog);
    
    const double dblMyTarget = std::ceil(125.9297462 * dblDigCount - 4.7601496
                                             * std::pow(dblDigCount, 2.0) + 0.0794196
                                             * std::pow(dblDigCount, 3.0) - 0.0001794
                                             * std::pow(dblDigCount, 4.0) - 1197);
    
    std::size_t myTarget = static_cast<std::size_t>(dblMyTarget);
    
    while (LimB < myTarget) {
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        fudge1 += 0.001;
    }
    
    const std::vector<int> facBase = GetPrimesQuadRes(myNum, LimB, fudge1, sqrLogLog, myTarget);
    
    // rawCoef <- round(unname(lm(MSize ~ poly(DigSize, 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //    intercept           x^1           x^2           x^3           x^4
    // -237.5850816    26.7743460    -1.1188164     0.0200060    -0.0000899
    //
    // Note that the smallest that dblDigCount can be in order for dblLenB to be
    // valid is 22. This is no problem as we factor numbers less than 23 digits
    // with Pollard's Rho algorithm.
    
    const double dblLenB = std::ceil(26.7743460 * dblDigCount - 1.1188164
                                         * std::pow(dblDigCount, 2.0) + 0.0200060
                                         * std::pow(dblDigCount, 3.0) - 0.0000899
                                         * std::pow(dblDigCount, 4.0) - 237.5850816);
    
    const std::size_t LenB = static_cast<std::size_t>(dblLenB) * 1000;
    const int DoubleLenB = 2 * LenB + 1;
    
    const std::size_t facSize = facBase.size();
    const int vecMaxSize = std::min(static_cast<int>((1 + facBase.back() / L1Cache) * L1Cache), DoubleLenB);
    const std::vector<std::size_t> SieveDist = SetSieveDist(facBase, myNum);
    
    // CounterSize <- c(10, 7, 5, 3.5, 2.25, 1.25, 1, 0.9)
    // rawCoef <- round(unname(lm(CounterSize ~ poly(DigSize[seq_along(CounterSize)], 4, raw = TRUE))$coefficients), 7)
    // names(rawCoef) <- c("intercept", "x^1", "x^2", "x^3", "x^4")
    // rawCoef
    //  intercept        x^1        x^2        x^3        x^4 
    // 40.1734307 -2.2578629  0.0571943 -0.0007403  0.0000039
    
    double dblFacMultiplier = std::max(std::ceil(-2.2578629 * dblDigCount + 0.0571943
                                                     * std::pow(dblDigCount, 2.0) - 0.0007403
                                                     * std::pow(dblDigCount, 3.0) + 0.0000039
                                                     * std::pow(dblDigCount, 4.0) + 40.1734307), 1.0);
    
    const std::size_t mpzContainerSize = static_cast<double>(facSize) * dblFacMultiplier;
    std::size_t mpzFSzLIMIT = facSize * 10;
    
    // This array will be passed to solutionSeach.
    std::vector<mpz_class> mpzFacBase;
    mpzFacBase.reserve(mpzFSzLIMIT);
    
    mpz_class Temp;
    Temp = myNum * 2;
    Temp = sqrt(Temp);
    Temp *= static_cast<unsigned long int>(LenB);
    
    const double fudge2 = (digCount < 45) ? 1.410 :
                          (digCount < 50) ? 1.440 :
                          (digCount < 60) ? 1.510 : 
                          (digCount < 65) ? 1.515 : 1.540;
    
    const int theCut = std::ceil(100.0 * fudge2 *
                                 static_cast<double>(mpz_sizeinbase(Temp.get_mpz_t(), 10)));
    
    std::vector<int> LnFB(facSize);
    
    for (std::size_t i = 0; i < facSize; ++i) {
        mpzFacBase.push_back(facBase[i]);
        LnFB[i] = std::floor(100.0 * std::log(static_cast<double>(facBase[i])));
    }
    
    Temp = sqrt(myNum);
    Temp -= static_cast<unsigned long int>(LenB);
    const int minPrime = static_cast<int>(mpz_sizeinbase(Temp.get_mpz_t(), 10) * 2);
    
    const auto it = std::find_if(facBase.cbegin(), facBase.cend(),
                                 [minPrime](int f) {return f > minPrime;});
    
    const std::size_t strt = std::distance(facBase.cbegin(), it) + 1u;
    
    mpz_class LowBound;
    LowBound = -1 * static_cast<int>(LenB);
    
    mpz_class NextPrime;
    mpz_mul_2exp(NextPrime.get_mpz_t(), myNum.get_mpz_t(), 1);
    NextPrime = sqrt(NextPrime);
    NextPrime /= static_cast<unsigned long int>(LenB);
    NextPrime = sqrt(NextPrime);
    
    if (cmp(NextPrime, facBase.back()) < 0)
        NextPrime = facBase.back();
    
    Polynomial myPoly(mpzContainerSize, facSize, bShowStats, myNum);
    bool xtraTime = true;
    
    if (IsParallel) {
        myPoly.FactorParallel(SieveDist, facBase, LnFB, mpzFacBase, NextPrime,
                              LowBound, myNum, theCut, DoubleLenB, vecMaxSize,
                              strt, checkPoint0, nThreads);

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
        myPoly.FactorSerial(SieveDist, facBase, LnFB, mpzFacBase, NextPrime,
                            LowBound, myNum, theCut, DoubleLenB, vecMaxSize,
                            strt, checkPoint0);

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
