#include "Polynomial.h"
#include <RcppThread.h>

void QuadraticSieve(mpz_t mpzNum, mpz_t *const factors,
                    std::size_t nThreads, bool bShowStats) {
    
    mpz_class myNum(mpzNum);
    const auto checkPoint0 = std::chrono::steady_clock::now();
    const std::size_t digCount = mpz_sizeinbase(myNum.get_mpz_t(), 10);
    const std::size_t bits = mpz_sizeinbase(myNum.get_mpz_t(), 2);
    
    const double lognum = bits / std::log2(std::exp(1.0));
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
    
    const std::size_t myTarget = static_cast<std::size_t>(dblMyTarget);
    
    while (LimB < myTarget) {
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        fudge1 += 0.001;
    }
    
    const std::vector<std::size_t> facBase = GetPrimesQuadRes(myNum.get_mpz_t(), LimB,
                                                              fudge1, sqrLogLog, myTarget);
    
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
    const int DoubleLenB = 2 * LenB + 1;
    
    const std::size_t facSize = facBase.size();
    const int vecMaxSize = std::min(static_cast<int>((1 + facBase.back() / L1Cache) * L1Cache), DoubleLenB);
    const std::vector<std::size_t> SieveDist = SetSieveDist(facBase, myNum.get_mpz_t());
    
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
    std::size_t mpzFSzLIMIT = facSize * 10;
    
    // This array will be passed to solutionSeach.
    std::vector<mpz_class> mpzFacBase;
    mpzFacBase.reserve(mpzFSzLIMIT);
    
    mpz_class Temp;
    Temp = myNum * 2;
    Temp = sqrt(Temp);
    Temp *= LenB;
    
    const double fudge2 = (digCount < 45) ? 1.410 :
                          (digCount < 50) ? 1.440 :
                          (digCount < 60) ? 1.500 : 
                          (digCount < 65) ? 1.515 : 1.530;
    
    const int theCut = std::ceil(100.0 * fudge2 *
                                 static_cast<double>(mpz_sizeinbase(Temp.get_mpz_t(), 10)));
    
    std::vector<int> LnFB(facSize);
    
    for (std::size_t i = 0; i < facSize; ++i) {
        mpzFacBase.push_back(facBase[i]);
        LnFB[i] = std::floor(100.0 * std::log(static_cast<double>(facBase[i])));
    }
    
    Temp = sqrt(myNum);
    Temp -= LenB;
    const std::size_t minPrime = static_cast<std::size_t>(mpz_sizeinbase(Temp.get_mpz_t(), 10) * 2);
    
    const auto it = std::find_if(facBase.cbegin(), facBase.cend(),
                                 [minPrime](std::size_t f) {return f > minPrime;});
    
    const std::size_t strt = std::distance(facBase.cbegin(), it) + 1u;
    
    mpz_class LowBound;
    LowBound = -1 * static_cast<int>(LenB);
    
    mpz_class NextPrime;
    mpz_mul_2exp(NextPrime.get_mpz_t(), myNum.get_mpz_t(), 1);
    NextPrime = sqrt(NextPrime);
    NextPrime /= LenB;
    NextPrime = sqrt(NextPrime);
    
    if (cmp(NextPrime, facBase.back()) < 0)
        NextPrime = facBase.back();
    
    if (nThreads > 1) {
        std::size_t setPolys = 29500;
        
        std::vector<Polynomial> vecPoly(nThreads,
                                        Polynomial(mpzContainerSize / nThreads,
                                                   facSize, setPolys / nThreads, bShowStats));
        
        for (std::size_t i = 0; i < setPolys; ++i) {
            for (bool LegendreTest = true; LegendreTest; ) {
                mpz_nextprime(NextPrime.get_mpz_t(), NextPrime.get_mpz_t());
                
                if (mpz_legendre(myNum.get_mpz_t(), NextPrime.get_mpz_t()) == 1)
                    LegendreTest = false;
            }
            
            mpzFacBase.push_back(NextPrime);
        }
        
        // RcppThread::ThreadPool pool(nThreads);
        std::size_t step = setPolys / nThreads;
        std::vector<std::thread> myThreads;
        
        for (std::size_t i = 0, mySize = facSize; i < nThreads; ++i, mySize += step) {
            vecPoly[i].SetMpzFacSize(mySize);

            myThreads.emplace_back(&Polynomial::SievePolys, vecPoly[i], std::cref(SieveDist),
                                   std::cref(facBase), std::cref(LnFB), std::cref(mpzFacBase),
                                   LowBound, myNum, theCut, DoubleLenB, vecMaxSize, strt);
        }
        
        for (auto &thr: myThreads)
            thr.join();
        
        mpz_set_ui(factors[0], 11u);
        mpz_set_ui(factors[1], 13u);
        
        // SievePolys(const std::vector<std::size_t> &SieveDist,
        //            const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
        //            const std::vector<mpz_class> &mpzFacBase,
        //            mpz_class LowBound, mpz_class myNum, int theCut, int DoubleLenB,
        //            int vecMaxSize, std::size_t strt)
        
        // while (mpz_cmp_ui(factors[0], 0) == 0) {
        //     myPoly.FactorFinish(SieveDist, facBase, LnFB, mpzFacBase, NextPrime,
        //                         LowBound, myNum, theCut, DoubleLenB, vecMaxSize,
        //                         strt, checkPoint0);
        //     
        //     myPoly.GetSolution(mpzFacBase, facBase, factors,
        //                        myNum.get_mpz_t(), nThreads, checkPoint0);
        //     NextPrime = mpzFacBase.back();
        // }
    } else {
        Polynomial myPoly(mpzContainerSize, facSize, bShowStats, myNum);
        
        while (mpz_cmp_ui(factors[0], 0) == 0) {
            myPoly.FactorFinish(SieveDist, facBase, LnFB, mpzFacBase, NextPrime,
                                LowBound, myNum, theCut, DoubleLenB, vecMaxSize,
                                strt, checkPoint0);
            
            myPoly.GetSolution(mpzFacBase, facBase, factors,
                               myNum.get_mpz_t(), nThreads, checkPoint0);
            NextPrime = mpzFacBase.back();
        }
    }
}
