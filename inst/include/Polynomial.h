#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "Cpp14MakeUnique.h"
#include "SieveUtils.h"
#include "StatsUtils.h"

using typeTimePoint = std::chrono::time_point<std::chrono::steady_clock>;

class SieveInd {
private:
    int ind_1;
    int ind_2;
    
public:
    bool IsDivisible(std::uint32_t myPrime, std::uint32_t ind) const {
        return !((ind_1 + ind) % myPrime && (ind_2 + ind) % myPrime);
    };
    
    void InitialSet(int temp, int q, int myMin, int myMax, int myPrime) {
        ind_1 = (temp) ? (myMin > q) ? myMin - q : myPrime + myMin - q : temp;
        ind_2 = (temp) ? (myMax > q) ? myMax - q : myPrime + myMax - q :
            (q == myMin) ? (myMax - myMin) : myPrime - (myMax - myMin);
    };
    
    void SmallSieve(std::vector<logType> &myLogs, int vecMaxSize,
                    int myPrime, logType LnFB) {
        for (int j = ind_1; j < vecMaxSize; j += myPrime)
            myLogs[j] += LnFB;
        
        for (int j = ind_2; j < vecMaxSize; j += myPrime)
            myLogs[j] += LnFB;
        
        ind_1 = ((ind_1 - vecMaxSize) % myPrime) + myPrime;
        ind_2 = ((ind_2 - vecMaxSize) % myPrime) + myPrime;
    };
    
    void LargeSieve(std::vector<logType> &myLogs, int vecMaxSize,
                    int myPrime, logType LnFB) {
        if (ind_1 < vecMaxSize) {
            myLogs[ind_1] += LnFB;
            ind_1 = ((ind_1 - vecMaxSize) % myPrime) + myPrime;
        } else {
            ind_1 -= vecMaxSize;
        }
        
        if (ind_2 < vecMaxSize) {
            myLogs[ind_2] += LnFB;
            ind_2 = ((ind_2 - vecMaxSize) % myPrime) + myPrime;
        } else {
            ind_2 -= vecMaxSize;
        }
    };
};

class Polynomial {
private:
    vec2dint powsOfSmooths;
    vec2dint powsOfPartials;
    
    std::vector<SieveInd> myStart;
    std::vector<logType> myLogs;
    hash64vec partFactorsMap;
    hash64mpz partIntvlMap;
    
    std::vector<mpz_class> smoothInterval;
    std::vector<std::uint64_t> largeCoFactors;
    std::vector<mpz_class> partialInterval;
    
    std::size_t nPolys;
    std::size_t nPartial;
    std::size_t nSmooth;
    int mpzFacSize;
    
    std::size_t SaPThresh;
    const std::size_t facSize;
    bool bShowStats;
    
    void MergeMaster(vec2dint &powsOfSmoothsBig, vec2dint &powsOfPartialsBig,
                     hash64vec &partFactorsMapBig, hash64mpz &partIntvlMapBig,
                     std::vector<mpz_class> &smoothIntervalBig,
                     std::vector<std::uint64_t> &largeCoFactorsBig, 
                     std::vector<mpz_class> &partialIntervalBig);
    
    void SetMpzFacSize(int _mpzFacSize) {mpzFacSize = _mpzFacSize;}
    
    void SievePolys(const std::vector<std::size_t> &SieveDist,
                    const std::vector<int> &facBase,
                    const std::vector<logType> &LnFB,
                    const std::vector<mpz_class> &mpzFacBase,
                    const mpz_class &myNum, int LowBound,
                    logType theCut, int TwiceLenB, int vecMaxSize,
                    std::size_t strt, std::size_t vecMaxStrt,
                    std::size_t polyLimit);
    
    void SinglePoly(const std::vector<std::size_t> &SieveDist,
                    const std::vector<int> &facBase,
                    const std::vector<logType> &LnFB,
                    const mpz_class &NextPrime, const mpz_class &myNum,
                    int LowBound, logType theCut, int TwiceLenB,
                    int vecMaxSize, std::size_t strt, std::size_t vecMaxStrt);
    
    void SieveListsInit(const std::vector<int> &facBase,
                        const std::vector<logType> &LnFB,
                        const std::vector<std::size_t> &SieveDist,
                        const mpz_class &firstSqrDiff, const mpz_class &VarA,
                        const mpz_class &VarB, std::size_t strt,
                        int LowBound, int vecMaxSize);
    
public:
    
    Polynomial(std::size_t _facSize, std::size_t _vecMaxSize);
    Polynomial(std::size_t _facSize, std::size_t _vecMaxSize,
               bool _bShowStats, const mpz_class &myNum);
    
    void InitialParSieve(const std::vector<std::size_t> &SieveDist,
                         const std::vector<int> &facBase,
                         const std::vector<logType> &LnFB,
                         std::vector<mpz_class> &mpzFacBase,
                         mpz_class &NextPrime, const mpz_class &myNum,
                         int LowBound, logType theCut, int TwiceLenB,
                         int vecMaxSize, std::size_t strt,
                         std::size_t vecMaxStrt, typeTimePoint checkPoint0);
    
    void FactorParallel(const std::vector<std::size_t> &SieveDist,
                        const std::vector<int> &facBase, 
                        const std::vector<logType> &LnFB,
                        std::vector<mpz_class> &mpzFacBase,
                        mpz_class &NextPrime, const mpz_class &myNum,
                        int LowBound, logType theCut, int TwiceLenB,
                        int vecMaxSize, std::size_t strt,
                        std::size_t vecMaxStrt, typeTimePoint checkPoint0,
                        std::size_t nThreads);
    
    void FactorSerial(const std::vector<std::size_t> &SieveDist,
                      const std::vector<int> &facBase,
                      const std::vector<logType> &LnFB,
                      std::vector<mpz_class> &mpzFacBase,
                      mpz_class &NextPrime, const mpz_class &myNum,
                      int LowBound, logType theCut, int TwiceLenB,
                      int vecMaxSize, std::size_t strt,
                      std::size_t vecMaxStrt, typeTimePoint checkPoint0);
    
    void GetSolution(const std::vector<mpz_class> &mpzFacBase, 
                     const std::vector<int> &facBase, std::vector<mpz_class> &factors,
                     const mpz_class &mpzNum, std::size_t nThreads,
                     typeTimePoint checkPoint0);
    
    bool ContinueToSolution() {return (nSmooth + nPartial) > SaPThresh;};
    void MakeStatsFalse() {bShowStats = false;};
};

#endif
