#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "Cpp14MakeUnique.h"
#include "MultPoly.h"
#include "StatsUtils.h"

using typeTimePoint = std::chrono::time_point<std::chrono::steady_clock>;

class Polynomial {
private:
    vec2dint powsOfSmooths;
    vec2dint powsOfPartials;
    
    std::vector<int> myStart;
    hash64vec partFactorsMap;
    hash64mpz partIntvlMap;
    
    std::vector<mpz_class> smoothInterval;
    std::vector<uint64_t> largeCoFactors;
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
                     std::vector<uint64_t> &largeCoFactorsBig, 
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
    
public:
    
    Polynomial(std::size_t _facSize);
    Polynomial(std::size_t _facSize, bool _bShowStats, const mpz_class &myNum);
    
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
