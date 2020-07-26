#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "SieveUtils.h"
#include <chrono>

class Polynomial {
private:
    vec2dint powsOfSmooths;
    vec2dint powsOfPartials;
    
    std::vector<std::size_t> coFactorIndexVec;
    std::vector<std::size_t> myStart;
    
    hash64vec partFactorsMap;
    hash64mpz partIntvlMap;
    hash64size_t keepingTrack;
    
    std::vector<mpz_class> smoothInterval;
    std::vector<mpz_class> largeCoFactors;
    std::vector<mpz_class> partialInterval;
    
    std::size_t nPolys;
    std::size_t nPartial;
    std::size_t nSmooth;
    std::size_t coFactorInd;
    std::size_t mpzFacSize;
    
    std::size_t SaPThresh;
    const std::size_t polyLimit;
    const bool bShowStats;
    
public:
    
    // SaPThresh: Smooth + Partial Threshold
    Polynomial(std::size_t _mpzContainerSize,
               std::size_t _facSize, bool _bShowStats, mpz_class myNum);
    
    Polynomial(std::size_t _mpzContainerSize, std::size_t _facSize,
               std::size_t _polyLimit, bool _bShowStats);
    
    void SievePolys(const std::vector<std::size_t> &SieveDist,
                    const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
                    const std::vector<mpz_class> &mpzFacBase,
                    mpz_class LowBound, mpz_class myNum, int theCut, int DoubleLenB,
                    int vecMaxSize, std::size_t strt);
    
    void FactorFinish(const std::vector<std::size_t> &SieveDist,
                      const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
                      std::vector<mpz_class> &mpzFacBase, mpz_class NextPrime,
                      mpz_class LowBound, mpz_class myNum,
                      int theCut, int DoubleLenB, int vecMaxSize, std::size_t strt,
                      std::chrono::time_point<std::chrono::steady_clock> checkPoint0);
    
    void GetSolution(const std::vector<mpz_class> &mpzFacBase, 
                     const std::vector<std::size_t> &facBase, mpz_t *const factors,
                     mpz_t mpzNum, std::size_t nThreads,
                     std::chrono::time_point<std::chrono::steady_clock> checkPoint0);
    
    void SetMpzFacSize(std::size_t _mpzFacSize) {mpzFacSize = _mpzFacSize;}
};

#endif
