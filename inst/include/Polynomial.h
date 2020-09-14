#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "Cpp14MakeUnique.h"
#include "SieveUtils.h"
#include <chrono>

using typeTimePoint = std::chrono::time_point<std::chrono::steady_clock>;

class Polynomial {
private:
    vec2dint powsOfSmooths;
    vec2dint powsOfPartials;
    
    std::vector<std::size_t> coFactorIndexVec;
    std::vector<int> myStart;
    
    hash64vec partFactorsMap;
    hash64mpz partIntvlMap;
    hash64size_t keepingTrack;
    
    std::vector<mpz_class> smoothInterval;
    std::vector<double> largeCoFactors;
    std::vector<mpz_class> partialInterval;
    
    std::size_t nPolys;
    std::size_t nPartial;
    std::size_t nSmooth;
    std::size_t coFactorInd;
    int mpzFacSize;
    
    std::size_t SaPThresh;
    const std::size_t facSize;
    const bool bShowStats;
    
    void MergeMaster(vec2dint &powsOfSmoothsBig, vec2dint &powsOfPartialsBig,
                     std::vector<std::size_t> &coFactorIndexVecBig,
                     hash64vec &partFactorsMapBig, hash64mpz &partIntvlMapBig,
                     hash64size_t &keepingTrackBig, std::vector<mpz_class> &smoothIntervalBig,
                     std::vector<double> &largeCoFactorsBig, 
                     std::vector<mpz_class> &partialIntervalBig);
    
    void SetMpzFacSize(int _mpzFacSize) {mpzFacSize = _mpzFacSize;}
    
    std::size_t GetSmoothSize() const {return smoothInterval.size();};
    
    void SievePolys(const std::vector<std::size_t> &SieveDist,
                    const std::vector<int> &facBase, const std::vector<int> &LnFB, 
                    const std::vector<mpz_class> &mpzFacBase,
                    const mpz_class &LowBound, const mpz_class &myNum,
                    int theCut, int DoubleLenB, int vecMaxSize,
                    std::size_t strt, std::size_t polyLimit);
    
public:
    
    // Polynomial(Polynomial const&) = delete;
    // Polynomial& operator=(Polynomial const&) = delete;
    // Polynomial(Polynomial&&) = delete;
    // //! copy constructor (forbidden)
    // RMonitor(RMonitor const&) = delete;
    // //! move constructor (forbidden)
    // RMonitor(RMonitor&&) = delete;
    // //! copy assignment (forbidden)
    // RMonitor& operator=(RMonitor const&) = delete;
    // //! move assignment (forbidden)
    // RMonitor& operator=(RMonitor &&) = delete;
    
    Polynomial(std::size_t _facSize);
    
    // SaPThresh: Smooth + Partial Threshold
    Polynomial(std::size_t _mpzContainerSize,
               std::size_t _facSize, bool _bShowStats, const mpz_class &myNum);
    
    void InitialParSieve(const std::vector<std::size_t> &SieveDist,
                         const std::vector<int> &facBase, const std::vector<int> &LnFB,
                         std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                         const mpz_class &LowBound, const mpz_class &myNum, int theCut,
                         int DoubleLenB, int vecMaxSize, std::size_t strt,
                         typeTimePoint checkPoint0);
    
    void FactorParallel(const std::vector<std::size_t> &SieveDist,
                        const std::vector<int> &facBase, const std::vector<int> &LnFB,
                        std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                        const mpz_class &LowBound, const mpz_class &myNum, int theCut,
                        int DoubleLenB, int vecMaxSize, std::size_t strt,
                        typeTimePoint checkPoint0, std::size_t nThreads);
    
    void FactorSerial(const std::vector<std::size_t> &SieveDist,
                      const std::vector<int> &facBase, const std::vector<int> &LnFB,
                      std::vector<mpz_class> &mpzFacBase, mpz_class &NextPrime,
                      const mpz_class &LowBound, const mpz_class &myNum, int theCut,
                      int DoubleLenB, int vecMaxSize, std::size_t strt,
                      typeTimePoint checkPoint0);
    
    void GetSolution(const std::vector<mpz_class> &mpzFacBase, 
                     const std::vector<int> &facBase, std::vector<mpz_class> &factors,
                     const mpz_class &mpzNum, std::size_t nThreads,
                     typeTimePoint checkPoint0);
};

#endif
