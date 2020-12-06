#ifndef MULTI_POLY_H
#define MULTI_POLY_H

#include "SieveUtils.h"

namespace MPQS {
    
    class SieveIndex {
    private:
        int ind_1;
        int ind_2;
        int offset;
        
    public:
        bool IsDivisible(std::uint32_t myPrime, std::uint32_t ind) const;
        void InitialSet(int temp, int q, int myMin,
                        int myMax, int myPrime, int vecMaxSize);
      
        void SmallSieve(std::vector<logType> &myLogs, int vecMaxSize,
                        int myPrime, logType LnFB);
        
        void LargeSieve(std::vector<logType> &myLogs,
                        int vecMaxSize, int myPrime, logType LnFB);
    };
    
    void SinglePoly(const std::vector<std::size_t> &SieveDist,
                    const std::vector<int> &facBase,
                    const std::vector<logType> &LnFB, vec2dint &powsOfSmooths,
                    vec2dint &powsOfPartials, std::vector<SieveIndex> &myStart,
                    hash64vec &partFactorsMap, hash64mpz &partIntvlMap,
                    std::vector<mpz_class> &smoothInterval,
                    std::vector<std::uint64_t> &largeCoFactors,
                    std::vector<mpz_class> &partialInterval,
                    const mpz_class &NextPrime, const mpz_class &myNum,
                    int LowBound, logType theCut, int TwiceLenB, int mpzFacSize,
                    int vecMaxSize, std::size_t strt, std::size_t vecMaxStrt);
}

#endif
