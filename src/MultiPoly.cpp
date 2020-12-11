#include "MultPoly.h"

namespace MPQS {

    bool SieveIndex::IsDivisible(std::uint32_t myPrime, std::uint32_t ind) const {
        return !((ind_1 + ind) % myPrime && (ind_2 + ind) % myPrime);
    }
    
    void SieveIndex::InitialSet(int temp, int q, int myMin,
                                int myMax, int myPrime, int vecMaxSize) {
        
        ind_1 = (temp) ? ((myMin > q) ? myMin - q : myPrime + myMin - q) : temp;
        ind_2 = (temp) ? ((myMax > q) ? myMax - q : myPrime + myMax - q) :
            (q == myMin) ? (myMax - myMin) : myPrime - (myMax - myMin);
        
        int next_ind_1 = ((ind_1 - vecMaxSize) % myPrime) + myPrime;
        offset = (next_ind_1 > ind_1) ? next_ind_1 - ind_1 : myPrime - ind_1 + next_ind_1;
    }
    
    void SieveIndex::SmallSieve(std::vector<logType> &myLogs, int vecMaxSize,
                                int myPrime, logType LnFB) {
        
        for (int j = ind_1; j < vecMaxSize; j += myPrime)
            myLogs[j] += LnFB;
        
        for (int j = ind_2; j < vecMaxSize; j += myPrime)
            myLogs[j] += LnFB;
        
        ind_1 = (ind_1 + offset >= myPrime) ? ind_1 + offset - myPrime : ind_1 + offset;
        ind_2 = (ind_2 + offset >= myPrime) ? ind_2 + offset - myPrime : ind_2 + offset;
    }

    void SieveIndex::LargeSieve(std::vector<logType> &myLogs,
                                int vecMaxSize, int myPrime, logType LnFB) {
        
        if (ind_1 < vecMaxSize) {
            myLogs[ind_1] += LnFB;
            ind_1 += (myPrime - vecMaxSize);
        } else {
            ind_1 -= vecMaxSize;
        }
        
        if (ind_2 < vecMaxSize) {
            myLogs[ind_2] += LnFB;
            ind_2 += (myPrime - vecMaxSize);
        } else {
            ind_2 -= vecMaxSize;
        }
    }
    
    void SieveListsInit(const std::vector<int> &facBase,
                        const std::vector<logType> &LnFB,
                        const std::vector<std::size_t> &SieveDist,
                        std::vector<logType> &myLogs, std::vector<SieveIndex> &myStart,
                        const mpz_class &firstSqrDiff, const mpz_class &VarA,
                        const mpz_class &VarB, std::size_t strt,
                        int LowBound, int vecMaxSize) {
        
        mpz_class Temp;
        
        for (std::size_t i = strt, facSize = facBase.size(); i < facSize; ++i) {
            auto&& myPrime = facBase[i];
            Temp = VarA % myPrime;
            const int AUtil = int_invert(Temp.get_ui(), myPrime);
            
            mpz_ui_sub(Temp.get_mpz_t(), SieveDist[i], VarB.get_mpz_t());
            Temp *= AUtil;
            mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), myPrime);
            int myMin = Temp.get_si();
            
            mpz_ui_sub(Temp.get_mpz_t(), myPrime - SieveDist[i], VarB.get_mpz_t());
            Temp *= AUtil;
            mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), myPrime);
            int myMax = Temp.get_si();
            
            const int q = (LowBound % myPrime) + myPrime;
            mpz_mod_ui(Temp.get_mpz_t(), firstSqrDiff.get_mpz_t(), myPrime);
            if (myMin > myMax) {std::swap(myMin, myMax);}
            
            myStart[i].InitialSet(Temp.get_si(), q, myMin,
                                  myMax, myPrime, vecMaxSize);
            
            if (myPrime < vecMaxSize) {
                myStart[i].SmallSieve(myLogs, vecMaxSize, myPrime, LnFB[i]);
            } else {
                myStart[i].LargeSieve(myLogs, vecMaxSize, myPrime, LnFB[i]);
            }
        }
    }
    
    void SinglePoly(const std::vector<std::size_t> &SieveDist,
                    const std::vector<int> &facBase,
                    const std::vector<logType> &LnFB, vec2dint &powsOfSmooths,
                    vec2dint &powsOfPartials, std::vector<SieveIndex> &myStart,
                    hash64vec &partFactorsMap, hash64mpz &partIntvlMap,
                    std::vector<mpz_class> &smoothInterval,
                    std::vector<uint64_t> &largeCoFactors,
                    std::vector<mpz_class> &partialInterval,
                    const mpz_class &NextPrime, const mpz_class &myNum,
                    int LowBound, logType theCut, int TwiceLenB, int mpzFacSize,
                    int vecMaxSize, std::size_t strt, std::size_t vecMaxStrt) {
        
        mpz_class VarA, VarB, VarC, IntVal;
        TonelliShanksC(myNum, NextPrime, VarC);
        
        IntVal = VarC * 2u;
        mpz_invert(IntVal.get_mpz_t(),
                   IntVal.get_mpz_t(),
                   NextPrime.get_mpz_t());
        
        VarA = NextPrime * NextPrime;
        VarB = (IntVal * (myNum - VarC * VarC) + VarC) % VarA;
        VarC = (VarB * VarB - myNum) / VarA;
        
        IntVal = LowBound * (VarA * LowBound) + VarB * 2 * LowBound + VarC;
        std::vector<logType> myLogs(vecMaxSize);
        
        SieveListsInit(facBase, LnFB, SieveDist, myLogs ,myStart,
                       IntVal, VarA, VarB, strt, LowBound, vecMaxSize);
        
        for (int chunk = 0; chunk < TwiceLenB; chunk += vecMaxSize) {
            std::vector<int> largeLogs;
            
            for (int i = 0; i < vecMaxSize; ++i)
                if (myLogs[i] > theCut)
                    largeLogs.push_back(i + chunk);
                
            for (const auto lrgLog: largeLogs) {
                std::vector<int> primeIndexVec = {mpzFacSize, mpzFacSize};
                const int myIntVal = LowBound + lrgLog;
                IntVal = (VarA * myIntVal) * myIntVal + (VarB * myIntVal) * 2 + VarC;
                
                // If Negative, we push zero (i.e. the index referring to -1)
                if (sgn(IntVal) < 0) {
                    IntVal = abs(IntVal);
                    primeIndexVec.push_back(0);
                }
                
                for (std::size_t j = 0; j < strt; ++j) {
                    while (mpz_divisible_ui_p(IntVal.get_mpz_t(), facBase[j])) {
                        IntVal /= facBase[j];
                        primeIndexVec.push_back(j + 1);
                    }
                }
                
                for (std::uint32_t j = strt, facSize = facBase.size(),
                     ind = vecMaxSize - (lrgLog - chunk); j < facSize; ++j) {
                    if (myStart[j].IsDivisible(facBase[j], ind)) {
                        do {
                            IntVal /= facBase[j];
                            primeIndexVec.push_back(j + 1);
                        } while (mpz_divisible_ui_p(IntVal.get_mpz_t(), facBase[j]));
                    }
                }
                
                if (cmp(IntVal, 1u) == 0) {
                    // Found a smooth number
                    smoothInterval.push_back(VarA * myIntVal + VarB);
                    powsOfSmooths.push_back(primeIndexVec);
                } else if (cmp(IntVal, Significand53) < 0) {
                    const uint64_t myKey = static_cast<uint64_t>(IntVal.get_d());
                    auto&& pFacIt = partFactorsMap.find(myKey);
                    
                    if (pFacIt != partFactorsMap.end()) {
                        largeCoFactors.push_back(myKey);
                        primeIndexVec.insert(primeIndexVec.begin(),
                                             pFacIt->second.cbegin(),
                                             pFacIt->second.cend());
                        
                        powsOfPartials.push_back(primeIndexVec);
                        auto&& intervalIt = partIntvlMap.find(myKey);
                        partialInterval.push_back((VarA * myIntVal + VarB) *
                                                  intervalIt->second);
                        
                        partFactorsMap.erase(pFacIt);
                        partIntvlMap.erase(intervalIt);
                    } else {
                        partFactorsMap[myKey] = primeIndexVec;
                        partIntvlMap[myKey] = VarA * myIntVal + VarB;
                    }
                }
            }
            
            if (chunk + vecMaxSize < TwiceLenB) {
                std::fill(myLogs.begin(), myLogs.end(),
                          static_cast<logType>(0));
                
                for (std::size_t i = strt; i < vecMaxStrt; ++i)
                    myStart[i].SmallSieve(myLogs, vecMaxSize, facBase[i], LnFB[i]);
                
                for (std::size_t i = vecMaxStrt, facSize = facBase.size(); i < facSize; ++i)
                    myStart[i].LargeSieve(myLogs, vecMaxSize, facBase[i], LnFB[i]);
            }
        }
    }
}
