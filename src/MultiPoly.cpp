#include "MultPoly.h"

namespace MPQS {
    
    void SieveListsInit(const std::vector<int> &facBase, const std::vector<logType> &LnFB,
                        const std::vector<std::size_t> &SieveDist,
                        std::vector<logType> &myLogs, std::vector<int> &myStart,
                        const mpz_class &firstSqrDiff, const mpz_class &VarA,
                        const mpz_class &VarB, std::size_t strt, int LowBound) {
        
        mpz_class Temp;
        const int vecMaxSize = myLogs.size();
        
        for (std::size_t i = strt, facSize = facBase.size(); i < facSize; ++i) {
            Temp = VarA % facBase[i];
            const int AUtil = int_invert(Temp.get_si(), facBase[i]);
            
            mpz_ui_sub(Temp.get_mpz_t(), SieveDist[i], VarB.get_mpz_t());
            Temp *= AUtil;
            mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), facBase[i]);
            int myMin = Temp.get_si();
            
            mpz_ui_sub(Temp.get_mpz_t(), facBase[i] - SieveDist[i], VarB.get_mpz_t());
            Temp *= AUtil;
            mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), facBase[i]);
            int myMax = Temp.get_si();
            
            const int q = (LowBound % facBase[i]) + facBase[i];
            mpz_mod_ui(Temp.get_mpz_t(), firstSqrDiff.get_mpz_t(), facBase[i]);
            myStart[i * 2] = Temp.get_si();
            
            if (myMin > myMax)
                std::swap(myMin, myMax);
            
            if (myStart[i * 2] == 0) {
                myStart[i * 2 + 1] = (q == myMin) ? (myMax - myMin) : facBase[i] - (myMax - myMin);
            } else {
                myStart[i * 2] = (myMin > q) ? myMin - q : facBase[i] + myMin - q;
                myStart[i * 2 + 1] = (myMax > q) ? myMax - q : facBase[i] + myMax - q;
            }
            
            if (facBase[i] < vecMaxSize) {
                for (int row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row) {
                    for (int j = myStart[row]; j < vecMaxSize; j += myPrime)
                        myLogs[j] += LnFB[i];
                    
                    myStart[row] = ((myStart[row] - vecMaxSize) % myPrime) + myPrime;
                }
            } else {
                for (int row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row) {
                    if (myStart[row] < vecMaxSize) {
                        myLogs[myStart[row]] += LnFB[i];
                        myStart[row] = ((myStart[row] - vecMaxSize) % myPrime) + myPrime;
                    } else if (myStart[row] < (2 * vecMaxSize)) {
                        myStart[row] %= vecMaxSize;
                    } else {
                        myStart[row] -= vecMaxSize;
                    }
                }
            }
        }
    }
    
    void SinglePoly(const std::vector<std::size_t> &SieveDist,
                    const std::vector<int> &facBase,
                    const std::vector<logType> &LnFB, vec2dint &powsOfSmooths,
                    vec2dint &powsOfPartials, std::vector<int> &myStart,
                    hash64vec &partFactorsMap, hash64mpz &partIntvlMap,
                    std::vector<mpz_class> &smoothInterval,
                    std::vector<uint64_t> &largeCoFactors,
                    std::vector<mpz_class> &partialInterval,
                    const mpz_class &NextPrime, const mpz_class &myNum,
                    int LowBound, logType theCut, int TwiceLenB, int mpzFacSize,
                    int vecMaxSize, std::size_t strt, std::size_t vecMaxStrt) {
        
        mpz_class VarA, VarB, VarC, Temp, IntVal;
        TonelliShanksC(myNum, NextPrime, VarC);
        
        IntVal = VarC * 2u;
        mpz_invert(IntVal.get_mpz_t(), IntVal.get_mpz_t(), NextPrime.get_mpz_t());
        
        VarA = NextPrime * NextPrime;
        VarB = (IntVal * (myNum - VarC * VarC) + VarC) % VarA;
        VarC = (VarB * VarB - myNum) / VarA;
        
        IntVal = LowBound * (VarA * LowBound) + VarB * 2 * LowBound + VarC;
        std::vector<logType> myLogs(vecMaxSize);
        
        SieveListsInit(facBase, LnFB, SieveDist, myLogs,
                       myStart, IntVal, VarA, VarB, strt, LowBound);
        
        for (int chunk = 0; chunk < TwiceLenB; chunk += vecMaxSize) {
            std::vector<int> largeLogs;
            
            for (int i = 0, myLogSize = myLogs.size(); i < myLogSize; ++i)
                if (myLogs[i] > theCut)
                    largeLogs.push_back(i + chunk);
                
            for (const auto lrgLog: largeLogs) {
                std::vector<int> primeIndexVec;
                const int myIntVal = LowBound + lrgLog;
                IntVal = (VarA * myIntVal) * myIntVal + (VarB * myIntVal) * 2 + VarC;
                
                // Add the index referring to A^2.. (i.e. add it twice)
                primeIndexVec.insert(primeIndexVec.end(), 2, mpzFacSize);
                
                // If Negative, we push zero (i.e. the index referring to -1)
                if (sgn(IntVal) < 0) {
                    IntVal = abs(IntVal);
                    primeIndexVec.push_back(0);
                }
                
                for (int j = 0, facSize = facBase.size(); j < facSize; ++j) {
                    while (mpz_divisible_ui_p(IntVal.get_mpz_t(), facBase[j])) {
                        IntVal /= facBase[j];
                        primeIndexVec.push_back(j + 1);
                    }
                }
                
                Temp = VarA * myIntVal + VarB;
                
                if (cmp(IntVal, 1u) == 0) {
                    // Found a smooth number
                    smoothInterval.push_back(Temp);
                    powsOfSmooths.push_back(primeIndexVec);
                } else if (cmp(IntVal, Significand53) < 0) {
                    const uint64_t myKey = static_cast<uint64_t>(IntVal.get_d());
                    const auto pFacIt = partFactorsMap.find(myKey);
                    
                    if (pFacIt != partFactorsMap.end()) {
                        largeCoFactors.push_back(myKey);
                        primeIndexVec.insert(primeIndexVec.begin(),
                                             pFacIt->second.cbegin(), pFacIt->second.cend());
                        
                        powsOfPartials.push_back(primeIndexVec);
                        const auto intervalIt = partIntvlMap.find(myKey);
                        partialInterval.push_back(Temp * intervalIt->second);
                        
                        partFactorsMap.erase(pFacIt);
                        partIntvlMap.erase(intervalIt);
                    } else {
                        partFactorsMap[myKey] = primeIndexVec;
                        partIntvlMap[myKey] = Temp;
                    }
                }
            }
            
            if (chunk + 2 * vecMaxSize < TwiceLenB) {
                std::fill(myLogs.begin(), myLogs.end(), 0);
                
                for (std::size_t i = strt; i < vecMaxStrt; ++i) {
                    for (int row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row) {
                        for (int j = myStart[row]; j < vecMaxSize; j += myPrime)
                            myLogs[j] += LnFB[i];
                        
                        myStart[row] = ((myStart[row] - vecMaxSize) % myPrime) + myPrime;
                    }
                }
                
                for (std::size_t i = vecMaxStrt, facSize = facBase.size(); i < facSize; ++i) {
                    for (int row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row) {
                        if (myStart[row] < vecMaxSize) {
                            myLogs[myStart[row]] += LnFB[i];
                            myStart[row] = ((myStart[row] - vecMaxSize) % myPrime) + myPrime;
                        } else if (myStart[row] < (2 * vecMaxSize)) {
                            myStart[row] %= vecMaxSize;
                        } else {
                            myStart[row] -= vecMaxSize;
                        }
                    }
                }
            } else if (chunk + vecMaxSize < LowBound) {
                myLogs.resize(LowBound % vecMaxSize);
                std::fill(myLogs.begin(), myLogs.end(), 0);
                
                for (std::size_t i = strt, facSize = facBase.size(), logSize = myLogs.size(); i < facSize; ++i)
                    for (std::size_t row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row)
                        for (std::size_t j = myStart[row]; j < logSize; j += myPrime)
                            myLogs[j] += LnFB[i];
            }
        }
    }
}
