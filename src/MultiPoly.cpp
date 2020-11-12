#include "MultPoly.h"
#include <deque>

namespace MPQS {

    struct bucketType {
        int startInd;
        int FBInd;
    };

    void SieveListsInit(const std::vector<int> &facBase,
                        const std::vector<logType> &LnFB,
                        const std::vector<std::size_t> &SieveDist,
                        std::deque<std::vector<bucketType>> &Buckets,
                        std::vector<logType> &myLogs, std::vector<int> &myStart,
                        const mpz_class &firstSqrDiff, const mpz_class &VarA,
                        const mpz_class &VarB, std::size_t strt, int LowBound) {
        
        mpz_class Temp;
        const int vecMaxSize = myLogs.size();
        
        for (std::size_t i = strt, facSize = facBase.size(); i < facSize; ++i) {
            auto&& myPrime = facBase[i];
            Temp = VarA % myPrime;
            const int AUtil = int_invert(Temp.get_si(), myPrime);
            
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
            int myStart0 = Temp.get_si();
            int myStart1 = 0;
            
            if (myStart0 == 0) {
                myStart1 = (q == myMin) ? (myMax - myMin) : myPrime - (myMax - myMin);
            } else {
                myStart0 = (myMin > q) ? myMin - q : myPrime + myMin - q;
                myStart1 = (myMax > q) ? myMax - q : myPrime + myMax - q;
            }
            
            if (myPrime < vecMaxSize) {
                myStart[i * 2] = myStart0;
                myStart[i * 2 + 1] = myStart1;
                
                for (int row = i * 2; row <= (i * 2 + 1); ++row) {
                    for (int j = myStart[row]; j < vecMaxSize; j += myPrime)
                        myLogs[j] += LnFB[i];
                    
                    myStart[row] = ((myStart[row] - vecMaxSize) % myPrime) + myPrime;
                }
            } else {
                if (myStart0 < vecMaxSize) {
                    myLogs[myStart0] += LnFB[i];
                    myStart0 += myPrime;
                }
                
                if (myStart1 < vecMaxSize) {
                    myLogs[myStart1] += LnFB[i];
                    myStart1 += myPrime;
                }
                
                myStart0 -= vecMaxSize;
                myStart1 -= vecMaxSize;
                
                Buckets[myStart0 / vecMaxSize].push_back(
                        {.startInd = myStart0 % vecMaxSize,
                         .FBInd = static_cast<int>(i)});
                        
                Buckets[myStart1 / vecMaxSize].push_back(
                        {.startInd = myStart1 % vecMaxSize,
                         .FBInd = static_cast<int>(i)});
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
        
        mpz_class VarA, VarB, VarC, IntVal;
        TonelliShanksC(myNum, NextPrime, VarC);
        
        IntVal = VarC * 2u;
        mpz_invert(IntVal.get_mpz_t(), IntVal.get_mpz_t(), NextPrime.get_mpz_t());
        
        VarA = NextPrime * NextPrime;
        VarB = (IntVal * (myNum - VarC * VarC) + VarC) % VarA;
        VarC = (VarB * VarB - myNum) / VarA;
        
        IntVal = LowBound * (VarA * LowBound) + VarB * 2 * LowBound + VarC;
        
        std::vector<logType> myLogs(vecMaxSize);
        std::deque<std::vector<bucketType>> Buckets(2 + facBase.back()/ vecMaxSize, 
                                                        std::vector<bucketType>());
        
        SieveListsInit(facBase, LnFB, SieveDist, Buckets, myLogs,
                       myStart, IntVal, VarA, VarB, strt, LowBound);
        
        for (int chunk = 0; chunk < TwiceLenB; chunk += vecMaxSize) {
            std::vector<int> largeLogs;
            
            for (int i = 0; i < vecMaxSize; ++i)
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
                                             pFacIt->second.cbegin(), pFacIt->second.cend());
                        
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
            
            if (chunk + 2 * vecMaxSize < TwiceLenB) {
                std::fill(myLogs.begin(), myLogs.end(), 0);
                
                for (std::size_t i = strt; i < vecMaxStrt; ++i) {
                    for (int row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row) {
                        for (int j = myStart[row]; j < vecMaxSize; j += myPrime)
                            myLogs[j] += LnFB[i];
                        
                        myStart[row] = ((myStart[row] - vecMaxSize) % myPrime) + myPrime;
                    }
                }
                
                for (const auto &v: Buckets.front()) {
                    myLogs[v.startInd] += LnFB[v.FBInd];
                    int newInd = (v.startInd + facBase[v.FBInd]);
                    
                    Buckets[newInd / vecMaxSize].push_back(
                            {.startInd = newInd % vecMaxSize,
                             .FBInd = v.FBInd});
                }
                
                Buckets.pop_front();
                Buckets.push_back(std::vector<bucketType>());
                
            } else if (chunk + vecMaxSize < TwiceLenB) {
                std::fill(myLogs.begin(), myLogs.end(), 0);
                
                for (std::size_t i = strt; i < vecMaxStrt; ++i)
                    for (int row = i * 2, myPrime = facBase[i]; row <= (i * 2 + 1); ++row)
                        for (int j = myStart[row]; j < vecMaxSize; j += myPrime)
                            myLogs[j] += LnFB[i];
                
                for (const auto &v: Buckets.front())
                    myLogs[v.startInd] += LnFB[v.FBInd];
            }
        }
    }
}
