#include "SieveUtils.h"

// Getting quadratic residues. See tonellishanks.cc for more details
std::vector<std::size_t> SetSieveDist(const std::vector<int> &facBase,
                                      const mpz_class &myNum) {
    
    const std::size_t facSize = facBase.size();
    std::vector<std::size_t> SieveDist(facSize * 2, 0u);
    mpz_class p, TS_1;
    
    for (std::size_t i = 1; i < facSize; ++i) {
        p = facBase[i];
        TonelliShanksC(myNum, p, TS_1);
        SieveDist[i * 2] = TS_1.get_ui();
        TS_1 = p - TS_1;
        SieveDist[i * 2 + 1] = TS_1.get_ui();
    }
    
    return SieveDist;
}

std::vector<int> GetPrimesQuadRes(const mpz_class &myN, double LimB, double fudge1,
                                  double sqrLogLog, std::size_t myTarget) {
    
    const std::size_t uN = LimB;
    std::vector<char> primes(uN + 1, 1);
    std::vector<int> myps;
    
    myps.reserve(LimB * 2.0 / std::log(LimB));
    const std::size_t fsqr = std::floor(std::sqrt(LimB));
    
    for (std::size_t j = 4; j <= uN; j += 2)
        primes[j] = 0;
    
    for (std::size_t lastP = 3; lastP <= fsqr;) {
        for (std::size_t j = lastP * lastP; j <= uN; j += 2 * lastP)
            primes[j] = 0;
        
        std::size_t k = lastP + 2;
        std::size_t ind = 2;
        
        while (!primes[k]) {
            k += 2;
            ind += 2;
        }
        
        lastP += ind;
    }
    
    myps.push_back(2u);
    mpz_class currP, nextP;
    
    for (int j = 3; j <= static_cast<int>(uN); j += 2) {
        if (primes[j]) {
            currP = j;
            
            if (mpz_legendre(myN.get_mpz_t(), currP.get_mpz_t()) == 1)
                myps.push_back(j);
        }
    }
    
    while (myps.size() < myTarget) {
        fudge1 += 0.005;
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        
        currP = myps.back();
        mpz_nextprime(nextP.get_mpz_t(),
                      currP.get_mpz_t());
        
        while (cmp(nextP, LimB) < 0) {
            currP = nextP;
            mpz_nextprime(nextP.get_mpz_t(),
                          currP.get_mpz_t());
            
            if (mpz_legendre(myN.get_mpz_t(), currP.get_mpz_t()) == 1)
                myps.push_back(currP.get_si());
        }
    }
    
    // Ensure that the facBase is utilized to most efficiently
    // based off the size of vecMaxSize (The size of each segment)
    if (myps.back() > (4 *L1Cache)) {
        double myDec = std::fmod(static_cast<double>(myps.back()) / static_cast<double>(L1Cache), 1.0);
        
        if (myDec > 0.5) {
            const int biggerTarget = ((myps.back() + L1Cache - 1) / L1Cache) * L1Cache;
            
            while (myps.back() < biggerTarget) {
                currP = nextP;
                mpz_nextprime(nextP.get_mpz_t(),
                              currP.get_mpz_t());
                
                if (mpz_legendre(myN.get_mpz_t(), currP.get_mpz_t()) == 1)
                    myps.push_back(currP.get_si());
            }
            
            myps.pop_back();
        } else {
            const int smallerTarget = (myps.back() / L1Cache) * L1Cache;
            
            while (myps.back() > smallerTarget) {
                myps.pop_back();
            }
        }
    }
    
    return myps;
}

void SieveListsInit(const std::vector<int> &FBase, const std::vector<int> &LnFB,
                    const std::vector<std::size_t> &SieveDist, std::vector<int> &myLogs,
                    std::vector<int> &myStart, const mpz_class &firstSqrDiff, const mpz_class &VarA,
                    const mpz_class &VarB, const mpz_class &LowBound, std::size_t strt) {
    
    mpz_class Temp, AUtil;
    
    for (std::size_t i = strt, row = strt * 2,
         vecMaxSize = myLogs.size(), facSize = FBase.size(); i < facSize; ++i, row += 2) {
        
        Temp = static_cast<unsigned long int>(FBase[i]);
        mpz_invert(AUtil.get_mpz_t(), VarA.get_mpz_t(), Temp.get_mpz_t());
        
        mpz_ui_sub(Temp.get_mpz_t(), SieveDist[row], VarB.get_mpz_t());
        Temp *= AUtil;
        mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), FBase[i]);
        int myMin = Temp.get_si();
        
        mpz_ui_sub(Temp.get_mpz_t(), SieveDist[row + 1], VarB.get_mpz_t());
        Temp *= AUtil;
        mpz_mod_ui(Temp.get_mpz_t(), Temp.get_mpz_t(), FBase[i]);
        int myMax = Temp.get_si();
        
        mpz_mod_ui(Temp.get_mpz_t(), LowBound.get_mpz_t(), FBase[i]);
        const int q = Temp.get_si();
        
        mpz_mod_ui(Temp.get_mpz_t(), firstSqrDiff.get_mpz_t(), FBase[i]);
        int myStart0 = Temp.get_si();
        int myStart1 = 0;
        
        if (myMin > myMax)
            std::swap(myMin, myMax);
        
        if (myStart0 == 0) {
            myStart1 = (q == myMin) ? (myMax - myMin) : FBase[i] - (myMax - myMin);
        } else {
            myStart0 = (myMin > q) ? myMin - q : FBase[i] + myMin - q;
            myStart1 = (myMax > q) ? myMax - q : FBase[i] + myMax - q;
        }
        
        for (std::size_t j = myStart0; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (std::size_t j = myStart1; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        myStart[row] = ((myStart0 - static_cast<int>(vecMaxSize)) % FBase[i]) + FBase[i];
        myStart[row + 1] = ((myStart1 - static_cast<int>(vecMaxSize)) % FBase[i]) + FBase[i];
    }
}

void SinglePoly(const std::vector<std::size_t> &SieveDist,
                const std::vector<int> &facBase, const std::vector<int> &LnFB,
                vec2dint &powsOfSmooths, vec2dint &powsOfPartials,
                std::vector<int> &myStart, hash64vec &partFactorsMap,
                hash64mpz &partIntvlMap, std::vector<mpz_class> &smoothInterval,
                std::vector<uint64_t> &largeCoFactors,
                std::vector<mpz_class> &partialInterval,
                const mpz_class &NextPrime, const mpz_class &LowBound,
                const mpz_class &myNum, int theCut,int DoubleLenB,
                int mpzFacSize, int vecMaxSize, std::size_t strt) {
    
    mpz_class VarA, VarB, VarC, AUtil, Temp, IntVal;
    TonelliShanksC(myNum, NextPrime, VarC);
    
    Temp = VarC * 2u;
    mpz_invert(Temp.get_mpz_t(), Temp.get_mpz_t(), NextPrime.get_mpz_t());
    mpz_pow_ui(VarB.get_mpz_t(), VarC.get_mpz_t(), 2u);
    
    VarB = myNum - VarB;
    VarB *= Temp;
    VarB += VarC;
    
    mpz_pow_ui(VarA.get_mpz_t(), NextPrime.get_mpz_t(), 2u);
    VarB %= VarA;
    
    mpz_pow_ui(VarC.get_mpz_t(), VarB.get_mpz_t(), 2u);
    VarC -= myNum;
    VarC /= VarA;
    
    const int intLowBound = LowBound.get_si();
    
    Temp = VarB * 2 * intLowBound;
    Temp += VarC;
    AUtil = VarA * LowBound;
    AUtil *= LowBound;
    IntVal = AUtil + Temp;
    
    std::vector<int> myLogs(vecMaxSize);
    SieveListsInit(facBase, LnFB, SieveDist, myLogs,
                   myStart, IntVal, VarA, VarB, LowBound, strt);
    
    for (int chunk = 0; chunk < DoubleLenB; chunk += vecMaxSize) {
        std::vector<int> largeLogs;
        
        for (int i = 0, myLogSize = myLogs.size(); i < myLogSize; ++i)
            if (myLogs[i] > theCut)
                largeLogs.push_back(i + chunk);
            
        for (const auto lrgLog: largeLogs) {
            std::vector<int> primeIndexVec;
            const int myIntVal = intLowBound + lrgLog;
            
            Temp = VarB * myIntVal;
            Temp *= 2;
            Temp += VarC;
            AUtil = VarA * myIntVal;
            AUtil *= myIntVal;
            IntVal = AUtil + Temp;

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
        
        if (chunk + 2 * vecMaxSize < DoubleLenB) {
            std::fill(myLogs.begin(), myLogs.end(), 0);
            
            for (std::size_t i = strt, row = strt * 2, facSize = facBase.size(); i < facSize; ++i, row += 2) {
                for (int j = myStart[row]; j < vecMaxSize; j += facBase[i])
                    myLogs[j] += LnFB[i];
                
                for (int j = myStart[row + 1]; j < vecMaxSize; j += facBase[i])
                    myLogs[j] += LnFB[i];
                
                myStart[row] = ((myStart[row] - vecMaxSize) % facBase[i]) + facBase[i];
                myStart[row + 1] = ((myStart[row + 1] - vecMaxSize) % facBase[i]) + facBase[i];
            }
        } else if (chunk + vecMaxSize < DoubleLenB) {
            myLogs.resize(DoubleLenB % vecMaxSize);
            std::fill(myLogs.begin(), myLogs.end(), 0);
            
            for (std::size_t i = strt, facSize = facBase.size(), logSize = myLogs.size(); i < facSize; ++i) {
                for (std::size_t j = myStart[i * 2]; j < logSize; j += facBase[i])
                    myLogs[j] += LnFB[i];
                
                for (std::size_t j = myStart[i * 2 + 1]; j < logSize; j += facBase[i])
                    myLogs[j] += LnFB[i];
            }
        }
    }
}
