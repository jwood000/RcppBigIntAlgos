#include "SieveUtils.h"

// Getting quadratic residues. See tonellishanks.cc for more details
std::vector<std::size_t> SetSieveDist(const std::vector<int> &facBase,
                                      const mpz_class &myNum) {
    
    const std::size_t facSize = facBase.size();
    std::vector<std::size_t> SieveDist(facSize * 2, 0u);
    SieveDist[0] = SieveDist[1] = 1;
    
    mpz_class p, TS_1, TS_2;
    
    for (std::size_t i = 1; i < facSize; ++i) {
        p = facBase[i];
        TonelliShanksC(myNum, p, TS_1);
        
        SieveDist[i * 2] = TS_1.get_ui();
        
        TS_2 = p - TS_1;
        SieveDist[i * 2 + 1] = TS_2.get_ui();
    }
    
    return SieveDist;
}

std::vector<std::uint8_t> MyIntToBit(std::size_t x, std::size_t dig) {
    
    std::vector<std::uint8_t> binaryVec(dig);
    
    for (std::size_t i = 0; x > 0; ++i) {
        binaryVec[i] = x % 2;
        x >>= 1;
    }
    
    return binaryVec;
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
    
    mpz_t test, jmpz, temp;
    mpz_init(test);
    mpz_init(jmpz);
    mpz_init(temp);
    
    myps.push_back(2u);
    
    for (std::size_t j = 3; j <= uN; j += 2) {
        if (primes[j]) {
            mpz_set_si(jmpz, j);
            mpz_set_si(temp, j);
            mpz_sub_ui(temp, jmpz, 1u);
            mpz_div_2exp(temp, temp, 1);
            mpz_powm(test, myN.get_mpz_t(), temp, jmpz);
            
            if (mpz_cmp_ui(test, 1u) == 0)
                myps.push_back(j);
        }
    }
    
    mpz_clear(jmpz); mpz_clear(temp);
    mpz_clear(test);
    
    mpz_t currP, nextP, resTest, CP1;
    mpz_init(currP); mpz_init(nextP);
    mpz_init(CP1); mpz_init(resTest);
    
    while (myps.size() < myTarget) {
        fudge1 += 0.005;
        LimB = std::exp((0.5 + fudge1) * sqrLogLog);
        mpz_set_ui(currP, myps.back());
        mpz_nextprime(nextP, currP);
        
        while (mpz_cmp_ui(nextP, LimB) < 0) {
            mpz_set(currP, nextP);
            mpz_nextprime(nextP, currP);
            mpz_sub_ui(CP1, currP, 1);
            mpz_div_2exp(CP1, CP1, 1);
            mpz_powm(resTest, myN.get_mpz_t(), CP1, currP);
            
            if (mpz_cmp_ui(resTest, 1) == 0)
                myps.push_back(mpz_get_si(currP));
        }
    }
    
    mpz_clear(currP); mpz_clear(nextP);
    mpz_clear(CP1); mpz_clear(resTest);
    return myps;
}

void SieveListsInit(const std::vector<int> &FBase, const std::vector<int> &LnFB,
                    const std::vector<std::size_t> &SieveDist, std::vector<int> &myLogs,
                    std::vector<int> &myStart, const mpz_class &firstSqrDiff, const mpz_class &VarA,
                    const mpz_class &VarB, const mpz_class &LowBound, std::size_t strt) {
    
    mpz_class Temp, AUtil;
    const int vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, row = strt * 2, facSize = FBase.size(); i < facSize; ++i, row += 2) {
        
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
            myStart0 = (myMin > q) ? myMin - q :
                       (sgn(LowBound) < 0) ? FBase[i] + myMin - q :
                       FBase[i] - ((myMax + q) % FBase[i]);
            
            myStart1 = (myMax > q) ? myMax - q :
                       (sgn(LowBound) < 0) ? FBase[i] + myMax - q :
                       FBase[i] - ((myMin + q) % FBase[i]);
        }
        
        for (int j = myStart0; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (int j = myStart1; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        myStart[row] = ((myStart0 - vecMaxSize) % FBase[i]) + FBase[i];
        myStart[row + 1] = ((myStart1 - vecMaxSize) % FBase[i]) + FBase[i];
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
                const mpz_class &myNum, std::size_t &nPartial,
                std::size_t &nSmooth, int theCut,int DoubleLenB,
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
            
        for (auto lrgLog: largeLogs) {
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
                ++nSmooth;
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
                    ++nPartial;
                } else {
                    partFactorsMap[myKey] = primeIndexVec;
                    partIntvlMap[myKey] = Temp;
                }
            }
        }
        
        if (chunk + 2 * vecMaxSize < DoubleLenB) {
            std::fill(myLogs.begin(), myLogs.end(), 0);
            
            for (std::size_t i = strt, facSize = facBase.size(); i < facSize; ++i) {
                for (int j = myStart[i * 2]; j < vecMaxSize; j += facBase[i])
                    myLogs[j] += LnFB[i];
                
                for (int j = myStart[i * 2 + 1]; j < vecMaxSize; j += facBase[i])
                    myLogs[j] += LnFB[i];
                
                myStart[i * 2] = ((myStart[i * 2] - vecMaxSize) % facBase[i]) + facBase[i];
                myStart[i * 2 + 1] = ((myStart[i * 2 + 1] - vecMaxSize) % facBase[i]) + facBase[i];
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
