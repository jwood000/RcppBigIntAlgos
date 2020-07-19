#include "SieveUtils.h"

// Getting quadratic residues. See tonellishanks.cc for more details
std::vector<std::size_t> SetSieveDist(const std::vector<std::size_t> &facBase, mpz_t myNum) {
    
    const std::size_t facSize = facBase.size();
    std::vector<std::size_t> SieveDist(facSize * 2, 0u);
    SieveDist[0] = SieveDist[1] = 1;
    
    mpz_t TS_1, TS_2;
    mpz_init(TS_1); mpz_init(TS_2);
    
    mpz_t p;
    mpz_init(p);
    
    for (std::size_t i = 1, row = 2; i < facSize; ++i, row += 2) {
        mpz_set_ui(p, facBase[i]);
        TonelliShanksC(myNum, p, TS_1);
        
        SieveDist[row] = mpz_get_ui(TS_1);
        
        mpz_sub(TS_2, p, TS_1);
        SieveDist[row + 1] = mpz_get_ui(TS_2);
    }
    
    mpz_clear(TS_1); mpz_clear(TS_2);
    mpz_clear(p);
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

std::vector<std::size_t> GetPrimesQuadRes(mpz_t myN, double LimB, double fudge1,
                                          double sqrLogLog, std::size_t myTarget) {
    
    const std::size_t uN = LimB;
    std::vector<char> primes(uN + 1, 1);
    std::vector<std::size_t> myps;
    
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
            mpz_sub_ui(temp, jmpz, 1);
            mpz_div_2exp(temp, temp, 1);
            mpz_powm(test,myN,temp,jmpz);
            
            if (mpz_cmp_ui(test, 1) == 0)
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
            mpz_powm(resTest, myN, CP1, currP);
            
            if (mpz_cmp_ui(resTest, 1) == 0)
                myps.push_back(mpz_get_ui(currP));
        }
    }
    
    mpz_clear(currP); mpz_clear(nextP);
    mpz_clear(CP1); mpz_clear(resTest);
    return myps;
}

void SieveListsInit(const std::vector<std::size_t> &FBase, const std::vector<int> &LnFB,
                    const std::vector<std::size_t> &SieveDist, std::vector<int> &myLogs,
                    std::vector<std::size_t> &myStart, mpz_t firstSqrDiff, mpz_t VarA,
                    mpz_t VarB, mpz_t LowBound, std::size_t strt) {
    
    mpz_t Temp, AUtil;
    mpz_init(Temp); mpz_init(AUtil);
    const std::size_t vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, row = strt * 2,
         facSize = FBase.size(); i < facSize; ++i, row += 2) {
        
        mpz_set_ui(Temp, FBase[i]);
        mpz_invert(AUtil, VarA, Temp);
        
        mpz_ui_sub(Temp, SieveDist[row], VarB);
        mpz_mul(Temp, Temp, AUtil);
        mpz_mod_ui(Temp, Temp, FBase[i]);
        std::int64_t myMin = mpz_get_si(Temp);
        
        mpz_ui_sub(Temp, SieveDist[row + 1], VarB);
        mpz_mul(Temp, Temp, AUtil);
        mpz_mod_ui(Temp, Temp, FBase[i]);
        std::int64_t myMax = mpz_get_si(Temp);
        
        mpz_mod_ui(Temp, LowBound, FBase[i]);
        const std::int64_t q = mpz_get_si(Temp);
        
        mpz_mod_ui(Temp, firstSqrDiff, FBase[i]);
        std::int64_t myStart0 = mpz_get_si(Temp);
        std::int64_t myStart1 = 0;
        
        if (myMin > myMax)
            std::swap(myMin, myMax);
        
        if (myStart0 == 0) {
            myStart1 = (q == myMin) ? (myMax - myMin) : FBase[i] - (myMax - myMin);
        } else {
            const std::int64_t signedFB = FBase[i];
             
            myStart0 = (myMin > q) ? myMin - q :
                       (mpz_sgn(LowBound) < 0) ? -1 * (q - signedFB - myMin) :
                       signedFB - ((myMax + q) % signedFB);
            
            myStart1 = (myMax > q) ? myMax - q :
                       (mpz_sgn(LowBound) < 0) ? -1 * (q - signedFB - myMax) :
                       signedFB - ((myMin + q) % signedFB);
        }
        
        for (std::size_t j = myStart0; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (std::size_t j = myStart1; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        myStart[row] = (FBase[i] * ((vecMaxSize - myStart0 + FBase[i] - 1u) / FBase[i]) + myStart0) % vecMaxSize;
        myStart[row + 1] = (FBase[i] * ((vecMaxSize - myStart1 + FBase[i] - 1u) / FBase[i]) + myStart1) % vecMaxSize;
    }
    
    mpz_clear(Temp); mpz_clear(AUtil);
}

void SieveLists(const std::vector<std::size_t> &FBase, const std::vector<int> &LnFB,
                std::vector<std::size_t> &myStart, std::vector<int> &myLogs,
                std::size_t strt) {
    
    std::fill(myLogs.begin(), myLogs.end(), 0);
    const std::size_t vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, row = strt * 2, facSize = FBase.size(); i < facSize; ++i, row += 2) {
        for (std::size_t j = myStart[row]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (std::size_t j = myStart[row + 1]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        myStart[row] = (FBase[i] * ((vecMaxSize - myStart[row] + FBase[i] - 1u) / FBase[i]) + myStart[row]) % vecMaxSize;
        myStart[row + 1] = (FBase[i] * ((vecMaxSize - myStart[row + 1] + FBase[i] - 1u) / FBase[i]) + myStart[row + 1]) % vecMaxSize;
    }
}

void SieveListsFinal(const std::vector<std::size_t> &FBase, const std::vector<int> &LnFB,
                     const std::vector<std::size_t> &myStart, std::vector<int> &myLogs,
                     std::size_t strt) {
    
    std::fill(myLogs.begin(), myLogs.end(), 0);
    const std::size_t vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, row = strt * 2, facSize = FBase.size(); i < facSize; ++i, row += 2) {
        for (std::size_t j = myStart[row]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (std::size_t j = myStart[row + 1]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
    }
}

void SinglePoly(const std::vector<std::size_t> &SieveDist,
                const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
                vec2dint &powsOfSmooths, vec2dint &powsOfPartials,
                std::vector<std::size_t> &coFactorIndexVec,
                std::vector<std::size_t> &myStart, hash64vec &partFactorsMap,
                hash64mpz_t &partIntvlMap, hash64size_t &keepingTrack,
                mpz_t *const smoothInterval, mpz_t *const largeCoFactors,
                mpz_t *const partialInterval, mpz_t NextPrime, mpz_t LowBound,
                mpz_t myNum, std::size_t &nPartial, std::size_t &nSmooth,
                std::size_t &coFactorInd, int theCut, int DoubleLenB,
                std::size_t mpzFacSize, int vecMaxSize, std::size_t strt) {
    
    mpz_t VarA, VarB, VarC, AUtil, Temp, IntVal;
    mpz_init(VarA); mpz_init(VarB); mpz_init(VarC);
    mpz_init(AUtil); mpz_init(Temp); mpz_init(IntVal);
    
    TonelliShanksC(myNum, NextPrime, VarC);
    mpz_mul_2exp(Temp, VarC, 1);
    mpz_invert(Temp, Temp, NextPrime);
    mpz_pow_ui(VarB, VarC, 2u);
    
    mpz_sub(VarB, myNum, VarB);
    mpz_mul(VarB, VarB, Temp);
    mpz_add(VarB, VarB, VarC);
    
    mpz_pow_ui(VarA, NextPrime, 2u);
    mpz_mod(VarB, VarB, VarA);
    
    mpz_pow_ui(VarC, VarB, 2u);
    mpz_sub(VarC, VarC, myNum);
    mpz_divexact(VarC, VarC, VarA);
    
    const int intLowBound = mpz_get_si(LowBound);
    
    mpz_mul_si(Temp, VarB, 2 * intLowBound);
    mpz_add(Temp, Temp, VarC);
    mpz_mul(AUtil, VarA, LowBound);
    mpz_mul(AUtil, AUtil, LowBound);
    mpz_add(IntVal, AUtil, Temp);
    
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

            mpz_mul_si(Temp, VarB, 2 * myIntVal);
            mpz_add(Temp, Temp, VarC);
            mpz_mul_si(AUtil, VarA, myIntVal);
            mpz_mul_si(AUtil, AUtil, myIntVal);
            mpz_add(IntVal, AUtil, Temp);

            // Add the index referring to A^2.. (i.e. add it twice)
            primeIndexVec.insert(primeIndexVec.end(), 2, static_cast<int>(mpzFacSize));

            // If Negative, we push zero (i.e. the index referring to -1)
            if (mpz_sgn(IntVal) < 0) {
                mpz_abs(IntVal, IntVal);
                primeIndexVec.push_back(0);
            }

            for (int j = 0, facSize = facBase.size(); j < facSize; ++j) {
                while (mpz_divisible_ui_p(IntVal, facBase[j])) {
                    mpz_divexact_ui(IntVal, IntVal, facBase[j]);
                    primeIndexVec.push_back(j + 1);
                }
            }

            mpz_mul_si(Temp, VarA, myIntVal);

            if (mpz_cmp_ui(IntVal, 1u) == 0) {
                // Found a smooth number
                mpz_add(smoothInterval[nSmooth], Temp, VarB);
                powsOfSmooths.push_back(primeIndexVec);
                ++nSmooth;
            } else if (mpz_cmp_d(IntVal, Significand53) < 0) {
                const uint64_t myKey = static_cast<uint64_t>(mpz_get_d(IntVal));
                const auto pFacIt = partFactorsMap.find(myKey);

                if (pFacIt != partFactorsMap.end()) {
                    const auto trackIt = keepingTrack.find(myKey);

                    if (trackIt != keepingTrack.end()) {
                        coFactorIndexVec.push_back(trackIt->second);
                    } else {
                        keepingTrack[myKey] = coFactorInd;
                        mpz_set(largeCoFactors[coFactorInd], IntVal);
                        coFactorIndexVec.push_back(coFactorInd);
                        ++coFactorInd;
                    }

                    primeIndexVec.insert(primeIndexVec.begin(),
                                         pFacIt->second.cbegin(), pFacIt->second.cend());

                    powsOfPartials.push_back(primeIndexVec);
                    const auto intervalIt = partIntvlMap.find(myKey);

                    mpz_add(Temp, Temp, VarB);
                    mpz_mul(partialInterval[nPartial],
                            Temp, intervalIt->second);

                    partFactorsMap.erase(pFacIt);
                    partIntvlMap.erase(intervalIt);
                    ++nPartial;
                } else {
                    partFactorsMap[myKey] = primeIndexVec;
                    mpz_init(partIntvlMap[myKey]);
                    mpz_add(partIntvlMap[myKey], Temp, VarB);
                }
            }
        }
        
        if (chunk + 2 * vecMaxSize < DoubleLenB) {
            SieveLists(facBase, LnFB, myStart, myLogs, strt);
        } else if (chunk + vecMaxSize < DoubleLenB) {
            myLogs.resize(DoubleLenB % vecMaxSize);
            SieveListsFinal(facBase, LnFB, myStart, myLogs, strt);
        }
    }
    
    mpz_clear(VarA); mpz_clear(VarB); mpz_clear(VarC);
    mpz_clear(AUtil); mpz_clear(Temp); mpz_clear(IntVal);
}
