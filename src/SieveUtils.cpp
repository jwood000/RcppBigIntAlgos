#include "SieveUtils.h"

// Getting quadratic residues. See tonellishanks.cc for more details
std::vector<std::size_t> SetSieveDist(const std::vector<std::size_t> &facBase, mpz_t myNum) {
    
    const std::size_t facSize = facBase.size();
    std::vector<std::size_t> SieveDist(facSize * 2, 0u);
    SieveDist[0] = SieveDist[1] = 1;
    
    mpz_class p, TS_1, TS_2;
    
    for (std::size_t i = 1; i < facSize; ++i) {
        p = facBase[i];
        TonelliShanksC(myNum, p.get_mpz_t(), TS_1.get_mpz_t());
        
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
            mpz_sub_ui(temp, jmpz, 1u);
            mpz_div_2exp(temp, temp, 1);
            mpz_powm(test, myN, temp, jmpz);
            
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
                    std::vector<int> &myStart, mpz_t firstSqrDiff, mpz_t VarA,
                    mpz_t VarB, mpz_t LowBound, std::size_t strt) {
    
    mpz_class Temp, AUtil;
    const int vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, row = strt * 2, facSize = FBase.size(); i < facSize; ++i, row += 2) {
        
        Temp = FBase[i];
        mpz_invert(AUtil.get_mpz_t(), VarA, Temp.get_mpz_t());
        
        mpz_ui_sub(Temp.get_mpz_t(), SieveDist[row], VarB);
        Temp *= AUtil;
        Temp %= FBase[i];
        int myMin = Temp.get_si();
        
        mpz_ui_sub(Temp.get_mpz_t(), SieveDist[row + 1], VarB);
        Temp *= AUtil;
        Temp %= FBase[i];
        int myMax = Temp.get_si();
        
        mpz_mod_ui(Temp.get_mpz_t(), LowBound, FBase[i]);
        const int q = Temp.get_si();
        
        mpz_mod_ui(Temp.get_mpz_t(), firstSqrDiff, FBase[i]);
        int myStart0 = Temp.get_si();
        int myStart1 = 0;
        
        if (myMin > myMax)
            std::swap(myMin, myMax);
        
        if (myStart0 == 0) {
            myStart1 = (q == myMin) ? (myMax - myMin) : static_cast<int>(FBase[i]) - (myMax - myMin);
        } else {
            myStart0 = (myMin > q) ? myMin - q :
                       (mpz_sgn(LowBound) < 0) ? static_cast<int>(FBase[i]) + myMin - q :
                       static_cast<int>(FBase[i]) - static_cast<int>((myMax + q) % static_cast<int>(FBase[i]));
            
            myStart1 = (myMax > q) ? myMax - q :
                       (mpz_sgn(LowBound) < 0) ? static_cast<int>(FBase[i]) + myMax - q :
                       static_cast<int>(FBase[i]) - static_cast<int>((myMin + q) % static_cast<int>(FBase[i]));
        }
        
        for (int j = myStart0; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (int j = myStart1; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        myStart[row] = ((myStart0 - vecMaxSize) % static_cast<int>(FBase[i])) + static_cast<int>(FBase[i]);
        myStart[row + 1] = ((myStart1 - vecMaxSize) % static_cast<int>(FBase[i])) + static_cast<int>(FBase[i]);
    }
}

void SieveLists(const std::vector<std::size_t> &FBase, const std::vector<int> &LnFB,
                std::vector<int> &myStart, std::vector<int> &myLogs, std::size_t strt) {
    
    std::fill(myLogs.begin(), myLogs.end(), 0);
    const int vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, row = strt * 2, facSize = FBase.size(); i < facSize; ++i, row += 2) {
        for (int j = myStart[row]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (int j = myStart[row + 1]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        myStart[row] = ((myStart[row] - vecMaxSize) % static_cast<int>(FBase[i])) + static_cast<int>(FBase[i]);
        myStart[row + 1] = ((myStart[row + 1] - vecMaxSize) % static_cast<int>(FBase[i])) + static_cast<int>(FBase[i]);
    }
}

void SieveListsFinal(const std::vector<std::size_t> &FBase, const std::vector<int> &LnFB,
                     const std::vector<int> &myStart, std::vector<int> &myLogs, std::size_t strt) {
    
    std::fill(myLogs.begin(), myLogs.end(), 0);
    const std::size_t vecMaxSize = myLogs.size();
    
    for (std::size_t i = strt, facSize = FBase.size(); i < facSize; ++i) {
        for (std::size_t j = myStart[i * 2]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
        
        for (std::size_t j = myStart[i * 2 + 1]; j < vecMaxSize; j += FBase[i])
            myLogs[j] += LnFB[i];
    }
}

void SinglePoly(const std::vector<std::size_t> &SieveDist,
                const std::vector<std::size_t> &facBase, const std::vector<int> &LnFB,
                vec2dint &powsOfSmooths, vec2dint &powsOfPartials,
                std::vector<std::size_t> &coFactorIndexVec, std::vector<int> &myStart,
                hash64vec &partFactorsMap, hash64mpz &partIntvlMap, hash64size_t &keepingTrack,
                std::vector<mpz_class> &smoothInterval, std::vector<mpz_class> &largeCoFactors,
                std::vector<mpz_class> &partialInterval, mpz_class NextPrime,
                mpz_class LowBound, mpz_class myNum, std::size_t &nPartial,
                std::size_t &nSmooth, std::size_t &coFactorInd, int theCut,
                int DoubleLenB, std::size_t mpzFacSize, int vecMaxSize, std::size_t strt) {
    
    mpz_class VarA, VarB, VarC, AUtil, Temp, IntVal;
    TonelliShanksC(myNum.get_mpz_t(), NextPrime.get_mpz_t(), VarC.get_mpz_t());
    
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
    
    SieveListsInit(facBase, LnFB, SieveDist, myLogs, myStart, IntVal.get_mpz_t(),
                   VarA.get_mpz_t(), VarB.get_mpz_t(), LowBound.get_mpz_t(), strt);
    
    for (int chunk = 0; chunk < DoubleLenB; chunk += vecMaxSize) {
        std::vector<int> largeLogs;
        
        for (int i = 0, myLogSize = myLogs.size(); i < myLogSize; ++i)
            if (myLogs[i] > theCut)
                largeLogs.push_back(i + chunk);
            
        for (auto lrgLog: largeLogs) {
            std::vector<int> primeIndexVec;
            const int myIntVal = intLowBound + lrgLog;
            
            Temp = VarB * 2 * myIntVal;
            Temp += VarC;
            AUtil = VarA * myIntVal;
            AUtil *= myIntVal;
            IntVal = AUtil + Temp;

            // Add the index referring to A^2.. (i.e. add it twice)
            primeIndexVec.insert(primeIndexVec.end(), 2, static_cast<int>(mpzFacSize));

            // If Negative, we push zero (i.e. the index referring to -1)
            if (mpz_sgn(IntVal.get_mpz_t()) < 0) {
                mpz_abs(IntVal.get_mpz_t(), IntVal.get_mpz_t());
                primeIndexVec.push_back(0);
            }

            for (int j = 0, facSize = facBase.size(); j < facSize; ++j) {
                while (mpz_divisible_ui_p(IntVal.get_mpz_t(), facBase[j])) {
                    IntVal /= facBase[j];
                    primeIndexVec.push_back(j + 1);
                }
            }

            Temp = VarA * myIntVal;
            
            if (cmp(IntVal, 1u) == 0) {
                // Found a smooth number
                smoothInterval.push_back(Temp + VarB);
                powsOfSmooths.push_back(primeIndexVec);
                ++nSmooth;
            } else if (mpz_cmp_d(IntVal.get_mpz_t(), Significand53) < 0) {
                const uint64_t myKey = static_cast<uint64_t>(IntVal.get_d());
                const auto pFacIt = partFactorsMap.find(myKey);

                if (pFacIt != partFactorsMap.end()) {
                    const auto trackIt = keepingTrack.find(myKey);

                    if (trackIt != keepingTrack.end()) {
                        coFactorIndexVec.push_back(trackIt->second);
                    } else {
                        keepingTrack[myKey] = coFactorInd;
                        largeCoFactors.push_back(IntVal);
                        coFactorIndexVec.push_back(coFactorInd);
                        ++coFactorInd;
                    }

                    primeIndexVec.insert(primeIndexVec.begin(),
                                         pFacIt->second.cbegin(), pFacIt->second.cend());

                    powsOfPartials.push_back(primeIndexVec);
                    const auto intervalIt = partIntvlMap.find(myKey);

                    Temp += VarB;
                    partialInterval.push_back(Temp * intervalIt->second);

                    partFactorsMap.erase(pFacIt);
                    partIntvlMap.erase(intervalIt);
                    ++nPartial;
                } else {
                    partFactorsMap[myKey] = primeIndexVec;
                    partIntvlMap[myKey] = Temp + VarB;
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
}
