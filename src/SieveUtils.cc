#include "SieveUtils.h"

// Getting quadratic residues. See tonellishanks.cc for more
// details. The array "TS" was used here to make the code
// more concise and since everything will be stored in
// SieveDist, TS can easily be cleared from memory when done.
std::vector<std::size_t> setSieveDist(mpz_t myNum, mpz_t *const TS,
                                      const std::vector<std::size_t> &facBase,
                                      std::size_t facSize) {
    
    std::vector<std::size_t> SieveDist(facSize * 2, 0u);
    SieveDist[0] = SieveDist[1] = 1;
    
    for (std::size_t i = 1, row = 2,
         myAns1 = 0, myAns2 = 0; i < facSize; ++i, row += 2) {
        
        TonelliShanksC(myNum, facBase[i],
                       myAns1, myAns2, TS);
        
        SieveDist[row] = myAns1;
        SieveDist[row + 1] = myAns2;
    }
    
    return SieveDist;
}

std::vector<std::size_t> outersect(std::vector<std::size_t> &x, std::vector<std::size_t> &y) {
    
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    
    const std::size_t lenX = x.size();
    const std::size_t lenY = y.size();
    std::vector<std::size_t> v(lenX + lenY);
    
    auto it = std::set_difference(x.begin(), x.end(), y.begin(), y.end(), v.begin());
    auto it2 = std::set_difference(y.begin(), y.end(), x.begin(), x.end(), it);
    v.resize(it2 - v.begin());
    
    return v;
}

// xi = di * 2 ^ ex  ==> log(xi) = log(di) + ex * log(2) :
double log_mpz_t(mpz_t x) {
    signed long int ex;
    const double di = mpz_get_d_2exp(&ex, x);
    return std::log(di) + M_LN2 * static_cast<double>(ex);
}

uint64_t makeKey(mpz_t x) {
    const double myLog = log_mpz_t(x);
    uint64_t quotient = static_cast<uint64_t>(Significand53 / myLog);
    return quotient;
}

std::vector<uint8_t> myIntToBit(std::size_t x, std::size_t dig) {
    
    std::vector<uint8_t> binaryVec(dig);
    
    for (std::size_t i = 0; x > 0; ++i) {
        binaryVec[i] = x % 2;
        x >>= 1;
    }
    
    return binaryVec;
}

std::vector<std::size_t> getPrimesQuadRes(mpz_t myN, double n) {
    
    const std::size_t uN = n;
    std::vector<char> primes(uN + 1, 1);
    std::vector<std::size_t> myps;
    
    myps.reserve(n * 2.0 / std::log(n));
    const std::size_t fsqr = std::floor(std::sqrt(n));
    
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
    
    for (std::size_t j = 3; j <= n; j += 2) {
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
    
    mpz_clear(jmpz); mpz_clear(temp); mpz_clear(test);
    return myps;
}

void sieveLists(std::size_t facSize, const std::vector<std::size_t> &FBase,
                std::size_t LenB2, mpz_t *const sqrDiff,
                const std::vector<double> &LnFB,
                std::vector<double> &myLogs,
                std::vector<bool> &indexDiv,
                std::size_t minPrime,
                const std::vector<std::size_t> &polySieveD,
                mpz_t lowerBound) {
    
    std::fill(myLogs.begin(), myLogs.end(), 0.0);
    const std::size_t indexLimit = LenB2 * facSize;
    const std::size_t strt = (mpz_even_p(sqrDiff[0]) != 0) ? 0 : 1;
    
    for (std::size_t j = strt * facSize, rowJump = 2 * facSize; j < indexLimit; j += rowJump)
        indexDiv[j] = true;
    
    mpz_t modTest;
    mpz_init(modTest);
    
    for (std::size_t i = 1, row = 2; i < facSize; ++i, row += 2) {
        const std::size_t rowJump = FBase[i] * facSize;
        
        mpz_mod_ui(modTest, lowerBound, FBase[i]);
        int64_t q = mpz_get_si(modTest);
        
        mpz_mod_ui(modTest, sqrDiff[0], FBase[i]);
        int64_t myStart0 = mpz_get_si(modTest);
        
        int64_t myStart1 = 0;
        int64_t myMin = 0;
        int64_t myMax = 0;
        
        if (polySieveD[row] > polySieveD[row + 1]) {
            myMax = polySieveD[row];
            myMin = polySieveD[row + 1];
        } else {
            myMin = polySieveD[row];
            myMax = polySieveD[row + 1];
        }
        
        if (myStart0 == 0) {
            myStart0 = 0;
            
            for (std::size_t j = 1; j < LenB2; ++j) {
                mpz_mod_ui(modTest, sqrDiff[j], FBase[i]);
                
                if (mpz_cmp_ui(modTest, 0) == 0) {
                    myStart1 = j;
                    break;
                }
            }
        } else {
            if (myMin > q) {
                myStart0 = myMin - q;
            } else if (mpz_sgn(lowerBound) < 0) {
                myStart0 = -1 * (q - FBase[i] - myMin);
            } else {
                myStart0 = FBase[i] - ((myMax + q) % FBase[i]);
            }
            
            if (myMax > q) {
                myStart1 = myMax - q;
            } else if (mpz_sgn(lowerBound) < 0) {
                myStart1 = -1 * (q - FBase[i] - myMax);
            } else {
                myStart1 = FBase[i] - ((myMin + q) % FBase[i]);
            }
        }
        
        for (std::size_t j = myStart0 * facSize + i; j < indexLimit; j += rowJump)
            indexDiv[j] = true;
        
        if (FBase[i] > minPrime)
            for (std::size_t j = myStart0; j < LenB2; j += FBase[i])
                myLogs[j] += LnFB[i];
        
        for (std::size_t j = myStart1 * facSize + i; j < indexLimit; j += rowJump)
            indexDiv[j] = true;
        
        if (FBase[i] > minPrime)
            for (std::size_t j = myStart1; j < LenB2; j += FBase[i])
                myLogs[j] += LnFB[i];
    }
    
    mpz_clear(modTest);
}
