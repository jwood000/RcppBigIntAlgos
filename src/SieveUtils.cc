#include "SieveUtils.h"

// Getting quadratic residues. See tonellishanks.cc for more
// details. The array "TS" was used here to make the code
// more concise and since everything will be stored in
// SieveDist, TS can easily be cleared from memory when done.
void setSieveDist(mpz_t myNum, const v1d &facBase,
                  std::size_t facSize, v2d &SieveDist) {
    
    mpz_t TS[13];
    
    for (std::size_t i = 0; i < 13; ++i)
        mpz_init(TS[i]);
    
    int pow2, iter1, iter2;
    mpz_set_ui(TS[12], 2);
    
    for (std::size_t i = 1; i < facSize; ++i) {
        mpz_set_ui(TS[5], facBase[i]);
        mpz_set_ui(TS[0], facBase[i]);
        mpz_sub_ui(TS[0], TS[0], 1);
        mpz_set(TS[1], TS[0]);
        mpz_set_ui(TS[7], 2);
        pow2 = static_cast<int>(mpz_scan1(TS[1], 0));
        mpz_div_2exp(TS[1], TS[1], pow2);
        
        if (pow2 == 1) {
            mpz_add_ui (TS[4], TS[5], 1);
            mpz_div_2exp (TS[4], TS[4], 2);
            mpz_powm (TS[2], myNum, TS[4], TS[5]);
            mpz_neg (TS[4], TS[2]);
            mpz_mod (TS[3], TS[4], TS[5]);
        } else {
            mpz_div_2exp (TS[4], TS[0], 1);
            mpz_powm (TS[6], TS[7], TS[4], TS[5]);
            while (mpz_cmp_ui(TS[6], 1) == 0) {
                mpz_add_ui(TS[7], TS[7], 1);
                mpz_powm (TS[6], TS[7], TS[4], TS[5]);
            }
            
            mpz_add_ui(TS[4], TS[1], 1);
            mpz_div_2exp(TS[4], TS[4], 1);
            mpz_powm(TS[10], myNum, TS[4], TS[5]);
            mpz_powm(TS[8], myNum, TS[1], TS[5]);
            mpz_powm(TS[9], TS[7], TS[1], TS[5]);
            
            iter1 = pow2;
            iter2 = 1;
            mpz_mod(TS[11], TS[8], TS[5]);
            
            while ((mpz_cmp_ui(TS[11], 1) != 0) && (iter2 != 0)) {
                iter2 = 0;
                mpz_mod(TS[11], TS[8], TS[5]);
                
                while (mpz_cmp_ui(TS[11], 1) != 0) {
                    ++iter2;
                    mpz_pow_ui(TS[4], TS[12], iter2);
                    mpz_powm(TS[11], TS[8], TS[4], TS[5]);
                }
                
                if (iter2 != 0) {
                    mpz_pow_ui(TS[4], TS[12], iter1 - iter2 - 1);
                    mpz_powm(TS[4], TS[9], TS[4], TS[5]);
                    mpz_mul(TS[4], TS[4], TS[10]);
                    mpz_mod(TS[10], TS[4], TS[5]);
                    
                    mpz_pow_ui(TS[4], TS[12], iter1 - iter2);
                    mpz_powm(TS[9], TS[9], TS[4], TS[5]);
                    
                    mpz_mul(TS[4], TS[8], TS[9]);
                    mpz_mod(TS[8], TS[4], TS[5]);
                    iter1 = iter2;
                }
                
                mpz_set_ui(TS[11], 0);
            }
            
            mpz_set(TS[2], TS[10]);
            mpz_sub(TS[4], TS[5], TS[10]);
            mpz_mod(TS[3], TS[4], TS[5]);
        }
        
        SieveDist[i][0] = static_cast<int64_t>(mpz_get_si(TS[2]));
        SieveDist[i][1] = static_cast<int64_t>(mpz_get_si(TS[3]));
    }
    
    // Finished generating residues.. now free memory
    for (std::size_t i = 0; i < 13; ++i)
        mpz_clear(TS[i]);
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

v1d getPrimesQuadRes(mpz_t myN, double n) {
    
    std::vector<char> primes(n + 1, 1);
    v1d myps;
    
    myps.reserve(n * 2.0 / std::log(n));
    const int fsqr = std::floor(std::sqrt(n));
    
    for (std::size_t j = 4; j <= static_cast<std::size_t>(n); j += 2)
        primes[j] = 0;
    
    for (int lastP = 3; lastP <= fsqr;) {
        for (int j = lastP * lastP; j <= n; j += 2 * lastP)
            primes[j] = 0;
        
        int k = lastP + 2;
        int ind = 2;
        
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
    
    myps.push_back(2);
    
    for (int64_t j = 3; j <= n; j += 2) {
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

v3d sieveLists(int64_t facLim, const v1d &FBase, int64_t vecLen, mpz_t *const sqrDiff) {
    
    v3d outList(facLim, v2d(2, v1d()));
    mpz_t modTest;
    mpz_init(modTest);
    
    for (int64_t i = 1, vStrt1 = 0, vStrt2 = 0; i < facLim; ++i) {
        const uint64_t uiFB = FBase[i];
        const int64_t siFB = FBase[i];
        
        for (int64_t j = 0; j < vecLen; ++j) {
            mpz_mod_ui(modTest, sqrDiff[j], uiFB);
            
            if (mpz_cmp_ui(modTest, 0) == 0) {
                vStrt1 = j;
                break;
            }
        }
        
        for (int64_t j = vStrt1; j < vecLen; j += siFB)
            outList[i][0].push_back(j);
        
        for (int64_t j = vStrt1 + 1; j < vecLen; ++j) {
            mpz_mod_ui(modTest, sqrDiff[j], uiFB);
            
            if (mpz_cmp_ui(modTest, 0) == 0) {
                vStrt2 = j;
                break;
            }
        }
        
        for (int64_t j = vStrt2; j < vecLen; j += siFB)
            outList[i][1].push_back(j);
    }
    
    mpz_clear(modTest);
    return outList;
}
