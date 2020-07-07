#include "RSAFactorUtils.h"
#include "GrowMPZArray.h"

// Max number of iterations in the main loop
constexpr std::size_t POLLARD_RHO_REPS = 100000u;

std::size_t GetPower(mpz_t nmpz) {
    
    mpz_t testRoot;
    mpz_init(testRoot);
    std::size_t p = 2;
    std::size_t myPow = 1;
    
    for (std::size_t i = 0; i < primesDiffPR.size(); ) {
        mpz_root(testRoot, nmpz, p);
        mpz_pow_ui(testRoot, testRoot, p);
        
        if (mpz_cmp(testRoot, nmpz) == 0) {
            std::size_t powPrime = p;
            bool keepGoing = true;
            
            // Determine if root is a power of a prime
            while (keepGoing) {
                powPrime *= p;
                mpz_root(testRoot, nmpz, powPrime);
                mpz_pow_ui(testRoot, testRoot, powPrime);
                
                if (mpz_cmp(testRoot, nmpz) != 0)
                    keepGoing = false;
            }
            
            powPrime /= p;
            myPow *= powPrime;
            mpz_root(nmpz, nmpz, powPrime);
        }
        
        p += primesDiffPR[i++];
        if (!mpz_perfect_power_p(nmpz))
            break;
    }
    
    // This means the powers involved are very large
    if (mpz_perfect_power_p(nmpz)) {
        mpz_t myNextP;
        mpz_init(myNextP);
        mpz_init_set_ui(myNextP, p);
        
        for (;;) {
            mpz_root(testRoot, nmpz, p);
            mpz_pow_ui(testRoot, testRoot, p);
            
            if (mpz_cmp(testRoot, nmpz) == 0) {
                std::size_t powPrime = p;
                bool keepGoing = true;
                
                // Determine if root is a power of a prime
                while (keepGoing) {
                    powPrime *= p;
                    mpz_root(testRoot, nmpz, powPrime);
                    mpz_pow_ui(testRoot, testRoot, powPrime);
                    
                    if (mpz_cmp(testRoot, nmpz) != 0)
                        keepGoing = false;
                }
                
                powPrime /= p;
                myPow *= powPrime;
                mpz_root(nmpz, nmpz, powPrime);
            }
            
            mpz_nextprime(myNextP, myNextP);
            p = mpz_get_ui(myNextP);
            
            if (!mpz_perfect_power_p(nmpz))
                break;
        }
        
        mpz_clear(myNextP);
    }
    
    mpz_clear(testRoot);
    return myPow;
}

int PollardRhoWithConstraint(mpz_t n, std::size_t a, mpz_t *const factors, std::size_t& numPs,
                             std::vector<std::size_t>& myLens, std::size_t myLimit,
                             std::size_t powMultiplier, std::size_t arrayMax,
                             std::vector<std::size_t>& xtraRecFacs) {
    mpz_t x, z, y, P;
    mpz_t t, t2;
    int returnVal = 0;
    
    mpz_init(t);
    mpz_init(t2);
    mpz_init_set_si(y, 2);
    mpz_init_set_si(x, 2);
    mpz_init_set_si(z, 2);
    mpz_init_set_ui(P, 1);
    
    std::size_t k = 1;
    std::size_t l = 1;
    std::size_t count = 0;
    
    while (mpz_cmp_ui(n, 1) != 0) {
        for (;;) {
            do {
                mpz_mul(t, x, x);
                mpz_mod(x, t, n);
                mpz_add_ui(x, x, a);
                
                mpz_sub(t, z, x);
                mpz_mul(t2, P, t);
                mpz_mod(P, t2, n);
                
                if (k % 32 == 1) {
                    mpz_gcd(t, P, n);
                    
                    if (mpz_cmp_ui(t, 1) != 0)
                        goto factor_found;
                    
                    mpz_set(y, x);
                }
                
                ++count;
            } while (--k != 0 && count < myLimit);
            
            if (count == myLimit)
                goto myReturn;
            
            mpz_set(z, x);
            k = l;
            l = 2 * l;
            
            for (std::size_t i = 0; i < k; ++i) {
                mpz_mul(t, x, x);
                mpz_mod(x, t, n);
                mpz_add_ui(x, x, a);
            }
            
            mpz_set(y, x);
        }
        
    factor_found:
        do {
            mpz_mul(t, y, y);
            mpz_mod(y, t, n);
            mpz_add_ui(y, y, a);
            
            mpz_sub(t, z, y);
            mpz_gcd(t, t, n);
        } while (mpz_cmp_ui(t, 1) == 0);
        
        mpz_divexact(n, n, t);	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p(t, MR_REPS) == 0) {
            returnVal = PollardRhoWithConstraint(t, a + 1, factors, numPs,
                                                  myLens, myLimit, powMultiplier, 
                                                  arrayMax, xtraRecFacs);
            if (returnVal == 1) {
                int ind = numPs - 1;
                
                while (mpz_probab_prime_p(factors[ind], MR_REPS) == 0)
                    --ind;
                
                xtraRecFacs.push_back(ind);
                mpz_mul(factors[ind], factors[ind], t);
                goto myReturn;
            }
        } else {
            mpz_set(factors[numPs], t);
            myLens.push_back(powMultiplier);
            
            while (mpz_divisible_p(n, t)) {
                mpz_divexact(n, n, t);
                myLens[numPs] += powMultiplier;
            }
            
            ++numPs;
            
            if (numPs == arrayMax) {
                returnVal = 1;
                goto myReturn;
            }
        }
        
        if (mpz_probab_prime_p(n, MR_REPS) != 0) {
            mpz_set(factors[numPs], n);
            mpz_set_ui(n, 1);
            myLens.push_back(powMultiplier);
            ++numPs;
            
            if (numPs == arrayMax)
                returnVal = 1;
            
            break;
        }

        mpz_mod(x, x, n);
        mpz_mod(z, z, n);
        mpz_mod(y, y, n);
    }
    
myReturn:
    mpz_clear(P);
    mpz_clear(t2);
    mpz_clear(t);
    mpz_clear(z);
    mpz_clear(x);
    mpz_clear(y);
    return returnVal;
}

void GetBigPrimeFacs(mpz_t n, mpz_t *const factors,
                     mpz_t *const result, std::size_t& numPs,
                     std::vector<std::size_t>& myLens, 
                     std::size_t nThreads, bool bShowStats,
                     std::size_t powMaster, std::size_t arrayMax,
                     std::vector<std::size_t>& extraRecursionFacs) {
    
    if (mpz_sizeinbase(n, 10) < 24) {
        PollardRhoWithConstraint(n, 1, factors, numPs, myLens, 10000000,
                                 powMaster, arrayMax, extraRecursionFacs);
    } else {
        QuadraticSieve(n, result, nThreads, bShowStats);
        
        for (std::size_t i = 0; i < 2; ++i) {
            const std::size_t myPow = ((mpz_perfect_power_p(result[i])) ?
                                        GetPower(result[i]) : 1) * powMaster;
            
            if (mpz_probab_prime_p(result[i], MR_REPS) == 0) {
                mpz_t recurseTemp[2];
                mpz_init(recurseTemp[0]); mpz_init(recurseTemp[1]);
                
                GetBigPrimeFacs(result[i], factors, recurseTemp,
                                numPs, myLens, nThreads, bShowStats,
                                myPow, arrayMax, extraRecursionFacs);
                
                mpz_clear(recurseTemp[0]);
                mpz_clear(recurseTemp[1]);
            } else {
                mpz_divexact(n, n, result[i]);
                mpz_set(factors[numPs], result[i]);
                myLens.push_back(myPow);

                while (mpz_divisible_p (n, result[i]))
                    mpz_divexact (n, n, result[i]);

                ++numPs;

                // This should not happen
                if (numPs == arrayMax)
                    Rcpp::stop("Too many prime factors!!");
            }
        }
    }
}

void QuadSieveHelper(mpz_t nmpz, std::unique_ptr<mpz_t[]> &factors, std::size_t &arrayMax,
                     std::size_t &numUni, std::vector<std::size_t> &lengths,
                     std::size_t nThreads, bool bShowStats) {
    
    std::vector<std::size_t> extraRecursionFacs;
    
    // First we test for small factors.
    int increaseSize = TrialDivision(nmpz, factors.get(), numUni, lengths, arrayMax);
    
    while (increaseSize) {
        const std::size_t oldMax = arrayMax;
        arrayMax <<= 1;
        Grow(factors, oldMax, arrayMax);
        increaseSize = TrialDivision(nmpz, factors.get(), numUni, lengths, arrayMax);
    }
    
    if (mpz_cmp_ui(nmpz, 1) != 0) {
        // We now test for larger primes using pollard's rho
        // algorithm, but constrain it to a limited number of checks
        increaseSize = PollardRhoWithConstraint(nmpz, 1, factors.get(), numUni,
                                                lengths, POLLARD_RHO_REPS,
                                                1, arrayMax, extraRecursionFacs);
        while (increaseSize) {
            const std::size_t oldMax = arrayMax;
            arrayMax <<= 1;
            Grow(factors, oldMax, arrayMax);
            increaseSize = PollardRhoWithConstraint(nmpz, 1, factors.get(), numUni,
                                                    lengths, POLLARD_RHO_REPS,
                                                    1, arrayMax, extraRecursionFacs);
        }

        // extraRecursionFacs are factors that are a result of the pollarRho algo
        // terminating early because of the limitations on the size of the factors
        // array. As a result, we are left with an extra factor, say f, that can't
        // be placed anywhere in the current array. To remedy this, we find an index
        // of the factors array that contains a prime factor, say pj, and set this
        // index to the product (pj * f). We also take note of the index by pushing
        // it to the extraRecursionFacs vector. Below, we take this info and fully
        // factorize these partially factored numbers and add them to our final array.
        const std::size_t eRFSize = extraRecursionFacs.size();

        if (eRFSize > 0) {
            const std::size_t oldMax = arrayMax;
            arrayMax += (eRFSize * mpzChunkBig);
            Grow(factors, oldMax, arrayMax);
            
            auto tempFacs = FromCpp14::make_unique<mpz_t[]>(mpzChunkBig);

            for (std::size_t i = 0; i < mpzChunkBig; ++i)
                mpz_init(tempFacs[i]);

            mpz_t tempNum;
            mpz_init(tempNum);

            for (const auto xtraFacs: extraRecursionFacs) {
                mpz_set(tempNum, factors[xtraFacs]);
                std::size_t tempUni = 0;
                std::vector<std::size_t> tempLens, dummyVec;

                int temp = PollardRhoWithConstraint(tempNum, 1, tempFacs.get(), tempUni,
                                                    tempLens, 10000000, 1,
                                                    100000, dummyVec);
                if (temp == 0) {
                    mpz_set(factors[xtraFacs], tempFacs[0]);

                    for (std::size_t i = 1; i < tempUni; ++i) {
                        mpz_set(factors[numUni], tempFacs[i]);
                        lengths.push_back(tempLens[i]);
                        ++numUni;
                    }
                } else {
                    Rcpp::stop("Too many prime factors!!");
                }
            }

            for (std::size_t i = 0; i < mpzChunkBig; ++i)
                mpz_clear(tempFacs[i]);

            tempFacs.reset();
            mpz_clear(tempNum);
        }

        // If there is less than 60% of mpzChunkBig, then increase array
        // size to ensure that the functions below have enough space
        if ((100.0 * static_cast<double>(arrayMax - numUni) /
                            static_cast<double>(mpzChunkBig)) < 60.0) {
            std::size_t oldMax = arrayMax;
            arrayMax <<= 1;
            Grow(factors, oldMax, arrayMax);
        }
        
        if (mpz_cmp_ui(nmpz, 1) != 0) {
            auto result = FromCpp14::make_unique<mpz_t[]>(2);
            mpz_init(result[0]);
            mpz_init(result[1]);
            
            // Shield quadratic sieve from perfect powers
            const std::size_t myPow = (mpz_perfect_power_p(nmpz)) ? GetPower(nmpz) : 1;
            
            if (mpz_probab_prime_p(nmpz, MR_REPS) != 0) {
                mpz_set(factors[numUni], nmpz);
                lengths.push_back(myPow);
                ++numUni;
            } else {
                GetBigPrimeFacs(nmpz, factors.get(), result.get(), numUni, lengths,
                                nThreads, bShowStats, myPow, arrayMax, extraRecursionFacs);
                mpz_clear(result[0]);
                mpz_clear(result[1]);
            }
            
            result.reset();
        }
    }
}
