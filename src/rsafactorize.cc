/*! 
 *  \file rsafactorize.cc
 *  \brief C function that transfers input from R to 
 *          quadraticSieve function for factoring large
 *            numbers and returning result to R console
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:33:33 EDT jwood000>
 *
 *  \author Joseph Wood
 *
 *  \note Licence: GPL (>=) 2  
 */

#include "quadraticsieve.h"
#include "rsafactorize.h"
#include "importExportMPZ.h"
#include <vector>
#include "Rgmp.h"

unsigned long int mpzChunkBig = 50;

static unsigned char primes_diff[] = {
    #define P(a,b,c) a,
    #include "primes.h"
    #undef P
};

#define PRIMES_PTAB_ENTRIES (sizeof(primes_diff) / sizeof(primes_diff[0]))

/* Number of Miller-Rabin tests to run when not proving primality. */
#define MR_REPS 25

// Max number of iterations in the main loop
#define POLLARD_RHO_REPS 100000

unsigned long int getPower(mpz_t nmpz) {
    mpz_t testRoot;
    mpz_init(testRoot);
    unsigned long int p = 2;
    unsigned long int myPow = 1;
    
    for (std::size_t i = 0; i < PRIMES_PTAB_ENTRIES; ) {
        mpz_root(testRoot, nmpz, p);
        mpz_pow_ui(testRoot, testRoot, p);
        if (mpz_cmp(testRoot, nmpz) == 0) {
            unsigned long int powPrime = p;
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
        
        p += primes_diff[i++];
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
                unsigned long int powPrime = p;
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
    }
    
    return myPow;
}

int trialDivision (mpz_t t, int numPrimes,
                    mpz_t factors[], unsigned int& numPs,
                    std::vector<unsigned int>& myLens, 
                    unsigned int arrayMax) {
    mpz_t q;
    unsigned long int p;
    int i;
    
    mpz_init (q);
    p = mpz_scan1 (t, 0);
    mpz_div_2exp (t, t, p);
    
    if (p) {
        mpz_set_ui(factors[numPs], 2);
        myLens.push_back(p);
        numPs++;
    }
    
    p = 3;
    for (i = 1; i < numPrimes;) {
        if (!mpz_divisible_ui_p(t, p)) {
            p += primes_diff[i++];
            if (mpz_cmp_ui (t, p * p) < 0)
                break;
        } else {
            mpz_tdiv_q_ui (t, t, p);
            mpz_set_ui(factors[numPs], p);
            myLens.push_back(1);
            
            while (mpz_divisible_ui_p(t, p)) {
                mpz_tdiv_q_ui(t, t, p);
                myLens[numPs]++;
            }
            
            numPs++;
            if (numPs == arrayMax)
                return 1;
            
            p += primes_diff[i++];
            if (mpz_cmp_ui (t, p * p) < 0)
                break;
        }
    }
    
    mpz_clear (q);
    return 0;
}

int pollardRhoWithConstraint (mpz_t n, unsigned long int a,
                               mpz_t factors[], unsigned int& numPs,
                               std::vector<unsigned int>& myLens,
                               unsigned long int myLimit,
                               unsigned long int powMultiplier,
                               unsigned int arrayMax,
                               std::vector<unsigned int>& extraRecursionFacs) {
    mpz_t x, z, y, P;
    mpz_t t, t2;
    unsigned long  k, l, i;
    int returnVal = 0;
    
    mpz_init (t);
    mpz_init (t2);
    mpz_init_set_si (y, 2);
    mpz_init_set_si (x, 2);
    mpz_init_set_si (z, 2);
    mpz_init_set_ui (P, 1);
    k = 1;
    l = 1;
    unsigned long int count = 0;
    
    while (mpz_cmp_ui (n, 1) != 0) {
        for (;;) {
            do {
                mpz_mul (t, x, x);
                mpz_mod (x, t, n);
                mpz_add_ui (x, x, a);
                
                mpz_sub (t, z, x);
                mpz_mul (t2, P, t);
                mpz_mod (P, t2, n);
                
                if (k % 32 == 1) {
                    mpz_gcd (t, P, n);
                    if (mpz_cmp_ui (t, 1) != 0)
                        goto factor_found;
                    mpz_set (y, x);
                }
                count++;
            }
            while (--k != 0 && count < myLimit);
            
            if (count == myLimit)
                goto myReturn;
            
            mpz_set (z, x);
            k = l;
            l = 2 * l;
            for (i = 0; i < k; i++) {
                mpz_mul (t, x, x);
                mpz_mod (x, t, n);
                mpz_add_ui (x, x, a);
            }
            mpz_set (y, x);
        }
        
    factor_found:
        do {
            mpz_mul (t, y, y);
            mpz_mod (y, t, n);
            mpz_add_ui (y, y, a);
            
            mpz_sub (t, z, y);
            mpz_gcd (t, t, n);
        }
        while (mpz_cmp_ui (t, 1) == 0);
        
        mpz_divexact (n, n, t);	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p (t, MR_REPS) == 0) {
            returnVal = pollardRhoWithConstraint (t, a + 1, factors, numPs,
                                                  myLens, myLimit, powMultiplier, 
                                                  arrayMax, extraRecursionFacs);
            if (returnVal == 1) {
                int ind = numPs - 1;
                while (mpz_probab_prime_p (factors[ind], MR_REPS) == 0) {ind--;}
                extraRecursionFacs.push_back(ind);
                mpz_mul(factors[ind], factors[ind], t);
                goto myReturn;
            }
        } else {
            mpz_set(factors[numPs], t);
            myLens.push_back(powMultiplier);
            
            while (mpz_divisible_p (n, t)) {
                mpz_divexact (n, n, t);
                myLens[numPs] += powMultiplier;
            }
            
            numPs++;
            if (numPs == arrayMax) {
                returnVal = 1;
                goto myReturn;
            }
        }
        
        if (mpz_probab_prime_p (n, MR_REPS) != 0) {
            mpz_set(factors[numPs], n);
            mpz_set_ui(n, 1);
            myLens.push_back(powMultiplier);
            numPs++;
            if (numPs == arrayMax)
                returnVal = 1;
            break;
        }

        mpz_mod (x, x, n);
        mpz_mod (z, z, n);
        mpz_mod (y, y, n);
    }
    
myReturn:
    mpz_clear (P);
    mpz_clear (t2);
    mpz_clear (t);
    mpz_clear (z);
    mpz_clear (x);
    mpz_clear (y);
    return returnVal;
}

void getBigPrimeFacs(mpz_t n, mpz_t factors[],
                    mpz_t result[], unsigned int& numPs,
                    std::vector<unsigned int>& myLens,
                    unsigned long int powMaster,
                    unsigned int arrayMax,
                    std::vector<unsigned int>& extraRecursionFacs) {
    
    if (mpz_sizeinbase(n, 10) < 24) {
        pollardRhoWithConstraint(n, 1, factors, numPs, myLens, 10000000,
                                powMaster, arrayMax, extraRecursionFacs);
    } else {
        quadraticSieve (n, 0.0, 0.0, 0, result);
        
        for (std::size_t i = 0; i < 2; i++) {
            unsigned long int myPow = 1;
            if (mpz_perfect_power_p(result[i]))
                myPow = getPower(result[i]);

            myPow *= powMaster;
            if (mpz_probab_prime_p(result[i], MR_REPS) == 0) {
                mpz_t recurseTemp[2];
                mpz_init(recurseTemp[0]); mpz_init(recurseTemp[1]);
                getBigPrimeFacs(result[i], factors, recurseTemp,
                                numPs, myLens, myPow, arrayMax, extraRecursionFacs);
                mpz_clear(recurseTemp[0]); mpz_clear(recurseTemp[1]);
            } else {
                mpz_divexact(n, n, result[i]);
                mpz_set(factors[numPs], result[i]);
                myLens.push_back(myPow);

                while (mpz_divisible_p (n, result[i]))
                    mpz_divexact (n, n, result[i]);

                numPs++;
                
                // This should not happen
                if (numPs == arrayMax)
                    error("Too many prime factors!!");
            }
        }
    }
}

SEXP QuadraticSieveContainer (SEXP Rn) {
    
    unsigned int vSize;
    
    switch (TYPEOF(Rn)) {
        case RAWSXP: {
            const char* raw = (char*)RAW(Rn);
            vSize = ((int*)raw)[0];
            break;
        }
        default:
            vSize = LENGTH(Rn);
    }
    
    mpz_t myVec[1];
    mpz_init(myVec[0]);
    
    if (vSize > 1)
        error(_("Can only factor one number at a time"));
    
    // This is from the importExportMPZ header
    createMPZArray(Rn, myVec, 1);
    mpz_t nmpz;
    mpz_init_set(nmpz, myVec[0]);
    mpz_t *result;
    result = (mpz_t *) malloc(2 * sizeof(mpz_t));
    mpz_init(result[0]); mpz_init(result[1]);

    std::vector<unsigned int> lengths, extraRecursionFacs;
    unsigned int arrayMax = mpzChunkBig, numUni = 0;
    unsigned long int myPow = 1;
    
    mpz_t *factors;
    factors = (mpz_t *) malloc(mpzChunkBig * sizeof(mpz_t));
    for (std::size_t i = 0; i < mpzChunkBig; i++)
        mpz_init(factors[i]);

    // First we test for small factors.
    int increaseSize = trialDivision(nmpz, PRIMES_PTAB_ENTRIES,
                            factors, numUni, lengths, arrayMax);
    
    while (increaseSize) {
        arrayMax += mpzChunkBig;
        factors = (mpz_t *) realloc(factors, arrayMax * sizeof factors[0]);
        for (std::size_t i = (arrayMax - mpzChunkBig); i < arrayMax; i++)
            mpz_init(factors[i]);
        
        increaseSize = trialDivision(nmpz, PRIMES_PTAB_ENTRIES,
                                     factors, numUni, lengths, arrayMax);
    }
    
    if (mpz_cmp_ui(nmpz, 1) != 0) {
        // We now test for larger primes using pollard's rho
        // algorithm, but constrain it to a limited number of checks
        increaseSize = pollardRhoWithConstraint(nmpz, 1, factors, numUni,
                                                lengths, POLLARD_RHO_REPS,
                                                1, arrayMax, extraRecursionFacs);
        while (increaseSize) {
            arrayMax += mpzChunkBig;
            factors = (mpz_t *) realloc(factors, arrayMax * sizeof factors[0]);
            for (std::size_t i = (arrayMax - mpzChunkBig); i < arrayMax; i++)
                mpz_init(factors[i]);
            
            increaseSize = pollardRhoWithConstraint(nmpz, 1, factors, numUni,
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
        unsigned int eRFSize = extraRecursionFacs.size();
        
        if (eRFSize > 0) {
            arrayMax += (eRFSize * mpzChunkBig);
            factors = (mpz_t *) realloc(factors, arrayMax * sizeof factors[0]);
            for (std::size_t i = (arrayMax - (eRFSize * mpzChunkBig)); i < arrayMax; i++)
                mpz_init(factors[i]);

            mpz_t tempNum;
            mpz_t *tempFacs;
            tempFacs = (mpz_t *) malloc(mpzChunkBig * sizeof(mpz_t));
            for (std::size_t i = 0; i < mpzChunkBig; i++)
                mpz_init(tempFacs[i]);
            
            mpz_init(tempNum);
            std::vector<unsigned int>::iterator it, itEnd = extraRecursionFacs.end();
            for (it = extraRecursionFacs.begin(); it < itEnd; it++) {
                mpz_set(tempNum, factors[*it]);
                unsigned int tempUni = 0;
                std::vector<unsigned int> tempLens, dummyVec;
                int temp = pollardRhoWithConstraint(tempNum, 1, tempFacs, tempUni,
                                                    tempLens, 10000000, 1,
                                                    100000, dummyVec);
                if (temp == 0) {
                    mpz_set(factors[*it], tempFacs[0]);
                    for (std::size_t i = 1; i < tempUni; i++) {
                        mpz_set(factors[numUni], tempFacs[i]);
                        lengths.push_back(tempLens[i]);
                        numUni++;
                    }
                } else {
                    error("Too many prime factors!!");
                }
            }
        }
        
        // If there is less than 60% of mpzChunkBig, then increase array
        // size to ensure that the functions below have enough space
        if ((100 * (arrayMax - numUni) / mpzChunkBig) < 60) {
            arrayMax += mpzChunkBig;
            factors = (mpz_t *) realloc(factors, arrayMax * sizeof factors[0]);
            for (std::size_t i = arrayMax - mpzChunkBig; i < arrayMax; i++)
                mpz_init(factors[i]);
        }
        
        if (mpz_cmp_ui(nmpz, 1) != 0) {
            // Protect quadratic sieve from perfect powers
            if (mpz_perfect_power_p(nmpz))
                myPow = getPower(nmpz);

            if (mpz_probab_prime_p (nmpz, MR_REPS) != 0) {
                mpz_set(factors[numUni], nmpz);
                lengths.push_back(1);
                numUni++;
            } else {
                getBigPrimeFacs(nmpz, factors, result, numUni,
                                lengths, myPow, arrayMax, extraRecursionFacs);
            }
        }
    }
    
    // Sort the prime factors as well as order the
    // lengths vector by the order of the factors array
    quickSort(factors, 0, numUni - 1, lengths);

    unsigned long int totalNum = 0;
    for (std::size_t i = 0; i < numUni; i++)
            totalNum += lengths[i];

    unsigned long int intSize = sizeof(int); // starting with vector-size-header
    unsigned long int numb = 8 * intSize;
    unsigned long int tempSize, size = intSize;
    std::vector<unsigned long int> mySizes(totalNum);
    unsigned long int count = 0;

    for (std::size_t i = 0; i < numUni; i++) { // adding each bigint's needed size
        for (std::size_t j = 0; j < lengths[i]; j++, count++) {
            tempSize = intSize * (2 + (mpz_sizeinbase(factors[i], 2) + numb - 1) / numb);
            size += tempSize;
            mySizes[count] = tempSize;
        }
    }

    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size));
    char* r = (char*)(RAW(ans));
    ((int*)(r))[0] = totalNum; // first int is vector-size-header

    // current position in pos[] (starting after vector-size-header)
    unsigned long int pos = intSize;
    count = 0;

    for (std::size_t i = 0; i < numUni; i++)
        for (std::size_t j = 0; j < lengths[i]; j++, count++)
            pos += myRaw(&r[pos], factors[i], mySizes[count]);

    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
    UNPROTECT(1);
    return(ans);
}