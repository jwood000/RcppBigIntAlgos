#ifndef POLLARDRHO_GMP_R
#define POLLARDRHO_GMP_R

#include <Algos.h>

/**
 * Get Prime Factorization
 * t: number to factorize
 * factors [out]: the list of factors
 *
 * Note: this is adapted from demo "factorize.c" file from gmplib
 */
void getPrimeFactors(mpz_t t, mpz_t factors[], unsigned long int& numPs,
                     std::vector<unsigned long int>& myLens);

void factor_using_division(mpz_t t, int numPrimes,
                           mpz_t factors[], unsigned long int& numPs,
                           std::vector<unsigned long int>& myLens);

#endif
