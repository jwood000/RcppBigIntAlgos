#ifndef POLLARDRHO_GMP_R
#define POLLARDRHO_GMP_R 1
#include "bigvec.h"

/**
 * Get Prime Factorization
 * t: number to factorize
 * factors [out]: the list of factors
 *
 * Note: this is adapted from demo "factorize.c" file from gmplib
 */
void getPrimeFactors (mpz_t t, bigvec & factors);

#endif
