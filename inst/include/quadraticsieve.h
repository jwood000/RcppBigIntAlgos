#ifndef QUADRATICSIEVE_GMP_R
#define QUADRATICSIEVE_GMP_R 1
#include "bigvec.h"

/**
 * Factorize a prime number
 * t: number to factorize
 * factors [out]: the list of factors
 *
 */

void quadraticSieve (mpz_t myNum, double fudge1,
                     double fudge2,
                     unsigned long int LenB, bigvec & factors);

#endif
