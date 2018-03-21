#ifndef QUADRATICSIEVE_GMP_R
#define QUADRATICSIEVE_GMP_R 1

#include "Rgmp.h"

/**
 * Factor large number into two smaller numbers (possibly prime)
 * t: number to factorize
 * fudge1: number that alters the size of the prime factor base
 * fudge2: number that adjusts the lower bound for the log sum
 * LenB: Size of the sieving interval
 * factors [out]: the list of factors
 *
 * If fudge1, fudge2, or LenB are left blank, they will be determined
 * in the algorithm based off of the literature regarding the QS
 */

void quadraticSieve (mpz_t myNum, double fudge1,
                     double fudge2,
                     unsigned long int LenB, mpz_t factors[]);

#endif
