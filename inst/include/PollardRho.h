#ifndef POLLARDRHO_GMP_R
#define POLLARDRHO_GMP_R

#include <Rcpp.h>
#include <gmp.h>
#include "Primes.h"

/**
 * Get Prime Factorization
 * t: number to factorize
 * factors [out]: the list of factors
 *
 * Note: this is adapted from demo "factorize.c" file from gmplib
 */
void GetPrimeFactors(mpz_t t, mpz_t *const factors, std::size_t &numPs,
                     std::vector<std::size_t> &myLens);

int TrialDivision(mpz_t t, mpz_t *const factors, std::size_t& numPs,
                  std::vector<std::size_t>& myLens, 
                  std::size_t arrayMax);

#endif
