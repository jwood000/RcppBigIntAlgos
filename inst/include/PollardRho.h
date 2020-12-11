#ifndef POLLARDRHO_GMP_R
#define POLLARDRHO_GMP_R

#include <Rcpp.h>
#include "GmpxxCopy.h"
#include "Primes.h"

void GetPrimeFactors(mpz_class &t, std::vector<mpz_class> &factors,
                     std::vector<std::size_t> &myLens);

void TrialDivision(mpz_class &t, std::vector<mpz_class> &factors,
                   std::vector<std::size_t> &myLens);

#endif
