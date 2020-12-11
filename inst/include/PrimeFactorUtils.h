#ifndef PRIME_FACTOR_UTILS_H
#define PRIME_FACTOR_UTILS_H

#include "QuadraticSieve.h"
#include "ImportExportMPZ.h"
#include "PollardRho.h"

void QuadSieveHelper(mpz_class &nMpz, std::vector<mpz_class> &factors,
                     std::vector<std::size_t> &lengths, std::size_t nThreads,
                     bool bShowStats, bool bSkipPR, bool bSkipECM);

SEXP PrimeFactorizeHuge(mpz_class &nMpz, std::size_t nThreads,
                        bool bShowStats, bool bSkipPR, bool bSkipECM);

#endif
