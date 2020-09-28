#ifndef RSA_FACTOR_UTILS_H
#define RSA_FACTOR_UTILS_H

#include "QuadraticSieve.h"
#include "ImportExportMPZ.h"
#include "PollardRho.h"

void QuadSieveHelper(mpz_class &nMpz, std::vector<mpz_class> &factors,
                     std::vector<std::size_t> &lengths, std::size_t nThreads,
                     bool bShowStats, bool bSkipExtPR);

#endif
