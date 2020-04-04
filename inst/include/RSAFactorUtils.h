#ifndef RSA_FACTOR_UTILS_H
#define RSA_FACTOR_UTILS_H

#include "QuadraticSieve.h"
#include "ImportExportMPZ.h"
#include "Cpp14MakeUnique.h"
#include "PollardRho.h"

void QuadSieveHelper(mpz_t nmpz, mpz_t factors[], std::size_t &arrayMax,
                     std::size_t &numUni, std::vector<std::size_t> &lengths,
                     std::size_t nThreads, bool bShowStats);

#endif
