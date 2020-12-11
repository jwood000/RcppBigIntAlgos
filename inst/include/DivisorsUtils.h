#ifndef DIVISORS_UTILS_H
#define DIVISORS_UTILS_H

#include "PollardRho.h"
#include "ImportExportMPZ.h"

SEXP FactorNum(mpz_class &val, std::size_t nThreads,
               bool bShowStats, bool bSkipPR, bool bSkipECM);

#endif
