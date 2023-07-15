#ifndef DIVISORS_UTILS_H
#define DIVISORS_UTILS_H

#include "CppConvert/GmpConvert.h"
#include "PollardRho.h"

SEXP FactorNum(mpz_class &val, std::size_t nThreads,
               bool bShowStats, bool bSkipPR, bool bSkipECM);

#endif
