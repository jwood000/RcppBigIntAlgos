#ifndef FACTOR_UTILS_H
#define FACTOR_UTILS_H

#include "PollardRho.h"
#include "ImportExportMPZ.h"
#include "Cpp14MakeUnique.h"
#include <memory>

SEXP FactorNum(mpz_t val, std::unique_ptr<mpz_t[]> &primeFacs);

#endif
