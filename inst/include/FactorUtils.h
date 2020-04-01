#ifndef FACTOR_UTILS_H
#define FACTOR_UTILS_H

#include "PollardRho.h"
#include "ImportExportMPZ.h"
#include "Cpp14MakeUnique.h"

SEXP factorNum(mpz_t val, mpz_t primeFacs[]);

#endif
