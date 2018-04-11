/*! \file rsafactorize.h
 *  \brief header for rsafactorize functions set
 *
 *  \version 1
 *
 *  \date Created: 2017 (jwood000)
 *  \date Last modified: Time-stamp: <2017-10-06 11:45:17 EDT jwood000>
 *
 *
 *  \note Licence: GPL
 */

#ifndef GMP_R_RSAFACTORIZE_HEADER_
#define GMP_R_RSAFACTORIZE_HEADER_ 1

#include "Rgmp.h"
#include <inttypes.h>

extern "C" {
    SEXP QuadraticSieveContainer (SEXP n);
}

void quadraticSieve (mpz_t a, double e1, double e2, int64_t lb, mpz_t factors[]);

#endif
