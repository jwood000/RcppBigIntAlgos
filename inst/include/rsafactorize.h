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

#ifndef GMP_R_RSAFACTORIZE_H
#define GMP_R_RSAFACTORIZE_H

#include <Algos.h>

extern "C" {
    SEXP QuadraticSieveContainer(SEXP n);
}

void QuadraticSieve(mpz_t a, mpz_t *const factors);

#endif
