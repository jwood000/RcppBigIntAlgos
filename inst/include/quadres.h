/*! \file quadres.h
 *  \brief header for quadres function
 *
 *  \version 1
 *
 *  \date Created: 2017
 *  \date Last modified: Time-stamp: <2017-10-06 11:45:17 EDT jwood000>
 *
 *  \note Licence: GPL (>=) 2
 */

#ifndef GMP_R_QUADRES_HEADER_
#define GMP_R_QUADRES_HEADER_ 1

#include "bigintegerR.h"

extern "C" {
    SEXP QuadraticResidueContainer (SEXP n, SEXP p);
}

void TonelliShanksC (mpz_t a, mpz_t p, bigvec & quadRes);

#endif
