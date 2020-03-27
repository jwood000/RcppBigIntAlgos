/*! \file factorization.h
 *  \brief header for factorization functions set
 *
 *  \version 1
 *
 *  \date Created: 2005 (antoine)
 *  \date Last modified: Time-stamp: <2018-03-15 15:13:00 EDT jwood000>
 *
 *
 *  \note Licence: GPL
 */

#ifndef GMP_R_FACTORIZATION_H
#define GMP_R_FACTORIZATION_H

#include <Algos.h>

extern "C" {
    SEXP getDivisorsC(SEXP Rv, SEXP RNamed);
}

/** \brief Function that call an algorithm for factorization
 */
void getPrimeFactors(mpz_t t, mpz_t factors[], std::size_t &numPs,
                     std::vector<std::size_t> &myLens);

#endif
