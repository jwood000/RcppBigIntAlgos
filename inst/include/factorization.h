/*! \file factorization.h
 *  \brief header for factorization functions set
 *
 *  \version 1
 *
 *  \date Created: 2005 (antoine)
 *  \date Last modified: Time-stamp: <2017-10-06 11:45:17 EDT jwood000>
 *
 *
 *  \note Licence: GPL
 */

#ifndef GMP_R_FACTORIZATION_HEADER_
#define GMP_R_FACTORIZATION_HEADER_ 1

#include "bigintegerR.h"

extern "C" {
    SEXP getDivisorsC (SEXP n);
}

/** \brief Function used to test factorization with small numbers
 */
void factor_using_division (mpz_t t, int numPrimes,  bigvec & result) ;

/** \brief Function used for factorization
 */
void factor_using_division_2kp (mpz_t t, unsigned int limit, unsigned long p,  bigvec & result) ;

/** \brief Pollard Rho method for factorization
 */
void factor_using_pollard_rho (mpz_t n, int a_int, unsigned long p, bigvec & result);

/** \brief Function that call an algorithm for factorization
 */
void getPrimeFactors (mpz_t t, bigvec & result);

#endif
