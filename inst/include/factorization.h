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

/** \brief Function used to test factorization with small numbers
 */
void factor_using_division(mpz_t t, int numPrimes,
                           mpz_t factors[], unsigned int& numPs,
                           std::vector<unsigned int>& myLens) ;

/** \brief Pollard Rho method for factorization
 */
void factor_using_pollard_rho(mpz_t n, unsigned long a,
                              mpz_t factors[], unsigned int& numPs,
                              std::vector<unsigned int>& myLens);

/** \brief Function that call an algorithm for factorization
 */
void getPrimeFactors(mpz_t t, mpz_t factors[], unsigned int& numPs,
                     std::vector<unsigned int>& myLens);

#endif
