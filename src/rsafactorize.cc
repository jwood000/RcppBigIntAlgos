/*! 
 *  \file rsafactorize.cc
 *  \brief C function that transfers input from R to 
 *          quadraticSieve function for factoring large
 *            numbers and returning result to R console
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:33:33 EDT jwood000>
 *
 *  \author Joseph Wood
 *
 *  \note Licence: GPL (>=) 2  
 */

#include "Rgmp.h"
#include "quadraticsieve.h"
#include "rsafactorize.h"

typedef std::vector<signed long int> v1d;
typedef std::vector<v1d> v2d;

SEXP QuadraticSieveContainer (SEXP n) {
    bigvec nBig = bigintegerR::create_bignum(n);
    mpz_t nmpz;
    mpz_init_set(nmpz, nBig[0].value.getValue());
    bigvec result;
    
    // factor_using_division(nmpz, result);
    
    quadraticSieve (nmpz, 0.0, 0.0, 0, result);
    return bigintegerR::create_SEXP(result);
}