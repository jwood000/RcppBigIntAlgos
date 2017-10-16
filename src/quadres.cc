/*! 
 *  \file quadRes.cc
 *  \brief C function that transfers input from R to 
 *          TonelliShanksC function for obtaining quadratic residues
 *            and returning result to R console
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:08:27 EDT jwood000>
 *  
 *  \author Joseph Wood
 *  
 */

#include "Rgmp.h"
#include "tonellishanks.h"
#include "quadres.h"

SEXP QuadraticResidueContainer (SEXP n, SEXP p) {
    mpz_t nmpz, pmpz, test;
    bigvec v1 = bigintegerR::create_bignum(n);
    bigvec v2 = bigintegerR::create_bignum(p);
    mpz_init(nmpz);
    mpz_init(pmpz);
    mpz_init(test);
    mpz_set(nmpz,v1[0].value.getValueTemp());
    mpz_set(pmpz,v2[0].value.getValueTemp());
    bigvec result;
    
    // Test is Legendre symbol is 1. If not exit,
    // otherwise get quadratic residues. First,
    // p must be prime.
    unsigned long int probTest = mpz_probab_prime_p(pmpz, 40);
    if (probTest == 0)
        error(_("p must be prime!!"));
    
    mpz_sub_ui(test, pmpz, 1);
    mpz_div_2exp(test, test, 1);
    mpz_powm(test, nmpz, test, pmpz);
    if (mpz_cmp_ui(test, 1) != 0)
        error(_("n must be a quadratic residue modulo p!!"));
    
    TonelliShanksC(nmpz, pmpz, result);
    
    return bigintegerR::create_SEXP(result);
}
