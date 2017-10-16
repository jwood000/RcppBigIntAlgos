/*! 
 *  \file factorization.cc
 *  \brief C functions used for integer factorization
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:33:33 EDT jwood000>
 *
 *  \author Joseph Wood. The first 20 lines (although slightly modified)
 *       are attributed to Antoine Lucas and help from Immanuel Scholz.
 *       See factor.cc from the R gmp package for more details.
 *          Original C code from libgmp.
 *
 *  \note Licence: GPL (>=) 2
 */

#include "Rgmp.h"
#include "pollardrho.h"
#include "factorization.h"

SEXP getDivisorsC (SEXP n) {
    bigvec v = bigintegerR::create_bignum(n);
    bigvec primeFacs, myFacs;
    if(v.size() > 0) {
        if (mpz_cmp_ui(v[0].value.getValueTemp(), 1) == 0) {
            myFacs.push_back(1);
        } else {
            mpz_t val;
            mpz_init(val);
            mpz_t_sentry val_s(val);
            mpz_set(val,v[0].value.getValueTemp());
            
            int sgn = mpz_sgn(val);
            if(sgn == 0)
                error(_("Cannot factorize 0"));
            if (sgn<0) {
                mpz_abs(val,val);
                primeFacs.value.push_back(biginteger(-1));
            }
            getPrimeFactors (val,primeFacs);
            
            std::vector<int> lengths;
            std::vector<biginteger>::iterator it;
            biginteger prev = primeFacs.value[0];
            
            unsigned long int i, j, k, n = primeFacs.size(), numUni = 0;
            mpz_t bigFacs[n];
            lengths.reserve(n);
            for (i = 0; i < n; i++) {mpz_init(bigFacs[i]);}
            
            mpz_set(bigFacs[0], primeFacs.value[0].getValue());
            lengths.push_back(1);
            k = 1;
            
            for(it = primeFacs.value.begin() + 1; it < primeFacs.value.end(); it++) {
                if (prev == *it) {
                    lengths[numUni]++;
                } else {
                    numUni++;
                    prev = *it;
                    lengths.push_back(1);
                    mpz_set(bigFacs[numUni], primeFacs[k].value.getValue());
                }
                k++;
            }
            
            unsigned long int ind, facSize = 1, numFacs = 1;
            for (i = 0; i <= numUni; i++) {numFacs *= (lengths[i]+1);}
            
            myFacs.value.reserve(numFacs);
            mpz_t temp;
            mpz_init(temp);
            
            for (i = 0; i <= lengths[0]; ++i) {
                mpz_pow_ui(temp, bigFacs[0], i);
                myFacs.push_back(temp);
            }
            
            if (numUni > 0) {
                for (j = 1; j <= numUni; j++) {
                    facSize *= (lengths[j-1] + 1);
                    for (i = 1; i <= lengths[j]; i++) {
                        ind = i*facSize;
                        for (k = 0; k < facSize; k++) {
                            mpz_pow_ui(temp, bigFacs[j], i);
                            mpz_mul(temp, temp, myFacs[k].value.getValue());
                            myFacs.push_back(temp);
                        }
                    }
                }
            }
            
            std::sort(myFacs.value.begin(), myFacs.value.end());
        }
    }
    return bigintegerR::create_SEXP(myFacs);
}