/*! 
 *  \file tonelliShanks.cc
 *  \brief C function for obtaining quadratic residues via Tonelli Shanks algorithm
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *
 *  \author Joseph Wood (C adaptation of algorithm obtained from the following 
 *      journal article: Brown, Ezra B. “Square Roots from 1; 24, 51, 10 to Dan Shanks.” 
 *                          The College Mathematics Journal, vol. 30, 1999, pp. 82–95.)
 *      URL: https://www.math.vt.edu/people/brown/doc/sqrts.pdf
 *
 *  \note Licence: GPL (>=) 2
*/

#include <cstdlib>
#include <cstdint>
#include <gmp.h>

void TonelliShanksC(mpz_t myNum, mpz_t p, mpz_t *const TS) {
    
    mpz_sub_ui(TS[0], p, 1);
    mpz_set_ui(TS[3], 2);
    mpz_set_ui(TS[8], 1u);
    
    const int pow2 = mpz_scan1(TS[0], 0);
    
    if (pow2 == 1) {
        mpz_set(TS[2], p);
        mpz_add_ui(TS[2], TS[2], 1);
        mpz_div_2exp(TS[2], TS[2], 2);
        mpz_powm(TS[1], myNum, TS[2], p);
        mpz_neg(TS[2], TS[1]);
        mpz_mod(TS[2], TS[2], p);
    } else {
        mpz_div_2exp(TS[9], TS[0], 1);
        mpz_powm(TS[2], TS[3], TS[9], p);
        
        while (mpz_cmp_ui(TS[2], 1) == 0) {
            mpz_add_ui(TS[3], TS[3], 1);
            mpz_powm(TS[2], TS[3], TS[9], p);
        }
        
        
        mpz_div_2exp(TS[9], TS[0], pow2);
        mpz_add_ui(TS[0], TS[9], 1);
        mpz_div_2exp(TS[0], TS[0], 1);
        
        mpz_powm(TS[4], myNum, TS[9], p);
        mpz_powm(TS[5], TS[3], TS[9], p);
        mpz_powm(TS[6], myNum, TS[0], p);
        mpz_mod(TS[7], TS[4], p);
        
        int r = pow2;
        int m = mpz_cmp_ui(TS[7], 1);
        
        while (m) {
            m = 0;
            mpz_mod(TS[7], TS[4], p);
            
            while (mpz_cmp_ui(TS[7], 1)) {
                ++m;
                mpz_mul_2exp(TS[1], TS[8], m);
                mpz_powm(TS[7], TS[4], TS[1], p);
            }
            
            if (m) {
                mpz_mul_2exp(TS[1], TS[8], r - m - 1);
                mpz_powm(TS[1], TS[5], TS[1], p);
                mpz_mul(TS[1], TS[1], TS[6]);
                mpz_mod(TS[6], TS[1], p);
                
                mpz_mul_2exp(TS[1], TS[8], r - m);
                mpz_powm(TS[5], TS[5], TS[1], p);
                
                mpz_mul(TS[1], TS[4], TS[5]);
                mpz_mod(TS[4], TS[1], p);
                r = m;
            }
        }
        
        mpz_set(TS[1], TS[6]);
        mpz_sub(TS[2], p, TS[1]);
    }
}
