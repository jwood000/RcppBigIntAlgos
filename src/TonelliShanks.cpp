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
#include <gmp.h>

void TonelliShanksC(mpz_t myNum, mpz_t p, mpz_t TS_1) {
    
    mpz_t TS_0;
    mpz_init(TS_0);
    
    mpz_sub_ui(TS_0, p, 1u);
    const int pow2 = mpz_scan1(TS_0, 0);
    
    if (pow2 == 1) {
        mpz_add_ui(TS_1, p, 1);
        mpz_div_2exp(TS_1, TS_1, 2);
        mpz_powm(TS_1, myNum, TS_1, p);
    } else if (pow2 == 2) {
        mpz_t TS_2, TS_3, TS_4;
        mpz_init(TS_2); mpz_init(TS_3); mpz_init(TS_4);
        
        mpz_mod(TS_3, myNum, p);
        mpz_mul_2exp(TS_4, TS_3, 1);
        mpz_sub_ui(TS_2, p, 5u);
        mpz_div_2exp(TS_2, TS_2, 3);
        
        mpz_powm(TS_2, TS_4, TS_2, p);
        mpz_mul(TS_1, TS_3, TS_2);
        mpz_mul(TS_1, TS_1, TS_2);
        mpz_mul_2exp(TS_1, TS_1, 1);
        
        mpz_mod(TS_1, TS_1, p);
        mpz_sub_ui(TS_1, TS_1, 1u);
        mpz_mul(TS_1, TS_1, TS_2);
        mpz_mul(TS_1, TS_1, TS_3);
        mpz_mod(TS_1, TS_1, p);
        
        mpz_clear(TS_2); mpz_clear(TS_3); mpz_clear(TS_4);
    } else {
        mpz_t TS_3, TS_4, TS_5, TS_6, TS_7, TS_8, TS_9;
        
        mpz_init(TS_3); mpz_init(TS_4); mpz_init(TS_5);
        mpz_init(TS_6); mpz_init(TS_7); mpz_init(TS_8);
        mpz_init(TS_9);
        
        mpz_set_ui(TS_3, 2u);
        mpz_set_ui(TS_8, 1u);
        mpz_div_2exp(TS_9, TS_0, 1);
        mpz_powm(TS_1, TS_3, TS_9, p);
        
        while (mpz_cmp_ui(TS_1, 1) == 0) {
            mpz_add_ui(TS_3, TS_3, 1);
            mpz_powm(TS_1, TS_3, TS_9, p);
        }
        
        mpz_div_2exp(TS_9, TS_0, pow2);
        mpz_add_ui(TS_0, TS_9, 1);
        mpz_div_2exp(TS_0, TS_0, 1);
        
        mpz_powm(TS_4, myNum, TS_9, p);
        mpz_powm(TS_5, TS_3, TS_9, p);
        mpz_powm(TS_6, myNum, TS_0, p);
        mpz_mod(TS_7, TS_4, p);
        
        int r = pow2;
        int m = mpz_cmp_ui(TS_7, 1);
        
        while (m) {
            m = 0;
            mpz_mod(TS_7, TS_4, p);
            
            while (mpz_cmp_ui(TS_7, 1u)) {
                ++m;
                mpz_mul_2exp(TS_1, TS_8, m);
                mpz_powm(TS_7, TS_4, TS_1, p);
            }
            
            if (m) {
                mpz_mul_2exp(TS_1, TS_8, r - m - 1);
                mpz_powm(TS_1, TS_5, TS_1, p);
                mpz_mul(TS_1, TS_1, TS_6);
                mpz_mod(TS_6, TS_1, p);
                
                mpz_mul_2exp(TS_1, TS_8, r - m);
                mpz_powm(TS_5, TS_5, TS_1, p);
                
                mpz_mul(TS_1, TS_4, TS_5);
                mpz_mod(TS_4, TS_1, p);
                r = m;
            }
        }
        
        mpz_set(TS_1, TS_6);
        
        mpz_clear(TS_3); mpz_clear(TS_4); mpz_clear(TS_5);
        mpz_clear(TS_6); mpz_clear(TS_7); mpz_clear(TS_8);
        mpz_clear(TS_9);
    }
    
    mpz_clear(TS_0);
}
