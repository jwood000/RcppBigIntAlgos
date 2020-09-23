/*! 
 *  \file tonelliShanks.cc
 *  \brief C function for obtaining quadratic residues via Tonelli Shanks algorithm
 *
 *  \version 2
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
#include "GmpxxCopy.h"

void TonelliShanksC(const mpz_class &myNum, const mpz_class &p, mpz_class &TS_1) {
    
    mpz_class TS_0 = p - 1;
    const int pow2 = mpz_scan1(TS_0.get_mpz_t(), 0);
    
    if (pow2 == 1) {
        TS_1 = p + 1;
        TS_1 /= 4;
        mpz_powm(TS_1.get_mpz_t(), myNum.get_mpz_t(),
                 TS_1.get_mpz_t(), p.get_mpz_t());
    } else if (pow2 == 2) {
        mpz_class TS_2, TS_3, TS_4;
        
        TS_3 = myNum % p;
        TS_4 = TS_3 * 2;
        TS_2 = p - 5;
        TS_2 /= 8;
        
        mpz_powm(TS_2.get_mpz_t(), TS_4.get_mpz_t(),
                 TS_2.get_mpz_t(), p.get_mpz_t());
        TS_1 = TS_3 * TS_2;
        TS_1 *= (TS_2 * 2);
        
        TS_1 %= p;
        --TS_1;
        
        TS_1 *= (TS_2 * TS_3);
        TS_1 %= p;
    } else {
        mpz_class TS_3, TS_4, TS_5, TS_6, TS_7, TS_8, TS_9;
        
        TS_3 = 2;
        TS_8 = 1;
        TS_9 = TS_0 / 2;
        mpz_powm(TS_1.get_mpz_t(), TS_3.get_mpz_t(),
                 TS_9.get_mpz_t(), p.get_mpz_t());
        
        while (cmp(TS_1, 1) == 0) {
            ++TS_3;
            mpz_powm(TS_1.get_mpz_t(), TS_3.get_mpz_t(),
                     TS_9.get_mpz_t(), p.get_mpz_t());
        }
        
        mpz_div_2exp(TS_9.get_mpz_t(), TS_0.get_mpz_t(), pow2);
        TS_0 = TS_9 + 1;
        TS_0 /= 2;
        
        mpz_powm(TS_4.get_mpz_t(), myNum.get_mpz_t(),
                 TS_9.get_mpz_t(), p.get_mpz_t());
        mpz_powm(TS_5.get_mpz_t(), TS_3.get_mpz_t(),
                 TS_9.get_mpz_t(), p.get_mpz_t());
        mpz_powm(TS_6.get_mpz_t(), myNum.get_mpz_t(),
                 TS_0.get_mpz_t(), p.get_mpz_t());
        TS_7 = TS_4 % p;
        
        int r = pow2;
        int m = cmp(TS_7, 1);
        
        while (m) {
            m = 0;
            TS_7 = TS_4 % p;
            
            while (cmp(TS_7, 1)) {
                ++m;
                mpz_mul_2exp(TS_1.get_mpz_t(), TS_8.get_mpz_t(), m);
                mpz_powm(TS_7.get_mpz_t(), TS_4.get_mpz_t(),
                         TS_1.get_mpz_t(), p.get_mpz_t());
            }
            
            if (m) {
                mpz_mul_2exp(TS_1.get_mpz_t(), TS_8.get_mpz_t(), r - m - 1);
                mpz_powm(TS_1.get_mpz_t(), TS_5.get_mpz_t(),
                         TS_1.get_mpz_t(), p.get_mpz_t());
                TS_1 *= TS_6;
                TS_6 = TS_1 % p;
                
                mpz_mul_2exp(TS_1.get_mpz_t(), TS_8.get_mpz_t(), r - m);
                mpz_powm(TS_5.get_mpz_t(), TS_5.get_mpz_t(),
                         TS_1.get_mpz_t(), p.get_mpz_t());
                
                TS_1 = TS_4 * TS_5;
                TS_4 = TS_1 % p;
                r = m;
            }
        }
        
        TS_1 = TS_6;
    }
}
