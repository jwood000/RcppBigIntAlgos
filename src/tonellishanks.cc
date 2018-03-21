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
#include <cstdio>
#include <cstring>
// This is needed as cinttypes is C++11
#include <inttypes.h>
#include "tonellishanks.h"

void TonelliShanksC (mpz_t a, mpz_t p, mpz_t quadRes[]) {
    mpz_t P1, s, myAns1, myAns2, temp;
    mpz_t Legendre2, n, b, g, x, Test, big2;

    mpz_init_set(P1, p);
    mpz_init(temp); mpz_init_set_ui(n, 2);
    mpz_init(Legendre2); mpz_sub_ui(P1, P1, 1);
    mpz_init_set(s, P1); mpz_init(myAns1); mpz_init(myAns2);
    mpz_init(x); mpz_init(b); mpz_init(g);
    mpz_init(Test); mpz_init_set_ui(big2, 2);

    unsigned long int j = 0, r, m = 1;

    j = mpz_scan1 (s, 0);
    mpz_div_2exp (s, s, j);

    if (j == 1) {
        mpz_add_ui (temp, p, 1);
        mpz_div_2exp (temp, temp, 2);
        mpz_powm (myAns1, a, temp, p);
        mpz_neg (temp, myAns1);
        mpz_mod (myAns2, temp, p);
    } else {
        mpz_div_2exp (temp, P1, 1);
        mpz_powm (Legendre2, n, temp, p);
        while (mpz_cmp_ui(Legendre2, 1) == 0) {
            mpz_add_ui(n, n, 1);
            mpz_powm (Legendre2, n, temp, p);
        }

        mpz_add_ui(temp, s, 1);
        mpz_div_2exp(temp, temp, 1);
        mpz_powm(x, a, temp, p);
        mpz_powm(b, a, s, p);
        mpz_powm(g, n, s, p);

        r = j;
        m = 1;
        mpz_mod(Test, b, p);

        while ((mpz_cmp_ui(Test, 1) != 0) && (m != 0)) {
            m = 0;
            mpz_mod(Test, b, p);
            while (mpz_cmp_ui(Test, 1) != 0) {
                m++;
                mpz_pow_ui(temp, big2, m);
                mpz_powm(Test, b, temp, p);
            }
            if (m != 0) {
                mpz_pow_ui(temp, big2, r-m-1);
                mpz_powm(temp, g, temp, p);
                mpz_mul(temp, temp, x);
                mpz_mod(x, temp, p);

                mpz_pow_ui(temp, big2, r-m);
                mpz_powm(g, g, temp, p);

                mpz_mul(temp, b, g);
                mpz_mod(b, temp, p);
                r = m;
            }
            mpz_set_ui(Test, 0);
        }
        mpz_set(myAns1, x);
        mpz_sub(temp, p, x);
        mpz_mod(myAns2, temp, p);
    }

    mpz_set(quadRes[0], myAns1);
    mpz_set(quadRes[1], myAns1);
}
