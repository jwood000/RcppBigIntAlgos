/* Factoring with Pollard's rho method.

Copyright 1995, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2009, 2012
Free Software Foundation, Inc.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see http://www.gnu.org/licenses/.  */

// This is needed as cinttypes is C++11
#include <inttypes.h>
#include <math.h>
#include "PollardRho.h"

static unsigned char primes_diff[] = {
    #define P(a,b,c) a,
    #include "Primes.h"
    #undef P
};

#define PRIMES_PTAB_ENTRIES (sizeof(primes_diff) / sizeof(primes_diff[0]))

/* Number of Miller-Rabin tests to run when not proving primality. */
#define MR_REPS 25

void factor_using_division (mpz_t t, int numPrimes,
                            mpz_t factors[], unsigned long int& numPs,
                            std::vector<unsigned long int>& myLens) {
    mpz_t q;
    unsigned long int p;
    int i;
    
    mpz_init (q);
    p = mpz_scan1 (t, 0);
    mpz_div_2exp (t, t, p);
    
    if (p) {
        mpz_set_ui(factors[numPs], 2);
        myLens.push_back(p);
        numPs++;
    }

    p = 3;
    for (i = 1; i < numPrimes;) {
        if (!mpz_divisible_ui_p(t, p)) {
            p += primes_diff[i++];
            if (mpz_cmp_ui (t, p * p) < 0)
                break;
        } else {
            mpz_tdiv_q_ui (t, t, p);
            mpz_set_ui(factors[numPs], p);
            myLens.push_back(1);
            
            while (mpz_divisible_ui_p(t, p)) {
                mpz_tdiv_q_ui(t, t, p);
                myLens[numPs]++;
            }
            
            numPs++;
            if (numPs == 50)
                Rf_error(_("Too many prime factors. Result will contain over one quadrillion (10^15) factors!!"));
            
            p += primes_diff[i++];
            if (mpz_cmp_ui (t, p * p) < 0)
                break;
        }
    }
    
    mpz_clear (q);
}

void factor_using_pollard_rho (mpz_t n, unsigned long a,
                               mpz_t factors[], unsigned long int& numPs,
                               std::vector<unsigned long int>& myLens) {
    mpz_t x, z, y, P;
    mpz_t t, t2;
    unsigned long  k, l, i;
    
    mpz_init (t);
    mpz_init (t2);
    mpz_init_set_si (y, 2);
    mpz_init_set_si (x, 2);
    mpz_init_set_si (z, 2);
    mpz_init_set_ui (P, 1);
    k = 1;
    l = 1;

    while (mpz_cmp_ui (n, 1) != 0) {
        for (;;) {
            do {
                mpz_mul (t, x, x);
                mpz_mod (x, t, n);
                mpz_add_ui (x, x, a);
                
                mpz_sub (t, z, x);
                mpz_mul (t2, P, t);
                mpz_mod (P, t2, n);
                
                if (k % 32 == 1) {
                    mpz_gcd (t, P, n);
                    if (mpz_cmp_ui (t, 1) != 0)
                        goto factor_found;
                    mpz_set (y, x);
                }
	        }
	        while (--k != 0);

            mpz_set (z, x);
            k = l;
            l = 2 * l;
            for (i = 0; i < k; i++) {
                mpz_mul (t, x, x);
                mpz_mod (x, t, n);
                mpz_add_ui (x, x, a);
            }
            mpz_set (y, x);
        }

        factor_found:
        do {
            mpz_mul (t, y, y);
            mpz_mod (y, t, n);
            mpz_add_ui (y, y, a);
            
            mpz_sub (t, z, y);
            mpz_gcd (t, t, n);
        }
        while (mpz_cmp_ui (t, 1) == 0);

        mpz_divexact (n, n, t);	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p (t, MR_REPS) == 0) {
            factor_using_pollard_rho (t, a + 1, factors, numPs, myLens);
        } else {
            mpz_set(factors[numPs], t);
            myLens.push_back(1);

            while (mpz_divisible_p (n, t)) {
                mpz_divexact (n, n, t);
                myLens[numPs]++;
            }

            numPs++;
            if (numPs == 50)
                Rf_error(_("Too many prime factors. Result will contain over one quadrillion (10^15) factors!!"));
        }

        if (mpz_probab_prime_p (n, MR_REPS) != 0) {
            mpz_set(factors[numPs], n);
            myLens.push_back(1);
            numPs++;
            break;
	    }

        mpz_mod (x, x, n);
        mpz_mod (z, z, n);
        mpz_mod (y, y, n);
    }

    mpz_clear (P);
    mpz_clear (t2);
    mpz_clear (t);
    mpz_clear (z);
    mpz_clear (x);
    mpz_clear (y);
}

void getPrimeFactors (mpz_t t, mpz_t factors[], unsigned long int& numPs,
                      std::vector<unsigned long int>& myLens) {
    
    if (mpz_sgn (t) != 0) {
        factor_using_division (t, PRIMES_PTAB_ENTRIES, factors, numPs, myLens);
        if (mpz_cmp_ui (t, 1) != 0) {
    	    if (mpz_probab_prime_p (t, MR_REPS) != 0) {
                mpz_set(factors[numPs], t);
    	        myLens.push_back(1);
    	        numPs++;
    	    } else{
    	        factor_using_pollard_rho (t, 1, factors, numPs, myLens);
            }
        }
    }
}
