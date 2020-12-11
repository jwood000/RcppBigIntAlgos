// Factoring with Pollard's rho method.
// 
// N.B. Below is the original message of the file. There have been minor
// style edits and slight modifications for the RcppBigIntAlgos library
// 
// Copyright 1995, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2009, 2012
// Free Software Foundation, Inc.
// 
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see http://www.gnu.org/licenses/. 

#include "PollardRho.h"

void TrialDivision(mpz_class &t, std::vector<mpz_class> &factors,
                   std::vector<std::size_t>& myLens) {
    
    unsigned long int p = mpz_scan1(t.get_mpz_t(), 0);
    mpz_div_2exp(t.get_mpz_t(), t.get_mpz_t(), p);
    
    if (p) {
        factors.push_back(2);
        myLens.push_back(p);
    }
    
    p = 3;
    
    for (std::size_t i = 1; i < primesDiffPR.size();) {
        if (!mpz_divisible_ui_p(t.get_mpz_t(), p)) {
            p += primesDiffPR[i++];
            
            if (cmp(t, p * p) < 0)
                break;
        } else {
            t /= p;
            factors.push_back(p);
            myLens.push_back(1);
            
            while (mpz_divisible_ui_p(t.get_mpz_t(), p)) {
                t /= p;
                ++myLens.back();
            }
            
            p += primesDiffPR[i++];
            
            if (cmp(t, p * p) < 0)
                break;
        }
    }
}

void PollardRho(mpz_class &n, unsigned long int a,
                std::vector<mpz_class> &factors,
                std::vector<std::size_t>& myLens) {
    
    mpz_class x, z, y, p, t;
    x = y = z = 2;
    p = 1;
    
    std::size_t k = 1u;
    std::size_t q = 1u;

    while (n != 1) {
        for (;;) {
            do {
                x = (x * x) % n + a;
                t = z - x;
                mpz_mod(t.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
                p *= t;
                p %= n;
                
                if (k % 32 == 1) {
                    t = gcd(p, n);
                    
                    if (cmp(t, 1) != 0)
                        goto factor_found;
                    
                    y = x;
                }
	        } while (--k != 0);

            z = x;
            k = q;
            q <<= 1;
            
            for (std::size_t i = 0; i < k; ++i)
                x = (x * x) % n + a;
            
            y = x;
        }

        factor_found:
        do {
            y = (y * y) % n + a;
            t = gcd(z - y, n);
        } while (t == 1);

        n /= t;	/* divide by t, before t is overwritten */

        if (mpz_probab_prime_p(t.get_mpz_t(), MR_REPS) == 0) {
            PollardRho(t, a + 1, factors, myLens);
        } else {
            factors.push_back(t);
            myLens.push_back(1);

            while (mpz_divisible_p(n.get_mpz_t(), t.get_mpz_t())) {
                n /= t;
                ++myLens.back();
            }
        }
        
        if (mpz_probab_prime_p(n.get_mpz_t(), MR_REPS) != 0) {
            factors.push_back(n);
            myLens.push_back(1);
            break;
	    }

        x %= n;
        z %= n;
        y %= n;
    }
}

void GetPrimeFactors(mpz_class &t, std::vector<mpz_class> &factors,
                     std::vector<std::size_t> &myLens) {
    
    if (sgn(t) != 0) {
        TrialDivision(t, factors, myLens);
        
        if (cmp(t, 1) != 0) {
    	    if (mpz_probab_prime_p(t.get_mpz_t(), MR_REPS) != 0) {
    	        factors.push_back(t);
    	        myLens.push_back(1);
    	    } else{
    	        PollardRho(t, 1, factors, myLens);
            }
        }
    }
}
