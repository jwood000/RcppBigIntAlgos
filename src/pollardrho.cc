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

#include <cstdlib>
#include <cstdio>
#include <cstring>
// This is needed as cinttypes is C++11
#include <inttypes.h>
#include <math.h>
#include "Rgmp.h"
#include "pollardrho.h"

static unsigned char primes_diff[] = {
    #define P(a,b,c) a,
    #include "primes.h"
    #undef P
};

#define PRIMES_PTAB_ENTRIES (sizeof(primes_diff) / sizeof(primes_diff[0]))

/* Prove primality or run probabilistic tests.  */
int flag_prove_primality = 1;

/* Number of Miller-Rabin tests to run when not proving primality. */
#define MR_REPS 25

void
factor_using_division (mpz_t t, bigvec & factors)
{
  mpz_t q;
  unsigned long int p;
  int i;

  mpz_init (q);

  p = mpz_scan1 (t, 0);
  mpz_div_2exp (t, t, p);
  while (p)
    {
      factors.push_back(2);
      --p;
    }

  p = 3;
  for (i = 1; i < PRIMES_PTAB_ENTRIES;)
    {
      if (! mpz_divisible_ui_p (t, p))
	{
	  p += primes_diff[i++];
	  if (mpz_cmp_ui (t, p * p) < 0)
	    break;
	}
      else
	{
	  mpz_tdiv_q_ui (t, t, p);
	  factors.push_back( p);
	}
    }

  mpz_clear (q);
}

static int
mp_millerrabin (mpz_srcptr n, mpz_srcptr nm1, mpz_ptr x, mpz_ptr y,
		mpz_srcptr q, unsigned long int k)
{
  unsigned long int i;

  mpz_powm (y, x, q, n);

  if (mpz_cmp_ui (y, 1) == 0 || mpz_cmp (y, nm1) == 0)
    return 1;

  for (i = 1; i < k; i++)
    {
      mpz_powm_ui (y, y, 2, n);
      if (mpz_cmp (y, nm1) == 0)
	return 1;
      if (mpz_cmp_ui (y, 1) == 0)
	return 0;
    }
  return 0;
}

int
mp_prime_p (mpz_t n)
{
  int k, r, is_prime;
  mpz_t q, a, nm1, tmp;

  bigvec factors;

  if (mpz_cmp_ui (n, 1) <= 0)
    return 0;

  /* We have already casted out small primes. */
  if (mpz_cmp_ui (n, (long) FIRST_OMITTED_PRIME * FIRST_OMITTED_PRIME) < 0)
    return 1;

  mpz_init (q);
  mpz_init(a);
  mpz_init( nm1);
  mpz_init(tmp);

  /* Precomputation for Miller-Rabin.  */
  mpz_sub_ui (nm1, n, 1);

  /* Find q and k, where q is odd and n = 1 + 2**k * q.  */
  k = mpz_scan1 (nm1, 0);
  mpz_tdiv_q_2exp (q, nm1, k);

  mpz_set_ui (a, 2);

  /* Perform a Miller-Rabin test, finds most composites quickly.  */
  if (!mp_millerrabin (n, nm1, a, tmp, q, k))
    {
      is_prime = 0;
      goto ret2;
    }

  if (flag_prove_primality)
    {
      /* Factor n-1 for Lucas.  */
      mpz_set (tmp, nm1);
      getPrimeFactors (tmp, factors);
    }

  /* Loop until Lucas proves our number prime, or Miller-Rabin proves our
     number composite.  */
  for (r = 0; r < PRIMES_PTAB_ENTRIES; r++)
    {
      int i;

      if (flag_prove_primality)
	{
	  is_prime = 1;
	  for (i = 0; i < factors.size() && is_prime; i++)
	    {
	      mpz_divexact (tmp, nm1, factors[i].value.getValue());
	      mpz_powm (tmp, a, tmp, n);
	      is_prime = mpz_cmp_ui (tmp, 1) != 0;
	    }
	}
      else
	{
	  /* After enough Miller-Rabin runs, be content. */
	  is_prime = (r == MR_REPS - 1);
	}

      if (is_prime)
	goto ret1;

      mpz_add_ui (a, a, primes_diff[r]);	/* Establish new base.  */

      if (!mp_millerrabin (n, nm1, a, tmp, q, k))
	{
	  is_prime = 0;
	  goto ret1;
	}
    }

  error( "Lucas prime test failure.  This should not happen\n");
 ret1:
  if (flag_prove_primality)
    factors.resize(0);

 ret2:
  mpz_clear (q);
  mpz_clear (a);
  mpz_clear(nm1);
  mpz_clear(tmp);

  return is_prime;
}

void
factor_using_pollard_rho (mpz_t n, unsigned long a, bigvec & factors)
{
  mpz_t x, z, y, P;
  mpz_t t, t2;
  unsigned long  k, l, i;

  mpz_init (t);
  mpz_init(t2);
  mpz_init_set_si (y, 2);
  mpz_init_set_si (x, 2);
  mpz_init_set_si (z, 2);
  mpz_init_set_ui (P, 1);
  k = 1;
  l = 1;

  while (mpz_cmp_ui (n, 1) != 0)
    {
      for (;;)
	{
	  do
	    {
	      mpz_mul (t, x, x);
	      mpz_mod (x, t, n);
	      mpz_add_ui (x, x, a);

	      mpz_sub (t, z, x);
	      mpz_mul (t2, P, t);
	      mpz_mod (P, t2, n);

	      if (k % 32 == 1)
		{
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
	  for (i = 0; i < k; i++)
	    {
	      mpz_mul (t, x, x);
	      mpz_mod (x, t, n);
	      mpz_add_ui (x, x, a);
	    }
	  mpz_set (y, x);
	}

    factor_found:
      do
	{
	  mpz_mul (t, y, y);
	  mpz_mod (y, t, n);
	  mpz_add_ui (y, y, a);

	  mpz_sub (t, z, y);
	  mpz_gcd (t, t, n);
	}
      while (mpz_cmp_ui (t, 1) == 0);

      mpz_divexact (n, n, t);	/* divide by t, before t is overwritten */

      if (!mp_prime_p (t))
	{
	  factor_using_pollard_rho (t, a + 1, factors);
	}
      else
	{
	  factors.push_back( t );
	}

      if (mp_prime_p (n))
	{
	  factors.push_back( n);
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

void getPrimeFactors (mpz_t t, bigvec & factors) {
  if (mpz_sgn (t) != 0) {
      factor_using_division (t, factors);
        if (mpz_cmp_ui (t, 1) != 0) {
    	    if (mp_prime_p (t))
    	        factors.push_back( t);
    	    else
    	        factor_using_pollard_rho (t, 1, factors);
	    }
    }
    std::sort(factors.value.begin(), factors.value.end());
}