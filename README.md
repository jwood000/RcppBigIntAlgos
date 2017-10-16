# bigIntegerAlgos

Overview
---------
bigIntegerAlgos uses the C library GMP (GNU Multiple Precision Arithmetic) for efficiently
factoring big integers. For very large integers, prime factorization is carried out by a variant of the quadratic sieve algorithm that implements multiple polynomials. For smaller integers, the Pollard's rho algorithm is used (original code from https://gmplib.org/... this is the same
algorithm found in the R gmp package (https://cran.r-project.org/web/packages/gmp/gmp.pdf) called by the function `factorize`). Finally, one can quickly obtain a complete factorization of given number `n` via `getDivisors`.

Usage
-----
``` r
## get all divisors of a given integer.
getDivisors(1000)
Big Integer ('bigz') object of length 16:
 [1] 1    2    4    5    8    10   20   25   40   50   100  125  200  250  500  1000
 
```
