![](http://cranlogs.r-pkg.org/badges/bigIntegersAlgos?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/bigIntegersAlgos?color=brightgreen)

# bigIntegerAlgos

Overview
---------
bigIntegerAlgos uses the C library GMP (GNU Multiple Precision Arithmetic) for efficiently
factoring big integers. For very large integers, prime factorization is carried out by a variant of the quadratic sieve algorithm that implements multiple polynomials. For smaller integers, a constrained version of the Pollard's rho algorithm is used (original code from https://gmplib.org/... this is the same algorithm found in the R gmp package (https://cran.r-project.org/web/packages/gmp/gmp.pdf) called by the function `factorize`). Finally, one can quickly obtain a complete factorization of given number `n` via `divisorsBig`.

Installation
------------

``` r
install.packages("bigIntegerAlgos")

## Or install the development version
devtools::install_github("jwood000/bigIntegerAlgos")
```

Usage
-----
First, we take a look at `divisorsBig`. It is vectorized and can also return a named list.
``` r
## Get all divisors of a given number:
divisorsBig(1000)
Big Integer ('bigz') object of length 16:
 [1] 1    2    4    5    8    10   20   25   40   50   100  125  200  250  500  1000
 
 
 ## Or, get all divisors of a vector:
divisorsBig(urand.bigz(nb = 2, size = 100, seed = 42), namedList = TRUE)
Seed initialisation
$`153675943236425922379228498617`
Big Integer ('bigz') object of length 16:
 [1] 1                              3                             
 [3] 7                              9                             
 [5] 21                             27                            
 [7] 63                             189                           
 [9] 813100228764158319466817453    2439300686292474958400452359  
[11] 5691701601349108236267722171   7317902058877424875201357077  
[13] 17075104804047324708803166513  21953706176632274625604071231 
[15] 51225314412141974126409499539  153675943236425922379228498617

$`261352009818227569107309994396`
Big Integer ('bigz') object of length 12:
 [1] 1                              2                             
 [3] 4                              155861                        
 [5] 311722                         623444                        
 [7] 419206873140534785974859       838413746281069571949718      
 [9] 1676827492562139143899436      65338002454556892276827498599 
[11] 130676004909113784553654997198 261352009818227569107309994396
```

It is very efficient as well. It is equipped with a modified merge sort algorithm that significantly outperforms the `std::sort`/`bigvec` (the class utilized in the `R gmp` package) combination.

```r
hugeNumber <- pow.bigz(2, 100) * pow.bigz(3, 100) * pow.bigz(5, 100)
system.time(overOneMillion <- divisorsBig(hugeNumber))
   user  system elapsed 
  0.637   0.043   0.682 
length(overOneMillion)
[1] 1030301

## Output is in ascending order
tail(overOneMillion)
Big Integer ('bigz') object of length 6:
[1] 858962534553352218394101882942702121170179203335000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 
[2] 1030755041464022662072922259531242545404215044002000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
[3] 1288443801830028327591152824414053181755268805002500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
[4] 1717925069106704436788203765885404242340358406670000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
[5] 2576887603660056655182305648828106363510537610005000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
[6] 5153775207320113310364611297656212727021075220010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
```
 Another benefit is that it will return correct orderings on extremely large numbers when compared to sorting large vectors in `base R`.  Typically in `base R` you must execute the following: `order(asNumeric(myVectorHere))`. When the numbers get large enough, precision is lost which leads to incorrect orderings. Observe:
 
```r
set.seed(101)
testBaseSort <- do.call(c, lapply(sample(100), function(x) add.bigz(pow.bigz(10,80), x)))
testBaseSort <- testBaseSort[order(asNumeric(testBaseSort))]
myDiff <- do.call(c, lapply(1:99, function(x) sub.bigz(testBaseSort[x+1], testBaseSort[x])))

## Should return integer(0) as the difference should always be positive
which(myDiff < 0)
 [1]  1  3  4  7  9 11 14 17 19 22 24 25 26 28 31 32 33 36 37 38 40 42 45 47 48
[26] 50 51 54 57 58 59 63 64 65 66 69 70 72 75 78 81 82 85 87 89 91 93 94 97 98

## N.B. The first and second elements are incorrect order (among others)
head(testBaseSort)
Big Integer ('bigz') object of length 6:
[1] 100000000000000000000000000000000000000000000000000000000000000000000000000000038
[2] 100000000000000000000000000000000000000000000000000000000000000000000000000000005
[3] 100000000000000000000000000000000000000000000000000000000000000000000000000000070
[4] 100000000000000000000000000000000000000000000000000000000000000000000000000000064
[5] 100000000000000000000000000000000000000000000000000000000000000000000000000000024
[6] 100000000000000000000000000000000000000000000000000000000000000000000000000000029

```

The Quadratic Sieve
--------------
The function `quadraticSieve` implements the multiple polynomial quadratic sieve algorithm. Currently, `quadraticSieve` can comfortably factor numbers with less than 50 digits (~160 bits).

```r
## Generate large semi-primes
semiPrime120bits <- prod(nextprime(urand.bigz(2,60,42)))
semiPrime130bits <- prod(nextprime(urand.bigz(2,65,1)))
semiPrime140bits <- prod(nextprime(urand.bigz(2,70,42)))

## Using factorize from gmp package which implements pollard's rho algorithm
## We did not test the 140 bit semi-prime as the 130 bit took a very long time

##**************gmp::factorize*********************
system.time(print(factorize(semiPrime120bits)))
Big Integer ('bigz') object of length 2:
[1] 638300143449131711  1021796573707617139
   user  system elapsed 
126.603   0.052 126.694

system.time(print(factorize(semiPrime130bits)))
Big Integer ('bigz') object of length 2:
[1] 14334377958732970351 29368224335577838231
    user   system  elapsed 
1513.055    1.455 1517.524

##**************quadraticSieve*********************
## quadraticSieve is much faster and scales better
system.time(print(quadraticSieve(semiPrime120bits)))
Big Integer ('bigz') object of length 2:
[1] 638300143449131711  1021796573707617139
   user  system elapsed 
  5.055   0.013   5.069 
  
system.time(print(quadraticSieve(semiPrime130bits)))
Big Integer ('bigz') object of length 2:
[1] 14334377958732970351 29368224335577838231
   user  system elapsed 
  7.570   0.010   7.583

system.time(print(quadraticSieve(semiPrime140bits)))
Big Integer ('bigz') object of length 2:
[1] 143600566714698156857  1131320166687668315849
   user  system elapsed 
 14.308   0.029  14.338
```

It can also be used as a general prime factoring function:

```r
quadraticSieve(urand.bigz(1,50,1))
Seed initialisation
Big Integer ('bigz') object of length 5:
[1] 5       31      307     2441    4702723
```

However `gmp::factorize` is more suitable for numbers smaller than 60 bits and should be used in such cases.

Current Research:
-----
Improvements are being made to the section of the quadratic sieve algorithm that selects smooth numbers. Right now, we are only selecting **_M-smooth_** numbers where **_M_** is the largest prime in our sieving base. A potential efficiency gain is to also keep track of mostly **_M-smooth_** numbers (i.e. numbers that are almost completely factored by our sieving base but have one factor that contains a prime number greater than **_M_**). This way, we can obtain more smooth numbers by essentially eliminating these larger factors when we find a pair of them (N.B. square terms come out in the wash, so the large factors won't influence the outcome).

We are also working on efficiently integrating `divisorsBig` with `quadraticSieve` as currently, `divisorsBig` utilizes `gmp::factorize`.

Contact
----
I welcome any and all feedback. If you would like to report a bug, have a question, or have suggestions for possible improvements, please contact me here: jwood000@gmail.com
