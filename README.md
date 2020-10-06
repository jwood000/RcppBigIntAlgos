[![CRAN status](https://www.r-pkg.org/badges/version/RcppBigIntAlgos)](https://cran.r-project.org/package=RcppBigIntAlgos)
[![Travis build status](https://travis-ci.com/jwood000/RcppBigIntAlgos.svg?branch=master)](https://travis-ci.com/jwood000/RcppBigIntAlgos)
![](http://cranlogs.r-pkg.org/badges/RcppBigIntAlgos?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/RcppBigIntAlgos?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/RcppBigIntAlgos)](https://cran.r-project.org/package=RcppBigIntAlgos)
[![Coverage status](https://codecov.io/gh/jwood000/RcppBigIntAlgos/branch/master/graph/badge.svg)](https://codecov.io/github/jwood000/RcppBigIntAlgos?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3be4c3c9e3554125b8cc0e13decaf95c)](https://app.codacy.com/manual/jwood000/RcppBigIntAlgos?utm_source=github.com&utm_medium=referral&utm_content=jwood000/RcppBigIntAlgos&utm_campaign=Badge_Grade)

# RcppBigIntAlgos

## Overview

RcppBigIntAlgos uses the C library GMP (GNU Multiple Precision Arithmetic) for efficiently
factoring big integers. Links to `RcppThread` for factoring in parallel. For very large integers, prime factorization is carried out by a variant of the quadratic sieve algorithm that implements multiple polynomials. For smaller integers, a constrained version of the Pollard's rho algorithm is used (original code from <https://gmplib.org/>... this is the same algorithm found in the [R gmp package](<https://CRAN.R-project.org/package=gmp>) called by the function `factorize`). Finally, one can quickly obtain a complete factorization of a given number `n` via `divisorsBig`.

## Installation

``` r
install.packages("RcppBigIntAlgos")

## Or install the development version
devtools::install_github("jwood000/RcppBigIntAlgos")
```

## Usage

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

### Efficiency

It is very efficient as well. It is equipped with a modified merge sort algorithm that significantly outperforms the `std::sort`/`bigvec` (the class utilized in the `R gmp` package) combination.

```r
hugeNumber <- pow.bigz(2, 100) * pow.bigz(3, 100) * pow.bigz(5, 100)
system.time(overOneMillion <- divisorsBig(hugeNumber))
   user  system elapsed 
  0.364   0.029   0.390
  
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

### Correct Ordering

 Another benefit is that it will return correct orderings on extremely large numbers when compared to sorting large vectors in `base R`.  Typically in `base R` you must execute the following: `order(asNumeric(myVectorHere))`. When the numbers get large enough, precision is lost which leads to incorrect orderings. Observe:
 
```r
set.seed(101)
testBaseSort <- do.call(c, lapply(sample(100), function(x) add.bigz(pow.bigz(10,80), x)))
testBaseSort <- testBaseSort[order(asNumeric(testBaseSort))]
myDiff <- do.call(c, lapply(1:99, function(x) sub.bigz(testBaseSort[x+1], testBaseSort[x])))

## Should return integer(0) as the difference should always be positive
## NOTE that the result will be unpredictable because of lack of precision
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

## The Quadratic Sieve

The function `quadraticSieve` implements the multiple polynomial quadratic sieve algorithm. Currently, `quadraticSieve` can comfortably factor numbers with less than 70 digits (~230 bits) on most standard personal computers. If you have access to powerful computers with many cores, factoring 100+ digit semiprimes in less than a day is not out of the question.

```r
## Generate large semi-primes
semiPrime120bits <- prod(nextprime(urand.bigz(2, 60, 42)))
semiPrime130bits <- prod(nextprime(urand.bigz(2, 65, 1)))
semiPrime140bits <- prod(nextprime(urand.bigz(2, 70, 42)))

## The 120 bit number is 36 digits
nchar(as.character(semiPrime120bits))
[1] 36

## The 130 bit number is 39 digits
nchar(as.character(semiPrime130bits))
[1] 39

## The 140 bit number is 42 digits
nchar(as.character(semiPrime140bits))
[1] 42

## Using factorize from gmp package which implements pollard's rho algorithm
##**************gmp::factorize*********************
system.time(print(factorize(semiPrime120bits)))
Big Integer ('bigz') object of length 2:
[1] 638300143449131711  1021796573707617139
   user  system elapsed 
125.117   0.139 125.113

system.time(print(factorize(semiPrime130bits)))
Big Integer ('bigz') object of length 2:
[1] 14334377958732970351 29368224335577838231
    user   system  elapsed 
1437.246    0.309 1437.505

system.time(print(factorize(semiPrime140bits)))
Big Integer ('bigz') object of length 2:
[1] 143600566714698156857  1131320166687668315849
    user   system  elapsed 
2239.374    0.299 2239.641


##**************quadraticSieve*********************
## quadraticSieve is much faster and scales better
system.time(print(quadraticSieve(semiPrime120bits)))
Big Integer ('bigz') object of length 2:
[1] 638300143449131711  1021796573707617139
   user  system elapsed 
  0.086   0.001   0.085
  
system.time(print(quadraticSieve(semiPrime130bits)))
Big Integer ('bigz') object of length 2:
[1] 14334377958732970351 29368224335577838231
   user  system elapsed 
  0.102   0.000   0.103

system.time(print(quadraticSieve(semiPrime140bits)))
Big Integer ('bigz') object of length 2:
[1] 143600566714698156857  1131320166687668315849
   user  system elapsed 
  0.173   0.000   0.174
```

### Using Multiple Threads

As of version `0.3.0`, we can utilize multiple threads with the help of [RcppThread](https://github.com/tnagler/RcppThread). For example, we factor the largest [Cunnaningham Most Wanted](<https://www.lehigh.edu/~bad0/msg06332.html>) number from the first edition released in 1983 in under a minute and [RSA-79](<https://members.loria.fr/PZimmermann/records/rsa.html>) can be factored in under 6 minutes on my machine. I obtained the best performance when `nThreads = stdThreadMax() / 2`. When the number of threads was maximized, there was a decrease in efficiency probably due to pollution of the cache.

```r
quadraticSieve(mostWanted1983, nThreads=4, skipExtPolRho=TRUE, showStats=TRUE)

Summary Statistics for Factoring:
    11111111111111111111111111111111111111111111111111111111111111111111111

|  Pollard Rho Time  |
|--------------------|
|        52ms        |

|      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
|--------------------|----------|-------------|------------|------------|
|      34s 765ms     |   100%   |    15919    |    4511    |    4461    |

|  Mat Algebra Time  |    Mat Dimension   |
|--------------------|--------------------|
|      5s 123ms      |     8836 x 8972    |

|     Total Time     |
|--------------------|
|      40s 133ms     |

Big Integer ('bigz') object of length 2:
[1] 241573142393627673576957439049            45994811347886846310221728895223034301839


## ***************************************************************************


quadraticSieve(rsa79, showStats=TRUE, nThreads=4, skipExtPolRho=TRUE)

Summary Statistics for Factoring:
    7293469445285646172092483905177589838606665884410340391954917800303813280275279

|  Pollard Rho Time  |
|--------------------|
|        66ms        |

|      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
|--------------------|----------|-------------|------------|------------|
|    5m 19s 965ms    |   100%   |    100725   |    5581    |    7166    |

|  Mat Algebra Time  |    Mat Dimension   |
|--------------------|--------------------|
|      12s 694ms     |    12605 x 12747   |

|     Total Time     |
|--------------------|
|    5m 33s 179ms    |

Big Integer ('bigz') object of length 2:
[1] 848184382919488993608481009313734808977  8598919753958678882400042972133646037727
```

### Factor More Than Just Semiprimes

If you encounter a number that is a product of multiple large primes, the algorithm will recursively factor the number into two numbers until every part is prime.

```r
threePrime195bits <- prod(nextprime(urand.bigz(3, 65, 97)))

quadraticSieve(threePrime195bits, showStats = TRUE)

Summary Statistics for Factoring:
    6634573213431810791169420577087478977215298519759798575509

|  Pollard Rho Time  |
|--------------------|
|        51ms        |

|      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
|--------------------|----------|-------------|------------|------------|
|      11s 273ms     |   100%   |     2914    |    1704    |    2098    |

|  Mat Algebra Time  |    Mat Dimension   |
|--------------------|--------------------|
|        470ms       |     3745 x 3802    |


Summary Statistics for Factoring:
    369498233670465681342232176125551121921

|      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
|--------------------|----------|-------------|------------|------------|
|        115ms       |   100%   |      63     |     597    |     246    |

|  Mat Algebra Time  |    Mat Dimension   |
|--------------------|--------------------|
|        18ms        |      801 x 843     |

|     Total Time     |
|--------------------|
|      11s 997ms     |

Big Integer ('bigz') object of length 3:
[1] 11281626468262639417 17955629036507943829 32752213052784053513
```

### General Prime Factoring

It can also be used as a general prime factoring function:

```r
quadraticSieve(urand.bigz(1, 50, 1))
Seed initialisation
Big Integer ('bigz') object of length 5:
[1] 5       31      307     2441    4702723
```

However `gmp::factorize` is more suitable for numbers smaller than 70 bits (about 22 decimal digits) and should be used in such cases.

## Safely Interrupt Execution in **`quadraticSieve`**

If you want to interrupt a command which will take a long time, hit Ctrl + c, or esc if using RStudio, to stop execution.

```r
## User hits Ctrl + c
## system.time(quadraticSieve(prod(nextprime(urand.bigz(2, 100, 42)))))
## Seed default initialisation
## Seed initialisation
## 
##  Error in QuadraticSieveContainer(n) : C++ call interrupted by the user.
##  
## Timing stopped at: 1.623 0.102 1.726
```

## Acknowledgments and Resources

  * Credit to [primo](<https://codegolf.stackexchange.com/users/4098/primo>) (Mike Tryczak) and his excellent answer to [Fastest semiprime factorization](<https://codegolf.stackexchange.com/a/9088/52987>).
  
  * [Factoring large numbers with quadratic sieve](<https://blogs.msdn.microsoft.com/devdev/2006/06/19/factoring-large-numbers-with-quadratic-sieve/>) on MSDN Archive.
  
  * A really nice concise example is given here: [Factorization of _n = 87463_ with the Quadratic Sieve](<https://www.math.colostate.edu/~hulpke/lectures/m400c/quadsievex.pdf>)
  
  * [Smooth numbers and the quadratic sieve](<http://library.msri.org/books/Book44/files/03carl.pdf>) by Carl Pomerance
  
  * [Integer Factorization using the Quadratic Sieve
](<http://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf>) by Chad Seibert


## Current Research

Currenlty, our main focus is on implementing our sieve in a parallel fashion.

## Contact

I welcome any and all feedback. If you would like to report a bug, have a question, or have suggestions for possible improvements, please file an [issue](<https://github.com/jwood000/RcppBigIntAlgos/issues>).
