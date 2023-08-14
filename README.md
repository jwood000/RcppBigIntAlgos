[![CRAN status](<https://www.r-pkg.org/badges/version/RcppBigIntAlgos>)](<https://cran.r-project.org/package=RcppBigIntAlgos>)
[![R build status](<https://github.com/jwood000/RcppBigIntAlgos/workflows/R-CMD-check/badge.svg>)](<https://github.com/jwood000/RcppBigIntAlgos/actions>)
![](<http://cranlogs.r-pkg.org/badges/RcppBigIntAlgos?color=orange>)
![](<http://cranlogs.r-pkg.org/badges/grand-total/RcppBigIntAlgos?color=brightgreen>)
[![Dependencies](<https://tinyverse.netlify.com/badge/RcppBigIntAlgos>)](<https://cran.r-project.org/package=RcppBigIntAlgos>)
[![Coverage status](<https://codecov.io/github/jwood000/RcppBigIntAlgos/branch/main/graph/badge.svg?token=LMO4DH4OtN>)](<https://app.codecov.io/github/jwood000/RcppBigIntAlgos>)
[![Codacy Badge](<https://app.codacy.com/project/badge/Grade/3be4c3c9e3554125b8cc0e13decaf95c>)](<https://app.codacy.com/gh/jwood000/RcppBigIntAlgos/dashboard?utm_source=github.com&utm_medium=referral&utm_content=jwood000/RcppBigIntAlgos&utm_campaign=Badge_Grade>)

# RcppBigIntAlgos

## Overview

`RcppBigIntAlgos` is an `R` package for efficiently factoring large integers. It is multithreaded and uses the C library GMP (GNU Multiple Precision Arithmetic). For very large integers, prime factorization is carried out by a variant of the quadratic sieve algorithm that implements multiple polynomials. For smaller integers, a simple elliptic curve algorithm is attempted followed by a constrained version of the Pollard’s rho algorithm (original code from <https://gmplib.org/>… this is the same algorithm found in the [R gmp package](<https://CRAN.R-project.org/package=gmp>) called by the function `factorize`). Finally, one can quickly obtain a complete factorization of a given number `n` via `divisorsBig`.

## Installation

``` r
## install.packages("RcppBigIntAlgos")

## Or install the development version
## devtools::install_github("jwood000/RcppBigIntAlgos")
```

## The Quadratic Sieve

The function `quadraticSieve` implements the multiple polynomial quadratic sieve algorithm. Currently, `quadraticSieve` can comfortably factor numbers with less than 70 digits (~230 bits) on most standard personal computers. If you have access to powerful computers with many cores, factoring 100+ digit semiprimes in less than a day is not out of the question. Note, the function `primeFactorizeBig(n, skipECM = T, skipPolRho = T)` is the same as `quadraticSieve(n)`.

``` r
library(gmp)
#> 
#> Attaching package: 'gmp'
#> The following objects are masked from 'package:base':
#> 
#>     %*%, apply, crossprod, matrix, tcrossprod
library(RcppBigIntAlgos)

## Generate large semi-primes
semiPrime120bits <- prod(nextprime(urand.bigz(2, 60, 42)))
#> Seed default initialisation
#> Seed initialisation
semiPrime130bits <- prod(nextprime(urand.bigz(2, 65, 1)))
#> Seed initialisation
semiPrime140bits <- prod(nextprime(urand.bigz(2, 70, 42)))
#> Seed initialisation

## The 120 bit number is 36 digits
nchar(as.character(semiPrime120bits))
#> [1] 36

## The 130 bit number is 39 digits
nchar(as.character(semiPrime130bits))
#> [1] 39

## The 140 bit number is 42 digits
nchar(as.character(semiPrime140bits))
#> [1] 42

## Using factorize from gmp package which implements pollard's rho algorithm
##**************gmp::factorize*********************
system.time(print(factorize(semiPrime120bits)))
#> Big Integer ('bigz') object of length 2:
#> [1] 638300143449131711  1021796573707617139
#>    user  system elapsed 
#> 100.788   0.291 101.112 

system.time(print(factorize(semiPrime130bits)))
#> Big Integer ('bigz') object of length 2:
#> [1] 14334377958732970351 29368224335577838231
#>      user   system  elapsed 
#>  1253.060    2.612 1256.065 
 
system.time(print(factorize(semiPrime140bits)))
#> Big Integer ('bigz') object of length 2:
#> [1] 143600566714698156857  1131320166687668315849
#>     user   system  elapsed 
#> 1673.666    3.154 1676.838 


##**************quadraticSieve*********************
## quadraticSieve is much faster and scales better
system.time(print(quadraticSieve(semiPrime120bits)))
#> Big Integer ('bigz') object of length 2:
#> [1] 638300143449131711  1021796573707617139
#>    user  system elapsed 
#>   0.041   0.000   0.042

system.time(print(quadraticSieve(semiPrime130bits)))
#> Big Integer ('bigz') object of length 2:
#> [1] 14334377958732970351 29368224335577838231
#>    user  system elapsed 
#>   0.051   0.001   0.052

system.time(print(quadraticSieve(semiPrime140bits)))
#> Big Integer ('bigz') object of length 2:
#> [1] 143600566714698156857  1131320166687668315849
#>    user  system elapsed 
#>   0.077   0.001   0.078
```
### Using Multiple Threads

As of version `0.3.0`, we can utilize multiple threads. Below, are a few examples:

1.  The largest [Cunnaningham Most Wanted](<https://www.lehigh.edu/~bad0/msg06332.html>) number from the first edition released in 1983 in less than 14 seconds.
2.  [RSA-79](<https://members.loria.fr/PZimmermann/records/rsa.html>) under 2 minutes.
3.  A 300-bit (91-digits) semiprime in 1 hour.
4.  [RSA-99](<https://members.loria.fr/PZimmermann/records/rsa.html>) under 5 hours.
5.  [RSA-100](<https://en.wikipedia.org/wiki/RSA-100>) under 10 hours.

Below are my machine specs and R version info:

``` r
## MacBook Air (2022)
## Processor: Apple M2
## Memory; 24 GB

sessionInfo()
#> R version 4.3.1 (2023-06-16)
#> Platform: aarch64-apple-darwin20 (64-bit)
#> Running under: macOS Ventura 13.4.1
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] RcppBigIntAlgos_1.1.0    gmp_0.7-2
#> 
#> loaded via a namespace (and not attached):
#> [1] compiler_4.3.1

## Maximum number of available threads
stdThreadMax()
#> [1] 8
```

#### Most Wanted 1983

```r
mostWanted1983 <- as.bigz(div.bigz(sub.bigz(pow.bigz(10, 71), 1), 9))
quadraticSieve(mostWanted1983, showStats = TRUE, nThreads = 8)
#> 
#> Summary Statistics for Factoring:
#>     11111111111111111111111111111111111111111111111111111111111111111111111
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |      10s 940ms     |   100%   |    16183    |    4089    |    4345    |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |      2s 107ms      |     8301 x 8434    |
#> 
#> |     Total Time     |
#> |--------------------|
#> |      13s 173ms     |
#> 
#> Big Integer ('bigz') object of length 2:
#> [1] 241573142393627673576957439049           
#> [2] 45994811347886846310221728895223034301839
```

#### RSA-79

```r
rsa79 <- as.bigz("7293469445285646172092483905177589838606665884410340391954917800303813280275279")
quadraticSieve(rsa79, showStats = TRUE, nThreads = 8)
#> 
#> Summary Statistics for Factoring:
#>     7293469445285646172092483905177589838606665884410340391954917800303813280275279
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |     1m 35s 89ms    |   100%   |    91221    |    5651    |    7096    |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |      5s 175ms      |    12625 x 12747   |
#> 
#> |     Total Time     |
#> |--------------------|
#> |    1m 40s 558ms    |
#> 
#> Big Integer ('bigz') object of length 2:
#> [1] 848184382919488993608481009313734808977 
#> [2] 8598919753958678882400042972133646037727
```

### Random 300-bit Semiprime

```r
semi_prime_300_bit <- prod(nextprime(urand.bigz(nb = 2, size = 150, seed = 42)))
#> Seed default initialisation
#> Seed initialisation

nchar(as.character(semi_prime_300_bit))
#> [1] 91

quadraticSieve(semi_prime_300_bit, showStats=TRUE, nThreads=8)
#> 
#> Summary Statistics for Factoring:
#>     1598678911004402782180963020655649301157676037614983547537229086778878619660017752047003277
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |   1h 3m 4s 604ms   |   100%   |   2447868   |    6835    |    12143   |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |      13s 202ms     |    18892 x 18978   |
#> 
#> |     Total Time     |
#> |--------------------|
#> |   1h 3m 18s 733ms  |
#> 
#> Big Integer ('bigz') object of length 2:
#> [1] 1262920060924323380524693068285713630353017593
#> [2] 1265859146963220756974764147023611562195937589
```

#### RSA-99

```r
rsa99 <- "256724393281137036243618548169692747168133997830674574560564321074494892576105743931776484232708881"
 
quadraticSieve(rsa99, showStats=TRUE, nThreads=8)
#> 
#> Summary Statistics for Factoring:
#>     256724393281137036243618548169692747168133997830674574560564321074494892576105743931776484232708881
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |  4h 41m 28s 410ms  |   100%   |   7351308   |    9203    |    15846   |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |      31s 359ms     |    24929 x 25049   |
#> 
#> |     Total Time     |
#> |--------------------|
#> |   4h 42m 1s 55ms   |
#> 
#> Big Integer ('bigz') object of length 2:
#> [1] 4868376167980921239824329271069101142472222111193 
#> [2] 52733064254484107837300974402288603361507691060217
```

#### RSA-100

```r
rsa100 <- "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139"

quadraticSieve(rsa100, showStats = TRUE, nThreads = 8)
#> 
#> Summary Statistics for Factoring:
#>     1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |   9h 43m 9s 13ms   |   100%   |   14999534  |    9230    |    16665   |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |      34s 414ms     |    25799 x 25895   |
#> 
#> |     Total Time     |
#> |--------------------|
#> |  9h 43m 45s 242ms  |
#> 
#> Big Integer ('bigz') object of length 2:
#> [1] 37975227936943673922808872755445627854565536638199
#> [2] 40094690950920881030683735292761468389214899724061
```

## **`primeFactorizeBig`**

As of version `1.0.0`, we can take advantage of the power of Lenstra’s elliptic curve method. This method is particularly useful for quickly finding smaller prime factors of very large composite numbers. It is automatically utilized in the vectorized prime factorization function `primeFactorizeBig`. This function should be preferred in most situations and is identical to `quadraticSieve` when both `skipECM` and `skipPolRho` are set to `TRUE`.

It is optimized for factoring big and small numbers by dynamically using different algorithms based off of the input. It takes cares to not spend too much time in any of the methods and avoids wastefully switching to the quadratic sieve when the number is very large.

For example, using the defaults on `mostWanted1983` above only adds a few seconds.

```r
primeFactorizeBig(mostWanted1983, showStats = TRUE, nThreads = 8)
#> 
#> Summary Statistics for Factoring:
#>     11111111111111111111111111111111111111111111111111111111111111111111111
#> 
#> |  Pollard Rho Time  |
#> |--------------------|
#> |        90ms        |
#> 
#> |  Lenstra ECM Time  |  Number of Curves  |
#> |--------------------|--------------------|
#> |      2s 129ms      |        4181        |
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |      10s 637ms     |   100%   |    16183    |    4089    |    4345    |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |       2s 73ms      |     8301 x 8434    |
#> 
#> |     Total Time     |
#> |--------------------|
#> |       15s 7ms      |
#> 
#> Big Integer ('bigz') object of length 2:
#> [1] 241573142393627673576957439049           
#> [2] 45994811347886846310221728895223034301839
```


### The Lentra’s power

Performing a few iterations of the Pollard’s rho algorithm followed by the quadratic sieve can prove to be a terrible solution when the input size is incredibly large. More than likely, your constrained Pollard’s rho algorithm won’t find any prime factors, and eventually will pass an enormous composite number to the quadratic sieve. The quadratic sieve isn’t optimized for finding *smallish* factors in large composites. It will treat the input similarly to a semiprime of the same size. For example, the number `ecmExpo1` below is 428 digits and if passed to the quadratic sieve algorithm, it could take years to factor as not even the most sophisticated semiprime algorithms can easily factor numbers of this size (See <https://en.wikipedia.org/wiki/RSA_Factoring_Challenge>). However, with the ECM, this is light work.

``` r
ecmExpo1 <- pow.bigz(prod(nextprime(urand.bigz(10, 47, 13))), 3) *
            nextprime(urand.bigz(1, 47, 123))
#> Seed initialisation
#> Seed initialisation

system.time(primeFactorizeBig(ecmExpo1, skipPolRho = TRUE, nThreads = 8))
#>    user  system elapsed 
#>  10.768   0.024   2.840

ecmExpo2 <- pow.bigz(prod(nextprime(urand.bigz(10, 40, 42))), 3) *
            pow.bigz(prod(nextprime(urand.bigz(10, 45, 42))), 5) *
            nextprime(urand.bigz(1, 80, 123)) *
            nextprime(urand.bigz(1, 90, 123))
#> Seed initialisation
#> Seed initialisation
#> Seed initialisation
#> Seed initialisation

ecm2 <- primeFactorizeBig(ecmExpo2, nThreads = 8, showStats = TRUE)
#> 
#> Summary Statistics for Factoring:
#>     215590243996472403826986190959172968924582673399860040722565967228161406366421528160220759715319487748182459438846599849086479492406371823155743767776724831907786824847947911279318978691074995051134928403920242678621702805852530836549215396392496547012632510254095452820560897056877208839289455859596425985161415486708317667524364034985007498553043810662487802224636724842739870398138687583364019516727138808818929150447675931521160760469848585726320428719422689614460288917155690252498420348768667307236940996727406088066077174182630522224530495122361513443818343671556622023320404680434511134662288782134897895034832405624062328931996299632052035105693354967480119521413889950571647670948619957219258020882926218179673783234886719528724327466781424377562596208305773866403070803284721517038640418360492984835882655310830736014378303508492181368455405610072971016964625940147583833553074577246736323251622391837090678894432642161330574148459946482734705537584444739753686240736188950594612747228733363204312823260997792251300802876878725604381077823
#> 
#> |  Pollard Rho Time  |
#> |--------------------|
#> |      13s 163ms     |
#> 
#> |  Lenstra ECM Time  |  Number of Curves  |
#> |--------------------|--------------------|
#> |      15s 256ms     |        24565       |
#> 
#> |      MPQS Time     | Complete | Polynomials |   Smooths  |  Partials  |
#> |--------------------|----------|-------------|------------|------------|
#> |        419ms       |   100%   |     1482    |    1013    |    1450    |
#> 
#> |  Mat Algebra Time  |    Mat Dimension   |
#> |--------------------|--------------------|
#> |        60ms        |     2425 x 2463    |
#> 
#> |     Total Time     |
#> |--------------------|
#> |      28s 648ms     |

head(ecm2)
#> Big Integer ('bigz') object of length 6:
#> [1] 66642484459  66642484459  66642484459  385217678221 385217678221
#> [6] 385217678221

tail(ecm2)
#> Big Integer ('bigz') object of length 6:
#> [1] 34936543827863               34936543827863              
#> [3] 34936543827863               34936543827863              
#> [5] 1130548241045557299883517    1051687085486158310119550449

length(ecm2)
#> [1] 82
```

### Safely Interrupt Execution

If you want to interrupt a command which will take a long time, hit Ctrl + c, or esc if using RStudio, to stop execution. When you utilize multiple threads with a very large number (e.g. 90 digit semiprime), you will be able to interrupt execution once every ~30 seconds.

``` r
## User hits Ctrl + c
## system.time(quadraticSieve(prod(nextprime(urand.bigz(2, 100, 42)))))
## Seed initialisation
## Error in PrimeFactorization(n, FALSE, showStats, TRUE, TRUE, nThreads,  :
##   C++ call interrupted by the user.
## Timing stopped at: 2.164 0.039 2.167
```

## Complete Factorization with **`divisorsBig`**

This function generates the complete factorization for many (possibly large) numbers. It is vectorized and can also return a named list.

``` r
## Get all divisors of a given number:
divisorsBig(1000)
#> Big Integer ('bigz') object of length 16:
#>  [1] 1    2    4    5    8    10   20   25   40   50   100  125  200  250  500 
#> [16] 1000

## Or, get all divisors of a vector:
divisorsBig(urand.bigz(nb = 2, size = 100, seed = 42), namedList = TRUE)
#> Seed initialisation
#> $`153675943236425922379228498617`
#> Big Integer ('bigz') object of length 16:
#>  [1] 1                              3                             
#>  [3] 7                              9                             
#>  [5] 21                             27                            
#>  [7] 63                             189                           
#>  [9] 813100228764158319466817453    2439300686292474958400452359  
#> [11] 5691701601349108236267722171   7317902058877424875201357077  
#> [13] 17075104804047324708803166513  21953706176632274625604071231 
#> [15] 51225314412141974126409499539  153675943236425922379228498617
#> 
#> $`261352009818227569107309994396`
#> Big Integer ('bigz') object of length 12:
#>  [1] 1                              2                             
#>  [3] 4                              155861                        
#>  [5] 311722                         623444                        
#>  [7] 419206873140534785974859       838413746281069571949718      
#>  [9] 1676827492562139143899436      65338002454556892276827498599 
#> [11] 130676004909113784553654997198 261352009818227569107309994396
```

### Efficiency

It is very efficient as well. It is equipped with a modified merge sort algorithm that significantly outperforms the `std::sort`/`bigvec` (the class utilized in the `R gmp` package) combination.

``` r
hugeNumber <- pow.bigz(2, 100) * pow.bigz(3, 100) * pow.bigz(5, 100)
system.time(overOneMillion <- divisorsBig(hugeNumber))
#>    user  system elapsed 
#>   0.117   0.008   0.126

length(overOneMillion)
#> [1] 1030301

## Output is in ascending order
tail(overOneMillion)
#> Big Integer ('bigz') object of length 6:
#> [1] 858962534553352218394101882942702121170179203335000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 
#> [2] 1030755041464022662072922259531242545404215044002000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#> [3] 1288443801830028327591152824414053181755268805002500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#> [4] 1717925069106704436788203765885404242340358406670000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#> [5] 2576887603660056655182305648828106363510537610005000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#> [6] 5153775207320113310364611297656212727021075220010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
```

### Correct Ordering

Another benefit is that it will return correct orderings on extremely large numbers when compared to sorting large vectors in `base R`. Typically in `base R` you must execute the following: `order(asNumeric(myVectorHere))`. When the numbers get large enough, precision is lost which leads to incorrect orderings. Observe:

``` r
set.seed(101)
testBaseSort <- do.call(
    c, lapply(sample(100), \(x) add.bigz(pow.bigz(10,80), x))
)
testBaseSort <- testBaseSort[order(asNumeric(testBaseSort))]
myDiff <- do.call(
    c, lapply(1:99, \(x) sub.bigz(testBaseSort[x+1], testBaseSort[x]))
)

## Should return integer(0) as the difference should always be positive
## NOTE that the result will be unpredictable because of lack of precision
which(myDiff < 0)
#>  [1]  1  2  4  5  7  8  9 11 14 16 17 19 21 22 24 26 28 29 31 33 37 38 40 41 43
#> [26] 45 47 48 50 51 52 54 55 56 58 59 61 63 65 67 68 69 71 74 76 79 80 83 85 88
#> [51] 89 91 93 94 96 98 99

## N.B. The first and second elements are incorrect order (among others)
head(testBaseSort)
#> Big Integer ('bigz') object of length 6:
#> [1] 100000000000000000000000000000000000000000000000000000000000000000000000000000073
#> [2] 100000000000000000000000000000000000000000000000000000000000000000000000000000057
#> [3] 100000000000000000000000000000000000000000000000000000000000000000000000000000046
#> [4] 100000000000000000000000000000000000000000000000000000000000000000000000000000095
#> [5] 100000000000000000000000000000000000000000000000000000000000000000000000000000081
#> [6] 100000000000000000000000000000000000000000000000000000000000000000000000000000058
```

## Acknowledgments and Resources

  * Credit to [primo](<https://codegolf.stackexchange.com/users/4098/primo>) (Mike Tryczak) and his excellent answer to [Fastest semiprime factorization](<https://codegolf.stackexchange.com/a/9088/52987>).

  * [Factoring large numbers with quadratic sieve](<https://blogs.msdn.microsoft.com/devdev/2006/06/19/factoring-large-numbers-with-quadratic-sieve/>) on MSDN Archive.

  * A really nice concise example is given here: [Factorization of *n = 87463* with the Quadratic Sieve](<https://www.math.colostate.edu/~hulpke/lectures/m400c/quadsievex.pdf>)

  * [Smooth numbers and the quadratic sieve](<http://library.msri.org/books/Book44/files/03carl.pdf>) by Carl Pomerance

  * [Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware](<https://www.iacr.org/archive/ches2006/10/10.pdf>) by Gaj K. et al.

  * [Integer Factorization using the Quadratic Sieve](<https://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf>) by Chad Seibert

  * In the stackoverflow question and answer [What is the most efficient factoring algorithm for quadratic sieve extraction phase?](<https://stackoverflow.com/q/63541365/4408538>) by [Ilya Gazman](<https://github.com/gazman-sdk>), an efficient method for checking divisibility is sketched out that utilizes built-in types. You can see more on a video Ilya put on youtube: [E15: Quadratic Sieve Running on Java - Receiving](<https://youtu.be/sXg_WrCUX-Q>). While `mpz_divisible_ui_p` is very efficient, we found better performance using this method.

  * [R Function for returning ALL factors](<https://stackoverflow.com/a/49742904/4408538>)

  * [Issues factoring large number that is 99 digits long](<https://stackoverflow.com/a/66128627/4408538>)


## Current Research

Currently, our main focus for version `1.2.0` will be implementing the self initializing quadratic sieve.

## Contact

If you would like to report a bug, have a question, or have suggestions for possible improvements, please file an [issue](<https://github.com/jwood000/RcppBigIntAlgos/issues>).
