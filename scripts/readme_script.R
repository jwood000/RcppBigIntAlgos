reprex::reprex({
    #' [![CRAN status](https://www.r-pkg.org/badges/version/RcppBigIntAlgos)](https://cran.r-project.org/package=RcppBigIntAlgos)
    #' [![R build status](<https://github.com/jwood000/RcppBigIntAlgos/workflows/R-CMD-check/badge.svg>)](<https://github.com/jwood000/RcppBigIntAlgos/actions>)
    #' ![](http://cranlogs.r-pkg.org/badges/RcppBigIntAlgos?color=orange)
    #' ![](http://cranlogs.r-pkg.org/badges/grand-total/RcppBigIntAlgos?color=brightgreen)
    #' [![Dependencies](https://tinyverse.netlify.com/badge/RcppBigIntAlgos)](https://cran.r-project.org/package=RcppBigIntAlgos)
    #' [![Coverage status](<https://codecov.io/github/jwood000/RcppBigIntAlgos/branch/main/graph/badge.svg?token=LMO4DH4OtN>)](<https://app.codecov.io/github/jwood000/RcppBigIntAlgos>)
    #' [![Codacy Badge](<https://app.codacy.com/project/badge/Grade/3be4c3c9e3554125b8cc0e13decaf95c>)](<https://app.codacy.com/gh/jwood000/RcppBigIntAlgos/dashboard?utm_source=github.com&utm_medium=referral&utm_content=jwood000/RcppBigIntAlgos&utm_campaign=Badge_Grade>)
    #'
    #' # RcppBigIntAlgos
    #'
    #' ## Overview
    #'
    #' `RcppBigIntAlgos` is an `R` package for efficiently factoring large integers. It is multithreaded and uses the C library GMP (GNU Multiple Precision Arithmetic). For very large integers, prime factorization is carried out by a variant of the quadratic sieve algorithm that implements multiple polynomials. For smaller integers, a simple elliptic curve algorithm is attempted followed by a constrained version of the Pollard's rho algorithm (original code from <https://gmplib.org/>... this is the same algorithm found in the [R gmp package](<https://CRAN.R-project.org/package=gmp>) called by the function `factorize`). Finally, one can quickly obtain a complete factorization of a given number `n` via `divisorsBig`.
    #'
    #' ## Installation
    #'

    ## install.packages("RcppBigIntAlgos")

    ## Or install the development version
    ## devtools::install_github("jwood000/RcppBigIntAlgos")

    #'
    #' ## The Quadratic Sieve
    #'
    #' The function `quadraticSieve` implements the multiple polynomial quadratic sieve algorithm. Currently, `quadraticSieve` can comfortably factor numbers with less than 70 digits (~230 bits) on most standard personal computers. If you have access to powerful computers with many cores, factoring 100+ digit semiprimes in less than a day is not out of the question. Note, the function `primeFactorizeBig(n, skipECM = T, skipPolRho = T)` is the same as `quadraticSieve(n)`.
    #'

    library(gmp)
    library(RcppBigIntAlgos)

    ## Generate large semi-primes
    semiPrime120bits <- prod(nextprime(urand.bigz(2, 60, 42)))
    semiPrime130bits <- prod(nextprime(urand.bigz(2, 65, 1)))
    semiPrime140bits <- prod(nextprime(urand.bigz(2, 70, 42)))

    ## The 120 bit number is 36 digits
    nchar(as.character(semiPrime120bits))

    ## The 130 bit number is 39 digits
    nchar(as.character(semiPrime130bits))

    ## The 140 bit number is 42 digits
    nchar(as.character(semiPrime140bits))

    ## Using factorize from gmp package which implements pollard's rho algorithm
    ##**************gmp::factorize*********************
    # system.time(print(factorize(semiPrime120bits)))
    #
    # system.time(print(factorize(semiPrime130bits)))
    #
    # system.time(print(factorize(semiPrime140bits)))

    ##**************quadraticSieve*********************
    ## quadraticSieve is much faster and scales better
    system.time(print(quadraticSieve(semiPrime120bits)))

    system.time(print(quadraticSieve(semiPrime130bits)))

    system.time(print(quadraticSieve(semiPrime140bits)))

    #'
    #' ### Using Multiple Threads
    #'
    #' As of version `0.3.0`, we can utilize multiple threads. Below, are a few examples:
    #'
    #'  1. The largest [Cunnaningham Most Wanted](<https://www.lehigh.edu/~bad0/msg06332.html>) number from the first edition released in 1983 in less than 25 seconds.
    #'  2. [RSA-79](<https://members.loria.fr/PZimmermann/records/rsa.html>) under 2 minutes.
    #'  3. A 300-bit (91-digits) semiprime in 1 hour.
    #'  4. [RSA-99](<https://members.loria.fr/PZimmermann/records/rsa.html>) under 5 hours.
    #'  5. [RSA-100](<https://en.wikipedia.org/wiki/RSA-100>) under 10 hours.
    #'
    #' Below are my machine specs and R version info:

    ## MacBook Air (2022)
    ## Processor: Apple M2
    ## Memory; 24 GB

    sessionInfo()

    ## Maximum number of available threads
    stdThreadMax()

    #'
    #' #### Most Wanted 1983
    #'

    mostWanted1983 <- as.bigz(div.bigz(sub.bigz(pow.bigz(10, 71), 1), 9))
    quadraticSieve(mostWanted1983, showStats = TRUE, nThreads = 8)

    #'
    #' #### RSA-79
    #'

    rsa79 <- as.bigz("7293469445285646172092483905177589838606665884410340391954917800303813280275279")
    quadraticSieve(rsa79, showStats = TRUE, nThreads = 8)

    #'
    #' #### RSA-99
    #'

    # rsa99 <- "256724393281137036243618548169692747168133997830674574560564321074494892576105743931776484232708881"
    # quadraticSieve(rsa99, showStats = TRUE, nThreads = 8)

    #'
    #' #### RSA-100
    #'

    # rsa100 <- "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139"
    # quadraticSieve(rsa100, showStats = TRUE, nThreads = 8)

    #'
    #' ## **`primeFactorizeBig`**
    #'
    #' As of version `1.0.0`, we can take advantage of the power of Lenstra's elliptic curve method. This method is particularly useful for quickly finding smaller prime factors of very large composite numbers. It is automatically utilized in the vectorized prime factorization function `primeFactorizeBig`. This function should be preferred in most situations and is identical to `quadraticSieve` when both `skipECM` and `skipPolRho` are set to `TRUE`.
    #'
    #' It is optimized for factoring big and small numbers by dynamically using different algorithms based off of the input. It takes cares to not spend too much time in any of the methods and avoids wastefully switching to the quadratic sieve when the number is very large.
    #'
    #' For example, using the defaults on `mostWanted1983` above only adds a few seconds.
    #'

    primeFactorizeBig(mostWanted1983, showStats = TRUE, nThreads = 8)

    #'
    #' ### The Lentra's power
    #'
    #' Performing a few iterations of the Pollard's rho algorithm followed by the quadratic sieve can prove to be a terrible solution when the input size is incredibly large. More than likely, your constrained Pollard's rho algorithm won't find any prime factors, and eventually will pass an enormous composite number to the quadratic sieve. The quadratic sieve isn't optimized for finding _smallish_ factors in large composites. It will treat the input similarly to a semiprime of the same size. For example, the number `ecmExpo1` below is 428 digits and if passed to the quadratic sieve algorithm, it could take years to factor as not even the most sophisticated semiprime algorithms can easily factor numbers of this size (See https://en.wikipedia.org/wiki/RSA_Factoring_Challenge). However, with the ECM, this is light work.
    #'

    ecmExpo1 <- pow.bigz(prod(nextprime(urand.bigz(10, 47, 13))), 3) * nextprime(urand.bigz(1, 47, 123))

    system.time(primeFactorizeBig(ecmExpo1, skipPolRho = TRUE, nThreads = 8))

    ecmExpo2 <- pow.bigz(prod(nextprime(urand.bigz(10, 40, 42))), 3) * pow.bigz(prod(nextprime(urand.bigz(10, 45, 42))), 5) * nextprime(urand.bigz(1, 80, 123)) * nextprime(urand.bigz(1, 90, 123))

    ecm2 <- primeFactorizeBig(ecmExpo2, nThreads = 8, showStats = TRUE)

    head(ecm2)

    tail(ecm2)

    length(ecm2)

    #'
    #' ### Safely Interrupt Execution
    #'
    #' If you want to interrupt a command which will take a long time, hit Ctrl + c, or esc if using RStudio, to stop execution. When you utilize multiple threads with a very large number (e.g. 90 digit semiprime), you will be able to interrupt execution once every ~30 seconds.
    #'

    ## User hits Ctrl + c
    ## system.time(quadraticSieve(prod(nextprime(urand.bigz(2, 100, 42)))))
    ## Seed initialisation
    ## Error in PrimeFactorization(n, FALSE, showStats, TRUE, TRUE, nThreads,  :
    ##   C++ call interrupted by the user.
    ## Timing stopped at: 2.164 0.039 2.167

    #'
    #' ## Complete Factorization with **`divisorsBig`**
    #'
    #' This function generates the complete factorization for many (possibly large) numbers. It is vectorized and can also return a named list.
    #'

    ## Get all divisors of a given number:
    divisorsBig(1000)

    ## Or, get all divisors of a vector:
    divisorsBig(urand.bigz(nb = 2, size = 100, seed = 42), namedList = TRUE)

    #'
    #' ### Efficiency
    #'
    #' It is very efficient as well. It is equipped with a modified merge sort algorithm that significantly outperforms the `std::sort`/`bigvec` (the class utilized in the `R gmp` package) combination.
    #'

    hugeNumber <- pow.bigz(2, 100) * pow.bigz(3, 100) * pow.bigz(5, 100)
    system.time(overOneMillion <- divisorsBig(hugeNumber))

    length(overOneMillion)

    ## Output is in ascending order
    tail(overOneMillion)

    #'
    #' ### Correct Ordering
    #'
    #' Another benefit is that it will return correct orderings on extremely large numbers when compared to sorting large vectors in `base R`.  Typically in `base R` you must execute the following: `order(asNumeric(myVectorHere))`. When the numbers get large enough, precision is lost which leads to incorrect orderings. Observe:
    #'

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

    ## N.B. The first and second elements are incorrect order (among others)
    head(testBaseSort)

    #'
    #' ## Acknowledgments and Resources
    #'
    #' * Credit to [primo](<https://codegolf.stackexchange.com/users/4098/primo>) (Mike Tryczak) and his excellent answer to [Fastest semiprime factorization](<https://codegolf.stackexchange.com/a/9088/52987>).
    #'
    #' * [Factoring large numbers with quadratic sieve](<https://blogs.msdn.microsoft.com/devdev/2006/06/19/factoring-large-numbers-with-quadratic-sieve/>) on MSDN Archive.
    #'
    #' * A really nice concise example is given here: [Factorization of _n = 87463_ with the Quadratic Sieve](<https://www.math.colostate.edu/~hulpke/lectures/m400c/quadsievex.pdf>)
    #'
    #' * [Smooth numbers and the quadratic sieve](<http://library.msri.org/books/Book44/files/03carl.pdf>) by Carl Pomerance
    #'
    #' * [Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware](<https://www.iacr.org/archive/ches2006/10/10.pdf>) by Gaj K. et al.
    #'
    #' * [Integer Factorization using the Quadratic Sieve](<http://micsymposium.org/mics_2011_proceedings/mics2011_submission_28.pdf>) by Chad Seibert
    #'
    #' * In the stackoverflow question and answer [What is the most efficient factoring algorithm for quadratic sieve extraction phase?](https://stackoverflow.com/q/63541365/4408538) by  [Ilya Gazman](https://github.com/gazman-sdk), an efficient method for checking divisibility is sketched out that utilizes built-in types. You can see more on a video Ilya put on youtube: [E15: Quadratic Sieve Running on Java - Receiving](https://youtu.be/sXg_WrCUX-Q). While `mpz_divisible_ui_p` is very efficient, we found better performance using this method.
    #'
    #' * [R Function for returning ALL factors](<https://stackoverflow.com/a/49742904/4408538>)
    #'
    #' * [Issues factoring large number that is 99 digits long](<https://stackoverflow.com/a/66128627/4408538>)
    #'
    #' ## Current Research
    #'
    #' Currently, our main focus for version `1.2.0` will be implementing the self initializing quadratic sieve.
    #'
    #' ## Contact
    #'
    #' If you would like to report a bug, have a question, or have suggestions for possible improvements, please file an [issue](<https://github.com/jwood000/RcppBigIntAlgos/issues>).
    #'

}, advertise = FALSE, venue = "r", html_preview = FALSE, wd = ".")
