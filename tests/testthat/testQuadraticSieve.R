context("testing quadraticSieve")

test_that("quadraticSieve generates correct numbers", {
    options(scipen = 999)
    expect_equal(asNumeric(quadraticSieve(100)), c(2, 2, 5, 5))
    expect_equal(asNumeric(quadraticSieve(-100)), c(-1, 2, 2, 5, 5))
    expect_equal(asNumeric(quadraticSieve(2)), 2)
    expect_equal(asNumeric(quadraticSieve(1)), numeric(0))
    expect_equal(asNumeric(quadraticSieve(1000)), c(2,2,2,5,5,5))
    
    hugeNumSmallPrimes = prod.bigz(rep(c(2, 3, 5, 7, 11, 13), times = 95:100))
    expect_equal(quadraticSieve(hugeNumSmallPrimes),
                 factorize(hugeNumSmallPrimes))
    
    ## Test perfect powers
    expect_equal(quadraticSieve(pow.bigz(nextprime("1234567890987654321"), 10)),
                 rep(nextprime("1234567890987654321"), 10))
    
    ## Test semi-primes
    testNums <- lapply(42:65, function(x) {
                            prod(nextprime(urand.bigz(2, x, x)))
                        })
    
    gmpFactorize <- lapply(42:65, function(x) {
                            facs <- nextprime(urand.bigz(2, x, x))
                            facs[order(asNumeric(facs))]
                        })
    
    quadSieveFacs <- lapply(testNums, quadraticSieve)
    expect_equal(gmpFactorize, quadSieveFacs);
    
    prime500 <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                  59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                  127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
                  191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
                  257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
                  331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
                  401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
                  467, 479, 487, 491, 499)
    
    ## ******************* Gmpxx Note **********************
    ## We are now using gmpxx, thus we are able to take advantage
    ## of std::vector<mpz_class> which means we no longer have to
    ## worry with the complex resizing issues we had before. We
    ## now only include these examples for posterity
    
    ## We are ensuring that the resizing of the factors array is correct
    expect_equal(quadraticSieve(prod.bigz(prime500)),
                 factorize(prod.bigz(prime500)))
    
    prime5K58 <- c(5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081,
                   5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171,
                   5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273,
                   5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381,
                   5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441,
                   5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519,
                   5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623,
                   5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689,
                   5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783,
                   5791)
    
    ## See Gmpxx Note above
    ## Again, we are ensuring that the resizing of the factors array is
    ## correct (It can occur in multiple places in quadraticSieve)
    a <- quadraticSieve(prod.bigz(c(prime500, prime5K58)))
    b <- factorize(prod.bigz(c(prime500, prime5K58)))
    
    ## factorize does not sort the output when the prime numbers are greater
    ## than 4001. This is actually when the pollard rho algorithm is engaged
    b <- b[order(asNumeric(b))]
    expect_equal(a, b)
})

test_that("quadraticSieve generates correct numbers with multiple threads", {
    test1 <- nextprime(urand.bigz(2, 82, 42))
    expect_equal(prod(quadraticSieve(prod(test1), nThreads = 2)), prod(test1))
    
    test2 <- nextprime(urand.bigz(2, 83, 42))
    expect_equal(prod(quadraticSieve(prod(test2), nThreads = 2)), prod(test2))
    
    test3 <- nextprime(urand.bigz(2, 84, 42))
    expect_equal(prod(quadraticSieve(prod(test3), nThreads = 2)), prod(test3))
    
    test4 <- nextprime(urand.bigz(2, 85, 42))
    expect_equal(prod(quadraticSieve(prod(test4), nThreads = 2)), prod(test4))
    
    test5 <- nextprime(urand.bigz(2, 86, 42))
    expect_equal(prod(quadraticSieve(prod(test5), nThreads = 2)), prod(test5))
    
    test6 <- nextprime(urand.bigz(2, 90, 42))
    expect_equal(prod(quadraticSieve(prod(test6), nThreads = 2)), prod(test6))
})

test_that("quadraticSieve produces appropriate error messages", {
    expect_error(quadraticSieve(1:10), "Can only factor one number at a time")
    expect_error(quadraticSieve(0), "Cannot factorize 0")
    expect_error(quadraticSieve(1234567, nThreads = "9"),
                 "This type is not supported! No conversion possible for nThreads")
    expect_error(quadraticSieve(1234567, nThreads = 3.5),
                 "nThreads must be a whole number")
    expect_error(quadraticSieve(1234567, showStats = "T"),
                 "Only logical values are supported for showStats")
    expect_error(quadraticSieve(1234567, skipExtPolRho = "T"),
                 "Only logical values are supported for skipExtPolRho")
})