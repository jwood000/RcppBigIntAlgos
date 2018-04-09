context("testing quadraticSieve")

test_that("quadraticSieve generates correct numbers", {
    options(scipen = 999)
    expect_equal(asNumeric(quadraticSieve(100)), c(2, 2, 5, 5))
    expect_equal(asNumeric(quadraticSieve(2)), 2)
    expect_equal(asNumeric(quadraticSieve(1)), numeric(0))
    expect_equal(asNumeric(quadraticSieve(1000)), c(2,2,2,5,5,5))
    
    ## Test semi-primes
    testNums <- lapply(1:5, function(x) {
                            prod(nextprime(urand.bigz(2,42,x)))
                        })
    
    gmpFactorize <- lapply(testNums, function(x) {
                            facs <- factorize(x)
                            facs[order(asNumeric(facs))]
                        })
    
    quadSieveFacs <- lapply(testNums, quadraticSieve)
    
    expect_equal(gmpFactorize, quadSieveFacs);
})

test_that("quadraticSieve produces appropriate error messages", {
    expect_error(quadraticSieve(1:10), "Can only factor one number at a time")
    expect_error(quadraticSieve(-2), "Can only factor positive numbers")
    expect_error(quadraticSieve(0), "Can only factor positive numbers")
})