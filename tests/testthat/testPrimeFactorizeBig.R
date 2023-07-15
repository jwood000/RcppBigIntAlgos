context("testing primeFactorizeBig")

test_that("primeFactorizeBig generates correct numbers", {
    options(scipen = 999)
    expect_equal(length(primeFactorizeBig(1:100)), 100)
    expect_equal(primeFactorizeBig(2), 2)
    expect_equal(length(primeFactorizeBig(1)), 0)
    expect_equal(primeFactorizeBig(-10), c(-1, 2, 5))
    expect_equal(primeFactorizeBig(1000), c(2, 2, 2, 5, 5, 5))
    expect_equal(primeFactorizeBig(1000, TRUE), c(2, 2, 2, 5, 5, 5))

    n5 <- gmp::as.bigz("94968915845307373740134800567566911")
    out <- c(gmp::as.bigz("216366620575959221"),
             gmp::as.bigz("438925910071081891"))
    expect_equal(primeFactorizeBig(n5), out)

    ## Test Names
    expect_equal(
        as.integer(names(primeFactorizeBig(100, namedList = TRUE))),
        integer(0)
    )
    expect_equal(
        as.numeric(names(primeFactorizeBig((10^12):(10^12 + 100),
                                           namedList = TRUE))),
        (10^12):(10^12 + 100)
    )
})

test_that("primeFactorizeBig generates correct numbers with elliptic curve", {
    capture.output(test1 <- nextprime(urand.bigz(10, 40, 123)),
                   file = nullfile())
    expect_equal(primeFactorizeBig(prod(test1), nThreads = 2,
                                   skipECM = FALSE, skipPolRho = TRUE),
                 test1[order(asNumeric(test1))])

    capture.output(test2 <- nextprime(urand.bigz(15, 40, 321)),
                   file = nullfile())
    expect_equal(primeFactorizeBig(prod(test2), nThreads = 2,
                                   skipECM = FALSE, skipPolRho = TRUE),
                 test2[order(asNumeric(test2))])
})

test_that("primeFactorizeBig generates correct numbers with multiple threads", {
    capture.output(test <- nextprime(urand.bigz(2, 82, 8191)),
                   file = nullfile())
    expect_equal(primeFactorizeBig(prod(test), nThreads = 2,
                                   skipECM = TRUE, skipPolRho = TRUE),
                 test[order(asNumeric(test))])
})

test_that("primeFactorizeBig produces appropriate error messages", {
    expect_error(primeFactorizeBig(0), "Cannot factorize 0")
    expect_error(primeFactorizeBig(1234567, skipPolRho = "T"),
                 "Only logical values are supported for skipPolRho")
    expect_error(primeFactorizeBig(1234567, skipECM = "T"),
                 "Only logical values are supported for skipECM")
})
