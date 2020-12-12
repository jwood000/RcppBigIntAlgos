context("testing divisorsBig")

test_that("divisorsBig generates correct numbers", {
    options(scipen = 999)
    expect_equal(length(divisorsBig(1:100)), 100)
    expect_equal(divisorsBig(2), c(1, 2))
    expect_equal(divisorsBig(1), 1)
    expect_equal(divisorsBig(-10), c(-10, -5, -2, -1, 1, 2, 5, 10))
    expect_equal(divisorsBig(1000), c(1,2,4,5,8,10,20,
                                       25,40,50,100,125,
                                       200,250,500,1000))
    expect_equal(divisorsBig(1000, TRUE), c(1,2,4,5,8,10,20,
                                             25,40,50,100,125,
                                             200,250,500,1000))
    
    ## Test Names
    expect_equal(as.integer(names(divisorsBig(100, namedList = TRUE))), integer(0))
    expect_equal(as.numeric(names(divisorsBig((10^12):(10^12 + 100),
                                                 namedList = TRUE))), (10^12):(10^12 + 100))
})

test_that("divisorsBig generates correct numbers with multiple threads", {
    test1 <- nextprime(urand.bigz(2, 82, 8191))
    expect_equal(prod(divisorsBig(prod(test1), nThreads = 2,
                                  skipECM = TRUE, skipPolRho = TRUE)), prod(test1)^2)

    test2 <- nextprime(urand.bigz(2, 83, 8191))
    expect_equal(prod(divisorsBig(prod(test2), nThreads = 2,
                                  skipECM = TRUE, skipPolRho = TRUE)), prod(test2)^2)

    test3 <- nextprime(urand.bigz(2, 84, 8191))
    expect_equal(prod(divisorsBig(prod(test3), nThreads = 2,
                                  skipECM = TRUE, skipPolRho = TRUE)), prod(test3)^2)

    test4 <- nextprime(urand.bigz(2, 85, 8191))
    expect_equal(prod(divisorsBig(prod(test4), nThreads = 2,
                                  skipECM = TRUE, skipPolRho = TRUE)), prod(test4)^2)

    test5 <- nextprime(urand.bigz(2, 86, 8191))
    expect_equal(prod(divisorsBig(prod(test5), nThreads = 2,
                                  skipECM = TRUE, skipPolRho = TRUE)), prod(test5)^2)
})

test_that("divisorsBig produces appropriate error messages", {
    expect_error(divisorsBig(0), "Cannot factorize 0")
    expect_error(divisorsBig(1234567, skipPolRho = "T"),
                 "Only logical values are supported for skipPolRho")
    expect_error(divisorsBig(1234567, skipECM = "T"),
                 "Only logical values are supported for skipECM")
    expect_error(divisorsBig(1234567, nThreads = "9"),
                 "This type is not supported! No conversion possible for nThreads")
    expect_error(divisorsBig(1234567, nThreads = 3.5),
                 "nThreads must be a whole number")
})
