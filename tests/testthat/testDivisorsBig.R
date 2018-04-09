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

test_that("divisorsBig produces appropriate error messages", {
    expect_error(divisorsBig(0), "Cannot factorize 0")
})
