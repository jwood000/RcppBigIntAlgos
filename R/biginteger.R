divisorsBig <- function(v, namedList = FALSE) .Call("getDivisorsC", v, namedList, PACKAGE = "bigIntegerAlgos")

quadraticSieve <- function(n) .Call("QuadraticSieveContainer", n, PACKAGE = "bigIntegerAlgos")

quadResSquareRoot <- function(n, p) .Call("QuadraticResidueContainer", n, p, PACKAGE = "bigIntegerAlgos")
