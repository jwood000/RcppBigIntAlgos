# algos

Overview
---------
A collection of optimized functions implemented in C++, Rcpp, and the C library GMP (GNU Multiple Precision Arithmetic) for solving problems in combinatorics and computational mathematics.

Usage
-----
``` r
## Find all 3-tuples without repetition of the numbers c(1, 2, 3, 4).
comboGeneral(4, 3)
     [,1] [,2] [,3]
[1,]   1    2    3
[2,]   1    2    4
[3,]   1    3    4
[4,]   2    3    4


## Find all 3-tuples without repetition of c(2, 3, 5, 7, 11), such that the product is less than 130.
comboGeneral(5, 3, v = c(2, 3, 5, 7, 11), 
             constraintFun = "prod", 
             comparisonFun = "<", 
             limitConstraints = 130)
     [,1] [,2] [,3]
[1,]    2    3    5
[2,]    2    3    7
[3,]    2    3   11
[4,]    2    5    7
[5,]    2    5   11
[6,]    3    5    7


## Find all 3-tuples with repetition of c(2, 3, 5, 7, 11), such that the product is greater than 170
comboGeneral(5, 3, v = c(2, 3, 5, 7, 11), 
             constraintFun = "prod", 
             comparisonFun = ">",
             repetition = TRUE,
             limitConstraints = 170)
      [,1] [,2] [,3]
 [1,]   11   11   11
 [2,]   11   11    7
 [3,]   11   11    5
 [4,]   11   11    3
 [5,]   11   11    2
 [6,]   11    7    7
 [7,]   11    7    5
 [8,]   11    7    3
 [9,]   11    5    5
[10,]    7    7    7
[11,]    7    7    5
[12,]    7    5    5


## get all divisors of a given integer.
getDivisors(1000)
Big Integer ('bigz') object of length 16:
 [1] 1    2    4    5    8    10   20   25   40   50   100  125  200  250  500  1000
 
 
 ## get prime factorization for every number from 1 to n
 primeFactorizationList(5)
[[1]]
integer(0)

[[2]]
[1] 2

[[3]]
[1] 3

[[4]]
[1] 2 2

[[5]]
[1] 5
```
