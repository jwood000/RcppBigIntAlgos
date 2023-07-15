#ifndef QUADRATICSIEVE_GMP_H
#define QUADRATICSIEVE_GMP_H

#include "CppConvert.h"

// Factoring with Multiple Polynomial Quadratic Sieve.
// 
// In addition to the references in the man file, the links below are very helpful:
// - 1:
//         URL: http://www.cs.virginia.edu/crab/QFS_Simple.pdf
//         author: Eric Landquist
//         date: December 14, 2001
//         title: The Quadratic Sieve Factoring Algorithm
// - 2:
//         URL: https://blogs.msdn.microsoft.com/devdev/2006/06/19/factoring-large-numbers-with-quadratic-sieve/
//         author: MSDN Archive
//         date: June 19, 2006
//         title: Factoring large numbers with quadratic sieve
// - 3:
//        URL: http://www.math.colostate.edu/~hulpke/lectures/m400c/quadsievex.pdf
//
// Factor large number into two smaller numbers (possibly prime)
// t: number to factorize
// fudge1: number that alters the size of the prime factor base
// fudge2: number that adjusts the lower bound for the log sum
// LenB: Size of the sieving interval
// factors [out]: the list of factors
// If fudge1, fudge2, or LenB are left blank, they will be determined
// in the algorithm based off of the literature regarding the QS

void QuadraticSieve(const mpz_class &myNum, std::vector<mpz_class> &factors,
                    std::size_t nThreads, bool bShowStats);

#endif
