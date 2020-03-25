#ifndef SIEVE_UTILS_R
#define SIEVE_UTILS_R

#include <algorithm>
#include <Algos.h>
#include <iostream>

typedef std::vector<int64_t> v1d;
typedef std::vector<v1d> v2d;
typedef std::vector<v2d> v3d;

constexpr std::size_t hundredThousand = 100000;
constexpr double Significand53 = 9007199254740991.0;

void setSieveDist(mpz_t myNum, const v1d &facBase,
                  std::size_t facSize, v2d &SieveDist);

std::vector<std::size_t> outersect(std::vector<std::size_t> &x,
                                   std::vector<std::size_t> &y);

uint64_t makeKey(mpz_t x);

std::vector<uint8_t> myIntToBit(std::size_t x, std::size_t dig);

v1d getPrimesQuadRes(mpz_t myN, double n);

v3d sieveLists(int64_t facLim, const v1d &FBase, int64_t vecLen, mpz_t *const sqrDiff);

#endif
