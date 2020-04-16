#ifndef SIEVE_UTILS_R
#define SIEVE_UTILS_R

#include "TonelliShanks.h"
#include <Rcpp.h>
#include <gmp.h>

constexpr std::size_t hundredThousand = 100000;
constexpr double Significand53 = 9007199254740991.0;

std::vector<std::size_t> setSieveDist(mpz_t myNum, mpz_t *const TS,
                                      const std::vector<std::size_t> &facBase,
                                      std::size_t facSize);

std::vector<std::size_t> outersect(std::vector<std::size_t> &x,
                                   std::vector<std::size_t> &y);

uint64_t makeKey(mpz_t x);

std::vector<uint8_t> myIntToBit(std::size_t x, std::size_t dig);

std::vector<std::size_t> getPrimesQuadRes(mpz_t myN, double LimB, double fudge1,
                                          double sqrLogLog, std::size_t myTarget);

void sieveLists(std::size_t facSize, const std::vector<std::size_t> &FBase,
                std::size_t LenB2, mpz_t *const sqrDiff,
                const std::vector<double> &LnFB,
                std::vector<double> &myLogs,
                std::size_t minPrime,
                const std::vector<std::size_t> &polySieveD,
                mpz_t lowerBound);

#endif
