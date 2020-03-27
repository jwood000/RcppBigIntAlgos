#ifndef SIEVE_UTILS_R
#define SIEVE_UTILS_R

#include <Rcpp.h>
#include <gmp.h>
#include <chrono>

// Used for checking whether user has interrupted computation
constexpr auto timeout = std::chrono::milliseconds(1000);

constexpr std::size_t hundredThousand = 100000;
constexpr double Significand53 = 9007199254740991.0;

void setSieveDist(mpz_t myNum, const std::vector<int64_t> &facBase,
                  std::size_t facSize, std::vector<std::size_t> &SieveDist);

std::vector<std::size_t> outersect(std::vector<std::size_t> &x,
                                   std::vector<std::size_t> &y);

uint64_t makeKey(mpz_t x);

std::vector<uint8_t> myIntToBit(std::size_t x, std::size_t dig);

std::vector<int64_t> getPrimesQuadRes(mpz_t myN, double n);

void sieveLists(std::size_t facSize, const std::vector<int64_t> &FBase,
                std::size_t LenB2, mpz_t *const sqrDiff,
                const std::vector<double> &LnFB,
                std::vector<double> &myLogs,
                std::vector<bool> &indexDiv,
                int64_t minPrime,
                const std::vector<std::size_t> &polySieveD,
                mpz_t lowerBound);

#endif
