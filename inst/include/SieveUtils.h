#ifndef SIEVE_UTILS_H
#define SIEVE_UTILS_H

#include "TonelliShanks.h"
#include <unordered_map>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <numeric>
#include <cmath>

constexpr double dblLogMult = 100.0;
constexpr double Significand53 = 9007199254740991.0;
constexpr int L1Cache = 32768;

using logType = std::uint16_t;
using vec2dint = std::vector<std::vector<int>>;
using hash64vec = std::unordered_map<std::uint64_t, std::vector<int>>;
using hash64mpz = std::unordered_map<std::uint64_t, mpz_class>;

int int_invert(unsigned int n, unsigned int p);

std::vector<std::size_t> GetSieveDist(const std::vector<int> &facBase,
                                      const mpz_class &myNum);

std::vector<int> GetPrimesQuadRes(const mpz_class &myN, double LimB, double fudge1,
                                  double sqrLogLog, std::size_t myTarget);

#endif
