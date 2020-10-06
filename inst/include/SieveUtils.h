#ifndef SIEVE_UTILS_R
#define SIEVE_UTILS_R

#include "TonelliShanks.h"
#include <unordered_map>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <cmath>

constexpr double Significand53 = 9007199254740991.0;
constexpr int L1Cache = 32768;

using vec2dint = std::vector<std::vector<int>>;
using hash64vec = std::unordered_map<std::uint64_t, std::vector<int>>;
using hash64mpz = std::unordered_map<std::uint64_t, mpz_class>;

std::vector<std::size_t> SetSieveDist(const std::vector<int> &facBase,
                                      const mpz_class &myNum);

std::vector<int> GetPrimesQuadRes(const mpz_class &myN, double LimB, double fudge1,
                                  double sqrLogLog, std::size_t myTarget);

void SinglePoly(const std::vector<std::size_t> &SieveDist,
                const std::vector<int> &facBase, const std::vector<int> &LnFB,
                vec2dint &powsOfSmooths, vec2dint &powsOfPartials,
                std::vector<int> &myStart, hash64vec &partFactorsMap,
                hash64mpz &partIntvlMap, std::vector<mpz_class> &smoothInterval,
                std::vector<uint64_t> &largeCoFactors,
                std::vector<mpz_class> &partialInterval,
                const mpz_class &NextPrime, const mpz_class &LowBound,
                const mpz_class &myNum, int theCut,int DoubleLenB,
                int mpzFacSize, int vecMaxSize, std::size_t strt);

#endif
