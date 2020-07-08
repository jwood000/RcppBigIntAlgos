#ifndef SIEVE_UTILS_R
#define SIEVE_UTILS_R

#include "TonelliShanks.h"
#include <Rcpp.h>
#include <gmp.h>

constexpr std::size_t hundredThousand = 100000;
constexpr double Significand53 = 9007199254740991.0;

using vec2dint = std::vector<std::vector<int>>;
using hash64vec = std::unordered_map<std::uint64_t, std::vector<int>>;
using hash64mpz_t = std::unordered_map<std::uint64_t, mpz_t>;
using hash64size_t = std::unordered_map<std::uint64_t, std::size_t>;

std::vector<std::size_t> setSieveDist(mpz_t myNum, mpz_t *const TS,
                                      const std::vector<std::size_t> &facBase);

std::vector<std::size_t> outersect(std::vector<std::size_t> &x,
                                   std::vector<std::size_t> &y);

std::vector<std::uint8_t> myIntToBit(std::size_t x, std::size_t dig);

std::vector<std::size_t> getPrimesQuadRes(mpz_t myN, double LimB, double fudge1,
                                          double sqrLogLog, std::size_t myTarget);

void SinglePoly(std::vector<std::size_t> &polySieveD, mpz_t *const smoothInterval,
                const std::vector<std::size_t> &SieveDist, mpz_t *const TS,
                const std::vector<std::size_t> &facBase, mpz_t *const mpzFacBase,
                const std::vector<double> &LnFB, mpz_t *const largeCoFactors,
                std::vector<std::size_t> &myStart, mpz_t *const partialInterval,
                vec2dint &powsOfSmooths, vec2dint &powsOfPartials,
                hash64vec &partFactorsMap, hash64mpz_t &partIntvlMap,
                hash64size_t &keepingTrack, std::vector<std::size_t> &coFactorIndexVec,
                std::size_t &nPartial, std::size_t &nSmooth, std::size_t &coFactorInd,
                mpz_t intVal, mpz_t myNum, mpz_t Atemp, mpz_t Btemp, mpz_t temp,
                mpz_t A, mpz_t B, mpz_t C, mpz_t Atemp2, mpz_t lowBound,
                double theCut, std::size_t DoubleLenB, std::size_t mpzFacSize,
                std::size_t vecMaxSize, std::size_t strt);

#endif
