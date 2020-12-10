#ifndef LENSTRA_ECM_H
#define LENSTRA_ECM_H

#include "ImportExportMPZ.h"
#include <limits>
#include <map>

constexpr unsigned long int MaxCurveBound = std::numeric_limits<unsigned long int>::max();

static const std::map<unsigned long int, int> CurveLookup = {{50, 8}, {55, 9},
                                                             {60, 10}, {65, 11},
                                                             {70, 12}, {77, 13},
                                                             {84, 14}, {91, 15},
                                                             {100, 16}, 
                                                             {MaxCurveBound, 17}};

unsigned long int GetMaxCurves(std::size_t maxLoopIter);
std::vector<unsigned long int> GenerateNPrimes(std::size_t limit);

bool LenstraECM(const mpz_class &n, std::size_t maxLoopIter,
                const std::vector<unsigned long int> &primes,
                std::vector<mpz_class> &factors,
                std::size_t &numCurves, std::size_t nThreads);
  
#endif
