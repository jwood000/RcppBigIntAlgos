#ifndef SOLUTION_SEARCH_H
#define SOLUTION_SEARCH_H

#include "SieveUtils.h"
#include "ReduceMatrix.h"

void SolutionSearch(const std::vector<std::uint8_t> &mat, std::size_t matNRows,
                    std::size_t matNCols, mpz_t n, const std::vector<mpz_class> &mpzFacBase,
                    const std::vector<mpz_class> &testInterval, mpz_t *const factors,
                    std::size_t nThreads);

#endif
