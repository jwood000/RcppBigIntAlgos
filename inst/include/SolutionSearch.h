#ifndef SOLUTION_SEARCH_H
#define SOLUTION_SEARCH_H

#include "SieveUtils.h"
#include "ReduceMatrix.h"

void SolutionSearch(const std::vector<std::uint8_t> &mat, std::size_t matNRows,
                    std::size_t matNCols, mpz_t n, mpz_t *const mpzFacBase,
                    mpz_t *const test, mpz_t *const factors,
                    std::size_t nThreads);

#endif
