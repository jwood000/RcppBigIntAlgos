#ifndef SOLUTION_SEARCH_H
#define SOLUTION_SEARCH_H

#include "ReduceMatrix.h"

void SolutionSearch(const std::vector<std::uint8_t> &mat, std::size_t matNRows,
                    std::size_t matNCols, const mpz_class &myNum,
                    const std::vector<mpz_class> &mpzFacBase,
                    const std::vector<mpz_class> &testInterval,
                    std::vector<mpz_class> &factors,
                    std::size_t nThreads, bool bShowStats);

#endif
