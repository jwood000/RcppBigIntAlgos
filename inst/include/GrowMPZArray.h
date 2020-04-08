#ifndef GROW_MPZ_ARRAY_H
#define GROW_MPZ_ARRAY_H

#include <cstdlib>
#include <memory>
#include <gmp.h>

void Grow(std::unique_ptr<mpz_t[]> &myArray,
          std::size_t oldMax, std::size_t newMax);

#endif
