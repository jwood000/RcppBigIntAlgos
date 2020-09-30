#ifndef SOLUTION_UTILS_H
#define SOLUTION_UTILS_H

#include <vector>
#include <numeric>
#include <cmath>
#include "GmpxxCopy.h"
#include <bitset>

constexpr unsigned long int oneThousand = 1000u;
constexpr std::uint64_t one64 = 1;

std::vector<std::uint8_t> MyIntToBit(std::size_t x, std::size_t dig);
constexpr size_t wordSize = sizeof(std::bitset<1>) * 8;

#endif
