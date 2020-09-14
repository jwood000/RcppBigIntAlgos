#ifndef IMPORTEXPORTMPZ_GMP_R
#define IMPORTEXPORTMPZ_GMP_R

#include <Rcpp.h>
#include "GmpxxCopy.h"

constexpr std::size_t intSize = sizeof(int);
constexpr std::size_t numb = 8 * intSize;

void CreateMPZVector(SEXP v, std::vector<mpz_class> &myVec, std::size_t sizevec);
void convertMpzClass(SEXP v, mpz_class &result);
int myRaw(char* raw, mpz_t value, std::size_t totals);

void QuickSort(std::vector<mpz_class> &arr, int left, int right,
               std::vector<std::size_t>& lens);

#endif
