#ifndef IMPORTEXPORTMPZ_GMP_R
#define IMPORTEXPORTMPZ_GMP_R

#include <Rcpp.h>
#include <gmp.h>

constexpr std::size_t intSize = sizeof(int);
constexpr std::size_t numb = 8 * intSize;

/**
 * Functions for importing/exporting and converting SEXPs to mpz_t
 * 
 */
void createMPZArray(SEXP v, mpz_t myVec[], std::size_t sizevec);

int myRaw(char* raw, mpz_t value, std::size_t totals);

void quickSort(mpz_t arr[], int left, int right, std::vector<std::size_t>& lens);

#endif
