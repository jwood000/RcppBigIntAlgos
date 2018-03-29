#ifndef IMPORTEXPORTMPZ_GMP_R
#define IMPORTEXPORTMPZ_GMP_R 1

#include <vector>
#include "Rgmp.h"

typedef void (*gmp_binary)(mpz_t, const mpz_t, const mpz_t);

/**
 * Functions for importing/exporting and converting SEXPs to mpz_t
 * 
 */
void createMPZArray (SEXP v, mpz_t myVec[], unsigned int sizevec);

int myRaw (char* raw, mpz_t value, unsigned long int totals);

void quickSort(mpz_t arr[], int left, int right, std::vector<unsigned int>& lens);

#endif
