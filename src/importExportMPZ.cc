/* This file contains utility functions that
 * are used for converting to and from type mpz_t,
 * as well as sorting an array of type mpz_t.
 * 
 * createMPZArray and myRaw are slightly modified versions
 * of "bigvec create_vector(const SEXP & param)" and 
 * "int biginteger::as_raw(char* raw) const", respectively,
 * from the source files bigintegerR.cc/ biginteger.cc from
 * the R gmp package.
 * 
 * The quickSort function is based off of the quicksort
 * algorithm in C++ found here:
 *      http://www.algolist.net/Algorithms/Sorting/Quicksort
 */

#include "ImportExportMPZ.h"
#include <algorithm>

void createMPZArray(SEXP v, mpz_t *const myVec, std::size_t sizevec) {
    
    switch (TYPEOF(v)) {
        case RAWSXP: {
            // deserialise the vector. first int is the size.
            const char* raw = (char*) RAW(v);
            int pos = intSize; // position in raw[]. Starting after header.
            
            for (std::size_t i = 0; i < sizevec; i++) {
                const int* r = (int*) (&raw[pos]);
                
                if (r[0] > 0) {
                    mpz_import(myVec[i], r[0], 1, intSize, 0, 0, (void*) &(r[2]));
                    
                    if(r[1] == -1)
                        mpz_neg(myVec[i], myVec[i]);
                } else {
                    mpz_set_si(myVec[i], 0);
                }
                
                pos += intSize * (2 + (mpz_sizeinbase(myVec[i], 2) + numb - 1) / numb);
            }
            
            break;
        }
        case REALSXP: {
            double* myDbl = REAL(v);
            
            for (std::size_t j = 0; j < sizevec; j++) {
                /// New:   numeric '+- Inf'  give  +- "Large" instead of NA
                const double dj = myDbl[j];
                
                if(R_FINITE(dj) || ISNAN(dj)) {
                    mpz_set_d(myVec[j], dj);
                } else { 
                    // dj is +- Inf : use LARGE ( = +- 2 ^ 80000 -- arbitrarily )
                    // FIXME: Keep 'LARGE' a static const; initialized only once
                    mpz_ui_pow_ui(myVec[j], 2u, 8000u);
                    
                    if (dj == R_NegInf)
                        mpz_neg(myVec[j], myVec[j]);
                }
            }
            
            break;
        }
        case INTSXP:
        case LGLSXP: {
            int* myInt = INTEGER(v);
            
            for (std::size_t j = 0; j < sizevec; j++)
                mpz_set_si(myVec[j], myInt[j]);
            
            break;
        }
        case STRSXP: {
            for (std::size_t i = 0; i < sizevec; i++) {
                if (STRING_ELT(v,i) == NA_STRING) {
                    mpz_set_si(myVec[i], 0);
                } else {
                    mpz_set_str(myVec[i], CHAR(STRING_ELT(v,i)), 10);
                }
            }
            
            break;
        }
        default:
            // no longer: can be fatal later! return bigvec();
            Rf_error(_("only logical, numeric or character (atomic) vectors can be coerced to 'bigz'"));
    }
}

int myRaw(char* raw, mpz_t value, std::size_t totals) {
    memset(raw, 0, totals);
    
    int* r = (int*) raw;
    r[0] = totals / intSize - 2;
    
    r[1] = static_cast<int>(mpz_sgn(value));
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);
    
    return totals;
}

void quickSort(mpz_t *const arr, int left, int right,
               std::vector<std::size_t> &lens) {
    
    int i = left, j = right, mid;
    mpz_t pivot;
    mpz_init(pivot);
    
    mid = (left + right) / 2;
    mpz_set(pivot, arr[mid]);
    
    /* partition */
    while (i <= j) {
        while (mpz_cmp(arr[i], pivot) < 0)
            ++i;
        
        while (j >= 0 && mpz_cmp(arr[j], pivot) > 0)
            --j;
        
        if (i <= j) {
            mpz_swap(arr[i], arr[j]);
            std::swap(lens[i], lens[j]);
            ++i;
            --j;
        }
    }
    
    mpz_clear(pivot);
    
    /* recursion */
    if (left < j)
        quickSort(arr, left, j, lens);
    
    if (i < right)
        quickSort(arr, i, right, lens);
}
