/*! 
 *  \file rsafactorize.cc
 *  \brief C function that transfers input from R to 
 *          quadraticSieve function for factoring large
 *            numbers and returning result to R console
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:33:33 EDT jwood000>
 *
 *  \author Joseph Wood
 *
 *  \note Licence: GPL (>=) 2  
 */

#include "quadraticsieve.h"
#include "rsafactorize.h"
#include "importExportMPZ.h"
#include <vector>
#include "Rgmp.h"

SEXP QuadraticSieveContainer (SEXP Rn) {
    
    unsigned int vSize;
    
    switch (TYPEOF(Rn)) {
        case RAWSXP: {
            const char* raw = (char*)RAW(Rn);
            vSize = ((int*)raw)[0];
            break;
        }
        default:
            vSize = LENGTH(Rn);
    }
    
    mpz_t myVec[1];
    mpz_init(myVec[0]);
    
    if (vSize > 1)
        error(_("Can only factor one number at a time"));
    
    createMPZArray(Rn, myVec, 1);
    mpz_t nmpz;
    mpz_init_set(nmpz, myVec[0]);
    mpz_t *result;
    result = (mpz_t *) malloc(2 * sizeof(mpz_t));
    mpz_init(result[0]); mpz_init(result[1]);
    
    // factor_using_division(nmpz, result);
    
    quadraticSieve (nmpz, 0.0, 0.0, 0, result);

    unsigned long int intSize = sizeof(int); // starting with vector-size-header
    unsigned long int numb = 8 * intSize;
    unsigned long int tempSize, size = intSize;
    std::vector<unsigned long int> mySizes(2);

    for (int i = 0; i < 2; i++) { // adding each bigint's needed size
        tempSize = intSize * (2 + (mpz_sizeinbase(result[i], 2) + numb - 1) / numb);
        size += tempSize;
        mySizes[i] = tempSize;
    }

    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size));
    char* r = (char*)(RAW(ans));
    ((int*)(r))[0] = 2; // first int is vector-size-header

    // current position in rPos[] (starting after vector-size-header)
    unsigned long int pos = intSize;
    for (int i = 0; i < 2; i++)
        pos += myRaw(&r[pos], result[i], mySizes[i]);

    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
    UNPROTECT(1);
    return(ans);
}