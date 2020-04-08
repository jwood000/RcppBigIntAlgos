/*! 
 *  \file RSAFactorize.cc
 *  
 *  \version 2
 *
 *  \date Created: 2017-10-06
 *  \date Last modified: Time-stamp: <2020-03-25 12:30:00 EDT jwood000>
 *
 *  \author Joseph Wood
 *
 *  \note Licence: GPL (>=) 2  
 */

#include "RSAFactorUtils.h"
#include "CleanConvert.h"

// [[Rcpp::export]]
SEXP QuadraticSieveContainer(SEXP Rn, SEXP RShowStats,
                             SEXP RNumThreads, int maxThreads) {
    
    std::size_t vSize;
    
    switch (TYPEOF(Rn)) {
        case RAWSXP: {
            const char* raw = (char*) RAW(Rn);
            vSize = ((int*) raw)[0];
            break;
        }
        default:
            vSize = LENGTH(Rn);
    }
    
    mpz_t myVec[1];
    mpz_init(myVec[0]);
    
    if (vSize > 1)
        Rcpp::stop("Can only factor one number at a time");
    
    // This is from the importExportMPZ header
    CreateMPZArray(Rn, myVec, 1);
    mpz_t nmpz;
    mpz_init_set(nmpz, myVec[0]);
    mpz_clear(myVec[0]);
    
    if (mpz_sgn(nmpz) == 0)
        Rcpp::stop("Cannot factorize 0");
    
    const std::size_t IsNegative = (mpz_sgn(nmpz) < 0) ? 1 : 0;
    
    if (IsNegative)
        mpz_abs(nmpz, nmpz);

    std::size_t arrayMax = mpzChunkBig;
    std::vector<std::size_t> lengths;
    std::size_t numUni = 0;
    
    auto factors = FromCpp14::make_unique<mpz_t[]>(mpzChunkBig);
    
    for (std::size_t i = 0; i < mpzChunkBig; ++i)
        mpz_init(factors[i]);
    
    int nThreads = 1;
    const bool bShowStats = convertLogical(RShowStats, "bShowStats");
    
    if (!Rf_isNull(RNumThreads))
        convertInt(RNumThreads, nThreads, "nThreads");
    
    QuadSieveHelper(nmpz, factors, arrayMax,
                    numUni, lengths, nThreads, bShowStats);
    
    // Sort the prime factors as well as order the
    // lengths vector by the order of the factors array
    QuickSort(factors.get(), 0, numUni - 1, lengths);
    
    const std::size_t totalNum = std::accumulate(lengths.cbegin(),
                                                 lengths.cend(), 0u) + IsNegative;
    std::size_t size = intSize;
    std::vector<std::size_t> mySizes(totalNum);
    
    mpz_t negOne;
    mpz_init(negOne);
    mpz_set_si(negOne, -1);
    
    if (IsNegative) {
        const std::size_t tempSize = intSize * (2 + (mpz_sizeinbase(negOne, 2) + numb - 1) / numb);
        size += tempSize;
        mySizes[0] = tempSize;
    }
    
    for (std::size_t i = 0, count = IsNegative; i < numUni; ++i) { // adding each bigint's needed size
        for (std::size_t j = 0; j < lengths[i]; ++j, ++count) {
            const std::size_t tempSize = intSize * (2 + (mpz_sizeinbase(factors[i], 2) + numb - 1) / numb);
            size += tempSize;
            mySizes[count] = tempSize;
        }
    }
    
    Rcpp::RawVector ans(size);
    char* r = (char*) (RAW(ans));
    ((int*) (r))[0] = totalNum; // first int is vector-size-header
    
    // current position in pos[] (starting after vector-size-header)
    std::size_t pos = intSize;
    
    if (IsNegative)
        pos += myRaw(&r[pos], negOne, mySizes.front());
    
    for (std::size_t i = 0, count = IsNegative; i < numUni; ++i)
        for (std::size_t j = 0; j < lengths[i]; ++j, ++count)
            pos += myRaw(&r[pos], factors[i], mySizes[count]);
    
    for (std::size_t i = 0; i < arrayMax; ++i)
        mpz_clear(factors[i]);
    
    mpz_clear(negOne);
    factors.reset();
    ans.attr("class") = Rcpp::CharacterVector::create("bigz");
    return ans;
}
