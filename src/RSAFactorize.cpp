/*! 
 *  \file RSAFactorize.cc
 *  
 *  \version 2
 *
 *  \date Created: 2017-10-06
 *  \date Last modified: Time-stamp: <2020-09-13 jwood000>
 *
 *  \author Joseph Wood
 *
 *  \note Licence: GPL (>=) 2  
 */

#include "RSAFactorUtils.h"
#include "CleanConvert.h"

// [[Rcpp::export]]
SEXP QuadraticSieveContainer(SEXP Rn, SEXP RShowStats, SEXP RNumThreads,
                             int maxThreads, SEXP RSkipExtPR) {
    
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
    
    if (vSize > 1)
        Rcpp::stop("Can only factor one number at a time");
    
    mpz_class nMpz;
    convertMpzClass(Rn, nMpz);
    
    if (sgn(nMpz) == 0)
        Rcpp::stop("Cannot factorize 0");
    
    const std::size_t IsNegative = (sgn(nMpz) < 0) ? 1 : 0;
    
    if (IsNegative)
        nMpz = abs(nMpz);
    
    if (cmp(nMpz, 1) == 0) {
        if (IsNegative) {
            mpz_class mpzNegOne = -1;
            Rcpp::RawVector myNegOne(intSize * 4);
            
            char* r = (char*) (RAW(myNegOne));
            ((int*) (r))[0] = 1;
            
            myRaw(&r[intSize], mpzNegOne.get_mpz_t(), intSize * 3);
            myNegOne.attr("class") = Rcpp::CharacterVector::create("bigz");
            return myNegOne;
        } else {
            Rcpp::RawVector noPrimeFacs(intSize);
            char* r = (char*) (RAW(noPrimeFacs));
            ((int*) (r))[0] = 0;
            noPrimeFacs.attr("class") = Rcpp::CharacterVector::create("bigz");
            return noPrimeFacs;
        }
    }
    
    std::vector<std::size_t> lengths;
    std::vector<mpz_class> factors;
    
    int nThreads = 1;
    const bool bShowStats = convertLogical(RShowStats, "showStats");
    const bool bSkipExtPR = convertLogical(RSkipExtPR, "skipExtPolRho");
    
    if (!Rf_isNull(RNumThreads))
        convertInt(RNumThreads, nThreads, "nThreads");
    
    QuadSieveHelper(nMpz, factors, lengths,
                    nThreads, bShowStats, bSkipExtPR);
    
    // Sort the prime factors as well as order the
    // lengths vector by the order of the factors array
    QuickSort(factors, 0, factors.size() - 1, lengths);
    
    const std::size_t totalNum = std::accumulate(lengths.cbegin(),
                                                 lengths.cend(), 0u) + IsNegative;
    std::size_t size = intSize;
    std::vector<std::size_t> mySizes(totalNum);
    
    mpz_class negOne = -1;
    
    if (IsNegative) {
        const std::size_t tempSize = intSize * (2 + (mpz_sizeinbase(negOne.get_mpz_t(), 2) + numb - 1) / numb);
        size += tempSize;
        mySizes[0] = tempSize;
    }
    
    for (std::size_t i = 0, count = IsNegative; i < factors.size(); ++i) { // adding each bigint's needed size
        for (std::size_t j = 0; j < lengths[i]; ++j, ++count) {
            const std::size_t tempSize = intSize * (2 + (mpz_sizeinbase(factors[i].get_mpz_t(), 2) + numb - 1) / numb);
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
        pos += myRaw(&r[pos], negOne.get_mpz_t(), mySizes.front());
    
    for (std::size_t i = 0, count = IsNegative; i < factors.size(); ++i)
        for (std::size_t j = 0; j < lengths[i]; ++j, ++count)
            pos += myRaw(&r[pos], factors[i].get_mpz_t(), mySizes[count]);
    
    ans.attr("class") = Rcpp::CharacterVector::create("bigz");
    return ans;
}
