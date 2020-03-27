/*! 
 *  \file factorization.cc
 *  \brief C functions used for integer factorization.
 *      - myMergeSort is a modified merge sort algorithm. The classic merge
 *        sort creates an empty vector of the type you are sorting and fills
 *        in the ordered elements from two smaller sorted vectors. This
 *        would be expensive with type mpz_t. For this reason, we only keep
 *        track of the indices and never actually swap any of the indices in
 *        the mpz_t array. An integer vector of indices is returned and is
 *        used in writing out the output at the very bottom.  This method
 *        shows great efficiency gains over the naive method of using the
 *        class bigvec (from the R gmp package)/std::sort combination.
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2018-04-07 12:20:00 EDT jwood000>
 *
 *  \author Joseph Wood. Original C code from libgmp.
 *       See factor.cc from the R gmp package for more details.
 *
 *  \note Licence: GPL (>=) 2
 */

#include "PollardRho.h"
#include "ImportExportMPZ.h"
#include "Cpp14MakeUnique.h"

std::vector<std::size_t> myMergeSort(mpz_t *const arr, 
                                     std::vector<std::size_t> indPass,
                                     std::size_t numSecs, 
                                     std::size_t secSize) {
    
    std::size_t x, y;
    std::size_t count = 0;
    const std::size_t totalSize  = numSecs * secSize;
    
    std::size_t tempSize = totalSize;
    std::vector<std::size_t> leftOver, myInd(totalSize);
    
    for (std::size_t i = 0; i < secSize; ++i)
        myInd[i] = indPass[i];
    
    for (std::size_t i = secSize; i < totalSize; ++i)
        myInd[i] = i;
    
    leftOver.reserve(numSecs / 2);
    
    while (numSecs > 1) {
        if (numSecs % 2 != 0) {
            leftOver.push_back(tempSize);
            tempSize -= secSize;
        }
        
        std::vector<std::size_t> endPoints(numSecs, secSize);
        
        for (std::size_t i = 1; i < numSecs; ++i) {
            x = endPoints[i - 1] + secSize;
            endPoints[i] = x;
        }
        
        secSize *= 2;
        std::size_t k = numSecs / 2;
        
        std::size_t lim = (k < 2) ? k : 2;
        std::size_t left = x = 0;
        
        std::vector<std::size_t> tempInd, defaultInd(secSize, 0);
        
        for (std::size_t i = 0; i < lim; ++i) {
            count = 0;
            std::size_t twoI = 2 * i;
            y = endPoints[twoI];
            tempInd = defaultInd;
                     
            while (x < endPoints[twoI] && y < endPoints[twoI + 1]) {
                if (mpz_cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
                    tempInd[count] = myInd[x];
                    ++x;
                } else {
                    tempInd[count] = myInd[y];
                    ++y;
                }
                
                ++count;
            }
            
            for (std::size_t j = 0; j < count; ++j)
                myInd[left + j] = tempInd[j];
            
            x = left = endPoints[2 * i + 1];
        }
        
        if (k > 2) {
            for (std::size_t i = 2; i < k; ++i) {
                x = endPoints[2*i - 1];
                y = x - secSize;
                
                for (std::size_t j = 0; j < count; ++j)
                    myInd[x + j] = tempInd[j] + y;
            }
        }
        
        numSecs /= 2;
    }
    
    const int LOSize = leftOver.size();

    if (LOSize > 0) {
        for (int j = LOSize - 1; j >= 0; --j) {
            x = count = 0;
            y = tempSize;
            std::vector<std::size_t> tempInd(leftOver[j]);

            while (x < tempSize && y < leftOver[j]) {
                if (mpz_cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
                    tempInd[count] = myInd[x];
                    ++x;
                } else {
                    tempInd[count] = myInd[y];
                    ++y;
                }
                
                ++count;
            }
            
            for (std::size_t r = 0; r < count; ++r)
                myInd[r] = tempInd[r];
            
            tempSize = leftOver[j];
        }
    }

    return myInd;
}

SEXP factorNum(mpz_t val, mpz_t primeFacs[]) {
    
    if (mpz_cmp_ui(val, 1) == 0) {
        mpz_t mpzOne;
        mpz_init_set_si(mpzOne, 1);
        
        std::size_t oneSize = intSize * 3;
        std::size_t totalSize = intSize;
        
        totalSize += oneSize;
        Rcpp::RawVector myFacs(totalSize);
        
        char* r = (char*) (RAW(myFacs));
        ((int*) (r))[0] = 1;
        
        std::size_t pos = intSize;
        pos += myRaw(&r[pos], mpzOne, oneSize);
        
        myFacs.attr("class") = Rcpp::CharacterVector::create("bigz");
        mpz_clear(mpzOne);
        return myFacs;
    } else {
        const int sgn = mpz_sgn(val);
        
        std::vector<std::size_t> lengths;
        std::size_t numUni = 0;
        bool isNegative = false;
        
        if (sgn == 0)
            Rcpp::stop("Cannot factorize 0");
        
        if (sgn < 0) {
            mpz_abs(val, val);
            isNegative = true;
        }
        
        getPrimeFactors(val, primeFacs, numUni, lengths);
        quickSort(primeFacs, 0, numUni - 1, lengths);
        
        std::vector<std::size_t> myIndex(lengths[0] + 1);
        std::size_t ind, facSize = 1, numFacs = 1;
        
        for (std::size_t i = 0; i < numUni; ++i)
            numFacs *= (lengths[i] + 1);
        
        auto myMPZ = FromCpp14::make_unique<mpz_t[]>(numFacs);
        
        for (std::size_t i = 0; i < numFacs; ++i)
            mpz_init(myMPZ[i]);
        
        mpz_t temp, myPow;
        mpz_init(temp);
        mpz_init(myPow);
        
        for (std::size_t i = 0; i <= lengths[0]; ++i) {
            mpz_pow_ui(temp, primeFacs[0], i);
            mpz_set(myMPZ[i], temp);
            myIndex[i] = i;
        }
        
        if (numUni > 0) {
            for (std::size_t j = 1; j < numUni; ++j) {
                facSize *= (lengths[j - 1] + 1);
                
                for (std::size_t i = 1; i <= lengths[j]; ++i) {
                    ind = i*facSize;
                    mpz_pow_ui(myPow, primeFacs[j], i);
                    
                    for (std::size_t k = 0; k < facSize; ++k) {
                        mpz_mul(temp, myPow, myMPZ[myIndex[k]]);
                        mpz_set(myMPZ[ind + k], temp);
                    }
                }
                
                myIndex = myMergeSort(myMPZ.get(), myIndex, lengths[j] + 1, facSize);
            }
        }
        
        std::size_t tempSize, size = intSize;
        std::vector<std::size_t> mySizes(numFacs);
        
        for (std::size_t i = 0; i < numFacs; ++i) { // adding each bigint's needed size
            tempSize = intSize * (2 + (mpz_sizeinbase(myMPZ[i],2) + numb - 1) / numb);
            size += tempSize;
            mySizes[i] = tempSize;
        }
        
        if (!isNegative) {
            Rcpp::RawVector ansPos = Rcpp::no_init_vector(size);
            
            char* rPos = (char*) (RAW(ansPos));
            ((int*) (rPos))[0] = numFacs; // first int is vector-size-header
            
            // current position in rPos[] (starting after vector-size-header)
            std::size_t posPos = intSize;
            
            for (std::size_t i = 0; i < numFacs; ++i)
                posPos += myRaw(&rPos[posPos], myMPZ[myIndex[i]], mySizes[myIndex[i]]);
            
            ansPos.attr("class") = Rcpp::CharacterVector::create("bigz");
            return(ansPos);
        } else {
            size *= 2; // double size as every element will have a negative counterpart
            size -= intSize; // Remove superfluous initializing size
            Rcpp::RawVector ansNeg = Rcpp::no_init_vector(size);
            
            char* rNeg = (char*) (RAW(ansNeg));
            ((int*) (rNeg))[0] = 2 * numFacs; // first int is vector-size-header
            
            // current position in rNeg[] (starting after vector-size-header)
            std::size_t posNeg = intSize;
            
            // First write out negative numbers in reverse "myIndex" order
            for (int i = numFacs - 1; i >= 0; --i) {
                mpz_neg(temp, myMPZ[myIndex[i]]);
                posNeg += myRaw(&rNeg[posNeg], temp, mySizes[myIndex[i]]);
            }
            
            for (std::size_t i = 0; i < numFacs; ++i)
                posNeg += myRaw(&rNeg[posNeg], myMPZ[myIndex[i]], mySizes[myIndex[i]]);
            
            ansNeg.attr("class") = Rcpp::CharacterVector::create("bigz");
            return(ansNeg);
        }
    }
}

// [[Rcpp::export]]
SEXP getDivisorsC(SEXP Rv, SEXP RNamed) {
    
    std::size_t vSize = 0;
    
    switch (TYPEOF(Rv)) {
        case RAWSXP: {
            const char* raw = (char*)RAW(Rv);
            vSize = ((int*)raw)[0];
            break;
        }
        default:
            vSize = LENGTH(Rv);
    }
    
    auto myVec = FromCpp14::make_unique<mpz_t[]>(vSize);
    auto primeFacs = FromCpp14::make_unique<mpz_t[]>(mpzChunkBig);
    
    for (std::size_t i = 0; i < vSize; ++i)
        mpz_init(myVec[i]);
    
    createMPZArray(Rv, myVec.get(), vSize);
    
    for (std::size_t i = 0; i < mpzChunkBig; ++i)
        mpz_init(primeFacs[i]);
    
    if (vSize > 0) {
        if (vSize == 1) {
            return factorNum(myVec[0], primeFacs.get());
        } else {
            Rcpp::List res(vSize);
            bool isNamed = false;
            
            if (!Rf_isNull(RNamed)) {
                if (TYPEOF(RNamed) == LGLSXP) {
                    double dblInp = Rcpp::as<double>(RNamed);
                    
                    if (Rcpp::NumericVector::is_na(dblInp) || std::isnan(dblInp))
                        Rcpp::stop("namedList cannot be NA or NaN");
                    
                    isNamed = Rcpp::as<bool>(RNamed);
                } else {
                    Rcpp::stop("Only logical values are supported for namedList");
                }
            }
            
            mpz_t val;
            mpz_init(val);
            
            if (isNamed) {
                constexpr int base10 = 10;
                Rcpp::CharacterVector myNames(vSize);

                for (std::size_t i = 0; i < vSize; ++i) {
                    auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(myVec[i], base10) + 2);
                    mpz_get_str(buffer.get(), base10, myVec[i]);
                    myNames[i] = Rf_mkChar(buffer.get());
                }
                
                for (std::size_t i = 0; i < vSize; ++i)
                    res[i] = factorNum(myVec[i], primeFacs.get());
                
                res.attr("names") = myNames;
            } else {
                for (std::size_t i = 0; i < vSize; ++i)
                    res[i] = factorNum(myVec[i], primeFacs.get());
            }

            return res;
        }
    }
    
    Rcpp::IntegerVector resTrivial(1);
    return resTrivial;
}
