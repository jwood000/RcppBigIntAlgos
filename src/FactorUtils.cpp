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

#include "FactorUtils.h"
#include "RSAFactorUtils.h"

std::vector<std::size_t> myMergeSort(const std::vector<mpz_class> &arr, 
                                     const std::vector<std::size_t> &indPass,
                                     std::size_t numSecs, 
                                     std::size_t secSize) {
    
    std::size_t x, y;
    std::size_t count = 0;
    const std::size_t totalSize  = numSecs * secSize;
    
    std::size_t tempSize = totalSize;
    std::vector<std::size_t> leftOver;
    std::vector<std::size_t> myInd(totalSize);
    
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
        
        std::vector<std::size_t> tempInd(secSize);
        
        for (std::size_t i = 0; i < lim; ++i) {
            count = 0;
            y = endPoints[2 * i];
            std::fill(tempInd.begin(), tempInd.end(), 0);
                     
            while (x < endPoints[2 * i] && y < endPoints[2 * i + 1]) {
                if (cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
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
                x = endPoints[2 * i - 1];
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
                if (cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
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

SEXP FactorNum(mpz_class &val) {
    
    if (cmp(val, 1) == 0) {
        mpz_class mpzOne = 1;
        Rcpp::RawVector myFacs(intSize * 4);
        
        char* r = (char*) (RAW(myFacs));
        ((int*) (r))[0] = 1;
        
        myRaw(&r[intSize], mpzOne.get_mpz_t(), intSize * 3);
        myFacs.attr("class") = Rcpp::CharacterVector::create("bigz");
        return myFacs;
    } else {
        const int mySgn = sgn(val);
        bool isNegative = false;
        
        std::vector<std::size_t> lengths;
        std::vector<mpz_class> primeFacs;
        
        if (mySgn == 0)
            Rcpp::stop("Cannot factorize 0");
        
        if (mySgn < 0) {
            val = abs(val);
            isNegative = true;
        }
        
        if (mpz_sizeinbase(val.get_mpz_t(), 10) > 23) {
            std::size_t nThreads = 1;
            QuadSieveHelper(val, primeFacs, lengths, nThreads, false);
        } else {
            GetPrimeFactors(val, primeFacs, lengths);
        }
        
        QuickSort(primeFacs, 0, lengths.size() - 1, lengths);
        
        std::vector<std::size_t> myIndex(lengths[0] + 1);
        std::size_t facSize = 1;
        std::size_t numFacs = 1;
        
        for (std::size_t i = 0; i < lengths.size(); ++i)
            numFacs *= (lengths[i] + 1);
        
        std::vector<mpz_class> myMPZ(numFacs);
        mpz_class temp, myPow;
        
        for (std::size_t i = 0; i <= lengths[0]; ++i) {
            mpz_pow_ui(temp.get_mpz_t(), primeFacs[0].get_mpz_t(), i);
            myMPZ[i] = temp;
            myIndex[i] = i;
        }
        
        for (std::size_t j = 1; j < lengths.size(); ++j) {
            facSize *= (lengths[j - 1] + 1);
            
            for (std::size_t i = 1; i <= lengths[j]; ++i) {
                const std::size_t ind = i*facSize;
                mpz_pow_ui(myPow.get_mpz_t(), primeFacs[j].get_mpz_t(), i);
                
                for (std::size_t k = 0; k < facSize; ++k)
                    myMPZ[ind + k] = myPow * myMPZ[myIndex[k]];
            }
            
            myIndex = myMergeSort(myMPZ, myIndex, lengths[j] + 1, facSize);
        }
        
        std::size_t size = intSize;
        std::vector<std::size_t> mySizes(numFacs);
        
        for (std::size_t i = 0; i < numFacs; ++i) { // adding each bigint's needed size
            const std::size_t tempSize = intSize * (2 + (mpz_sizeinbase(myMPZ[i].get_mpz_t(), 2) + numb - 1) / numb);
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
                posPos += myRaw(&rPos[posPos], myMPZ[myIndex[i]].get_mpz_t(), mySizes[myIndex[i]]);
            
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
                mpz_neg(temp.get_mpz_t(), myMPZ[myIndex[i]].get_mpz_t());
                posNeg += myRaw(&rNeg[posNeg], temp.get_mpz_t(), mySizes[myIndex[i]]);
            }
            
            for (std::size_t i = 0; i < numFacs; ++i)
                posNeg += myRaw(&rNeg[posNeg], myMPZ[myIndex[i]].get_mpz_t(), mySizes[myIndex[i]]);
            
            ansNeg.attr("class") = Rcpp::CharacterVector::create("bigz");
            return(ansNeg);
        }
    }
}
