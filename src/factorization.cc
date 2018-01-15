/*! 
 *  \file factorization.cc
 *  \brief C functions used for integer factorization
 *
 *  \version 1
 *
 *  \date Created: 10/06/17
 *  \date Last modified: Time-stamp: <2017-10-06 12:33:33 EDT jwood000>
 *
 *  \author Joseph Wood. The first 20 lines (although slightly modified)
 *       are attributed to Antoine Lucas and help from Immanuel Scholz.
 *       See factor.cc from the R gmp package for more details.
 *          Original C code from libgmp.
 *
 *  \note Licence: GPL (>=) 2
 */

#include "Rgmp.h"
#include "pollardrho.h"
#include "factorization.h"

unsigned long int intSize = sizeof(int); // starting with vector-size-header
unsigned long int numb = 8*intSize;

int myRaw(char* raw, mpz_t value, unsigned long int totals) {
    memset(raw, 0, totals);
    
    int* r = (int*)raw;
    r[0] = totals/intSize - 2;
    
    r[1] = (int) mpz_sgn(value);
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);
    
    return totals;
}

std::vector<unsigned int> myMergeSort(mpz_t arr[], std::vector<unsigned int> indPass,
                 unsigned int numSecs, unsigned int secSize) {
    
    unsigned int i, j, k, x, y, left, count, lim, tempSize, totalSize;
    tempSize = totalSize = numSecs * secSize;
    std::vector<unsigned int> leftOver, myInd(totalSize);
    for (i = 0; i < secSize; i++) {myInd[i] = indPass[i];}
    for (i = secSize; i < totalSize; i++) {myInd[i] = i;}
    
    while (numSecs > 1) {

        if (numSecs % 2 != 0) {
            leftOver.push_back(tempSize);
            tempSize -= secSize;
        }
        
        std::vector<unsigned int> endPoints(numSecs, secSize);
        
        for (i = 1; i < numSecs; i++) {
            x = endPoints[i-1] + secSize;
            endPoints[i] = x;
        }
        
        k = std::floor(numSecs/2);
        if (k < 2) {lim = k;} else {lim = 2;}
        left = x = 0;
        std::vector<unsigned int> tempInd, defaultInd(secSize*2, 0);
        
        for (i = 0; i < lim; i++) {
            count = 0;
            y = endPoints[2*i];
            tempInd = defaultInd;
                     
            while (x < endPoints[2*i] && y < endPoints[2*i + 1]) {
                if (mpz_cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
                    tempInd[count] = myInd[x];
                    x++;
                } else {
                    tempInd[count] = myInd[y];
                    y++;
                }
                count++;
            }
            
            for (j = 0; j < count; j++) {myInd[left + j] = tempInd[j];}
            x = left = endPoints[2*i + 1];
        }
        
        if (k > 2) {
            for (i = 2; i < k; i++) {
                x = endPoints[2*i - 1];
                y = x - secSize*2;
                for (j = 0; j < count; j++) {
                    myInd[x + j] = tempInd[j] + y;
                }
            }
        }

        secSize *= 2;
        numSecs = std::floor(numSecs/2);
    }

    if (leftOver.size() > 0) {
        std::reverse(leftOver.begin(), leftOver.end());

        for (j = 0; j < leftOver.size(); j++) {
            x = count = 0;
            y = tempSize;
            std::vector<unsigned int> tempInd(leftOver[j]);

            while (x < tempSize && y < leftOver[j]) {
                if (mpz_cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
                    tempInd[count] = myInd[x];
                    x++;
                } else {
                    tempInd[count] = myInd[y];
                    y++;
                }
                count++;
            }
            
            for (k = 0; k < count; k++) {myInd[k] = tempInd[k];}
            tempSize = leftOver[j];
        }

    }

    return myInd;
}

SEXP getDivisorsC (SEXP n) {
    bigvec v = bigintegerR::create_bignum(n);
    bigvec primeFacs;
    
    if(v.size() > 0) {
        if (mpz_cmp_ui(v[0].value.getValueTemp(), 1) == 0) {
            bigvec myFacs;
            myFacs.push_back(1);
            return bigintegerR::create_SEXP(myFacs);
        } else {
            mpz_t val;
            mpz_init(val);
            mpz_t_sentry val_s(val);
            mpz_set(val,v[0].value.getValueTemp());
            
            int sgn = mpz_sgn(val);
            if(sgn == 0)
                error(_("Cannot factorize 0"));
            if (sgn<0) {
                mpz_abs(val,val);
                primeFacs.value.push_back(biginteger(-1));
            }
            
            getPrimeFactors (val,primeFacs);
            
            std::vector<int> lengths;
            std::vector<biginteger>::iterator it, primeEnd;
            primeEnd = primeFacs.value.end();
            biginteger prev = primeFacs.value[0];
            
            unsigned long int i, j, k, n = primeFacs.size(), numUni = 0;
            mpz_t bigFacs[n];
            lengths.reserve(n);
            for (i = 0; i < n; i++) {mpz_init(bigFacs[i]);}
            
            mpz_set(bigFacs[0], primeFacs.value[0].getValue());
            lengths.push_back(1);
            k = 1;
            
            for(it = primeFacs.value.begin() + 1; it < primeEnd; it++) {
                if (prev == *it) {
                    lengths[numUni]++;
                } else {
                    numUni++;
                    prev = *it;
                    lengths.push_back(1);
                    mpz_set(bigFacs[numUni], primeFacs[k].value.getValue());
                }
                k++;
            }
            
            unsigned long int ind, facSize = 1, numFacs = 1;
            for (i = 0; i <= numUni; i++) {numFacs *= (lengths[i]+1);}
            
            mpz_t *myMPZ;
            
            myMPZ = (mpz_t *) malloc(numFacs * sizeof(mpz_t));
            for (i = 0; i < numFacs; i++) {mpz_init(myMPZ[i]);}
            mpz_t temp, myPow;
            mpz_init(temp);
            mpz_init(myPow);
            std::vector<unsigned int> myIndex;
            myIndex.reserve(lengths[0] + 1);
            
            for (i = 0; i <= lengths[0]; ++i) {
                mpz_pow_ui(temp, bigFacs[0], i);
                mpz_set(myMPZ[i], temp);
                myIndex.push_back(i);
            }
            
            if (numUni > 0) {
                for (j = 1; j <= numUni; j++) {
                    facSize *= (lengths[j-1] + 1);
                    for (i = 1; i <= lengths[j]; i++) {
                        ind = i*facSize;
                        mpz_pow_ui(myPow, bigFacs[j], i);
                        for (k = 0; k < facSize; k++) {
                            mpz_mul(temp, myPow, myMPZ[myIndex[k]]);
                            mpz_set(myMPZ[ind + k], temp);
                        }
                    }
                    myIndex = myMergeSort(myMPZ, myIndex, lengths[j] + 1, facSize);
                }
            }
            
            unsigned long int tempSize, size = intSize;
            std::vector<unsigned long int> mySizes(numFacs);
                
            for (i = 0; i < numFacs; ++i) { // adding each bigint's needed size
                tempSize = intSize * (2 + (mpz_sizeinbase(myMPZ[i],2)+numb-1) / numb);
                size += tempSize;
                mySizes[i] = tempSize;
            }
                
            SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size));
            char* r = (char*)(RAW(ans));
            ((int*)(r))[0] = numFacs; // first int is vector-size-header
            
            unsigned long int pos = intSize; // current position in r[] (starting after vector-size-header)
            for (i = 0; i < numFacs; ++i)
                pos += myRaw(&r[pos], myMPZ[myIndex[i]], mySizes[myIndex[i]]);
            
            Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));
            UNPROTECT(1);
            return(ans);
        }
    }
}