// Description:
//     myMergeSort is a modified merge sort algorithm. The classic merge
//     sort creates an empty vector of the type you are sorting and fills
//     in the ordered elements from two smaller sorted vectors. This
//     would be expensive with type mpz_t. For this reason, we only keep
//     track of the indices and never actually swap any of the indices in
//     the mpz_t array. An integer vector of indices is returned and is
//     used in writing out the output at the very bottom.  This method
//     shows great efficiency gains over the naive method of using the
//     class bigvec (from the R gmp package)/std::sort combination.
//
// Date Created: 10/06/17
// Date Last modified: Time-stamp: <2020-09-14 jwood000>
//
// Author Joseph Wood. Original C code from libgmp.
//    See factor.cc from the R gmp package for more details.
//
// Note Licence: GPL (>=) 2

#include "cpp11/raws.hpp"
#include "DivisorsUtils.h"
#include "PrimeFactorUtils.h"
#include <memory>

std::vector<int> myMergeSort(mpz_t *const arr, const std::vector<int> &indPass,
                             std::size_t numSecs, std::size_t secSize) {

    std::size_t count = 0;
    const std::size_t totalSize  = numSecs * secSize;

    int tempSize = totalSize;
    std::vector<int> leftOver;
    std::vector<int> myInd(totalSize);

    for (std::size_t i = 0; i < secSize; ++i) {
        myInd[i] = indPass[i];
    }

    for (std::size_t i = secSize; i < totalSize; ++i) {
        myInd[i] = i;
    }

    leftOver.reserve(numSecs / 2);

    while (numSecs > 1) {
        if (numSecs % 2 != 0) {
            leftOver.push_back(tempSize);
            tempSize -= static_cast<int>(secSize);
        }

        std::vector<int> endPoints(numSecs);

        for (std::size_t i = 0; i < numSecs; ++i) {
            endPoints[i] = secSize * (i + 1);
        }

        secSize *= 2;
        std::size_t k = numSecs / 2;

        int lim = (k < 2) ? k : 2;
        std::vector<int> tempInd(secSize);

        for (int x = 0, i = 0, left = 0; i < lim; ++i) {
            count = 0;
            int y = endPoints[2 * i];
            std::fill(tempInd.begin(), tempInd.end(), 0);

            while (x < endPoints[2 * i] && y < endPoints[2 * i + 1]) {
                if (mpz_cmp(arr[myInd[x]], arr[myInd[y]]) < 0) {
                    tempInd[count] = myInd[x];
                    ++x;
                } else {
                    tempInd[count] = myInd[y];
                    ++y;
                }

                ++count;
            }

            for (std::size_t j = 0; j < count; ++j) {
                myInd[left + j] = tempInd[j];
            }

            x = left = endPoints[2 * i + 1];
        }

        if (k > 2) {
            for (std::size_t i = 2; i < k; ++i) {
                const int x = endPoints[2 * i - 1];
                const int y = x - secSize;

                for (std::size_t j = 0; j < count; ++j) {
                    myInd[x + j] = tempInd[j] + y;
                }
            }
        }

        numSecs /= 2;
    }

    const int LOSize = leftOver.size();

    if (LOSize > 0) {
        for (int j = LOSize - 1; j >= 0; --j) {
            int x = count = 0;
            int y = tempSize;
            std::vector<int> tempInd(leftOver[j]);

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

            for (std::size_t r = 0; r < count; ++r) {
                myInd[r] = tempInd[r];
            }

            tempSize = leftOver[j];
        }
    }

    return myInd;
}

SEXP FactorNum(mpz_class &val, std::size_t nThreads,
               bool bShowStats, bool bSkipPR, bool bSkipECM) {

    if (cmp(val, 1) == 0) {
        mpz_class mpzOne = 1;
        cpp11::writable::raws myFacs(intSize * 4);

        char* r = (char*) (RAW(myFacs));
        ((int*) (r))[0] = 1;

        CppConvert::rawExport(&r[intSize], mpzOne.get_mpz_t(), intSize * 3);
        Rf_setAttrib(myFacs, R_ClassSymbol, Rf_mkString("bigz"));
        return myFacs;
    } else {
        const int mySgn = sgn(val);
        bool isNegative = false;

        std::vector<std::size_t> lengths;
        std::vector<mpz_class> primeFacs;

        if (mySgn == 0) {
            cpp11::stop("Cannot factorize 0");
        }

        if (mySgn < 0) {
            val = abs(val);
            isNegative = true;
        }

        if (mpz_sizeinbase(val.get_mpz_t(), 10) > 23) {
            QuadSieveHelper(val, primeFacs, lengths, nThreads,
                            bShowStats, bSkipPR, bSkipECM);
        } else {
            GetPrimeFactors(val, primeFacs, lengths);
        }

        CppConvert::QuickSort(primeFacs, 0, lengths.size() - 1, lengths);
        std::vector<int> myIndex(lengths[0] + 1);
        std::size_t numFacs = 1;

        for (std::size_t i = 0; i < lengths.size(); ++i) {
            numFacs *= (lengths[i] + 1);
        }

        auto myMPZ = std::make_unique<mpz_t[]>(numFacs);

        for (std::size_t i = 0; i < numFacs; ++i) {
            mpz_init(myMPZ[i]);
        }

        mpz_t myPow;
        mpz_init(myPow);

        for (std::size_t i = 0; i <= lengths[0]; ++i) {
            mpz_pow_ui(myMPZ[i], primeFacs[0].get_mpz_t(), i);
            myIndex[i] = i;
        }

        for (std::size_t j = 1, facSize = 1; j < lengths.size(); ++j) {
            facSize *= (lengths[j - 1] + 1);

            for (std::size_t i = 1; i <= lengths[j]; ++i) {
                mpz_pow_ui(myPow, primeFacs[j].get_mpz_t(), i);

                for (std::size_t k = 0, ind = i * facSize; k < facSize; ++k) {
                    mpz_mul(myMPZ[ind + k], myPow, myMPZ[myIndex[k]]);
                }
            }

            myIndex = myMergeSort(
                myMPZ.get(), myIndex, lengths[j] + 1, facSize
            );
        }

        std::size_t size = intSize;
        std::vector<std::size_t> mySizes(numFacs);

        // adding each bigint's needed size
        for (std::size_t i = 0; i < numFacs; ++i) {
            const std::size_t tempSize = intSize *
                (2 + (mpz_sizeinbase(myMPZ[i], 2) + numb - 1) / numb);
            size += tempSize;
            mySizes[i] = tempSize;
        }

        if (!isNegative) {
            cpp11::writable::raws ansPos(size);

            char* rPos = (char*) (RAW(ansPos));
            ((int*) (rPos))[0] = numFacs; // first int is vector-size-header

            // current position in rPos[] (starting after vector-size-header)
            std::size_t posPos = intSize;

            for (std::size_t i = 0; i < numFacs; ++i) {
                posPos += CppConvert::rawExport(
                    &rPos[posPos], myMPZ[myIndex[i]], mySizes[myIndex[i]]
                );
            }

            Rf_setAttrib(ansPos, R_ClassSymbol, Rf_mkString("bigz"));
            return(ansPos);
        } else {
            // double size as every element will have a negative counterpart
            size *= 2;
            size -= intSize; // Remove superfluous initializing size
            cpp11::writable::raws ansNeg(size);
            char* rNeg = (char*) (RAW(ansNeg));

            // first int is vector-size-header
            ((int*) (rNeg))[0] = 2 * numFacs;

            // current position in rNeg[] (starting after vector-size-header)
            std::size_t posNeg = intSize;
            mpz_t myNegative;
            mpz_init(myNegative);

            // First write out negative numbers in reverse "myIndex" order
            for (int i = numFacs - 1; i >= 0; --i) {
                mpz_neg(myNegative, myMPZ[myIndex[i]]);
                posNeg += CppConvert::rawExport(
                    &rNeg[posNeg], myNegative, mySizes[myIndex[i]]
                );
            }

            for (std::size_t i = 0; i < numFacs; ++i) {
                posNeg += CppConvert::rawExport(
                    &rNeg[posNeg], myMPZ[myIndex[i]], mySizes[myIndex[i]]
                );
            }

            Rf_setAttrib(ansNeg, R_ClassSymbol, Rf_mkString("bigz"));
            return(ansNeg);
        }
    }
}
