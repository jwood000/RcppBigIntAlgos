#include "cpp11/strings.hpp"
#include "cpp11/list.hpp"

#include "PrimeFactorUtils.h"
#include "CppConvert.h"

[[cpp11::register]]
SEXP PrimeFactorization(SEXP Rv, SEXP RNamed, SEXP RShowStats,
                        SEXP RSkipExtPR, SEXP RSkipECM,
                        SEXP RNumThreads, int maxThreads) {

    std::size_t vSize = 0;

    switch (TYPEOF(Rv)) {
        case RAWSXP: {
            const char* raw = (char*)RAW(Rv);
            vSize = ((int*)raw)[0];
            break;
        } default: {
            vSize = LENGTH(Rv);
        }
    }

    int nThreads = 1;
    const bool bShowStats = CppConvert::convertFlag(RShowStats, "showStats");
    const bool bSkipPR = CppConvert::convertFlag(RSkipExtPR, "skipPolRho");
    const bool bSkipECM = CppConvert::convertFlag(RSkipECM, "skipECM");

    if (!Rf_isNull(RNumThreads)) {
        CppConvert::convertPrimitive(
            RNumThreads, nThreads, VecType::Integer, "nThreads"
        );
    }

    if (vSize > 0) {
        if (vSize == 1) {
            mpz_class myNum;
            CppConvert::convertMpzClass(Rv, myNum, "n", true);

            cpp11::sexp res = PrimeFactorizeHuge(
                myNum, nThreads, bShowStats, bSkipPR, bSkipECM
            );

            return res;
        } else {
            std::vector<mpz_class> myVec(vSize);
            CppConvert::convertMPZVector(Rv, myVec, vSize, "v", true);
            cpp11::writable::list res(vSize);
            const bool isNamed = CppConvert::convertFlag(RNamed, "namedList");

            if (isNamed) {
                cpp11::writable::strings myNames(vSize);

                for (std::size_t i = 0; i < vSize; ++i)
                    myNames[i] = myVec[i].get_str();

                for (std::size_t i = 0; i < vSize; ++i) {
                    res[i] = PrimeFactorizeHuge(
                        myVec[i], nThreads, bShowStats, bSkipPR, bSkipECM
                    );
                }

                res.attr("names") = myNames;
            } else {
                for (std::size_t i = 0; i < vSize; ++i) {
                    res[i] = PrimeFactorizeHuge(
                        myVec[i], nThreads, bShowStats, bSkipPR, bSkipECM
                    );
                }
            }

            return res;
        }
    }

    cpp11::writable::integers resTrivial(1);
    return resTrivial;
}
