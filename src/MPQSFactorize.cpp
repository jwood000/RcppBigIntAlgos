#include "PrimeFactorUtils.h"
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
    
    if (vSize > 1)
        Rcpp::stop("Can only factor one number at a time");
    
    mpz_class nMpz;
    convertMpzClass(Rn, nMpz);
    
    int nThreads = 1;
    const bool bShowStats = convertLogical(RShowStats, "showStats");
    
    if (!Rf_isNull(RNumThreads))
        convertInt(RNumThreads, nThreads, "nThreads");
    
    return PrimeFactorizeHuge(nMpz, nThreads, bShowStats, true, true);
}
