#include "FactorUtils.h"

// [[Rcpp::export]]
SEXP GetDivisorsC(SEXP Rv, SEXP RNamed, SEXP RNumThreads, int maxThreads) {
    
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
    
    CreateMPZArray(Rv, myVec.get(), vSize);
    
    for (std::size_t i = 0; i < mpzChunkBig; ++i)
        mpz_init(primeFacs[i]);
    
    if (vSize > 0) {
        if (vSize == 1) {
            return FactorNum(myVec[0], primeFacs);
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
            
            if (isNamed) {
                constexpr int base10 = 10;
                Rcpp::CharacterVector myNames(vSize);

                for (std::size_t i = 0; i < vSize; ++i) {
                    auto buffer = FromCpp14::make_unique<char[]>(mpz_sizeinbase(myVec[i], base10) + 2);
                    mpz_get_str(buffer.get(), base10, myVec[i]);
                    myNames[i] = Rf_mkChar(buffer.get());
                }
                
                for (std::size_t i = 0; i < vSize; ++i)
                    res[i] = FactorNum(myVec[i], primeFacs);
                
                res.attr("names") = myNames;
            } else {
                for (std::size_t i = 0; i < vSize; ++i)
                    res[i] = FactorNum(myVec[i], primeFacs);
            }

            return res;
        }
    }
    
    for (std::size_t i = 0; i < mpzChunkBig; ++i)
        mpz_clear(primeFacs[i]);
    
    for (std::size_t i = 0; i < vSize; ++i)
        mpz_clear(myVec[i]);
    
    myVec.reset();
    primeFacs.reset();
    Rcpp::IntegerVector resTrivial(1);
    return resTrivial;
}
