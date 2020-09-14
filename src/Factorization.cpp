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
    
    if (vSize > 0) {
        if (vSize == 1) {
            mpz_class myNum;
            convertMpzClass(Rv, myNum);
            return FactorNum(myNum);
        } else {
            std::vector<mpz_class> myVec(vSize);
            CreateMPZVector(Rv, myVec, vSize);
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
                Rcpp::CharacterVector myNames(vSize);

                for (std::size_t i = 0; i < vSize; ++i)
                    myNames[i] = myVec[i].get_str();
                
                for (std::size_t i = 0; i < vSize; ++i)
                    res[i] = FactorNum(myVec[i]);
                
                res.attr("names") = myNames;
            } else {
                for (std::size_t i = 0; i < vSize; ++i)
                    res[i] = FactorNum(myVec[i]);
            }
            
            return res;
        }
    }
    
    Rcpp::IntegerVector resTrivial(1);
    return resTrivial;
}
