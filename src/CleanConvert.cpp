#include <Rcpp.h>

bool convertLogical(SEXP boolInput, const std::string &nameOfBool) {
    bool result = false;
    
    if (!Rf_isNull(boolInput)) {
        if (TYPEOF(boolInput) == LGLSXP) {
            double dblInp = Rcpp::as<double>(boolInput);
            
            if (Rcpp::NumericVector::is_na(dblInp) || std::isnan(dblInp))
                Rcpp::stop(nameOfBool + " cannot be NA or NaN");
            
            result = Rcpp::as<bool>(boolInput);
        } else {
            Rcpp::stop("Only logical values are supported for " + nameOfBool);
        }
    }
    
    return result;
}

void convertInt(SEXP input, int &result, const std::string &nameOfObject) {
    
    const int maxType = std::numeric_limits<int>::max();
    
    switch(TYPEOF(input)) {
        case REALSXP:
        case INTSXP: {
            const double dblInp = Rcpp::as<double>(input);
            
            if (Rcpp::NumericVector::is_na(dblInp) || std::isnan(dblInp))
                Rcpp::stop(nameOfObject + " cannot be NA or NaN");
            
            if (dblInp < 1)
                Rcpp::stop(nameOfObject + " must be a positive whole number");
            
            if (dblInp > maxType)
                Rcpp::stop("The abs value of " + nameOfObject + " must be less than"
                               " or equal to " + std::to_string(maxType));
            
            if (static_cast<int64_t>(dblInp) != dblInp)
                Rcpp::stop(nameOfObject + " must be a whole number");
                
            result = Rcpp::as<int>(input);
            break;
        }
        default:
            Rcpp::stop("This type is not supported! No conversion possible for " + nameOfObject);
    }
}