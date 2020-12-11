// This file contains utility functions that
// are used for converting to and from type mpz_t,
// as well as sorting an array of type mpz_t.
// 
// createMPZArray and myRaw are slightly modified versions
// of "bigvec create_vector(const SEXP & param)" and 
// "int biginteger::as_raw(char* raw) const", respectively,
// from the source files bigintegerR.cc/ biginteger.cc from
// the R gmp package.
// 
// The quickSort function is based off of the quicksort
// algorithm in C++ found here:
//      http://www.algolist.net/Algorithms/Sorting/Quicksort

#include "ImportExportMPZ.h"

void CreateMPZVector(SEXP v, std::vector<mpz_class> &myVec, std::size_t sizevec) {
    
    switch (TYPEOF(v)) {
        case RAWSXP: {
            // deserialise the vector. first int is the size.
            const char* raw = (char*) RAW(v);
            int pos = intSize; // position in raw[]. Starting after header.
            
            for (std::size_t i = 0; i < sizevec; ++i) {
                const int* r = (int*) (&raw[pos]);
                
                if (r[0] > 0) {
                    mpz_import(myVec[i].get_mpz_t(), r[0], 1, intSize, 0, 0, (void*) &(r[2]));
                    
                    if(r[1] == -1)
                        mpz_neg(myVec[i].get_mpz_t(), myVec[i].get_mpz_t());
                } else {
                    myVec[i] = 0;
                }
                
                pos += intSize * (2 + (mpz_sizeinbase(myVec[i].get_mpz_t(), 2) + numb - 1) / numb);
            }
            
            break;
        }
        case REALSXP: {
            const std::vector<double> myDbl = Rcpp::as<std::vector<double>>(v);
            constexpr double Sig53 = 9007199254740991.0;
            
            for (std::size_t j = 0; j < sizevec; ++j) {
                if (Rcpp::NumericVector::is_na(myDbl[j]) || std::isnan(myDbl[j]))
                    Rcpp::stop("Elements in v cannot be NA or NaN");
                
                if (myDbl[j] > Sig53) {
                    Rcpp::stop("Number is too large for double precision. "
                                   "Consider wrapping v with gmp::as.bigz or "
                                   "as.character (e.g. gmp::as.bigz(v))");
                }
                
                if (static_cast<int64_t>(myDbl[j]) != myDbl[j])
                    Rcpp::stop("Elements in v must be whole numbers");
                
                myVec[j] = myDbl[j];
            }
            
            break;
        }
        case INTSXP:
        case LGLSXP: {
            const std::vector<int> myInt = Rcpp::as<std::vector<int>>(v);
            const std::vector<double> dblVec = Rcpp::as<std::vector<double>>(v);
            
            for (std::size_t j = 0; j < sizevec; ++j) {
                if (Rcpp::NumericVector::is_na(dblVec[j]) || std::isnan(dblVec[j]))
                    Rcpp::stop("Elements in v cannot be NA or NaN");
                
                myVec[j] = myInt[j];
            }
            
            break;
        }
        case STRSXP: {
            for (std::size_t i = 0; i < sizevec; ++i) {
                if (STRING_ELT(v, i) == NA_STRING) {
                    Rcpp::stop("Elements in v  cannot be NA or NaN");
                } else {
                    myVec[i].set_str(CHAR(STRING_ELT(v, i)), 10);
                }
            }
            
            break;
        }
        default:
            // no longer: can be fatal later! return bigvec();
            Rcpp::stop("only logical, numeric or character (atomic) vectors can be coerced to 'bigz'");
    }
}

void convertMpzClass(SEXP v, mpz_class &result) {
    
    switch (TYPEOF(v)) {
        case RAWSXP: {
            // deserialise the vector. first int is the size.
            const char* raw = (char*) RAW(v);
            const int* r = (int*) (&raw[intSize]);
            
            if (r[0] > 0) {
                mpz_import(result.get_mpz_t(), r[0], 1, intSize, 0, 0, (void*) &(r[2]));
                
                if(r[1] == -1)
                    mpz_neg(result.get_mpz_t(), result.get_mpz_t());
            } else {
                result = 0;
            }
            
            break;
        }
        case REALSXP: {
            const double myDbl = Rcpp::as<double>(v);
            constexpr double Sig53 = 9007199254740991.0;
            
            if (Rcpp::NumericVector::is_na(myDbl) || std::isnan(myDbl))
                Rcpp::stop("Elements in v cannot be NA or NaN");
            
            if (myDbl > Sig53) {
                Rcpp::stop("Number is too large for double precision. "
                               "Consider wrapping v with gmp::as.bigz or "
                               "as.character (e.g. gmp::as.bigz(v))");
            }
            
            if (static_cast<int64_t>(myDbl) != myDbl)
                Rcpp::stop("Elements in v must be whole numbers");
            
            result = myDbl;
            break;
        }
        case INTSXP:
        case LGLSXP: {
            const int myInt = Rcpp::as<int>(v);
            const double dblVec = Rcpp::as<double>(v);
            
            if (Rcpp::NumericVector::is_na(dblVec) || std::isnan(dblVec))
                Rcpp::stop("Elements in v cannot be NA or NaN");
            
            result = myInt;
            break;
        }
        case STRSXP: {
            if (v == NA_STRING) {
                Rcpp::stop("Elements in v cannot be NA or NaN");
            } else {
                result = Rcpp::as<std::string>(v);
            }
        
            break;
        }
        default:
            // no longer: can be fatal later! return bigvec();
            Rcpp::stop("only logical, numeric or character (atomic) vectors can be coerced to 'bigz'");
    }
}

int myRaw(char* raw, mpz_t value, std::size_t totals) {
    memset(raw, 0, totals);
    
    int* r = (int*) raw;
    r[0] = totals / intSize - 2;
    
    r[1] = static_cast<int>(mpz_sgn(value));
    mpz_export(&r[2], 0, 1, intSize, 0, 0, value);
    
    return totals;
}

void QuickSort(std::vector<mpz_class> &arr, int left, int right,
               std::vector<std::size_t> &lens) {
    
    int i = left;
    int j = right;
    mpz_class pivot;
    
    int mid = (left + right) / 2;
    pivot = arr[mid];
    
    /* partition */
    while (i <= j) {
        while (cmp(arr[i], pivot) < 0)
            ++i;
        
        while (j >= 0 && cmp(arr[j], pivot) > 0)
            --j;
        
        if (i <= j && j >= 0) {
            swap(arr[i], arr[j]);
            std::swap(lens[i], lens[j]);
            ++i;
            --j;
        }
    }
    
    /* recursion */
    if (left < j)
        QuickSort(arr, left, j, lens);
    
    if (i < right)
        QuickSort(arr, i, right, lens);
}
