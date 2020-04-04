#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#include <Rcpp.h>

bool convertLogical(SEXP boolInput, const std::string &nameOfBool);
void convertInt(SEXP input, int &result, const std::string &nameOfObject);

#endif
