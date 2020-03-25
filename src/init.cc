#include "Factorization.h"
#include "RSAFactorize.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
// Factorization.h :
    CALLDEF(getDivisorsC, 2),
// RSAFactorize.h :
    CALLDEF(QuadraticSieveContainer, 1),
    {NULL, NULL, 0}
};

extern "C"
void R_init_gmp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

