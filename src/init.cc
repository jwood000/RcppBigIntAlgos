// This ensures registration -- see also useDynLib(...) in  ../NAMESPACE

#include <R.h>
#include <Rinternals.h>

// include those that have an  extern "C" { .... } :
#include "factorization.h"
#include "rsafactorize.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
// factorization.h :
    CALLDEF(getDivisorsC, 2),
// rsafactorize.h :
    CALLDEF(QuadraticSieveContainer, 1),
    {NULL, NULL, 0}
};

extern "C"
void R_init_gmp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

