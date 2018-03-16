// This ensures registration -- see also useDynLib(...) in  ../NAMESPACE

#include <R.h>
#include <Rinternals.h>

// include those that have an  extern "C" { .... } :
#include "bigintegerR.h"
#include "factorization.h"
#include "quadres.h"
#include "rsafactorize.h"
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
// bigintegerR.h :
  CALLDEF(R_gmp_get_version, 0),
  CALLDEF(biginteger_add, 2),
  CALLDEF(biginteger_sub, 2),
  CALLDEF(biginteger_mul, 2),
  CALLDEF(biginteger_div, 2),
  CALLDEF(biginteger_divq, 2),
  CALLDEF(biginteger_mod, 2),
  CALLDEF(biginteger_pow, 2),
  CALLDEF(biginteger_inv, 2),
  CALLDEF(biginteger_gcd, 2),
  CALLDEF(biginteger_lcm, 2),
  CALLDEF(biginteger_get_at, 2),
  CALLDEF(biginteger_set_at, 3),
  CALLDEF(biginteger_as, 2),
  CALLDEF(biginteger_as_character, 2),
  CALLDEF(biginteger_as_numeric, 1),
  CALLDEF(biginteger_as_integer, 1),
  CALLDEF(biginteger_length, 1),
  CALLDEF(biginteger_setlength, 2),
  CALLDEF(biginteger_is_na, 1),
  CALLDEF(biginteger_sgn, 1),
  CALLDEF(biginteger_lt, 2),
  CALLDEF(biginteger_gt, 2),
  CALLDEF(biginteger_lte, 2),
  CALLDEF(biginteger_gte, 2),
  CALLDEF(biginteger_eq, 2),
  CALLDEF(biginteger_neq, 2),
  CALLDEF(biginteger_c, 1),
  CALLDEF(biginteger_cbind, 1),
  CALLDEF(biginteger_rep, 2),
  CALLDEF(biginteger_is_prime, 2),
  CALLDEF(biginteger_nextprime, 1),
  CALLDEF(biginteger_abs, 1),
  CALLDEF(biginteger_gcdex, 2),
  CALLDEF(biginteger_rand_u, 4),
  CALLDEF(biginteger_sizeinbase, 2),
  CALLDEF(bigI_frexp, 1),
  CALLDEF(bigI_choose, 2),
  CALLDEF(bigI_factorial, 1),
  CALLDEF(bigI_fibnum, 1),
  CALLDEF(bigI_fibnum2, 1),
  CALLDEF(bigI_lucnum, 1),
  CALLDEF(bigI_lucnum2, 1),
  CALLDEF(biginteger_max, 2),
  CALLDEF(biginteger_min, 2),
  CALLDEF(biginteger_cumsum, 1),
  CALLDEF(biginteger_sum, 1),
  CALLDEF(biginteger_prod, 1),
  CALLDEF(biginteger_powm, 3),
  CALLDEF(biginteger_log2, 1),
  CALLDEF(biginteger_log, 1),

// factorization.h :
  CALLDEF(getDivisorsC, 2),

// quadres.h :
  CALLDEF(QuadraticResidueContainer, 2),
  
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

