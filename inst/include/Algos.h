#ifndef ALGOS_H
#define ALGOS_H

#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <gmp.h>
#include <vector>
#include <cstdint>

#define USE_RINTERNALS
#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

#undef PROTECT
#undef UNPROTECT
#define PROTECT(x) Rf_protect(x)
#define UNPROTECT(x) Rf_unprotect(x)
#undef Length
#define Length(x) Rf_length(x)

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("main", String)
#else
#define _(String) (String)
#endif

#endif
