/************************************************************/
/*! \file bigintegerR.cc
 *  \brief C function to interface R and libgmp with big integer values
 *
 *  \version 1
 *
 *  \date Created: 27/10/04
 *  \date Last modified: $Id: bigintegerR.cc,v 1.37 2014-09-16 07:33:43 mmaechler Exp $
 *
 *  \author Immanuel Scholz (help from A. Lucas)
 *
 *  \note Licence: GPL
 */

#include "Rgmp.h"

#include "bigintegerR.h"

#include <vector>
#include <algorithm>

using namespace std;

#include <Rmath.h>

/* Globals variables */

static gmp_randstate_t seed_state;
static int seed_init=0;

namespace bigintegerR
{
  // \brief create a vector of bigvecs, all without a modulus.
  bigvec create_vector(const SEXP & param) {
    switch (TYPEOF(param)) {
    case NILSXP:
	return bigvec(); // = bigz(0)
    case RAWSXP:
      {
	// deserialise the vector. first int is the size.
	bigvec v;
	const char* raw = (char*)RAW(param);
	int pos = sizeof(int); // position in raw[]. Starting after header.
	int sizevec = ((int*)raw)[0];
	//std::cout << "nb element a lire " << sizevec << std::endl;
	v.value.resize(sizevec);
	for (int i = 0; i < sizevec; ++i) {
	  v.value[i] = biginteger(&raw[pos]);
	  pos += v.value[i].raw_size(); // increment number of bytes read.
	}
	return v;
      }
    case REALSXP:
      {
	double* d = REAL(param);
	//bigvec v(d,d+LENGTH(param));
	bigvec v;
	v.value.resize(LENGTH(param));
	for (int j = 0; j < LENGTH(param); ++j) {
	    /// New:   numeric '+- Inf'  give  +- "Large" instead of NA
	    double dj = d[j];
	    if(R_FINITE(dj) || ISNAN(dj))
		v.value[j] = dj;
	    else { // dj is +- Inf : use LARGE ( =   +- 2 ^ 80000 -- arbitrarily )
		mpz_t LARGE;
		mpz_init(LARGE);
		// FIXME: Keep 'LARGE' a static const; initialized only once
		mpz_ui_pow_ui (LARGE, (unsigned long int) 2, (unsigned long int) 8000);
		if(dj == R_PosInf)
		    v.value[j] = LARGE;
		else if(dj == R_NegInf) {
		    mpz_t neg_L;
		    mpz_init(neg_L);
		    mpz_neg(neg_L, LARGE);
		    v.value[j] = neg_L;
		    mpz_clear(neg_L);
		}
		else// should never happen
		    v.value[j] = dj;

		mpz_clear(LARGE);
	    }
	}
	return v;
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);
	//bigvec v(i,i+LENGTH(param));
	bigvec v;
	v.value.resize(LENGTH(param));

	for (int j = 0; j < LENGTH(param); ++j)
	    v.value[j] = i[j];

	return v;
      }
    case STRSXP:
      {
	bigvec v;
	v.value.resize(LENGTH(param));
	for (int i = 0; i < LENGTH(param); ++i) {
	  if (STRING_ELT(param,i) == NA_STRING)
	    v.value[i]= biginteger();
	  else
	    v.value[i]=biginteger(std::string(CHAR(STRING_ELT(param,i))));
	}
	return v;
      }
    default:
	// no longer: can be fatal later! return bigvec();
	error(_("only logical, numeric or character (atomic) vectors can be coerced to 'bigz'"));
    }
  }

  bigvec create_bignum(const SEXP & param) {
      SEXP
	  modAttr = Rf_getAttrib(param, Rf_mkString("mod")),
	  dimAttr = Rf_getAttrib(param, Rf_mkString("nrow"));

    // try to catch biz-nrow dimension value
    //std::cout << "import value" << std::endl;
    bigvec v = bigintegerR::create_vector(param);

    if (TYPEOF(dimAttr) == INTSXP)
	v.nrow = INTEGER(dimAttr)[0];
    else {
	// catch to get std matrix dimensions value
	dimAttr = Rf_getAttrib(param, Rf_mkString("dim"));
	v.nrow = (TYPEOF(dimAttr) == INTSXP) ? INTEGER(dimAttr)[0] : -1;// -1: want support 0-row
    }

    if (TYPEOF(modAttr) != NILSXP)
      {
	//std::cout << "import value" << std::endl;
	v.modulus = bigintegerR::create_vector(modAttr).value;
      }
    return v;

  }

  std::vector<int> create_int(const SEXP & param) {
    switch (TYPEOF(param)) {
    case REALSXP:
      {
	double* d = REAL(param);
	// copy vector manually to avoid stupid conversion warning in STL-code :-/
	vector<int> v;
	v.reserve(LENGTH(param));
	for (int i = 0; i < LENGTH(param); ++i)
	  v.push_back(static_cast<int>(d[i]));
	return v;
	//return vector<int>(d, d+LENGTH(param));
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);

	return std::vector<int>(i, i+LENGTH(param));
      }
    default:
      return std::vector<int>();
    }
  }


  SEXP create_SEXP(const std::vector<biginteger>& v)
  {
    unsigned int i;
    int size = sizeof(int); // starting with vector-size-header
    for (i = 0; i < v.size(); ++i)
      size += v[i].raw_size(); // adding each bigint's needed size
    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, size));
    // Rprintf("    o create_SEXP(vect<biginteger>): size=%d, v.size()=%d\n", size, v.size());
    char* r = (char*)(RAW(ans));
    ((int*)(r))[0] = v.size(); // first int is vector-size-header
    int pos = sizeof(int); // current position in r[] (starting after vector-size-header)
    for (i = 0; i < v.size(); ++i)
      pos += v[i].as_raw(&r[pos]);
    UNPROTECT(1);
    return(ans);
  }

  SEXP create_SEXP(const bigvec& v) {

    SEXP ans = PROTECT(create_SEXP(v.value));
    // set the class attribute to "bigz"
    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));

    // Rprintf("   o create_SEXP(<bigvec>): v.nrow=%d", v.nrow);
    // set the dim attribute
    if(v.nrow >= 0) // {
      Rf_setAttrib(ans, Rf_mkString("nrow"), Rf_ScalarInteger((int) v.nrow));
      // Rprintf(" *SET*\n");
      // } else Rprintf(" no set\n");
    // set the mod attribute
    if(v.modulus.size() > 0) {
      SEXP mod = PROTECT(create_SEXP(v.modulus)); // and set *its* class
      Rf_setAttrib(mod, R_ClassSymbol, Rf_mkString("bigz"));
      Rf_setAttrib(ans, Rf_mkString("mod"), mod);
      UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
  }

  /**
   * \brief Main function of doing a binary operation on bigintegers.
   * It calls a function argument for doing the correct thing.
   * This could also be written as a class functor (template)
   * to save one function call, but then code bloat will happen.
   */
  SEXP biginteger_binary_operation(const SEXP& a,const SEXP& b, biginteger_binary_fn f)
  {
    bigvec va = bigintegerR::create_bignum(a);
    bigvec vb = bigintegerR::create_bignum(b), result;
    int size = (va.value.empty() || vb.value.empty()) ? 0 : max(va.value.size(), vb.value.size());
    result.value.reserve(size);
    for (int i = 0; i < size; ++i)
	result.push_back(f(va[i%va.size()], vb[i%vb.size()]));
    
    // Don't need matrix - Joseph Wood 10/2/17
    // result.nrow = matrixz::checkDims(va.nrow,vb.nrow);
    // Rprintf(" o bigI_b_op(.); size=%d -> nrow=%d\n", size, result.nrow);
    return bigintegerR::create_SEXP(result);
  }


  SEXP biginteger_logical_binary_operation(const SEXP & a,const SEXP & b, biginteger_logical_binary_fn f)
  {
    bigvec va = bigintegerR::create_bignum(a);
    bigvec vb = bigintegerR::create_bignum(b), result;
    int size = (va.value.empty() || vb.value.empty()) ? 0 : max(va.value.size(), vb.value.size());
    //	int sizemod = max(va.modulus.size(), vb.modulus.size());
    SEXP ans = PROTECT(Rf_allocVector(LGLSXP, size));
    int *r = LOGICAL(ans);
    /* TODO: this kind of situation 5 == (5 %% 17)*/
    for (int i = 0; i < size; ++i) {
      biginteger am = va.value[i % va.value.size()];
      biginteger bm = vb.value[i % vb.value.size()];
      if (am.isNA() || bm.isNA())
	r[i] = NA_LOGICAL;
      else
	r[i] = f(am, bm) ? 1 : 0;
    }
    
    // Don't need matrix - Joseph Wood 10/2/17
    // int nrow = matrixz::checkDims(va.nrow,vb.nrow) ;

    // Add dimension parameter when available
//     if(nrow >= 0)
//       {
// 	SEXP dimVal;
// 	PROTECT(dimVal = Rf_allocVector(INTSXP, 2));
// 	INTEGER(dimVal)[0] = (int) nrow;
// 	INTEGER(dimVal)[1] = (int) size / nrow;
// 	Rf_setAttrib(ans, Rf_mkString("dim"), dimVal);
// 	UNPROTECT(1);
//       }

    UNPROTECT(1);
    return ans;
  }

  bool lt(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) < 0;}
  bool gt(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) > 0;}
  bool lte(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) <= 0;}
  bool gte(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) >= 0;}
  bool eq(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) == 0;}
  bool neq(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) != 0;}

}

/* End of namespace bigintegerR*/


SEXP R_gmp_get_version() {
    return Rf_mkString(gmp_version);
}

SEXP biginteger_add (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator+);}
SEXP biginteger_sub (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator-);}
SEXP biginteger_mul (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator*);}
SEXP biginteger_divq(SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator/);}
SEXP biginteger_mod (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator%);}

SEXP biginteger_div (SEXP a, SEXP b) { // called from  "/.bigz" == div.bigz
    bigvec A = bigintegerR::create_bignum(a),
	B = bigintegerR::create_bignum(b);
    // Note: a or b may be simple numbers (e.g. '1') => work with (A,B)
    int len_m_a = A.modulus.size(),
	len_m_b = B.modulus.size();
    // Don't need bigrational - Joseph Wood 10/2/17
    // if(len_m_a == 0 && len_m_b == 0) // deal with important case quickly:
	// return bigrational_div(a, b);
    if(len_m_a == 0) { // and  len_m_b > 0:
	// should work directly using b's "mod" --> compute  a * b^(-1)
    }
    else if(len_m_b == 0) { // and  len_m_a > 0:
	// should use a's "mod" for b: div_via_inv() need's  b's modulus
	B.modulus = A.modulus;
	return bigintegerR::biginteger_binary_operation(a,
							bigintegerR::create_SEXP(B),
							div_via_inv);
    }
    else { // len_m_a > 0 and  len_m_b > 0:
	bool same_mod = true;// are the two mods the "same" (after recycling)?
	int m = (len_m_a < len_m_b) ? len_m_b : len_m_a; // = max(l..a, l..b)
	for(int i = 0; i < m; i++)
	    if(A.modulus[i % len_m_a] != B.modulus[i % len_m_b])  {
		same_mod = false; break;
	    }
	if(same_mod) {
	    // compute   a * b^(-1) ... should work w/o more
	} else {
	    // use *rational* a/b  (not considering 'mod' anymore):
	    // Don't need bigrational - Joseph Wood 10/2/17
	    // return bigrational_div(a, b);
	}
    }
    return bigintegerR::biginteger_binary_operation(a,b, div_via_inv);
}

SEXP biginteger_pow (SEXP a, SEXP b) {
  bigvec v = bigintegerR::create_bignum(a),
    exp = bigintegerR::create_bignum(b);
  if(v.modulus.size() == 0) { /* has no modulus: now, if any b < 0, the
				 result must be (non-integer) bigrational */
    bool use_rat = FALSE;
    for (unsigned int i = 0; i < exp.value.size(); ++i) {
      if(mpz_sgn(exp.value[i].getValueTemp()) < 0) {
	use_rat = TRUE;
	break;
      }
    }
    // Don't need bigrational - Joseph Wood 10/2/17
    // if (use_rat) { // a ^ b  with some b negative --> rational result
    //   // 1)  a := as.bigq(a, 1)
    //   SEXP aq = bigrational_as(a, Rf_ScalarInteger(1));
    //   // 2)  result =  <bigq a> ^ b:
    //   return bigrational_pow(aq, b);
    // }
  }
  // else, either, a has a modulus, or (no modulus *and*  exp >= 0) :
  return bigintegerR::biginteger_binary_operation(a,b, pow); // -> pow() in ./bigmod.cc
}
SEXP biginteger_inv (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,inv);}
SEXP biginteger_gcd (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,gcd);}
SEXP biginteger_lcm (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,lcm);}
SEXP biginteger_as (SEXP a, SEXP mod){return bigintegerR::biginteger_binary_operation(a,mod,set_modulus);}
//								set_modulus :  -> ./bigmod.cc

SEXP biginteger_lt  (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::lt);}
SEXP biginteger_gt  (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::gt);}
SEXP biginteger_lte (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::lte);}
SEXP biginteger_gte (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::gte);}
SEXP biginteger_eq  (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::eq);}
SEXP biginteger_neq (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::neq);}



SEXP biginteger_as_character(SEXP a, SEXP b)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans;
  int base = INTEGER(AS_INTEGER(b))[0];
  if (base < 2 || base > 36)
    error(_("select a base between 2 and 36"));

  PROTECT(ans = Rf_allocVector(STRSXP, v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    SET_STRING_ELT(ans, i, Rf_mkChar(v.str(i,base).c_str()));
  // matrix part
  if(v.nrow >= 0)
    {
      SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
      INTEGER(dim)[0] = v.nrow;
      INTEGER(dim)[1] = v.value.size() / v.nrow;
      Rf_setAttrib(ans, Rf_mkString("dim"), dim);
      UNPROTECT(1);
    }

  UNPROTECT(1);
  return ans;
}

SEXP biginteger_as_numeric(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
  double *r = REAL(ans);
  for (unsigned int i = 0; i < v.size(); ++i)
    r[i] = v.value[i].isNA() ? NA_REAL : v.value[i].as_double();
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_as_integer(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(INTSXP,v.size()));
  int *r = INTEGER(ans);
  for (unsigned int i = 0; i < v.size(); ++i) {
    if(v.value[i].isNA()) {
      r[i] = NA_INTEGER;
    }
    else if(!mpz_fits_slong_p(v.value[i].getValueTemp())) {
      Rf_warning("NAs introduced by coercion from big integer");
      r[i] = NA_INTEGER;
    } else {
      r[i] = mpz_get_si(v.value[i].getValueTemp());
    }
  }
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_get_at(SEXP a, SEXP i)
{
  //result = a [i]

  bigvec va = bigintegerR::create_bignum(a);
  return(bigintegerR::create_SEXP(bigintegerR::biginteger_get_at_C(va,i)));

}

// also called from  matrix_get_at_z(.)  in ./extract_matrix.cc :
bigvec bigintegerR::biginteger_get_at_C(bigvec va, SEXP ind)
{
  vector<int> v_ind = bigintegerR::create_int(ind);
  bigvec result;
  // logical: ind = true/false
  if (TYPEOF(ind) == LGLSXP) {
    for (unsigned int i = 0; i < va.size(); ++i)
      if (v_ind[i%v_ind.size()])
	{
	  //std::cout << "cas LOGIC "<< std::endl;
	  result.push_back(va[i]);
	}
    return result;
  }
  else {
    std::remove(v_ind.begin(), v_ind.end(), 0); // remove all zeroes from ind
    if (v_ind.empty())
      return bigvec();

    // case: a[-ind]
    if (v_ind[0] < 0) {
      //std::cout << "cas ngatif" << std::cout;
      for (vector<int>::iterator it = v_ind.begin(); it != v_ind.end(); ++it)
	if (*it > 0)
	  error(_("only 0's may mix with negative subscripts"));
	else if (-(*it)-1 >= (int)va.size())
	  error(_("subscript out of bounds"));

      // TODO: This is optimized for large va.size and small v_ind.size.
      // Maybe add a condition to use a different approach for large v_ind's
      result.value.reserve(va.size()-v_ind.size());
      for (int i = 0; i < (int)va.size(); ++i)
	if (find(v_ind.begin(), v_ind.end(), -i-1) == v_ind.end())
	  {
	    result.push_back(va[i]);
	  }
    }
    else {
      // standard case: a[ind] with ind: integers
      result.value.reserve(v_ind.size());
      for (vector<int>::iterator it = v_ind.begin(); it != v_ind.end(); ++it) {
	if (*it <= 0)
	  error(_("only 0's may mix with negative subscripts"));
	if (*it <= (int)va.size())
	  {
	    //std::cout << "on sort " << va.value[(*it)-1].str(10) << std::endl;
	    result.push_back(va[(*it)-1]);
	  }
	else
	  result.push_back(bigmod()); // NA for out of range's
      }
    }
  }
  return (result);
}

SEXP biginteger_set_at(SEXP src, SEXP idx, SEXP value)
{
  // return = ( src[idx] <- value )

  bigvec result = bigintegerR::create_bignum(src);
  bigvec vvalue = bigintegerR::create_bignum(value);
  vector<int> vidx = bigintegerR::create_int(idx);

  if(vvalue.size() == 0) {
      if(result.size() == 0)
	  return bigintegerR::create_SEXP(result);
      else
	  error(_("replacement has length zero"));
  }
  //case: logicals
  if (TYPEOF(idx) == LGLSXP) {
    int pos = 0;
    for (unsigned int i = 0; i < result.size(); ++i)
      if (vidx[i%vidx.size()])
	result.set(i, vvalue[pos++ % vvalue.size()]);
    return bigintegerR::create_SEXP(result);
  }
  else {
    std::remove(vidx.begin(), vidx.end(), 0); // remove all zeroes
    if (vidx.empty())
      return bigintegerR::create_SEXP(result);
    // return = (src[-idx] <- value)
    if (vidx[0] < 0) {
      for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it)
	if (*it > 0)
	  error(_("only 0's may mix with negative subscripts"));
	else if (-(*it)-1 >= (int)result.size())
	  error(_("subscript out of bounds"));
      int pos = 0;
      for (int i = 0; i < (int)result.size(); ++i)
	if (find(vidx.begin(), vidx.end(), -i-1) == vidx.end())
	  result.set(i, vvalue[pos++%vvalue.size()]);
    }
    //standard case: return = (src[idx] <- value) with idx: positive integer
    else {
      // finding maximum to resize vector if needed
      int maximum = INT_MIN;
      for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it)
	maximum = max(maximum, *it);
      if (maximum > (int)result.size())
	result.resize(maximum);
      int pos = 0;
      for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it) {
	if (*it < 0)
	  error(_("only 0's may mix with negative subscripts"));
	result.set((*it)-1,vvalue[pos++%vvalue.size()]);
      }
    }
  }
  return bigintegerR::create_SEXP(result);
}

SEXP biginteger_length(SEXP a)
{
  return Rf_ScalarInteger(bigintegerR::create_bignum(a).size());
}

SEXP biginteger_setlength(SEXP vec, SEXP value)
{
  int len = 0;
  switch (TYPEOF(value)) {
  case INTSXP:
  case LGLSXP:
    if (LENGTH(value) != 1)
      error(_("invalid second argument"));
    len = *INTEGER(value);
    if (len < 0)
      error(_("vector size cannot be negative"));
    else if (len == NA_INTEGER)
      error(_("vector size cannot be NA"));
    break;
  case REALSXP:
    if (LENGTH(value) != 1)
      error(_("invalid second argument"));
    len = (int)*REAL(value);
    if (len < 0)
      error(_("vector size cannot be negative"));
    else if (! (R_FINITE (len ) ))
      error(_("vector size cannot be NA, NaN of Inf"));
    break;
  case STRSXP:
    // dunno why R spits out this strange error on "length(foo) <- -1"
    // but I always follow the holy standard ;-)
    error(_("negative length vectors are not allowed"));
  default:
    error(_("invalid second argument"));
  }
  bigvec v =bigintegerR::create_bignum(vec);
  v.resize(len);
  return bigintegerR::create_SEXP(v);
}

SEXP biginteger_is_na(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(LGLSXP, v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    LOGICAL(ans)[i] = v[i].value.isNA();
  UNPROTECT(1);
  return ans;
}


SEXP biginteger_sgn(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans = PROTECT(Rf_allocVector(INTSXP, v.size()));
  int *r = INTEGER(ans);
  for (unsigned int i = 0; i < v.size(); ++i)
    r[i] = mpz_sgn(v[i].value.getValueTemp());
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_c(SEXP args)
{
  // if(TYPEOF(args) != VECSXP) error(_("should be a list"));
  bigvec result;
  for(int i=0; i < LENGTH(args); i++) {
    bigvec v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
    for(unsigned int j=0; j < v.size(); j++)
      result.push_back(v[j]);
    v.clear();
  }
  return bigintegerR::create_SEXP(result);
}

SEXP biginteger_cbind(SEXP args)
{
  // if(TYPEOF(args) != VECSXP) error(_("should be a list"));
  bigvec result = bigintegerR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow <= 0)
    result.nrow = result.size();

  for(int i = 1; i < LENGTH(args);i++)
    {
      bigvec v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
      for(unsigned int j=0; j< v.size() ; j++)
	result.push_back(v[j]);
      v.clear();
    }

  return bigintegerR::create_SEXP(result);
}

SEXP biginteger_rep(SEXP x, SEXP times)
{
  bigvec v = bigintegerR::create_bignum(x),
    result;
  int rep = INTEGER(AS_INTEGER(times))[0];

  result.value.reserve(v.size()*rep);
  for(int i = 0 ; i < rep ; i++)
    for(unsigned int j = 0 ; j < v.size() ; j++)
      result.push_back(v[j]);

  return bigintegerR::create_SEXP(result);
}


SEXP biginteger_is_prime(SEXP a, SEXP reps)
{
  bigvec v = bigintegerR::create_bignum(a);
  vector<int> vb = bigintegerR::create_int(reps);
  unsigned int i;
  SEXP ans = PROTECT(Rf_allocVector(INTSXP, v.size()));
  int *r = INTEGER(ans);
  if(v.size() == vb.size())
    for (i = 0; i < v.size(); ++i)
      r[i] = v[i].value.isprime(vb[i]);
  else
    for (i = 0; i < v.size(); ++i)
      r[i] = v[i].value.isprime(vb[0]);
  UNPROTECT(1);
  return ans;
}


SEXP biginteger_nextprime(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a),
    result;
  result.value.reserve(v.size());

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);

  for (unsigned int i = 0; i < v.size(); ++i) {
    mpz_nextprime(val,v[i].value.getValueTemp());
    result.push_back(bigmod(val));
  }
  return bigintegerR::create_SEXP(result);
}

SEXP biginteger_abs(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a),
    result;
  result.value.reserve(v.size());

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);

  for (unsigned int i = 0; i < v.size(); ++i)
    {
      mpz_abs(val,v[i].value.getValueTemp());
      result.push_back(bigmod(val));

      // TODO: understand why following lines don't work.
      //result.push_back(bigmod());
      //result[i].value.setValue(val);
    }

  result.modulus = v.modulus;
  return bigintegerR::create_SEXP(result);
}


/** @brief Bezoult coefficients: compute g,s and t as as + bt = g
 *  @param a BigInteger
 *  @param b BigInteger
 */
SEXP biginteger_gcdex(SEXP a, SEXP b)
{
  bigvec
    va = bigintegerR::create_bignum(a),
    vb = bigintegerR::create_bignum(b),
    result;

  if(va.size() != vb.size())
    return bigintegerR::create_SEXP(bigvec());

  result.value.reserve(3*va.size());

  mpz_t g;
  mpz_t s;
  mpz_t t;
  mpz_init(g);
  mpz_init(s);
  mpz_init(t);
  mpz_t_sentry val_g(g);
  mpz_t_sentry val_s(s);
  mpz_t_sentry val_t(t);

  for(unsigned int i=0; i < va.size(); i++)
    {
      mpz_gcdext (g,s,t,va[i].value.getValueTemp(),vb[i].value.getValueTemp());
      result.value.push_back(biginteger(g)); // Hem... not very elegant !
      result.value.push_back(biginteger(s));
      result.value.push_back(biginteger(t));
      /*      result[i*3].value.setValue(g);
      result[i*3+1].value.setValue(s);
      result[i*3+2].value.setValue(t);*/

    }
  return bigintegerR::create_SEXP(result);
}

/** @brief Random number generation
    \note If seed is not initialised: generation of a new seed
    @param nb  Integer: number of number to generate
    @param length Integer number will be of length 2^length
    @param newseed Integer, seed initialisation (if exists)
    @param ok Integer 1: seed generation 0 not
*/
SEXP biginteger_rand_u (SEXP nb, SEXP length, SEXP newseed, SEXP ok)
{
  mpz_t    bz;
  int i,flag,len,size;
  bigvec result;

  //extern int seed_init;
  //extern gmp_randstate_t seed_state;


  /* store input data into appropriate mode */
  bigvec va = bigintegerR::create_bignum(newseed);
  PROTECT (ok = AS_INTEGER(ok));
  PROTECT (length = AS_INTEGER(length));
  PROTECT (nb = AS_INTEGER(nb));
  flag = INTEGER(ok)[0];
  len = INTEGER(length)[0];
  size = INTEGER(nb)[0];
  UNPROTECT(3);

  result.value.reserve(size);

  /* Random seed initialisation */

  if(seed_init==0)
    {
      gmp_randinit_default(seed_state);
      Rprintf("Seed default initialisation\n");
    }
  if(flag == 1)
    {
      gmp_randseed(seed_state,va[0].value.getValueTemp());
      Rprintf("Seed initialisation\n");
    }

  seed_init = 1;

  mpz_init (bz);
  mpz_t_sentry val_s(bz);

  for(i= 0; i<size; i++)
    {
      /*  Random number generation  */
      mpz_urandomb(bz,seed_state,len);
      result.push_back(bigmod(bz));
    }
  return bigintegerR::create_SEXP(result);
}


/** @brief biginteger_sizeinbase return
 *  @param x BigInteger
 *  @param base BigInteger
 */
SEXP biginteger_sizeinbase(SEXP x, SEXP base)
{
  bigvec vx = bigintegerR::create_bignum(x);
  int basesize= INTEGER(AS_INTEGER(base))[0];
  SEXP ans = PROTECT(Rf_allocVector(INTSXP,vx.size()));
  int *r = INTEGER(ans);
  for(unsigned int i=0; i < vx.size(); i++)
    r[i] = mpz_sizeinbase(vx[i].value.getValueTemp(), basesize);
  UNPROTECT(1);
  return ans;
}


/** @brief bigI_factorial returns n!
 *  @param n non-negative integer vector
 */
SEXP bigI_factorial(SEXP n)
{
  bigvec result;
  int *nn = INTEGER(AS_INTEGER(n)), size = Length(n);
  result.value.resize(size);
  for (int i = 0; i < size; ++i) {
    result.value[i].NA(false);
    if(nn[i] != NA_INTEGER && nn[i] >= 0) {
      mpz_fac_ui(result.value[i].getValue(), (unsigned long int)nn[i]);
    }
  }
  return bigintegerR::create_SEXP(result);
} // bigI_factorial

/** @brief bigI_choose(n, k) returns binomial coefficient (n \choose k)
 *  @param n integer, either R "integer" (non-negative), or a "bigz"
 *  @param k non-negative integer
 */
SEXP bigI_choose(SEXP n, SEXP k)
{
  bigvec result, n_ = bigintegerR::create_bignum(n);
  int *kk = INTEGER(AS_INTEGER(k)), n_k = Length(k);
  int size = (n_.value.empty() || n_k == 0) ? 0 :
    // else:  max(n_.value.size(), n_k)
    (((int)n_.value.size() <= n_k) ? n_k : n_.value.size());

  result.value.resize(size);
  for (int i = 0; i < size; ++i) {
    result.value[i].NA(false);
    int ik_i = kk[i % n_k];
    // check if k in range:
    if(ik_i != NA_INTEGER && ik_i >= 0) {
      unsigned long int k_i = (unsigned long int)ik_i;
      /* void mpz_bin_ui (mpz_t ROP, mpz_t N, unsigned long int K) */
      mpz_bin_ui(result.value[i].getValue(),
		 n_.value[i % n_.value.size()].getValueTemp(), k_i);
    }
  }
  return bigintegerR::create_SEXP(result);
}

/** @brief fibnum return nth Fibonacci number
 *  @param n integer
 */
SEXP bigI_fibnum(SEXP n)
{
  bigvec result;
  if(Length(n) > 0)
    {
      int nn = INTEGER(AS_INTEGER(n))[0];
      unsigned long int num = nn;
      if(nn < 0 || nn == NA_INTEGER)
	  error(_("argument must be non-negative"));
      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);

      mpz_fib_ui(val,num);
      result.push_back(bigmod(val));
      //      result[0].value.setValue(val);
    }
  // else
  //   error(_("argument must not be an empty list"));

  return bigintegerR::create_SEXP(result);

}

/** @brief fibnum2 return nth and n-1th Fibonacci number
 *  @param n integer
 */
SEXP bigI_fibnum2(SEXP n)
{
  bigvec result;
  if(Length(n) > 0)
    {
      int nn = INTEGER(AS_INTEGER(n))[0];
      unsigned long int num = nn;
      if(nn < 0 || nn == NA_INTEGER)
	  error(_("argument must be non-negative"));
      result.value.reserve(1);
      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);
      mpz_t val2;
      mpz_init(val2);
      mpz_t_sentry val_s2(val2);

      mpz_fib2_ui(val,val2, num);
      result.push_back(bigmod(val2));
      result.push_back(bigmod(val));
    }
  else
    error(_("argument must not be an empty list"));

  return bigintegerR::create_SEXP(result);

}

/** @brief lucnum return nth lucas number
 *  @param n integer
 */
SEXP bigI_lucnum(SEXP n)
{
  bigvec result;
  if(Length(n) > 0)
    {
      int nn = INTEGER(AS_INTEGER(n))[0];
      unsigned long int num = nn;
      if(nn < 0 || nn == NA_INTEGER)
	  error(_("argument must be non-negative"));

      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);

      mpz_lucnum_ui(val,num);
      result.push_back(bigmod(val));
    }
  // else
  //   error(_("argument must not be an empty list"));

  return bigintegerR::create_SEXP(result);

}

/** @brief lucnum2 return nth and n-1th lucas number
 *  @param n integer
 */
SEXP bigI_lucnum2(SEXP n)
{
  bigvec result;

  if(Length(n) > 0) {
    int nn = INTEGER(AS_INTEGER(n))[0];
    unsigned long int num = nn;
    if(nn < 0 || nn == NA_INTEGER)
      error(_("argument must be non-negative"));
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);
    mpz_t val2;
    mpz_init(val2);
    mpz_t_sentry val_s2(val2);

    mpz_lucnum2_ui(val,val2,num);
    result.push_back(bigmod(val2));
    result.push_back(bigmod(val));
  }
  else
    error(_("argument must not be an empty list"));

  return bigintegerR::create_SEXP(result);
}



//Return max

SEXP biginteger_max(SEXP a, SEXP narm)
{
  bigvec result;
  bigvec va = bigintegerR::create_bignum(a);

  if( ! va.size())
    return bigintegerR::create_SEXP(result);

  unsigned int maximum = 0;
  int na_remove = Rf_asInteger(narm);

  for(unsigned int i = 1 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() && !na_remove)
	return(bigintegerR::create_SEXP(result));
      else
	if(!(va.value[i] <  va.value[maximum] ))
	  maximum = i; // if va.value[maximum = 0] is NA => false for the "<" => maximum changed = good
    }

  result.push_back(va.value[maximum]);


  // now the modulus !
  if(va.modulus.size() == 1)
    result.modulus.push_back(va.modulus[0]);

  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());
      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i)
	{
	  if(modulus != va.modulus[i]) // if one is different: no modulus
	    return(bigintegerR::create_SEXP(result));
	}
      result.modulus.push_back(modulus);
    }

  return(bigintegerR::create_SEXP(result));
}


// Return min
SEXP biginteger_min(SEXP a, SEXP narm)
{
  bigvec result;
  bigvec va = bigintegerR::create_bignum(a);

  if( ! va.size())
    return bigintegerR::create_SEXP(result);

  unsigned int minimum = 0;
  int na_remove = Rf_asInteger(narm);

  for(unsigned int i = 1 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() && !na_remove)
	return bigintegerR::create_SEXP(result);
      else
	if(!(va.value[i] >  va.value[minimum] ))
	  minimum = i;
    }

  result.push_back(va.value[minimum]);

  // now the modulus !
  if(va.modulus.size() == 1)
    result.modulus.push_back(va.modulus[0]);

  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());
      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i)
	{
	  if(modulus != va.modulus[i]) // if one is different: no modulus
	    return bigintegerR::create_SEXP(result);
	}
      result.modulus.push_back(modulus);
    }

  return bigintegerR::create_SEXP(result);
}


SEXP biginteger_cumsum(SEXP a)
{
  bigvec result, va = bigintegerR::create_bignum(a);

  result.value.resize(va.value.size());

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);

  bool hasmodulus = true;

  // first the modulus !
  if(va.modulus.size() > 1) {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());

      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i) {
	    if(modulus != va.modulus[i]) { // if one is different: no modulus
		hasmodulus = false;
		break;
	    }
      }
      if(hasmodulus)
	result.modulus.push_back(modulus);
  }
  else if(va.modulus.size() == 1) {
      result.modulus.push_back(va.modulus[0]);
      hasmodulus = true;
  }
  else hasmodulus = false;

  for(unsigned int i = 0 ; i < va.size(); ++i)
    {
      {
	if(va.value[i].isNA() )
	  {
	    break; // all last values are NA.
	  }

	mpz_add(val,val,va.value[i].getValueTemp());

	if(hasmodulus)
	  mpz_mod(val,val,va.modulus[0].getValueTemp() );

	result.value[i].setValue(val);
      }
    }

  return(bigintegerR::create_SEXP(result));
}



SEXP biginteger_sum(SEXP a)
{
  bigvec result, va = bigintegerR::create_bignum(a);

  result.value.resize(1);

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);

  bool hasmodulus = true;

  // first the modulus !
  if(va.modulus.size() > 1) {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());

      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i) {
	    if(modulus != va.modulus[i]) { // if one is different: no modulus
		hasmodulus = false;
		break;
	    }
      }
      if(hasmodulus)
	result.modulus.push_back(modulus);
  }
  else if(va.modulus.size() == 1) {
      result.modulus.push_back(va.modulus[0]);
      hasmodulus = true;
  }
  else
      hasmodulus = false;





  for(unsigned int i = 0 ; i < va.size(); ++i)
    {
      {
	if(va.value[i].isNA() )
	  {
	    break; // all last values are NA.
	  }

	mpz_add(val,val,va.value[i].getValueTemp());

	if(hasmodulus)
	  mpz_mod(val,val,va.modulus[0].getValueTemp() );
      }
    }

  result.value[0].setValue(val);

  return(bigintegerR::create_SEXP(result));
}


SEXP biginteger_prod(SEXP a)
{
  bigvec result;
  bigvec va = bigintegerR::create_bignum(a);

  result.value.resize(1);

  mpz_t val;
  mpz_init(val);
  mpz_set_ui(val,1);
  mpz_t_sentry val_s(val);

  bool hasmodulus = true;

  // first the modulus !
  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());

      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i) {
	if(modulus != va.modulus[i]) { // if one is different: no modulus
	  hasmodulus = false;
	  break;
	}
      }
      if(hasmodulus)
	result.modulus.push_back(modulus);

    }
  else
    hasmodulus = false;

  if(va.modulus.size() == 1)
    {
      result.modulus.push_back(va.modulus[0]);
      hasmodulus = true;
    }

  for(unsigned int i = 0 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() )
	{
	  return (bigintegerR::create_SEXP(result));
	}

      mpz_mul(val,val,va.value[i].getValueTemp());

      if(hasmodulus)
	mpz_mod(val,val,va.modulus[0].getValueTemp() );
    }

  result.value[0].setValue(val);

  return(bigintegerR::create_SEXP(result));
}

// return x ^ y [n]
SEXP biginteger_powm(SEXP x, SEXP y, SEXP n)
{
  bigvec result;

  bigvec vx = bigintegerR::create_bignum(x);
  bigvec vy = bigintegerR::create_bignum(y);
  bigvec vn = bigintegerR::create_bignum(n);

  result.value.resize(vx.value.size());

  for (unsigned int i = 0 ; i < vx.value.size(); i++)
  {

    result.value[i].NA(false);
    // check if n != 0
    if(mpz_sgn(vn.value[i % vn.value.size()].getValueTemp()) != 0)
      mpz_powm(result.value[i].getValue(),
	       vx.value[i].getValueTemp(),
	       vy.value[i % vy.value.size()].getValueTemp(),
	       vn.value[i % vn.value.size()].getValueTemp());

  }

  return bigintegerR::create_SEXP(result);
} // ..._powm()


// TODO:  A version that only returns  'ex'  {i.e. the binary precision}
// ----   GMP manual suggests that     size_t mpz_sizeinbase (mpz_t OP, int BASE)
//        i.e.,  mpz_sizeinbase (OP, 2)  would give that

// (d, ex) where  x = d * 2^ex, and  0.5 <= |d| < 1
SEXP bigI_frexp(SEXP x)
{
// double mpz_get_d_2exp (signed long int *exp, mpz t op )

// Convert op to a double, truncating if necessary (ie. rounding towards zero), and returning
// the exponent separately.
// The return value is in the range 0.5 ≤ |d| < 1 and the exponent is stored to *exp . d ∗ 2exp is
// the (truncated) op value. If op is zero, the return is 0.0 and 0 is stored to *exp .
// This is similar to the standard C frexp function.

  const char *nms[] = {"d", "exp", ""};
  SEXP ans, d_R, exp_R;
  bigvec vx = bigintegerR::create_bignum(x);
  int n = vx.value.size();

  PROTECT(ans = Rf_mkNamed(VECSXP, nms)); // =~= R  list(d = . , exp= .)
  d_R   = Rf_allocVector(REALSXP, n); SET_VECTOR_ELT(ans, 0, d_R);
  exp_R = Rf_allocVector(INTSXP,  n); SET_VECTOR_ELT(ans, 1, exp_R);
  double *d_ = REAL(d_R);
  int  *exp_ = INTEGER(exp_R);
  for (int i = 0 ; i < n; i++) {
    signed long int ex;
    d_[i] = mpz_get_d_2exp(&ex, vx.value[i].getValueTemp());
    if(abs(ex) < INT_MAX)
      exp_[i] = (int) ex;
    else
      error(_("exponent too large to fit into an integer"));
  }
  UNPROTECT(1);
  return ans;
} // bigI_frexp()


SEXP biginteger_log2(SEXP x)
{
  bigvec v = bigintegerR::create_bignum(x);
  SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
  double *r = REAL(ans);

  for (unsigned int i = 0; i < v.size(); ++i) {
    signed long int ex;
    double di = mpz_get_d_2exp(&ex, v.value[i].getValueTemp());
    // xi = di * 2 ^ ex  ==> log2(xi) = log2(di) + ex :
    r[i] = log(di) / M_LN2 + (double)ex;
  }
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_log(SEXP x)
{
  bigvec v = bigintegerR::create_bignum(x);
  SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
  double *r = REAL(ans);

  for (unsigned int i = 0; i < v.size(); ++i) {
    signed long int ex;
    double di = mpz_get_d_2exp(&ex, v.value[i].getValueTemp());
    // xi = di * 2 ^ ex  ==> log(xi) = log(di) + ex*log(2) :
    r[i] = log(di) + M_LN2*(double)ex;
  }
  UNPROTECT(1);
  return ans;
}
