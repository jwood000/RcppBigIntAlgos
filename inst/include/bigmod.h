/*! \file bigmod.h
 *  \brief Description of class bigmod
 *
 *  \date Created: 22/05/06
 *  \date Last modified: Time-stamp: <2014-07-10 08:44:04 antoine>
 *
 *  \author Immanuel Scholz
 *
 *  \note Licence: GPL
 */


#ifndef BIGMOD_HEADER_
#define BIGMOD_HEADER_ 1

#include "biginteger.h"

typedef void (*gmp_binary)(mpz_t, const mpz_t, const mpz_t);

extern "C" {
  /**
   * division
   * result = a / b
   */
  void integer_div(mpz_t result,const mpz_t a, const mpz_t b);
}


/**
 * \brief class for bigmod values. Represent any integer in Z/nZ
 *
 * Represents two biginteger: a value and a modulus. These both are used
 * to operate arithmetic functions on it. If the modulus is NA, no modulus
 * to the operation result is applied. If the value is NA, the result is always NA.
 */
class bigmod {
 public:
  /** \brief  Value of our bigmod */
  biginteger value;
  /** \brief  modulus of our bigmod representation: value %% modulus */
  biginteger modulus;

  /** \brief creator
   */
  bigmod(const biginteger& value_ = biginteger(),
	 const biginteger& modulus_ = biginteger()) :
    value(value_),modulus(modulus_) {}

  /** \brief copy operator  */
  bigmod(const bigmod & rhs) :
    value(rhs.value),
    modulus(rhs.modulus){}

  /**
   * \brief  Return as a human readible string
   */
  std::string str(int b) const;

  /** \brief assignement operator */
  bigmod & operator= (const bigmod& rhs);

  /** \brief return sign (-1 if negative, 0 if 0; +1 if positive)
   */
  int sgn() const
    {
      return(mpz_sgn(value.getValueTemp()));
    }

  bigmod inv () const;


};


/** \brief comparison operator
 */
bool operator!= (const bigmod& rhs, const bigmod& lhs);

/** \brief comparison operator
 */
bool operator== (const bigmod& rhs, const bigmod& lhs);



/**
 * \brief Add two bigmods together.
 *
 * If only one has a modulus set, the result will have this
 * modulus. If both bigmods disagree with the modulus, the result will not have
 * a modulus set. If none modulus for either bigmod is set, the result will not
 * have a modulus as well.
 */
bigmod operator+(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Subtract two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator-(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Multiply two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator*(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Divide two bigmods   a / b  :=  a * b^(-1)
 */
bigmod div_via_inv(const bigmod& a, const bigmod& b);

/**
 * \brief Divide two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator/(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Calculate the modulus (remainder) of two bigmods.
 *
 * The resulting bigmod will have set the intern modulus to
 * the value of lhs, no matter what rhs.modulus or lhs.modulus
 * was before, except if rhs and lhs has both no modulus set,
 * in which case the resulting modulus will be unset too.
 */
bigmod operator%(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Return the power of "exp" to the base of "base" (return = base^exp).
 *
 * If both moduli are unset or unequal, this may EAT your memory alive,
 * since then the infinite "pow" is used instead of the modulus "powm".
 * You  may not try to pow a value this way with an exponent that does
 * not fit into a long value.
 *
 * For other modulus description, see operator+(bigmod, bigmod)
 */
bigmod pow(const bigmod& base, const bigmod& exp);

/**
 * \brief Return the modulo inverse to x mod m. (return = x^-1 % m)
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod inv(const bigmod& x, const bigmod& m);

/**
 * \brief Return a bigmod with value (x % m) and the intern modulus set to m.
 * Intern modulus settings of x and m are ignored.
 *
 * Do not confuse this with operator%(bigmod, bigmod).
 */
bigmod set_modulus(const bigmod& x, const bigmod& m);


biginteger get_modulus(const bigmod& b1, const bigmod& b2);
/**
 * \brief Return the greatest common divisor of both parameters
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod gcd(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief  Return the least common multiply of both parameter.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod lcm(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief function used to make any binary operation between
 * two bigmod that return a bigmod (addition substraction... )
 */
bigmod create_bigmod(const bigmod& lhs, const bigmod& rhs, gmp_binary f,
		     bool zeroRhsAllowed = true) ;

#endif
