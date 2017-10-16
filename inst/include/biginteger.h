/*! \file biginteger.h
 *  \brief Description of class biginteger
 *
 *  \date Created: 2004
 *  \date Last modified: Time-stamp: <2010-04-10 19:09:32 antoine>
 *
 *  \author Immanuel Scholz
 *
 *  \note Licence: GPL
 */

#ifndef biginteger_HEADER
#define biginteger_HEADER 1

#include <string>
#include <vector>

#include "Rgmp.h"

/** \mainpage Gmp package for R language.
 *
 * Theses pages made to help developpers to enter into the C++ code related to
 * the gmp package.
 *
 * \section sec_intro Introduction
 *
 * The gmp R package uses directly the gmp C/C++ library (http://gmplib.org)
 * to provide powerful computation of big integers and rational within R.
 * This package is more than a simple wrapper as it allows to handle
 *    - easy computation in Z/nZ (it is an option of bigz class)
 *    - matrix computation
 *
 */


/** \brief mpz struct
 *
 * Use this to clear mpz_t structs at end-of-function automatically
 */
struct mpz_t_sentry {
  /** \brief value */
  mpz_t& value;
  /** \brief constructor */
  mpz_t_sentry(mpz_t& v): value(v) {}
  /** \brief detructor (and clear mpz) */
  ~mpz_t_sentry() {mpz_clear(value);}
};


/** \brief Class biginteger
 *
 * A big integer. Actually a wrapper for mpz_t to work with plus
 * some special stuff.
 *
 * The biginteger special state "NA" means, no value is assigned.
 * This does not mean, the internal state is not constructed, but
 * the value explicit is "not available".
 */





class biginteger
{
 private:
  /**
   * The actual integer value.
   */
  mpz_t value;

  /**
   * True, if the value is "NA".
   */
  bool na;

 public:
  /**
   * Construct a "NA" biginteger.
   */
  biginteger() : na(true) {mpz_init(value);}

  /**
   * Construct a biginteger from a raw expression.
   */
  biginteger(const char* raw);

  /**
   * Create a biginteger from a value. Remember to free the
   * parameter's mpz_t if you allocated them by yourself -
   * biginteger will copy the value.
   */
  biginteger(const mpz_t value_);

  /**
   * Construct a biginteger from a int value.
   */
  biginteger(const int value_) : na(false) {
    if(value_ ==  NA_INTEGER)
      {mpz_init(value); na = true  ;}
    else
      mpz_init_set_si(value, value_);}

  /**
   * Construct a biginteger from a long value.
   */
  biginteger(const long int value_) : na(false) {
    if(value_ ==  NA_INTEGER)
      {mpz_init(value); na = true  ;}
    else
      mpz_init_set_si(value, value_);}

  /**
   * Construct a biginteger from a unsigned long value.
   */
  biginteger(const unsigned long int value_) : na(false) {
    if(value_ == (unsigned long int) NA_INTEGER)
      {mpz_init(value); na = true  ;}
    else
      mpz_init_set_ui(value, value_);}

  /**
   * Construct a biginteger from a double value.
   */
  biginteger(const double value_) : na(false) {
    if( R_FINITE(value_) )
      mpz_init_set_d(value, value_);
    else
      {mpz_init(value); na = true  ;}
  }

  /**
   * Construct a biginteger from a string value.
   */
  biginteger(const std::string& value_) : na(false)
    {
      /* mpz_init.. return -1 when error, 0: ok */
      if(mpz_init_set_str(value, value_.c_str(), 0))
	{
	  mpz_set_si(value, 0);
	  na=true;
	}
      /*	if(mpz_init_set_str(value, value_.c_str(), 0) == -1)
		Rf_error("Not a valid number");    */
    }

  /**
   *  Copy constructor (mpz_t aren't standard-copyable)
   */
  biginteger(const biginteger& rhs) : na(rhs.na)
    {
      mpz_init_set(value, rhs.getValueTemp());
    }


  /**
   * Free the owned mpz_t structs
   */
  virtual ~biginteger() {mpz_clear(value);}


  /**
   * Set the biginteger to state "NA".
   */
  void setValue() {mpz_set_si(value, 0); na = true;}

  /**
   * Set the biginteger to a specific value.
   */
  void setValue(const mpz_t & value_ ) {
    mpz_set(value, value_); na = false;
  }

  /** \brief set value from an integer
   */
  void setValue(int value_) {
    if(value_ == NA_INTEGER)
      {mpz_set_ui(value, 0); na = true  ;}
    else
      {
	mpz_set_si(value, value_);
	na = false;
      }
  }
  /** \brief set value from a long integer
   */
  void setValue(long int value_) {
    if(value_ == NA_INTEGER)
      {mpz_set_ui(value, 0); na = true  ;}
    else
      {
	mpz_set_si(value, value_);
	na = false;
      }
  }

  /** \brief set value from an unsigned int
   */
  void setValue(unsigned long int value_) {
    if((int)value_ == NA_INTEGER)
      {mpz_set_ui(value, 0); na = true  ;}
    else
      {
	mpz_set_ui(value, value_);
	na = false;
      }
  }

  /** \brief set value from a float
   */
  void setValue(double value_) {
    if(R_FINITE (value_) )
      {mpz_set_d(value, value_); na = false;}
    else
      {mpz_set_ui(value, 0); na = true  ;}
  }


  /** \brief set value from a biginteger
   */
  void setValue(const biginteger  & value_) {
    setValue(value_.getValueTemp());
    na = value_.isNA();

  }
  /**
   * For const-purposes, return the value. Remember, that the return value
   * only lives as long as this class live, so do not call getValueTemp on
   * temporary objects.
   */
  const mpz_t& getValueTemp() const {return value;}


  /** \brief accessor on value
   */
  mpz_t & getValue()
  {
    return value;
  }

  /**
   * Return true, if the value is NA.
   */
  bool isNA() const {return na;}

  /**
   * set NA value
   */
  void NA(bool value_p)  {na = value_p;}

  /**
   * Return true, if the value is 0.
   */
  int sgn() const {return mpz_sgn(value);}

  /**
   *  Convert the biginteger into a standard string.
   */
  std::string str(int b) const;

  /**
   * Convert the biginteger into a long value (cut off the msb's if it don't
   * fit).
   */
  long as_long() const {return mpz_get_si(value);}

  /**
   * \brief Convert the biginteger into a double value
   * (you may loose precision)
   */
  double as_double() const {return mpz_get_d(value);}

  /**
   * Convert the biginteger to a raw memory block. Obtain the size needed
   * from biginteger_raw_size() first and make sure, the buffer provided is
   * large enough to hold the data.
   *
   * Also remember, that the modulus is not saved this way. To obtain a
   * modulus raw byte use get_modulus().as_raw(void*).
   *
   * @return number of bytes used (same as raw_size())
   */
  int as_raw(char* raw) const;

  /**
   * Return the number of bytes needed to store this biginteger in a
   * continous memory block.
   */
  size_t raw_size() const;

  /**
   * Swap values with the argument
   */
  void swap(biginteger& other);

  /**
   * Test prime numbers
   */
  int isprime(int reps){return  mpz_probab_prime_p(value,reps);}


  /** \brief overload affectation operator
   */
  biginteger & operator= (const biginteger& rhs);

};




/** \brief comparison operator
 */
bool operator!= (const biginteger& rhs, const biginteger& lhs);


/** \brief comparison operator
 */
bool operator> (const biginteger& rhs, const biginteger& lhs);

/** \brief comparison operator
 */
bool operator< (const biginteger& rhs, const biginteger& lhs);

/** \brief standard operator */
biginteger operator* (const biginteger& rhs, const biginteger& lhs);

/** \brief standard operator */
biginteger operator- (const biginteger& rhs, const biginteger& lhs);

/** \brief standard operator */
biginteger operator% (const biginteger& rhs, const biginteger& lhs);

/** \brief function used by rational that should export
 *   mpz value
 */
int as_raw(char* raw,mpz_t value,bool na);




#endif
