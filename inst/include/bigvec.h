/*! \file bigvec.h
 *  \brief bigvec class definition
 *
 *  \version 1
 *
 *  \date Created: 2005
 *  \date Last modified: Time-stamp: <2013-06-08 22:27:37 antoine>
 *
 *
 *  \note Licence: GPL
 */


#ifndef BIGVEC_HEADER_
#define BIGVEC_HEADER_ 1

#include "bigmod.h"


/** \brief class bigvec
 *
 * It a class composed of 2 vectors, (value & modulus) that
 * can be of different size and a nrow
 * parameter (for matrix support)
 */
class bigvec {
 public:
  /** \brief value */
  std::vector<biginteger> value;
  /** \brief modulus */
  std::vector<biginteger> modulus;
  /** \brief optional parameter used with matrix -- set to -1 for non-matrix */
  int nrow ;

  /** \brief initialize value to size i
   */
  bigvec(unsigned int i = 0);

  /**
   * \brief copy constructor
   */
  bigvec(const bigvec & vecteur);

  /**
   * \brief construct a bigmod at indice i
   *
   * It gets values from value[i] & modulus[i % modulus.size()]
   *
   * \note should  not used for assignement
   */
  bigmod  operator[] (unsigned int i) const;


  /**
   * \brief assign a value at indice i
   */
  void set(unsigned int i,const bigmod & val);

  /**
   * \brief extend our vectors.
   *
   * This function will check if modulus should be added.
   * Modulus can be set "globaly" i.e. one modulus for the
   * whole vector/matrix.
   * Or set by row (a constant modulus for each row)
   * Or set by cell (as many modulus as value)
   */
  void push_back(const bigmod &i);

  /** 
   * insert int value
   */
  void push_back(int value_p);

  /**
   * Insert Big Integer value
   */
  //  void push_back(mpz_t & value_p);
  void push_back(const __mpz_struct* value_p);

  /**
   * \brief return size of vector value
   */
  unsigned int size() const ;

  /**
   * \brief extend vector value.
   */
  void resize(unsigned int i);

  /**
   * \brief clear all
   */
  void clear();

  /**
   * \brief Return as a human readible string
   */
  std::string str(int i, int b) const;

  /**
   * \brief assignement operator
   */
  bigvec & operator= (const bigvec & rhs);


  /**
   * \brief substract lambda * line j to line i
   */
  void subLine(unsigned int i,unsigned int j,const bigvec & lambda);

  /**
   * \brief multiply line i by lambda
   */
  void mulLine(unsigned int i,const bigvec & lambda);

  /** \brief print matrix to stdout
   *
   * use for debug purpose
   */
  void print();



};
//


/** \brief comparison operator
 */
bool operator!= (const bigvec& rhs, const bigvec& lhs);



#endif
