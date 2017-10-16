
#include "bigvec.h"

/** \brief constructor
 *
 */
bigvec::bigvec(unsigned int i) :
  value(i),
  modulus(0),
  nrow(-1)
{
}


bigvec::bigvec(const bigvec & vecteur) :
  value(vecteur.value.size()),
  modulus(vecteur.modulus.size()),
  nrow(vecteur.nrow)
{
  //  *this = vecteur;
  value.resize(vecteur.value.size());
  modulus.resize(vecteur.modulus.size());
  std::vector<biginteger>::const_iterator jt=vecteur.modulus.begin();
  std::vector<biginteger>::iterator it = modulus.begin();
  while(it != modulus.end())
    {
      *it = *jt;
      ++it;
      ++jt;
    }
  jt = vecteur.value.begin();
  for(it=value.begin(); it != value.end(); ++it)
    {
      *it= *jt;
      ++jt;
    }
}


//
std::string bigvec::str(int i,int b) const
{
    if (value[i].isNA())
      return ("NA");

    std::string s; // sstream seems to collide with libgmp :-(
    bool onemodulus = modulus.size()>0;
    if(onemodulus)
      onemodulus = !modulus[i%modulus.size()].isNA() ;
    if (onemodulus)
	s = "(";
    s += value[i].str(b);
    if (onemodulus) {
	s += " %% ";
	s += modulus[i%modulus.size()].str(b);
	s += ")";
    }

    return s;
}


bigmod bigvec::operator[] (unsigned int i) const
{
  if(modulus.size()>0)
    return(bigmod(value[i], modulus[i%modulus.size()]));
  return(bigmod(value[i]));
}

void bigvec::set(unsigned int i,const bigmod & val)
{

  value[i] = val.value;

  if(!val.modulus.isNA())
    {
      if(modulus.size()> i)
	{
	  modulus[i] = val.modulus ;
	  return;
	}

      int nrow_mod = nrow;
      if(nrow<1)
	nrow_mod = 1;
      if( (modulus.size() ==  (unsigned int)nrow_mod ) || (modulus.size() == 1) )
	{
	  // check "row-mod" or "global" modulus
	  if(!(val.modulus != modulus[i % modulus.size()]))
	    return;
	}

      // we add "previous"
      nrow_mod = modulus.size();
      for(unsigned int j = nrow_mod; j< i;++j)
	modulus.push_back(modulus[j % nrow_mod]);
      modulus.push_back(val.modulus);

    }
}

void bigvec::push_back(const bigmod & number)
{
  int nrow_mod = (nrow < 0) ? 1 : nrow;

  value.push_back(number.value);

  if((!number.modulus.isNA()) || (modulus.size()>0) )
    {
      // pathologic case: value.size >0 ; modulus.size =0
      // we assume previous mod where NA
      if((modulus.size() ==0) && (value.size()>0))
	{
	  modulus.resize(value.size()-1);
	  modulus.push_back( number.modulus);
	  return;
	}

      // standard cas
      if((modulus.size() != 1 ) && (static_cast<int>(modulus.size()) != nrow_mod) )
	{
	  modulus.push_back(number.modulus);
	  return;
	}

      // pathologic case:
      //  value modulus [nrow=4]
      //  2     2  push_back: ok
      //  2     2  push_back: nothing
      //  2     1  push_back: shoud add previous modulus then 1
      // note nrow_mod is either 1 ither nrow (when nrow >1)
      nrow_mod = modulus.size();
      if(modulus[(value.size() -1)% nrow_mod ] != number.modulus)
	{
	  // we add "previous"
	  for(unsigned int i = nrow_mod; i< value.size()-1;i++)
	    modulus.push_back(modulus[i % nrow_mod]);
	  modulus.push_back(number.modulus);
	  }
    }
}

/** 
 * insert int value
 */
void bigvec::push_back(int value_p)
{
  value.push_back(biginteger(value_p));
}

/**
 * Insert Big Integer value
 */
void bigvec::push_back(const __mpz_struct * value_p)
{
  value.push_back(biginteger(value_p));
}


// return size of value
unsigned int bigvec::size() const
{
  return(value.size());
}

// hummm. to avoid !
void bigvec::resize(unsigned int i)
{

  value.resize(i);
  if(i < modulus.size())
      modulus.resize(i);
}

// clear all
void bigvec::clear()
{
  value.clear();
  modulus.clear();
  nrow = -1;
}


// assignment operator
bigvec & bigvec::operator= (const bigvec & rhs)
{
  if(this != &rhs)
    {
      value.resize(rhs.value.size());
      modulus.resize(rhs.modulus.size());
      std::vector<biginteger>::const_iterator jt=rhs.modulus.begin();
      std::vector<biginteger>::iterator it = modulus.begin();
      while(it != modulus.end())
	{
	  *it = *jt;
	  ++it;
	  ++jt;
	}
      jt = rhs.value.begin();
      for(it=value.begin(); it != value.end(); ++it)
	{
	  *it= *jt;
	  ++jt;
	}
      nrow = rhs.nrow;
    }
  return(*this);
}



// Comparison operator
bool operator!=(const bigvec & rhs, const bigvec& lhs)
{
  std::vector<biginteger>::const_iterator it = rhs.value.begin();
  std::vector<biginteger>::const_iterator jt = lhs.value.begin();

  if( (rhs.value.size() != lhs.value.size()) || \
      (rhs.nrow != lhs.nrow )  )
    return(false);

  // check value
  while(it != rhs.value.end())
    {
      if(*it != *jt)
	return(false);
      it++; jt++;
    }
  for(unsigned int i=0; i < std::max( rhs.modulus.size() ,lhs.modulus.size() ); ++i)
    if(rhs.modulus[i % rhs.modulus.size()] != \
       lhs.modulus[i % lhs.modulus.size()] )
      return(false);
  return(true);
}



//
// \brief substract lambda[0] * line j to line i
//
void bigvec::subLine(unsigned int i,unsigned int j,const bigvec & lambda)
{
  if(nrow <= 0)
    error(_("Need matrix with at least one row to do this operation"));

  unsigned int k, n = value.size() / (unsigned int) nrow;
  if(modulus.size() != 1)
    {
      for(k=0; k < n; ++k)
	value[i + k*nrow] =  value[i + k*nrow] - ( value[j + k*nrow] * lambda.value[0] ) ;
    }
  else
     for(k=0; k < n ; ++k)
       {
	 value[i + k*nrow] =  value[i + k*nrow] - ( value[j + k*nrow] * lambda.value[0] ) ;
	 value[i + k*nrow] =  value[i + k*nrow] % modulus[0];
       }
}

/*
 * \brief multiply line i by lambda
 */
void bigvec::mulLine(unsigned int i,const bigvec & lambda)
{
  if(nrow <= 0)
    error(_("Need matrix with at least one row to do this operation"));

  unsigned int k;
  // n number of columns
  unsigned int n = value.size() / (unsigned int) nrow;
  if(modulus.size() != 1)
    {
      for(k=0; k < n; ++k)
	value[i + k*nrow] =  value[i + k*nrow]  * lambda.value[0]  ;
    }
  else
     for(k=0; k < n ; ++k)
       {
	 value[i + k*nrow] =  value[i + k*nrow] * lambda.value [0] ;
	 value[i + k*nrow] =  value[i + k*nrow] % modulus[0];
       }
}


// never used
void bigvec::print()
{
  if(nrow > 0) {
    for(int i=0; i < nrow; ++i)
      {
	for(int j=0; j < (value.size() / nrow); ++j)
	  Rprintf("%s\t", value[i+j* nrow].str(10).c_str() );
	Rprintf("\n");
      }
  }
  else {
    for(unsigned int i=0; i < value.size(); ++i)
      Rprintf("%s\t", value[i].str(10).c_str() );
    Rprintf("\n");
  }
}
