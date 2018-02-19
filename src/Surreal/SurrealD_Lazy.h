#ifndef SURREALD_LAZY_H
#define SURREALD_LAZY_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <cassert>
#include <iostream>

#include <algorithm>

#include "always_inline.h"

class SurrealDTypeBase {};

//A class used to distinguish surreal expressions
template< class Derived >
struct SurrealDType : SurrealDTypeBase
{
  //A convenient method for casting to the derived type
  ALWAYS_INLINE const Derived& cast() const { return static_cast<const Derived&>(*this); }

  //A simple way to call value without having to case first
  ALWAYS_INLINE double value() const { return cast().value(); }
};

namespace SurrealDExp
{
template<class L, class R > class OpMul;
class OpMul_impl;
}

//----------------------------------------------------------------------------//
// SurrealD:  value, N derivatives
//
// Operators with Return Value Optimization (RVO)
//
//----------------------------------------------------------------------------//

class SurrealD : public SurrealDType< SurrealD >
{
public:
  //The default constructor is intentionally not included here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  SurrealD( const SurrealD& z );
  explicit SurrealD( int n=0 );
  SurrealD( const double& v0, const double d0[], int n );
  SurrealD( const double& v0, const double& d0, int n );
  SurrealD( const double& v0, int n );
  SurrealD( const double& v0 );
  template<class Expr>
  ALWAYS_INLINE
  SurrealD( const SurrealDType<Expr>& r ) : v_(0), d_(NULL), N_(0) { operator=(r); }
  ALWAYS_INLINE ~SurrealD();

#if 0     // these are not needed and will be removed
  operator double() const { return v_; }
  operator int() const { return int(v_); }
#endif

  ALWAYS_INLINE int size() const { return N_; }

  // value accessor operators
  ALWAYS_INLINE       double& value()       { return v_; }
  ALWAYS_INLINE const double& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE double& deriv( int i=0 )       { assert(d_ != NULL); return d_[i]; }
  ALWAYS_INLINE double  deriv( int i=0 ) const { return N_ > 0 ? d_[i] : 0; }

  // assignment
  SurrealD& operator=( const SurrealD& );
  SurrealD& operator=( const double& );

  template<class Expr> SurrealD& operator= ( const SurrealDType<Expr>& );
  template<class Expr> SurrealD& operator+=( const SurrealDType<Expr>& );
  template<class Expr> SurrealD& operator-=( const SurrealDType<Expr>& );


  // unary operators; no side effects
  const SurrealD& operator+() const;
  const SurrealDExp::OpMul<SurrealD, double>  operator-() const;

  // binary accumulation operators
  SurrealD& operator+=( const SurrealD& );
  SurrealD& operator+=( const double& );
  SurrealD& operator+=( const int& );
  SurrealD& operator-=( const SurrealD& );
  SurrealD& operator-=( const double& );
  SurrealD& operator-=( const int& );
  SurrealD& operator*=( const SurrealD& );
  SurrealD& operator*=( const double& );
  SurrealD& operator*=( const int& );
  SurrealD& operator/=( const SurrealD& );
  SurrealD& operator/=( const double& );
  SurrealD& operator/=( const int& );

  // misc functions <cmath>
  friend SurrealDExp::OpMul<SurrealD, double> abs( const SurrealD& );
  friend SurrealDExp::OpMul<SurrealD, double> fabs( const SurrealD& );

  // classification functions <cmath>
  friend bool isfinite( const SurrealD& );
  friend bool isinf( const SurrealD& );
  friend bool isnan( const SurrealD& );

  // input
  friend std::istream& operator>>( std::istream&, SurrealD& );

  void dump( int indentSize=0, std::ostream& out = std::cout ) const;

private:
  double v_;          // value
  double *d_;         // derivative array
  int N_;           // size of derivative array
};


// constructors

ALWAYS_INLINE
SurrealD::SurrealD( const SurrealD& z ) : v_(z.v_), d_(NULL), N_(z.N_)
{
  if (N_ > 0)
  {
    d_ = new double[N_];
    for (int i = 0; i < N_; i++)
      d_[i] = z.d_[i];
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( int n ) : v_(n), d_(NULL), N_(0)
{
//  if (N_ > 0)
//  {
//    d_ = new double[N_];
//    for (int i = 0; i < N_; i++)
//      d_[i] = 0;
//  }
}

ALWAYS_INLINE
SurrealD::SurrealD( const double& v0, const double d0[], int n ) : v_(v0), d_(NULL), N_(n)
{
  assert( N_ > 0 );

  d_ = new double[N_];
  for (int i = 0; i < N_; i++)
    d_[i] = d0[i];
}

ALWAYS_INLINE
SurrealD::SurrealD( const double& v0, const double& d0, int n ) : v_(v0), d_(NULL), N_(n)
{
  assert( N_ >= 0 );

  if (N_ > 0)
  {
    d_ = new double[N_];
    for (int i = 0; i < N_; i++)
      d_[i] = d0;
  }
}


ALWAYS_INLINE
SurrealD::SurrealD( const double& v0, int n ) : v_(v0), d_(NULL), N_(n)
{
  assert( N_ >= 0 );

  if (N_ > 0)
  {
    d_ = new double[N_];
    for (int i = 0; i < N_; i++)
      d_[i] = 0;
  }
}

ALWAYS_INLINE
SurrealD::SurrealD( const double& v0 ) : v_(v0), d_(NULL), N_(0) {}

ALWAYS_INLINE
SurrealD::~SurrealD()
{
  delete [] d_;
}


// assignment

ALWAYS_INLINE SurrealD&
SurrealD::operator=( const SurrealD& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new double[N_];
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it's a scalar
    (*this) = z.v_;
    return *this;
  }

  v_ = z.v_;
  for (int i = 0; i < N_; i++)
    d_[i] = z.d_[i];
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator=( const double& r )
{
  v_ = r;
  delete [] d_;
  d_ = NULL;
  N_ = 0;
//  if (N_ > 0)
//  {
//    for (int i = 0; i < N_; i++)
//      d_[i] = 0;
//  }
  return *this;
}


// Lazy assignment and unary accumulation

template< class Expr >
ALWAYS_INLINE SurrealD&
SurrealD::operator=( const SurrealDType<Expr>& r )
{
  const Expr& Tree = static_cast<const Expr&>(r);
  const int size = Tree.size();

  if (size == 0)
  {
    //Size is zero so the tree is a scalar expression
    (*this) = Tree.value();
    return *this;
  }
  else if ( N_ == 0 )
  {
    N_ = Tree.size();
    d_ = new double[N_];
  }
  else
    assert( N_ == Tree.size() );

  //This allows the compiler to do more optimizations as it knows now that the loop
  //count cannot change during the execution of the loop
  const int tmp = N_;

  //Set the derivatives
  for (int i = 0; i < tmp; ++i)
    d_[i] = Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ = Tree.value();

  return *this;
}

template< class Expr >
ALWAYS_INLINE SurrealD&
SurrealD::operator+=( const SurrealDType<Expr>& r )
{
  const Expr& Tree = static_cast<const Expr&>(r);
  const int size = Tree.size();

  if ( size == 0 )
  {
    //Size is zero so the tree is a scalar expression
    (*this) += Tree.value();
    return *this;
  }
  else if ( N_ == 0 )
  {
    N_ = Tree.size();
    d_ = new double[N_](); //Initialize to zero
  }
  else
    assert( N_ == Tree.size() );

  //This allows the compiler to do more optimizations as it knows now that the loop
  //count cannot change during the execution of the loop
  const int tmp = N_;

  //Add to the derivatives
  for (int i = 0; i < tmp; ++i)
    d_[i] += Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ += Tree.value();

  return *this;
}

template< class Expr >
ALWAYS_INLINE SurrealD&
SurrealD::operator-=( const SurrealDType<Expr>& r )
{
  const Expr& Tree = static_cast<const Expr&>(r);
  const int size = Tree.size();

  if ( size == 0 )
  {
    //Size is zero so the tree is a scalar expression
    (*this) -= Tree.value();
    return *this;
  }
  else if ( N_ == 0 )
  {
    N_ = Tree.size();
    d_ = new double[N_](); //Initialize to zero
  }
  else
    assert( N_ == Tree.size() );

  //This allows the compiler to do more optimizations as it knows now that the loop
  //count cannot change during the execution of the loop
  const int tmp = N_;

  //Add to the derivatives
  for (int i = 0; i < tmp; ++i)
    d_[i] -= Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ -= Tree.value();

  return *this;
}

// binary accumulation operators

ALWAYS_INLINE SurrealD&
SurrealD::operator+=( const SurrealD& z )
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new double[N_]();    // NOTE: d_ is value-initialized here
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) += z.v_;
    return *this;
  }

  v_ += z.v_;
  for (int i = 0; i < N_; i++)
    d_[i] += z.d_[i];

  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator+=( const double& r )
{
  v_ += r;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator+=( const int& r )
{
  v_ += r;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator-=( const SurrealD& z )
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new double[N_]();    // NOTE: d_ is value-initialized here
  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) -= z.v_;
    return *this;
  }

  v_ -= z.v_;
  for (int i = 0; i < N_; i++)
    d_[i] -= z.d_[i];

  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator-=( const double& r )
{
  v_ -= r;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator-=( const int& r )
{
  v_ -= r;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator*=( const SurrealD& z )
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new double[N_]();    // NOTE: d_ is value-initialized here

  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) *= z.v_;
    return *this;
  }

  for (int i = 0; i < N_; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;

  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator*=( const double& r )
{
  for (int i = 0; i < N_; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator*=( const int& r )
{
  for (int i = 0; i < N_; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator/=( const SurrealD& z)
{
  if ((N_ == 0) && (z.N_ > 0))
  {
    N_ = z.N_;
    d_ = new double[N_]();    // NOTE: d_ is value-initialized here

  }
  else if ( z.N_ == 0 )
  {
    //z has no derivatives, so it is a scalar
    (*this) /= z.v_;
    return *this;
  }


  double tmp = 1./(z.v_*z.v_);
  for (int i = 0; i < N_; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;

  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator/=( const double& r )
{
  double tmp = 1./r;
  for (int i = 0; i < N_; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}

ALWAYS_INLINE SurrealD&
SurrealD::operator/=( const int& r )
{
  double tmp = 1./r;
  for (int i = 0; i < N_; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}

namespace SurrealDExp
{

// Binary Lazy expressions


// Addition and Subtraction

template<class L, class R>
class OpAdd : public SurrealDType< OpAdd<L,R> >
{
public:
  OpAdd(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE double value() const { return Ll.value() + Rr.value(); }
  ALWAYS_INLINE double deriv(const int& i) const { return Ll.deriv(i) + Rr.deriv(i); }

  ALWAYS_INLINE const OpAdd&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return std::max(Ll.size(),Rr.size()); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R>
ALWAYS_INLINE SurrealDExp::OpAdd<L,R>
operator+(const SurrealDType<L>& Ll, const SurrealDType<R>& Rr)
{
  return SurrealDExp::OpAdd<L,R>( static_cast<const L&>(Ll), static_cast<const R&>(Rr) );
}

namespace SurrealDExp
{

template<class L, class R>
class OpSub : public SurrealDType< OpSub<L,R> >
{
public:
  OpSub(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE double value() const { return Ll.value() - Rr.value(); }
  ALWAYS_INLINE double deriv(const int& i) const { return Ll.deriv(i) - Rr.deriv(i); }

  ALWAYS_INLINE const OpSub&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return std::max(Ll.size(),Rr.size()); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R>
ALWAYS_INLINE SurrealDExp::OpSub<L, R>
operator-(const SurrealDType<L>& Ll, const SurrealDType<R>& Rr)
{
  return SurrealDExp::OpSub<L, R>( static_cast<const L&>(Ll), static_cast<const R&>(Rr) );
}

//Addition and Subtraction with scalar quantities

namespace SurrealDExp
{

template<class Expr>
class OpScalar : public SurrealDType< OpScalar<Expr> >
{
public:

  OpScalar(const Expr& e, const double esgn, const double s) : e(e), esgn(esgn), s(s) {}

  ALWAYS_INLINE double value() const { return esgn*e.value() + s; }
  ALWAYS_INLINE double deriv(const int& i) const { return esgn*e.deriv(i); }

  ALWAYS_INLINE const OpScalar&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const double esgn;
  const double s;
};

}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpScalar<Expr>
operator+(const SurrealDType<Expr>& e, const double& s)
{
  return SurrealDExp::OpScalar<Expr>( static_cast<const Expr&>(e), 1, s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpScalar<Expr>
operator+(const double& s, const SurrealDType<Expr>& e)
{
  return SurrealDExp::OpScalar<Expr>( static_cast<const Expr&>(e), 1, s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpScalar<Expr>
operator-(const SurrealDType<Expr>& e, const double& s)
{
  return SurrealDExp::OpScalar<Expr>( static_cast<const Expr&>(e), 1, -s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpScalar<Expr>
operator-(const double& s, const SurrealDType<Expr>& e)
{
  return SurrealDExp::OpScalar<Expr>( static_cast<const Expr&>(e), -1, s );
}


//Multiplication with Surreals

namespace SurrealDExp
{

template<class ExprL, class ExprR>
class OpMul : public SurrealDType< OpMul<ExprL, ExprR> >
{
public:
  OpMul(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), e1_val(eL.value()), e2_val(eR.value()) {}

  ALWAYS_INLINE double value() const { return e1_val*e2_val; }
  ALWAYS_INLINE double deriv(const int& i) const { return e1_val*eR.deriv(i) + eL.deriv(i)*e2_val; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return std::max(eL.size(),eR.size()); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const double e1_val, e2_val;
};

template<class ExprL>
class OpMul<ExprL, double> : public SurrealDType< OpMul<ExprL, double> >
{
public:
  OpMul(const ExprL& e, const double& s) : e(e), s(s) {}

  ALWAYS_INLINE double value() const { return e.value()*s; }
  ALWAYS_INLINE double deriv(const int& i) const { return e.deriv(i)*s; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }

  const ExprL& e;
  const double s;
};
}

//=============================================================================
template<class ExprL, class ExprR>
ALWAYS_INLINE SurrealDExp::OpMul<ExprL, ExprR>
operator*(const SurrealDType<ExprL>& z1, const SurrealDType<ExprR>& z2)
{
  return SurrealDExp::OpMul<ExprL, ExprR>( static_cast<const ExprL&>(z1), static_cast<const ExprR&>(z2) );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpMul<Expr, double>
operator*(const SurrealDType<Expr>& e, const double s)
{
  return SurrealDExp::OpMul<Expr, double>( static_cast<const Expr&>(e), s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpMul<Expr, double>
operator/(const SurrealDType<Expr>& e, const double s)
{
  return SurrealDExp::OpMul<Expr, double>( static_cast<const Expr&>(e), double(1)/s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpMul<Expr, double>
operator*(const double s, const SurrealDType<Expr>& e)
{
  return SurrealDExp::OpMul<Expr, double>( static_cast<const Expr&>(e), s );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr>
ALWAYS_INLINE SurrealDExp::OpMul<Expr, double>
operator*(const SurrealDExp::OpMul<Expr, double>& MulScal, const double s)
{
  return SurrealDExp::OpMul<Expr, double>( MulScal.e, MulScal.s*s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpMul<Expr, double>
operator/(const SurrealDExp::OpMul<Expr, double>& MulScal, const double s)
{
  return SurrealDExp::OpMul<Expr, double>( MulScal.e, MulScal.s/s );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpMul<Expr, double>
operator*(const double s, const SurrealDExp::OpMul<Expr, double>& MulScal)
{
  return SurrealDExp::OpMul<Expr, double>( MulScal.e, MulScal.s*s );
}


//=============================================================================
//Division with Surreals

namespace SurrealDExp
{

template<class ExprL, class ExprR>
class OpDiv : public SurrealDType< OpDiv<ExprL, ExprR> >
{
public:
  OpDiv(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), e1_val(eL.value()), e2_val(eR.value())
                                        , vali(1/(e2_val*e2_val)) {}

  ALWAYS_INLINE double value() const { return e1_val/e2_val; }
  ALWAYS_INLINE double deriv(const int& i) const { return (e2_val*eL.deriv(i) - eR.deriv(i)*e1_val)*vali; }

  ALWAYS_INLINE const OpDiv&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return std::max(eL.size(),eR.size()); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const double e1_val, e2_val, vali;
};

}

template< class ExprL, class ExprR >
ALWAYS_INLINE SurrealDExp::OpDiv<ExprL, ExprR>
operator/(const SurrealDType<ExprL>& eL, const SurrealDType<ExprR>& eR)
{
  return SurrealDExp::OpDiv<ExprL, ExprR>( static_cast<const ExprL&>(eL), static_cast<const ExprR&>(eR) );
}


namespace SurrealDExp
{

template<class Expr>
class OpDivScalarNumerator : public SurrealDType< OpDivScalarNumerator<Expr> >
{
public:
  OpDivScalarNumerator(const Expr& e, const double& s) : e(e), s(s), e_val(e.value()), se_val2i(s/(e_val*e_val)) {}

  ALWAYS_INLINE double value() const { return s/e_val; }
  ALWAYS_INLINE double deriv(const int& i) const { return -se_val2i*e.deriv(i); }

  ALWAYS_INLINE const OpDivScalarNumerator&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const double s;
  const double e_val, se_val2i;
};

}

template<class Expr>
ALWAYS_INLINE SurrealDExp::OpDivScalarNumerator<Expr>
operator/(const double s, const SurrealDType<Expr>& e)
{
  return SurrealDExp::OpDivScalarNumerator<Expr>( static_cast<const Expr&>(e), s );
}

// unary operators; no side effects

ALWAYS_INLINE const SurrealD&
SurrealD::operator+() const
{
  return *this;
}

ALWAYS_INLINE const SurrealDExp::OpMul<SurrealD, double>
SurrealD::operator-() const
{
  return SurrealDExp::OpMul<SurrealD, double>(*this, -1);
}

template< class Expr >
ALWAYS_INLINE const SurrealDExp::OpMul<Expr, double>
operator-(SurrealDType<Expr> const& e)
{
  return SurrealDExp::OpMul<Expr, double>( static_cast<const Expr&>(e), -1 );
}


// relational operators

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator==( const SurrealDType<ExprL>& lhs, const SurrealDType<ExprR>& rhs )
{
  return lhs.value() == rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator==( const SurrealDType<Expr>& lhs, const double& rhs )
{
  return lhs.value() == rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator==( const double& lhs, const SurrealDType<Expr>& rhs )
{
  return lhs == rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator!=( const SurrealDType<ExprL>& lhs, const SurrealDType<ExprR>& rhs )
{
  return lhs.value() != rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator!=( const SurrealDType<Expr>& lhs, const double& rhs )
{
  return lhs.value() != rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator!=( const double& lhs, const SurrealDType<Expr>& rhs )
{
  return lhs != rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator>( const SurrealDType<ExprL>& lhs, const SurrealDType<ExprR>& rhs )
{
  return lhs.value() > rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator>( const SurrealDType<Expr>& lhs, const double& rhs )
{
  return lhs.value() > rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator>( const double& lhs, const SurrealDType<Expr>& rhs )
{
  return lhs > rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator<( const SurrealDType<ExprL>& lhs, const SurrealDType<ExprR>& rhs )
{
  return lhs.value() < rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator<( const SurrealDType<Expr>& lhs, const double& rhs )
{
  return lhs.value() < rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator<( const double& lhs, const SurrealDType<Expr>& rhs )
{
  return lhs < rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator>=( const SurrealDType<ExprL>& lhs, const SurrealDType<ExprR>& rhs )
{
  return lhs.value() >= rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator>=( const SurrealDType<Expr>& lhs, const double& rhs )
{
  return lhs.value() >= rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator>=( const double& lhs, const SurrealDType<Expr>& rhs )
{
  return lhs >= rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator<=( const SurrealDType<ExprL>& lhs, const SurrealDType<ExprR>& rhs )
{
  return lhs.value() <= rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator<=( const SurrealDType<Expr>& lhs, const double& rhs )
{
  return lhs.value() <= rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator<=( const double& lhs, const SurrealDType<Expr>& rhs )
{
  return lhs <= rhs.value();
}



//Functions for Surreals
#define SURREALD_FUNC1( NAME, FUNC, DERIV ) \
namespace SurrealDExp \
{  \
template<class Expr> \
class Surreal_ ## NAME : public SurrealDType< Surreal_ ## NAME<Expr> > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  Surreal_ ## NAME(const Expr& e) : e(e), z(e.value()), der(DERIV) {} \
  \
  ALWAYS_INLINE double value() const { return FUNC; } \
  ALWAYS_INLINE double deriv(const int& i) const { return der*e.deriv(i); } \
  \
  ALWAYS_INLINE const Surreal_ ## NAME& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return e.size(); } \
private: \
  const Expr& e; \
  const double z, der; \
}; \
} \
\
template<class Expr> \
ALWAYS_INLINE SurrealDExp::Surreal_ ## NAME<Expr> \
NAME(const SurrealDType<Expr>& z) { return SurrealDExp::Surreal_ ## NAME<Expr>( static_cast<const Expr&>(z) ); }


#define SURREALD_FUNC2( NAME, FUNC, DERIV ) \
namespace SurrealDExp \
{  \
template<class ExprL, class ExprR> \
class Surreal_ ## NAME : public SurrealDType< Surreal_ ## NAME<ExprL, ExprR> > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  Surreal_ ## NAME(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), z1(eL.value()), z2(eR.value()), \
                                                                 der(DERIV) {} \
  \
  ALWAYS_INLINE double value() const { return FUNC; } \
  ALWAYS_INLINE double deriv(const int& i) const { return der*(z2*eL.deriv(i) - z1*eR.deriv(i)); } \
  \
  ALWAYS_INLINE const Surreal_ ## NAME& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return std::max(eL.size(),eR.size()); } \
private: \
  const ExprL& eL; \
  const ExprR& eR; \
  const double z1, z2, der; \
}; \
  \
} \
\
template<class ExprL, class ExprR> \
ALWAYS_INLINE SurrealDExp::Surreal_ ## NAME<ExprL, ExprR> \
NAME(const SurrealDType<ExprL>& z1, const SurrealDType<ExprR>& z2) \
{ return SurrealDExp::Surreal_ ## NAME<ExprL, ExprR>( static_cast<const ExprL&>(z1), static_cast<const ExprR&>(z2) ); } \
\
ALWAYS_INLINE SurrealDExp::Surreal_ ## NAME<SurrealD, SurrealD> \
NAME(const SurrealD& z1, const SurrealD& z2) \
{ return SurrealDExp::Surreal_ ## NAME<SurrealD, SurrealD>( z1, z2 ); }

// trig functions <cmath>

SURREALD_FUNC1( cos, std::cos(z), -std::sin(z) )
SURREALD_FUNC1( sin, std::sin(z),  std::cos(z) )
SURREALD_FUNC1( tan, std::tan(z),  double(1)/(std::cos(z)*std::cos(z)) )
SURREALD_FUNC1( acos, std::acos(z), -double(1)/std::sqrt(1. - z*z) )
SURREALD_FUNC1( asin, std::asin(z),  double(1)/std::sqrt(1. - z*z) )
SURREALD_FUNC1( atan, std::atan(z),  double(1)/(1. + z*z) )

SURREALD_FUNC2( atan2, std::atan2(z1, z2),  double(1)/(z1*z1 + z2*z2) )

// hyperbolic functions <cmath>

SURREALD_FUNC1( cosh, std::cosh(z), std::sinh(z) )
SURREALD_FUNC1( sinh, std::sinh(z), std::cosh(z) )
SURREALD_FUNC1( tanh, std::tanh(z), double(1)/(std::cosh(z)*std::cosh(z)) )

// exp and log functions <cmath>

SURREALD_FUNC1( exp, std::exp(z), std::exp(z) )
SURREALD_FUNC1( log, std::log(z), double(1)/z )
SURREALD_FUNC1( log10, std::log10(z), double(1)/(z*std::log(10.)) )

// power functions <cmath>

namespace SurrealDExp
{

template<class ExprL, class ExprR>
class Surreal_pow : public SurrealDType< Surreal_pow<ExprL, ExprR> >
{ /*This is for functions when the argument is an expression*/
public:
  Surreal_pow(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), a(eL.value()), b(eR.value()),
                                                powab(std::pow(a,b)),
                                                tmp1( b*std::pow(a, b - 1) ),
                                                tmp2( powab*std::log(a) ) {}

  ALWAYS_INLINE double value() const { return powab; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp1*eL.deriv(i) + tmp2*eR.deriv(i); }

  ALWAYS_INLINE const Surreal_pow&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return std::max(eL.size(),eR.size()); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const double a, b, powab, tmp1, tmp2;
};

template<class ExprL>
class Surreal_pow<ExprL, double> : public SurrealDType< Surreal_pow<ExprL, double> >
{ /*This is optimized when the argument is SurrealView and double*/
public:
  Surreal_pow(const ExprL& eL, const double& b) : eL(eL), a(eL.value()),
                                               powab(std::pow(a,b)),
                                               tmp1( b*std::pow(a, b - 1) ) {}

  ALWAYS_INLINE double value() const { return powab; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp1*eL.deriv(i); }

  ALWAYS_INLINE const Surreal_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const double a, powab, tmp1;
};


template<class ExprR>
class Surreal_pow<double, ExprR> : public SurrealDType< Surreal_pow<double, ExprR> >
{ /*This is optimized when the argument is a double and SurrealView*/
public:
  Surreal_pow(const double& a, const ExprR& eR) : eR(eR), b(eR.value()),
                                               powab(std::pow(a,b)),
                                               tmp2( powab*std::log(a) ) {}

  ALWAYS_INLINE double value() const { return powab; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp2*eR.deriv(i); }

  ALWAYS_INLINE const Surreal_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eR.size(); }
private:
  const ExprR& eR;
  const double b, powab, tmp2;

};

}

template<class ExprL, class ExprR>
ALWAYS_INLINE SurrealDExp::Surreal_pow<ExprL, ExprR>
pow(const SurrealDType<ExprL>& a, const SurrealDType<ExprR>& b)
{
  return SurrealDExp::Surreal_pow<ExprL, ExprR>( static_cast<const ExprL&>(a), static_cast<const ExprR&>(b) );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::Surreal_pow<Expr, double>
pow(const SurrealDType<Expr>& a, const double b )
{
  return SurrealDExp::Surreal_pow<Expr, double>( static_cast<const Expr&>(a), b );
}

template<class Expr>
ALWAYS_INLINE SurrealDExp::Surreal_pow<double, Expr>
pow(const double a, const SurrealDType<Expr>& b)
{
  return SurrealDExp::Surreal_pow<double, Expr>( a, static_cast<const Expr&>(b) );
}


namespace SurrealDExp
{

template<class Expr>
class Surreal_sqrt : public SurrealDType< Surreal_sqrt<Expr> >
{ /*This is optimized when the argument is an Expression*/
public:
  ALWAYS_INLINE
  Surreal_sqrt(const Expr& e) : e(e), sqrtv( sqrt(e.value()) ), tmp( sqrtv == 0 ? 0 : 0.5/sqrtv ) {}

  ALWAYS_INLINE double value() const { return sqrtv; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp*e.deriv(i); }

  ALWAYS_INLINE const Surreal_sqrt
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const double sqrtv, tmp;
};

}

template<class Expr>
ALWAYS_INLINE SurrealDExp::Surreal_sqrt<Expr>
sqrt(const SurrealDType<Expr>& z)
{
  return SurrealDExp::Surreal_sqrt<Expr>( static_cast<const Expr&>(z) );
}


// rounding functions <cmath>

SURREALD_FUNC1( ceil, std::ceil(z), 0 )
SURREALD_FUNC1( floor, std::floor(z), 0 )

// misc functions <cmath>

ALWAYS_INLINE SurrealDExp::OpMul<SurrealD, double>
abs( const SurrealD& z )
{
  return (z.v_ < 0) ?
         SurrealDExp::OpMul<SurrealD, double>( z, -1 ) :
         SurrealDExp::OpMul<SurrealD, double>( z,  1 );
}

ALWAYS_INLINE SurrealDExp::OpMul<SurrealD, double>
fabs( const SurrealD& z )
{
  return (z.v_ < 0) ?
         SurrealDExp::OpMul<SurrealD, double>( z, -1 ) :
         SurrealDExp::OpMul<SurrealD, double>( z,  1 );
}

//Clean up macro definitions
#undef SURREALD_FUNC1
#undef SURREALD_FUNC2


template<class Expr>
std::ostream&
operator<<( std::ostream& os, const SurrealDType<Expr>& ztype )
{
  const Expr& z = static_cast<const Expr&>(ztype);
  const int N = z.size();
  os << '(' << z.value();
  if (N > 0)
  {
    os << ';';
    for (int i = 0; i < N - 1; i++)
      os << z.deriv(i) << ',';
    os << z.deriv(N - 1);
  }
  os << ')';
  return os;
}



#endif // SURREALD_LAZY_H
