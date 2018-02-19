#ifndef SURREALD_TRAD_H
#define SURREALD_TRAD_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <cassert>

#include "always_inline.h"

//----------------------------------------------------------------------------//
// SurrealD:  value, N derivatives
//
// Operators with traditional implementation
//
//----------------------------------------------------------------------------//

class SurrealD
{
public:
  //The default constructor is intentionally not included here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  SurrealD( const SurrealD& z );
  explicit SurrealD( int n=0 );
  SurrealD( const double& v0, const double d0[], const int n );
  SurrealD( const double& v0, const double& d0, const int n );
  SurrealD( const int& v0, const int& d0, const int n );
  SurrealD( const double& v0 );
  SurrealD( const double& v0, const int n );
  ~SurrealD();

  int size() const { return N_; }

  // value accessor operators
  ALWAYS_INLINE       double& value()       { return v_; }
  ALWAYS_INLINE const double& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE double& deriv( int i=0 )       { assert(N_ > 0); return d_[i]; }
  ALWAYS_INLINE double  deriv( int i=0 ) const { return N_ > 0 ? d_[i] : 0; }

  // assignment
  SurrealD& operator=( const SurrealD& );
  SurrealD& operator=( const double& );

  // unary operators; no side effects
  const SurrealD& operator+() const;
  const SurrealD  operator-() const;

  // binary accumulation operators
  SurrealD& operator+=( const SurrealD& );
  SurrealD& operator+=( const double& );
  SurrealD& operator-=( const SurrealD& );
  SurrealD& operator-=( const double& );
  SurrealD& operator*=( const SurrealD& );
  SurrealD& operator*=( const double& );
  SurrealD& operator/=( const SurrealD& );
  SurrealD& operator/=( const double& );

  // binary operators
  friend SurrealD operator+( const SurrealD&, const SurrealD& );
  friend SurrealD operator+( const SurrealD&, const double& );
  friend SurrealD operator+( const double&, const SurrealD& );
  friend SurrealD operator-( const SurrealD&, const SurrealD& );
  friend SurrealD operator-( const SurrealD&, const double& );
  friend SurrealD operator-( const double&, const SurrealD& );
  friend SurrealD operator*( const SurrealD&, const SurrealD& );
  friend SurrealD operator*( const SurrealD&, const double& );
  friend SurrealD operator*( const double&, const SurrealD& );
  friend SurrealD operator/( const SurrealD&, const SurrealD& );
  friend SurrealD operator/( const SurrealD&, const double& );
  friend SurrealD operator/( const double&, const SurrealD& );

  // relational operators
  friend bool operator==( const SurrealD&, const SurrealD& );
  friend bool operator==( const SurrealD&, const double& );
  friend bool operator==( const double&, const SurrealD& );
  friend bool operator!=( const SurrealD&, const SurrealD& );
  friend bool operator!=( const SurrealD&, const double& );
  friend bool operator!=( const double&, const SurrealD& );
  friend bool operator>( const SurrealD&, const SurrealD& );
  friend bool operator>( const SurrealD&, const double& );
  friend bool operator>( const double&, const SurrealD& );
  friend bool operator<( const SurrealD&, const SurrealD& );
  friend bool operator<( const SurrealD&, const double& );
  friend bool operator<( const double&, const SurrealD& );
  friend bool operator>=( const SurrealD&, const SurrealD& );
  friend bool operator>=( const SurrealD&, const double& );
  friend bool operator>=( const double&, const SurrealD& );
  friend bool operator<=( const SurrealD&, const SurrealD& );
  friend bool operator<=( const SurrealD&, const double& );
  friend bool operator<=( const double&, const SurrealD& );

  // trig functions <cmath>
  friend SurrealD cos( const SurrealD& );
  friend SurrealD sin( const SurrealD& );
  friend SurrealD tan( const SurrealD& );
  friend SurrealD acos( const SurrealD& );
  friend SurrealD asin( const SurrealD& );
  friend SurrealD atan( const SurrealD& );
  friend SurrealD atan2( const SurrealD&, const SurrealD& );

  // hyperbolic functions <cmath>
  friend SurrealD cosh( const SurrealD& );
  friend SurrealD sinh( const SurrealD& );
  friend SurrealD tanh( const SurrealD& );

  // exp and log functions <cmath>
  friend SurrealD exp( const SurrealD& );
  friend SurrealD log( const SurrealD& );
  friend SurrealD log10( const SurrealD& );

  // power functions <cmath>
  friend SurrealD pow( const SurrealD&, const SurrealD& );
  friend SurrealD pow( const SurrealD&, const double& );
  friend SurrealD pow( const double&, const SurrealD& );
  friend SurrealD sqrt( const SurrealD& );

  // rounding functions <cmath>
  friend SurrealD ceil( const SurrealD& );
  friend SurrealD floor( const SurrealD& );

  // misc functions <cmath>
  friend SurrealD abs( const SurrealD& );
  friend SurrealD fabs( const SurrealD& );

  // classification functions <cmath>
  friend bool isfinite( const SurrealD& );
  friend bool isinf( const SurrealD& );
  friend bool isnan( const SurrealD& );

  // input/output
  friend std::istream& operator>>( std::istream&, SurrealD& );
  friend std::ostream& operator<<( std::ostream&, const SurrealD& );

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
SurrealD::SurrealD( int n ) : v_(n), d_(NULL), N_(0) {}

ALWAYS_INLINE
SurrealD::SurrealD( const double& v0, const double d0[], const int n ) : v_(v0), d_(NULL), N_(n)
{
  assert( N_ > 0 );

  d_ = new double[N_];
  for (int i = 0; i < N_; i++)
    d_[i] = d0[i];
}

ALWAYS_INLINE
SurrealD::SurrealD( const double& v0, const double& d0, const int n ) : v_(v0), d_(NULL), N_(n)
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
SurrealD::SurrealD( const int& v0, const int& d0, const int n ) : v_(v0), d_(NULL), N_(n)
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
SurrealD::SurrealD( const double& v0 ) : v_(v0), d_(NULL), N_(0) {}

ALWAYS_INLINE
SurrealD::SurrealD( const double& v0, const int n ) : v_(v0), d_(NULL), N_(n)
{
  assert( N_ > 0 );

  d_ = new double[N_];
  for ( int i = 0; i < N_; i++ ) d_[i] = 0;
}

ALWAYS_INLINE
SurrealD::~SurrealD()
{
  delete [] d_;
}


// assignment

ALWAYS_INLINE  SurrealD&
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

ALWAYS_INLINE  SurrealD&
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


// unary operators; no side effects

ALWAYS_INLINE  const SurrealD&
SurrealD::operator+() const
{
  return *this;
}

ALWAYS_INLINE  const SurrealD
SurrealD::operator-() const
{
  if (N_ == 0)
  {
    return SurrealD( -v_ );
  }
  else
  {
    SurrealD c(-v_, N_);
    for (int i = 0; i < N_; i++)
      c.d_[i] = -d_[i];
    return c;
  }
}


// binary accumulation operators

ALWAYS_INLINE  SurrealD&
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

ALWAYS_INLINE  SurrealD&
SurrealD::operator+=( const double& r )
{
  v_ += r;
  return *this;
}

ALWAYS_INLINE  SurrealD&
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

ALWAYS_INLINE  SurrealD&
SurrealD::operator-=( const double& r )
{
  v_ -= r;
  return *this;
}

ALWAYS_INLINE  SurrealD&
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

ALWAYS_INLINE  SurrealD&
SurrealD::operator*=( const double& r )
{
  for (int i = 0; i < N_; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

ALWAYS_INLINE  SurrealD&
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

ALWAYS_INLINE  SurrealD&
SurrealD::operator/=( const double& r )
{
  double tmp = 1./r;
  for (int i = 0; i < N_; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}


// binary operators

ALWAYS_INLINE  SurrealD
operator+( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ + b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
    return SurrealD(a.v_ + b.v_, a.d_, a.N_);
  else if (a.N_ == 0 && b.N_ > 0)
    return SurrealD(a.v_ + b.v_, b.d_, b.N_);

  assert( a.N_ == b.N_ );

  const int N = b.N_;
  SurrealD c(a.v_ + b.v_, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] + b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator+( const SurrealD& a, const double& b )
{
  if (a.N_ == 0 )
    return SurrealD(a.v_ + b);

  return SurrealD(a.v_ + b, a.d_, a.N_);
}

ALWAYS_INLINE  SurrealD
operator+( const double& a, const SurrealD& b )
{
  if (b.N_ == 0 )
    return SurrealD(a + b.v_);

  return SurrealD(a + b.v_, b.d_, b.N_);
}

ALWAYS_INLINE  SurrealD
operator-( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ - b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
    return SurrealD(a.v_ - b.v_, a.d_, a.N_);
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const int N = b.N_;
    SurrealD c(a.v_ - b.v_, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = -b.d_[i];
    return c;
  }

  assert( a.N_ == b.N_ );

  const int N = b.N_;
  SurrealD c(a.v_ - b.v_, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] - b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator-( const SurrealD& a, const double& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ - b);

  return SurrealD(a.v_ - b, a.d_, a.N_);
}

ALWAYS_INLINE  SurrealD
operator-( const double& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a - b.v_);

  const int N = b.N_;
  SurrealD c(a - b.v_, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = -b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator*( const SurrealD& a, const SurrealD& b )
{
  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ * b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    const int N = a.N_;
    SurrealD c(a.v_ * b.v_, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = a.d_[i]*b.v_;
    return c;
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const int N = b.N_;
    SurrealD c(a.v_ * b.v_, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = a.v_*b.d_[i];
    return c;
  }

  const int N = b.N_;
  SurrealD c(a.v_ * b.v_, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*b.v_ + a.v_*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator*( const SurrealD& a, const double& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ * b);

  const int N = a.N_;
  SurrealD c(a.v_ * b, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*b;
  return c;
}

ALWAYS_INLINE  SurrealD
operator*( const double& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a * b.v_);

  const int N = b.N_;
  SurrealD c(a * b.v_, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = a*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
operator/( const SurrealD& a, const SurrealD& b )
{
  double tmp = 1./(b.v_*b.v_);

  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(a.v_ / b.v_);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    const int N = a.N_;
    SurrealD c(a.v_ / b.v_, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = (b.v_*a.d_[i])*tmp;
    return c;
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    const int N = b.N_;
    SurrealD c(a.v_ / b.v_, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = (-a.v_*b.d_[i])*tmp;
    return c;
  }

  assert( a.N_ == b.N_ );

  const int N = b.N_;
  SurrealD c(a.v_ / b.v_, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = (b.v_*a.d_[i] - a.v_*b.d_[i])*tmp;
  return c;
}

ALWAYS_INLINE  SurrealD operator/( const SurrealD& a, const double& b )
{
  if (a.N_ == 0)
    return SurrealD(a.v_ / b);

  double tmp = 1./(b);
  const int N = a.N_;
  SurrealD c(a.v_ / b, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*tmp;
  return c;
}

ALWAYS_INLINE  SurrealD
operator/( const double& a, const SurrealD& b )
{
  if (b.N_ == 0)
    return SurrealD(a / b.v_);

  const int N = b.N_;
  double tmpv = a/(b.v_);
  double tmpd = -1./(b.v_);
  SurrealD c(tmpv, N);
  for (int i = 0; i < N; i++)
    c.d_[i] = tmpv*tmpd*b.d_[i];
  return c;
}


// relational operators

ALWAYS_INLINE  bool
operator==( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ == rhs.v_;
}

ALWAYS_INLINE  bool
operator==( const SurrealD& lhs, const double& rhs )
{
  return lhs.v_ == rhs;
}

ALWAYS_INLINE  bool
operator==( const double& lhs, const SurrealD& rhs )
{
  return lhs == rhs.v_;
}

ALWAYS_INLINE  bool
operator!=( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ != rhs.v_;
}

ALWAYS_INLINE  bool
operator!=( const SurrealD& lhs, const double& rhs )
{
  return lhs.v_ != rhs;
}

ALWAYS_INLINE  bool
operator!=( const double& lhs, const SurrealD& rhs )
{
  return lhs != rhs.v_;
}

ALWAYS_INLINE  bool
operator>( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ > rhs.v_;
}

ALWAYS_INLINE  bool
operator>( const SurrealD& lhs, const double& rhs )
{
  return lhs.v_ > rhs;
}

ALWAYS_INLINE  bool
operator>( const double& lhs, const SurrealD& rhs )
{
  return lhs > rhs.v_;
}

ALWAYS_INLINE  bool
operator<( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ < rhs.v_;
}

ALWAYS_INLINE  bool
operator<( const SurrealD& lhs, const double& rhs )
{
  return lhs.v_ < rhs;
}

ALWAYS_INLINE  bool
operator<( const double& lhs, const SurrealD& rhs )
{
  return lhs < rhs.v_;
}

ALWAYS_INLINE  bool
operator>=( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ >= rhs.v_;
}

ALWAYS_INLINE  bool
operator>=( const SurrealD& lhs, const double& rhs )
{
  return lhs.v_ >= rhs;
}

ALWAYS_INLINE  bool
operator>=( const double& lhs, const SurrealD& rhs )
{
  return lhs >= rhs.v_;
}

ALWAYS_INLINE  bool
operator<=( const SurrealD& lhs, const SurrealD& rhs )
{
  return lhs.v_ <= rhs.v_;
}

ALWAYS_INLINE  bool
operator<=( const SurrealD& lhs, const double& rhs )
{
  return lhs.v_ <= rhs;
}

ALWAYS_INLINE  bool
operator<=( const double& lhs, const SurrealD& rhs )
{
  return lhs <= rhs.v_;
}


//Macros for functions

#define SURREALD_FUNC1( NAME, FUNC, DERIV ) \
ALWAYS_INLINE  SurrealD \
NAME( const SurrealD& z ) \
{ \
  if ( z.N_ == 0 ) \
    return SurrealD(FUNC); \
  \
  const int N = z.N_; \
  double tmp = DERIV; \
  SurrealD c(FUNC, N); \
  for (int i = 0; i < N; i++) \
    c.d_[i] = tmp*z.d_[i]; \
  return c; \
}

#define SURREALD_FUNC2( NAME, FUNC, DERIV ) \
ALWAYS_INLINE  SurrealD \
NAME( const SurrealD& z1, const SurrealD& z2) \
{ \
  if ( z1.N_ == 0 && z2.N_ == 0 ) \
    return SurrealD(FUNC); \
  else if ( z1.N_ > 0 && z2.N_ == 0 ) \
  { \
    const int N = z1.N_; \
    double tmp = DERIV; \
    SurrealD c(FUNC, N); \
    for (int i = 0; i < N; i++) \
      c.d_[i] = tmp*(z2.v_*z1.d_[i]); \
    return c; \
  } \
  else if ( z1.N_ == 0 && z2.N_ > 0 ) \
  { \
    const int N = z2.N_; \
    double tmp = DERIV; \
    SurrealD c(FUNC, N); \
    for (int i = 0; i < N; i++) \
      c.d_[i] = tmp*(-z1.v_*z2.d_[i]); \
    return c; \
  } \
  \
  assert( z1.N_ == z2.N_ ); \
  \
  const int N = z1.N_; \
  double tmp = DERIV; \
  SurrealD c(FUNC, N); \
  for (int i = 0; i < N; i++) \
    c.d_[i] = tmp*(z2.v_*z1.d_[i] - z1.v_*z2.d_[i]); \
  return c; \
}

// trig functions <cmath>

SURREALD_FUNC1( cos, std::cos(z.v_), -std::sin(z.v_) )
SURREALD_FUNC1( sin, std::sin(z.v_),  std::cos(z.v_) )
SURREALD_FUNC1( tan, std::tan(z.v_),  double(1)/(std::cos(z.v_)*std::cos(z.v_)) )
SURREALD_FUNC1( acos, std::acos(z.v_), -double(1)/std::sqrt(1 - z.v_*z.v_) )
SURREALD_FUNC1( asin, std::asin(z.v_),  double(1)/std::sqrt(1 - z.v_*z.v_) )
SURREALD_FUNC1( atan, std::atan(z.v_),  double(1)/(1 + z.v_*z.v_) )

SURREALD_FUNC2( atan2, std::atan2(z1.v_, z2.v_),  double(1)/(z1.v_*z1.v_ + z2.v_*z2.v_) )

// hyperbolic functions <cmath>

SURREALD_FUNC1( cosh, std::cosh(z.v_), std::sinh(z.v_) )
SURREALD_FUNC1( sinh, std::sinh(z.v_), std::cosh(z.v_) )
SURREALD_FUNC1( tanh, std::tanh(z.v_), double(1)/(std::cosh(z.v_)*std::cosh(z.v_)) )

// exp and log functions <cmath>

SURREALD_FUNC1( exp, std::exp(z.v_), std::exp(z.v_) )
SURREALD_FUNC1( log, std::log(z.v_), double(1)/z.v_ )
SURREALD_FUNC1( log10, std::log10(z.v_), double(1)/(z.v_*std::log(10.)) )

// power functions <cmath>

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const SurrealD& b)
{
  double powab=std::pow(a.v_,b.v_);

  if (a.N_ == 0 && b.N_ == 0)
    return SurrealD(powab);
  else if (a.N_ > 0 && b.N_ == 0)
  {
    double tmp1 = b.v_*std::pow(a.v_, b.v_ - 1);

    const int N = a.N_;
    SurrealD c(powab, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = tmp1*a.d_[i];
    return c;
  }
  else if (a.N_ == 0 && b.N_ > 0)
  {
    double tmp2 = powab*std::log(a.v_);

    const int N = b.N_;
    SurrealD c(powab, N);
    for (int i = 0; i < N; i++)
      c.d_[i] = tmp2*b.d_[i];
    return c;
  }

  assert(a.N_ == b.N_);
  // many sticky points were derivative is undefined or infinite
  // badness if 0 <= b < 1 and a == 0
  double tmp1 = b.v_*std::pow(a.v_, b.v_ - 1);
  double tmp2 = powab*std::log(a.v_);

  SurrealD c(powab, a.N_);
  for (int i = 0; i < a.N_; i++)
    c.d_[i] = tmp1*a.d_[i] + tmp2*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
pow( const SurrealD& a, const double& b)
{
  double powab=std::pow(a.v_,b);

  if (a.N_ == 0)
    return SurrealD(powab);

  double tmp = b*std::pow(a.v_, b-1);
  SurrealD c(powab, a.N_);
  for (int i = 0; i < a.N_; i++)
    c.d_[i] = tmp*a.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
pow( const double& a, const SurrealD& b)
{
  double powab=std::pow(a, b.v_);

  if (b.N_ == 0)
    return SurrealD(powab);

  double tmp = powab*std::log(a);
  SurrealD c(powab, b.N_);
  for (int i = 0; i < b.N_; i++)
    c.d_[i] = tmp*b.d_[i];
  return c;
}

ALWAYS_INLINE  SurrealD
sqrt( const SurrealD& z )
{
  double sqrtv=std::sqrt(z.v_);
  if (sqrtv == 0)
  {
    return SurrealD(0., 0., z.N_);
  }
  else
  {
    double tmp = 0.5/sqrtv;
    SurrealD c(sqrtv, z.N_);
    for (int i = 0; i < z.N_; i++)
      c.d_[i] = tmp*z.d_[i];
    return c;
  }
}


// rounding functions <cmath>

ALWAYS_INLINE  SurrealD
ceil( const SurrealD& z )
{
  if (z.N_ == 0)
    return SurrealD(std::ceil(z.v_));

  return SurrealD(std::ceil(z.v_), 0., z.N_);
}

ALWAYS_INLINE  SurrealD
floor( const SurrealD& z )
{
  if (z.N_ == 0)
    return SurrealD(std::floor(z.v_));

  return SurrealD(std::floor(z.v_), 0., z.N_);
}


// misc functions <cmath>

ALWAYS_INLINE  SurrealD
abs( const SurrealD& z )
{
  return (z.v_ < 0) ? -z : z;
}

ALWAYS_INLINE  SurrealD
fabs( const SurrealD& z )
{
  return (z.v_ < 0) ? -z : z;
}

//Clean up macro definitions
#undef SURREALD_FUNC1
#undef SURREALD_FUNC2


// output format: (v;d0,d1,d2,...,dN)

inline std::ostream&
operator<<( std::ostream& os, const SurrealD& z )
{
  os << '(' << z.v_;
  if (z.N_ > 0)
  {
    os << ';';
    for (int i = 0; i < z.N_ - 1; i++)
      os << z.d_[i] << ',';
    os << z.d_[z.N_ - 1];
  }
  os << ')';
  return os;
}


// debug dump of private data
inline void
SurrealD::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "SurrealD: v_ = " << v_;
  if (N_ == 0)
  {
    if ( d_ == NULL )
      out << "  d_ = " << 0;
    else
      out << "  d_ = " << d_;
    out << "  N_ = " << N_ << std::endl;
  }
  else
  {
    out << "  d_[" << N_ << "] = (";
    for (int n = 0; n < N_-1; n++)
      out << d_[n] << ",";
    out << d_[N_-1] << ")" << std::endl;
  }
}

#endif // SURREALD_TRAD_H
