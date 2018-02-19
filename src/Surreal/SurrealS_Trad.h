#ifndef SURREALS_TRAD_H
#define SURREALS_TRAD_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <cassert>

#include <string>

#include "always_inline.h"

//----------------------------------------------------------------------------//
// SurrealS<N>:  value, N derivatives
//
// Traditional implementation of Operators
//
// statically defined derivative array
//----------------------------------------------------------------------------//

template<int N>
class SurrealS
{
public:
  //The default constructor is intentionally empty here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  ALWAYS_INLINE SurrealS() {}
  SurrealS( const SurrealS& z );
  SurrealS( const double& v0, const double d0[], const int n );
  SurrealS( const double& v0, const double& d0 );
  SurrealS( const double& v0 );

  ALWAYS_INLINE ~SurrealS() {}

#if 0
  operator double() const { return v_; }
  operator int() const { return int(v_); }
#endif

  ALWAYS_INLINE int size() const { return N; }

  // value accessor operators
  ALWAYS_INLINE double& value()       { return v_; }
  ALWAYS_INLINE double  value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE double& deriv( int i=0 )       { return d_[i]; }
  ALWAYS_INLINE double  deriv( int i=0 ) const { return d_[i]; }

  // assignment
  SurrealS& operator=( const SurrealS& );
  SurrealS& operator=( const double& );

  // unary operators; no side effects
  const SurrealS& operator+() const;
  const SurrealS  operator-() const;

  // binary accumulation operators
  SurrealS& operator+=( const SurrealS& );
  SurrealS& operator+=( const double& );
  SurrealS& operator-=( const SurrealS& );
  SurrealS& operator-=( const double& );
  SurrealS& operator*=( const SurrealS& );
  SurrealS& operator*=( const double& );
  SurrealS& operator/=( const SurrealS& );
  SurrealS& operator/=( const double& );

  // binary operators
  template<int M> friend SurrealS<M> operator+( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator+( const SurrealS<M>&, const double& );
  template<int M> friend SurrealS<M> operator+( const double&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator-( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator-( const SurrealS<M>&, const double& );
  template<int M> friend SurrealS<M> operator-( const double&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator*( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator*( const SurrealS<M>&, const double& );
  template<int M> friend SurrealS<M> operator*( const double&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator/( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> operator/( const SurrealS<M>&, const double& );
  template<int M> friend SurrealS<M> operator/( const double&, const SurrealS<M>& );

  // relational operators
  template<int M> friend bool operator==( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator==( const SurrealS<M>&, const double& );
  template<int M> friend bool operator==( const double&, const SurrealS<M>& );
  template<int M> friend bool operator!=( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator!=( const SurrealS<M>&, const double& );
  template<int M> friend bool operator!=( const double&, const SurrealS<M>& );
  template<int M> friend bool operator>( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator>( const SurrealS<M>&, const double& );
  template<int M> friend bool operator>( const double&, const SurrealS<M>& );
  template<int M> friend bool operator<( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator<( const SurrealS<M>&, const double& );
  template<int M> friend bool operator<( const double&, const SurrealS<M>& );
  template<int M> friend bool operator>=( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator>=( const SurrealS<M>&, const double& );
  template<int M> friend bool operator>=( const double&, const SurrealS<M>& );
  template<int M> friend bool operator<=( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend bool operator<=( const SurrealS<M>&, const double& );
  template<int M> friend bool operator<=( const double&, const SurrealS<M>& );

  // trig functions <cmath>
  template<int M> friend SurrealS<M> cos( const SurrealS<M>& );
  template<int M> friend SurrealS<M> sin( const SurrealS<M>& );
  template<int M> friend SurrealS<M> tan( const SurrealS<M>& );
  template<int M> friend SurrealS<M> acos( const SurrealS<M>& );
  template<int M> friend SurrealS<M> asin( const SurrealS<M>& );
  template<int M> friend SurrealS<M> atan( const SurrealS<M>& );
  template<int M> friend SurrealS<M> atan2( const SurrealS<M>&, const SurrealS<M>& );

  // hyperbolic functions <cmath>
  template<int M> friend SurrealS<M> cosh( const SurrealS<M>& );
  template<int M> friend SurrealS<M> sinh( const SurrealS<M>& );
  template<int M> friend SurrealS<M> tanh( const SurrealS<M>& );

  // exp and log functions <cmath>
  template<int M> friend SurrealS<M> exp( const SurrealS<M>& );
  template<int M> friend SurrealS<M> log( const SurrealS<M>& );
  template<int M> friend SurrealS<M> log10( const SurrealS<M>& );

  // power functions <cmath>
  template<int M> friend SurrealS<M> pow( const SurrealS<M>&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> pow( const SurrealS<M>&, const double& );
  template<int M> friend SurrealS<M> pow( const SurrealS<M>&, const int& );
  template<int M> friend SurrealS<M> pow( const double&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> pow( const int&, const SurrealS<M>& );
  template<int M> friend SurrealS<M> sqrt( const SurrealS<M>& );

  // rounding functions <cmath>
  template<int M> friend SurrealS<M> ceil( const SurrealS<M>& );
  template<int M> friend SurrealS<M> floor( const SurrealS<M>& );

  // misc functions <cmath>
  template<int M> friend SurrealS<M> abs( const SurrealS<M>& );
  template<int M> friend SurrealS<M> fabs( const SurrealS<M>& );

  // classification functions <cmath>
  template<int M> friend bool isfinite( const SurrealS<M>& );
  template<int M> friend bool isinf( const SurrealS<M>& );
  template<int M> friend bool isnan( const SurrealS<M>& );

  // input/output
  template<int M> friend std::istream& operator>>( std::istream&, SurrealS<M>& );
  template<int M> friend std::ostream& operator<<( std::ostream&, const SurrealS<M>& );

  void dump( int indentSize=0 ) const;

private:
  double v_;      // value
  double d_[N];   // derivative array
};


// constructors

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const SurrealS& z ) : v_(z.v_)
{
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const double& v0, const double d0[], const int n ) : v_(v0)
{
  assert( n == N );
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const double& v0, const double& d0 ) : v_(v0)
{
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const double& v0 ) : v_(v0)
{
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}


// assignment

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const SurrealS<N>& z )
{
  //Do nothing if assigning self to self
  if ( &z == this ) return *this;

  v_ = z.v_;
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const double& r )
{
  v_ = r;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
  return *this;
}


// unary operators; no side effects

template<int N>
ALWAYS_INLINE const SurrealS<N>&
SurrealS<N>::operator+() const
{
  return *this;
}

template<int N>
ALWAYS_INLINE const SurrealS<N>
SurrealS<N>::operator-() const
{
  SurrealS<N> c;
  c.v_ = -v_;
  for (int i = 0; i < N; i++)
    c.d_[i] = -d_[i];

  return c;
}


// binary accumulation operators

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator+=( const SurrealS<N>& z )
{
  v_ += z.v_;
  for (int i = 0; i < N; i++)
    d_[i] += z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator+=( const double& r )
{
  v_ += r;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator-=( const SurrealS<N>& z )
{
  v_ -= z.v_;
  for (int i = 0; i < N; i++)
    d_[i] -= z.d_[i];
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator-=( const double& r )
{
  v_ -= r;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator*=( const SurrealS<N>& z )
{
  for (int i = 0; i < N; i++)
    d_[i] = v_*z.d_[i] + d_[i]*z.v_;
  v_ *= z.v_;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator*=( const double& r )
{
  for (int i = 0; i < N; i++)
    d_[i] *= r;
  v_ *= r;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator/=( const SurrealS<N>& z)
{
  double tmp = 1./(z.v_*z.v_);
  for (int i = 0; i < N; i++)
    d_[i] = (z.v_*d_[i] - v_*z.d_[i])*tmp;
  v_ /= z.v_;
  return *this;
}

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator/=( const double& r )
{
  double tmp = 1./r;
  for (int i = 0; i < N; i++)
    d_[i] *= tmp;
  v_ *= tmp;
  return *this;
}


// debug dump of private data
template<int N>
void
SurrealS<N>::dump( int indentSize ) const
{
  std::string indent(indentSize, ' ');
  std::cout << indent << "SurrealS<" << N << ">: v_ = " << v_;
  std::cout << "  d_[" << N << "] = (";
  for (int n = 0; n < N-1; n++)
    std::cout << d_[n] << ",";
  std::cout << d_[N-1] << ")" << std::endl;
}


// binary operators

template<int N>
ALWAYS_INLINE SurrealS<N>
operator+( const SurrealS<N>& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] + b.d_[i];
  c.v_ = a.v_ + b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator+( const SurrealS<N>& a, const double& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i];
  c.v_ = a.v_ + b;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator+( const double& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = b.d_[i];
  c.v_ = a + b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator-( const SurrealS<N>& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i] - b.d_[i];
  c.v_ = a.v_ - b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator-( const SurrealS<N>& a, const double& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i];
  c.v_ = a.v_ - b;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator-( const double& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = -b.d_[i];
  c.v_ = a - b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator*( const SurrealS<N>& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.v_*b.d_[i] + a.d_[i]*b.v_;
  c.v_ = a.v_*b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator*( const SurrealS<N>& a, const double& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*b;
  c.v_ = a.v_*b;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator*( const double& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = a*b.d_[i];
  c.v_ = a*b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator/( const SurrealS<N>& a, const SurrealS<N>& b )
{
  SurrealS<N> c;
  double tmp = 1./(b.v_*b.v_);
  for (int i = 0; i < N; i++)
    c.d_[i] = (b.v_*a.d_[i] - a.v_*b.d_[i])*tmp;
  c.v_ = a.v_/b.v_;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N> operator/( const SurrealS<N>& a, const double& b )
{
  SurrealS<N> c;
  double tmp = 1./b;
  for (int i = 0; i < N; i++)
    c.d_[i] = a.d_[i]*tmp;
  c.v_ = a.v_/b;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
operator/( const double& a, const SurrealS<N>& b )
{
  double tmpv = a/(b.v_);
  double tmpd = -1./(b.v_);
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmpv*tmpd*b.d_[i];
  c.v_ = tmpv;
  return c;
}


// relational operators

template<int N>
ALWAYS_INLINE bool
operator==( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ == rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator==( const SurrealS<N>& lhs, const double& rhs )
{
  return lhs.v_ == rhs;
}

template<int N>
ALWAYS_INLINE bool
operator==( const double& lhs, const SurrealS<N>& rhs )
{
  return lhs == rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator!=( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ != rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator!=( const SurrealS<N>& lhs, const double& rhs )
{
  return lhs.v_ != rhs;
}

template<int N>
ALWAYS_INLINE bool
operator!=( const double& lhs, const SurrealS<N>& rhs )
{
  return lhs != rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator>( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ > rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator>( const SurrealS<N>& lhs, const double& rhs )
{
  return lhs.v_ > rhs;
}

template<int N>
ALWAYS_INLINE bool
operator>( const double& lhs, const SurrealS<N>& rhs )
{
  return lhs > rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator<( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ < rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator<( const SurrealS<N>& lhs, const double& rhs )
{
  return lhs.v_ < rhs;
}

template<int N>
ALWAYS_INLINE bool
operator<( const double& lhs, const SurrealS<N>& rhs )
{
  return lhs < rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator>=( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ >= rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator>=( const SurrealS<N>& lhs, const double& rhs )
{
  return lhs.v_ >= rhs;
}

template<int N>
ALWAYS_INLINE bool
operator>=( const double& lhs, const SurrealS<N>& rhs )
{
  return lhs >= rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator<=( const SurrealS<N>& lhs, const SurrealS<N>& rhs )
{
  return lhs.v_ <= rhs.v_;
}

template<int N>
ALWAYS_INLINE bool
operator<=( const SurrealS<N>& lhs, const double& rhs )
{
  return lhs.v_ <= rhs;
}

template<int N>
ALWAYS_INLINE bool
operator<=( const double& lhs, const SurrealS<N>& rhs )
{
  return lhs <= rhs.v_;
}

//Macros for functions

#define SURREALS_FUNC1( NAME, FUNC, DERIV ) \
template<int N> \
ALWAYS_INLINE SurrealS<N> \
NAME( const SurrealS<N>& z ) \
{ \
  double tmp = DERIV; \
  SurrealS<N> c; \
  for (int i = 0; i < N; i++) \
    c.d_[i] = tmp*z.d_[i]; \
  c.v_ = FUNC; \
  return c; \
}

#define SURREALS_FUNC2( NAME, FUNC, DERIV ) \
template<int N> \
ALWAYS_INLINE SurrealS<N> \
NAME( const SurrealS<N>& z1, const SurrealS<N>& z2) \
{ \
  double tmp = DERIV; \
  SurrealS<N> c; \
  for (int i = 0; i < N; i++) \
    c.d_[i] = tmp*(z2.v_*z1.d_[i] - z1.v_*z2.d_[i]); \
  c.v_ = FUNC; \
  return c; \
}

// trig functions <cmath>

SURREALS_FUNC1( cos, std::cos(z.v_), -std::sin(z.v_) )
SURREALS_FUNC1( sin, std::sin(z.v_),  std::cos(z.v_) )
SURREALS_FUNC1( tan, std::tan(z.v_),  double(1)/(std::cos(z.v_)*std::cos(z.v_)) )
SURREALS_FUNC1( acos, std::acos(z.v_), -double(1)/std::sqrt(1 - z.v_*z.v_) )
SURREALS_FUNC1( asin, std::asin(z.v_),  double(1)/std::sqrt(1 - z.v_*z.v_) )
SURREALS_FUNC1( atan, std::atan(z.v_),  double(1)/(1 + z.v_*z.v_) )

SURREALS_FUNC2( atan2, std::atan2(z1.v_, z2.v_),  double(1)/(z1.v_*z1.v_ + z2.v_*z2.v_) )

// hyperbolic functions <cmath>

SURREALS_FUNC1( cosh, std::cosh(z.v_), std::sinh(z.v_) )
SURREALS_FUNC1( sinh, std::sinh(z.v_), std::cosh(z.v_) )
SURREALS_FUNC1( tanh, std::tanh(z.v_), double(1)/(std::cosh(z.v_)*std::cosh(z.v_)) )

// exp and log functions <cmath>

SURREALS_FUNC1( exp, std::exp(z.v_), std::exp(z.v_) )
SURREALS_FUNC1( log, std::log(z.v_), double(1)/z.v_ )
SURREALS_FUNC1( log10, std::log10(z.v_), double(1)/(z.v_*std::log(10.)) )

// power functions <cmath>

template<int N>
ALWAYS_INLINE SurrealS<N>
pow( const SurrealS<N>& a, const SurrealS<N>& b)
{
  // many sticky points were derivative is undefined or infinite
  // badness if 0 <= b < 1 and a == 0
  double powab=std::pow(a.v_,b.v_);
  double tmp1 = b.v_*std::pow(a.v_, b.v_ - 1);
  double tmp2 = powab*std::log(a.v_);
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp1*a.d_[i] + tmp2*b.d_[i];
  c.v_ = powab;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
pow( const SurrealS<N>& a, const double& b)
{
  double tmp = b*std::pow(a.v_, b-1);
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*a.d_[i];
  c.v_ = std::pow(a.v_,b);
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
pow( const SurrealS<N>& a, const int& b)
{
  double tmp = static_cast<double>(b)*std::pow(a.v_, static_cast<double>(b - 1));
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*a.d_[i];
  c.v_ = std::pow(a.v_, static_cast<double>(b));
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
pow( const double& a, const SurrealS<N>& b)
{
  double powab=std::pow(a, b.v_);
  double tmp = powab*std::log(a);
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*b.d_[i];
  c.v_ = powab;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
pow( const int& a, const SurrealS<N>& b)
{
  double powab=std::pow(static_cast<double>(a), b.v_);
  double tmp = powab*std::log(static_cast<double>(a));
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = tmp*b.d_[i];
  c.v_ = powab;
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
sqrt( const SurrealS<N>& z )
{
  double sqrtv=std::sqrt(z.v_);
  if (sqrtv == 0)
  {
    return SurrealS<N>(0, 0);
  }
  else
  {
    double tmp = 0.5/sqrtv;
    SurrealS<N> c;
    for (int i = 0; i < N; i++)
      c.d_[i] = tmp*z.d_[i];
    c.v_ = sqrtv;
    return c;
  }
}


// rounding functions <cmath>

template<int N>
ALWAYS_INLINE SurrealS<N>
ceil( const SurrealS<N>& z )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = 0;
  c.v_ = std::ceil(z.v_);
  return c;
}

template<int N>
ALWAYS_INLINE SurrealS<N>
floor( const SurrealS<N>& z )
{
  SurrealS<N> c;
  for (int i = 0; i < N; i++)
    c.d_[i] = 0;
  c.v_ = std::floor(z.v_);
  return c;
}


// misc functions <cmath>

template<int N>
ALWAYS_INLINE SurrealS<N>
abs( const SurrealS<N>& z )
{
  return (z.v_ < 0) ? -z : SurrealS<N>(z);
}

template<int N>
ALWAYS_INLINE SurrealS<N>
fabs( const SurrealS<N>& z )
{
  return (z.v_ < 0) ? -z : SurrealS<N>(z);
}


// I/O

template<int N>
std::istream&
operator>>( std::istream& is, SurrealS<N>& z )
{
  double v = 0;
  double d[10] = {0};
  char c = 0;
  int n = 0;

  is >> c;
  if (c == '(')
  {
    is >> v;

    is >> c;
    bool done = false;
    while (! done)
    {
      if (c != ')') is.clear(std::ios::badbit);
      if (c == ',')
      {
        is >> d[n]; n++;
      }
      else if (c == ')')
      {
        done = true;
      }
    }
  }
  else
  {
    is.putback(c);
    is >> v;
  }

  if (is) z = SurrealS<N>(v, d, n);
  return is;
}

template<int N>
std::ostream&
operator<<( std::ostream& os, const SurrealS<N>& z )
{
  os << '(' << z.value() << ';';
  for (int i = 0; i < N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(N - 1) << ')';
  return os;
}


//Clean up macro definitions
#undef SURREALS_FUNC1
#undef SURREALS_FUNC2


#endif // SURREALS_TRAD_H
