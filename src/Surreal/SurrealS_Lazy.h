#ifndef SURREALS_LAZY_H
#define SURREALS_LAZY_H

//  overloaded derivative operator
//  ref: derivify.h (Google it)

#include <cmath>
#include <iostream>
#include <cassert>

#include <string>

#include "always_inline.h"

class SurrealSTypeBase {};

template< class Derived >
struct SurrealSType : SurrealSTypeBase
{
  //A convenient method for casting to the derived type
  ALWAYS_INLINE const Derived& cast() const { return static_cast<const Derived&>(*this); }

  //A simple way to call value without having to case first
  ALWAYS_INLINE double value() const { return cast().value(); }
};

//Forward declarations
template<int N>
class SurrealS;

//#define SURREALS_LOOP_UNROLL

namespace SurrealSExpr
{
template<class L, class R > class OpMul;
class OpMul_impl;
}


//----------------------------------------------------------------------------//
// SurrealS:  value, N derivatives
//
// Operators with Lazy Expressions
//
// statically defined derivative array
//----------------------------------------------------------------------------//
template<int N_>
class SurrealS : public SurrealSType< SurrealS<N_> >
{
public:
  static const int N = N_;

  //The default constructor is intentionally empty here. This means Surreal is
  //not initialized when declared, which is consistent with regular numbers. This also
  //improves performance.
  ALWAYS_INLINE SurrealS() {}
  SurrealS( const SurrealS& z );
  SurrealS( const double v0 );
  SurrealS( const double v0, const double d0[], int n );
  SurrealS( const double v0, const double& d0 );
  template<class Expr>
  ALWAYS_INLINE
  SurrealS( const SurrealSType<Expr>& r ) : v_(0) { operator=(r); }
  ALWAYS_INLINE ~SurrealS() {}

  ALWAYS_INLINE int size() const { return N; }

  // value accessor operators
  ALWAYS_INLINE       double& value()       { return v_; }
  ALWAYS_INLINE const double& value() const { return v_; }

  // derivative accessor operators
  ALWAYS_INLINE       double& deriv( int i=0 )       { return d_[i]; }
  ALWAYS_INLINE const double& deriv( int i=0 ) const { return d_[i]; }

  // assignment
  SurrealS& operator=( const SurrealS& );
  SurrealS& operator=( const double& );

  template<class Expr> SurrealS& operator= ( const SurrealSType<Expr>& );
  template<class Expr> SurrealS& operator+=( const SurrealSType<Expr>& );
  template<class Expr> SurrealS& operator-=( const SurrealSType<Expr>& );

  // unary operators; no side effects
  const SurrealS& operator+() const;

  // binary accumulation operators
  SurrealS& operator+=( const double& );
  SurrealS& operator-=( const double& );
  SurrealS& operator*=( const SurrealS& );
  SurrealS& operator*=( const double& );
  SurrealS& operator/=( const SurrealS& );
  SurrealS& operator/=( const double& );

#if 0
  // classification functions <cmath>
  friend bool isfinite( const SurrealS& );
  friend bool isinf( const SurrealS& );
  friend bool isnan( const SurrealS& );
#endif

  // input/output
  template<int M> friend std::istream& operator>>( std::istream&, SurrealS<M>& );

protected:
  double d_[N];   // derivative array
  double v_;      // value
};

//Constructors

template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const SurrealS& z )
{
  v_ = z.v_;
  for (int i = 0; i < N; i++)
    d_[i] = z.d_[i];
}
template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const double v0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = 0;
}
template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const double v0, const double d0[], int n )
{
  assert( n == N );
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = d0[i];
}
template<int N>
ALWAYS_INLINE
SurrealS<N>::SurrealS( const double v0, const double& d0 )
{
  v_ = v0;
  for (int i = 0; i < N; i++)
    d_[i] = d0;
}


namespace SurrealSExpr
{

// Lazy expressions


// Addition and Subtraction

template<class L, class R>
class OpAdd : public SurrealSType< OpAdd<L,R> >
{
public:
  static const int N = L::N;
  BOOST_MPL_ASSERT_RELATION( L::N, ==, R::N );

  ALWAYS_INLINE
  OpAdd(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE double value() const { return Ll.value() + Rr.value(); }
  ALWAYS_INLINE double deriv(const int& i) const { return Ll.deriv(i) + Rr.deriv(i); }

  ALWAYS_INLINE const OpAdd&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return Ll.size(); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R>
ALWAYS_INLINE SurrealSExpr::OpAdd<L,R>
operator+(const SurrealSType<L>& Ll, const SurrealSType<R>& Rr)
{
  return SurrealSExpr::OpAdd<L,R>( Ll.cast(), Rr.cast() );
}

namespace SurrealSExpr
{

template<class L, class R>
class OpSub : public SurrealSType< OpSub<L,R> >
{
public:
  static const int N = L::N;
  BOOST_MPL_ASSERT_RELATION( L::N, ==, R::N );

  ALWAYS_INLINE
  OpSub(const L& Ll, const R& Rr) : Ll(Ll), Rr(Rr) {}

  ALWAYS_INLINE double value() const { return Ll.value() - Rr.value(); }
  ALWAYS_INLINE double deriv(const int& i) const { return Ll.deriv(i) - Rr.deriv(i); }

  ALWAYS_INLINE const OpSub&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return Ll.size(); }
private:
  const L& Ll;
  const R& Rr;
};

}

template<class L, class R>
ALWAYS_INLINE SurrealSExpr::OpSub<L, R>
operator-(const SurrealSType<L>& Ll, const SurrealSType<R>& Rr)
{
  return SurrealSExpr::OpSub<L, R>( Ll.cast(), Rr.cast() );
}

//Addition and Subtraction with scalar quantities

namespace SurrealSExpr
{

template<class Expr>
class OpScalar : public SurrealSType< OpScalar<Expr> >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
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
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator+(const SurrealSType<Expr>& e, const double& s)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator+(const double& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), 1, s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator-(const SurrealSType<Expr>& e, const double& s)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), 1, -s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpScalar<Expr>
operator-(const double& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpScalar<Expr>( e.cast(), -1, s );
}


//Multiplication with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class OpMul : public SurrealSType< OpMul<ExprL, ExprR> >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpMul(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value()) {}

  ALWAYS_INLINE double value() const { return eL_val*eR_val; }
  ALWAYS_INLINE double deriv(const int& i) const { return eL_val*eR.deriv(i) + eL.deriv(i)*eR_val; }

  ALWAYS_INLINE const OpMul&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const double eL_val, eR_val;
};

template<class ExprL>
class OpMul<ExprL, double> : public SurrealSType< OpMul<ExprL, double> >
{
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  OpMul(const ExprL& e, const double s) : e(e), s(s) {}

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
ALWAYS_INLINE SurrealSExpr::OpMul<ExprL, ExprR>
operator*(const SurrealSType<ExprL>& z1, const SurrealSType<ExprR>& z2)
{
  return SurrealSExpr::OpMul<ExprL, ExprR>( z1.cast(), z2.cast() );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
operator*(const SurrealSType<Expr>& e, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double>( e.cast(), s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
operator/(const SurrealSType<Expr>& e, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double>( e.cast(), double(1)/s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
operator*(const double& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpMul<Expr, double>( e.cast(), s );
}

//=============================================================================
// This is a special case when multiplies scalars are multiplying from two sides, i.e. B = 2*A*2;
// This reduces the complexity of the expression tree and hence reduces code bloat
template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
operator*(const SurrealSExpr::OpMul<Expr, double>& MulScal, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double>( MulScal.e, MulScal.s*s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
operator/(const SurrealSExpr::OpMul<Expr, double>& MulScal, const double& s)
{
  return SurrealSExpr::OpMul<Expr, double>( MulScal.e, MulScal.s/s );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
operator*(const double& s, const SurrealSExpr::OpMul<Expr, double>& MulScal)
{
  return SurrealSExpr::OpMul<Expr, double>( MulScal.e, MulScal.s*s );
}


//=============================================================================
//Division with SurrealSs

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class OpDiv : public SurrealSType< OpDiv<ExprL, ExprR> >
{
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION( ExprL::N, ==, ExprR::N );

  ALWAYS_INLINE
  OpDiv(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), eL_val(eL.value()), eR_val(eR.value())
                                          , vali(1/(eR_val*eR_val)) {}

  ALWAYS_INLINE double value() const { return eL_val/eR_val; }
  ALWAYS_INLINE double deriv(const int& i) const { return (eR_val*eL.deriv(i) - eR.deriv(i)*eL_val)*vali; }

  ALWAYS_INLINE const OpDiv&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const double eL_val, eR_val, vali;
};

}

template< class ExprL, class ExprR >
ALWAYS_INLINE SurrealSExpr::OpDiv<ExprL, ExprR>
operator/(const SurrealSType<ExprL>& eL, const SurrealSType<ExprR>& eR)
{
  return SurrealSExpr::OpDiv<ExprL, ExprR>( eL.cast(), eR.cast() );
}


namespace SurrealSExpr
{

template<class Expr>
class OpDivScalarNumerator : public SurrealSType< OpDivScalarNumerator<Expr> >
{
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
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
ALWAYS_INLINE SurrealSExpr::OpDivScalarNumerator<Expr>
operator/(const double& s, const SurrealSType<Expr>& e)
{
  return SurrealSExpr::OpDivScalarNumerator<Expr>( e.cast(), s );
}


// assignment

template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const SurrealS& z )
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

template<int N>
template< class Expr >
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator=( const SurrealSType<Expr>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for (int i = 0; i < N; ++i)
    d_[i] = Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ = Tree.value();

  return *this;
}

template<int N>
template< class Expr >
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator+=( const SurrealSType<Expr>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for (int i = 0; i < N; ++i)
    d_[i] += Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ += Tree.value();

  return *this;
}

template<int N>
template< class Expr >
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator-=( const SurrealSType<Expr>& r )
{
  const Expr& Tree = r.cast();

  BOOST_MPL_ASSERT_RELATION( N, ==, Expr::N );

  for (int i = 0; i < N; ++i)
    d_[i] -= Tree.deriv(i);

  //Value must be set last as it might be used in the derivative calculation
  v_ -= Tree.value();

  return *this;
}


// unary operators; no side effects

template<int N>
ALWAYS_INLINE const SurrealS<N>&
SurrealS<N>::operator+() const
{
  return *this;
}

template< class Expr >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, double>
operator-(SurrealSType<Expr> const& e)
{
  return SurrealSExpr::OpMul<Expr, double>( e.cast(), -1 );
}

template< class Expr >
ALWAYS_INLINE const SurrealSExpr::OpMul<Expr, double>
operator-(SurrealSExpr::OpMul<Expr, double> const& Mul)
{
  return SurrealSExpr::OpMul<Expr, double>( Mul.e, -1*Mul.s );
}

// binary accumulation operators


template<int N>
ALWAYS_INLINE SurrealS<N>&
SurrealS<N>::operator+=( const double& r )
{
  v_ += r;
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
SurrealS<N>::operator*=( const SurrealS& z )
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
SurrealS<N>::operator/=( const SurrealS& z)
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

// relational operators

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator==( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() == rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator==( const SurrealSType<Expr>& lhs, const double& rhs )
{
  return lhs.value() == rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator==( const double& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs == rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator!=( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() != rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator!=( const SurrealSType<Expr>& lhs, const double& rhs )
{
  return lhs.value() != rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator!=( const double& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs != rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator>( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() > rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator>( const SurrealSType<Expr>& lhs, const double& rhs )
{
  return lhs.value() > rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator>( const double& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs > rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator<( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() < rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator<( const SurrealSType<Expr>& lhs, const double& rhs )
{
  return lhs.value() < rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator<( const double& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs < rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator>=( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() >= rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator>=( const SurrealSType<Expr>& lhs, const double& rhs )
{
  return lhs.value() >= rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator>=( const double& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs >= rhs.value();
}

template<class ExprL, class ExprR>
ALWAYS_INLINE bool
operator<=( const SurrealSType<ExprL>& lhs, const SurrealSType<ExprR>& rhs )
{
  return lhs.value() <= rhs.value();
}

template<class Expr>
ALWAYS_INLINE bool
operator<=( const SurrealSType<Expr>& lhs, const double& rhs )
{
  return lhs.value() <= rhs;
}

template<class Expr>
ALWAYS_INLINE bool
operator<=( const double& lhs, const SurrealSType<Expr>& rhs )
{
  return lhs <= rhs.value();
}


//Functions for SurrealSs
#define SURREALS_FUNC1( NAME, FUNC, DERIV ) \
namespace SurrealSExpr \
{  \
template<class Expr> \
class BOOST_PP_CAT(SurrealS_, NAME) : public SurrealSType< SurrealS_ ## NAME<Expr> > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = Expr::N; \
  \
  ALWAYS_INLINE \
  BOOST_PP_CAT(SurrealS_, NAME)(const Expr& e) : e(e), z(e.value()), der(DERIV) {} \
  \
  ALWAYS_INLINE double value() const { return FUNC; } \
  ALWAYS_INLINE double deriv(const int& i) const { return der*e.deriv(i); } \
  \
  ALWAYS_INLINE const BOOST_PP_CAT(SurrealS_, NAME)& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return e.size(); } \
private: \
  const Expr& e; \
  const double z, der; \
}; \
} \
\
template<class Expr> \
ALWAYS_INLINE SurrealSExpr::SurrealS_ ## NAME<Expr> \
NAME(const SurrealSType<Expr>& z) { return SurrealSExpr::SurrealS_ ## NAME<Expr>( z.cast() ); }


#define SURREALS_FUNC2( NAME, FUNC, DERIV ) \
namespace SurrealSExpr \
{  \
template<class ExprL, class ExprR> \
class BOOST_PP_CAT(SurrealS_, NAME) : public SurrealSType< SurrealS_ ## NAME<ExprL, ExprR> > \
{ /*This is for functions when the argument is an expression*/ \
public: \
  static const int N = ExprL::N; \
  BOOST_MPL_ASSERT_RELATION(ExprL::N, ==, ExprR::N ); \
  \
  ALWAYS_INLINE \
  BOOST_PP_CAT(SurrealS_, NAME)(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), z1(eL.value()), z2(eR.value()), \
                                                                 der(DERIV) {} \
  \
  ALWAYS_INLINE double value() const { return FUNC; } \
  ALWAYS_INLINE double deriv(const int& i) const { return der*(z2*eL.deriv(i) - z1*eR.deriv(i)); } \
  \
  ALWAYS_INLINE const BOOST_PP_CAT(SurrealS_, NAME)& \
  operator+() const { return *this; } \
  ALWAYS_INLINE int size() const { return eL.size(); } \
private: \
  const ExprL& eL; \
  const ExprR& eR; \
  const double z1, z2, der; \
}; \
  \
} \
\
template<class ExprL, class ExprR> \
ALWAYS_INLINE SurrealSExpr::SurrealS_ ## NAME<ExprL, ExprR> \
NAME(const SurrealSType<ExprL>& z1, const SurrealSType<ExprR>& z2) \
{ return SurrealSExpr::SurrealS_ ## NAME<ExprL, ExprR>( z1.cast(), z2.cast() ); }

// trig functions <cmath>

SURREALS_FUNC1( cos, std::cos(z), -std::sin(z) )
SURREALS_FUNC1( sin, std::sin(z),  std::cos(z) )
SURREALS_FUNC1( tan, std::tan(z),  double(1)/(std::cos(z)*std::cos(z)) )
SURREALS_FUNC1( acos, std::acos(z), -double(1)/std::sqrt(1 - z*z) )
SURREALS_FUNC1( asin, std::asin(z),  double(1)/std::sqrt(1 - z*z) )
SURREALS_FUNC1( atan, std::atan(z),  double(1)/(1 + z*z) )

SURREALS_FUNC2( atan2, std::atan2(z1, z2),  double(1)/(z1*z1 + z2*z2) )

// hyperbolic functions <cmath>

SURREALS_FUNC1( cosh, std::cosh(z), std::sinh(z) )
SURREALS_FUNC1( sinh, std::sinh(z), std::cosh(z) )
SURREALS_FUNC1( tanh, std::tanh(z), double(1)/(std::cosh(z)*std::cosh(z)) )

// exp and log functions <cmath>

SURREALS_FUNC1( exp, std::exp(z), std::exp(z) )
SURREALS_FUNC1( log, std::log(z), double(1)/z )
SURREALS_FUNC1( log10, std::log10(z), double(1)/(z*std::log(10.)) )

// power functions <cmath>

namespace SurrealSExpr
{

template<class ExprL, class ExprR>
class SurrealS_pow : public SurrealSType< SurrealS_pow<ExprL, ExprR> >
{ /*This is for functions when the argument is an expression*/
public:
  static const int N = ExprL::N;
  BOOST_MPL_ASSERT_RELATION(ExprL::N, ==, ExprR::N);

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const ExprR& eR) : eL(eL), eR(eR), a(eL.value()), b(eR.value()),
                                                   powab(std::pow(a,b)),
                                                   tmp1( b*std::pow(a, b - 1) ),
                                                   tmp2( powab*std::log(a) ) {}

  ALWAYS_INLINE double value() const { return powab; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp1*eL.deriv(i) + tmp2*eR.deriv(i); }

  ALWAYS_INLINE const SurrealS_pow&
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const ExprR& eR;
  const double a, b, powab, tmp1, tmp2;
};

template<class ExprL>
class SurrealS_pow<ExprL, double> : public SurrealSType< SurrealS_pow<ExprL, double> >
{ /*This is optimized when the argument is SurrealS and double*/
public:
  static const int N = ExprL::N;

  ALWAYS_INLINE
  SurrealS_pow(const ExprL& eL, const double& b) : eL(eL), a(eL.value()),
                                                 powab(std::pow(a,b)),
                                                 tmp1( b*std::pow(a, b - 1) ) {}

  ALWAYS_INLINE double value() const { return powab; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp1*eL.deriv(i); }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eL.size(); }
private:
  const ExprL& eL;
  const double a, powab, tmp1;
};


template<class ExprR>
class SurrealS_pow<double, ExprR> : public SurrealSType< SurrealS_pow<double, ExprR> >
{ /*This is optimized when the argument is a double and SurrealS*/
public:
  static const int N = ExprR::N;

  ALWAYS_INLINE
  SurrealS_pow(const double& a, const ExprR& eR) : eR(eR), b(eR.value()),
                                               powab(std::pow(a,b)),
                                               tmp2( powab*std::log(a) ) {}

  ALWAYS_INLINE double value() const { return powab; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp2*eR.deriv(i); }

  template<int I>
  ALWAYS_INLINE double get_deriv() const { return tmp2*eR.template get_deriv<I>(); }

  ALWAYS_INLINE const SurrealS_pow
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return eR.size(); }
private:
  const ExprR& eR;
  const double b, powab, tmp2;

};

}

template<class ExprL, class ExprR>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<ExprL, ExprR>
pow(const SurrealSType<ExprL>& a, const SurrealSType<ExprR>& b)
{
  return SurrealSExpr::SurrealS_pow<ExprL, ExprR>( a.cast(), b.cast() );
}

template<class Expr, typename T>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<Expr, double>
pow(const SurrealSType<Expr>& a, const double& b )
{
  return SurrealSExpr::SurrealS_pow<Expr, double>( a.cast(), b );
}

template<class Expr, typename T>
ALWAYS_INLINE SurrealSExpr::SurrealS_pow<double, Expr>
pow(const double& a, const SurrealSType<Expr>& b)
{
  return SurrealSExpr::SurrealS_pow<double, Expr>( a, b.cast() );
}


namespace SurrealSExpr
{

template<class Expr>
class SurrealS_sqrt : public SurrealSType< SurrealS_sqrt<Expr> >
{ /*This is optimized when the argument is an Expression*/
public:
  static const int N = Expr::N;

  ALWAYS_INLINE
  SurrealS_sqrt(const Expr& e) : e(e), sqrtv( sqrt(e.value()) ), tmp( sqrtv == 0 ? 0 : 0.5/sqrtv ) {}

  ALWAYS_INLINE double value() const { return sqrtv; }
  ALWAYS_INLINE double deriv(const int& i) const { return tmp*e.deriv(i); }

  template<int I>
  ALWAYS_INLINE double get_deriv() const { return tmp*e.template get_deriv<I>(); }

  ALWAYS_INLINE const SurrealS_sqrt
  operator+() const { return *this; }
  ALWAYS_INLINE int size() const { return e.size(); }
private:
  const Expr& e;
  const double sqrtv, tmp;
};

}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::SurrealS_sqrt<Expr>
sqrt(const SurrealSType<Expr>& z)
{
  return SurrealSExpr::SurrealS_sqrt<Expr>( z.cast() );
}


// rounding functions <cmath>

SURREALS_FUNC1( ceil, std::ceil(z), 0 )
SURREALS_FUNC1( floor, std::floor(z), 0 )

// misc functions <cmath>

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
abs( const SurrealSType<Expr>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<Expr, double>( z.cast(), -1 ) :
         SurrealSExpr::OpMul<Expr, double>( z.cast(),  1 );
}

template<class Expr>
ALWAYS_INLINE SurrealSExpr::OpMul<Expr, double>
fabs( const SurrealSType<Expr>& z )
{
  return (z.value() < 0) ?
         SurrealSExpr::OpMul<Expr, double>( z.cast(), -1 ) :
         SurrealSExpr::OpMul<Expr, double>( z.cast(),  1 );
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


template<class Expr>
std::ostream&
operator<<( std::ostream& os, const SurrealSType<Expr>& ztype )
{
  const Expr& z = ztype.cast();
  os << '(' << z.value() << ';';
  for (int i = 0; i < Expr::N - 1; i++)
    os << z.deriv(i) << ',';
  os << z.deriv(Expr::N - 1) << ')';
  return os;
}

//Clean up macro definitions
#undef SURREALS_FUNC1
#undef SURREALS_FUNC2

#endif // SURREALS_LAZY_H
