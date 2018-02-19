/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Spline Fitting Functions w/ Derivatives
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egads.h"
#include "Surreal/SurrealS.h"
#include "Surreal/SurrealD.h"

#ifdef WIN32
#define DllExport   __declspec( dllexport )
#else
#define DllExport
#endif

#define NITER    10000
#define RELAX    0.15
#define MAXDEG   21
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)          (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))



/* to deal with tolerance & print statements */

static double value(double val)
{
  return val;
}

template <int N>
static double value(SurrealS<N> valS)
{
  return valS.value();
}


template<class T, class T2>
static int
FindSpan(int nKnots, int degree, T2 u, T *U)
{
  int n, low, mid, high;

  if (u <= U[degree]) return degree;
  n = nKnots - degree - 1;
  if (u >= U[n]) return n-1;

  low  = degree;
  high = n;
  mid  = (low+high)/2;
  while ((u < U[mid]) || (u >= U[mid+1])) {
    if (u < U[mid]) {
      high = mid;
    } else {
      low  = mid;
    }
    mid = (low+high)/2;
  }

  return mid;
}


template<class T, class T2>
static void
BasisFuns(int span, int degree, T2 u, T *U, T *N)
{
  int j, r;
  T   saved, temp, left[MAXDEG], right[MAXDEG];

  N[0] = 1.0;
  for (j = 1; j <= degree; j++) {
    left[j]  = u - U[span+1-j];
    right[j] = U[span+j] - u;
    saved    = 0.0;
    for (r = 0; r < j; r++) {
      temp  = N[r]/(right[r+1]+left[j-r]);
      N[r]  = saved + right[r+1]*temp;
      saved = left[j-r]*temp;
    }
    N[j] = saved;
  }
}


template<class T, class T2>
static void
DersBasisFuns(int i, int p, T2 u, T *knot, int der, T **ders)
{
  int j, k, j1, j2, r, s1, s2, rk, pk;
  T   d, saved, temp, ndu[MAXDEG+1][MAXDEG+1];
  T   a[2][MAXDEG+1], left[MAXDEG+1], right[MAXDEG+1];

  ndu[0][0] = 1.0;
  for (j = 1; j <= p; j++) {
    left[j]  = u - knot[i+1-j];
    right[j] = knot[i+j] - u;
    saved = 0.0;
    for (r = 0; r < j; r++) {
      ndu[j][r] = right[r+1] + left[j-r];
      temp      = ndu[r][j-1]/ndu[j][r];
      ndu[r][j] = saved + right[r+1]*temp;
      saved     = left[j-r]*temp;
    }
    ndu[j][j] = saved;
  }

  for (j = 0; j <= p; j++) ders[0][j] = ndu[j][p];  /* basis function */

  /* compute derivatives */
  for (r = 0; r <= p; r++ ) {
    s1      = 0;
    s2      = 1;
    a[0][0] = 1.0;
    /* compute k'th derivative */
    for (k = 1; k <= der; k++ ) {
      d  = 0.0;
      rk = r - k;
      pk = p - k;
      if (r >= k) {
	a[s2][0] = a[s1][0]/ndu[pk+1][rk];
	d        = a[s2][0]*ndu[rk][pk];
      }
      j1 = rk >= -1  ? 1   : -rk;
      j2 = (r-1<=pk) ? k-1 : p-r;
      for (j = j1; j <= j2; j++) {
	a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
	d       +=  a[s2][j]*ndu[rk+j][pk];
      }
      if (r <= pk) {
	a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
	d       +=  a[s2][k]*ndu[r][pk];
      }
      ders[k][r] = d;
      /* switch rows */
      j  = s1;
      s1 = s2;
      s2 = j;
    }
  }

  r = p;
  for (k = 1; k <= der; k++) {
    for (j = 0; j <= p; j++ ) ders[k][j] *= r;
    r *= p - k;
  }

}


/* ************************** Curve Functions ****************************** */

template<class T, class T2>
static int
EG_spline1dEval_impl(int *ivec, T *data, T2 t, T *point)
{
  int i, j, degree, nKnots, span;
  T   N[MAXDEG], *CP;

  point[0] = point[1] = point[2] = 0.0;
  degree   = ivec[1];
  nKnots   = ivec[3];
  CP       = data + ivec[3];
  if (ivec[0] != 0) {
    printf(" EG_spline1dEval: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degree >= MAXDEG) {
    printf(" EG_spline1dEval: degree %d >= %d!\n", degree, MAXDEG);
    return EGADS_CONSTERR;
  }

  span = FindSpan(nKnots, degree, t, data);
  BasisFuns(span, degree, t, data, N);

  for (j = 0; j <= degree;  j++) {
    i         = span-degree+j;
    point[0] += N[j]*CP[3*i  ];
    point[1] += N[j]*CP[3*i+1];
    point[2] += N[j]*CP[3*i+2];
  }

  return EGADS_SUCCESS;
}


extern "C"
int
EG_spline1dEval(int *ivec, double *data, double t, double *point)
{
  return EG_spline1dEval_impl(ivec, data, t, point);
}


extern "C"
int
EG_spline1dEval_psens(int *ivec, const double *data, const double *rdata_xyz,
                       const double *xyz_dot, double t,
                       double *point, double *point_dot)
{
  SurrealS<1>* dataS = NULL;
  SurrealS<1>  pointS[3];

  point[0]     = point[1]     = point[2]     = 0.0;
  point_dot[0] = point_dot[1] = point_dot[2] = 0.0;
  int icp      = ivec[2];
  int iknot    = ivec[3];
  int imax     = icp-2;

  dataS = new SurrealS<1>[iknot+3*icp];
  if (dataS == NULL) return EGADS_MALLOC;

  for (int i = 0; i < iknot+3*icp; i++) {
    dataS[i] = data[i];
    double deriv = 0;
    for (int j = 0; j < 3*imax; j++) {
      deriv += rdata_xyz[j + i*3*imax]*xyz_dot[j];
    }
    dataS[i].deriv() = deriv;
  }

  int stat = EG_spline1dEval_impl(ivec, dataS, t, pointS);
  delete [] dataS;
  if (stat != EGADS_SUCCESS) return stat;

  for (int i = 0; i < 3; i++) {
    point[i]     = pointS[i].value();
    point_dot[i] = pointS[i].deriv();
  }
  return EGADS_SUCCESS;
}


int
EG_spline1d_rdata_dot(int *ivec, const double *rdata, const double *rdata_xyz,
                      const SurrealS<1> *xyz_dot, SurrealS<1> **rdata_dot)
{
  int icp      = ivec[2];
  int iknot    = ivec[3];
  int imax     = icp-2;

  (*rdata_dot) = new SurrealS<1>[iknot+3*icp];
  if ((*rdata_dot) == NULL) return EGADS_MALLOC;

  for (int i = 0; i < iknot+3*icp; i++) {
    (*rdata_dot)[i] = rdata[i];
    double deriv = 0;
    for (int j = 0; j < 3*imax; j++) {
      deriv += rdata_xyz[j + i*3*imax]*xyz_dot[j].deriv();
    }
    (*rdata_dot)[i].deriv() = deriv;
  }

  return EGADS_SUCCESS;
}


template<int N>
int
EG_spline1dEval(int *ivec, SurrealS<N> *data, double t, SurrealS<N> *point)
{
  return EG_spline1dEval_impl(ivec, data, t, point);
}

// Create explicit instantiations of the function
template DllExport int EG_spline1dEval<1>(int *, SurrealS<1> *, double,
                                          SurrealS<1> *);


template<class T, class T2>
static int
EG_spline1dDeriv_impl(int *ivec, T *data, int der, T2 t, T *deriv)
{
  int i, j, k, degree, nKnots, span, dt;
  T   Nders[MAXDEG+1][MAXDEG+1], *Nder[MAXDEG+1], *CP;

  for (k = 0; k <= der; k++) deriv[3*k  ] = deriv[3*k+1] = deriv[3*k+2] = 0.0;

  degree = ivec[1];
  nKnots = ivec[3];
  CP     = data + ivec[3];
  dt     = MIN(der, degree);
  if (ivec[0] != 0) {
    printf(" EG_spline1dDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degree >= MAXDEG) {
    printf(" EG_spline1dDeriv: degree %d >= %d!\n", degree, MAXDEG);
    return EGADS_CONSTERR;
  }
  for (i = 0; i <= degree; i++) Nder[i] = &Nders[i][0];

  span = FindSpan(nKnots, degree, t, data);
  DersBasisFuns(span, degree, t, data, dt, Nder);

  for (k = 0; k <= dt; k++)
    for (j = 0; j <= degree; j++) {
      i             = span-degree+j;
      deriv[3*k  ] += Nders[k][j]*CP[3*i  ];
      deriv[3*k+1] += Nders[k][j]*CP[3*i+1];
      deriv[3*k+2] += Nders[k][j]*CP[3*i+2];
    }
  return EGADS_SUCCESS;
}


extern "C"
int
EG_spline1dDeriv(int *ivec, double *data, int der, double t, double *deriv)
{
  return EG_spline1dDeriv_impl(ivec, data, der, t, deriv);
}


/*
 ************************************************************************
 *                                                                      *
 *   EG_spline1dFit - create 1d cubic spline from input data            *
 *                                                                      *
 ************************************************************************
 */

template<class T>
static int
EG_spline1dFit_impl(int endx, int imaxx, const T *xyz, const double *kn,
                    double tol, int *header, T *rdata)
{
    int i, kk, iknot, icp, iter, endc, imax, *mdata;
    T   du, dx, dy, dz, dxyzmax, rj[3], u21, u20, data[9];
    T   d2xdt2L, d2ydt2L, d2zdt2L, d2sdt2L, d2xdt2R, d2ydt2R, d2zdt2R, d2sdt2R;
    T   *cp, *knots, *cps;

    endc = endx;
    imax = imaxx;
    if (imax < 0)                   imax = -imax;
    if (imax < 2)                   return EGADS_DEGEN;
    if ((endx <  0) || (endx >  2)) return EGADS_RANGERR;
    if ((imax == 2) && (endc == 2)) endc = 1;

    /* indices associated with various arrays:
              xyz       knot       cp
                         0*
                         1*
               0         2*         0
               .         3*         1       (for end condition)
               1         4          2
               2         5          3

                         :

             imax-2    imax+1     imax-1
               .       imax+2*    imax      (for end condition)
             imax-1    imax+3*    imax+1
                       imax+4*
                       imax+5*

        *note: there are 4 repeated knots at beginning and
                         4 repeated knots at end */

    mdata = (int *) EG_alloc(imax*sizeof(int));
    if (mdata == NULL) return EGADS_MALLOC;
    icp   = imax + 2;
    iknot = imax + 6;
    knots =  rdata;
    cps   = &rdata[iknot];
    cp    = &rdata[iknot+3*icp];
  
    /* create spline curve */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;

    /* knots */
    if (kn == NULL) {
      kk          = 0;
      knots[kk++] = 0.0;
      knots[kk++] = 0.0;
      knots[kk++] = 0.0;

      knots[kk++] = 0.0;
      if (imaxx > 0) {
        /* arc-length spaced */
        for (i = 1; i < imax; i++) {
          knots[kk] = knots[kk-1] +
                      sqrt((xyz[3*i  ]-xyz[3*i-3])*(xyz[3*i  ]-xyz[3*i-3]) +
                           (xyz[3*i+1]-xyz[3*i-2])*(xyz[3*i+1]-xyz[3*i-2]) +
                           (xyz[3*i+2]-xyz[3*i-1])*(xyz[3*i+2]-xyz[3*i-1]));
          kk++;
        }
      } else {
        /* equally spaced */
        for (i = 1; i < imax; i++) knots[kk++] = i;
      }

      knots[kk] = knots[kk-1]; kk++;
      knots[kk] = knots[kk-1]; kk++;
      knots[kk] = knots[kk-1]; kk++;

      /* normalize */
      for (i = 0; i < kk; i++) knots[i] /= knots[kk-1];
    } else {
      kk          = 0;
      knots[kk++] = 0.0;
      knots[kk++] = 0.0;
      knots[kk++] = 0.0;
      for (i = 0; i < iknot; i++, kk++) knots[kk] = kn[i];
      knots[kk] = knots[kk-1]; kk++;
      knots[kk] = knots[kk-1]; kk++;
      knots[kk] = knots[kk-1]; kk++;
      for (i = 0; i < kk; i++) knots[i] /= knots[kk-1];
    }
  
    /* look for knots with multiplicity 2 or 3 */
    mdata[0] = mdata[imax-1] = 1;
    for (i = 1; i < imax-1; i++)
        if ((knots[i+3] == knots[i+4]) && (knots[i+3] == knots[i+5])) {
            if (knots[i+3] == knots[i+6]) {
                EG_free(mdata);
                return EGADS_DEGEN;
            }
#ifdef DEBUG
            printf("repeated knot at i=%d, %d, and %d\n", i, i+1, i+2);
#endif
            mdata[i++] = -3;
            mdata[i++] =  1;
            mdata[i  ] = +3;
        } else if (knots[i+3] == knots[i+4]) {
#ifdef DEBUG
            printf("repeated knot at i=%d and %d\n", i, i+1);
#endif
            mdata[i++] = -2;
            mdata[i  ] = +2;
        } else {
            mdata[i  ] =  1;
        }
#ifdef DEBUG
    printf("after pass 1\n");
    for (i = 0; i < imax; i++)
        printf("mdata[%2d]=%2d\n", i, mdata[i]);
#endif
  
    /* for all multiplicity 2 knots, determine where the flat spot is */
    for (i = 1; i < imax-1; i++)
        if ((mdata[i] == -2) && (mdata[i+1] == +2)) {
            if (i < 2) {
                /* flat on left */
                mdata[i+1] = 1;
            } else if (i > imax-4) {
                /* flat on right */
                mdata[i  ] = 1;
            } else if ((knots[i+2] == knots[i+1]) && (knots[i+6] == knots[i+5])) {
                /* flat on both sides */
                printf("EG_spline1dFit: Not enough on either side of %d and %d\n",
                       i, i+1);
                printf("                This is probably an error!\n");
                mdata[i  ] = 1;
                mdata[i+1] = 1;
            } else if (knots[i+2] == knots[i+1]) {
                /* flat on left because of another multiplicity 2 too close */
#ifdef DEBUG
                printf("flat on left at %d and %d due to close data\n", i, i+1);
#endif
                mdata[i+1] = 1;
            } else if (knots[i+6] == knots[i+5]) {
                /* flat on right because of another multiplicity 2 too close */
#ifdef DEBUG
                printf("flat on right at %d and %d due to close data\n", i, i+1);
#endif
                mdata[i  ] = 1;
            } else {
                /* find second derivatives on both sides */
                d2xdt2L = ((xyz[3*i   ]-xyz[3*i-3])/(knots[i+3]-knots[i+2]) -
                           (xyz[3*i- 3]-xyz[3*i-6])/(knots[i+2]-knots[i+1])) /
                                                    (knots[i+3]-knots[i+1]);
                d2ydt2L = ((xyz[3*i+ 1]-xyz[3*i-2])/(knots[i+3]-knots[i+2]) -
                           (xyz[3*i- 2]-xyz[3*i-5])/(knots[i+2]-knots[i+1])) /
                                                    (knots[i+3]-knots[i+1]);
                d2zdt2L = ((xyz[3*i+ 2]-xyz[3*i-1])/(knots[i+3]-knots[i+2]) -
                           (xyz[3*i- 1]-xyz[3*i-4])/(knots[i+2]-knots[i+1])) /
                                                    (knots[i+3]-knots[i+1]);
        
                d2xdt2R = ((xyz[3*i+ 9]-xyz[3*i+6])/(knots[i+6]-knots[i+5]) -
                           (xyz[3*i+ 6]-xyz[3*i+3])/(knots[i+5]-knots[i+4])) /
                                                    (knots[i+6]-knots[i+4]);
                d2ydt2R = ((xyz[3*i+10]-xyz[3*i+7])/(knots[i+6]-knots[i+5]) -
                           (xyz[3*i+ 7]-xyz[3*i+4])/(knots[i+5]-knots[i+4])) /
                                                    (knots[i+6]-knots[i+4]);
                d2zdt2R = ((xyz[3*i+11]-xyz[3*i+8])/(knots[i+6]-knots[i+5]) -
                           (xyz[3*i+ 8]-xyz[3*i+5])/(knots[i+5]-knots[i+4])) /
                                                    (knots[i+6]-knots[i+4]);
        
                d2sdt2L = d2xdt2L*d2xdt2L + d2ydt2L*d2ydt2L + d2zdt2L*d2zdt2L;
                d2sdt2R = d2xdt2R*d2xdt2R + d2ydt2R*d2ydt2R + d2zdt2R*d2zdt2R;
                if (d2sdt2L > d2sdt2R) {
                    /* flat on right because left curvature is smaller */
#ifdef DEBUG
                    printf("flat on right\n");
#endif
                    mdata[i  ] = 1;
                } else {
                    /* flat on left because right curvture is smaller  */
#ifdef DEBUG
                    printf("flat on left\n");
#endif
                    mdata[i+1] = 1;
                }
            }
            i++;
        }
#ifdef DEBUG
    printf("after pass 2\n");
    for (i = 0; i < imax; i++)
        printf("mdata[%2d]=%2d\n", i, mdata[i]);
#endif

    /* initial control point */
    kk        = 0;
    cps[kk++] = xyz[0];
    cps[kk++] = xyz[1];
    cps[kk++] = xyz[2];

    /* initial interior control point (for slope) */
    cps[kk++] = (3 * xyz[0] + xyz[3]) / 4;
    cps[kk++] = (3 * xyz[1] + xyz[4]) / 4;
    cps[kk++] = (3 * xyz[2] + xyz[5]) / 4;

    /* interior control points */
    for (i = 1; i < imax-1; i++) {
        cps[kk++] = xyz[3*i  ];
        cps[kk++] = xyz[3*i+1];
        cps[kk++] = xyz[3*i+2];
    }

    /* penultimate interior control point (for slope) */
    cps[kk++] = (3 * xyz[3*(imax-1)  ] + xyz[3*(imax-2)  ]) / 4;
    cps[kk++] = (3 * xyz[3*(imax-1)+1] + xyz[3*(imax-2)+1]) / 4;
    cps[kk++] = (3 * xyz[3*(imax-1)+2] + xyz[3*(imax-2)+2]) / 4;

    /* final control point */
    cps[kk++] = xyz[3*(imax-1)  ];
    cps[kk++] = xyz[3*(imax-1)+1];
    cps[kk++] = xyz[3*(imax-1)+2];

    /* iterate to have knot evaluations match data points */
    for (iter = 0; iter < NITER; iter++) {
        dxyzmax = 0.0;

        /* condition at beginning */
        EG_spline1dDeriv_impl(header, rdata, 2, knots[3], data);
        du = knots[4] - knots[3];
        if (endc == 0) {
            /* natural end */
            dx = du * du * data[6];
            dy = du * du * data[7];
            dz = du * du * data[8];
        } else if (endc == 1) {
            /* FD slope */
            dx = xyz[3] - xyz[0] - du * data[3];
            dy = xyz[4] - xyz[1] - du * data[4];
            dz = xyz[5] - xyz[2] - du * data[5];
        } else {
            /* quadratic fit */
            u20    = knots[5] - knots[3];
            u21    = knots[5] - knots[4];
            rj[0]  = xyz[3]*u20*u20 - xyz[0]*u21*u21 - xyz[6]*du*du;
            rj[1]  = xyz[4]*u20*u20 - xyz[1]*u21*u21 - xyz[7]*du*du;
            rj[2]  = xyz[5]*u20*u20 - xyz[2]*u21*u21 - xyz[8]*du*du;
            rj[0] /= 2.0*u21*du;
            rj[1] /= 2.0*u21*du;
            rj[2] /= 2.0*u21*du;
            dx     = rj[0] - xyz[0] - 0.5 * u20 * data[3];
            dy     = rj[1] - xyz[1] - 0.5 * u20 * data[4];
            dz     = rj[2] - xyz[2] - 0.5 * u20 * data[5];
        }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3] = cps[3] + RELAX * dx;
        cp[4] = cps[4] + RELAX * dy;
        cp[5] = cps[5] + RELAX * dz;

        /* match interior spline points */
        for (i = 1; i < imax-1; i++) {
          
            if (mdata[i] == +2) {
                /* multiplicity 2 with flat spot on right
                   note: i is the second repeated data point */
                EG_spline1dDeriv_impl(header, rdata, 1, knots[i+2], data);
                du =  knots[i+4] - knots[i+3];
                dx = (xyz[3*i+3] - xyz[3*i  ]) - du * data[3];
                dy = (xyz[3*i+4] - xyz[3*i+1]) - du * data[4];
                dz = (xyz[3*i+5] - xyz[3*i+2]) - du * data[5];
            
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            
                cp[3*i+3] = cps[3*i+3] + RELAX * dx;
                cp[3*i+4] = cps[3*i+4] + RELAX * dy;
                cp[3*i+5] = cps[3*i+5] + RELAX * dz;
                continue;
            }
          
            if (mdata[i] == -2) {
                /* multiplicity 2 with flat spot on left
                   note: i is the first repeated data point */
                EG_spline1dDeriv_impl(header, rdata, 1, knots[i+4], data);
                du =  knots[i+3] - knots[i+2];
                dx = (xyz[3*i-3] - xyz[3*i  ]) + du * data[3];
                dy = (xyz[3*i-2] - xyz[3*i+1]) + du * data[4];
                dz = (xyz[3*i-1] - xyz[3*i+2]) + du * data[5];
            
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            
                cp[3*i+3] = cps[3*i+3] + RELAX * dx;
                cp[3*i+4] = cps[3*i+4] + RELAX * dy;
                cp[3*i+5] = cps[3*i+5] + RELAX * dz;
                continue;
            }

            EG_spline1dEval_impl(header, rdata, knots[i+3], data);
            dx = xyz[3*i  ] - data[0];
            dy = xyz[3*i+1] - data[1];
            dz = xyz[3*i+2] - data[2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*i+3] = cps[3*i+3] + dx;
            cp[3*i+4] = cps[3*i+4] + dy;
            cp[3*i+5] = cps[3*i+5] + dz;
        }

        /* condition at end */
        EG_spline1dDeriv_impl(header, rdata, 2, knots[imax+2], data);
        du = knots[imax+2] - knots[imax+1];
        if (endc == 0) {
            /* natural end */
            dx = du * du * data[6];
            dy = du * du * data[7];
            dz = du * du * data[8];
        } else if (endc == 1) {
            /* FD slope */
            dx = xyz[3*(imax-2)  ] - xyz[3*(imax-1)  ] + du * data[3];
            dy = xyz[3*(imax-2)+1] - xyz[3*(imax-1)+1] + du * data[4];
            dz = xyz[3*(imax-2)+2] - xyz[3*(imax-1)+2] + du * data[5];
        } else {
            /* quadratic fit */
            u20    = knots[imax+2] - knots[imax];
            u21    = knots[imax+1] - knots[imax];
            rj[0]  = xyz[3*(imax-2)  ]*u20*u20 - xyz[3*(imax-1)  ]*u21*u21 -
                     xyz[3*(imax-3)  ]*du*du;
            rj[1]  = xyz[3*(imax-2)+1]*u20*u20 - xyz[3*(imax-1)+1]*u21*u21 -
                     xyz[3*(imax-3)+1]*du*du;
            rj[2]  = xyz[3*(imax-2)+2]*u20*u20 - xyz[3*(imax-1)+2]*u21*u21 -
                     xyz[3*(imax-3)+2]*du*du;
            rj[0] /= 2.0*u21*du;
            rj[1] /= 2.0*u21*du;
            rj[2] /= 2.0*u21*du;
            dx     = rj[0] - xyz[3*(imax-1)  ] + 0.5 * u20 * data[3];
            dy     = rj[1] - xyz[3*(imax-1)+1] + 0.5 * u20 * data[4];
            dz     = rj[2] - xyz[3*(imax-1)+2] + 0.5 * u20 * data[5];
        }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*imax  ] = cps[3*imax  ] + RELAX * dx;
        cp[3*imax+1] = cps[3*imax+1] + RELAX * dy;
        cp[3*imax+2] = cps[3*imax+2] + RELAX * dz;

        /* update the control points */
        for (i = 0; i < imax; i++) {
          cps[3*i+3] = cp[3*i+3];
          cps[3*i+4] = cp[3*i+4];
          cps[3*i+5] = cp[3*i+5];
        }

        /* convergence check */
        if (dxyzmax < tol) break;

    }
    if (dxyzmax >= tol)
        printf(" Warning: Not Converged (EG_spline1dFit)!\n");

    EG_free(mdata);
    return EGADS_SUCCESS;
}


extern "C"
int
EG_spline1dFit(int endx, int imaxx, const double *xyz, const double *kn,
               double tol, int *ivec, double **rdata)
{
  *rdata   = NULL;
  int imax = imaxx;
  if (imax < 0)                   imax = -imax;
  if (imax < 2)                   return EGADS_DEGEN;
  if ((endx <  0) || (endx >  2)) return EGADS_RANGERR;

  int icp      = imax + 2;
  int iknot    = imax + 6;
  double *rvec = (double *) EG_alloc((iknot+6*icp)*sizeof(double));
  if (rvec == NULL) return EGADS_MALLOC;

  int stat = EG_spline1dFit_impl(endx, imaxx, xyz, kn, tol, ivec, rvec);
  if (stat == EGADS_SUCCESS) {
    *rdata = (double *) EG_alloc((iknot+3*icp)*sizeof(double));
    if (*rdata == NULL) {
      EG_free(rvec);
      return EGADS_MALLOC;
    }
    for (int i = 0; i < iknot+3*icp; i++) (*rdata)[i] = rvec[i];
  }
  EG_free(rvec);

  return stat;
}


extern "C"
int
EG_spline1dFit_psens(int endx, int imaxx, const double *xyz, const double *kn,
                     double tol, int *ivec, double **rdata, double **rdata_xyz)
{
  SurrealD *xyzS = NULL, *rdataS = NULL;

  int imax = imaxx;
  if (imax == 0)                  return EGADS_DEGEN;
  if (imax < 0)                   imax = -imax;
  if (imax < 2)                   return EGADS_DEGEN;
  if ((endx <  0) || (endx >  2)) return EGADS_RANGERR;

  int icp   = imax + 2;
  int iknot = imax + 6;
  xyzS      = new SurrealD[3*imax];
  if (xyzS == NULL) return EGADS_MALLOC;
  rdataS    = new SurrealD[iknot+6*icp];
  if (rdataS == NULL) {
    delete [] xyzS;
    return EGADS_MALLOC;
  }

#ifndef __clang_analyzer__
  for (int i = 0; i < 3*imax; i++) {
    xyzS[i]          = SurrealD(xyz[i],3*imax);
    xyzS[i].deriv(i) = 1;
  }
#endif

  int stat = EG_spline1dFit_impl(endx, imaxx, xyzS, kn, tol, ivec, rdataS);
  delete [] xyzS;
  if (stat == EGADS_SUCCESS) {
    *rdata = (double *) EG_alloc((iknot+3*icp)*sizeof(double));
    if (*rdata == NULL) {
      delete [] rdataS;
      return EGADS_MALLOC;
    }

    *rdata_xyz = (double *) EG_alloc((iknot+3*icp)*3*imax*sizeof(double));
    if (*rdata_xyz == NULL) {
      EG_free(*rdata);
      delete [] rdataS;
      return EGADS_MALLOC;
    }

    for (int i = 0; i < iknot+3*icp; i++) {
      (*rdata)[i] = rdataS[i].value();
      if (rdataS[i].size() == 0) {
        for (int j = 0; j < 3*imax; j++)
          (*rdata_xyz)[i*3*imax+j] = 0;
      } else {
        for (int j = 0; j < 3*imax; j++)
          (*rdata_xyz)[i*3*imax+j] = rdataS[i].deriv(j);
      }
    }
  }
  delete [] rdataS;

  return stat;
}


template<int N>
int
EG_spline1dFit(int endx, int imaxx, const SurrealS<N> *xyz, const double *kn,
               double tol, int *ivec, SurrealS<N> **rdata)
{
  *rdata   = NULL;
  int imax = imaxx;
  if (imax < 0)                   imax = -imax;
  if (imax < 2)                   return EGADS_DEGEN;
  if ((endx <  0) || (endx >  2)) return EGADS_RANGERR;

  int icp      = imax + 2;
  int iknot    = imax + 6;
  SurrealS<N> *rvec = new SurrealS<N>[iknot+6*icp];
  if (rvec == NULL) return EGADS_MALLOC;

  int stat = EG_spline1dFit_impl(endx, imaxx, xyz, kn, tol, ivec, rvec);
  if (stat == EGADS_SUCCESS) {
    *rdata = new SurrealS<N>[iknot+3*icp];
    if (*rdata == NULL) {
      delete [] rvec;
      return EGADS_MALLOC;
    }
    for (int i = 0; i < iknot+3*icp; i++) (*rdata)[i] = rvec[i];
  }
  delete [] rvec;

  return stat;
}

// Create explicit instantiations of the function
template DllExport int EG_spline1dFit<1>(int, int, const SurrealS<1> *,
                                         const double *, double, int *,
                                         SurrealS<1> **);


extern "C"
int
EG_spline1d(egObject *context, int endx, int imaxx, const double *xyz,
            double tol, egObject **ecurv)
{
  int    stat, icp, iknot, imax, ivec[4];
  double *rvec;
  
  *ecurv = NULL;
  imax   = imaxx;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (imax < 0)                      imax = -imax;
  if (imax < 2)                      return EGADS_DEGEN;
  if ((endx <  0) || (endx >  2))    return EGADS_RANGERR;
  
  icp   = imax + 2;
  iknot = imax + 6;
  rvec  = (double *) EG_alloc((iknot+6*icp)*sizeof(double));
  if (rvec == NULL) return EGADS_MALLOC;
  
  stat = EG_spline1dFit_impl(endx, imaxx, xyz, NULL, tol, ivec, rvec);
  if (stat == EGADS_SUCCESS)
    stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL, ivec, rvec, ecurv);
  EG_free(rvec);
  
  return stat;
}


extern "C"
int
EG_spline1dPCrv(egObject *context, int endx, int imax, const double *xy,
                const double *knots, double tol, egObject **ecurv)
{
  int    i, stat, icp, iknot, ivec[4];
  double *rvec, *xyz, *cp;
  
  *ecurv = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (imax < 2)                      return EGADS_DEGEN;
  if ((endx <  0) || (endx >  2))    return EGADS_RANGERR;
  
  xyz = (double *) EG_alloc(3*imax*sizeof(double));
  if (xyz == NULL) return EGADS_MALLOC;
  for (i = 0; i < imax; i++) {
    xyz[3*i  ] = xy[2*i  ];
    xyz[3*i+1] = xy[2*i+1];
    xyz[3*i+2] = 0.0;
  }
  
  icp   = imax + 2;
  iknot = imax + 6;
  rvec  = (double *) EG_alloc((iknot+6*icp)*sizeof(double));
  if (rvec == NULL) {
    EG_free(xyz);
    return EGADS_MALLOC;
  }
  
  stat = EG_spline1dFit_impl(endx, imax, xyz, knots, tol, ivec, rvec);
  if (stat == EGADS_SUCCESS) {
    cp = &rvec[iknot];
    for (i = 1; i < imax; i++) {
      cp[2*i  ] = cp[3*i  ];
      cp[2*i+1] = cp[3*i+1];
    }
    stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, ivec, rvec, ecurv);
  }
  EG_free(rvec);
  EG_free(xyz);
  
  return stat;
}


/* ************************** Surface Functions **************************** */

template<class T, class T2>
static int
EG_spline2dEval_impl(int *ivec, T *data, T2 *uv, T *point)
{
  int i, j, l, degu, degv, nKu, nKv, nCPu, spanu, spanv;
  T   Nu[MAXDEG], Nv[MAXDEG], temp[4*MAXDEG], *Kv, *CP;

  point[0] = point[1] = point[2] = 0.0;
  degu     = ivec[1];
  nCPu     = ivec[2];
  nKu      = ivec[3];
  degv     = ivec[4];
  nKv      = ivec[6];
  Kv       = data + ivec[3];
  CP       = data + ivec[3] + ivec[6];
  if (ivec[0] != 0) {
    printf(" EG_spline2dEval: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degu >= MAXDEG) {
    printf(" EG_spline2dEval: degreeU %d >= %d!\n", degu, MAXDEG);
    return EGADS_CONSTERR;
  }
  if (degv >= MAXDEG) {
    printf(" EG_spline2dEval: degreeV %d >= %d!\n", degv, MAXDEG);
    return EGADS_CONSTERR;
  }

  spanu = FindSpan(nKu, degu, uv[0], data);
  BasisFuns(spanu, degu, uv[0], data, Nu);
  spanv = FindSpan(nKv, degv, uv[1], Kv);
  BasisFuns(spanv, degv, uv[1], Kv, Nv);

  for (l = 0; l <= degv; l++) {
    temp[3*l] = temp[3*l+1] = temp[3*l+2] = 0.0;
    for (j = 0; j <= degu; j++) {
      i            = spanu-degu+j + nCPu*(spanv-degv+l);
      temp[3*l  ] += Nu[j]*CP[3*i  ];
      temp[3*l+1] += Nu[j]*CP[3*i+1];
      temp[3*l+2] += Nu[j]*CP[3*i+2];
    }
  }

  for (l = 0; l <= degv; l++) {
    point[0] += Nv[l]*temp[3*l  ];
    point[1] += Nv[l]*temp[3*l+1];
    point[2] += Nv[l]*temp[3*l+2];
  }
  return EGADS_SUCCESS;
}


extern "C"
int
EG_spline2dEval(int *ivec, double *data, const double *uv, double *point)
{
  return EG_spline2dEval_impl(ivec, data, uv, point);
}


extern "C"
int
EG_spline2dEval_psens(int *ivec, const double *data, const double *data_dot,
                         const double *xyz_dot, const double *uv,
                         double *point, double *point_dot)
{
  SurrealS<1>* dataS = NULL;
  SurrealS<1>  pointS[3];

  point[0]     = point[1]     = point[2]     = 0.0;
  point_dot[0] = point_dot[1] = point_dot[2] = 0.0;
  int icp      =  ivec[2]*ivec[5];
  int iknot    =  ivec[3]+ivec[6];
  int imax     = (ivec[2]-2)*(ivec[5]-2);

  dataS = new SurrealS<1>[iknot+3*icp];
  if (dataS == NULL) return EGADS_MALLOC;

  for (int i = 0; i < iknot+3*icp; i++) {
    dataS[i] = data[i];
    double deriv = 0.0;
    for (int j = 0; j < 3*imax; j++) {
      deriv += data_dot[j + i*3*imax]*xyz_dot[j];
    }
    dataS[i].deriv() = deriv;
  }

  int stat = EG_spline2dEval_impl(ivec, dataS, uv, pointS);
  delete [] dataS;
  if (stat != EGADS_SUCCESS) return stat;

  for (int i = 0; i < 3; i++) {
    point[i]     = pointS[i].value();
    point_dot[i] = pointS[i].deriv();
  }
  return EGADS_SUCCESS;
}


template<int N>
int
EG_spline2dEval(int *ivec, SurrealS<N> *data, const double *uv,
                SurrealS<N> *point)
{
  return EG_spline2dEval_impl(ivec, data, uv, point);
}

// Create explicit instantiations of the function
template DllExport int EG_spline2dEval<1>(int *, SurrealS<1> *, const double *,
                                          SurrealS<1> *);


template<class T, class T2>
static int
EG_spline2dDeriv_impl(int *ivec, T *data, int der, const T2 *uv, T *deriv)
{
  int i, j, k, l, m, s, degu, degv, nKu, nKv, nCPu, spanu, spanv, du, dv;
  T   *Kv, *CP, *NderU[MAXDEG+1], *NderV[MAXDEG+1];
  T   Nu[MAXDEG+1][MAXDEG+1], Nv[MAXDEG+1][MAXDEG+1], temp[3*MAXDEG];

  degu = ivec[1];
  nCPu = ivec[2];
  nKu  = ivec[3];
  degv = ivec[4];
  nKv  = ivec[6];
  Kv   = data + ivec[3];
  CP   = data + ivec[3] + ivec[6];
  du   = MIN(der, degu);
  dv   = MIN(der, degv);
  for (m = l = 0; l <= der; l++)
    for (k = 0; k <= der-l; k++, m++)
      deriv[3*m  ] = deriv[3*m+1] = deriv[3*m+2] = 0.0;
  if (ivec[0] != 0) {
    printf(" EG_spline2dDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degu >= MAXDEG) {
    printf(" EG_spline2dDeriv: degreeU %d >= %d!\n", degu, MAXDEG);
    return EGADS_CONSTERR;
  }
  if (degv >= MAXDEG) {
    printf(" EG_spline2dDeriv: degreeV %d >= %d!\n", degv, MAXDEG);
    return EGADS_CONSTERR;
  }
  for (i = 0; i <= degu; i++) NderU[i] = &Nu[i][0];
  for (i = 0; i <= degv; i++) NderV[i] = &Nv[i][0];

  spanu = FindSpan(nKu, degu, uv[0], data);
  DersBasisFuns(spanu, degu, uv[0], data, du, NderU);
  spanv = FindSpan(nKv, degv, uv[1], Kv);
  DersBasisFuns(spanv, degv, uv[1], Kv, dv, NderV);

  for (m = l = 0; l <= dv; l++)
    for (k = 0; k <= der-l; k++, m++) {
      if (k > du) continue;
      for (s = 0; s <= degv; s++) {
        temp[3*s] = temp[3*s+1] = temp[3*s+2] = 0.0;
        for (j = 0; j <= degu; j++) {
          i            = spanu-degu+j + nCPu*(spanv-degv+s);
          temp[3*s  ] += Nu[k][j]*CP[3*i  ];
          temp[3*s+1] += Nu[k][j]*CP[3*i+1];
          temp[3*s+2] += Nu[k][j]*CP[3*i+2];
        }
      }

      for (s = 0; s <= degv;  s++) {
        deriv[3*m  ] += Nv[l][s]*temp[3*s  ];
        deriv[3*m+1] += Nv[l][s]*temp[3*s+1];
        deriv[3*m+2] += Nv[l][s]*temp[3*s+2];
      }
    }
  return EGADS_SUCCESS;
}


extern "C"
int
EG_spline2dDeriv(int *ivec, double *data, int der, const double *uv,
                 double *deriv)
{
  return EG_spline2dDeriv_impl(ivec, data, der, uv, deriv);
}


template<int N>
int
EG_spline2dDeriv(int *ivec, SurrealS<N> *data, int der,
                 const double *uv, SurrealS<N> *deriv)
{
  return EG_spline2dDeriv_impl(ivec, data, der, uv, deriv);
}

// Create explicit instantiations of the function
template DllExport int EG_spline2dDeriv<1>(int *, SurrealS<1> *, int,
                                           const double *, SurrealS<1> *);



template<class T, class T2>
static T
EG_cubicBasisFun(int nKnots, T *knots, int i, T2 u)
{
  int j, k;
  T   N[4], saved, Uleft, Uright, temp;
  
  if ((i < 0) || (i > nKnots-5)) return -1.0;
  
  if (((i == 0)        && (u == knots[0])) ||
      ((i == nKnots-5) && (u == knots[nKnots-1]))) return 1.0;
  
  if ((u < knots[i]) || (u >= knots[i+4])) return 0.0;
  
  for (j = 0; j <= 3; j++)
    if ((u >= knots[i+j]) && (u < knots[i+j+1])) {
      N[j] = 1.0;
    } else {
      N[j] = 0.0;
    }
  
  for (k = 1; k <= 3; k++) {
    if (N[0] == 0.0) {
      saved = 0.0;
    } else {
      saved = ((u-knots[i])*N[0])/(knots[i+k]-knots[i]);
    }
    for (j = 0; j < 3-k+1; j++) {
      Uleft  = knots[i+j  +1];
      Uright = knots[i+j+k+1];
      if (N[j+1] == 0.0) {
        N[j]  = saved;
        saved = 0.0;
      } else {
        temp  = N[j+1]/(Uright-Uleft);
        N[j]  = saved + (Uright-u)*temp;
        saved = (u-Uleft)*temp;
      }
    }
  }
  
  return N[0];
}


template<class T, class T2>
static void
EG_NzeroBSp(int nknots, T *knots, T2 t, int *nbasis, T *basis)
{
  int jj, k;
  
  jj = 3;
  while((knots[jj] < t) && (jj <= nknots-4)) jj++;
  
  for (k = 0; k <= 3; k++)
    basis[k] = EG_cubicBasisFun(nknots, knots, jj+k-4, t);
  
  *nbasis  = 4;
  if (basis[0] != 0.0) return;
  
  *nbasis  = 3;
  basis[0] = basis[1];
  basis[1] = basis[2];
  basis[2] = basis[3];
  basis[3] = 0.0;
}


template<class T>
static T
EG_getEllRad(const T *drnd, T *nell)
{
  T csth, snth;
  
  csth = drnd[0]*(drnd[1]*nell[0] + drnd[2]*nell[1] + drnd[3]*nell[2]);
  snth = drnd[4]*(drnd[5]*nell[0] + drnd[6]*nell[1] + drnd[7]*nell[2]);
  return sqrt(csth*csth + snth*snth);
}


/*
 ************************************************************************
 *                                                                      *
 *   EG_spline2dAppr - Approximate 2d cubic spline from input data      *
 *                                                                      *
 ************************************************************************
 */

template<class T>
static int
EG_spline2dAppr(int endc, int imax, int jmax, const T *xyz,
                /*@null@*/ const T   *uknot, /*@null@*/ const T *vknot,
                /*@null@*/ const int *vdata,
                /*@null@*/ const T   *wesT,  /*@null@*/ const T *easT,
                /*@null@*/ const T   *south, /*@null@*/       T *snor,
                /*@null@*/ const T   *north, /*@null@*/       T *nnor,
                double tol, int *header, T **rdata)
{
    int i, j, iknot, jknot, icp, jcp, iter, tanOK;
    int endi, endj, jj, kk, ms, mn;
    T   ns[3], nn[3], rs[3][3], rn[3][3], thet, q0, q1, q2, q3, x2[3];
    T   r, tt, du, dv, dx, dy, dz, dist, mmu, con, con2, dxyzmax, normnell;
    T   dD, eE, F, G, rj[3], u21, u20, norm[3], nell[3], x0[3], x1[3], t[3][3];
    T   basis[4], uv[2], eval[18], *rvec, *knotu, *knotv, *cp, *cpsav;
    double box[6], rsize;

#ifdef DEBUG
    printf(" In EG_splin2dAppr: endc = %d, imax = %d, jmax = %d  tol = %le\n",
           endc, imax, jmax, tol);
    printf("                    south ");
    if (south == NULL) {
      printf("= NULL");
    } else {
      printf("!= NULL");
    }
    printf("   snor ");
    if (snor == NULL) {
      printf("= NULL");
    } else {
      printf("!= NULL");
    }
    printf("   north ");
    if (north == NULL) {
      printf("= NULL");
    } else {
      printf("!= NULL");
    }
    printf("   nnor ");
    if (nnor == NULL) {
      printf("= NULL");
    } else {
      printf("!= NULL");
    }
    printf("\n");
    if (uknot != NULL)
        for (i = 0; i < imax; i+=7) {
            if (i == 0) {
                printf("  uKnot =");
            } else {
                printf("         ");
            }
            for (j = i; j < i+7; j++) {
              if (j >= imax) continue;
              printf(" %lf", value(uknot[j]));
            }
            printf("\n");
        }
    if (vknot != NULL)
        for (i = 0; i < jmax; i+=7) {
            if (i == 0) {
                printf("  vKnot =");
            } else {
                printf("         ");
            }
            for (j = i; j < i+7; j++) {
              if (j >= jmax) continue;
              printf(" %lf", value(vknot[j]));
            }
            printf("\n");
        }
    if (vdata != NULL)
        for (i = 0; i < jmax; i+=15) {
            if (i == 0) {
                printf("  vData =");
            } else {
                printf("         ");
            }
            for (j = i; j < i+15; j++) {
              if (j >= jmax) continue;
              printf(" %d", vdata[j]);
            }
            printf("\n");
        }
    for (j = 0; j < jmax; j++)
        for (i = 0; i < imax; i++)
            printf(" xyz[%d,%d] = %lf %lf %lf\n", j, i,
                   value(xyz[3*(j*imax+i)  ]), value(xyz[3*(j*imax+i)+1]),
                   value(xyz[3*(j*imax+i)+2]));
    if (wesT != NULL)
        for (j = 0; j < jmax; j++)
            printf(" wesT[%d] = %lf %lf %lf\n", j, value(wesT[3*j  ]),
                   value(wesT[3*j+1]), value(wesT[3*j+2]));
    if (easT != NULL)
        for (j = 0; j < jmax; j++)
            printf(" easT[%d] = %lf %lf %lf\n", j, value(easT[3*j  ]),
                   value(easT[3*j+1]), value(easT[3*j+2]));
#endif

    ms       = 1;
    rs[0][0] = rs[1][1] = rs[2][2] = 1.0;
    rs[0][1] = rs[0][2] = rs[1][0] = rs[1][2] = rs[2][0] = rs[2][1] = 0.0;
  
    mn       = 1;
    rn[0][0] = rn[1][1] = rn[2][2] = 1.0;
    rn[0][1] = rn[0][2] = rn[1][0] = rn[1][2] = rn[2][0] = rn[2][1] = 0.0;

    *rdata = NULL;
    endi   = endj = endc;
    if (imax == 2) endi = 1;
    if (jmax == 2) endj = 1;

    icp   = imax + 2;
    iknot = imax + 6;
    jcp   = jmax + 2;
    jknot = jmax + 6;
    rvec  = (T *) EG_alloc((iknot+jknot+3*icp*jcp)*sizeof(T));
    if (rvec == NULL) return EGADS_MALLOC;
    knotu =  rvec;
    knotv = &rvec[iknot      ];
    cpsav = &rvec[iknot+jknot];
    cp    = (T *) EG_alloc((3*icp*jcp)*sizeof(T));
    if (cp == NULL) {
      EG_free(rvec);
      return EGADS_MALLOC;
    }
  
    /* get relative size */
    box[0] = box[3] = value(xyz[0]);
    box[1] = box[4] = value(xyz[1]);
    box[2] = box[5] = value(xyz[2]);
    for (j = 0; j < jmax; j++)
        for (i = 0; i < imax; i++) {
            if (value(xyz[3*((i)+(j)*imax)  ]) < box[0])
                box[0] = value(xyz[3*((i)+(j)*imax)  ]);
            if (value(xyz[3*((i)+(j)*imax)  ]) > box[3])
                box[3] = value(xyz[3*((i)+(j)*imax)  ]);
            if (value(xyz[3*((i)+(j)*imax)+1]) < box[1])
                box[1] = value(xyz[3*((i)+(j)*imax)+1]);
            if (value(xyz[3*((i)+(j)*imax)+1]) > box[4])
                box[4] = value(xyz[3*((i)+(j)*imax)+1]);
            if (value(xyz[3*((i)+(j)*imax)+2]) < box[2])
                box[2] = value(xyz[3*((i)+(j)*imax)+2]);
            if (value(xyz[3*((i)+(j)*imax)+2]) > box[5])
                box[5] = value(xyz[3*((i)+(j)*imax)+2]);
        }
    rsize = sqrt((box[3]-box[0])*(box[3]-box[0]) +
                 (box[4]-box[1])*(box[4]-box[1]) +
                 (box[5]-box[2])*(box[5]-box[2]));

    /* create spline surface */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;
    header[4] = 3;
    header[5] = jcp;
    header[6] = jknot;

    /* knots in i-direction */
    kk          = 0;
    knotu[kk++] = 0.0;
    knotu[kk++] = 0.0;
    knotu[kk++] = 0.0;
    if (uknot == NULL) {
        /* arc-length spaced */
        for (i = 0; i < imax; i++) knotu[kk+i] = 0.0;
        dz = jmax;
        for (j = 0; j < jmax; j++) {
            dy = 0.0;
            for (i = 1; i < imax; i++) {
                dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
                           (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
                           (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
                           (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
                           (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
                           (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
            }
            if (dy < 1.e-13*rsize) {
                dz -= 1.0;
                continue;
            }
            dx = 0.0;
            for (i = 1; i < imax; i++) {
              dx += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
                         (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
                         (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
                         (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
                         (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
                         (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]))/dy;
              knotu[kk+i] += dx;
            }
       }
        for (iter = i = 1; i < imax; i++)
            if (knotu[kk+i] <= knotu[kk+i-1]) iter++;
        if (iter == 1) {
          for (i = 0; i < imax; i++) knotu[kk++] /= dz;
        } else {
            /* equally spaced */
            for (i = 0; i < imax; i++) {
                dx          = i;
                knotu[kk++] = dx/(imax-1);
            }
        }
    } else {
        /* knots in u set by input */
        for (i = 0; i < imax; i++) knotu[kk++] = uknot[i];
    }
    knotu[kk++] = 1.0;
    knotu[kk++] = 1.0;
    knotu[kk++] = 1.0;

    /* knots in j-direction -- allow for multiplicity */
    kk          = 0;
    knotv[kk++] = 0.0;
    knotv[kk++] = 0.0;
    knotv[kk++] = 0.0;
    if (vknot == NULL) {
        /* arc-length spaced */
        for (j = 0; j < jmax; j++) knotv[kk+j] = 0.0;
        dz = imax;
        for (i = 0; i < imax; i++) {
            dy = 0.0;
            for (j = 1; j < jmax; j++) {
                dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
                           (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
                           (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
                           (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
                           (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
                           (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
            }
            if (dy < 1.e-13*rsize) {
                dz -= 1.0;
                continue;
            }
            dx = 0.0;
            for (j = 1; j < jmax; j++) {
                dx += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
                           (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
                           (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
                           (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
                           (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
                           (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]))/dy;
                knotv[kk+j] += dx;
            }
        }
        for (j = 0; j < jmax; j++) knotv[kk++] /= dz;
    } else {
        /* knots in v set by input */
        for (j = 0; j < jmax; j++) knotv[kk++] = vknot[j];
    }
    knotv[kk++] = 1.0;
    knotv[kk++] = 1.0;
    knotv[kk++] = 1.0;
  
    /* map of IDs/indices for imax=8 and jmax=5 (used in comments below)

                    0  1  2  3  4  5  6  7  8  9 <-CPs
                                                    v
               4    nw O  n  n  n  n  n  n  P ne    6
               .    J  K  L  L  L  L  L  L  M  N    5
               3    w  H  *  *  *  *  *  *  I  e    4
               2    w  H  *  *  *  *  *  *  I  e    3
               1    w  H  *  *  *  *  *  *  I  e    2
               .    C  D  E  E  E  E  E  E  F  G    1
               0    sw A  s  s  s  s  s  s  B se    0
               ^
              xyz-> 0  .  1  2  3  4  5  6  .  7
     
       2 additional CPs per row/column for end condition control */
  
    kk = 0;

    /* southwest control point */
    cp[kk++] = xyz[3*((0)+(0)*imax)  ];
    cp[kk++] = xyz[3*((0)+(0)*imax)+1];
    cp[kk++] = xyz[3*((0)+(0)*imax)+2];

    /* point A */
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)  ] + xyz[3*((1)+(0)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)+1] + xyz[3*((1)+(0)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)+2] + xyz[3*((1)+(0)*imax)+2]) / 4;

    /* south control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = xyz[3*((i)+(0)*imax)  ];
        cp[kk++] = xyz[3*((i)+(0)*imax)+1];
        cp[kk++] = xyz[3*((i)+(0)*imax)+2];
    }

    /* point B */
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)  ] +
                    xyz[3*((imax-2)+(0)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)+1] +
                    xyz[3*((imax-2)+(0)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)+2] +
                    xyz[3*((imax-2)+(0)*imax)+2]) / 4;

    /* southeast control point */
    cp[kk++] = xyz[3*((imax-1)+(0)*imax)  ];
    cp[kk++] = xyz[3*((imax-1)+(0)*imax)+1];
    cp[kk++] = xyz[3*((imax-1)+(0)*imax)+2];

    /* point C */
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)  ] + xyz[3*((0)+(1)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)+1] + xyz[3*((0)+(1)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)+2] + xyz[3*((0)+(1)*imax)+2]) / 4;

    /* point D */
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)  ] + xyz[3*((1)+(1)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)+1] + xyz[3*((1)+(1)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(0)*imax)+2] + xyz[3*((1)+(1)*imax)+2]) / 4;

    /* points E */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = (3 * xyz[3*((i)+(0)*imax)  ] + xyz[3*((i)+(1)*imax)  ]) / 4;
        cp[kk++] = (3 * xyz[3*((i)+(0)*imax)+1] + xyz[3*((i)+(1)*imax)+1]) / 4;
        cp[kk++] = (3 * xyz[3*((i)+(0)*imax)+2] + xyz[3*((i)+(1)*imax)+2]) / 4;
    }

    /* point F */
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)  ] +
                    xyz[3*((imax-2)+(1)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)+1] +
                    xyz[3*((imax-2)+(1)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)+2] +
                    xyz[3*((imax-2)+(1)*imax)+2]) / 4;

    /* point G */
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)  ] +
                    xyz[3*((imax-1)+(1)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)+1] +
                    xyz[3*((imax-1)+(1)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(0)*imax)+2] +
                    xyz[3*((imax-1)+(1)*imax)+2]) / 4;

    /* loop through interior j lines */
    for (j = 1; j < jmax-1; j++) {

        /* west control point */
        cp[kk++] = xyz[3*((0)+(j)*imax)  ];
        cp[kk++] = xyz[3*((0)+(j)*imax)+1];
        cp[kk++] = xyz[3*((0)+(j)*imax)+2];

        /* point H */
        cp[kk++] = (3 * xyz[3*((0)+(j)*imax)  ] + xyz[3*((1)+(j)*imax)  ]) / 4;
        cp[kk++] = (3 * xyz[3*((0)+(j)*imax)+1] + xyz[3*((1)+(j)*imax)+1]) / 4;
        cp[kk++] = (3 * xyz[3*((0)+(j)*imax)+2] + xyz[3*((1)+(j)*imax)+2]) / 4;

        /* interior points */
        for (i = 1; i < imax-1; i++) {
            cp[kk++] = xyz[3*((i)+(j)*imax)  ];
            cp[kk++] = xyz[3*((i)+(j)*imax)+1];
            cp[kk++] = xyz[3*((i)+(j)*imax)+2];
        }

        /* point I */
        cp[kk++] = (3 * xyz[3*((imax-1)+(j)*imax)  ] +
                        xyz[3*((imax-2)+(j)*imax)  ]) / 4;
        cp[kk++] = (3 * xyz[3*((imax-1)+(j)*imax)+1] +
                        xyz[3*((imax-2)+(j)*imax)+1]) / 4;
        cp[kk++] = (3 * xyz[3*((imax-1)+(j)*imax)+2] +
                        xyz[3*((imax-2)+(j)*imax)+2]) / 4;

        /* east control point */
        cp[kk++] = xyz[3*((imax-1)+(j)*imax)  ];
        cp[kk++] = xyz[3*((imax-1)+(j)*imax)+1];
        cp[kk++] = xyz[3*((imax-1)+(j)*imax)+2];
    }

    /* point J */
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)  ] +
                    xyz[3*((0)+(jmax-2)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)+1] +
                    xyz[3*((0)+(jmax-2)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)+2] +
                    xyz[3*((0)+(jmax-2)*imax)+2]) / 4;

    /* point K */
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)  ] +
                    xyz[3*((1)+(jmax-2)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)+1] +
                    xyz[3*((1)+(jmax-2)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)+2] +
                    xyz[3*((1)+(jmax-2)*imax)+2]) / 4;

    /* points L */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = (3 * xyz[3*((i)+(jmax-1)*imax)  ] +
                        xyz[3*((i)+(jmax-2)*imax)  ]) / 4;
        cp[kk++] = (3 * xyz[3*((i)+(jmax-1)*imax)+1] +
                        xyz[3*((i)+(jmax-2)*imax)+1]) / 4;
        cp[kk++] = (3 * xyz[3*((i)+(jmax-1)*imax)+2] +
                        xyz[3*((i)+(jmax-2)*imax)+2]) / 4;
    }

    /* point M */
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)  ] +
                    xyz[3*((imax-2)+(jmax-2)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)+1] +
                    xyz[3*((imax-2)+(jmax-2)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)+2] +
                    xyz[3*((imax-2)+(jmax-2)*imax)+2]) / 4;

    /* point N */
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)  ] +
                    xyz[3*((imax-1)+(jmax-2)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)+1] +
                    xyz[3*((imax-1)+(jmax-2)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)+2] +
                    xyz[3*((imax-1)+(jmax-2)*imax)+2]) / 4;

    /* northwest control point */
    cp[kk++] = xyz[3*((0)+(jmax-1)*imax)  ];
    cp[kk++] = xyz[3*((0)+(jmax-1)*imax)+1];
    cp[kk++] = xyz[3*((0)+(jmax-1)*imax)+2];

    /* point O */
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)  ] +
                    xyz[3*((1)+(jmax-1)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)+1] +
                    xyz[3*((1)+(jmax-1)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((0)+(jmax-1)*imax)+2] +
                    xyz[3*((1)+(jmax-1)*imax)+2]) / 4;

    /* north control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = xyz[3*((i)+(jmax-1)*imax)  ];
        cp[kk++] = xyz[3*((i)+(jmax-1)*imax)+1];
        cp[kk++] = xyz[3*((i)+(jmax-1)*imax)+2];
    }

    /* point P */
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)  ] +
                    xyz[3*((imax-2)+(jmax-1)*imax)  ]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)+1] +
                    xyz[3*((imax-2)+(jmax-1)*imax)+1]) / 4;
    cp[kk++] = (3 * xyz[3*((imax-1)+(jmax-1)*imax)+2] +
                    xyz[3*((imax-2)+(jmax-1)*imax)+2]) / 4;

    /* northeast control point */
    cp[kk++] = xyz[3*((imax-1)+(jmax-1)*imax)  ];
    cp[kk++] = xyz[3*((imax-1)+(jmax-1)*imax)+1];
    cp[kk++] = xyz[3*((imax-1)+(jmax-1)*imax)+2];
  
    /* fix things up for degenerate treatments */
    if (south != NULL) {
        x0[0] = south[1];
        x0[1] = south[2];
        x0[2] = south[3];
        x1[0] = south[5];
        x1[1] = south[6];
        x1[2] = south[7];
        CROSS(ns, x0, x1);
        dist  = sqrt(DOT(ns,ns));
        x0[0] = xyz[3*(2*imax+0)  ] - xyz[3*(0*imax+0)  ];
        x0[1] = xyz[3*(2*imax+0)+1] - xyz[3*(0*imax+0)+1];
        x0[2] = xyz[3*(2*imax+0)+2] - xyz[3*(0*imax+0)+2];
        if (DOT(x0,ns) > 0.0) {
          ns[0] /=  dist;
          ns[1] /=  dist;
          ns[2] /=  dist;
        } else {
          ns[0] /= -dist;
          ns[1] /= -dist;
          ns[2] /= -dist;
        }
        ms   = 0;
        dist = ns[0];
        if (fabs(ns[1]) > fabs(dist)) {
          ms   = 1;
          dist = ns[1];
        }
        if (fabs(ns[2]) > fabs(dist)) ms = 2;
        thet = acos(ns[ms]);
#ifdef DEBUG
        printf("\n Nose plane = %lf %lf %lf   m = %d, thet = %lf\n\n",
               value(ns[0]), value(ns[1]), value(ns[2]), ms, value(thet));
#endif
        if (thet != 0.0) {
          x1[0]  = x1[1] = x1[2] = 0.0;
          x1[ms] = 1.0;
          CROSS(x0, x1, ns);
          dist = sqrt(DOT(x0, x0));
          if (dist != 0.0) {
            x0[0] /= dist;
            x0[1] /= dist;
            x0[2] /= dist;
#ifdef DEBUG
            printf(" Unit axis of rotation = %lf %lf %lf    %lf\n",
                   value(x0[0]), value(x0[1]), value(x0[2]), value(dist));
#endif
            dist   = sin(thet/2.0);
            q0     = cos(thet/2.0);
            q1     = x0[0]*dist;
            q2     = x0[1]*dist;
            q3     = x0[2]*dist;
#ifdef DEBUG
            printf(" Quaterion = %lf %lf %lf %lf\n",
                   value(q0), value(q1), value(q2), value(q3));
#endif
            rs[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
            rs[0][1] = 2.0*(q1*q2 - q0*q3);
            rs[0][2] = 2.0*(q1*q3 + q0*q2);
            rs[1][0] = 2.0*(q1*q2 + q0*q3);
            rs[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
            rs[1][2] = 2.0*(q2*q3 - q0*q1);
            rs[2][0] = 2.0*(q1*q3 - q0*q2);
            rs[2][1] = 2.0*(q2*q3 + q0*q1);
            rs[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
          }
        }
#ifdef DEBUG
        printf(" Rotation Matrix = %lf %lf %lf\n",
               value(rs[0][0]), value(rs[0][1]), value(rs[0][2]));
        printf("                   %lf %lf %lf\n",
               value(rs[1][0]), value(rs[1][1]), value(rs[1][2]));
        printf("                   %lf %lf %lf\n",
               value(rs[2][0]), value(rs[2][1]), value(rs[2][2]));
        printf("\n");
#endif
        mmu = 2.0*(knotv[4] - knotv[3])/(3.0*((knotv[5] - knotv[3])));
        for (i = -1; i < imax+1; i++) {
            x2[0]    = cp[3*(2*(imax+2)+i+1)  ] - cp[3*(0*(imax+2)+i+1)  ];
            x2[1]    = cp[3*(2*(imax+2)+i+1)+1] - cp[3*(0*(imax+2)+i+1)+1];
            x2[2]    = cp[3*(2*(imax+2)+i+1)+2] - cp[3*(0*(imax+2)+i+1)+2];
            x2[ms]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rs[0], x2);
            nell[1]  = DOT(rs[1], x2);
            nell[2]  = DOT(rs[2], x2);
            r        = EG_getEllRad(south, nell);
            x1[0]    = cp[3*(2*(imax+2)+i+1)  ] - cp[3*(0*(imax+2)+i+1)  ];
            x1[1]    = cp[3*(2*(imax+2)+i+1)+1] - cp[3*(0*(imax+2)+i+1)+1];
            x1[2]    = cp[3*(2*(imax+2)+i+1)+2] - cp[3*(0*(imax+2)+i+1)+2];
            CROSS(x2, nell, x1);
            tt       = sqrt(r*mmu*sqrt(DOT(x2, x2)));
#ifdef DEBUG
            printf(" south init %d: %lf %lf %lf\n", i+1,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]));
#endif
            cp[3*((i+1)+(1)*(imax+2))  ] = cp[3*((i+1)+(0)*(imax+2))  ] +
                                           tt*nell[0];
            cp[3*((i+1)+(1)*(imax+2))+1] = cp[3*((i+1)+(0)*(imax+2))+1] +
                                           tt*nell[1];
            cp[3*((i+1)+(1)*(imax+2))+2] = cp[3*((i+1)+(0)*(imax+2))+2] +
                                           tt*nell[2];
        }
    }
  
    if (north != NULL) {
        x0[0] = north[1];
        x0[1] = north[2];
        x0[2] = north[3];
        x1[0] = north[5];
        x1[1] = north[6];
        x1[2] = north[7];
        CROSS(nn, x0, x1);
        dist  = sqrt(DOT(nn,nn));
        x0[0] = xyz[3*((jmax-3)*imax+0)  ] - xyz[3*((jmax-1)*imax+0)  ];
        x0[1] = xyz[3*((jmax-3)*imax+0)+1] - xyz[3*((jmax-1)*imax+0)+1];
        x0[2] = xyz[3*((jmax-3)*imax+0)+2] - xyz[3*((jmax-1)*imax+0)+2];
        if (DOT(x0,nn) > 0.0) {
          nn[0] /=  dist;
          nn[1] /=  dist;
          nn[2] /=  dist;
        } else {
          nn[0] /= -dist;
          nn[1] /= -dist;
          nn[2] /= -dist;
        }
        mn   = 0;
        dist = nn[0];
        if (fabs(nn[1]) > fabs(dist)) {
          mn   = 1;
          dist = nn[1];
        }
        if (fabs(nn[2]) > fabs(dist)) mn = 2;
        thet = acos(nn[mn]);
#ifdef DEBUG
        printf("\n Nose plane = %lf %lf %lf   m = %d, thet = %lf\n\n",
               value(nn[0]), value(nn[1]), value(nn[2]), mn, value(thet));
#endif
        if (thet != 0.0) {
          x1[0]  = x1[1] = x1[2] = 0.0;
          x1[mn] = 1.0;
          CROSS(x0, x1, nn);
          dist = sqrt(DOT(x0, x0));
          if (dist != 0.0) {
            x0[0] /= dist;
            x0[1] /= dist;
            x0[2] /= dist;
#ifdef DEBUG
            printf(" Unit axis of rotation = %lf %lf %lf    %lf\n",
                   value(x0[0]), value(x0[1]), value(x0[2]), value(dist));
#endif
            dist   = sin(thet/2.0);
            q0     = cos(thet/2.0);
            q1     = x0[0]*dist;
            q2     = x0[1]*dist;
            q3     = x0[2]*dist;
#ifdef DEBUG
            printf(" Quaterion = %lf %lf %lf %lf\n",
                   value(q0), value(q1), value(q2), value(q3));
#endif
            rn[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
            rn[0][1] = 2.0*(q1*q2 - q0*q3);
            rn[0][2] = 2.0*(q1*q3 + q0*q2);
            rn[1][0] = 2.0*(q1*q2 + q0*q3);
            rn[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
            rn[1][2] = 2.0*(q2*q3 - q0*q1);
            rn[2][0] = 2.0*(q1*q3 - q0*q2);
            rn[2][1] = 2.0*(q2*q3 + q0*q1);
            rn[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
          }
        }
#ifdef DEBUG
        printf(" Rotation Matrix = %lf %lf %lf\n",
               value(rn[0][0]), value(rn[0][1]), value(rn[0][2]));
        printf("                   %lf %lf %lf\n",
               value(rn[1][0]), value(rn[1][1]), value(rn[1][2]));
        printf("                   %lf %lf %lf\n",
               value(rn[2][0]), value(rn[2][1]), value(rn[2][2]));
        printf("\n");
#endif
        mmu =  2.0*( knotv[jmax+1] - knotv[jmax+2])/
              (3.0*((knotv[jmax]   - knotv[jmax+2])));
        for (i = -1; i < imax+1; i++) {
            x2[0]    = cp[3*((jmax-1)*(imax+2)+i+1)  ] -
                       cp[3*((jmax+1)*(imax+2)+i+1)  ];
            x2[1]    = cp[3*((jmax-1)*(imax+2)+i+1)+1] -
                       cp[3*((jmax+1)*(imax+2)+i+1)+1];
            x2[2]    = cp[3*((jmax-1)*(imax+2)+i+1)+2] -
                       cp[3*((jmax+1)*(imax+2)+i+1)+2];
            x2[mn]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rn[0], x2);
            nell[1]  = DOT(rn[1], x2);
            nell[2]  = DOT(rn[2], x2);
            r        = EG_getEllRad(north, nell);
            x1[0]    = cp[3*((jmax-1)*(imax+2)+i+1)  ] -
                       cp[3*((jmax+1)*(imax+2)+i+1)  ];
            x1[1]    = cp[3*((jmax-1)*(imax+2)+i+1)+1] -
                       cp[3*((jmax+1)*(imax+2)+i+1)+1];
            x1[2]    = cp[3*((jmax-1)*(imax+2)+i+1)+2] -
                       cp[3*((jmax+1)*(imax+2)+i+1)+2];
            CROSS(x2, nell, x1);
            tt       = sqrt(r*mmu*sqrt(DOT(x2, x2)));
#ifdef DEBUG
            printf(" north init %d: %lf %lf %lf\n", i+1,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]));
#endif
            cp[3*((i+1)+(jmax)*(imax+2))  ] = cp[3*((i+1)+(jmax+1)*(imax+2))  ] +
                                              tt*nell[0];
            cp[3*((i+1)+(jmax)*(imax+2))+1] = cp[3*((i+1)+(jmax+1)*(imax+2))+1] +
                                              tt*nell[1];
            cp[3*((i+1)+(jmax)*(imax+2))+2] = cp[3*((i+1)+(jmax+1)*(imax+2))+2] +
                                              tt*nell[2];
        }
    }

    /* iterate to have knot evaluations match data points */
    for (iter = 0; iter < NITER; iter++) {

        dxyzmax = 0.0;
        for (i = 0; i < 3*icp*jcp; i++) cpsav[i] = cp[i];

        /* match interior spline points */
        for (j = 1; j < jmax-1; j++) {
            uv[1] = knotv[j+3];
            if (vdata != NULL) {
                if (vdata[j] == +2) {
                    /* multiplicity 2 with flat spot on right
                       note: j is the second repeated data point */
                    dv = knotv[j+4] - knotv[j+3];
                    for (i = 1; i < imax-1; i++) {
                        uv[0] = knotu[i+3];
                        EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
                        dx = xyz[3*((i)+(j+1)*imax)  ] -
                             xyz[3*((i)+(j  )*imax)  ] - dv*eval[ 9];
                        dy = xyz[3*((i)+(j+1)*imax)+1] -
                             xyz[3*((i)+(j  )*imax)+1] - dv*eval[10];
                        dz = xyz[3*((i)+(j+1)*imax)+2] -
                             xyz[3*((i)+(j  )*imax)+2] - dv*eval[11];
                        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                        cp[3*((i+1)+(j+1)*(imax+2))  ] += RELAX * dx;
                        cp[3*((i+1)+(j+1)*(imax+2))+1] += RELAX * dy;
                        cp[3*((i+1)+(j+1)*(imax+2))+2] += RELAX * dz;
                    }
                    continue;
                }
                if (vdata[j] == -2) {
                    /* multiplicity 2 with flat spot on left
                       note: j is the first repeated data point */
                    dv = knotv[j+3] - knotv[j+2];
                    for (i = 1; i < imax-1; i++) {
                        uv[0] = knotu[i+3];
                        EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
                        dx = xyz[3*((i)+(j-1)*imax)  ] -
                             xyz[3*((i)+(j  )*imax)  ] + dv*eval[ 9];
                        dy = xyz[3*((i)+(j-1)*imax)+1] -
                             xyz[3*((i)+(j  )*imax)+1] + dv*eval[10];
                        dz = xyz[3*((i)+(j-1)*imax)+2] -
                             xyz[3*((i)+(j  )*imax)+2] + dv*eval[11];
                        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                        cp[3*((i+1)+(j+1)*(imax+2))  ] += RELAX * dx;
                        cp[3*((i+1)+(j+1)*(imax+2))+1] += RELAX * dy;
                        cp[3*((i+1)+(j+1)*(imax+2))+2] += RELAX * dz;
                    }
                    continue;
                }
            }

            for (i = 1; i < imax-1; i++) {
                uv[0] = knotu[i+3];
                EG_spline2dEval_impl(header, rvec, uv, eval);
                dx = xyz[3*((i)+(j)*imax)  ] - eval[0];
                dy = xyz[3*((i)+(j)*imax)+1] - eval[1];
                dz = xyz[3*((i)+(j)*imax)+2] - eval[2];
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(j+1)*(imax+2))  ] += dx;
                cp[3*((i+1)+(j+1)*(imax+2))+1] += dy;
                cp[3*((i+1)+(j+1)*(imax+2))+2] += dz;
            }
        }

        /* point A */
        uv[0] = knotu[3];
        uv[1] = knotv[3];
        EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
        du = knotu[4] - knotu[3];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * eval[6];
            dy = du * du * eval[7];
            dz = du * du * eval[8];
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((1)+(0)*imax)  ] - xyz[3*((0)+(0)*imax)  ] - du*eval[3];
            dy = xyz[3*((1)+(0)*imax)+1] - xyz[3*((0)+(0)*imax)+1] - du*eval[4];
            dz = xyz[3*((1)+(0)*imax)+2] - xyz[3*((0)+(0)*imax)+2] - du*eval[5];
        } else {
            /* quadratic fit */
            u20    = knotu[5] - knotu[3];
            u21    = knotu[5] - knotu[4];
            rj[0]  = xyz[3*((1)+(0)*imax)  ]*u20*u20 -
                     xyz[3*((0)+(0)*imax)  ]*u21*u21 -
                     xyz[3*((2)+(0)*imax)  ]*du*du;
            rj[1]  = xyz[3*((1)+(0)*imax)+1]*u20*u20 -
                     xyz[3*((0)+(0)*imax)+1]*u21*u21 -
                     xyz[3*((2)+(0)*imax)+1]*du*du;
            rj[2]  = xyz[3*((1)+(0)*imax)+2]*u20*u20 -
                     xyz[3*((0)+(0)*imax)+2]*u21*u21 -
                     xyz[3*((2)+(0)*imax)+2]*du*du;
            rj[0] /= 2.0*u21*du;
            rj[1] /= 2.0*u21*du;
            rj[2] /= 2.0*u21*du;
            dx     = rj[0] - xyz[3*((0)+(0)*imax)  ] - 0.5 * u20 * eval[3];
            dy     = rj[1] - xyz[3*((0)+(0)*imax)+1] - 0.5 * u20 * eval[4];
            dz     = rj[2] - xyz[3*((0)+(0)*imax)+2] - 0.5 * u20 * eval[5];
        }
        /* match input tangent */
        if (wesT != NULL)
            if ((wesT[0] != 0.0) || (wesT[1] != 0.0) || (wesT[2] != 0.0)) {
                dx = (wesT[0] - eval[3])*du;
                dy = (wesT[1] - eval[4])*du;
                dz = (wesT[2] - eval[5])*du;
            }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((1)+(0)*(imax+2))  ] += RELAX * dx;
        cp[3*((1)+(0)*(imax+2))+1] += RELAX * dy;
        cp[3*((1)+(0)*(imax+2))+2] += RELAX * dz;
      
        /* point C */
        dv = knotv[4] - knotv[3];
        if (endj == 0) {
            /* d2/dv2 = 0 */
            dx = dv * dv * eval[15];
            dy = dv * dv * eval[16];
            dz = dv * dv * eval[17];
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((0)+(1)*imax)  ] - xyz[3*((0)+(0)*imax)  ] - dv*eval[ 9];
            dy = xyz[3*((0)+(1)*imax)+1] - xyz[3*((0)+(0)*imax)+1] - dv*eval[10];
            dz = xyz[3*((0)+(1)*imax)+2] - xyz[3*((0)+(0)*imax)+2] - dv*eval[11];
        } else {
            /* quadratic fit */
            u20    = knotv[5] - knotv[3];
            u21    = knotv[5] - knotv[4];
            rj[0]  = xyz[3*((0)+(1)*imax)  ]*u20*u20 -
                     xyz[3*((0)+(0)*imax)  ]*u21*u21 -
                     xyz[3*((0)+(2)*imax)  ]*dv*dv;
            rj[1]  = xyz[3*((0)+(1)*imax)+1]*u20*u20 -
                     xyz[3*((0)+(0)*imax)+1]*u21*u21 -
                     xyz[3*((0)+(2)*imax)+1]*dv*dv;
            rj[2]  = xyz[3*((0)+(1)*imax)+2]*u20*u20 -
                     xyz[3*((0)+(0)*imax)+2]*u21*u21 -
                     xyz[3*((0)+(2)*imax)+2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            dx     = rj[0] - xyz[3*((0)+(0)*imax)  ] - 0.5 * u20 * eval[ 9];
            dy     = rj[1] - xyz[3*((0)+(0)*imax)+1] - 0.5 * u20 * eval[10];
            dz     = rj[2] - xyz[3*((0)+(0)*imax)+2] - 0.5 * u20 * eval[11];
        }
        if (south != NULL) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
            x2[0]    = xyz[3*((0)+(1)*(imax))  ] - xyz[3*((0)+(0)*(imax))  ];
            x2[1]    = xyz[3*((0)+(1)*(imax))+1] - xyz[3*((0)+(0)*(imax))+1];
            x2[2]    = xyz[3*((0)+(1)*(imax))+2] - xyz[3*((0)+(0)*(imax))+2];
            x2[ms]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rs[0], x2);
            nell[1]  = DOT(rs[1], x2);
            nell[2]  = DOT(rs[2], x2);
            r        = EG_getEllRad(south, nell);
            normnell = 1.0;
            x1[0]    = cpsav[3*((0)+(2)*(imax+2))  ] -
                       cpsav[3*((0)+(0)*(imax+2))  ];
            x1[1]    = cpsav[3*((0)+(2)*(imax+2))+1] -
                       cpsav[3*((0)+(0)*(imax+2))+1];
            x1[2]    = cpsav[3*((0)+(2)*(imax+2))+2] -
                       cpsav[3*((0)+(0)*(imax+2))+2];
            tt       = sqrt(r*con2*fabs(DOT(x1,snor))/con)/normnell;
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", 0,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((0)+(0)*imax)  ] + tt*nell[0] -
                                           cp[3*((0)+(1)*(imax+2))  ];
            dy = xyz[3*((0)+(0)*imax)+1] + tt*nell[1] -
                                           cp[3*((0)+(1)*(imax+2))+1];
            dz = xyz[3*((0)+(0)*imax)+2] + tt*nell[2] -
                                           cp[3*((0)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(1)*(imax+2))  ] += dx;
            cp[3*((0)+(1)*(imax+2))+1] += dy;
            cp[3*((0)+(1)*(imax+2))+2] += dz;
        } else {
            if (snor != NULL) {
                dx = (snor[0] - eval[ 9])*dv;
                dy = (snor[1] - eval[10])*dv;
                dz = (snor[2] - eval[11])*dv;
            }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((0)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((0)+(1)*(imax+2))+2] += RELAX * dz;
        }
      
        /* point D [opposite sign] */
        dv = knotv[3] - knotv[4];
        if (endc == 0) {
            /* d2/du/dv = 0 */
            dx = du * dv * eval[12];
            dy = du * dv * eval[13];
            dz = du * dv * eval[14];
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((1)+(1)*imax)  ] - xyz[3*((0)+(1)*imax)  ]) -
                 (xyz[3*((1)+(0)*imax)  ] - xyz[3*((0)+(0)*imax)  ]) +
                 du * dv * eval[12];
            dy = (xyz[3*((1)+(1)*imax)+1] - xyz[3*((0)+(1)*imax)+1]) -
                 (xyz[3*((1)+(0)*imax)+1] - xyz[3*((0)+(0)*imax)+1]) +
                 du * dv * eval[13];
            dz = (xyz[3*((1)+(1)*imax)+2] - xyz[3*((0)+(1)*imax)+2]) -
                 (xyz[3*((1)+(0)*imax)+2] - xyz[3*((0)+(0)*imax)+2]) +
                 du * dv * eval[14];
        } else {
            /* quadratic fit */
            u20 = knotu[5] - knotu[3];
            u21 = knotu[5] - knotu[4];
            for (tanOK = j = 0; j < 3; j++)
                if (wesT != NULL)
                    if ((wesT[3*j  ] != 0.0) || (wesT[3*j+1] != 0.0) ||
                        (wesT[3*j+2] != 0.0)) tanOK++;
            for (j = 0; j < 3; j++) {
                if (tanOK == 3) {
                    t[j][0] = wesT[3*j  ]*du;
                    t[j][1] = wesT[3*j+1]*du;
                    t[j][2] = wesT[3*j+2]*du;
                    continue;
                }
                rj[0]   = xyz[3*((1)+(j)*imax)  ]*u20*u20 -
                          xyz[3*((0)+(j)*imax)  ]*u21*u21 -
                          xyz[3*((2)+(j)*imax)  ]*du*du;
                rj[1]   = xyz[3*((1)+(j)*imax)+1]*u20*u20 -
                          xyz[3*((0)+(j)*imax)+1]*u21*u21 -
                          xyz[3*((2)+(j)*imax)+1]*du*du;
                rj[2]   = xyz[3*((1)+(j)*imax)+2]*u20*u20 -
                          xyz[3*((0)+(j)*imax)+2]*u21*u21 -
                          xyz[3*((2)+(j)*imax)+2]*du*du;
                t[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((0)+(j)*imax)  ];
                t[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((0)+(j)*imax)+1];
                t[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((0)+(j)*imax)+2];
            }
            u20    = knotv[5] - knotv[3];
            u21    = knotv[5] - knotv[4];
            rj[0]  = t[1][0]*u20*u20 - t[0][0]*u21*u21 - t[2][0]*dv*dv;
            rj[1]  = t[1][1]*u20*u20 - t[0][1]*u21*u21 - t[2][1]*dv*dv;
            rj[2]  = t[1][2]*u20*u20 - t[0][2]*u21*u21 - t[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[5] - knotu[3]);
            if (tanOK == 3) u21 = du;
            dx     = -rj[0] - t[0][0] - 0.5*u20*u21*eval[12];
            dy     = -rj[1] - t[0][1] - 0.5*u20*u21*eval[13];
            dz     = -rj[2] - t[0][2] - 0.5*u20*u21*eval[14];
        }
        if (south != NULL) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
            dist  = 0.5* (knotu[3] + knotu[4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
            x2[0]    = 0.5*(xyz[3*((1)+(1)*(imax))  ] +
                            xyz[3*((0)+(1)*(imax))  ]) -
                       0.5*(xyz[3*((1)+(0)*(imax))  ] +
                            xyz[3*((0)+(0)*(imax))  ]);
            x2[1]    = 0.5*(xyz[3*((1)+(1)*(imax))+1] +
                            xyz[3*((0)+(1)*(imax))+1]) -
                       0.5*(xyz[3*((1)+(0)*(imax))+1] +
                            xyz[3*((0)+(0)*(imax))+1]);
            x2[2]    = 0.5*(xyz[3*((1)+(1)*(imax))+2] +
                            xyz[3*((0)+(1)*(imax))+2]) -
                       0.5*(xyz[3*((1)+(0)*(imax))+2] +
                            xyz[3*((0)+(0)*(imax))+2]);
            x2[ms]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rs[0], x2);
            nell[1]  = DOT(rs[1], x2);
            nell[2]  = DOT(rs[2], x2);
            r        = EG_getEllRad(south, nell);
            normnell = 1.0;
            norm[0]  = norm[1] = norm[2] = 0.0;
            for (kk = 0; kk < 4; kk++) {
                if (kk == 1) continue;
                x0[0]    = cpsav[3*((kk)+(1)*(imax+2))  ] -
                           cpsav[3*((kk)+(0)*(imax+2))  ];
                x0[1]    = cpsav[3*((kk)+(1)*(imax+2))+1] -
                           cpsav[3*((kk)+(0)*(imax+2))+1];
                x0[2]    = cpsav[3*((kk)+(1)*(imax+2))+2] -
                           cpsav[3*((kk)+(0)*(imax+2))+2];
                norm[0] += x0[0]*basis[kk];
                norm[1] += x0[1]*basis[kk];
                norm[2] += x0[2]*basis[kk];
            }
            dD = con*normnell*normnell*basis[1]*basis[1];
            eE = 2.0*con*basis[1]*DOT(nell, norm);
            F  = con*DOT(norm, norm);
            G  = 0.0;
            for (kk = 0; kk < 4; kk++) {
                x1[0] = cpsav[3*((kk)+(2)*(imax+2))  ] -
                        cpsav[3*((kk)+(0)*(imax+2))  ];
                x1[1] = cpsav[3*((kk)+(2)*(imax+2))+1] -
                        cpsav[3*((kk)+(0)*(imax+2))+1];
                x1[2] = cpsav[3*((kk)+(2)*(imax+2))+2] -
                        cpsav[3*((kk)+(0)*(imax+2))+2];
                G += fabs(DOT(x1, snor))*basis[kk];
            }
            G *= con2;
            tt = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", 1,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((0)+(0)*imax)  ] + tt*nell[0] -
                                           cp[3*((1)+(1)*(imax+2))  ];
            dy = xyz[3*((0)+(0)*imax)+1] + tt*nell[1] -
                                           cp[3*((1)+(1)*(imax+2))+1];
            dz = xyz[3*((0)+(0)*imax)+2] + tt*nell[2] -
                                           cp[3*((1)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(1)*(imax+2))  ] += dx;
            cp[3*((1)+(1)*(imax+2))+1] += dy;
            cp[3*((1)+(1)*(imax+2))+2] += dz;
        } else {
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((1)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((1)+(1)*(imax+2))+2] += RELAX * dz;
        }
      
        /* match south & E points */
        dv    = knotv[4] - knotv[3];
        con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
        con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
        uv[1] = knotv[3];
        for (i = 1; i < imax-1; i++) {
            uv[0] = knotu[i+3];
            EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
            dx = xyz[3*((i)+(0)*imax)  ] - eval[0];
            dy = xyz[3*((i)+(0)*imax)+1] - eval[1];
            dz = xyz[3*((i)+(0)*imax)+2] - eval[2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(0)*(imax+2))  ] += dx;
            cp[3*((i+1)+(0)*(imax+2))+1] += dy;
            cp[3*((i+1)+(0)*(imax+2))+2] += dz;

            if (endj == 0) {
                /* d2/dv2 = 0 */
                dx = dv * dv * eval[15];
                dy = dv * dv * eval[16];
                dz = dv * dv * eval[17];
            } else if (endj == 1) {
                /* match FD d/dv */
                dx = xyz[3*((i)+(1)*imax)  ] -
                     xyz[3*((i)+(0)*imax)  ] - dv * eval[ 9];
                dy = xyz[3*((i)+(1)*imax)+1] -
                     xyz[3*((i)+(0)*imax)+1] - dv * eval[10];
                dz = xyz[3*((i)+(1)*imax)+2] -
                     xyz[3*((i)+(0)*imax)+2] - dv * eval[11];
            } else {
                /* quadratic fit */
                u20    = knotv[5] - knotv[3];
                u21    = knotv[5] - knotv[4];
                rj[0]  = xyz[3*((i)+(1)*imax)  ]*u20*u20 -
                         xyz[3*((i)+(0)*imax)  ]*u21*u21 -
                         xyz[3*((i)+(2)*imax)  ]*dv*dv;
                rj[1]  = xyz[3*((i)+(1)*imax)+1]*u20*u20 -
                         xyz[3*((i)+(0)*imax)+1]*u21*u21 -
                         xyz[3*((i)+(2)*imax)+1]*dv*dv;
                rj[2]  = xyz[3*((i)+(1)*imax)+2]*u20*u20 -
                         xyz[3*((i)+(0)*imax)+2]*u21*u21 -
                         xyz[3*((i)+(2)*imax)+2]*dv*dv;
                rj[0] /= 2.0*u21*dv;
                rj[1] /= 2.0*u21*dv;
                rj[2] /= 2.0*u21*dv;
                dx     = rj[0] - xyz[3*((i)+(0)*imax)  ] - 0.5 * u20 * eval[ 9];
                dy     = rj[1] - xyz[3*((i)+(0)*imax)+1] - 0.5 * u20 * eval[10];
                dz     = rj[2] - xyz[3*((i)+(0)*imax)+2] - 0.5 * u20 * eval[11];
            }
            if (south != NULL) {
                EG_NzeroBSp(iknot, knotu, knotu[i+3], &kk, basis);
                x2[0]    = xyz[3*((i)+(1)*(imax))  ] - xyz[3*((i)+(0)*(imax))  ];
                x2[1]    = xyz[3*((i)+(1)*(imax))+1] - xyz[3*((i)+(0)*(imax))+1];
                x2[2]    = xyz[3*((i)+(1)*(imax))+2] - xyz[3*((i)+(0)*(imax))+2];
                x2[ms]   = 0.0;
                dist     = sqrt(DOT(x2, x2));
                x2[0]   /= dist;
                x2[1]   /= dist;
                x2[2]   /= dist;
                nell[0]  = DOT(rs[0], x2);
                nell[1]  = DOT(rs[1], x2);
                nell[2]  = DOT(rs[2], x2);
                r        = EG_getEllRad(south, nell);
                normnell = 1.0;
                x0[0]    = cpsav[3*((i)  +(1)*(imax+2))  ] -
                           cpsav[3*((i)  +(0)*(imax+2))  ];
                x0[1]    = cpsav[3*((i)  +(1)*(imax+2))+1] -
                           cpsav[3*((i)  +(0)*(imax+2))+1];
                x0[2]    = cpsav[3*((i)  +(1)*(imax+2))+2] -
                           cpsav[3*((i)  +(0)*(imax+2))+2];
                x1[0]    = cpsav[3*((i+2)+(1)*(imax+2))  ] -
                           cpsav[3*((i+2)+(0)*(imax+2))  ];
                x1[1]    = cpsav[3*((i+2)+(1)*(imax+2))+1] -
                           cpsav[3*((i+2)+(0)*(imax+2))+1];
                x1[2]    = cpsav[3*((i+2)+(1)*(imax+2))+2] -
                           cpsav[3*((i+2)+(0)*(imax+2))+2];
                norm[0]  = x0[0]*basis[0] + x1[0]*basis[2];
                norm[1]  = x0[1]*basis[0] + x1[1]*basis[2];
                norm[2]  = x0[2]*basis[0] + x1[2]*basis[2];
                dD = con*normnell*normnell*basis[1]*basis[1];
                eE = 2.0*con*basis[1]*DOT(nell, norm);
                F  = con*DOT(norm, norm);
                G  = 0.0;
                for (kk = 0; kk < 3; kk++) {
                    x1[0] = cpsav[3*((i+kk)+(2)*(imax+2))  ] -
                            cpsav[3*((i+kk)+(0)*(imax+2))  ];
                    x1[1] = cpsav[3*((i+kk)+(2)*(imax+2))+1] -
                            cpsav[3*((i+kk)+(0)*(imax+2))+1];
                    x1[2] = cpsav[3*((i+kk)+(2)*(imax+2))+2] -
                            cpsav[3*((i+kk)+(0)*(imax+2))+2];
                    G += fabs(DOT(x1, snor))*basis[kk];
                }
                G *= con2;
                tt = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
                printf(" south %d: %lf %lf %lf (iter = %d)\n", i+1,
                       value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                       iter);
#endif
                if (tt < 0.0) tt = 0.0;
                dx = xyz[3*((i)+(0)*imax)  ] + tt*nell[0] -
                                               cp[3*((i+1)+(1)*(imax+2))  ];
                dy = xyz[3*((i)+(0)*imax)+1] + tt*nell[1] -
                                               cp[3*((i+1)+(1)*(imax+2))+1];
                dz = xyz[3*((i)+(0)*imax)+2] + tt*nell[2] -
                                               cp[3*((i+1)+(1)*(imax+2))+2];
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(1)*(imax+2))  ] += dx;
                cp[3*((i+1)+(1)*(imax+2))+1] += dy;
                cp[3*((i+1)+(1)*(imax+2))+2] += dz;
            } else {
                if (snor != NULL) {
                    dx = (snor[0] - eval[ 9])*dv;
                    dy = (snor[1] - eval[10])*dv;
                    dz = (snor[2] - eval[11])*dv;
                }
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(1)*(imax+2))  ] += RELAX * dx;
                cp[3*((i+1)+(1)*(imax+2))+1] += RELAX * dy;
                cp[3*((i+1)+(1)*(imax+2))+2] += RELAX * dz;
            }
        }
      
        /* point B */
        uv[0] = knotu[imax+2];
        uv[1] = knotv[3];
        EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
        du = knotu[imax+2] - knotu[imax+1];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * eval[6];
            dy = du * du * eval[7];
            dz = du * du * eval[8];
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((imax-2)+(0)*imax)  ] -
                 xyz[3*((imax-1)+(0)*imax)  ] + du * eval[3];
            dy = xyz[3*((imax-2)+(0)*imax)+1] -
                 xyz[3*((imax-1)+(0)*imax)+1] + du * eval[4];
            dz = xyz[3*((imax-2)+(0)*imax)+2] -
                 xyz[3*((imax-1)+(0)*imax)+2] + du * eval[5];
        } else {
            /* quadratic fit */
            u20    = knotu[imax+2] - knotu[imax];
            u21    = knotu[imax+1] - knotu[imax];
            rj[0]  = xyz[3*((imax-2)+(0)*imax)  ]*u20*u20 -
                     xyz[3*((imax-1)+(0)*imax)  ]*u21*u21 -
                     xyz[3*((imax-3)+(0)*imax)  ]*du*du;
            rj[1]  = xyz[3*((imax-2)+(0)*imax)+1]*u20*u20 -
                     xyz[3*((imax-1)+(0)*imax)+1]*u21*u21 -
                     xyz[3*((imax-3)+(0)*imax)+1]*du*du;
            rj[2]  = xyz[3*((imax-2)+(0)*imax)+2]*u20*u20 -
                     xyz[3*((imax-1)+(0)*imax)+2]*u21*u21 -
                     xyz[3*((imax-3)+(0)*imax)+2]*du*du;
            rj[0] /= 2.0*u21*du;
            rj[1] /= 2.0*u21*du;
            rj[2] /= 2.0*u21*du;
            dx     = rj[0] - xyz[3*((imax-1)+(0)*imax)  ] + 0.5 * u20 * eval[3];
            dy     = rj[1] - xyz[3*((imax-1)+(0)*imax)+1] + 0.5 * u20 * eval[4];
            dz     = rj[2] - xyz[3*((imax-1)+(0)*imax)+2] + 0.5 * u20 * eval[5];
        }
        /* match input tangent */
        if (easT != NULL)
            if ((easT[0] != 0.0) || (easT[1] != 0.0) || (easT[2] != 0.0)) {
                dx = (easT[0] + eval[3])*du;
                dy = (easT[1] + eval[4])*du;
                dz = (easT[2] + eval[5])*du;
            }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax)+(0)*(imax+2))  ] += RELAX * dx;
        cp[3*((imax)+(0)*(imax+2))+1] += RELAX * dy;
        cp[3*((imax)+(0)*(imax+2))+2] += RELAX * dz;
      
        /* point G */
        dv = knotv[4] - knotv[3];
        if (endj == 0) {
            /* d2/dv2 = 0 */
            dx = dv * dv * eval[15];
            dy = dv * dv * eval[16];
            dz = dv * dv * eval[17];
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((imax-1)+(1)*imax)  ] -
                 xyz[3*((imax-1)+(0)*imax)  ] - dv * eval[ 9];
            dy = xyz[3*((imax-1)+(1)*imax)+1] -
                 xyz[3*((imax-1)+(0)*imax)+1] - dv * eval[10];
            dz = xyz[3*((imax-1)+(1)*imax)+2] -
                 xyz[3*((imax-1)+(0)*imax)+2] - dv * eval[11];
        } else {
            /* quadratic fit */
            u20    = knotv[5] - knotv[3];
            u21    = knotv[5] - knotv[4];
            rj[0]  = xyz[3*((imax-1)+(1)*imax)  ]*u20*u20 -
                     xyz[3*((imax-1)+(0)*imax)  ]*u21*u21 -
                     xyz[3*((imax-1)+(2)*imax)  ]*dv*dv;
            rj[1]  = xyz[3*((imax-1)+(1)*imax)+1]*u20*u20 -
                     xyz[3*((imax-1)+(0)*imax)+1]*u21*u21 -
                     xyz[3*((imax-1)+(2)*imax)+1]*dv*dv;
            rj[2]  = xyz[3*((imax-1)+(1)*imax)+2]*u20*u20 -
                     xyz[3*((imax-1)+(0)*imax)+2]*u21*u21 -
                     xyz[3*((imax-1)+(2)*imax)+2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            dx     = rj[0] - xyz[3*((imax-1)+(0)*imax)  ] - 0.5 * u20 * eval[ 9];
            dy     = rj[1] - xyz[3*((imax-1)+(0)*imax)+1] - 0.5 * u20 * eval[10];
            dz     = rj[2] - xyz[3*((imax-1)+(0)*imax)+2] - 0.5 * u20 * eval[11];
        }
        if (south != NULL) {
            con      = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2     = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
            x2[0]    = xyz[3*((imax-1)+(1)*(imax))  ] -
                       xyz[3*((imax-1)+(0)*(imax))  ];
            x2[1]    = xyz[3*((imax-1)+(1)*(imax))+1] -
                       xyz[3*((imax-1)+(0)*(imax))+1];
            x2[2]    = xyz[3*((imax-1)+(1)*(imax))+2] -
                       xyz[3*((imax-1)+(0)*(imax))+2];
            x2[ms]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rs[0], x2);
            nell[1]  = DOT(rs[1], x2);
            nell[2]  = DOT(rs[2], x2);
            r        = EG_getEllRad(south, nell);
            normnell = 1.0;
            x1[0]    = cpsav[3*((imax+1)+(2)*(imax+2))  ] -
                       cpsav[3*((imax+1)+(0)*(imax+2))  ];
            x1[1]    = cpsav[3*((imax+1)+(2)*(imax+2))+1] -
                       cpsav[3*((imax+1)+(0)*(imax+2))+1];
            x1[2]    = cpsav[3*((imax+1)+(2)*(imax+2))+2] -
                       cpsav[3*((imax+1)+(0)*(imax+2))+2];
            tt       = sqrt(r*con2*fabs(DOT(x1,snor))/con)/normnell;
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", imax+1,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((imax-1)+(0)*imax)  ] + tt*nell[0] -
                                                cp[3*((imax+1)+(1)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(0)*imax)+1] + tt*nell[1] -
                                                cp[3*((imax+1)+(1)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(0)*imax)+2] + tt*nell[2] -
                                                cp[3*((imax+1)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(1)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(1)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(1)*(imax+2))+2] += dz;
        } else {
            if (snor != NULL) {
                dx = (snor[0] - eval[ 9])*dv;
                dy = (snor[1] - eval[10])*dv;
                dz = (snor[2] - eval[11])*dv;
            }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax+1)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax+1)+(1)*(imax+2))+2] += RELAX * dz;
        }

        /* point F */
        if (endc == 0) {
            /* d2/du/dv = 0 */
            dx = du * dv * eval[12];
            dy = du * dv * eval[13];
            dz = du * dv * eval[14];
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((imax-2)+(1)*imax)  ] - xyz[3*((imax-1)+(1)*imax)  ]) -
                 (xyz[3*((imax-2)+(0)*imax)  ] - xyz[3*((imax-1)+(0)*imax)  ]) +
                 du * dv * eval[12];
            dy = (xyz[3*((imax-2)+(1)*imax)+1] - xyz[3*((imax-1)+(1)*imax)+1]) -
                 (xyz[3*((imax-2)+(0)*imax)+1] - xyz[3*((imax-1)+(0)*imax)+1]) +
                 du * dv * eval[13];
            dz = (xyz[3*((imax-2)+(1)*imax)+2] - xyz[3*((imax-1)+(1)*imax)+2]) -
                 (xyz[3*((imax-2)+(0)*imax)+2] - xyz[3*((imax-1)+(0)*imax)+2]) +
                 du * dv * eval[14];
        } else {
            /* quadratic fit */
            u20 = knotu[imax+2] - knotu[imax];
            u21 = knotu[imax+1] - knotu[imax];
            for (tanOK = j = 0; j < 3; j++)
                if (easT != NULL)
                    if ((easT[3*j  ] != 0.0) || (easT[3*j+1] != 0.0) ||
                        (easT[3*j+2] != 0.0)) tanOK++;
            for (j = 0; j < 3; j++) {
                if (tanOK == 3) {
                    t[j][0] = easT[3*j  ]*du;
                    t[j][1] = easT[3*j+1]*du;
                    t[j][2] = easT[3*j+2]*du;
                    continue;
                }
                rj[0]  = xyz[3*((imax-2)+(j)*imax)  ]*u20*u20 -
                         xyz[3*((imax-1)+(j)*imax)  ]*u21*u21 -
                         xyz[3*((imax-3)+(j)*imax)  ]*du*du;
                rj[1]  = xyz[3*((imax-2)+(j)*imax)+1]*u20*u20 -
                         xyz[3*((imax-1)+(j)*imax)+1]*u21*u21 -
                         xyz[3*((imax-3)+(j)*imax)+1]*du*du;
                rj[2]  = xyz[3*((imax-2)+(j)*imax)+2]*u20*u20 -
                         xyz[3*((imax-1)+(j)*imax)+2]*u21*u21 -
                         xyz[3*((imax-3)+(j)*imax)+2]*du*du;
                t[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((imax-1)+(j)*imax)  ];
                t[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((imax-1)+(j)*imax)+1];
                t[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((imax-1)+(j)*imax)+2];
            }
            u20    = knotv[5] - knotv[3];
            u21    = knotv[5] - knotv[4];
            rj[0]  = t[1][0]*u20*u20 - t[0][0]*u21*u21 - t[2][0]*dv*dv;
            rj[1]  = t[1][1]*u20*u20 - t[0][1]*u21*u21 - t[2][1]*dv*dv;
            rj[2]  = t[1][2]*u20*u20 - t[0][2]*u21*u21 - t[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[imax+2] - knotu[imax]);
            if (tanOK == 3) u21 = du;
            dx     = rj[0] - t[0][0] + 0.5*u20*u21*eval[12];
            dy     = rj[1] - t[0][1] + 0.5*u20*u21*eval[13];
            dz     = rj[2] - t[0][2] + 0.5*u20*u21*eval[14];
        }
        if (south != NULL) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
            dist  = 0.5*(knotu[iknot-5] + knotu[iknot-4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
            x2[0]    = 0.5*(xyz[3*((imax-2)+(1)*(imax))  ] +
                            xyz[3*((imax-1)+(1)*(imax))  ]) -
                       0.5*(xyz[3*((imax-2)+(0)*(imax))  ] +
                            xyz[3*((imax-1)+(0)*(imax))  ]);
            x2[1]    = 0.5*(xyz[3*((imax-2)+(1)*(imax))+1] +
                            xyz[3*((imax-1)+(1)*(imax))+1]) -
                       0.5*(xyz[3*((imax-2)+(0)*(imax))+1] +
                            xyz[3*((imax-1)+(0)*(imax))+1]);
            x2[2]    = 0.5*(xyz[3*((imax-2)+(1)*(imax))+2] +
                            xyz[3*((imax-1)+(1)*(imax))+2]) -
                       0.5*(xyz[3*((imax-2)+(0)*(imax))+2] +
                            xyz[3*((imax-1)+(0)*(imax))+2]);
            x2[ms]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rs[0], x2);
            nell[1]  = DOT(rs[1], x2);
            nell[2]  = DOT(rs[2], x2);
            r        = EG_getEllRad(south, nell);
            normnell = 1.0;
            norm[0]  = norm[1] = norm[2] = 0.0;
            for (kk = 0; kk < 4; kk++) {
                if (kk == 1) continue;
                x0[0]    = cpsav[3*((imax+1-kk)+(1)*(imax+2))  ] -
                           cpsav[3*((imax+1-kk)+(0)*(imax+2))  ];
                x0[1]    = cpsav[3*((imax+1-kk)+(1)*(imax+2))+1] -
                           cpsav[3*((imax+1-kk)+(0)*(imax+2))+1];
                x0[2]    = cpsav[3*((imax+1-kk)+(1)*(imax+2))+2] -
                           cpsav[3*((imax+1-kk)+(0)*(imax+2))+2];
                norm[0] += x0[0]*basis[3-kk];
                norm[1] += x0[1]*basis[3-kk];
                norm[2] += x0[2]*basis[3-kk];
            }
            dD = con*normnell*normnell*basis[2]*basis[2];
            eE = 2.0*con*basis[2]*DOT(nell, norm);
            F  = con*DOT(norm, norm);
            G  = 0.0;
            for (kk = 0; kk < 4; kk++) {
                x1[0] = cpsav[3*((imax+1-kk)+(2)*(imax+2))  ] -
                        cpsav[3*((imax+1-kk)+(0)*(imax+2))  ];
                x1[1] = cpsav[3*((imax+1-kk)+(2)*(imax+2))+1] -
                        cpsav[3*((imax+1-kk)+(0)*(imax+2))+1];
                x1[2] = cpsav[3*((imax+1-kk)+(2)*(imax+2))+2] -
                        cpsav[3*((imax+1-kk)+(0)*(imax+2))+2];
                G += fabs(DOT(x1, snor))*basis[3-kk];
            }
            G *= con2;
            tt = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", imax,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((imax-1)+(0)*imax)  ] + tt*nell[0] -
                                                cp[3*((imax)+(1)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(0)*imax)+1] + tt*nell[1] -
                                                cp[3*((imax)+(1)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(0)*imax)+2] + tt*nell[2] -
                                                cp[3*((imax)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(1)*(imax+2))  ] += dx;
            cp[3*((imax)+(1)*(imax+2))+1] += dy;
            cp[3*((imax)+(1)*(imax+2))+2] += dz;
        } else {
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax)+(1)*(imax+2))+2] += RELAX * dz;
        }
      
        /* point O */
        uv[0] = knotu[3];
        uv[1] = knotv[jmax+2];
        EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
        du = knotu[4] - knotu[3];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * eval[6];
            dy = du * du * eval[7];
            dz = du * du * eval[8];
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((1)+(jmax-1)*imax)  ] -
                 xyz[3*((0)+(jmax-1)*imax)  ] - du * eval[3];
            dy = xyz[3*((1)+(jmax-1)*imax)+1] -
                 xyz[3*((0)+(jmax-1)*imax)+1] - du * eval[4];
            dz = xyz[3*((1)+(jmax-1)*imax)+2] -
                 xyz[3*((0)+(jmax-1)*imax)+2] - du * eval[5];
        } else {
            /* quadratic fit */
            u20    = knotu[5] - knotu[3];
            u21    = knotu[5] - knotu[4];
            rj[0]  = xyz[3*((1)+(jmax-1)*imax)  ]*u20*u20 -
                     xyz[3*((0)+(jmax-1)*imax)  ]*u21*u21 -
                     xyz[3*((2)+(jmax-1)*imax)  ]*du*du;
            rj[1]  = xyz[3*((1)+(jmax-1)*imax)+1]*u20*u20 -
                     xyz[3*((0)+(jmax-1)*imax)+1]*u21*u21 -
                     xyz[3*((2)+(jmax-1)*imax)+1]*du*du;
            rj[2]  = xyz[3*((1)+(jmax-1)*imax)+2]*u20*u20 -
                     xyz[3*((0)+(jmax-1)*imax)+2]*u21*u21 -
                     xyz[3*((2)+(jmax-1)*imax)+2]*du*du;
            rj[0] /= 2.0*u21*du;
            rj[1] /= 2.0*u21*du;
            rj[2] /= 2.0*u21*du;
            dx     = rj[0] - xyz[3*((0)+(jmax-1)*imax)  ] - 0.5 * u20 * eval[3];
            dy     = rj[1] - xyz[3*((0)+(jmax-1)*imax)+1] - 0.5 * u20 * eval[4];
            dz     = rj[2] - xyz[3*((0)+(jmax-1)*imax)+2] - 0.5 * u20 * eval[5];
        }
        /* match input tangent */
        if (wesT != NULL)
            if ((wesT[3*(jmax-1)  ] != 0.0) || (wesT[3*(jmax-1)+1] != 0.0) ||
                (wesT[3*(jmax-1)+2] != 0.0)){
                dx = (wesT[3*(jmax-1)  ] - eval[3])*du;
                dy = (wesT[3*(jmax-1)+1] - eval[4])*du;
                dz = (wesT[3*(jmax-1)+2] - eval[5])*du;
            }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((1)+(jmax+1)*(imax+2))  ] += RELAX * dx;
        cp[3*((1)+(jmax+1)*(imax+2))+1] += RELAX * dy;
        cp[3*((1)+(jmax+1)*(imax+2))+2] += RELAX * dz;
      
        /* point J */
        dv = knotv[jmax+2] - knotv[jmax+1];
        if (endj == 0) {
            /* d2/dv2 = 0 */
            dx = dv * dv * eval[15];
            dy = dv * dv * eval[16];
            dz = dv * dv * eval[17];
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((0)+(jmax-2)*imax)  ] -
                 xyz[3*((0)+(jmax-1)*imax)  ] + dv * eval[ 9];
            dy = xyz[3*((0)+(jmax-2)*imax)+1] -
                 xyz[3*((0)+(jmax-1)*imax)+1] + dv * eval[10];
            dz = xyz[3*((0)+(jmax-2)*imax)+2] -
                 xyz[3*((0)+(jmax-1)*imax)+2] + dv * eval[11];
        } else {
            /* quadratic fit */
            u20    = knotv[jmax+2] - knotv[jmax];
            u21    = knotv[jmax+1] - knotv[jmax];
            rj[0]  = xyz[3*((0)+(jmax-2)*imax)  ]*u20*u20 -
                     xyz[3*((0)+(jmax-1)*imax)  ]*u21*u21 -
                     xyz[3*((0)+(jmax-3)*imax)  ]*dv*dv;
            rj[1]  = xyz[3*((0)+(jmax-2)*imax)+1]*u20*u20 -
                     xyz[3*((0)+(jmax-1)*imax)+1]*u21*u21 -
                     xyz[3*((0)+(jmax-3)*imax)+1]*dv*dv;
            rj[2]  = xyz[3*((0)+(jmax-2)*imax)+2]*u20*u20 -
                     xyz[3*((0)+(jmax-1)*imax)+2]*u21*u21 -
                     xyz[3*((0)+(jmax-3)*imax)+2]*dv*dv;
            rj[0]  /= 2.0*u21*dv;
            rj[1]  /= 2.0*u21*dv;
            rj[2]  /= 2.0*u21*dv;
            dx    = rj[0] - xyz[3*((0)+(jmax-1)*imax)  ] + 0.5 * u20 * eval[ 9];
            dy    = rj[1] - xyz[3*((0)+(jmax-1)*imax)+1] + 0.5 * u20 * eval[10];
            dz    = rj[2] - xyz[3*((0)+(jmax-1)*imax)+2] + 0.5 * u20 * eval[11];
        }
        if (north != NULL) {
            con      = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                            (knotv[jmax+1] - knotv[jmax+2]));
            con2     = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                            (knotv[jmax+1] - knotv[jmax+2]));
            x2[0]    = xyz[3*((0)+(jmax-2)*(imax))  ] -
                       xyz[3*((0)+(jmax-1)*(imax))  ];
            x2[1]    = xyz[3*((0)+(jmax-2)*(imax))+1] -
                       xyz[3*((0)+(jmax-1)*(imax))+1];
            x2[2]    = xyz[3*((0)+(jmax-2)*(imax))+2] -
                       xyz[3*((0)+(jmax-1)*(imax))+2];
            x2[mn]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rn[0], x2);
            nell[1]  = DOT(rn[1], x2);
            nell[2]  = DOT(rn[2], x2);
            r        = EG_getEllRad(north, nell);
            normnell = 1.0;
            x1[0]    = cpsav[3*((0)+(jmax-1)*(imax+2))  ] -
                       cpsav[3*((0)+(jmax+1)*(imax+2))  ];
            x1[1]    = cpsav[3*((0)+(jmax-1)*(imax+2))+1] -
                       cpsav[3*((0)+(jmax+1)*(imax+2))+1];
            x1[2]    = cpsav[3*((0)+(jmax-1)*(imax+2))+2] -
                       cpsav[3*((0)+(jmax+1)*(imax+2))+2];
            tt       = sqrt(r*con2*fabs(DOT(x1,nnor))/con)/normnell;
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", 0,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((0)+(jmax-1)*imax)  ] + tt*nell[0] -
                                                cp[3*((0)+(jmax)*(imax+2))  ];
            dy = xyz[3*((0)+(jmax-1)*imax)+1] + tt*nell[1] -
                                                cp[3*((0)+(jmax)*(imax+2))+1];
            dz = xyz[3*((0)+(jmax-1)*imax)+2] + tt*nell[2] -
                                                cp[3*((0)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(jmax)*(imax+2))  ] += dx;
            cp[3*((0)+(jmax)*(imax+2))+1] += dy;
            cp[3*((0)+(jmax)*(imax+2))+2] += dz;
        } else {
            if (nnor != NULL) {
                dx = (nnor[0] + eval[ 9])*dv;
                dy = (nnor[1] + eval[10])*dv;
                dz = (nnor[2] + eval[11])*dv;
            }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(jmax)*(imax+2))  ] += RELAX * dx;
            cp[3*((0)+(jmax)*(imax+2))+1] += RELAX * dy;
            cp[3*((0)+(jmax)*(imax+2))+2] += RELAX * dz;
        }

        /* point K */
        if (endc == 0) {
            /* d2/du/dv = 0 */
            dx = du * dv * eval[12];
            dy = du * dv * eval[13];
            dz = du * dv * eval[14];
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((1)+(jmax-2)*imax)  ] - xyz[3*((0)+(jmax-2)*imax)  ]) -
                 (xyz[3*((1)+(jmax-1)*imax)  ] - xyz[3*((0)+(jmax-1)*imax)  ]) +
                 du * dv * eval[12];
            dy = (xyz[3*((1)+(jmax-2)*imax)+1] - xyz[3*((0)+(jmax-2)*imax)+1]) -
                 (xyz[3*((1)+(jmax-1)*imax)+1] - xyz[3*((0)+(jmax-1)*imax)+1]) +
                 du * dv * eval[13];
            dz = (xyz[3*((1)+(jmax-2)*imax)+2] - xyz[3*((0)+(jmax-2)*imax)+2]) -
                 (xyz[3*((1)+(jmax-1)*imax)+2] - xyz[3*((0)+(jmax-1)*imax)+2]) +
                 du * dv * eval[14];
        } else {
            /* quadratic fit */
            u20 = knotu[5] - knotu[3];
            u21 = knotu[5] - knotu[4];
            for (tanOK = j = 0; j < 3; j++) {
                jj = jmax - j - 1;
                if (wesT != NULL)
                  if ((wesT[3*jj  ] != 0.0) || (wesT[3*jj+1] != 0.0) ||
                      (wesT[3*jj+2] != 0.0)) tanOK++;
            }
            for (j = 0; j < 3; j++) {
                jj = jmax - j - 1;
                if (tanOK == 3) {
                    t[j][0] = wesT[3*jj  ]*du;
                    t[j][1] = wesT[3*jj+1]*du;
                    t[j][2] = wesT[3*jj+2]*du;
                    continue;
                }
                rj[0]   = xyz[3*((1)+(jj)*imax)  ]*u20*u20 -
                          xyz[3*((0)+(jj)*imax)  ]*u21*u21 -
                          xyz[3*((2)+(jj)*imax)  ]*du*du;
                rj[1]   = xyz[3*((1)+(jj)*imax)+1]*u20*u20 -
                          xyz[3*((0)+(jj)*imax)+1]*u21*u21 -
                          xyz[3*((2)+(jj)*imax)+1]*du*du;
                rj[2]   = xyz[3*((1)+(jj)*imax)+2]*u20*u20 -
                          xyz[3*((0)+(jj)*imax)+2]*u21*u21 -
                          xyz[3*((2)+(jj)*imax)+2]*du*du;
                t[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((0)+(jj)*imax)  ];
                t[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((0)+(jj)*imax)+1];
                t[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((0)+(jj)*imax)+2];
            }
            u20    = knotv[jmax+2] - knotv[jmax];
            u21    = knotv[jmax+1] - knotv[jmax];
            rj[0]  = t[1][0]*u20*u20 - t[0][0]*u21*u21 - t[2][0]*dv*dv;
            rj[1]  = t[1][1]*u20*u20 - t[0][1]*u21*u21 - t[2][1]*dv*dv;
            rj[2]  = t[1][2]*u20*u20 - t[0][2]*u21*u21 - t[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[5] - knotu[3]);
            if (tanOK == 3) u21 = du;
            dx     = rj[0] - t[0][0] + 0.5*u20*u21*eval[12];
            dy     = rj[1] - t[0][1] + 0.5*u20*u21*eval[13];
            dz     = rj[2] - t[0][2] + 0.5*u20*u21*eval[14];
        }
        if (north != NULL) {
            con   = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
            con2  = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
            dist  = 0.5*(knotu[3] + knotu[4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
            x2[0]    = 0.5*(xyz[3*((1)+(jmax-2)*(imax))  ] +
                            xyz[3*((0)+(jmax-2)*(imax))  ]) -
                       0.5*(xyz[3*((1)+(jmax-1)*(imax))  ] +
                            xyz[3*((0)+(jmax-1)*(imax))  ]);
            x2[1]    = 0.5*(xyz[3*((1)+(jmax-2)*(imax))+1] +
                            xyz[3*((0)+(jmax-2)*(imax))+1]) -
                       0.5*(xyz[3*((1)+(jmax-1)*(imax))+1] +
                            xyz[3*((0)+(jmax-1)*(imax))+1]);
            x2[2]    = 0.5*(xyz[3*((1)+(jmax-2)*(imax))+2] +
                            xyz[3*((0)+(jmax-2)*(imax))+2]) -
                       0.5*(xyz[3*((1)+(jmax-1)*(imax))+2] +
                            xyz[3*((0)+(jmax-1)*(imax))+2]);
            x2[mn]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rn[0], x2);
            nell[1]  = DOT(rn[1], x2);
            nell[2]  = DOT(rn[2], x2);
            r        = EG_getEllRad(north, nell);
            normnell = 1.0;
            norm[0]  = norm[1] = norm[2] = 0.0;
            for (kk = 0; kk < 4; kk++) {
                if (kk == 1) continue;
                x0[0]    = cpsav[3*((kk)+(jmax  )*(imax+2))  ] -
                           cpsav[3*((kk)+(jmax+1)*(imax+2))  ];
                x0[1]    = cpsav[3*((kk)+(jmax  )*(imax+2))+1] -
                           cpsav[3*((kk)+(jmax+1)*(imax+2))+1];
                x0[2]    = cpsav[3*((kk)+(jmax  )*(imax+2))+2] -
                           cpsav[3*((kk)+(jmax+1)*(imax+2))+2];
                norm[0] += x0[0]*basis[kk];
                norm[1] += x0[1]*basis[kk];
                norm[2] += x0[2]*basis[kk];
            }
            dD = con*normnell*normnell*basis[1]*basis[1];
            eE = 2.0*con*basis[1]*DOT(nell, norm);
            F  = con*DOT(norm, norm);
            G  = 0.0;
            for (kk = 0; kk < 4; kk++) {
                x1[0] = cpsav[3*((kk)+(jmax-1)*(imax+2))  ] -
                        cpsav[3*((kk)+(jmax+1)*(imax+2))  ];
                x1[1] = cpsav[3*((kk)+(jmax-1)*(imax+2))+1] -
                        cpsav[3*((kk)+(jmax+1)*(imax+2))+1];
                x1[2] = cpsav[3*((kk)+(jmax-1)*(imax+2))+2] -
                        cpsav[3*((kk)+(jmax+1)*(imax+2))+2];
                G += fabs(DOT(x1, nnor))*basis[kk];
            }
            G *= con2;
            tt = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", 1,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((0)+(jmax-1)*imax)  ] + tt*nell[0] -
                                                cp[3*((1)+(jmax)*(imax+2))  ];
            dy = xyz[3*((0)+(jmax-1)*imax)+1] + tt*nell[1] -
                                                cp[3*((1)+(jmax)*(imax+2))+1];
            dz = xyz[3*((0)+(jmax-1)*imax)+2] + tt*nell[2] -
                                                cp[3*((1)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(jmax)*(imax+2))  ] += dx;
            cp[3*((1)+(jmax)*(imax+2))+1] += dy;
            cp[3*((1)+(jmax)*(imax+2))+2] += dz;
        } else {
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(jmax)*(imax+2))  ] += RELAX * dx;
            cp[3*((1)+(jmax)*(imax+2))+1] += RELAX * dy;
            cp[3*((1)+(jmax)*(imax+2))+2] += RELAX * dz;
        }

        /* match north & L points */
        dv    = knotv[jmax+2] - knotv[jmax+1];
        con   = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                     (knotv[jmax+1] - knotv[jmax+2]));
        con2  = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                     (knotv[jmax+1] - knotv[jmax+2]));
        uv[1] = knotv[jmax+2];
        for (i = 1; i < imax-1; i++) {
            uv[0] = knotu[i+3];
            EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
            dx = xyz[3*((i)+(jmax-1)*imax)  ] - eval[0];
            dy = xyz[3*((i)+(jmax-1)*imax)+1] - eval[1];
            dz = xyz[3*((i)+(jmax-1)*imax)+2] - eval[2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(jmax+1)*(imax+2))  ] += dx;
            cp[3*((i+1)+(jmax+1)*(imax+2))+1] += dy;
            cp[3*((i+1)+(jmax+1)*(imax+2))+2] += dz;

            if (endj == 0) {
                /* d2/dv2 = 0 */
                dx = dv * dv * eval[15];
                dy = dv * dv * eval[16];
                dz = dv * dv * eval[17];
            } else if (endj == 1) {
                /* match FD d/dv */
                dx = xyz[3*((i)+(jmax-2)*imax)  ] -
                     xyz[3*((i)+(jmax-1)*imax)  ] + dv * eval[ 9];
                dy = xyz[3*((i)+(jmax-2)*imax)+1] -
                     xyz[3*((i)+(jmax-1)*imax)+1] + dv * eval[10];
                dz = xyz[3*((i)+(jmax-2)*imax)+2] -
                     xyz[3*((i)+(jmax-1)*imax)+2] + dv * eval[11];
            } else {
                /* quadratic fit */
                u20    = knotv[jmax+2] - knotv[jmax];
                u21    = knotv[jmax+1] - knotv[jmax];
                rj[0]  = xyz[3*((i)+(jmax-2)*imax)  ]*u20*u20 -
                         xyz[3*((i)+(jmax-1)*imax)  ]*u21*u21 -
                         xyz[3*((i)+(jmax-3)*imax)  ]*dv*dv;
                rj[1]  = xyz[3*((i)+(jmax-2)*imax)+1]*u20*u20 -
                         xyz[3*((i)+(jmax-1)*imax)+1]*u21*u21 -
                         xyz[3*((i)+(jmax-3)*imax)+1]*dv*dv;
                rj[2]  = xyz[3*((i)+(jmax-2)*imax)+2]*u20*u20 -
                         xyz[3*((i)+(jmax-1)*imax)+2]*u21*u21 -
                         xyz[3*((i)+(jmax-3)*imax)+2]*dv*dv;
                rj[0] /= 2.0*u21*dv;
                rj[1] /= 2.0*u21*dv;
                rj[2] /= 2.0*u21*dv;
                dx     = rj[0] - xyz[3*((i)+(jmax-1)*imax)  ] + 0.5*u20*eval[ 9];
                dy     = rj[1] - xyz[3*((i)+(jmax-1)*imax)+1] + 0.5*u20*eval[10];
                dz     = rj[2] - xyz[3*((i)+(jmax-1)*imax)+2] + 0.5*u20*eval[11];
            }
            if (north != NULL) {
                EG_NzeroBSp(iknot, knotu, knotu[i+3], &kk, basis);
                x2[0]    = xyz[3*((i)+(jmax-2)*(imax))  ] -
                           xyz[3*((i)+(jmax-1)*(imax))  ];
                x2[1]    = xyz[3*((i)+(jmax-2)*(imax))+1] -
                           xyz[3*((i)+(jmax-1)*(imax))+1];
                x2[2]    = xyz[3*((i)+(jmax-2)*(imax))+2] -
                           xyz[3*((i)+(jmax-1)*(imax))+2];
                x2[mn]   = 0.0;
                dist     = sqrt(DOT(x2, x2));
                x2[0]   /= dist;
                x2[1]   /= dist;
                x2[2]   /= dist;
                nell[0]  = DOT(rn[0], x2);
                nell[1]  = DOT(rn[1], x2);
                nell[2]  = DOT(rn[2], x2);
                r        = EG_getEllRad(north, nell);
                normnell = 1.0;
                x0[0]    = cpsav[3*((i)  +(jmax  )*(imax+2))  ] -
                           cpsav[3*((i)  +(jmax+1)*(imax+2))  ];
                x0[1]    = cpsav[3*((i)  +(jmax  )*(imax+2))+1] -
                           cpsav[3*((i)  +(jmax+1)*(imax+2))+1];
                x0[2]    = cpsav[3*((i)  +(jmax  )*(imax+2))+2] -
                           cpsav[3*((i)  +(jmax+1)*(imax+2))+2];
                x1[0]    = cpsav[3*((i+2)+(jmax  )*(imax+2))  ] -
                           cpsav[3*((i+2)+(jmax+1)*(imax+2))  ];
                x1[1]    = cpsav[3*((i+2)+(jmax  )*(imax+2))+1] -
                           cpsav[3*((i+2)+(jmax+1)*(imax+2))+1];
                x1[2]    = cpsav[3*((i+2)+(jmax  )*(imax+2))+2] -
                           cpsav[3*((i+2)+(jmax+1)*(imax+2))+2];
                norm[0]  = x0[0]*basis[0] + x1[0]*basis[2];
                norm[1]  = x0[1]*basis[0] + x1[1]*basis[2];
                norm[2]  = x0[2]*basis[0] + x1[2]*basis[2];
                dD = con*normnell*normnell*basis[1]*basis[1];
                eE = 2.0*con*basis[1]*DOT(nell, norm);
                F  = con*DOT(norm, norm);
                G  = 0.0;
                for (kk = 0; kk < 3; kk++) {
                    x1[0] = cpsav[3*((i+kk)+(jmax-1)*(imax+2))  ] -
                            cpsav[3*((i+kk)+(jmax+1)*(imax+2))  ];
                    x1[1] = cpsav[3*((i+kk)+(jmax-1)*(imax+2))+1] -
                            cpsav[3*((i+kk)+(jmax+1)*(imax+2))+1];
                    x1[2] = cpsav[3*((i+kk)+(jmax-1)*(imax+2))+2] -
                            cpsav[3*((i+kk)+(jmax+1)*(imax+2))+2];
                    G += fabs(DOT(x1, nnor))*basis[kk];
                }
                G *= con2;
                tt = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
                printf(" north %d: %lf %lf %lf (iter = %d)\n", i+1,
                       value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                       iter);
#endif
                if (tt < 0.0) tt = 0.0;
                dx = xyz[3*((i)+(jmax-1)*imax)  ] + tt*nell[0] -
                                                cp[3*((i+1)+(jmax)*(imax+2))  ];
                dy = xyz[3*((i)+(jmax-1)*imax)+1] + tt*nell[1] -
                                                cp[3*((i+1)+(jmax)*(imax+2))+1];
                dz = xyz[3*((i)+(jmax-1)*imax)+2] + tt*nell[2] -
                                                cp[3*((i+1)+(jmax)*(imax+2))+2];
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(jmax)*(imax+2))  ] += dx;
                cp[3*((i+1)+(jmax)*(imax+2))+1] += dy;
                cp[3*((i+1)+(jmax)*(imax+2))+2] += dz;
            } else {
                if (nnor != NULL) {
                    dx = (nnor[0] + eval[ 9])*dv;
                    dy = (nnor[1] + eval[10])*dv;
                    dz = (nnor[2] + eval[11])*dv;
                }
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(jmax)*(imax+2))  ] += RELAX * dx;
                cp[3*((i+1)+(jmax)*(imax+2))+1] += RELAX * dy;
                cp[3*((i+1)+(jmax)*(imax+2))+2] += RELAX * dz;
            }
        }

        /* point P */
        uv[0] = knotu[imax+2];
        uv[1] = knotv[jmax+2];
        EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
        du = knotu[imax+2] - knotu[imax+1];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * eval[6];
            dy = du * du * eval[7];
            dz = du * du * eval[8];
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((imax-2)+(jmax-1)*imax)  ] -
                 xyz[3*((imax-1)+(jmax-1)*imax)  ] + du * eval[3];
            dy = xyz[3*((imax-2)+(jmax-1)*imax)+1] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+1] + du * eval[4];
            dz = xyz[3*((imax-2)+(jmax-1)*imax)+2] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+2] + du * eval[5];
        } else {
            /* quadratic fit */
            u20    = knotu[imax+2] - knotu[imax];
            u21    = knotu[imax+1] - knotu[imax];
            rj[0]  = xyz[3*((imax-2)+(jmax-1)*imax)  ]*u20*u20 -
                     xyz[3*((imax-1)+(jmax-1)*imax)  ]*u21*u21 -
                     xyz[3*((imax-3)+(jmax-1)*imax)  ]*du*du;
            rj[1]  = xyz[3*((imax-2)+(jmax-1)*imax)+1]*u20*u20 -
                     xyz[3*((imax-1)+(jmax-1)*imax)+1]*u21*u21 -
                     xyz[3*((imax-3)+(jmax-1)*imax)+1]*du*du;
            rj[2]  = xyz[3*((imax-2)+(jmax-1)*imax)+2]*u20*u20 -
                     xyz[3*((imax-1)+(jmax-1)*imax)+2]*u21*u21 -
                     xyz[3*((imax-3)+(jmax-1)*imax)+2]*du*du;
            rj[0] /= 2.0*u21*du;
            rj[1] /= 2.0*u21*du;
            rj[2] /= 2.0*u21*du;
            dx     = rj[0] - xyz[3*((imax-1)+(jmax-1)*imax)  ] + 0.5*u20*eval[3];
            dy     = rj[1] - xyz[3*((imax-1)+(jmax-1)*imax)+1] + 0.5*u20*eval[4];
            dz     = rj[2] - xyz[3*((imax-1)+(jmax-1)*imax)+2] + 0.5*u20*eval[5];
        }
        /* match input tangent */
        if (easT != NULL)
            if ((easT[3*(jmax-1)  ] != 0.0) || (easT[3*(jmax-1)+1] != 0.0) ||
                (easT[3*(jmax-1)+2] != 0.0)) {
                dx = (easT[3*(jmax-1)  ] + eval[3])*du;
                dy = (easT[3*(jmax-1)+1] + eval[4])*du;
                dz = (easT[3*(jmax-1)+2] + eval[5])*du;
            }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax)+(jmax+1)*(imax+2))  ] += RELAX * dx;
        cp[3*((imax)+(jmax+1)*(imax+2))+1] += RELAX * dy;
        cp[3*((imax)+(jmax+1)*(imax+2))+2] += RELAX * dz;
      
        /* point N */
        dv = knotv[jmax+2] - knotv[jmax+1];
        if (endj == 0) {
            /* d2/dv2 = 0 */
            dx = dv * dv * eval[15];
            dy = dv * dv * eval[16];
            dz = dv * dv * eval[17];
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((imax-1)+(jmax-2)*imax)  ] -
                 xyz[3*((imax-1)+(jmax-1)*imax)  ] + dv * eval[ 9];
            dy = xyz[3*((imax-1)+(jmax-2)*imax)+1] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+1] + dv * eval[10];
            dz = xyz[3*((imax-1)+(jmax-2)*imax)+2] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+2] + dv * eval[11];
        } else {
            /* quadratic fit */
            u20    = knotv[jmax+2] - knotv[jmax];
            u21    = knotv[jmax+1] - knotv[jmax];
            rj[0]  = xyz[3*((imax-1)+(jmax-2)*imax)  ]*u20*u20 -
                     xyz[3*((imax-1)+(jmax-1)*imax)  ]*u21*u21 -
                     xyz[3*((imax-1)+(jmax-3)*imax)  ]*dv*dv;
            rj[1]  = xyz[3*((imax-1)+(jmax-2)*imax)+1]*u20*u20 -
                     xyz[3*((imax-1)+(jmax-1)*imax)+1]*u21*u21 -
                     xyz[3*((imax-1)+(jmax-3)*imax)+1]*dv*dv;
            rj[2]  = xyz[3*((imax-1)+(jmax-2)*imax)+2]*u20*u20 -
                     xyz[3*((imax-1)+(jmax-1)*imax)+2]*u21*u21 -
                     xyz[3*((imax-1)+(jmax-3)*imax)+2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            dx     = rj[0] - xyz[3*((imax-1)+(jmax-1)*imax)  ] + 0.5*u20*eval[ 9];
            dy     = rj[1] - xyz[3*((imax-1)+(jmax-1)*imax)+1] + 0.5*u20*eval[10];
            dz     = rj[2] - xyz[3*((imax-1)+(jmax-1)*imax)+2] + 0.5*u20*eval[11];
        }
        if (north != NULL) {
            con      = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                            (knotv[jmax+1] - knotv[jmax+2]));
            con2     = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                            (knotv[jmax+1] - knotv[jmax+2]));
            x2[0]    = xyz[3*((imax-1)+(jmax-2)*(imax))  ] -
                       xyz[3*((imax-1)+(jmax-1)*(imax))  ];
            x2[1]    = xyz[3*((imax-1)+(jmax-2)*(imax))+1] -
                       xyz[3*((imax-1)+(jmax-1)*(imax))+1];
            x2[2]    = xyz[3*((imax-1)+(jmax-2)*(imax))+2] -
                       xyz[3*((imax-1)+(jmax-1)*(imax))+2];
            x2[mn]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rn[0], x2);
            nell[1]  = DOT(rn[1], x2);
            nell[2]  = DOT(rn[2], x2);
            r        = EG_getEllRad(north, nell);
            normnell = 1.0;
            x1[0]    = cpsav[3*((imax+1)+(jmax-1)*(imax+2))  ] -
                       cpsav[3*((imax+1)+(jmax+1)*(imax+2))  ];
            x1[1]    = cpsav[3*((imax+1)+(jmax-1)*(imax+2))+1] -
                       cpsav[3*((imax+1)+(jmax+1)*(imax+2))+1];
            x1[2]    = cpsav[3*((imax+1)+(jmax-1)*(imax+2))+2] -
                       cpsav[3*((imax+1)+(jmax+1)*(imax+2))+2];
            tt       = sqrt(r*con2*fabs(DOT(x1,nnor))/con)/normnell;
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", imax+1,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((imax-1)+(jmax-1)*imax)  ] + tt*nell[0] -
                                             cp[3*((imax+1)+(jmax)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(jmax-1)*imax)+1] + tt*nell[1] -
                                             cp[3*((imax+1)+(jmax)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(jmax-1)*imax)+2] + tt*nell[2] -
                                             cp[3*((imax+1)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(jmax)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(jmax)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(jmax)*(imax+2))+2] += dz;
        } else {
            if (nnor != NULL) {
                dx = (nnor[0] + eval[ 9])*dv;
                dy = (nnor[1] + eval[10])*dv;
                dz = (nnor[2] + eval[11])*dv;
            }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(jmax)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax+1)+(jmax)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax+1)+(jmax)*(imax+2))+2] += RELAX * dz;
        }

        /* point M (opposite sign) */
        dv = knotv[jmax+1] - knotv[jmax+2];
        if (endc == 0) {
            /* d2/du/dv = 0 */
            dx = du * dv * eval[12];
            dy = du * dv * eval[13];
            dz = du * dv * eval[14];
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((imax-2)+(jmax-2)*imax)  ] -
                  xyz[3*((imax-1)+(jmax-2)*imax)  ])  -
                 (xyz[3*((imax-2)+(jmax-1)*imax)  ] -
                  xyz[3*((imax-1)+(jmax-1)*imax)  ])  + du * dv * eval[12];
            dy = (xyz[3*((imax-2)+(jmax-2)*imax)+1] -
                  xyz[3*((imax-1)+(jmax-2)*imax)+1])  -
                 (xyz[3*((imax-2)+(jmax-1)*imax)+1] -
                  xyz[3*((imax-1)+(jmax-1)*imax)+1])  + du * dv * eval[13];
            dz = (xyz[3*((imax-2)+(jmax-2)*imax)+2] -
                  xyz[3*((imax-1)+(jmax-2)*imax)+2])  -
                 (xyz[3*((imax-2)+(jmax-1)*imax)+2] -
                  xyz[3*((imax-1)+(jmax-1)*imax)+2])  + du * dv * eval[14];
        } else {
            /* quadratic fit */
            u20 = knotu[imax+2] - knotu[imax];
            u21 = knotu[imax+1] - knotu[imax];
            for (tanOK = j = 0; j < 3; j++) {
                jj = jmax - j - 1;
                if (easT != NULL)
                    if ((easT[3*jj  ] != 0.0) || (easT[3*jj+1] != 0.0) ||
                        (easT[3*jj+2] != 0.0)) tanOK++;
            }
            for (j = 0; j < 3; j++) {
                jj = jmax - j - 1;
                if (tanOK == 3) {
                    t[j][0] = easT[3*jj  ]*du;
                    t[j][1] = easT[3*jj+1]*du;
                    t[j][2] = easT[3*jj+2]*du;
                    continue;
                }
                rj[0]   = xyz[3*((imax-2)+(jj)*imax)  ]*u20*u20 -
                          xyz[3*((imax-1)+(jj)*imax)  ]*u21*u21 -
                          xyz[3*((imax-3)+(jj)*imax)  ]*du*du;
                rj[1]   = xyz[3*((imax-2)+(jj)*imax)+1]*u20*u20 -
                          xyz[3*((imax-1)+(jj)*imax)+1]*u21*u21 -
                          xyz[3*((imax-3)+(jj)*imax)+1]*du*du;
                rj[2]   = xyz[3*((imax-2)+(jj)*imax)+2]*u20*u20 -
                          xyz[3*((imax-1)+(jj)*imax)+2]*u21*u21 -
                          xyz[3*((imax-3)+(jj)*imax)+2]*du*du;
                t[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((imax-1)+(jj)*imax)  ];
                t[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((imax-1)+(jj)*imax)+1];
                t[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((imax-1)+(jj)*imax)+2];
            }
            u20    = knotv[jmax+2] - knotv[jmax];
            u21    = knotv[jmax+1] - knotv[jmax];
            rj[0]  = t[1][0]*u20*u20 - t[0][0]*u21*u21 - t[2][0]*dv*dv;
            rj[1]  = t[1][1]*u20*u20 - t[0][1]*u21*u21 - t[2][1]*dv*dv;
            rj[2]  = t[1][2]*u20*u20 - t[0][2]*u21*u21 - t[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[imax+2] - knotu[imax]);
            if (tanOK == 3) u21 = du;
            dx     = -rj[0] - t[0][0] - 0.5*u20*u21*eval[12];
            dy     = -rj[1] - t[0][1] - 0.5*u20*u21*eval[13];
            dz     = -rj[2] - t[0][2] - 0.5*u20*u21*eval[14];
        }
        if (north != NULL) {
            con  = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                        (knotv[jmax+1] - knotv[jmax+2]));
            con2 = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                        (knotv[jmax+1] - knotv[jmax+2]));
            dist  = 0.5*(knotu[iknot-5] + knotu[iknot-4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
            x2[0]    = 0.5*(xyz[3*((imax-2)+(jmax-2)*(imax))  ] +
                            xyz[3*((imax-1)+(jmax-2)*(imax))  ]) -
                       0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))  ] +
                            xyz[3*((imax-1)+(jmax-1)*(imax))  ]);
            x2[1]    = 0.5*(xyz[3*((imax-2)+(jmax-2)*(imax))+1] +
                            xyz[3*((imax-1)+(jmax-2)*(imax))+1]) -
                       0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))+1] +
                            xyz[3*((imax-1)+(jmax-1)*(imax))+1]);
            x2[2]    = 0.5*(xyz[3*((imax-2)+(jmax-2)*(imax))+2] +
                            xyz[3*((imax-1)+(jmax-2)*(imax))+2]) -
                       0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))+2] +
                            xyz[3*((imax-1)+(jmax-1)*(imax))+2]);
            x2[mn]   = 0.0;
            dist     = sqrt(DOT(x2, x2));
            x2[0]   /= dist;
            x2[1]   /= dist;
            x2[2]   /= dist;
            nell[0]  = DOT(rn[0], x2);
            nell[1]  = DOT(rn[1], x2);
            nell[2]  = DOT(rn[2], x2);
            r        = EG_getEllRad(north, nell);
            normnell = 1.0;
            norm[0]  = norm[1] = norm[2] = 0.0;
            for (kk = 0; kk < 4; kk++) {
                if (kk == 1) continue;
                x0[0]    = cpsav[3*((imax+1-kk)+(jmax  )*(imax+2))  ] -
                           cpsav[3*((imax+1-kk)+(jmax+1)*(imax+2))  ];
                x0[1]    = cpsav[3*((imax+1-kk)+(jmax  )*(imax+2))+1] -
                           cpsav[3*((imax+1-kk)+(jmax+1)*(imax+2))+1];
                x0[2]    = cpsav[3*((imax+1-kk)+(jmax  )*(imax+2))+2] -
                           cpsav[3*((imax+1-kk)+(jmax+1)*(imax+2))+2];
                norm[0] += x0[0]*basis[3-kk];
                norm[1] += x0[1]*basis[3-kk];
                norm[2] += x0[2]*basis[3-kk];
            }
            dD = con*normnell*normnell*basis[2]*basis[2];
            eE = 2.0*con*basis[2]*DOT(nell, norm);
            F  = con*DOT(norm, norm);
            G  = 0.0;
            for (kk = 0; kk < 4; kk++) {
                x1[0] = cpsav[3*((imax+1-kk)+(jmax-1)*(imax+2))  ] -
                        cpsav[3*((imax+1-kk)+(jmax+1)*(imax+2))  ];
                x1[1] = cpsav[3*((imax+1-kk)+(jmax-1)*(imax+2))+1] -
                        cpsav[3*((imax+1-kk)+(jmax+1)*(imax+2))+1];
                x1[2] = cpsav[3*((imax+1-kk)+(jmax-1)*(imax+2))+2] -
                        cpsav[3*((imax+1-kk)+(jmax+1)*(imax+2))+2];
                G += fabs(DOT(x1, nnor))*basis[3-kk];
            }
            G *= con2;
            tt = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", imax,
                   value(tt*nell[0]), value(tt*nell[1]), value(tt*nell[2]),
                   iter);
#endif
            if (tt < 0.0) tt = 0.0;
            dx = xyz[3*((imax-1)+(jmax-1)*imax)  ] + tt*nell[0] -
                                               cp[3*((imax)+(jmax)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(jmax-1)*imax)+1] + tt*nell[1] -
                                               cp[3*((imax)+(jmax)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(jmax-1)*imax)+2] + tt*nell[2] -
                                               cp[3*((imax)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(jmax)*(imax+2))  ] += dx;
            cp[3*((imax)+(jmax)*(imax+2))+1] += dy;
            cp[3*((imax)+(jmax)*(imax+2))+2] += dz;
        } else {
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(jmax)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax)+(jmax)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax)+(jmax)*(imax+2))+2] += RELAX * dz;
        }

        /* match west & H points */
        du    = knotu[4] - knotu[3];
        uv[0] = knotu[3];
        for (j = 1; j < jmax-1; j++) {
            uv[1] = knotv[j+3];
            EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
            if (vdata != NULL) {
                if (vdata[j] == +2) {
                    dv = knotv[j+4] - knotv[j+3];
                    dx = xyz[3*((0)+(j+1)*imax)  ] - xyz[3*((0)+(j)*imax)  ] -
                         dv * eval[ 9];
                    dy = xyz[3*((0)+(j+1)*imax)+1] - xyz[3*((0)+(j)*imax)+1] -
                         dv * eval[10];
                    dz = xyz[3*((0)+(j+1)*imax)+2] - xyz[3*((0)+(j)*imax)+2] -
                         dv * eval[11];
                    if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                    if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                    if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                    cp[3*((0)+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((0)+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((0)+(j+1)*(imax+2))+2] += RELAX * dz;
                    cp[3*((1)+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((1)+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((1)+(j+1)*(imax+2))+2] += RELAX * dz;
                    continue;
                }
                if (vdata[j] == -2) {
                    dv = knotv[j+3] - knotv[j+2];
                    dx = xyz[3*((0)+(j-1)*imax)  ] - xyz[3*((0)+(j)*imax)  ] +
                         dv * eval[ 9];
                    dy = xyz[3*((0)+(j-1)*imax)+1] - xyz[3*((0)+(j)*imax)+1] +
                         dv * eval[10];
                    dz = xyz[3*((0)+(j-1)*imax)+2] - xyz[3*((0)+(j)*imax)+2] +
                         dv * eval[11];
                    if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                    if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                    if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                    cp[3*((0)+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((0)+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((0)+(j+1)*(imax+2))+2] += RELAX * dz;
                    cp[3*((1)+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((1)+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((1)+(j+1)*(imax+2))+2] += RELAX * dz;
                    continue;
                }
            }
            dx = xyz[3*((0)+(j)*imax)  ] - eval[0];
            dy = xyz[3*((0)+(j)*imax)+1] - eval[1];
            dz = xyz[3*((0)+(j)*imax)+2] - eval[2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(j+1)*(imax+2))  ] += dx;
            cp[3*((0)+(j+1)*(imax+2))+1] += dy;
            cp[3*((0)+(j+1)*(imax+2))+2] += dz;

            if (endi == 0) {
                /* d2/du2 = 0 */
                dx = du * du * eval[6];
                dy = du * du * eval[7];
                dz = du * du * eval[8];
            } else if (endi == 1) {
                /* match FD d/du */
                dx = xyz[3*((1)+(j)*imax)  ] -
                     xyz[3*((0)+(j)*imax)  ] - du * eval[3];
                dy = xyz[3*((1)+(j)*imax)+1] -
                     xyz[3*((0)+(j)*imax)+1] - du * eval[4];
                dz = xyz[3*((1)+(j)*imax)+2] -
                     xyz[3*((0)+(j)*imax)+2] - du * eval[5];
            } else {
                /* quadratic fit */
                u20    = knotu[5] - knotu[3];
                u21    = knotu[5] - knotu[4];
                rj[0]  = xyz[3*((1)+(j)*imax)  ]*u20*u20 -
                         xyz[3*((0)+(j)*imax)  ]*u21*u21 -
                         xyz[3*((2)+(j)*imax)  ]*du*du;
                rj[1]  = xyz[3*((1)+(j)*imax)+1]*u20*u20 -
                         xyz[3*((0)+(j)*imax)+1]*u21*u21 -
                         xyz[3*((2)+(j)*imax)+1]*du*du;
                rj[2]  = xyz[3*((1)+(j)*imax)+2]*u20*u20 -
                         xyz[3*((0)+(j)*imax)+2]*u21*u21 -
                         xyz[3*((2)+(j)*imax)+2]*du*du;
                rj[0] /= 2.0*u21*du;
                rj[1] /= 2.0*u21*du;
                rj[2] /= 2.0*u21*du;
                dx     = rj[0] - xyz[3*((0)+(j)*imax)  ] - 0.5 * u20 * eval[3];
                dy     = rj[1] - xyz[3*((0)+(j)*imax)+1] - 0.5 * u20 * eval[4];
                dz     = rj[2] - xyz[3*((0)+(j)*imax)+2] - 0.5 * u20 * eval[5];
            }
            /* match input tangent */
            if (wesT != NULL)
                if ((wesT[3*j  ] != 0.0) || (wesT[3*j+1] != 0.0) ||
                    (wesT[3*j+2] != 0.0)) {
                    dx = (wesT[3*j  ] - eval[3])*du;
                    dy = (wesT[3*j+1] - eval[4])*du;
                    dz = (wesT[3*j+2] - eval[5])*du;
                }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(j+1)*(imax+2))  ] += RELAX * dx;
            cp[3*((1)+(j+1)*(imax+2))+1] += RELAX * dy;
            cp[3*((1)+(j+1)*(imax+2))+2] += RELAX * dz;
        }

        /* match east & I points */
        du = knotu[imax+2] - knotu[imax+1];
        uv[0] = knotu[imax+2];
        for (j = 1; j < jmax-1; j++) {
            uv[1] = knotv[j+3];
            EG_spline2dDeriv_impl(header, rvec, 2, uv, eval);
            if (vdata != NULL) {
                if (vdata[j] == +2) {
                    dv = knotv[j+4] - knotv[j+3];
                    dx = xyz[3*((imax-1)+(j+1)*imax)  ] -
                         xyz[3*((imax-1)+(j  )*imax)  ] - dv * eval[ 9];
                    dy = xyz[3*((imax-1)+(j+1)*imax)+1] -
                         xyz[3*((imax-1)+(j  )*imax)+1] - dv * eval[10];
                    dz = xyz[3*((imax-1)+(j+1)*imax)+2] -
                         xyz[3*((imax-1)+(j  )*imax)+2] - dv * eval[11];
                    if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                    if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                    if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                    cp[3*((imax+1)+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((imax+1)+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((imax+1)+(j+1)*(imax+2))+2] += RELAX * dz;
                    cp[3*((imax  )+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((imax  )+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((imax  )+(j+1)*(imax+2))+2] += RELAX * dz;
                    continue;
                }
                if (vdata[j] == -2) {
                    dv = knotv[j+3] - knotv[j+2];
                    dx = xyz[3*((imax-1)+(j-1)*imax)  ] -
                         xyz[3*((imax-1)+(j  )*imax)  ] + dv * eval[ 9];
                    dy = xyz[3*((imax-1)+(j-1)*imax)+1] -
                         xyz[3*((imax-1)+(j  )*imax)+1] + dv * eval[10];
                    dz = xyz[3*((imax-1)+(j-1)*imax)+2] -
                         xyz[3*((imax-1)+(j  )*imax)+2] + dv * eval[11];
                    if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                    if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                    if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                    cp[3*((imax+1)+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((imax+1)+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((imax+1)+(j+1)*(imax+2))+2] += RELAX * dz;
                    cp[3*((imax  )+(j+1)*(imax+2))  ] += RELAX * dx;
                    cp[3*((imax  )+(j+1)*(imax+2))+1] += RELAX * dy;
                    cp[3*((imax  )+(j+1)*(imax+2))+2] += RELAX * dz;
                    continue;
                }
            }
            dx = xyz[3*((imax-1)+(j)*imax)  ] - eval[0];
            dy = xyz[3*((imax-1)+(j)*imax)+1] - eval[1];
            dz = xyz[3*((imax-1)+(j)*imax)+2] - eval[2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(j+1)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(j+1)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(j+1)*(imax+2))+2] += dz;
              
            if (endi == 0) {
                /* d2/du2 = 0 */
                dx = du * du * eval[6];
                dy = du * du * eval[7];
                dz = du * du * eval[8];
            } else if (endi == 1) {
                /* match FD d/du */
                dx = xyz[3*((imax-2)+(j)*imax)  ] -
                     xyz[3*((imax-1)+(j)*imax)  ] + du * eval[3];
                dy = xyz[3*((imax-2)+(j)*imax)+1] -
                     xyz[3*((imax-1)+(j)*imax)+1] + du * eval[4];
                dz = xyz[3*((imax-2)+(j)*imax)+2] -
                     xyz[3*((imax-1)+(j)*imax)+2] + du * eval[5];
            } else {
                /* quadratic fit */
                u20    = knotu[imax+2] - knotu[imax];
                u21    = knotu[imax+1] - knotu[imax];
                rj[0]  = xyz[3*((imax-2)+(j)*imax)  ]*u20*u20 -
                         xyz[3*((imax-1)+(j)*imax)  ]*u21*u21 -
                         xyz[3*((imax-3)+(j)*imax)  ]*du*du;
                rj[1]  = xyz[3*((imax-2)+(j)*imax)+1]*u20*u20 -
                         xyz[3*((imax-1)+(j)*imax)+1]*u21*u21 -
                         xyz[3*((imax-3)+(j)*imax)+1]*du*du;
                rj[2]  = xyz[3*((imax-2)+(j)*imax)+2]*u20*u20 -
                         xyz[3*((imax-1)+(j)*imax)+2]*u21*u21 -
                         xyz[3*((imax-3)+(j)*imax)+2]*du*du;
                rj[0] /= 2.0*u21*du;
                rj[1] /= 2.0*u21*du;
                rj[2] /= 2.0*u21*du;
                dx     = rj[0] - xyz[3*((imax-1)+(j)*imax)  ] + 0.5*u20*eval[3];
                dy     = rj[1] - xyz[3*((imax-1)+(j)*imax)+1] + 0.5*u20*eval[4];
                dz     = rj[2] - xyz[3*((imax-1)+(j)*imax)+2] + 0.5*u20*eval[5];
            }
            /* match input tangent */
            if (easT != NULL)
                if ((easT[3*j  ] != 0.0) || (easT[3*j+1] != 0.0) ||
                    (easT[3*j+2] != 0.0)) {
                    dx = (easT[3*j  ] + eval[3])*du;
                    dy = (easT[3*j+1] + eval[4])*du;
                    dz = (easT[3*j+2] + eval[5])*du;
                }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(j+1)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax)+(j+1)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax)+(j+1)*(imax+2))+2] += RELAX * dz;
        }

        /* convergence check */
        if (value(dxyzmax) < tol) break;

    }
#ifdef DEBUG
    printf(" EG_spline2d: Convergence at iteration %d (%le)!\n", iter,
           value(dxyzmax));
#endif
    if (value(dxyzmax) >= tol)
        printf(" EGADS Warning: Not Converged (EG_spline2d)!\n");
  
    for (i = 0; i < 3*icp*jcp; i++) cpsav[i] = cp[i];
    EG_free(cp);
  
    *rdata = rvec;
    return EGADS_SUCCESS;
}


template<int N>
int
EG_spline2dAprx(int endc, int imax, int jmax, const SurrealS<N> *xyz,
                const SurrealS<N> *uknot,     const SurrealS<N> *vknot,
                const int *vdata,
                const SurrealS<N> *wesT,      const SurrealS<N> *easT,
                const SurrealS<N> *south,           SurrealS<N> *snor,
                const SurrealS<N> *north,           SurrealS<N> *nnor,
                double tol, int *header,            SurrealS<N> **rdata)
{
  return EG_spline2dAppr(endc, imax, jmax, xyz, uknot, vknot, vdata, wesT, easT,
                         south, snor, north, nnor, tol, header, rdata);
}

// Create explicit instantiations of the function
template DllExport int EG_spline2dAprx<1>(int endc, int imax, int jmax,
                                          const SurrealS<1> *, const SurrealS<1> *,
                                          const SurrealS<1> *, const int *,
                                          const SurrealS<1> *, const SurrealS<1> *,
                                          const SurrealS<1> *,       SurrealS<1> *,
                                          const SurrealS<1> *,       SurrealS<1> *,
                                          double tol, int *header,   SurrealS<1> **);


extern "C"
int
EG_spline2dAppx(egObject *context,  int endc,   /*@null@*/ const double *uknot,
                /*@null@*/ const double *vknot, /*@null@*/ const int    *vdata,
                /*@null@*/ const double *wesT,  /*@null@*/ const double *easT,
                /*@null@*/ const double *south, /*@null@*/       double *snor,
                /*@null@*/ const double *north, /*@null@*/       double *nnor,
                int imax, int jmax, const double *xyz,           double tol,
                egObject **esurf)
{
  int    stat, header[7];
  double *rdata;
  
  stat = EG_spline2dAppr(endc, imax, jmax, xyz, uknot, vknot, vdata, wesT, easT,
                         south, snor, north, nnor, tol, header, &rdata);
  if (stat != EGADS_SUCCESS) return stat;
  
  stat = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rdata, esurf);
  EG_free(rdata);

  return stat;
}


extern "C"
int
EG_spline2dRaw(int endc, /*@null@*/ const double **drnd, int imax, int jmax,
               const double *xyz, double tol, int *header, double **rdata)
{
  int          i, j, k, stat;
  double       dist, dy, *norT, *souT, *easT, *wesT;
  double       x0[3], x1[3], nnor[3], snor[3], enor[3], wnor[3], *rot;
  const double *north, *south, *east, *west;

  if ((imax < 2) || (jmax < 2)) return EGADS_DEGEN;
  if ((endc < 0) || (endc > 2)) return EGADS_RANGERR;
  
  /* check for degenerate sides */
  north = south = east = west = NULL;
  norT  = souT  = easT = wesT = NULL;
  i  = 0;
  dy = 0.0;
  for (j = 1; j < jmax; j++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Imin (west) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[0] != NULL) {
        west  = drnd[0];
        x0[0] = drnd[0][1];
        x0[1] = drnd[0][2];
        x0[2] = drnd[0][3];
        x1[0] = drnd[0][5];
        x1[1] = drnd[0][6];
        x1[2] = drnd[0][7];
        CROSS(wnor, x0, x1);
        dist  = DOT(wnor, wnor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol) ||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          printf(" EGADS Error: BAD West Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        wnor[0] *= dist;
        wnor[1] *= dist;
        wnor[2] *= dist;
        wesT     = wnor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               wnor[0], wnor[1], wnor[2]);
#endif
      }
  }
  j  = 0;
  dy = 0.0;
  for (i = 1; i < imax; i++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Jmin (south) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[1] != NULL) {
        south = drnd[1];
        x0[0] = drnd[1][1];
        x0[1] = drnd[1][2];
        x0[2] = drnd[1][3];
        x1[0] = drnd[1][5];
        x1[1] = drnd[1][6];
        x1[2] = drnd[1][7];
        CROSS(snor, x0, x1);
        dist  = DOT(snor, snor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol) ||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          printf(" EGADS Error: BAD South Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        snor[0] *= dist;
        snor[1] *= dist;
        snor[2] *= dist;
        souT     = snor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               snor[0], snor[1], snor[2]);
#endif
      }
  }
  i  = imax-1;
  dy = 0.0;
  for (j = 1; j < jmax; j++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Imax (east) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[2] != NULL) {
        east  = drnd[2];
        x0[0] = drnd[2][1];
        x0[1] = drnd[2][2];
        x0[2] = drnd[2][3];
        x1[0] = drnd[2][5];
        x1[1] = drnd[2][6];
        x1[2] = drnd[2][7];
        CROSS(enor, x0, x1);
        dist  = DOT(enor, enor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol)||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          printf(" EGADS Error: BAD East Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        enor[0] *= dist;
        enor[1] *= dist;
        enor[2] *= dist;
        easT     = enor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               enor[0], enor[1], enor[2]);
#endif
      }
  }
  j  = jmax-1;
  dy = 0.0;
  for (i = 1; i < imax; i++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  spline2d: Jmax (north) degenerate!\n");
#endif
    if (drnd != NULL)
      if (drnd[3] != NULL) {
        north = drnd[3];
        x0[0] = drnd[3][1];
        x0[1] = drnd[3][2];
        x0[2] = drnd[3][3];
        x1[0] = drnd[3][5];
        x1[1] = drnd[3][6];
        x1[2] = drnd[3][7];
        CROSS(nnor, x0, x1);
        dist  = DOT(nnor, nnor);
        if ((dist == 0.0) || (DOT(x0, x1) > tol)||
            (fabs(1.0-DOT(x0, x0)) > tol) ||
            (fabs(1.0-DOT(x1, x1)) > tol)) {
          printf(" EGADS Error: BAD North Axes (EG_spline2d)!\n");
          return EGADS_NOTORTHO;
        }
        dist   = 1.0/sqrt(dist);
        nnor[0] *= dist;
        nnor[1] *= dist;
        nnor[2] *= dist;
        norT     = nnor;
#ifdef DEBUG
        printf("            with normal = %lf %lf %lf!\n",
               nnor[0], nnor[1], nnor[2]);
#endif
      }
  }
  if ((north != NULL) && ((east != NULL) || (west != NULL))) return EGADS_DEGEN;
  if ((south != NULL) && ((east != NULL) || (west != NULL))) return EGADS_DEGEN;
  
  /* call approx as is */
  if ((east == NULL) && (west == NULL))
    return EG_spline2dAppr(endc, imax, jmax, xyz, (double *) NULL,
                           (double *) NULL, NULL, west, east,
                           south, souT, north, norT, tol, header, rdata);
  
  /* rotate to get special treatment as north/south */
  rot = (double *) EG_alloc(3*imax*jmax*sizeof(double));
  if (rot == NULL) return EGADS_MALLOC;
  
  for (k = i = 0; i < imax; i++)
    for (j = jmax-1; j >= 0; j--, k++) {
      rot[3*k  ] = xyz[3*(i+j*imax)  ];
      rot[3*k+1] = xyz[3*(i+j*imax)+1];
      rot[3*k+2] = xyz[3*(i+j*imax)+2];
    }
#ifdef DEBUG
  printf(" spline2d: rotate -- south=west, north=east!\n");
#endif
  stat = EG_spline2dAppr(endc, jmax, imax, rot, (double *) NULL,
                         (double *) NULL, NULL, south, north,
                         west, wesT, east, easT, tol, header, rdata);
  
  EG_free(rot);
  return stat;
}
