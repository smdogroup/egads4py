/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             BSpline Fit Functions
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


/* NOTE: Most all of this has been deactivated!
 *       The code used can now be found in egadsSplineFit.cpp
 *       The only function still not "ifdef"ed out (OLDSPLINE) is EG_spline2dFit
 *           which is used for mapping BSpline surfaces when knots have changed
 *           and is currently not called due to and "ifdef" of MAPBSPLINE
 *           (which is not active) in egadsTopo.cpp & egadsGeom.cpp.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"


#define NITER    10000
#define RELAX    0.15

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define NOSECROSS
#define NEWNOSE


extern "C" int  EG_makeGeometry( egObject *context, int oclass, int mtype,
                                 /*@null@*/ egObject *refGeo, const int *ivec,
                                 const double *rvec, egObject **geom );

#ifdef OLDSPLINE
extern "C" int  EG_spline1d( egObject *context, int endc, int imax,
                             const double *xyz, double tol, egObject **ecrv );
extern "C" int  EG_spline2dAppx( egObject *context, int endc,
                                 /*@null@*/ const double *uknot,
                                 /*@null@*/ const double *vknot,
                                 /*@null@*/ /*unused@*/ const int *vdata,
                                 /*@null@*/ const double *wesT,
                                 /*@null@*/ const double *easT,
                                 /*@null@*/ const double *south,
                                 /*@null@*/       double *snor,
                                 /*@null@*/ const double *north,
                                 /*@null@*/       double *nnor,
                                 int imax, int jmax, const double *xyz,
                                 double tol, egObject **esurf );
#endif
extern "C" int  EG_spline2dFit( egObject *context, const double *crosT,
                                int imax, const double *uknot,
                                const double *souT, const double *norT,
                                int jmax, const double *vknot,
                                const double *wesT, const double *easT,
                                const double *xyz, double tol,
                                egObject **esurf );


#ifdef OLDSPLINE
/*
 ************************************************************************
 *                                                                      *
 *   EG_spline1d - create 1d cubic spline from input data               *
 *                                                                      *
 ************************************************************************
 */

int EG_spline1d(egObject *context, int endx, int imaxx, const double *xyz,
                double tol, egObject **ecurv)
{
    int    i, kk, iknot, icp, iter, header[4], status, endc, imax;
    double du, dx, dy, dz, dxyzmax, rj[3], u21, u20;
    double *rvec, *cp, *knots;
    gp_Pnt P0;
    gp_Vec V1, V2;

    *ecurv = NULL;
    endc   = endx;
    imax   = imaxx;
    if (imax    < 0)                   imax = -imax;
    if (context == NULL)               return EGADS_NULLOBJ;
    if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
    if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
    if (imax < 2)                      return EGADS_DEGEN;

    if ((endx < 0) || (endx > 2))      return EGADS_RANGERR;
    if ((imax == 2) && (endc == 2))    endc = 1;

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

    icp   = imax + 2;
    iknot = imax + 6;
    rvec  = (double *) EG_alloc((iknot+3*icp)*sizeof(double));
    if (rvec == NULL) return EGADS_MALLOC;
    knots =  rvec;
    cp    = &rvec[iknot];

    /* create spline curve */
    header[0] = 0;
    header[1] = 3;
    header[2] = icp;
    header[3] = iknot;

    /* knots */
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

    /* initial control point */
    kk       = 0;
    cp[kk++] = xyz[0];
    cp[kk++] = xyz[1];
    cp[kk++] = xyz[2];

    /* initial interior control point (for slope) */
    cp[kk++] = (3 * xyz[0] + xyz[3]) / 4;
    cp[kk++] = (3 * xyz[1] + xyz[4]) / 4;
    cp[kk++] = (3 * xyz[2] + xyz[5]) / 4;

    /* interior control points */
    for (i = 1; i < imax-1; i++) {
        cp[kk++] = xyz[3*i  ];
        cp[kk++] = xyz[3*i+1];
        cp[kk++] = xyz[3*i+2];
    }

    /* penultimate interior control point (for slope) */
    cp[kk++] = (3 * xyz[3*(imax-1)  ] + xyz[3*(imax-2)  ]) / 4;
    cp[kk++] = (3 * xyz[3*(imax-1)+1] + xyz[3*(imax-2)+1]) / 4;
    cp[kk++] = (3 * xyz[3*(imax-1)+2] + xyz[3*(imax-2)+2]) / 4;

    /* final control point */
    cp[kk++] = xyz[3*(imax-1)  ];
    cp[kk++] = xyz[3*(imax-1)+1];
    cp[kk++] = xyz[3*(imax-1)+2];

    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rvec,
                             ecurv);
    if (status != EGADS_SUCCESS) {
        EG_free(rvec);
        return status;
    }

    /* iterate to have knot evaluations match data points */
    for (iter = 0; iter < NITER; iter++) {
        egadsCurve         *pcurve = (egadsCurve *) (*ecurv)->blind;
        Handle(Geom_Curve) hCurve  = pcurve->handle;
        dxyzmax = 0.0;

        /* condition at beginning */
        hCurve->D2(knots[3], P0, V1, V2);
        du = knots[4] - knots[3];
        if (endc == 0) {
            /* natural end */
            dx = du * du * V2.X();
            dy = du * du * V2.Y();
            dz = du * du * V2.Z();
        } else if (endc == 1) {
            /* FD slope */
            dx = xyz[3] - xyz[0] - du * V1.X();
            dy = xyz[4] - xyz[1] - du * V1.Y();
            dz = xyz[5] - xyz[2] - du * V1.Z();
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
            dx     = rj[0] - xyz[0] - 0.5 * u20 * V1.X();
            dy     = rj[1] - xyz[1] - 0.5 * u20 * V1.Y();
            dz     = rj[2] - xyz[2] - 0.5 * u20 * V1.Z();
        }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3] += RELAX * dx;
        cp[4] += RELAX * dy;
        cp[5] += RELAX * dz;

        /* match interior spline points */
        for (i = 1; i < imax-1; i++) {
            hCurve->D0(knots[i+3], P0);
            dx = xyz[3*i  ] - P0.X();
            dy = xyz[3*i+1] - P0.Y();
            dz = xyz[3*i+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*i+3] += dx;
            cp[3*i+4] += dy;
            cp[3*i+5] += dz;
        }
      
        /* condition at end */
        hCurve->D2(knots[imax+2], P0, V1, V2);
        du = knots[imax+2] - knots[imax+1];
        if (endc == 0) {
            /* natural end */
            dx = du * du * V2.X();
            dy = du * du * V2.Y();
            dz = du * du * V2.Z();
        } else if (endc == 1) {
            /* FD slope */
            dx = xyz[3*(imax-2)  ] - xyz[3*(imax-1)  ] + du * V1.X();
            dy = xyz[3*(imax-2)+1] - xyz[3*(imax-1)+1] + du * V1.Y();
            dz = xyz[3*(imax-2)+2] - xyz[3*(imax-1)+2] + du * V1.Z();
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
            dx     = rj[0] - xyz[3*(imax-1)  ] + 0.5 * u20 * V1.X();
            dy     = rj[1] - xyz[3*(imax-1)+1] + 0.5 * u20 * V1.Y();
            dz     = rj[2] - xyz[3*(imax-1)+2] + 0.5 * u20 * V1.Z();
        }
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*imax  ] += RELAX * dx;
        cp[3*imax+1] += RELAX * dy;
        cp[3*imax+2] += RELAX * dz;

        /* convergence check */
        if (dxyzmax < tol) break;

        /* make the new curve (after deleting old one) */
        status = EG_deleteObject(*ecurv);
        if (status != EGADS_SUCCESS) {
            *ecurv = NULL;
            EG_free(rvec);
            return status;
        }
        status = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, rvec,
                                 ecurv);
        if (status != EGADS_SUCCESS) {
            EG_free(rvec);
            return status;
        }
    }
    if (dxyzmax >= tol)
        printf(" EGADS Warning: Not Converged (EG_spline1d)!\n");

    EG_free(rvec);

    return EGADS_SUCCESS;
}


static double
EG_cubicBasisFun(int nKnots, double *knots, int i, double u)
{
    int    j, k;
    double N[4], saved, Uleft, Uright, temp;
  
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


static void
EG_NzeroBSp(int nknots, double *knots, double t, int *nbasis, double *basis)
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


static double
EG_getEllRad(const double *drnd, double *nell)
{
#ifdef NEWNOSE
    double csth, snth;
  
    csth = drnd[0]*(drnd[1]*nell[0] + drnd[2]*nell[1] + drnd[3]*nell[2]);
    snth = drnd[4]*(drnd[5]*nell[0] + drnd[6]*nell[1] + drnd[7]*nell[2]);
    return sqrt(csth*csth + snth*snth);
#else
    double norm, x1, x2, thet, csth, snth;

    norm = sqrt(DOT(nell, nell));
    x1   = drnd[1]*nell[0]/norm + drnd[2]*nell[1]/norm + drnd[3]*nell[2]/norm;
    x2   = drnd[5]*nell[0]/norm + drnd[6]*nell[1]/norm + drnd[7]*nell[2]/norm;
    thet = atan2(x2, x1);
    csth = cos(thet);
    snth = sin(thet);
    return sqrt(drnd[0]*csth*drnd[0]*csth + drnd[4]*snth*drnd[4]*snth);
#endif
}


/*
 ************************************************************************
 *                                                                      *
 *   EG_spline2dAppx - Approximate 2d cubic spline from input data      *
 *                                                                      *
 ************************************************************************
 */

int
EG_spline2dAppx(egObject *context, int endc,
                /*@null@*/ const double *uknot, /*@null@*/ const double *vknot,
                /*@null@*/ /*unused@*/ const int *vdata,
                /*@null@*/ const double *wesT,  /*@null@*/ const double *easT,
                /*@null@*/ const double *south, /*@null@*/       double *snor,
                /*@null@*/ const double *north, /*@null@*/       double *nnor,
                int imax, int jmax, const double *xyz, double tol,
                egObject **esurf)
{
    int    i, j, iknot, jknot, icp, jcp, iter, header[7], status;
    int    endi, endj, jj, kk;
    double r, t, du, dv, dx, dy, dz, dist, mmu, con, con2, dxyzmax, normnell;
    double dD, eE, F, G, rj[3], u21, u20, norm[3], nell[3], x0[3], x1[3], T[3][3];
    double basis[4], *rvec, *knotu, *knotv, *cp, *cpsav;
    gp_Pnt P0;
    gp_Vec V1, V2, U1, U2, UV;
#ifdef NEWNOSE
    int    ms,    mn;
    double ns[3], nn[3], rs[3][3], rn[3][3], thet, q0, q1, q2, q3, x2[3];
  
    ms       = 1;
    rs[0][0] = rs[1][1] = rs[2][2] = 1.0;
    rs[0][1] = rs[0][2] = rs[1][0] = rs[1][2] = rs[2][0] = rs[2][1] = 0.0;
  
    mn       = 1;
    rn[0][0] = rn[1][1] = rn[2][2] = 1.0;
    rn[0][1] = rn[0][2] = rn[1][0] = rn[1][2] = rn[2][0] = rn[2][1] = 0.0;
#endif

    *esurf = NULL;
    endi   = endj = endc;
    if (imax == 2) endi = 1;
    if (jmax == 2) endj = 1;

    cpsav = NULL;
    icp   = imax + 2;
    iknot = imax + 6;
    jcp   = jmax + 2;
    jknot = jmax + 6;
    rvec  = (double *) EG_alloc((iknot+jknot+3*icp*jcp)*sizeof(double));
    if (rvec == NULL) return EGADS_MALLOC;
    knotu =  rvec;
    knotv = &rvec[iknot      ];
    cp    = &rvec[iknot+jknot];

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
            if (dy == 0.0) {
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
/*  for (i = 0; i < kk; i++) printf(" %d:  %lf\n", i, knotu[i]);  */

    /* knots in j-direction */
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
            if (dy == 0.0) {
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
        for (iter = j = 1; j < jmax; j++)
            if (knotv[kk+j] <= knotv[kk+j-1]) iter++;
        if (iter == 1) {
            for (j = 0; j < jmax; j++) knotv[kk++] /= dz;
        } else {
            /* equally spaced */
            for (j = 0; j < jmax; j++) {
                dy          = j;
                knotv[kk++] = dy/(jmax-1);
            }
        }
    } else {
        /* knots in v set by input */
        for (j = 0; j < jmax; j++) knotv[kk++] = vknot[j];
    }
    knotv[kk++] = 1.0;
    knotv[kk++] = 1.0;
    knotv[kk++] = 1.0;
/*  for (i = 0; i < jknot; i++) printf(" %d:  %lf\n", i, knotv[i]);  */
  
    if ((north != NULL) || (south != NULL)) {
        cpsav = (double *) EG_alloc((3*icp*jcp)*sizeof(double));
        if (cpsav == NULL) {
            EG_free(rvec);
            return EGADS_MALLOC;
        }
    }

/*
  {
    static int count = 0;
    char       filename[32];
    int        k;
    FILE       *fp;
    
    sprintf(filename, "X%d.dat", count);
    fp = fopen(filename, "w");
    for (j = 0; j < jmax; j++) {
      for (i = 0; i < imax; i++)
        fprintf(fp, " %lf", xyz[3*(j*imax+i)  ]);
      fprintf(fp, "\n");
    }
    fclose(fp);
    sprintf(filename, "Y%d.dat", count);
    fp = fopen(filename, "w");
    for (j = 0; j < jmax; j++) {
      for (i = 0; i < imax; i++)
        fprintf(fp, " %lf", xyz[3*(j*imax+i)+1]);
      fprintf(fp, "\n");
    }
    fclose(fp);
    sprintf(filename, "Z%d.dat", count);
    fp = fopen(filename, "w");
    for (j = 0; j < jmax; j++) {
      for (i = 0; i < imax; i++)
        fprintf(fp, " %lf", xyz[3*(j*imax+i)+2]);
      fprintf(fp, "\n");
    }
    fclose(fp);
    sprintf(filename, "param%d.dat", count);
    fp = fopen(filename, "w");
    for (i = header[3]+3; i < header[3]+header[6]-3; i++)
      fprintf(fp, " %lf", rvec[i]);
    fprintf(fp, "\n");
    for (k = 3; k < header[3]-3; k++)
      fprintf(fp, " %lf", rvec[k]);
    fprintf(fp, "\n");
    if (south != NULL) {
      fprintf(fp," %lf %lf %lf %lf\n", south[0], south[1], south[2], south[3]);
      fprintf(fp," %lf %lf %lf %lf\n", south[4], south[5], south[6], south[7]);
    }
    if (north != NULL) {
      fprintf(fp," %lf %lf %lf %lf\n", north[0], north[1], north[2], north[3]);
      fprintf(fp," %lf %lf %lf %lf\n", north[4], north[5], north[6], north[7]);
    }
    fclose(fp);
    if ((wesT != NULL) && (easT != NULL)) {
      sprintf(filename, "tangent%d.dat", count);
      fp = fopen(filename, "w");
      for (j = 0; j < jmax; j++) {
        fprintf(fp, " %lf %lf %lf   %lf %lf %lf\n",
                 wesT[3*j  ],  wesT[3*j+1],  wesT[3*j+2],
                -easT[3*j  ], -easT[3*j+1], -easT[3*j+2]);
      }
      fclose(fp);
    }
    count++;
  }
 */

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
#ifdef NEWNOSE
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
               ns[0], ns[1], ns[2], ms, thet);
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
                   x0[0], x0[1], x0[2], dist);
#endif
            dist   = sin(thet/2.0);
            q0     = cos(thet/2.0);
            q1     = x0[0]*dist;
            q2     = x0[1]*dist;
            q3     = x0[2]*dist;
#ifdef DEBUG
            printf(" Quaterion = %lf %lf %lf %lf\n", q0, q1, q2, q3);
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
        printf(" Rotation Matrix = %lf %lf %lf\n", rs[0][0], rs[0][1], rs[0][2]);
        printf("                   %lf %lf %lf\n", rs[1][0], rs[1][1], rs[1][2]);
        printf("                   %lf %lf %lf\n", rs[2][0], rs[2][1], rs[2][2]);
        printf("\n");
#endif
#endif
        mmu = 2.0*(knotv[4] - knotv[3])/(3.0*((knotv[5] - knotv[3])));
        for (i = -1; i < imax+1; i++) {
#ifdef NEWNOSE
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
            t        = sqrt(r*mmu*sqrt(DOT(x2, x2)));
#else
            x0[0] = cp[3*((i+1)+(2)*(imax+2))  ] - cp[3*((i+1)+(0)*(imax+2))  ];
            x0[1] = cp[3*((i+1)+(2)*(imax+2))+1] - cp[3*((i+1)+(0)*(imax+2))+1];
            x0[2] = cp[3*((i+1)+(2)*(imax+2))+2] - cp[3*((i+1)+(0)*(imax+2))+2];
            x1[0] = cp[3*((i+1)+(3)*(imax+2))  ] - cp[3*((i+1)+(0)*(imax+2))  ];
            x1[1] = cp[3*((i+1)+(3)*(imax+2))+1] - cp[3*((i+1)+(0)*(imax+2))+1];
            x1[2] = cp[3*((i+1)+(3)*(imax+2))+2] - cp[3*((i+1)+(0)*(imax+2))+2];
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, snor);
                nell[0] = x0[0] - dist*snor[0];
                nell[1] = x0[1] - dist*snor[1];
                nell[2] = x0[2] - dist*snor[2];
            } else {
                CROSS(nell, norm, snor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(south, nell);
            normnell = sqrt(DOT(nell, nell));
            CROSS(x1, nell, x0);
            dist = sqrt(DOT(x1, x1));
            t    = sqrt(r*mmu*dist/(normnell*normnell*normnell));
#endif
#ifdef DEBUG
            printf(" south init %d: %lf %lf %lf\n", i+1,
                   t*nell[0], t*nell[1], t*nell[2]);
#endif
            cp[3*((i+1)+(1)*(imax+2))  ] = cp[3*((i+1)+(0)*(imax+2))  ] +
                                           t*nell[0];
            cp[3*((i+1)+(1)*(imax+2))+1] = cp[3*((i+1)+(0)*(imax+2))+1] +
                                           t*nell[1];
            cp[3*((i+1)+(1)*(imax+2))+2] = cp[3*((i+1)+(0)*(imax+2))+2] +
                                           t*nell[2];
        }
    }
  
    if (north != NULL) {
#ifdef NEWNOSE
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
               nn[0], nn[1], nn[2], mn, thet);
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
                   x0[0], x0[1], x0[2], dist);
#endif
            dist   = sin(thet/2.0);
            q0     = cos(thet/2.0);
            q1     = x0[0]*dist;
            q2     = x0[1]*dist;
            q3     = x0[2]*dist;
#ifdef DEBUG
            printf(" Quaterion = %lf %lf %lf %lf\n", q0, q1, q2, q3);
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
        printf(" Rotation Matrix = %lf %lf %lf\n", rn[0][0], rn[0][1], rn[0][2]);
        printf("                   %lf %lf %lf\n", rn[1][0], rn[1][1], rn[1][2]);
        printf("                   %lf %lf %lf\n", rn[2][0], rn[2][1], rn[2][2]);
        printf("\n");
#endif
#endif
        mmu =  2.0*( knotv[jmax+1] - knotv[jmax+2])/
              (3.0*((knotv[jmax]   - knotv[jmax+2])));
        for (i = -1; i < imax+1; i++) {
#ifdef NEWNOSE
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
            t        = sqrt(r*mmu*sqrt(DOT(x2, x2)));
#else
            x0[0] = cp[3*((i+1)+(jmax-1)*(imax+2))  ] -
                    cp[3*((i+1)+(jmax+1)*(imax+2))  ];
            x0[1] = cp[3*((i+1)+(jmax-1)*(imax+2))+1] -
                    cp[3*((i+1)+(jmax+1)*(imax+2))+1];
            x0[2] = cp[3*((i+1)+(jmax-1)*(imax+2))+2] -
                    cp[3*((i+1)+(jmax+1)*(imax+2))+2];
            x1[0] = cp[3*((i+1)+(jmax-2)*(imax+2))  ] -
                    cp[3*((i+1)+(jmax+1)*(imax+2))  ];
            x1[1] = cp[3*((i+1)+(jmax-2)*(imax+2))+1] -
                    cp[3*((i+1)+(jmax+1)*(imax+2))+1];
            x1[2] = cp[3*((i+1)+(jmax-2)*(imax+2))+2] -
                    cp[3*((i+1)+(jmax+1)*(imax+2))+2];
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, nnor);
                nell[0] = x0[0] - dist*nnor[0];
                nell[1] = x0[1] - dist*nnor[1];
                nell[2] = x0[2] - dist*nnor[2];
            } else {
                CROSS(nell, norm, nnor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(north, nell);
            normnell = sqrt(DOT(nell, nell));
            CROSS(x1, nell, x0);
            dist = sqrt(DOT(x1, x1));
            t    = sqrt(r*mmu*dist/(normnell*normnell*normnell));
#endif
#ifdef DEBUG
            printf(" north init %d: %lf %lf %lf\n", i+1,
                   t*nell[0], t*nell[1], t*nell[2]);
#endif
            cp[3*((i+1)+(jmax)*(imax+2))  ] = cp[3*((i+1)+(jmax+1)*(imax+2))  ] +
                                              t*nell[0];
            cp[3*((i+1)+(jmax)*(imax+2))+1] = cp[3*((i+1)+(jmax+1)*(imax+2))+1] +
                                              t*nell[1];
            cp[3*((i+1)+(jmax)*(imax+2))+2] = cp[3*((i+1)+(jmax+1)*(imax+2))+2] +
                                              t*nell[2];
        }
    }

    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rvec,
                             esurf);
    if (status != EGADS_SUCCESS) {
        if (cpsav != NULL) EG_free(cpsav);
        EG_free(rvec);
        return status;
    }

    /* iterate to have knot evaluations match data points */
    for (iter = 0; iter < NITER; iter++) {
      
        egadsSurface *psurf        = (egadsSurface *) (*esurf)->blind;
        Handle(Geom_Surface) hSurf = psurf->handle;
        dxyzmax = 0.0;
        if (cpsav != NULL)
            for (i = 0; i < 3*icp*jcp; i++) cpsav[i] = cp[i];

        /* match interior spline points */
        for (j = 1; j < jmax-1; j++) {
            for (i = 1; i < imax-1; i++) {
                hSurf->D0(knotu[i+3], knotv[j+3], P0);
                dx = xyz[3*((i)+(j)*imax)  ] - P0.X();
                dy = xyz[3*((i)+(j)*imax)+1] - P0.Y();
                dz = xyz[3*((i)+(j)*imax)+2] - P0.Z();
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(j+1)*(imax+2))  ] += dx;
                cp[3*((i+1)+(j+1)*(imax+2))+1] += dy;
                cp[3*((i+1)+(j+1)*(imax+2))+2] += dz;
            }
        }

        /* point A */
        hSurf->D2(knotu[3], knotv[3], P0, U1, V1, U2, V2, UV);
        du = knotu[4] - knotu[3];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * U2.X();
            dy = du * du * U2.Y();
            dz = du * du * U2.Z();
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((1)+(0)*imax)  ] - xyz[3*((0)+(0)*imax)  ] - du*U1.X();
            dy = xyz[3*((1)+(0)*imax)+1] - xyz[3*((0)+(0)*imax)+1] - du*U1.Y();
            dz = xyz[3*((1)+(0)*imax)+2] - xyz[3*((0)+(0)*imax)+2] - du*U1.Z();
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
            dx     = rj[0] - xyz[3*((0)+(0)*imax)  ] - 0.5 * u20 * U1.X();
            dy     = rj[1] - xyz[3*((0)+(0)*imax)+1] - 0.5 * u20 * U1.Y();
            dz     = rj[2] - xyz[3*((0)+(0)*imax)+2] - 0.5 * u20 * U1.Z();
        }
        if (wesT != NULL) {
            /* match input tangent */
            dx = (wesT[0] - U1.X())*du;
            dy = (wesT[1] - U1.Y())*du;
            dz = (wesT[2] - U1.Z())*du;
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
            dx = dv * dv * V2.X();
            dy = dv * dv * V2.Y();
            dz = dv * dv * V2.Z();
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((0)+(1)*imax)  ] - xyz[3*((0)+(0)*imax)  ] - dv*V1.X();
            dy = xyz[3*((0)+(1)*imax)+1] - xyz[3*((0)+(0)*imax)+1] - dv*V1.Y();
            dz = xyz[3*((0)+(1)*imax)+2] - xyz[3*((0)+(0)*imax)+2] - dv*V1.Z();
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
            dx     = rj[0] - xyz[3*((0)+(0)*imax)  ] - 0.5 * u20 * V1.X();
            dy     = rj[1] - xyz[3*((0)+(0)*imax)+1] - 0.5 * u20 * V1.Y();
            dz     = rj[2] - xyz[3*((0)+(0)*imax)+2] - 0.5 * u20 * V1.Z();
        }
        if ((south != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
#ifdef NEWNOSE
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
#else
            x0[0] = xyz[3*((0)+(1)*(imax))  ] - xyz[3*((0)+(0)*(imax))  ];
            x0[1] = xyz[3*((0)+(1)*(imax))+1] - xyz[3*((0)+(0)*(imax))+1];
            x0[2] = xyz[3*((0)+(1)*(imax))+2] - xyz[3*((0)+(0)*(imax))+2];
            x1[0] = xyz[3*((0)+(2)*(imax))  ] - xyz[3*((0)+(0)*(imax))  ];
            x1[1] = xyz[3*((0)+(2)*(imax))+1] - xyz[3*((0)+(0)*(imax))+1];
            x1[2] = xyz[3*((0)+(2)*(imax))+2] - xyz[3*((0)+(0)*(imax))+2];
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
              dist    = DOT(x0, snor);
              nell[0] = x0[0] - dist*snor[0];
              nell[1] = x0[1] - dist*snor[1];
              nell[2] = x0[2] - dist*snor[2];
            } else {
              CROSS(nell, norm, snor);
            }
            if (DOT(x0, nell) < 0.0) {
              nell[0] = -nell[0];
              nell[1] = -nell[1];
              nell[2] = -nell[2];
            }
            r        = EG_getEllRad(south, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
            x1[0]    = cpsav[3*((0)+(2)*(imax+2))  ] -
                       cpsav[3*((0)+(0)*(imax+2))  ];
            x1[1]    = cpsav[3*((0)+(2)*(imax+2))+1] -
                       cpsav[3*((0)+(0)*(imax+2))+1];
            x1[2]    = cpsav[3*((0)+(2)*(imax+2))+2] -
                       cpsav[3*((0)+(0)*(imax+2))+2];
            t        = sqrt(r*con2*fabs(DOT(x1,snor))/con)/normnell;
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", 0,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((0)+(0)*imax)  ] + t*nell[0] -
                 cp[3*((0)+(1)*(imax+2))  ];
            dy = xyz[3*((0)+(0)*imax)+1] + t*nell[1] -
                 cp[3*((0)+(1)*(imax+2))+1];
            dz = xyz[3*((0)+(0)*imax)+2] + t*nell[2] -
                 cp[3*((0)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(1)*(imax+2))  ] += dx;
            cp[3*((0)+(1)*(imax+2))+1] += dy;
            cp[3*((0)+(1)*(imax+2))+2] += dz;
        } else {
            if (snor != NULL) {
                dx = (snor[0] - V1.X())*dv;
                dy = (snor[1] - V1.Y())*dv;
                dz = (snor[2] - V1.Z())*dv;
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
            dx = du * dv * UV.X();
            dy = du * dv * UV.Y();
            dz = du * dv * UV.Z();
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((1)+(1)*imax)  ] - xyz[3*((0)+(1)*imax)  ]) -
                 (xyz[3*((1)+(0)*imax)  ] - xyz[3*((0)+(0)*imax)  ]) +
                 du * dv * UV.X();
            dy = (xyz[3*((1)+(1)*imax)+1] - xyz[3*((0)+(1)*imax)+1]) -
                 (xyz[3*((1)+(0)*imax)+1] - xyz[3*((0)+(0)*imax)+1]) +
                 du * dv * UV.Y();
            dz = (xyz[3*((1)+(1)*imax)+2] - xyz[3*((0)+(1)*imax)+2]) -
                 (xyz[3*((1)+(0)*imax)+2] - xyz[3*((0)+(0)*imax)+2]) +
                 du * dv * UV.Z();
        } else {
            /* quadratic fit */
            u20 = knotu[5] - knotu[3];
            u21 = knotu[5] - knotu[4];
            for (j = 0; j < 3; j++) {
                if (wesT != NULL) {
                    T[j][0] = wesT[3*j  ]*du;
                    T[j][1] = wesT[3*j+1]*du;
                    T[j][2] = wesT[3*j+2]*du;
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
                T[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((0)+(j)*imax)  ];
                T[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((0)+(j)*imax)+1];
                T[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((0)+(j)*imax)+2];
/*              if (iter == 0) {
                    printf(" t%d = %lf %lf %lf   %lf %lf %lf\n", j, T[j][0],
                        T[j][1], T[j][2], wesT[3*j], wesT[3*j+1], wesT[3*j+2]);
                    printf("    = %lf %lf %lf\n",
                           xyz[3*((1)+(j)*imax)  ] - xyz[3*((0)+(j)*imax)  ],
                           xyz[3*((1)+(j)*imax)+1] - xyz[3*((0)+(j)*imax)+1],
                           xyz[3*((1)+(j)*imax)+2] - xyz[3*((0)+(j)*imax)+2]);
                }  */
            }
            u20    = knotv[5] - knotv[3];
            u21    = knotv[5] - knotv[4];
            rj[0]  = T[1][0]*u20*u20 - T[0][0]*u21*u21 - T[2][0]*dv*dv;
            rj[1]  = T[1][1]*u20*u20 - T[0][1]*u21*u21 - T[2][1]*dv*dv;
            rj[2]  = T[1][2]*u20*u20 - T[0][2]*u21*u21 - T[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[5] - knotu[3]);
            if (wesT != NULL) u21 = du;
            dx     = -rj[0] - T[0][0] - 0.5*u20*u21*UV.X();
            dy     = -rj[1] - T[0][1] - 0.5*u20*u21*UV.Y();
            dz     = -rj[2] - T[0][2] - 0.5*u20*u21*UV.Z();
        }
#ifdef NOSECROSS
        if ((south != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
            dist  = 0.5*(knotu[3] + knotu[4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
#ifdef NEWNOSE
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
#else
            x0[0] = 0.5*(xyz[3*((1)+(1)*(imax))  ] +
                         xyz[3*((0)+(1)*(imax))  ]) -
                    0.5*(xyz[3*((1)+(0)*(imax))  ] +
                         xyz[3*((0)+(0)*(imax))  ]);
            x0[1] = 0.5*(xyz[3*((1)+(1)*(imax))+1] +
                         xyz[3*((0)+(1)*(imax))+1]) -
                    0.5*(xyz[3*((1)+(0)*(imax))+1] +
                         xyz[3*((0)+(0)*(imax))+1]);
            x0[2] = 0.5*(xyz[3*((1)+(1)*(imax))+2] +
                         xyz[3*((0)+(1)*(imax))+2]) -
                    0.5*(xyz[3*((1)+(0)*(imax))+2] +
                         xyz[3*((0)+(0)*(imax))+2]);
            x1[0] = 0.5*(xyz[3*((1)+(2)*(imax))  ] +
                         xyz[3*((0)+(2)*(imax))  ]) -
                    0.5*(xyz[3*((1)+(0)*(imax))  ] +
                         xyz[3*((0)+(0)*(imax))  ]);
            x1[1] = 0.5*(xyz[3*((1)+(2)*(imax))+1] +
                         xyz[3*((0)+(2)*(imax))+1]) -
                    0.5*(xyz[3*((1)+(0)*(imax))+1] +
                         xyz[3*((0)+(0)*(imax))+1]);
            x1[2] = 0.5*(xyz[3*((1)+(2)*(imax))+2] +
                         xyz[3*((0)+(2)*(imax))+2]) -
                    0.5*(xyz[3*((1)+(0)*(imax))+2] +
                         xyz[3*((0)+(0)*(imax))+2]);
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, snor);
                nell[0] = x0[0] - dist*snor[0];
                nell[1] = x0[1] - dist*snor[1];
                nell[2] = x0[2] - dist*snor[2];
            } else {
                CROSS(nell, norm, snor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(south, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
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
            t  = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", 1,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((0)+(0)*imax)  ] + t*nell[0] -
                 cp[3*((1)+(1)*(imax+2))  ];
            dy = xyz[3*((0)+(0)*imax)+1] + t*nell[1] -
                 cp[3*((1)+(1)*(imax+2))+1];
            dz = xyz[3*((0)+(0)*imax)+2] + t*nell[2] -
                 cp[3*((1)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(1)*(imax+2))  ] += dx;
            cp[3*((1)+(1)*(imax+2))+1] += dy;
            cp[3*((1)+(1)*(imax+2))+2] += dz;
        } else {
#endif
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((1)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((1)+(1)*(imax+2))+2] += RELAX * dz;
#ifdef NOSECROSS
        }
#endif
      
        /* match south & E points */
        dv   = knotv[4] - knotv[3];
        con  = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
        con2 = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
        for (i = 1; i < imax-1; i++) {
            hSurf->D2(knotu[i+3], knotv[3], P0, U1, V1, U2, V2, UV);
            dx = xyz[3*((i)+(0)*imax)  ] - P0.X();
            dy = xyz[3*((i)+(0)*imax)+1] - P0.Y();
            dz = xyz[3*((i)+(0)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(0)*(imax+2))  ] += dx;
            cp[3*((i+1)+(0)*(imax+2))+1] += dy;
            cp[3*((i+1)+(0)*(imax+2))+2] += dz;

            if (endj == 0) {
                /* d2/dv2 = 0 */
                dx = dv * dv * V2.X();
                dy = dv * dv * V2.Y();
                dz = dv * dv * V2.Z();
            } else if (endj == 1) {
                /* match FD d/dv */
                dx = xyz[3*((i)+(1)*imax)  ] -
                     xyz[3*((i)+(0)*imax)  ] - dv * V1.X();
                dy = xyz[3*((i)+(1)*imax)+1] -
                     xyz[3*((i)+(0)*imax)+1] - dv * V1.Y();
                dz = xyz[3*((i)+(1)*imax)+2] -
                     xyz[3*((i)+(0)*imax)+2] - dv * V1.Z();
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
                dx     = rj[0] - xyz[3*((i)+(0)*imax)  ] - 0.5 * u20 * V1.X();
                dy     = rj[1] - xyz[3*((i)+(0)*imax)+1] - 0.5 * u20 * V1.Y();
                dz     = rj[2] - xyz[3*((i)+(0)*imax)+2] - 0.5 * u20 * V1.Z();
            }
            if ((south != NULL) && (cpsav != NULL)) {
                EG_NzeroBSp(iknot, knotu, knotu[i+3], &kk, basis);
#ifdef NEWNOSE
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
#else
                x0[0] = xyz[3*((i)+(1)*(imax))  ] - xyz[3*((i)+(0)*(imax))  ];
                x0[1] = xyz[3*((i)+(1)*(imax))+1] - xyz[3*((i)+(0)*(imax))+1];
                x0[2] = xyz[3*((i)+(1)*(imax))+2] - xyz[3*((i)+(0)*(imax))+2];
                x1[0] = xyz[3*((i)+(2)*(imax))  ] - xyz[3*((i)+(0)*(imax))  ];
                x1[1] = xyz[3*((i)+(2)*(imax))+1] - xyz[3*((i)+(0)*(imax))+1];
                x1[2] = xyz[3*((i)+(2)*(imax))+2] - xyz[3*((i)+(0)*(imax))+2];
                CROSS(norm, x0, x1);
                dist  = sqrt(DOT(norm, norm));
                if (dist < 0.000001) {
                    dist    = DOT(x0, snor);
                    nell[0] = x0[0] - dist*snor[0];
                    nell[1] = x0[1] - dist*snor[1];
                    nell[2] = x0[2] - dist*snor[2];
                } else {
                    CROSS(nell, norm, snor);
                }
                if (DOT(x0, nell) < 0.0) {
                    nell[0] = -nell[0];
                    nell[1] = -nell[1];
                    nell[2] = -nell[2];
                }
                r        = EG_getEllRad(south, nell);
                normnell = sqrt(DOT(nell, nell));
#endif
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
                t  = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
                printf(" south %d: %lf %lf %lf (iter = %d)\n", i+1,
                       t*nell[0], t*nell[1], t*nell[2], iter);
#endif
                if (t < 0.0) t = 0.0;
                dx = xyz[3*((i)+(0)*imax)  ] + t*nell[0] -
                     cp[3*((i+1)+(1)*(imax+2))  ];
                dy = xyz[3*((i)+(0)*imax)+1] + t*nell[1] -
                     cp[3*((i+1)+(1)*(imax+2))+1];
                dz = xyz[3*((i)+(0)*imax)+2] + t*nell[2] -
                     cp[3*((i+1)+(1)*(imax+2))+2];
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(1)*(imax+2))  ] += dx;
                cp[3*((i+1)+(1)*(imax+2))+1] += dy;
                cp[3*((i+1)+(1)*(imax+2))+2] += dz;
            } else {
                if (snor != NULL) {
                    dx = (snor[0] - V1.X())*dv;
                    dy = (snor[1] - V1.Y())*dv;
                    dz = (snor[2] - V1.Z())*dv;
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
        hSurf->D2(knotu[imax+2], knotv[3], P0, U1, V1, U2, V2, UV);
        du = knotu[imax+2] - knotu[imax+1];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * U2.X();
            dy = du * du * U2.Y();
            dz = du * du * U2.Z();
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((imax-2)+(0)*imax)  ] -
                 xyz[3*((imax-1)+(0)*imax)  ] + du * U1.X();
            dy = xyz[3*((imax-2)+(0)*imax)+1] -
                 xyz[3*((imax-1)+(0)*imax)+1] + du * U1.Y();
            dz = xyz[3*((imax-2)+(0)*imax)+2] -
                 xyz[3*((imax-1)+(0)*imax)+2] + du * U1.Z();
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
            dx     = rj[0] - xyz[3*((imax-1)+(0)*imax)  ] + 0.5 * u20 * U1.X();
            dy     = rj[1] - xyz[3*((imax-1)+(0)*imax)+1] + 0.5 * u20 * U1.Y();
            dz     = rj[2] - xyz[3*((imax-1)+(0)*imax)+2] + 0.5 * u20 * U1.Z();
        }
        if (easT != NULL) {
            /* match input tangent */
            dx = (easT[0] + U1.X())*du;
            dy = (easT[1] + U1.Y())*du;
            dz = (easT[2] + U1.Z())*du;
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
            dx = dv * dv * V2.X();
            dy = dv * dv * V2.Y();
            dz = dv * dv * V2.Z();
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((imax-1)+(1)*imax)  ] -
                 xyz[3*((imax-1)+(0)*imax)  ] - dv * V1.X();
            dy = xyz[3*((imax-1)+(1)*imax)+1] -
                 xyz[3*((imax-1)+(0)*imax)+1] - dv * V1.Y();
            dz = xyz[3*((imax-1)+(1)*imax)+2] -
                 xyz[3*((imax-1)+(0)*imax)+2] - dv * V1.Z();
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
            dx     = rj[0] - xyz[3*((imax-1)+(0)*imax)  ] - 0.5 * u20 * V1.X();
            dy     = rj[1] - xyz[3*((imax-1)+(0)*imax)+1] - 0.5 * u20 * V1.Y();
            dz     = rj[2] - xyz[3*((imax-1)+(0)*imax)+2] - 0.5 * u20 * V1.Z();
        }
        if ((south != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
#ifdef NEWNOSE
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
#else
            x0[0] = xyz[3*((imax-1)+(1)*(imax))  ] -
                    xyz[3*((imax-1)+(0)*(imax))  ];
            x0[1] = xyz[3*((imax-1)+(1)*(imax))+1] -
                    xyz[3*((imax-1)+(0)*(imax))+1];
            x0[2] = xyz[3*((imax-1)+(1)*(imax))+2] -
                    xyz[3*((imax-1)+(0)*(imax))+2];
            x1[0] = xyz[3*((imax-1)+(2)*(imax))  ] -
                    xyz[3*((imax-1)+(0)*(imax))  ];
            x1[1] = xyz[3*((imax-1)+(2)*(imax))+1] -
                    xyz[3*((imax-1)+(0)*(imax))+1];
            x1[2] = xyz[3*((imax-1)+(2)*(imax))+2] -
                    xyz[3*((imax-1)+(0)*(imax))+2];
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, snor);
                nell[0] = x0[0] - dist*snor[0];
                nell[1] = x0[1] - dist*snor[1];
                nell[2] = x0[2] - dist*snor[2];
            } else {
                CROSS(nell, norm, snor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(south, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
            x1[0]    = cpsav[3*((imax+1)+(2)*(imax+2))  ] -
                       cpsav[3*((imax+1)+(0)*(imax+2))  ];
            x1[1]    = cpsav[3*((imax+1)+(2)*(imax+2))+1] -
                       cpsav[3*((imax+1)+(0)*(imax+2))+1];
            x1[2]    = cpsav[3*((imax+1)+(2)*(imax+2))+2] -
                       cpsav[3*((imax+1)+(0)*(imax+2))+2];
            t        = sqrt(r*con2*fabs(DOT(x1,snor))/con)/normnell;
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", imax+1,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((imax-1)+(0)*imax)  ] + t*nell[0] -
                 cp[3*((imax+1)+(1)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(0)*imax)+1] + t*nell[1] -
                 cp[3*((imax+1)+(1)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(0)*imax)+2] + t*nell[2] -
                 cp[3*((imax+1)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(1)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(1)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(1)*(imax+2))+2] += dz;
        } else {
            if (snor != NULL) {
                dx = (snor[0] - V1.X())*dv;
                dy = (snor[1] - V1.Y())*dv;
                dz = (snor[2] - V1.Z())*dv;
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
            dx = du * dv * UV.X();
            dy = du * dv * UV.Y();
            dz = du * dv * UV.Z();
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((imax-2)+(1)*imax)  ] - xyz[3*((imax-1)+(1)*imax)  ]) -
                 (xyz[3*((imax-2)+(0)*imax)  ] - xyz[3*((imax-1)+(0)*imax)  ]) +
                 du * dv * UV.X();
            dy = (xyz[3*((imax-2)+(1)*imax)+1] - xyz[3*((imax-1)+(1)*imax)+1]) -
                 (xyz[3*((imax-2)+(0)*imax)+1] - xyz[3*((imax-1)+(0)*imax)+1]) +
                 du * dv * UV.Y();
            dz = (xyz[3*((imax-2)+(1)*imax)+2] - xyz[3*((imax-1)+(1)*imax)+2]) -
                 (xyz[3*((imax-2)+(0)*imax)+2] - xyz[3*((imax-1)+(0)*imax)+2]) +
                 du * dv * UV.Z();
        } else {
            /* quadratic fit */
            u20 = knotu[imax+2] - knotu[imax];
            u21 = knotu[imax+1] - knotu[imax];
            for (j = 0; j < 3; j++) {
                if (easT != NULL) {
                    T[j][0] = easT[3*j  ]*du;
                    T[j][1] = easT[3*j+1]*du;
                    T[j][2] = easT[3*j+2]*du;
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
                T[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((imax-1)+(j)*imax)  ];
                T[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((imax-1)+(j)*imax)+1];
                T[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((imax-1)+(j)*imax)+2];
/*              if (iter == 0) {
                  printf(" t%d = %lf %lf %lf   %lf %lf %lf\n", j, T[j][0],
                         T[j][1], T[j][2], easT[3*j], easT[3*j+1], easT[3*j+2]);
                  printf("    = %lf %lf %lf\n",
                      xyz[3*((imax-2)+(j)*imax)  ]-xyz[3*((imax-1)+(j)*imax)  ],
                      xyz[3*((imax-2)+(j)*imax)+1]-xyz[3*((imax-1)+(j)*imax)+1],
                      xyz[3*((imax-2)+(j)*imax)+2]-xyz[3*((imax-1)+(j)*imax)+2]);
                }  */
            }
            u20    = knotv[5] - knotv[3];
            u21    = knotv[5] - knotv[4];
            rj[0]  = T[1][0]*u20*u20 - T[0][0]*u21*u21 - T[2][0]*dv*dv;
            rj[1]  = T[1][1]*u20*u20 - T[0][1]*u21*u21 - T[2][1]*dv*dv;
            rj[2]  = T[1][2]*u20*u20 - T[0][2]*u21*u21 - T[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[imax+2] - knotu[imax]);
            if (easT != NULL) u21 = du;
            dx     = rj[0] - T[0][0] + 0.5*u20*u21*UV.X();
            dy     = rj[1] - T[0][1] + 0.5*u20*u21*UV.Y();
            dz     = rj[2] - T[0][2] + 0.5*u20*u21*UV.Z();
        }
#ifdef NOSECROSS
        if ((south != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[4] - knotv[3])*(knotv[4] - knotv[3]));
            con2  = 6.0/((knotv[5] - knotv[3])*(knotv[4] - knotv[3]));
            dist  = 0.5*(knotu[iknot-5] + knotu[iknot-4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
#ifdef NEWNOSE
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
#else
            x0[0] = 0.5*(xyz[3*((imax-2)+(1)*(imax))  ] +
                         xyz[3*((imax-1)+(1)*(imax))  ]) -
                    0.5*(xyz[3*((imax-2)+(0)*(imax))  ] +
                         xyz[3*((imax-1)+(0)*(imax))  ]);
            x0[1] = 0.5*(xyz[3*((imax-2)+(1)*(imax))+1] +
                         xyz[3*((imax-1)+(1)*(imax))+1]) -
                    0.5*(xyz[3*((imax-2)+(0)*(imax))+1] +
                         xyz[3*((imax-1)+(0)*(imax))+1]);
            x0[2] = 0.5*(xyz[3*((imax-2)+(1)*(imax))+2] +
                         xyz[3*((imax-1)+(1)*(imax))+2]) -
                    0.5*(xyz[3*((imax-2)+(0)*(imax))+2] +
                         xyz[3*((imax-1)+(0)*(imax))+2]);
            x1[0] = 0.5*(xyz[3*((imax-2)+(2)*(imax))  ] +
                         xyz[3*((imax-1)+(2)*(imax))  ]) -
                    0.5*(xyz[3*((imax-2)+(0)*(imax))  ] +
                         xyz[3*((imax-1)+(0)*(imax))  ]);
            x1[1] = 0.5*(xyz[3*((imax-2)+(2)*(imax))+1] +
                         xyz[3*((imax-1)+(2)*(imax))+1]) -
                    0.5*(xyz[3*((imax-2)+(0)*(imax))+1] +
                         xyz[3*((imax-1)+(0)*(imax))+1]);
            x1[2] = 0.5*(xyz[3*((imax-2)+(2)*(imax))+2] +
                         xyz[3*((imax-1)+(2)*(imax))+2]) -
                    0.5*(xyz[3*((imax-2)+(0)*(imax))+2] +
                         xyz[3*((imax-1)+(0)*(imax))+2]);
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, snor);
                nell[0] = x0[0] - dist*snor[0];
                nell[1] = x0[1] - dist*snor[1];
                nell[2] = x0[2] - dist*snor[2];
            } else {
                CROSS(nell, norm, snor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(south, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
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
            t  = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" south %d: %lf %lf %lf (iter = %d)\n", imax,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((imax-1)+(0)*imax)  ] + t*nell[0] -
                 cp[3*((imax)+(1)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(0)*imax)+1] + t*nell[1] -
                 cp[3*((imax)+(1)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(0)*imax)+2] + t*nell[2] -
                 cp[3*((imax)+(1)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(1)*(imax+2))  ] += dx;
            cp[3*((imax)+(1)*(imax+2))+1] += dy;
            cp[3*((imax)+(1)*(imax+2))+2] += dz;
        } else {
#endif
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax)+(1)*(imax+2))+2] += RELAX * dz;
#ifdef NOSECROSS
        }
#endif
      
        /* point O */
        hSurf->D2(knotu[3], knotv[jmax+2], P0, U1, V1, U2, V2, UV);
        du = knotu[4] - knotu[3];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * U2.X();
            dy = du * du * U2.Y();
            dz = du * du * U2.Z();
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((1)+(jmax-1)*imax)  ] -
                 xyz[3*((0)+(jmax-1)*imax)  ] - du * U1.X();
            dy = xyz[3*((1)+(jmax-1)*imax)+1] -
                 xyz[3*((0)+(jmax-1)*imax)+1] - du * U1.Y();
            dz = xyz[3*((1)+(jmax-1)*imax)+2] -
                 xyz[3*((0)+(jmax-1)*imax)+2] - du * U1.Z();
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
            dx     = rj[0] - xyz[3*((0)+(jmax-1)*imax)  ] - 0.5 * u20 * U1.X();
            dy     = rj[1] - xyz[3*((0)+(jmax-1)*imax)+1] - 0.5 * u20 * U1.Y();
            dz     = rj[2] - xyz[3*((0)+(jmax-1)*imax)+2] - 0.5 * u20 * U1.Z();
        }
        if (wesT != NULL) {
            /* match input tangent */
            dx = (wesT[3*(jmax-1)  ] - U1.X())*du;
            dy = (wesT[3*(jmax-1)+1] - U1.Y())*du;
            dz = (wesT[3*(jmax-1)+2] - U1.Z())*du;
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
            dx = dv * dv * V2.X();
            dy = dv * dv * V2.Y();
            dz = dv * dv * V2.Z();
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((0)+(jmax-2)*imax)  ] -
                 xyz[3*((0)+(jmax-1)*imax)  ] + dv * V1.X();
            dy = xyz[3*((0)+(jmax-2)*imax)+1] -
                 xyz[3*((0)+(jmax-1)*imax)+1] + dv * V1.Y();
            dz = xyz[3*((0)+(jmax-2)*imax)+2] -
                 xyz[3*((0)+(jmax-1)*imax)+2] + dv * V1.Z();
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
            dx    = rj[0] - xyz[3*((0)+(jmax-1)*imax)  ] + 0.5 * u20 * V1.X();
            dy    = rj[1] - xyz[3*((0)+(jmax-1)*imax)+1] + 0.5 * u20 * V1.Y();
            dz    = rj[2] - xyz[3*((0)+(jmax-1)*imax)+2] + 0.5 * u20 * V1.Z();
        }
        if ((north != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
            con2  = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
#ifdef NEWNOSE
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
#else
            x0[0] = xyz[3*((0)+(jmax-2)*(imax))  ] -
                    xyz[3*((0)+(jmax-1)*(imax))  ];
            x0[1] = xyz[3*((0)+(jmax-2)*(imax))+1] -
                    xyz[3*((0)+(jmax-1)*(imax))+1];
            x0[2] = xyz[3*((0)+(jmax-2)*(imax))+2] -
                    xyz[3*((0)+(jmax-1)*(imax))+2];
            x1[0] = xyz[3*((0)+(jmax-3)*(imax))  ] -
                    xyz[3*((0)+(jmax-1)*(imax))  ];
            x1[1] = xyz[3*((0)+(jmax-3)*(imax))+1] -
                    xyz[3*((0)+(jmax-1)*(imax))+1];
            x1[2] = xyz[3*((0)+(jmax-3)*(imax))+2] -
                    xyz[3*((0)+(jmax-1)*(imax))+2];
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, nnor);
                nell[0] = x0[0] - dist*nnor[0];
                nell[1] = x0[1] - dist*nnor[1];
                nell[2] = x0[2] - dist*nnor[2];
            } else {
                CROSS(nell, norm, nnor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(north, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
            x1[0]    = cpsav[3*((0)+(jmax-1)*(imax+2))  ] -
                       cpsav[3*((0)+(jmax+1)*(imax+2))  ];
            x1[1]    = cpsav[3*((0)+(jmax-1)*(imax+2))+1] -
                       cpsav[3*((0)+(jmax+1)*(imax+2))+1];
            x1[2]    = cpsav[3*((0)+(jmax-1)*(imax+2))+2] -
                       cpsav[3*((0)+(jmax+1)*(imax+2))+2];
            t        = sqrt(r*con2*fabs(DOT(x1,nnor))/con)/normnell;
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", 0,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((0)+(jmax-1)*imax)  ] + t*nell[0] -
                 cp[3*((0)+(jmax)*(imax+2))  ];
            dy = xyz[3*((0)+(jmax-1)*imax)+1] + t*nell[1] -
                 cp[3*((0)+(jmax)*(imax+2))+1];
            dz = xyz[3*((0)+(jmax-1)*imax)+2] + t*nell[2] -
                 cp[3*((0)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(jmax)*(imax+2))  ] += dx;
            cp[3*((0)+(jmax)*(imax+2))+1] += dy;
            cp[3*((0)+(jmax)*(imax+2))+2] += dz;
        } else {
            if (nnor != NULL) {
                dx = (nnor[0] + V1.X())*dv;
                dy = (nnor[1] + V1.Y())*dv;
                dz = (nnor[2] + V1.Z())*dv;
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
            dx = du * dv * UV.X();
            dy = du * dv * UV.Y();
            dz = du * dv * UV.Z();
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((1)+(jmax-2)*imax)  ] - xyz[3*((0)+(jmax-2)*imax)  ]) -
                 (xyz[3*((1)+(jmax-1)*imax)  ] - xyz[3*((0)+(jmax-1)*imax)  ]) +
                 du * dv * UV.X();
            dy = (xyz[3*((1)+(jmax-2)*imax)+1] - xyz[3*((0)+(jmax-2)*imax)+1]) -
                 (xyz[3*((1)+(jmax-1)*imax)+1] - xyz[3*((0)+(jmax-1)*imax)+1]) +
                 du * dv * UV.Y();
            dz = (xyz[3*((1)+(jmax-2)*imax)+2] - xyz[3*((0)+(jmax-2)*imax)+2]) -
                 (xyz[3*((1)+(jmax-1)*imax)+2] - xyz[3*((0)+(jmax-1)*imax)+2]) +
                 du * dv * UV.Z();
        } else {
            /* quadratic fit */
            u20 = knotu[5] - knotu[3];
            u21 = knotu[5] - knotu[4];
            for (j = 0; j < 3; j++) {
                jj = jmax - j - 1;
                if (wesT != NULL) {
                    T[j][0] = wesT[3*jj  ]*du;
                    T[j][1] = wesT[3*jj+1]*du;
                    T[j][2] = wesT[3*jj+2]*du;
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
                T[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((0)+(jj)*imax)  ];
                T[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((0)+(jj)*imax)+1];
                T[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((0)+(jj)*imax)+2];
/*              if (iter == 0) {
                  printf(" t%d = %lf %lf %lf   %lf %lf %lf\n", j, T[j][0],
                      T[j][1], T[j][2], wesT[3*jj], wesT[3*jj+1], wesT[3*jj+2]);
                  printf("    = %lf %lf %lf\n",
                         xyz[3*((1)+(jj)*imax)  ] - xyz[3*((0)+(jj)*imax)  ],
                         xyz[3*((1)+(jj)*imax)+1] - xyz[3*((0)+(jj)*imax)+1],
                         xyz[3*((1)+(jj)*imax)+2] - xyz[3*((0)+(jj)*imax)+2]);
                }  */
          }
          u20    = knotv[jmax+2] - knotv[jmax];
          u21    = knotv[jmax+1] - knotv[jmax];
          rj[0]  = T[1][0]*u20*u20 - T[0][0]*u21*u21 - T[2][0]*dv*dv;
          rj[1]  = T[1][1]*u20*u20 - T[0][1]*u21*u21 - T[2][1]*dv*dv;
          rj[2]  = T[1][2]*u20*u20 - T[0][2]*u21*u21 - T[2][2]*dv*dv;
          rj[0] /= 2.0*u21*dv;
          rj[1] /= 2.0*u21*dv;
          rj[2] /= 2.0*u21*dv;
          u21    = 0.5*(knotu[5] - knotu[3]);
          if (wesT != NULL) u21 = du;
          dx     = rj[0] - T[0][0] + 0.5*u20*u21*UV.X();
          dy     = rj[1] - T[0][1] + 0.5*u20*u21*UV.Y();
          dz     = rj[2] - T[0][2] + 0.5*u20*u21*UV.Z();
        }
#ifdef NOSECROSS
        if ((north != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
            con2  = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
            dist  = 0.5*(knotu[3] + knotu[4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
#ifdef NEWNOSE
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
#else
            x0[0] = 0.5*(xyz[3*((1)+(jmax-2)*(imax))  ] +
                         xyz[3*((0)+(jmax-2)*(imax))  ]) -
                    0.5*(xyz[3*((1)+(jmax-1)*(imax))  ] +
                         xyz[3*((0)+(jmax-1)*(imax))  ]);
            x0[1] = 0.5*(xyz[3*((1)+(jmax-2)*(imax))+1] +
                         xyz[3*((0)+(jmax-2)*(imax))+1]) -
                    0.5*(xyz[3*((1)+(jmax-1)*(imax))+1] +
                         xyz[3*((0)+(jmax-1)*(imax))+1]);
            x0[2] = 0.5*(xyz[3*((1)+(jmax-2)*(imax))+2] +
                         xyz[3*((0)+(jmax-2)*(imax))+2]) -
                    0.5*(xyz[3*((1)+(jmax-1)*(imax))+2] +
                         xyz[3*((0)+(jmax-1)*(imax))+2]);
            x1[0] = 0.5*(xyz[3*((1)+(jmax-3)*(imax))  ] +
                         xyz[3*((0)+(jmax-3)*(imax))  ]) -
                    0.5*(xyz[3*((1)+(jmax-1)*(imax))  ] +
                         xyz[3*((0)+(jmax-1)*(imax))  ]);
            x1[1] = 0.5*(xyz[3*((1)+(jmax-3)*(imax))+1] +
                         xyz[3*((0)+(jmax-3)*(imax))+1]) -
                    0.5*(xyz[3*((1)+(jmax-1)*(imax))+1] +
                         xyz[3*((0)+(jmax-1)*(imax))+1]);
            x1[2] = 0.5*(xyz[3*((1)+(jmax-3)*(imax))+2] +
                         xyz[3*((0)+(jmax-3)*(imax))+2]) -
                    0.5*(xyz[3*((1)+(jmax-1)*(imax))+2] +
                         xyz[3*((0)+(jmax-1)*(imax))+2]);
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, nnor);
                nell[0] = x0[0] - dist*nnor[0];
                nell[1] = x0[1] - dist*nnor[1];
                nell[2] = x0[2] - dist*nnor[2];
            } else {
                CROSS(nell, norm, nnor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(north, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
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
            t  = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", 1,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((0)+(jmax-1)*imax)  ] + t*nell[0] -
                 cp[3*((1)+(jmax)*(imax+2))  ];
            dy = xyz[3*((0)+(jmax-1)*imax)+1] + t*nell[1] -
                 cp[3*((1)+(jmax)*(imax+2))+1];
            dz = xyz[3*((0)+(jmax-1)*imax)+2] + t*nell[2] -
                 cp[3*((1)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(jmax)*(imax+2))  ] += dx;
            cp[3*((1)+(jmax)*(imax+2))+1] += dy;
            cp[3*((1)+(jmax)*(imax+2))+2] += dz;
        } else {
#endif
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(jmax)*(imax+2))  ] += RELAX * dx;
            cp[3*((1)+(jmax)*(imax+2))+1] += RELAX * dy;
            cp[3*((1)+(jmax)*(imax+2))+2] += RELAX * dz;
#ifdef NOSECROSS
        }
#endif

        /* match north & L points */
        dv   = knotv[jmax+2] - knotv[jmax+1];
        con  = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                    (knotv[jmax+1] - knotv[jmax+2]));
        con2 = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                    (knotv[jmax+1] - knotv[jmax+2]));
        for (i = 1; i < imax-1; i++) {
            hSurf->D2(knotu[i+3], knotv[jmax+2], P0, U1, V1, U2, V2, UV);
            dx = xyz[3*((i)+(jmax-1)*imax)  ] - P0.X();
            dy = xyz[3*((i)+(jmax-1)*imax)+1] - P0.Y();
            dz = xyz[3*((i)+(jmax-1)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(jmax+1)*(imax+2))  ] += dx;
            cp[3*((i+1)+(jmax+1)*(imax+2))+1] += dy;
            cp[3*((i+1)+(jmax+1)*(imax+2))+2] += dz;

            if (endj == 0) {
                /* d2/dv2 = 0 */
                dx = dv * dv * V2.X();
                dy = dv * dv * V2.Y();
                dz = dv * dv * V2.Z();
            } else if (endj == 1) {
                /* match FD d/dv */
                dx = xyz[3*((i)+(jmax-2)*imax)  ] -
                     xyz[3*((i)+(jmax-1)*imax)  ] + dv * V1.X();
                dy = xyz[3*((i)+(jmax-2)*imax)+1] -
                     xyz[3*((i)+(jmax-1)*imax)+1] + dv * V1.Y();
                dz = xyz[3*((i)+(jmax-2)*imax)+2] -
                     xyz[3*((i)+(jmax-1)*imax)+2] + dv * V1.Z();
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
                dx     = rj[0] - xyz[3*((i)+(jmax-1)*imax)  ] + 0.5*u20*V1.X();
                dy     = rj[1] - xyz[3*((i)+(jmax-1)*imax)+1] + 0.5*u20*V1.Y();
                dz     = rj[2] - xyz[3*((i)+(jmax-1)*imax)+2] + 0.5*u20*V1.Z();
            }
            if ((north != NULL) && (cpsav != NULL)) {
                EG_NzeroBSp(iknot, knotu, knotu[i+3], &kk, basis);
#ifdef NEWNOSE
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
#else
                x0[0] = xyz[3*((i)+(jmax-2)*(imax))  ] -
                        xyz[3*((i)+(jmax-1)*(imax))  ];
                x0[1] = xyz[3*((i)+(jmax-2)*(imax))+1] -
                        xyz[3*((i)+(jmax-1)*(imax))+1];
                x0[2] = xyz[3*((i)+(jmax-2)*(imax))+2] -
                        xyz[3*((i)+(jmax-1)*(imax))+2];
                x1[0] = xyz[3*((i)+(jmax-3)*(imax))  ] -
                        xyz[3*((i)+(jmax-1)*(imax))  ];
                x1[1] = xyz[3*((i)+(jmax-3)*(imax))+1] -
                        xyz[3*((i)+(jmax-1)*(imax))+1];
                x1[2] = xyz[3*((i)+(jmax-3)*(imax))+2] -
                        xyz[3*((i)+(jmax-1)*(imax))+2];
                CROSS(norm, x0, x1);
                dist  = sqrt(DOT(norm, norm));
                if (dist < 0.000001) {
                    dist    = DOT(x0, nnor);
                    nell[0] = x0[0] - dist*nnor[0];
                    nell[1] = x0[1] - dist*nnor[1];
                    nell[2] = x0[2] - dist*nnor[2];
                } else {
                    CROSS(nell, norm, nnor);
                }
                if (DOT(x0, nell) < 0.0) {
                    nell[0] = -nell[0];
                    nell[1] = -nell[1];
                    nell[2] = -nell[2];
                }
                r        = EG_getEllRad(north, nell);
                normnell = sqrt(DOT(nell, nell));
#endif
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
                t  = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
                printf(" north %d: %lf %lf %lf (iter = %d)\n", i+1,
                       t*nell[0], t*nell[1], t*nell[2], iter);
#endif
                if (t < 0.0) t = 0.0;
                dx = xyz[3*((i)+(jmax-1)*imax)  ] + t*nell[0] -
                     cp[3*((i+1)+(jmax)*(imax+2))  ];
                dy = xyz[3*((i)+(jmax-1)*imax)+1] + t*nell[1] -
                     cp[3*((i+1)+(jmax)*(imax+2))+1];
                dz = xyz[3*((i)+(jmax-1)*imax)+2] + t*nell[2] -
                     cp[3*((i+1)+(jmax)*(imax+2))+2];
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(jmax)*(imax+2))  ] += dx;
                cp[3*((i+1)+(jmax)*(imax+2))+1] += dy;
                cp[3*((i+1)+(jmax)*(imax+2))+2] += dz;
            } else {
                if (nnor != NULL) {
                    dx = (nnor[0] + V1.X())*dv;
                    dy = (nnor[1] + V1.Y())*dv;
                    dz = (nnor[2] + V1.Z())*dv;
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
        hSurf->D2(knotu[imax+2], knotv[jmax+2], P0, U1, V1, U2, V2, UV);
        du = knotu[imax+2] - knotu[imax+1];
        if (endi == 0) {
            /* d2/du2 = 0 */
            dx = du * du * U2.X();
            dy = du * du * U2.Y();
            dz = du * du * U2.Z();
        } else if (endi == 1) {
            /* match FD d/du */
            dx = xyz[3*((imax-2)+(jmax-1)*imax)  ] -
                 xyz[3*((imax-1)+(jmax-1)*imax)  ] + du * U1.X();
            dy = xyz[3*((imax-2)+(jmax-1)*imax)+1] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+1] + du * U1.Y();
            dz = xyz[3*((imax-2)+(jmax-1)*imax)+2] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+2] + du * U1.Z();
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
            dx     = rj[0] - xyz[3*((imax-1)+(jmax-1)*imax)  ] + 0.5*u20*U1.X();
            dy     = rj[1] - xyz[3*((imax-1)+(jmax-1)*imax)+1] + 0.5*u20*U1.Y();
            dz     = rj[2] - xyz[3*((imax-1)+(jmax-1)*imax)+2] + 0.5*u20*U1.Z();
        }
        if (easT != NULL) {
            /* match input tangent */
            dx = (easT[3*(jmax-1)  ] + U1.X())*du;
            dy = (easT[3*(jmax-1)+1] + U1.Y())*du;
            dz = (easT[3*(jmax-1)+2] + U1.Z())*du;
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
            dx = dv * dv * V2.X();
            dy = dv * dv * V2.Y();
            dz = dv * dv * V2.Z();
        } else if (endj == 1) {
            /* match FD d/dv */
            dx = xyz[3*((imax-1)+(jmax-2)*imax)  ] -
                 xyz[3*((imax-1)+(jmax-1)*imax)  ] + dv * V1.X();
            dy = xyz[3*((imax-1)+(jmax-2)*imax)+1] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+1] + dv * V1.Y();
            dz = xyz[3*((imax-1)+(jmax-2)*imax)+2] -
                 xyz[3*((imax-1)+(jmax-1)*imax)+2] + dv * V1.Z();
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
            dx     = rj[0] - xyz[3*((imax-1)+(jmax-1)*imax)  ] + 0.5*u20*V1.X();
            dy     = rj[1] - xyz[3*((imax-1)+(jmax-1)*imax)+1] + 0.5*u20*V1.Y();
            dz     = rj[2] - xyz[3*((imax-1)+(jmax-1)*imax)+2] + 0.5*u20*V1.Z();
        }
        if ((north != NULL) && (cpsav != NULL)) {
            con   = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
            con2  = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                         (knotv[jmax+1] - knotv[jmax+2]));
#ifdef NEWNOSE
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
#else
            x0[0] = xyz[3*((imax-1)+(jmax-2)*(imax))  ] -
                    xyz[3*((imax-1)+(jmax-1)*(imax))  ];
            x0[1] = xyz[3*((imax-1)+(jmax-2)*(imax))+1] -
                    xyz[3*((imax-1)+(jmax-1)*(imax))+1];
            x0[2] = xyz[3*((imax-1)+(jmax-2)*(imax))+2] -
                    xyz[3*((imax-1)+(jmax-1)*(imax))+2];
            x1[0] = xyz[3*((imax-1)+(jmax-3)*(imax))  ] -
                    xyz[3*((imax-1)+(jmax-1)*(imax))  ];
            x1[1] = xyz[3*((imax-1)+(jmax-3)*(imax))+1] -
                    xyz[3*((imax-1)+(jmax-1)*(imax))+1];
            x1[2] = xyz[3*((imax-1)+(jmax-3)*(imax))+2] -
                    xyz[3*((imax-1)+(jmax-1)*(imax))+2];
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, nnor);
                nell[0] = x0[0] - dist*nnor[0];
                nell[1] = x0[1] - dist*nnor[1];
                nell[2] = x0[2] - dist*nnor[2];
            } else {
                CROSS(nell, norm, nnor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(north, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
            x1[0]    = cpsav[3*((imax+1)+(jmax-1)*(imax+2))  ] -
                       cpsav[3*((imax+1)+(jmax+1)*(imax+2))  ];
            x1[1]    = cpsav[3*((imax+1)+(jmax-1)*(imax+2))+1] -
                       cpsav[3*((imax+1)+(jmax+1)*(imax+2))+1];
            x1[2]    = cpsav[3*((imax+1)+(jmax-1)*(imax+2))+2] -
                       cpsav[3*((imax+1)+(jmax+1)*(imax+2))+2];
            t        = sqrt(r*con2*fabs(DOT(x1,nnor))/con)/normnell;
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", imax+1,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((imax-1)+(jmax-1)*imax)  ] + t*nell[0] -
                 cp[3*((imax+1)+(jmax)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(jmax-1)*imax)+1] + t*nell[1] -
                 cp[3*((imax+1)+(jmax)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(jmax-1)*imax)+2] + t*nell[2] -
                 cp[3*((imax+1)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(jmax)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(jmax)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(jmax)*(imax+2))+2] += dz;
        } else {
            if (nnor != NULL) {
                dx = (nnor[0] + V1.X())*dv;
                dy = (nnor[1] + V1.Y())*dv;
                dz = (nnor[2] + V1.Z())*dv;
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
            dx = du * dv * UV.X();
            dy = du * dv * UV.Y();
            dz = du * dv * UV.Z();
        } else if ((endi == 1) || (endj == 1)) {
            /* match FD cross term */
            dx = (xyz[3*((imax-2)+(jmax-2)*imax)  ] -
                  xyz[3*((imax-1)+(jmax-2)*imax)  ])  -
                 (xyz[3*((imax-2)+(jmax-1)*imax)  ] -
                  xyz[3*((imax-1)+(jmax-1)*imax)  ])  + du * dv * UV.X();
            dy = (xyz[3*((imax-2)+(jmax-2)*imax)+1] -
                  xyz[3*((imax-1)+(jmax-2)*imax)+1])  -
                 (xyz[3*((imax-2)+(jmax-1)*imax)+1] -
                  xyz[3*((imax-1)+(jmax-1)*imax)+1])  + du * dv * UV.Y();
            dz = (xyz[3*((imax-2)+(jmax-2)*imax)+2] -
                  xyz[3*((imax-1)+(jmax-2)*imax)+2])  -
                 (xyz[3*((imax-2)+(jmax-1)*imax)+2] -
                  xyz[3*((imax-1)+(jmax-1)*imax)+2])  + du * dv * UV.Z();
        } else {
            /* quadratic fit */
            u20 = knotu[imax+2] - knotu[imax];
            u21 = knotu[imax+1] - knotu[imax];
            for (j = 0; j < 3; j++) {
                jj = jmax - j - 1;
                if (easT != NULL) {
                    T[j][0] = easT[3*jj  ]*du;
                    T[j][1] = easT[3*jj+1]*du;
                    T[j][2] = easT[3*jj+2]*du;
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
                T[j][0] = rj[0]/(2.0*u21*du) - xyz[3*((imax-1)+(jj)*imax)  ];
                T[j][1] = rj[1]/(2.0*u21*du) - xyz[3*((imax-1)+(jj)*imax)+1];
                T[j][2] = rj[2]/(2.0*u21*du) - xyz[3*((imax-1)+(jj)*imax)+2];
/*              if (iter == 0) {
                  printf(" t%d = %lf %lf %lf   %lf %lf %lf\n", j, T[j][0],
                      T[j][1], T[j][2], easT[3*jj], easT[3*jj+1], easT[3*jj+2]);
                  printf("    = %lf %lf %lf\n",
                    xyz[3*((imax-2)+(jj)*imax)  ]-xyz[3*((imax-1)+(jj)*imax)  ],
                    xyz[3*((imax-2)+(jj)*imax)+1]-xyz[3*((imax-1)+(jj)*imax)+1],
                    xyz[3*((imax-2)+(jj)*imax)+2]-xyz[3*((imax-1)+(jj)*imax)+2]);
                }  */
            }
            u20    = knotv[jmax+2] - knotv[jmax];
            u21    = knotv[jmax+1] - knotv[jmax];
            rj[0]  = T[1][0]*u20*u20 - T[0][0]*u21*u21 - T[2][0]*dv*dv;
            rj[1]  = T[1][1]*u20*u20 - T[0][1]*u21*u21 - T[2][1]*dv*dv;
            rj[2]  = T[1][2]*u20*u20 - T[0][2]*u21*u21 - T[2][2]*dv*dv;
            rj[0] /= 2.0*u21*dv;
            rj[1] /= 2.0*u21*dv;
            rj[2] /= 2.0*u21*dv;
            u21    = 0.5*(knotu[imax+2] - knotu[imax]);
            if (easT != NULL) u21 = du;
            dx     = -rj[0] - T[0][0] - 0.5*u20*u21*UV.X();
            dy     = -rj[1] - T[0][1] - 0.5*u20*u21*UV.Y();
            dz     = -rj[2] - T[0][2] - 0.5*u20*u21*UV.Z();
        }
#ifdef NOSECROSS
        if ((north != NULL) && (cpsav != NULL)) {
            con  = 9.0/((knotv[jmax+1] - knotv[jmax+2])*
                        (knotv[jmax+1] - knotv[jmax+2]));
            con2 = 6.0/((knotv[jmax  ] - knotv[jmax+2])*
                        (knotv[jmax+1] - knotv[jmax+2]));
            dist  = 0.5*(knotu[iknot-5] + knotu[iknot-4]);
            EG_NzeroBSp(iknot, knotu, dist, &kk, basis);
#ifdef NEWNOSE
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
#else
            x0[0] = 0.5*(xyz[3*((imax-2)+(jmax-2)*(imax))  ] +
                         xyz[3*((imax-1)+(jmax-2)*(imax))  ]) -
                    0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))  ] +
                         xyz[3*((imax-1)+(jmax-1)*(imax))  ]);
            x0[1] = 0.5*(xyz[3*((imax-2)+(jmax-2)*(imax))+1] +
                         xyz[3*((imax-1)+(jmax-2)*(imax))+1]) -
                    0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))+1] +
                         xyz[3*((imax-1)+(jmax-1)*(imax))+1]);
            x0[2] = 0.5*(xyz[3*((imax-2)+(jmax-2)*(imax))+2] +
                         xyz[3*((imax-1)+(jmax-2)*(imax))+2]) -
                    0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))+2] +
                         xyz[3*((imax-1)+(jmax-1)*(imax))+2]);
            x1[0] = 0.5*(xyz[3*((imax-2)+(jmax-3)*(imax))  ] +
                         xyz[3*((imax-1)+(jmax-3)*(imax))  ]) -
                    0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))  ] +
                         xyz[3*((imax-1)+(jmax-1)*(imax))  ]);
            x1[1] = 0.5*(xyz[3*((imax-2)+(jmax-3)*(imax))+1] +
                         xyz[3*((imax-1)+(jmax-3)*(imax))+1]) -
                    0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))+1] +
                         xyz[3*((imax-1)+(jmax-1)*(imax))+1]);
            x1[2] = 0.5*(xyz[3*((imax-2)+(jmax-3)*(imax))+2] +
                         xyz[3*((imax-1)+(jmax-3)*(imax))+2]) -
                    0.5*(xyz[3*((imax-2)+(jmax-1)*(imax))+2] +
                         xyz[3*((imax-1)+(jmax-1)*(imax))+2]);
            CROSS(norm, x0, x1);
            dist  = sqrt(DOT(norm, norm));
            if (dist < 0.000001) {
                dist    = DOT(x0, nnor);
                nell[0] = x0[0] - dist*nnor[0];
                nell[1] = x0[1] - dist*nnor[1];
                nell[2] = x0[2] - dist*nnor[2];
            } else {
                CROSS(nell, norm, nnor);
            }
            if (DOT(x0, nell) < 0.0) {
                nell[0] = -nell[0];
                nell[1] = -nell[1];
                nell[2] = -nell[2];
            }
            r        = EG_getEllRad(north, nell);
            normnell = sqrt(DOT(nell, nell));
#endif
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
            t  = (-eE + sqrt(eE*eE-4.0*dD*(F-r*G)))/(2.0*dD);
#ifdef DEBUG
            printf(" north %d: %lf %lf %lf (iter = %d)\n", imax,
                   t*nell[0], t*nell[1], t*nell[2], iter);
#endif
            if (t < 0.0) t = 0.0;
            dx = xyz[3*((imax-1)+(jmax-1)*imax)  ] + t*nell[0] -
                 cp[3*((imax)+(jmax)*(imax+2))  ];
            dy = xyz[3*((imax-1)+(jmax-1)*imax)+1] + t*nell[1] -
                 cp[3*((imax)+(jmax)*(imax+2))+1];
            dz = xyz[3*((imax-1)+(jmax-1)*imax)+2] + t*nell[2] -
                 cp[3*((imax)+(jmax)*(imax+2))+2];
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(jmax)*(imax+2))  ] += dx;
            cp[3*((imax)+(jmax)*(imax+2))+1] += dy;
            cp[3*((imax)+(jmax)*(imax+2))+2] += dz;
        } else {
#endif
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(jmax)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax)+(jmax)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax)+(jmax)*(imax+2))+2] += RELAX * dz;
#ifdef NOSECROSS
        }
#endif

        /* match west & H points */
        du = knotu[4] - knotu[3];
        for (j = 1; j < jmax-1; j++) {
            hSurf->D2(knotu[3], knotv[j+3], P0, U1, V1, U2, V2, UV);
            dx = xyz[3*((0)+(j)*imax)  ] - P0.X();
            dy = xyz[3*((0)+(j)*imax)+1] - P0.Y();
            dz = xyz[3*((0)+(j)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(j+1)*(imax+2))  ] += dx;
            cp[3*((0)+(j+1)*(imax+2))+1] += dy;
            cp[3*((0)+(j+1)*(imax+2))+2] += dz;

            if (endi == 0) {
                /* d2/du2 = 0 */
                dx = du * du * U2.X();
                dy = du * du * U2.Y();
                dz = du * du * U2.Z();
            } else if (endi == 1) {
                /* match FD d/du */
                dx = xyz[3*((1)+(j)*imax)  ] -
                     xyz[3*((0)+(j)*imax)  ] - du * U1.X();
                dy = xyz[3*((1)+(j)*imax)+1] -
                     xyz[3*((0)+(j)*imax)+1] - du * U1.Y();
                dz = xyz[3*((1)+(j)*imax)+2] -
                     xyz[3*((0)+(j)*imax)+2] - du * U1.Z();
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
                dx     = rj[0] - xyz[3*((0)+(j)*imax)  ] - 0.5 * u20 * U1.X();
                dy     = rj[1] - xyz[3*((0)+(j)*imax)+1] - 0.5 * u20 * U1.Y();
                dz     = rj[2] - xyz[3*((0)+(j)*imax)+2] - 0.5 * u20 * U1.Z();
            }
            if (wesT != NULL) {
                /* match input tangent */
                dx = (wesT[3*j  ] - U1.X())*du;
                dy = (wesT[3*j+1] - U1.Y())*du;
                dz = (wesT[3*j+2] - U1.Z())*du;
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
        for (j = 1; j < jmax-1; j++) {
            hSurf->D2(knotu[imax+2], knotv[j+3], P0, U1, V1, U2, V2, UV);
            dx = xyz[3*((imax-1)+(j)*imax)  ] - P0.X();
            dy = xyz[3*((imax-1)+(j)*imax)+1] - P0.Y();
            dz = xyz[3*((imax-1)+(j)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(j+1)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(j+1)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(j+1)*(imax+2))+2] += dz;

            if (endi == 0) {
                /* d2/du2 = 0 */
                dx = du * du * U2.X();
                dy = du * du * U2.Y();
                dz = du * du * U2.Z();
            } else if (endi == 1) {
                /* match FD d/du */
                dx = xyz[3*((imax-2)+(j)*imax)  ] -
                     xyz[3*((imax-1)+(j)*imax)  ] + du * U1.X();
                dy = xyz[3*((imax-2)+(j)*imax)+1] -
                     xyz[3*((imax-1)+(j)*imax)+1] + du * U1.Y();
                dz = xyz[3*((imax-2)+(j)*imax)+2] -
                     xyz[3*((imax-1)+(j)*imax)+2] + du * U1.Z();
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
                dx     = rj[0] - xyz[3*((imax-1)+(j)*imax)  ] + 0.5*u20*U1.X();
                dy     = rj[1] - xyz[3*((imax-1)+(j)*imax)+1] + 0.5*u20*U1.Y();
                dz     = rj[2] - xyz[3*((imax-1)+(j)*imax)+2] + 0.5*u20*U1.Z();
            }
            if (easT != NULL) {
                /* match input tangent */
                dx = (easT[3*j  ] + U1.X())*du;
                dy = (easT[3*j+1] + U1.Y())*du;
                dz = (easT[3*j+2] + U1.Z())*du;
            }
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(j+1)*(imax+2))  ] += RELAX * dx;
            cp[3*((imax)+(j+1)*(imax+2))+1] += RELAX * dy;
            cp[3*((imax)+(j+1)*(imax+2))+2] += RELAX * dz;
        }

        /* convergence check */
        if (dxyzmax < tol) break;

        /* make the new surface (after deleting old one) */
        status = EG_deleteObject(*esurf);
        if (status != EGADS_SUCCESS) {
            *esurf = NULL;
            if (cpsav != NULL) EG_free(cpsav);
            EG_free(rvec);
            return status;
        }
        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rvec,
                                 esurf);
        if (status != EGADS_SUCCESS) {
            if (cpsav != NULL) EG_free(cpsav);
            EG_free(rvec);
            return status;
        }
    }
#ifdef DEBUG
    printf(" EG_spline2d: Convergence at iteration %d (%le)!\n", iter, dxyzmax);
#endif
    if (dxyzmax >= tol)
        printf(" EGADS Warning: Not Converged (EG_spline2d)!\n");

    if (cpsav != NULL) EG_free(cpsav);
    EG_free(rvec);

    return EGADS_SUCCESS;
}
#endif


/*
 ************************************************************************
 *                                                                      *
 *   EG_spline2dFit - Fit a 2d cubic spline from input data             *
 *                                                                      *
 ************************************************************************
 */

int
EG_spline2dFit(egObject *context, const double *crosT, int imax,
               const double *uknot, const double *souT, const double *norT,
               int jmax, const double *vknot, const double *wesT,
               const double *easT, const double *xyz, double tol,
               egObject **esurf)
{
    int    i, j, kk, iknot, jknot, icp, jcp, iter, header[7], status;
    double du, dv, dx, dy, dz, dxyzmax, *rvec, *knotu, *knotv, *cp;
    gp_Pnt P0;
    gp_Vec V1, V2, U1, U2, UV;

    *esurf = NULL;
    icp    = imax + 2;
    iknot  = imax + 6;
    jcp    = jmax + 2;
    jknot  = jmax + 6;
    rvec   = (double *) EG_alloc((iknot+jknot+3*icp*jcp)*sizeof(double));
    if (rvec == NULL) return EGADS_MALLOC;
    knotu  =  rvec;
    knotv  = &rvec[iknot      ];
    cp     = &rvec[iknot+jknot];

    /* create spline surface */
    header[0] = 0;
    header[1] = 3;                      /* cubic only */
    header[2] = icp;
    header[3] = iknot;
    header[4] = 3;                      /* cubic only */
    header[5] = jcp;
    header[6] = jknot;

    /* knots in i-direction */
    kk          = 0;
    knotu[kk++] = uknot[0];
    knotu[kk++] = uknot[0];
    knotu[kk++] = uknot[0];
    for (i = 0; i < imax; i++) knotu[kk++] = uknot[i];
    knotu[kk++] = uknot[imax-1];
    knotu[kk++] = uknot[imax-1];
    knotu[kk++] = uknot[imax-1];

    /* knots in j-direction */
    kk          = 0;
    knotv[kk++] = vknot[0];
    knotv[kk++] = vknot[0];
    knotv[kk++] = vknot[0];
    for (j = 0; j < jmax; j++) knotv[kk++] = vknot[j];
    knotv[kk++] = vknot[jmax-1];
    knotv[kk++] = vknot[jmax-1];
    knotv[kk++] = vknot[jmax-1];

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
     
       crosT -> D, F, M, K

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
  
    /* make the original BSPLINE (based upon the assumed control points) */
    status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rvec,
                             esurf);
    if (status != EGADS_SUCCESS) {
        EG_free(rvec);
        return status;
    }

    /* iterate to have knot evaluations match data points */
    for (iter = 0; iter < NITER; iter++) {
      
        egadsSurface *psurf        = (egadsSurface *) (*esurf)->blind;
        Handle(Geom_Surface) hSurf = psurf->handle;
        dxyzmax = 0.0;

        /* match interior spline points */
        for (j = 1; j < jmax-1; j++) {
            for (i = 1; i < imax-1; i++) {
                hSurf->D0(knotu[i+3], knotv[j+3], P0);
                dx = xyz[3*((i)+(j)*imax)  ] - P0.X();
                dy = xyz[3*((i)+(j)*imax)+1] - P0.Y();
                dz = xyz[3*((i)+(j)*imax)+2] - P0.Z();
                if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
                if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
                if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
                cp[3*((i+1)+(j+1)*(imax+2))  ] += dx;
                cp[3*((i+1)+(j+1)*(imax+2))+1] += dy;
                cp[3*((i+1)+(j+1)*(imax+2))+2] += dz;
            }
        }

        /* point A */
        hSurf->D2(knotu[3], knotv[3], P0, U1, V1, U2, V2, UV);
        du = knotu[4] - knotu[3];
        dx = (wesT[0] - U1.X())*du;
        dy = (wesT[1] - U1.Y())*du;
        dz = (wesT[2] - U1.Z())*du;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((1)+(0)*(imax+2))  ] += RELAX * dx;
        cp[3*((1)+(0)*(imax+2))+1] += RELAX * dy;
        cp[3*((1)+(0)*(imax+2))+2] += RELAX * dz;
      
        /* point C */
        dv = knotv[4] - knotv[3];
        dx = (souT[0] - V1.X())*dv;
        dy = (souT[1] - V1.Y())*dv;
        dz = (souT[2] - V1.Z())*dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((0)+(1)*(imax+2))  ] += RELAX * dx;
        cp[3*((0)+(1)*(imax+2))+1] += RELAX * dy;
        cp[3*((0)+(1)*(imax+2))+2] += RELAX * dz;
      
        /* point D */
        dx = (crosT[0] - UV.X()) * du * dv;
        dy = (crosT[1] - UV.Y()) * du * dv;
        dz = (crosT[2] - UV.Z()) * du * dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((1)+(1)*(imax+2))  ] += RELAX * dx;
        cp[3*((1)+(1)*(imax+2))+1] += RELAX * dy;
        cp[3*((1)+(1)*(imax+2))+2] += RELAX * dz;
      
        /* match south & E points */
        dv = knotv[4] - knotv[3];
        for (i = 1; i < imax-1; i++) {
            hSurf->D1(knotu[i+3], knotv[3], P0, U1, V1);
  
            dx = xyz[3*((i)+(0)*imax)  ] - P0.X();
            dy = xyz[3*((i)+(0)*imax)+1] - P0.Y();
            dz = xyz[3*((i)+(0)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(0)*(imax+2))  ] += dx;
            cp[3*((i+1)+(0)*(imax+2))+1] += dy;
            cp[3*((i+1)+(0)*(imax+2))+2] += dz;

            dx = (souT[3*i  ] - V1.X())*dv;
            dy = (souT[3*i+1] - V1.Y())*dv;
            dz = (souT[3*i+2] - V1.Z())*dv;
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(1)*(imax+2))  ] += RELAX * dx;
            cp[3*((i+1)+(1)*(imax+2))+1] += RELAX * dy;
            cp[3*((i+1)+(1)*(imax+2))+2] += RELAX * dz;
        }
      
        /* point B */
        hSurf->D2(knotu[imax+2], knotv[3], P0, U1, V1, U2, V2, UV);
        du = knotu[imax+2] - knotu[imax+1];
        dx = (easT[0] - U1.X())*du;
        dy = (easT[1] - U1.Y())*du;
        dz = (easT[2] - U1.Z())*du;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax)+(0)*(imax+2))  ] -= RELAX * dx;
        cp[3*((imax)+(0)*(imax+2))+1] -= RELAX * dy;
        cp[3*((imax)+(0)*(imax+2))+2] -= RELAX * dz;
      
        /* point G */
        dv = knotv[4] - knotv[3];
        dx = (souT[3*(imax-1)  ] - V1.X())*dv;
        dy = (souT[3*(imax-1)+1] - V1.Y())*dv;
        dz = (souT[3*(imax-1)+2] - V1.Z())*dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax+1)+(1)*(imax+2))  ] += RELAX * dx;
        cp[3*((imax+1)+(1)*(imax+2))+1] += RELAX * dy;
        cp[3*((imax+1)+(1)*(imax+2))+2] += RELAX * dz;

        /* point F */
        dx = (crosT[3] - UV.X()) * du * dv;
        dy = (crosT[4] - UV.Y()) * du * dv;
        dz = (crosT[5] - UV.Z()) * du * dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax)+(1)*(imax+2))  ] -= RELAX * dx;
        cp[3*((imax)+(1)*(imax+2))+1] -= RELAX * dy;
        cp[3*((imax)+(1)*(imax+2))+2] -= RELAX * dz;
      
        /* point O */
        hSurf->D2(knotu[3], knotv[jmax+2], P0, U1, V1, U2, V2, UV);
        du = knotu[4] - knotu[3];
        dx = (wesT[3*(jmax-1)  ] - U1.X())*du;
        dy = (wesT[3*(jmax-1)+1] - U1.Y())*du;
        dz = (wesT[3*(jmax-1)+2] - U1.Z())*du;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((1)+(jmax+1)*(imax+2))  ] += RELAX * dx;
        cp[3*((1)+(jmax+1)*(imax+2))+1] += RELAX * dy;
        cp[3*((1)+(jmax+1)*(imax+2))+2] += RELAX * dz;
      
        /* point J */
        dv = knotv[jmax+2] - knotv[jmax+1];
        dx = (norT[0] - V1.X())*dv;
        dy = (norT[1] - V1.Y())*dv;
        dz = (norT[2] - V1.Z())*dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((0)+(jmax)*(imax+2))  ] -= RELAX * dx;
        cp[3*((0)+(jmax)*(imax+2))+1] -= RELAX * dy;
        cp[3*((0)+(jmax)*(imax+2))+2] -= RELAX * dz;

        /* point K */
        dx = (crosT[ 9] - UV.X()) * du * dv;
        dy = (crosT[10] - UV.Y()) * du * dv;
        dz = (crosT[11] - UV.Z()) * du * dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((1)+(jmax)*(imax+2))  ] -= RELAX * dx;
        cp[3*((1)+(jmax)*(imax+2))+1] -= RELAX * dy;
        cp[3*((1)+(jmax)*(imax+2))+2] -= RELAX * dz;

        /* match north & L points */
        dv = knotv[jmax+2] - knotv[jmax+1];
        for (i = 1; i < imax-1; i++) {
            hSurf->D1(knotu[i+3], knotv[jmax+2], P0, U1, V1);

            dx = xyz[3*((i)+(jmax-1)*imax)  ] - P0.X();
            dy = xyz[3*((i)+(jmax-1)*imax)+1] - P0.Y();
            dz = xyz[3*((i)+(jmax-1)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(jmax+1)*(imax+2))  ] += dx;
            cp[3*((i+1)+(jmax+1)*(imax+2))+1] += dy;
            cp[3*((i+1)+(jmax+1)*(imax+2))+2] += dz;

            dx = (norT[3*i  ] - V1.X())*dv;
            dy = (norT[3*i+1] - V1.Y())*dv;
            dz = (norT[3*i+2] - V1.Z())*dv;
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((i+1)+(jmax)*(imax+2))  ] -= RELAX * dx;
            cp[3*((i+1)+(jmax)*(imax+2))+1] -= RELAX * dy;
            cp[3*((i+1)+(jmax)*(imax+2))+2] -= RELAX * dz;
        }

        /* point P */
        hSurf->D2(knotu[imax+2], knotv[jmax+2], P0, U1, V1, U2, V2, UV);
        du = knotu[imax+2] - knotu[imax+1];
        dx = (easT[3*(jmax-1)  ] - U1.X())*du;
        dy = (easT[3*(jmax-1)+1] - U1.Y())*du;
        dz = (easT[3*(jmax-1)+2] - U1.Z())*du;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax)+(jmax+1)*(imax+2))  ] -= RELAX * dx;
        cp[3*((imax)+(jmax+1)*(imax+2))+1] -= RELAX * dy;
        cp[3*((imax)+(jmax+1)*(imax+2))+2] -= RELAX * dz;
      
        /* point N */
        dv = knotv[jmax+2] - knotv[jmax+1];
        dx = (norT[3*(imax-1)  ] - V1.X())*dv;
        dy = (norT[3*(imax-1)+1] - V1.Y())*dv;
        dz = (norT[3*(imax-1)+2] - V1.Z())*dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax+1)+(jmax)*(imax+2))  ] -= RELAX * dx;
        cp[3*((imax+1)+(jmax)*(imax+2))+1] -= RELAX * dy;
        cp[3*((imax+1)+(jmax)*(imax+2))+2] -= RELAX * dz;

        /* point M */
        dx = (crosT[6] - UV.X()) * du * dv;
        dy = (crosT[7] - UV.Y()) * du * dv;
        dz = (crosT[8] - UV.Z()) * du * dv;
        if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
        if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
        if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
        cp[3*((imax)+(jmax)*(imax+2))  ] += RELAX * dx;
        cp[3*((imax)+(jmax)*(imax+2))+1] += RELAX * dy;
        cp[3*((imax)+(jmax)*(imax+2))+2] += RELAX * dz;

        /* match west & H points */
        du = knotu[4] - knotu[3];
        for (j = 1; j < jmax-1; j++) {
            hSurf->D1(knotu[3], knotv[j+3], P0, U1, V1);

            dx = xyz[3*((0)+(j)*imax)  ] - P0.X();
            dy = xyz[3*((0)+(j)*imax)+1] - P0.Y();
            dz = xyz[3*((0)+(j)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((0)+(j+1)*(imax+2))  ] += dx;
            cp[3*((0)+(j+1)*(imax+2))+1] += dy;
            cp[3*((0)+(j+1)*(imax+2))+2] += dz;

            dx = (wesT[3*j  ] - U1.X())*du;
            dy = (wesT[3*j+1] - U1.Y())*du;
            dz = (wesT[3*j+2] - U1.Z())*du;
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((1)+(j+1)*(imax+2))  ] += RELAX * dx;
            cp[3*((1)+(j+1)*(imax+2))+1] += RELAX * dy;
            cp[3*((1)+(j+1)*(imax+2))+2] += RELAX * dz;
        }

        /* match east & I points */
        du = knotu[imax+2] - knotu[imax+1];
        for (j = 1; j < jmax-1; j++) {
            hSurf->D1(knotu[imax+2], knotv[j+3], P0, U1, V1);

            dx = xyz[3*((imax-1)+(j)*imax)  ] - P0.X();
            dy = xyz[3*((imax-1)+(j)*imax)+1] - P0.Y();
            dz = xyz[3*((imax-1)+(j)*imax)+2] - P0.Z();
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax+1)+(j+1)*(imax+2))  ] += dx;
            cp[3*((imax+1)+(j+1)*(imax+2))+1] += dy;
            cp[3*((imax+1)+(j+1)*(imax+2))+2] += dz;

            dx = (easT[3*j  ] - U1.X())*du;
            dy = (easT[3*j+1] - U1.Y())*du;
            dz = (easT[3*j+2] - U1.Z())*du;
            if (fabs(dx) > dxyzmax) dxyzmax = fabs(dx);
            if (fabs(dy) > dxyzmax) dxyzmax = fabs(dy);
            if (fabs(dz) > dxyzmax) dxyzmax = fabs(dz);
            cp[3*((imax)+(j+1)*(imax+2))  ] -= RELAX * dx;
            cp[3*((imax)+(j+1)*(imax+2))+1] -= RELAX * dy;
            cp[3*((imax)+(j+1)*(imax+2))+2] -= RELAX * dz;
        }

        /* convergence check */
        if (dxyzmax < tol) break;

        /* make the new surface (after deleting old one) */
        status = EG_deleteObject(*esurf);
        if (status != EGADS_SUCCESS) {
            *esurf = NULL;
            EG_free(rvec);
            return status;
        }
        status = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, header, rvec,
                                 esurf);
        if (status != EGADS_SUCCESS) {
            EG_free(rvec);
            return status;
        }
    }
#ifdef DEBUG
    printf(" EG_spline2dFit: Convergence at iteration %d (%le)!\n", iter, dxyzmax);
#endif
    if (dxyzmax >= tol)
        printf(" EGADS Warning: Not Converged (EG_spline2dFit)!\n");

    EG_free(rvec);

    return EGADS_SUCCESS;
}
