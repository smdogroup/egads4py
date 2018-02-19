/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Blend & Rule Functions
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

//#define SHARPER

#define PI               3.1415926535897931159979635

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


extern /*@kept@*/ /*@null@*/ egObject *EG_context( const ego object );
extern int EG_outLevel( const ego object );
extern int EG_isPlanar( const ego object );
extern int EG_getPlane( const ego object, ego *plane );
extern int EG_spline2dAppx( ego context,     int    endc,
                            /*@null@*/ const double *uknot,
                            /*@null@*/ const double *vknot,
                            /*@null@*/ const int    *vdata,
                            /*@null@*/ const double *wesT,
                            /*@null@*/ const double *easT,
                            /*@null@*/ const double *south,
                            /*@null@*/       double *snor,
                            /*@null@*/ const double *north,
                            /*@null@*/       double *nnor, int imax, int jmax,
                                       const double *xyz,  double tol,
                            ego *esurf );


typedef struct {
  int    ncp;                   /* number of control points */
  int    nave;                  /* number used to average positions */
  double *knots;                /* the knot positions */
} egSequ;


/* knot sequence utility functions */

int
EG_allocSeq(int nstripe, egSequ **sequ)
{
  int    i;
  egSequ *seq;
  
  *sequ = NULL;
  seq   = (egSequ *) EG_alloc(nstripe*sizeof(egSequ));
  if (seq == NULL) return EGADS_MALLOC;
  for (i = 0; i < nstripe; i++) {
    seq[i].ncp   = 0;
    seq[i].nave  = 0;
    seq[i].knots = NULL;
  }
  
  *sequ = seq;
  return EGADS_SUCCESS;
}


void
EG_freeSeq(int nstripe, /*@only@*/ egSequ *seq)
{
  int i;
  
  for (i = 0; i < nstripe; i++)
    if (seq[i].knots != NULL) EG_free(seq[i].knots);
  
  EG_free(seq);
}


int
EG_setSeq(int stripe, egSequ *seq, int num, /*@null@*/ int    *iinfo,
                                            /*@null@*/ double *rinfo)
{
  int    i, nk, n, ndeg;
  double dt, *nomulti;

  if ((iinfo == NULL) || (rinfo == NULL)) {
    /* not a spline */
    if (seq[stripe].ncp >= num) return EGADS_SUCCESS;
    if (seq[stripe].knots != NULL) EG_free(seq[stripe].knots);
    seq[stripe].knots = (double *) EG_alloc(num*sizeof(double));
    if (seq[stripe].knots == NULL) return EGADS_MALLOC;
    /* equi-spaced sampling */
    dt = 1.0/(num-1.0);
    for (i = 0; i < num; i++) seq[stripe].knots[i] = i*dt;
    seq[stripe].nave = 0;
    seq[stripe].ncp  = num;
  } else {
    ndeg = iinfo[1];
    nk   = iinfo[3] - 2*ndeg;
    nomulti = (double *) EG_alloc(nk*sizeof(double));
    if (nomulti == NULL) return EGADS_MALLOC;
    nomulti[0] = rinfo[ndeg];
    for (n = i = 1; i < nk; i++) {
      nomulti[n] = rinfo[i+ndeg];
      if (nomulti[n] != nomulti[n-1]) n++;
    }
/*  for (i = 0; i < n; i++) printf(" %d/%d: %lf\n", i+1, n, nomulti[i]);  */
    if (seq[stripe].ncp > n) {
      EG_free(nomulti);
      return EGADS_SUCCESS;
    }
    if ((seq[stripe].ncp < num) && (n < num)) {
      EG_free(nomulti);
      if (seq[stripe].knots != NULL) EG_free(seq[stripe].knots);
      seq[stripe].knots = (double *) EG_alloc(num*sizeof(double));
      if (seq[stripe].knots == NULL) return EGADS_MALLOC;
      /* equi-spaced sampling */
      dt = 1.0/(num-1.0);
      for (i = 0; i < num; i++) seq[stripe].knots[i] = i*dt;
      seq[stripe].nave = 0;
      seq[stripe].ncp  = num;
      return EGADS_SUCCESS;
    }
    /* same number of knots -- sum with previous */
    if (seq[stripe].ncp == n) {
      if (seq[stripe].nave == 0) {
        EG_free(nomulti);     /* keep equi-spaced */
        return EGADS_SUCCESS;
      }
      for (i = 0; i < n; i++) {
        seq[stripe].knots[i] *= seq[stripe].nave;
        seq[stripe].knots[i] += (nomulti[i  ] - nomulti[0])/
                                (nomulti[n-1] - nomulti[0]);
      }
      dt = seq[stripe].knots[n-1];
      for (i = 0; i < n; i++) seq[stripe].knots[i] /= dt;
      seq[stripe].nave++;
      EG_free(nomulti);
      return EGADS_SUCCESS;
    }
    /* greater sequence -- use it */
    if (seq[stripe].knots != NULL) EG_free(seq[stripe].knots);
    seq[stripe].knots = (double *) EG_alloc(n*sizeof(double));
    if (seq[stripe].knots == NULL) {
      EG_free(nomulti);
      return EGADS_MALLOC;
    }
    for (i = 0; i < n; i++)
      seq[stripe].knots[i] = nomulti[i] - nomulti[0];
    dt = seq[stripe].knots[n-1];
    for (i = 0; i < n; i++) seq[stripe].knots[i] /= dt;
    seq[stripe].nave = 1;
    seq[stripe].ncp  = n;
    EG_free(nomulti);
  }
  
  return EGADS_SUCCESS;
}


/*
 ************************************************************************
 *                                                                      *
 *   EG_spline2d - create 2d cubic spline from input data               *
 *                                                                      *
 ************************************************************************
 */

int EG_spline2d(ego context, int endc, /*@null@*/ const double **drnd,
                int imax, int jmax, const double *xyz, double tol, ego *esurf)
{
    int          i, j, k, stat, outLevel;
    double       dist, dy, *norT, *souT, *easT, *wesT;
    double       x0[3], x1[3], nnor[3], snor[3], enor[3], wnor[3], *rot;
    const double *north, *south, *east, *west;

    *esurf = NULL;
    if (context == NULL)               return EGADS_NULLOBJ;
    if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
    if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
    if ((imax < 2) || (jmax < 2))      return EGADS_DEGEN;
    if ((endc < 0) || (endc > 2))      return EGADS_RANGERR;
    outLevel = EG_outLevel(context);

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
                    if (outLevel > 0)
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
                    if (outLevel > 0)
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
                    if (outLevel > 0)
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
                    if (outLevel > 0)
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
    if ((north != NULL) && ((east != NULL) || (west != NULL)))
      return EGADS_DEGEN;
    if ((south != NULL) && ((east != NULL) || (west != NULL)))
      return EGADS_DEGEN;
  
    /* call approx as is */
    if ((east == NULL) && (west == NULL))
      return EG_spline2dAppx(context, endc, NULL, NULL, NULL, NULL, NULL,
                             south, souT, north, norT, imax, jmax, xyz, tol,
                             esurf);
        
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
    stat = EG_spline2dAppx(context, endc, NULL, NULL, NULL, NULL, NULL,
                           west, wesT, east, easT, jmax, imax, rot, tol,
                           esurf);

    EG_free(rot);
    return stat;
}
        

static int 
EG_loft2spline(ego context, const double *vknot, const double **drnd,
               egSequ *seq, int jmax, const double *xyz, double tol, ego *esurf)
{
    int          i, j, imax, outLevel, *vdata;
    double       dist, dy;
    double       x0[3], x1[3], nnor[3], snor[3], *norT, *souT;
    const double *north, *south;

    *esurf = NULL;
    if (context == NULL)               return EGADS_NULLOBJ;
    if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
    if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
    outLevel = EG_outLevel(context);
    imax     = seq->ncp;
    vdata    = (int *) &vknot[jmax];

    /* check for degenerate sides */
    north = south = NULL;
    norT  = souT  = NULL;
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
        if (outLevel > 0)
            printf(" EGADS Error: Imin (west) degenerate (EG_ loft2spline)!\n");
        return EGADS_DEGEN;
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
        printf("  loft2spline: Jmin (south) degenerate!\n");
#endif
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
                if (outLevel > 0)
                    printf(" EGADS Error: BAD South Axes (EG_loft2spline)!\n");
                return EGADS_NOTORTHO;
            }
            dist   = 1.0/sqrt(dist);
            snor[0] *= dist;
            snor[1] *= dist;
            snor[2] *= dist;
            souT     = snor;
#ifdef DEBUG
            printf("               with normal = %lf %lf %lf!\n",
                   snor[0], snor[1], snor[2]);
#endif
        }
    } else {
        if (drnd[1] != NULL)
            if (drnd[1][0] != 0.0) {
                snor[0] = drnd[1][1]*drnd[1][0];
                snor[1] = drnd[1][2]*drnd[1][0];
                snor[2] = drnd[1][3]*drnd[1][0];
                souT    = snor;
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
        if (outLevel > 0)
            printf(" EGADS Error: Imax (east) degenerate (EG_ loft2spline)!\n");
        return EGADS_DEGEN;
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
        printf("  loft2spline: Jmax (north) degenerate!\n");
#endif
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
                if (outLevel > 0)
                    printf(" EGADS Error: BAD North Axes (EG_spline2d)!\n");
                return EGADS_NOTORTHO;
            }
            dist   = 1.0/sqrt(dist);
            nnor[0] *= dist;
            nnor[1] *= dist;
            nnor[2] *= dist;
            norT     = nnor;
#ifdef DEBUG
            printf("               with normal = %lf %lf %lf!\n",
                   nnor[0], nnor[1], nnor[2]);
#endif
        }
    } else {
        if (drnd[3] != NULL)
            if (drnd[3][0] != 0.0) {
                nnor[0] = drnd[3][1]*drnd[3][0];
                nnor[1] = drnd[3][2]*drnd[3][0];
                nnor[2] = drnd[3][3]*drnd[3][0];
                norT    = nnor;
            }
    }

    if (imax == 3) {
        return EG_spline2dAppx(context, 1, seq->knots, vknot, vdata,
                               NULL, NULL, south, souT, north, norT,
                               imax, jmax, xyz, tol, esurf);
    } else {
        /* blend2xy (x and/or y is "b") fails with 2 */
        return EG_spline2dAppx(context, 1, seq->knots, vknot, vdata,
                               drnd[0], drnd[2], south, souT, north, norT,
                               imax, jmax, xyz, tol, esurf);
    }
}


static void
EG_checkDirs(ego edge, double t, double *pnt, double *dir2)
{
  int    stat;
  double dot, data[9], tan[3], dir[3];

  stat = EG_evaluate(edge, &t, data);
  if (stat != EGADS_SUCCESS) return;
  dir[0]   = (data[0]-pnt[0]);
  dir[1]   = (data[1]-pnt[1]);
  dir[2]   = (data[2]-pnt[2]);
  dot      = sqrt(DOT(dir, dir));
  if (dot != 0.0) {
    dir[0] /= dot;
    dir[1] /= dot;
    dir[2] /= dot;
  }
  tan[0]    = dir2[0];
  tan[1]    = dir2[1];
  tan[2]    = dir2[2];
  dot       = sqrt(DOT(dir2, dir2));
  if (dot != 0.0) {
    tan[0] /= dot;
    tan[1] /= dot;
    tan[2] /= dot;
  }
  if (DOT(dir, tan) < 0.0) {
//  printf(" checkDir flip = %lf!\n", DOT(dir, tan));
    dir2[0] = -dir2[0];
    dir2[1] = -dir2[1];
    dir2[2] = -dir2[2];
  }
}


static void
EG_cleanupTip(ego *objs, int filled)
{
  int    i, j, n, status, mtype, oclass, *senses;
  double range[4];
  ego    surf, geom, loop, pcrv[3], *loops, *edges;
  
  for (i = 0; i < filled; i++) {
    status = EG_getTopology(objs[2*i+1], &surf, &oclass, &mtype, range, &n,
                            &loops, &senses);
    if (status != EGADS_SUCCESS) continue;
    loop   = loops[0];
    status = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &n, &edges,
                            &senses);
    if (status != EGADS_SUCCESS) continue;
    if (n > 3) n = 3;
    for (j = 0; j < n; j++) pcrv[j] = edges[n+j];
    EG_deleteObject(objs[2*i+1]);
    EG_deleteObject(loop);
    EG_deleteObject(surf);
    for (j = 0; j < n; j++) EG_deleteObject(pcrv[j]);
  }
}


static int
EG_getPCurve(ego context, int *first, int sense, ego *pcurve)
{
  double line[4];
  
  if (*first == 0) {
    line[0] = 0.0;
    line[1] = 0.0;
    line[2] = 1.0;
    line[3] = 0.0;
  } else {
    line[1] = 1.0;
    line[3] = 0.0;
    if (*first == sense) {
      line[0] =  1.0;
      line[2] = -1.0;
    } else {
      line[0] =  0.0;
      line[2] =  1.0;
    }
  }
  
  *first = sense;
  return EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, line, pcurve);
}


static int
EG_wingTip(ego body, int nFace, double *ratios, ego *faces, int *bstrp,
           ego *newBody)
{
  int    i, j, k, m, n, nn, outLevel, oclass, mtype, ftype, nLoop, nEdge, nKnot;
  int    *info, *senses, *lsenses, status, open, sense, te = -1;
  double t, alen, mlen, tol, x0[3], x1[3], norm[3], lnrm[3], limits[4], uv[2];
  double *data, *uKnots, *split, *west, *east, results[18], ratio, frac = 0.1;
  ego    context, geom, rGeom, surf, newLoop, newFace, nBody;
  ego    edgpcv[8], neigbr[3], *objs, *loops, *edges, *nodes, *dum;
  
  outLevel = EG_outLevel(body);
  context  = EG_context(body);
  if (context == NULL) return EGADS_NOTCNTX;
  
  *newBody = NULL;
  /* do we have duplicate Faces? */
  for (i = 0; i < nFace-1; i++)
    for (j = i+1; j < nFace; j++)
      if (faces[i] == faces[j]) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d is the same as Face %d!\n", i+1, j+1);
        return EGADS_EXISTS;
      }

  /* get the storage for replacing the Faces in the Body */
  objs = (ego *) EG_alloc(2*nFace*sizeof(ego));
  if (objs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d egos\n", nFace);
    return EGADS_MALLOC;
  }
  split = (double *) EG_alloc(10*nFace*sizeof(double));
  if (split == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d splits\n", nFace);
    EG_free(objs);
    return EGADS_MALLOC;
  }
  
  /* look at each Face */
  for (i = 0; i < nFace; i++) {
    ratio  = ratios[i];
    status = EG_getTopology(faces[i], &geom, &oclass, &ftype, limits, &nLoop,
                            &loops, &senses);
    if (status != EGADS_SUCCESS) {
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTopology L = %d for Face %d (EG_wingTip)!\n",
               status, i+1);
      return status;
    }
    
    /* get the plane normal */
    status = EG_getGeometry(geom, &oclass, &mtype, &rGeom, &info, &data);
    if (status != EGADS_SUCCESS) {
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: EG_getGeometry = %d for Face %d (EG_wingTip)!\n",
               status, i+1);
      return status;
    }
    if (mtype != PLANE) {
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: Face %d is not Planar (EG_wingTip)!\n", i+1);
      return EGADS_GEOMERR;
    }
    x0[0]    = data[3];
    x0[1]    = data[4];
    x0[2]    = data[5];
    x1[0]    = data[6];
    x1[1]    = data[7];
    x1[2]    = data[8];
    CROSS(norm, x0, x1);
    alen     = sqrt(DOT(norm, norm));
    norm[0] /= alen*ftype;
    norm[1] /= alen*ftype;
    norm[2] /= alen*ftype;
#ifdef DEBUG
    printf(" Face %d: %d  norm = %lf %lf %lf\n",
           i+1, ftype, norm[0], norm[1], norm[2]);
#endif
    if (info != NULL) EG_free(info);
    if (data != NULL) EG_free(data);
    
    /* look at the loops */
    if (nLoop != 1) {
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: Face %d has %d Loops (EG_wingTip)!\n",
               i+1, nLoop);
      return EGADS_GEOMERR;
    }
    status = EG_getTopology(loops[0], &geom, &oclass, &mtype, NULL, &nEdge,
                            &edges, &lsenses);
    if (status != EGADS_SUCCESS) {
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTopology E = %d for Face %d (EG_wingTip)!\n",
               status, i+1);
      return status;
    }
    if ((nEdge != 2) && (nEdge != 3)) {
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: Face %d has %d Edges (EG_wingTip)!\n",
               i+1, nEdge);
      return EGADS_GEOMERR;
    }

    /* find the te Edge */
    if (nEdge == 3) {
      mlen = 1.e100;
      for (j = 0; j < nEdge; j++) {
        status = EG_getTopology(edges[j], &geom, &oclass, &mtype, limits, &n,
                                &nodes, &senses);
        if (status != EGADS_SUCCESS) {
          EG_cleanupTip(objs, i);
          EG_free(split);
          EG_free(objs);
          if (outLevel > 0)
            printf(" EGADS Error: EG_getTopology M = %d for Face %d (EG_wingTip)!\n",
                   status, i+1);
          return status;
        }
        status = EG_arcLength(edges[j], limits[0], limits[1], &alen);
        if (status != EGADS_SUCCESS) {
          EG_cleanupTip(objs, i);
          EG_free(split);
          EG_free(objs);
          if (outLevel > 0)
            printf(" EGADS Error: EG_arcLength = %d for Face %d (EG_wingTip)!\n",
                   status, i+1);
          return status;
        }
        if (alen > mlen) continue;
        mlen = alen;
        te   = j;
      }
    }
#ifdef DEBUG
    printf("  nEdge = %d   te = %d\n", nEdge, te+1);
#endif
    
    /* get the Knots */
    uKnots = NULL;
    tol    = 0.0;
    for (nKnot = j = 0; j < nEdge; j++) {
      edgpcv[j] = edges[j];
      neigbr[j] = NULL;
      /* get neighboring Face */
      status = EG_getBodyTopos(body, edges[j], FACE, &n, &dum);
      if (status == EGADS_SUCCESS) {
        if (n == 2) {
          neigbr[j] = dum[0];
          if (neigbr[j] == faces[i]) neigbr[j] = dum[1];
/*        printf(" %d/%d: %lx -- %lx %lx  %lx\n",
                 i,j, neigbr[j], dum[0], dum[1], faces[i]);  */
        }
        EG_free(dum);
      }
      status = EG_getTopology(edges[j], &geom, &oclass, &mtype, limits, &n,
                              &nodes, &senses);
      if (status != EGADS_SUCCESS) {
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTopology N = %d for Face %d (EG_wingTip)!\n",
                 status, i+1);
        return status;
      }
      EG_getTolerance(edges[j], &t);
      if (t > tol) tol = t;
      status = EG_getGeometry(geom, &oclass, &mtype, &rGeom, &info, &data);
      if (status != EGADS_SUCCESS) {
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: EG_getGeometry = %d for Face %d (EG_wingTip)!\n",
                 status, i+1);
        return status;
      }
#ifdef DEBUG
      printf("  Edge #%d:  %lf %lf  sense = %d:\n",
             j+1, limits[0], limits[1], lsenses[j]);
#endif
      if (j == te) {
        if (info != NULL) EG_free(info);
        if (data != NULL) EG_free(data);
        continue;
      }
      if ((mtype != BSPLINE) || (info == NULL) || (data == NULL)) {
        if (info != NULL) EG_free(info);
        if (data != NULL) EG_free(data);
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: Curve type = %d for Face %d (EG_wingTip)!\n",
                 mtype, i+1);
        return EGADS_GEOMERR;
      }
      if (info[1] != 3) {
        EG_free(info);
        EG_free(data);
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: nDeg = %d for Face %d (EG_wingTip)!\n",
                 info[1], i+1);
        return EGADS_GEOMERR;
      }
      if ((limits[0] != 0.0) || (limits[1] != 1.0)) {
        EG_free(info);
        EG_free(data);
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: Bad Range %lf %lf - Face %d (EG_wingTip)!\n",
                 limits[0], limits[1], i+1);
        return EGADS_GEOMERR;
      }
      if (nKnot == 0) {
        nKnot  = info[3]-6;
        uKnots = (double *) EG_alloc(nKnot*sizeof(double));
        if (uKnots == NULL) {
          EG_free(info);
          EG_free(data);
          EG_cleanupTip(objs, i);
          EG_free(split);
          EG_free(objs);
          if (outLevel > 0)
            printf(" EGADS Error: Knot Alloc %d Face %d (EG_wingTip)!\n",
                   info[3], i+1);
          return EGADS_MALLOC;
        }
        for (k = 0; k < nKnot; k++) uKnots[k] = data[k+3];
        sense = lsenses[j];
      } else {
        if (info[3]-6 != nKnot) {
          EG_free(uKnots);
          EG_free(info);
          EG_free(data);
          EG_cleanupTip(objs, i);
          EG_free(split);
          EG_free(objs);
          if (outLevel > 0)
            printf(" EGADS Error: nKnot = %d %d for Face %d (EG_wingTip)!\n",
                   nKnot, info[3]-6, i+1);
          return EGADS_GEOMERR;
        }
/*@-nullderef@*/
        if (lsenses[j] == sense) {
          for (k = 0; k < nKnot; k++) uKnots[nKnot-k-1] += 1.0 - data[k+3];
        } else {
          for (k = 0; k < nKnot; k++) uKnots[k] += data[k+3];
        }
        for (k = 0; k < nKnot; k++) uKnots[k] /= 2.0;
/*@+nullderef@*/
      }
      EG_free(info);
      EG_free(data);
    }
    
    /* make the grid of points - nKnots by 3/5 */
#ifdef SHARPER
    data = (double *) EG_alloc(5*3*nKnot*sizeof(double));
#else
    data = (double *) EG_alloc(7*3*nKnot*sizeof(double));
#endif
    if (data == NULL) {
      EG_free(uKnots);
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: %d x 3 Allocate for Face %d (EG_wingTip)!\n",
               nKnot, i+1);
      return EGADS_MALLOC;
    }
#ifdef SHARPER
    west = &data[3*3*nKnot];
    east = &data[4*3*nKnot];
#else
    west = &data[5*3*nKnot];
    east = &data[6*3*nKnot];
#endif
    for (n = j = 0; j < nEdge; j++) {
      if (j == te) continue;
      status = EG_getTopology(edges[j], &geom, &oclass, &mtype, limits, &nn,
                              &nodes, &senses);
      if (status != EGADS_SUCCESS) {
        EG_free(data);
        EG_free(uKnots);
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTopology N = %d for Face %d (EG_wingTip)!\n",
                 nKnot, i+1);
        return status;
      }
      for (k = 0; k < nKnot; k++) {
/*@-nullderef@*/
        t = uKnots[k];
        if ((n != 0) && (lsenses[j] == sense)) t = 1.0 - uKnots[k];
/*@+nullderef@*/
        if ((k == 0) || (k == nKnot-1)) {
          open = 0;
          if ((t > 0.5) && (nn > 1)) open = 1;
          status = EG_getTopology(nodes[open], &geom, &oclass, &mtype, results,
                                  &m, &dum, &senses);
        } else {
          status = EG_evaluate(edges[j], &t, results);
        }
        if (status != EGADS_SUCCESS) {
          EG_free(data);
          EG_free(uKnots);
          EG_cleanupTip(objs, i);
          EG_free(split);
          EG_free(objs);
          if (outLevel > 0)
            printf(" EGADS Error: %d/%d Eval = %d for Face %d (EG_wingTip)!\n",
                   k, nKnot, status, i+1);
          return status;
        }
        data[(n*nKnot+k)*3  ] = results[0];
        data[(n*nKnot+k)*3+1] = results[1];
        data[(n*nKnot+k)*3+2] = results[2];
        if (n == 0) {
          west[3*k  ] =  norm[0];
          west[3*k+1] =  norm[1];
          west[3*k+2] =  norm[2];
        } else {
          east[3*k  ] = -norm[0];
          east[3*k+1] = -norm[1];
          east[3*k+2] = -norm[2];
        }
        status = EG_getEdgeUV(neigbr[j], edges[j], -lsenses[j], t, uv);
        if (status != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: %d/%d getEdgeUV = %d for Face %d (EG_wingTip)!\n",
                    k, nKnot, status, i+1);
          continue;
        }
        status = EG_evaluate(neigbr[j], uv, results);
        if (status != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: %d/%d EvalF = %d for Face %d (EG_wingTip)!\n",
                    k, nKnot, status, i+1);
          continue;
        }
        x1[0] = results[6];
        x1[1] = results[7];
        x1[2] = results[8];
        alen  = sqrt(DOT(x1,x1));
        if (alen != 0.0) {
          x1[0] /= alen;
          x1[1] /= alen;
          x1[2] /= alen;
        }
        alen  = 1.0;
        if (n == 0) {
          if (DOT(x1, norm) < 0.0) alen = -1.0;
          west[3*k  ] = alen*x1[0];
          west[3*k+1] = alen*x1[1];
          west[3*k+2] = alen*x1[2];
        } else {
          if (DOT(x1, norm) > 0.0) alen = -1.0;
          east[3*k  ] = alen*x1[0];
          east[3*k+1] = alen*x1[1];
          east[3*k+2] = alen*x1[2];
        }
      }
#ifdef SHARPER
      n = 2;
#else
      n = 4;
#endif
    }
    for (open = m = k = 0; k < nKnot; k++) {
      x0[0] =      data[3*k  ] - data[(n*nKnot+k)*3  ];
      x0[1] =      data[3*k+1] - data[(n*nKnot+k)*3+1];
      x0[2] =      data[3*k+2] - data[(n*nKnot+k)*3+2];
      x1[0] = 0.5*(data[3*k  ] + data[(n*nKnot+k)*3  ]);
      x1[1] = 0.5*(data[3*k+1] + data[(n*nKnot+k)*3+1]);
      x1[2] = 0.5*(data[3*k+2] + data[(n*nKnot+k)*3+2]);
      alen  = sqrt(DOT(x0,x0));
      if (alen != 0.0) {
        x0[0] /= alen;
        x0[1] /= alen;
        x0[2] /= alen;
      }
      lnrm[0] = west[3*k  ] - east[3*k  ];
      lnrm[1] = west[3*k+1] - east[3*k+1];
      lnrm[2] = west[3*k+2] - east[3*k+2];
      t       = sqrt(DOT(lnrm,lnrm));
      if (t != 0.0) {
        lnrm[0] /= t;
        lnrm[1] /= t;
        lnrm[2] /= t;
      }
      if (DOT(lnrm,norm) < 0.0) {
        lnrm[0] = -lnrm[0];
        lnrm[1] = -lnrm[1];
        lnrm[2] = -lnrm[2];
      }
      if ((k >= nKnot/2) && (m == 0)) {
        split[10*i  ] = x1[0];
        split[10*i+1] = x1[1];
        split[10*i+2] = x1[2];
        split[10*i+3] = lnrm[0];
        split[10*i+4] = lnrm[1];
        split[10*i+5] = lnrm[2];
        split[10*i+6] = x0[0];
        split[10*i+7] = x0[1];
        split[10*i+8] = x0[2];
        split[10*i+9] = alen*ratio;
        m++;
      }
      if ((k == 0)       && (alen > tol)) open = -1;
      if ((k == nKnot-1) && (alen > tol)) open =  1;
      alen *= 0.5*ratio;
      if (n == 2) {
        if ((k == 0) || (k == nKnot-1)) alen = 0.0;
        data[(  nKnot+k)*3  ] = x1[0] + alen*lnrm[0];
        data[(  nKnot+k)*3+1] = x1[1] + alen*lnrm[1];
        data[(  nKnot+k)*3+2] = x1[2] + alen*lnrm[2];
      } else {
        if ((k == 0)       && (open != -1)) alen = 0.0;
        if ((k == nKnot-1) && (open !=  1)) alen = 0.0;
        mlen = alen;
        if (alen != 0.0) {
          mlen = sqrt(west[3*k  ]*west[3*k  ] + west[3*k+1]*west[3*k+1] +
                      west[3*k+2]*west[3*k+2]);
          mlen = frac*alen/mlen;
        }
        data[(  nKnot+k)*3  ] = data[3*k  ] + mlen*west[3*k  ];
        data[(  nKnot+k)*3+1] = data[3*k+1] + mlen*west[3*k+1];
        data[(  nKnot+k)*3+2] = data[3*k+2] + mlen*west[3*k+2];
        data[(2*nKnot+k)*3  ] = x1[0] + alen*lnrm[0];
        data[(2*nKnot+k)*3+1] = x1[1] + alen*lnrm[1];
        data[(2*nKnot+k)*3+2] = x1[2] + alen*lnrm[2];
        if (alen != 0.0) {
          mlen = sqrt(east[3*k  ]*east[3*k  ] + east[3*k+1]*east[3*k+1] +
                      east[3*k+2]*east[3*k+2]);
          mlen = -frac*alen/mlen;
        }
        data[(3*nKnot+k)*3  ] = data[(n*nKnot+k)*3  ] + mlen*east[3*k  ];
        data[(3*nKnot+k)*3+1] = data[(n*nKnot+k)*3+1] + mlen*east[3*k+1];
        data[(3*nKnot+k)*3+2] = data[(n*nKnot+k)*3+2] + mlen*east[3*k+2];
      }
    }
    
    /* correct points for open TE */
    if (open != 0) {
                      j = 6;
      if (nKnot < 24) j = 3;
      if (nKnot < 12) j = 2;
      for (m = 0; m < j; m++) {
        k = m;
        if (open == 1) k = nKnot - m - 1;
        x1[0] = 0.5*(data[3*k  ] + data[(n*nKnot+k)*3  ]);
        x1[1] = 0.5*(data[3*k+1] + data[(n*nKnot+k)*3+1]);
        x1[2] = 0.5*(data[3*k+2] + data[(n*nKnot+k)*3+2]);
        t     = 0.5*PI*m;
        alen  = sin(t/j);
        if (n == 2) {
          data[(  nKnot+k)*3  ] = (1.0-alen)*x1[0] + alen*data[(  nKnot+k)*3  ];
          data[(  nKnot+k)*3+1] = (1.0-alen)*x1[1] + alen*data[(  nKnot+k)*3+1];
          data[(  nKnot+k)*3+2] = (1.0-alen)*x1[2] + alen*data[(  nKnot+k)*3+2];
        } else {
          x0[0] = (1.0-frac)*data[3*k  ] + frac*data[(n*nKnot+k)*3  ];
          x0[1] = (1.0-frac)*data[3*k+1] + frac*data[(n*nKnot+k)*3+1];
          x0[2] = (1.0-frac)*data[3*k+2] + frac*data[(n*nKnot+k)*3+2];
          data[(  nKnot+k)*3  ] = (1.0-alen)*x0[0] + alen*data[(  nKnot+k)*3  ];
          data[(  nKnot+k)*3+1] = (1.0-alen)*x0[1] + alen*data[(  nKnot+k)*3+1];
          data[(  nKnot+k)*3+2] = (1.0-alen)*x0[2] + alen*data[(  nKnot+k)*3+2];
          data[(2*nKnot+k)*3  ] = (1.0-alen)*x1[0] + alen*data[(2*nKnot+k)*3  ];
          data[(2*nKnot+k)*3+1] = (1.0-alen)*x1[1] + alen*data[(2*nKnot+k)*3+1];
          data[(2*nKnot+k)*3+2] = (1.0-alen)*x1[2] + alen*data[(2*nKnot+k)*3+2];
          x0[0] = frac*data[3*k  ] + (1.0-frac)*data[(n*nKnot+k)*3  ];
          x0[1] = frac*data[3*k+1] + (1.0-frac)*data[(n*nKnot+k)*3+1];
          x0[2] = frac*data[3*k+2] + (1.0-frac)*data[(n*nKnot+k)*3+2];
          data[(3*nKnot+k)*3  ] = (1.0-alen)*x0[0] + alen*data[(3*nKnot+k)*3  ];
          data[(3*nKnot+k)*3+1] = (1.0-alen)*x0[1] + alen*data[(3*nKnot+k)*3+1];
          data[(3*nKnot+k)*3+2] = (1.0-alen)*x0[2] + alen*data[(3*nKnot+k)*3+2];
        }
      }
    }
    
    /* fit the points */
    status = EG_spline2dAppx(context, 0, uKnots, NULL, NULL, west, east, NULL,
                             NULL, NULL, NULL, nKnot, n+1, data, 1.e-8, &surf);
    if (status != EGADS_SUCCESS) {
      EG_free(data);
      EG_free(uKnots);
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: EG_spline2dAppx = %d for Face %d (EG_wingTip)!\n",
               status, i+1);
      return status;
    }
    EG_free(data);
    EG_free(uKnots);
    
    /* get the PCurves */
    for (k = j = 0; j < nEdge; j++) {
      if (j != te) {
        status = EG_getPCurve(context, &k, lsenses[j], &edgpcv[nEdge+j]);
      } else {
        status = EG_otherCurve(surf, edges[j], 0.0, &edgpcv[nEdge+j]);
      }
      if (status != EGADS_SUCCESS) {
        for (k = 0; k < j; k++) EG_deleteObject(edgpcv[nEdge+k]);
        EG_deleteObject(surf);
        EG_cleanupTip(objs, i);
        EG_free(split);
        EG_free(objs);
        if (outLevel > 0)
          printf(" EGADS Error: get PCurve = %d for Face %d (EG_wingTip)!\n",
                 status, i+1);
        return status;
      }
    }

    /* make the new Loop */
    status = EG_makeTopology(context, surf, LOOP, CLOSED, NULL, nEdge, edgpcv,
                             lsenses, &newLoop);
    if (status != EGADS_SUCCESS) {
      for (k = 0; k < nEdge; k++) EG_deleteObject(edgpcv[nEdge+k]);
      EG_deleteObject(surf);
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: EG_makeTopology L = %d for Face %d (EG_wingTip)\n",
               status, i+1);
      return status;
    }
    
    /* get Face orientation */
    uv[0]  = uv[1] = 0.5;
    status = EG_evaluate(surf, uv, results);
    if (status != EGADS_SUCCESS)
      printf(" EGADS Warning: EG_evaluate = %d (EG_wingTip)!", status);
    x0[0]  = results[3];
    x0[1]  = results[4];
    x0[2]  = results[5];
    x1[0]  = results[6];
    x1[1]  = results[7];
    x1[2]  = results[8];
    CROSS(results, x0, x1);
    alen   = sqrt(DOT(results, results));
    results[0] /= alen;
    results[1] /= alen;
    results[2] /= alen;
    mtype = SFORWARD;
    if (DOT(norm, results) < 0.0) mtype = SREVERSE;
#ifdef DEBUG
    printf("  newFace orientation = %lf  %d!\n", DOT(norm, results), mtype);
#endif
    
    /* make the new Face */
    status = EG_makeTopology(context, surf, FACE, mtype, NULL, 1, &newLoop,
                             &mtype, &newFace);
    if (status != EGADS_SUCCESS) {
      EG_deleteObject(newLoop);
      for (k = 0; k < nEdge; k++) EG_deleteObject(edgpcv[nEdge+k]);
      EG_deleteObject(surf);
      EG_cleanupTip(objs, i);
      EG_free(split);
      EG_free(objs);
      if (outLevel > 0)
        printf(" EGADS Error: EG_makeTopology F = %d for Face %d (EG_wingTip)\n",
               status, i+1);
      return status;
    }
    
    /* transfer the attributes */
    status = EG_attributeDup(faces[i], newFace);
    if (status != EGADS_SUCCESS)
      if (outLevel > 0)
        printf(" EGADS Warning: EG_attributeDup = %d for Face %d (EG_wingTip)\n",
               status, i+1);
    status = EG_attributeAdd(newFace, ".blendStrip", ATTRINT, 1, &bstrp[i],
                             NULL, NULL);
    if (status != EGADS_SUCCESS)
      if (outLevel > 0)
        printf(" EGADS Warning: Strip %d blendStrip = %d (EG_wingTip)!\n",
               k, status);
  
    objs[2*i  ] = faces[i];
    objs[2*i+1] = newFace;
  }
  
  /* replace the Faces and make a new Body */
  status = EG_replaceFaces(body, nFace, objs, &nBody);
  EG_cleanupTip(objs, nFace);
  EG_free(objs);
  if (status != EGADS_SUCCESS) {
    EG_free(split);
    if (outLevel > 0)
      printf(" EGADS Error: EG_replaceFaces = %d (EG_wingTip)!\n", status);
    return status;
  }
  *newBody = nBody;
  
  /* splice the wing tips (for better tessellations)
  for (i = 0; i < nFace; i++) {
    status = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, &split[10*i],
                             &geom);
    if (status != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: EG_makeGeometry = %d (EG_wingTip)!\n", status);
      continue;
    }
    limits[0] = -0.1*split[10*i+9];
    limits[1] =  1.1*split[10*i+9];
    limits[2] = -0.6*split[10*i+9];
    limits[3] =  0.6*split[10*i+9];
    status = EG_makeFace(geom, SFORWARD, limits, &newFace);
    if (status != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: EG_makeFace = %d (EG_wingTip)!\n", status);
      EG_deleteObject(geom);
      continue;
    }
    status = EG_intersection(*newBody, newFace, &n, &edges, &rGeom);
    if (status != EGADS_SUCCESS) {
      EG_deleteObject(newFace);
      EG_deleteObject(geom);
      if (outLevel > 0)
        printf(" EGADS Warning: EG_intersection = %d (EG_wingTip)!\n", status);
      continue;
    }
    if (n == 0) continue;
    status = EG_imprintBody(*newBody, n, edges, &nBody);
    EG_deleteObject(rGeom);
    EG_deleteObject(newFace);
    EG_deleteObject(geom);
    EG_free(edges);
    if (status != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: EG_imprintBody = %d (EG_wingTip)!\n", status);
      continue;
    }
    EG_deleteObject(*newBody);
    *newBody = nBody;
  }  */
  EG_free(split);
  
  return EGADS_SUCCESS;
}


int
EG_blend(int nsex, const ego *secs, /*@null@*/ double *rc1,
                                    /*@null@*/ double *rcN, ego *result)
{
  int    i, j, k, n, jj, nn, outLevel, stat, nstripe, oclass, mtype, npt, nchld;
  int    solid, closed, n0, n1, nnode, ncurv, nedge, nface, npcrv;
  int    *senses, *sen, *iinfo, plane[2], planar, lsenses[4] = {1, 1, -1, -1};
  int    nsec, begRC, endRC, rite, left, tips[2], bstrp[2], *vdata;
  double dx, dy, dz, data[18], uv[2], ts[2], t, dt, lims[4], v1[3], v2[3];
  double *rinfo, *xyzs, *t1, *tN, *vknot, ratio[2];
  double d2xdt2L, d2ydt2L, d2zdt2L, d2sdt2L, d2xdt2R, d2ydt2R, d2zdt2R, d2sdt2R;
  ego    context, loop, ref, geom, objs[8], shell, newbody;
  ego    *chldrn, *surfs, *pcurve, *curvs, *nodes, *edges, *loops, *faces, *dum;
  egSequ *ncp;
  const double *sides[4];
  
  nsec    = nsex;
  if (nsec < 0) nsec = -nsec;
  *result = NULL;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL)               return EGADS_NOTCNTX;
  
  /* look at the input and check to see if OK */
  begRC    = endRC  = -1;
  solid    = closed =  1;
  sides[0] = NULL;
  sides[1] = rc1;
  sides[2] = NULL;
  sides[3] = rcN;
  planar   = 0;
  for (nstripe = i = 0; i < nsec; i++) {
    if (secs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Section Object %d (EG_blend)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (secs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is not an EGO (EG_blend)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (secs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d has no data (EG_blend)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(secs[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Context Mismatch (EG_blend)!\n", i+1);
      return EGADS_MIXCNTX;
    }
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d getTopology = %d (EG_blend)!\n",
               i+1, stat);
      return stat;
    }
    if (oclass == NODE) {
      if ((i != 0) && (i != nsec-1)) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Node and not Bound (EG_blend)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      if (i == 0) {
        begRC = 0;
        if (rc1 != NULL) begRC = 1;
      }
      if (i == nsec-1) {
        endRC = 0;
        if (rcN != NULL) endRC = 1;
      }
      loop = NULL;
    } else if (oclass == FACE) {
      if ((i != 0) && (i != nsec-1)) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Face and not Bound (EG_blend)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      if (n != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Face with %d Loops (EG_blend)!\n",
                 i+1, n);
        return EGADS_TOPOERR;
      }
      loop = chldrn[0];
      if (ref->mtype != PLANE) planar = 1;
    } else if (oclass == BODY) {
      if (secs[i]->mtype != WIREBODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Not a WireBody (EG_blend)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      if ((i == 0) || (i == nsec-1)) solid = 0;
      loop = chldrn[0];
      if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
    } else if (oclass == LOOP) {
      if ((i == 0) || (i == nsec-1)) solid = 0;
      loop = secs[i];
      if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is Not a Loop (EG_blend)!\n", i+1);
      return EGADS_NOTTOPO;
    }

    if (loop == NULL) continue;
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &jj, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Loop getTopology = %d (EG_blend)!\n",
               i+1, stat);
      return stat;
    }
    if (mtype == OPEN) closed = solid = 0;
    for (n = j = 0; j < jj; j++) {
      stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, data, &k, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Edge %d getTopo = %d (EG_blend)!\n",
                 i+1, j+1, stat);
        return stat;
      }
      if (mtype != DEGENERATE) n++;
    }
    if (nstripe == 0) {
      nstripe = n;
    } else {
      if (n != nstripe) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d has %d Edges -- prev = %d (EG_blend)!\n",
                 i+1, n, nstripe);
        return EGADS_TOPOERR;
      }
    }
  }
  if (nstripe == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Edges found (EG_blend)!\n");
    return EGADS_TOPOERR;    
  }
  nface = nstripe;
  if (((begRC == 1) || (endRC == 1)) && (nsec <= 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 2 for Nose Treatment (EG_blend)!\n");
    return EGADS_GEOMERR;
  }
  if ((begRC == 1) && (endRC == 1) && (nsec <= 3)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 3 for 2Nose Treatment (EG_blend)!\n");
    return EGADS_GEOMERR;
  }
  
  /* set the knots in the loft direction */
  vknot = (double *) EG_alloc(nsec*(sizeof(double)+sizeof(int)));
  if (vknot == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Knots (EG_blend)!\n", nsec);
    return EGADS_MALLOC;
  }
  vdata = (int *) &vknot[nsec];
  xyzs  = (double *) EG_alloc(3*(3*nstripe+1-closed)*nsec*sizeof(double));
  if (xyzs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %dX%d Points (EG_blend)!\n",
             nstripe, nsec);
    EG_free(vknot);
    return EGADS_MALLOC;
  }
  for (npt = j = 0; j < nsec; j++) {
    stat = EG_getTopology(secs[j], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      for (k = 0; k < 3*nstripe+1-closed; k++, npt++) {
        xyzs[3*npt  ] = data[0];
        xyzs[3*npt+1] = data[1];
        xyzs[3*npt+2] = data[2];
      }
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      loop = chldrn[0];
    } else {
      loop = secs[j];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (n1 = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d GetTOPO = %d (EG_blend)!\n",
                 i+1, j+1, stat);
        EG_free(xyzs);
        EG_free(vknot);
        return stat;
      }
      if (mtype == DEGENERATE) continue;
      stat = EG_getRange(edges[jj], ts, &k);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_blend)!\n",
                 i+1, j+1, stat);
        EG_free(xyzs);
        EG_free(vknot);
        return stat;
      }
      n0 = 1;
      if ((closed == 0) && (n1 == 0)) {
        n1++;
        n0 = 0;
      }
      for (k = n0; k < 4; k++, npt++) {
        if (senses[jj] == 1) {
          t = ts[0] + k*(ts[1] - ts[0])/3.0;
        } else {
          t = ts[1] - k*(ts[1] - ts[0])/3.0;
        }
        stat = EG_evaluate(edges[jj], &t, data);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d eval = %d (EG_blend)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          EG_free(vknot);
          return stat;
        }
        xyzs[3*npt  ] =  data[0];
        xyzs[3*npt+1] =  data[1];
        xyzs[3*npt+2] =  data[2];
      }
    }
  }
  n0 = 3*nstripe+1-closed;
  /* arc-length spaced */
  for (j = 0; j < nsec; j++) vknot[j] = 0.0;
  dz = n0;
  for (i = 0; i < n0; i++) {
    dy = 0.0;
    for (j = 1; j < nsec; j++) {
      dy += sqrt((xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ])*
                 (xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ]) +
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1])*
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1]) +
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2])*
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2]));
    }
    if (dy == 0.0) {
      dz -= 1.0;
      continue;
    }
    dx = 0.0;
    for (j = 1; j < nsec; j++) {
      dx += sqrt((xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ])*
                 (xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ]) +
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1])*
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1]) +
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2])*
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2]))/dy;
      vknot[j] += dx;
    }
  }
  for (n = j = 1; j < nsec; j++)
    if (vknot[j] < vknot[j-1]) n++;
  if (n == 1) {
    for (j = 0; j < nsec; j++) vknot[j] /= dz;
  } else {
    /* equally spaced */
    for (jj = j = 0; j < nsec-1; j++)
      if (secs[j] != secs[j+1]) jj++;
    for (k = j = 0; j < nsec; j++) {
      dy       = k;
      vknot[j] = dy/jj;
      if (secs[j] != secs[j+1]) k++;
    }
  }

  vdata[0] = vdata[nsec-1] = 1;
  for (j = 1; j < nsec-1; j++)
    if ((vknot[j] == vknot[j+1]) && (vknot[j] == vknot[j+2])) {
      if (vknot[j] == vknot[j+3]) {
        printf("EG_blend: Multiplicity > 3 at %d\n", j);
        EG_free(xyzs);
        EG_free(vknot);
        return EGADS_DEGEN;
      }
#ifdef DEBUG
      printf("repeated vknot at j=%d, %d, and %d\n", j, j+1, j+2);
#endif
      vdata[j++] = -3;
      vdata[j++] =  1;
      vdata[j  ] = +3;
    } else if (vknot[j] == vknot[j+1]) {
#ifdef DEBUG
      printf("repeated vknot at j=%d and %d\n", j, j+1);
#endif
      vdata[j++] = -2;
      vdata[j  ] = +2;
    } else {
      vdata[j  ] =  1;
    }
#ifdef DEBUG
  printf("after pass 1\n");
  for (j = 0; j < nsec; j++)
    printf("vdata[%2d]=%2d   %lf\n", j, vdata[j], vknot[j]);
#endif
  
  /* for all multiplicity 2 knots, determine where the flat spot is */
  for (j = 1; j < nsec-1; j++)
    if ((vdata[j] == -2) && (vdata[j+1] == +2)) {
      if (j < 2) {
        /* flat on left */
        if (begRC == 1) {
          printf("EG_blend: C1 -- too close to nose %d\n", j);
          EG_free(xyzs);
          EG_free(vknot);
          return EGADS_DEGEN;
        }
        vdata[j+1] = 1;
      } else if (j > nsec-4) {
        /* flat on right */
        if (endRC == 1) {
          printf("EG_blend: C1 -- too close to rounded tail %d\n", j);
          EG_free(xyzs);
          EG_free(vknot);
          return EGADS_DEGEN;
        }
        vdata[j  ] = 1;
      } else if ((vknot[j-1] == vknot[j-2]) && (vknot[j+3] == vknot[j+2])) {
        /* flat on both sides */
        printf("EG_blend: C1 -- Too few sections on either side of %d and %d\n",
               j, j+1);
        EG_free(xyzs);
        EG_free(vknot);
        return EGADS_DEGEN;
      } else if (vknot[j-1] == vknot[j-2]) {
        /* flat on left because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on left at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j+1] = 1;
      } else if (vknot[j+3] == vknot[j+2]) {
        /* flat on right because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on right at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j  ] = 1;
      } else {
        /* find second derivatives on both sides */
        for (left = rite = i = 0; i < n0; i++) {
          d2xdt2L = ((xyzs[3*(i+(j  )*n0)  ]-
                      xyzs[3*(i+(j-1)*n0)  ])/(vknot[j  ]-vknot[j-1])
                  -  (xyzs[3*(i+(j-1)*n0)  ]-
                      xyzs[3*(i+(j-2)*n0)  ])/(vknot[j-1]-vknot[j-2]))
                  / (vknot[j  ]-vknot[j-2]);
          d2ydt2L = ((xyzs[3*(i+(j  )*n0)+1]-
                      xyzs[3*(i+(j-1)*n0)+1])/(vknot[j  ]-vknot[j-1])
                  -  (xyzs[3*(i+(j-1)*n0)+1]-
                      xyzs[3*(i+(j-2)*n0)+1])/(vknot[j-1]-vknot[j-2]))
                  / (vknot[j  ]-vknot[j-2]);
          d2zdt2L = ((xyzs[3*(i+(j  )*n0)+2]-
                      xyzs[3*(i+(j-1)*n0)+2])/(vknot[j  ]-vknot[j-1])
                  -  (xyzs[3*(i+(j-1)*n0)+2]-
                      xyzs[3*(i+(j-2)*n0)+2])/(vknot[j-1]-vknot[j-2]))
                  / (vknot[j  ]-vknot[j-2]);
          
          d2xdt2R = ((xyzs[3*(i+(j+3)*n0)  ]-
                      xyzs[3*(i+(j+2)*n0)  ])/(vknot[j+3]-vknot[j+2])
                  -  (xyzs[3*(i+(j+2)*n0)  ]-
                      xyzs[3*(i+(j+1)*n0)  ])/(vknot[j+2]-vknot[j+1]))
                  / (vknot[j+3]-vknot[j+1]);
          d2ydt2R = ((xyzs[3*(i+(j+3)*n0)+1]-
                      xyzs[3*(i+(j+2)*n0)+1])/(vknot[j+3]-vknot[j+2])
                  -  (xyzs[3*(i+(j+2)*n0)+1]-
                      xyzs[3*(i+(j+1)*n0)+1])/(vknot[j+2]-vknot[j+1]))
                  / (vknot[j+3]-vknot[j+1]);
          d2zdt2R = ((xyzs[3*(i+(j+3)*n0)+2]-
                      xyzs[3*(i+(j+2)*n0)+2])/(vknot[j+3]-vknot[j+2])
                  -  (xyzs[3*(i+(j+2)*n0)+2]-
                      xyzs[3*(i+(j+1)*n0)+2])/(vknot[j+2]-vknot[j+1]))
                  / (vknot[j+3]-vknot[j+1]);
          
          d2sdt2L = d2xdt2L*d2xdt2L + d2ydt2L*d2ydt2L + d2zdt2L*d2zdt2L;
          d2sdt2R = d2xdt2R*d2xdt2R + d2ydt2R*d2ydt2R + d2zdt2R*d2zdt2R;
          if (d2sdt2L > d2sdt2R) {
            rite++;
          } else {
            left++;
          }
        }
        if (rite > left) {
          /* flat on right because left curvature is smaller */
#ifdef DEBUG
          printf("flat on right %d %d\n", left, rite);
#endif
          vdata[j  ] = 1;
        } else {
          /* flat on left because right curvture is smaller  */
#ifdef DEBUG
          printf("flat on left %d %d\n", left, rite);
#endif
          vdata[j+1] = 1;
        }
      }
      j++;
    }
#ifdef DEBUG
  printf("after pass 2\n");
  for (j = 0; j < nsec; j++)
    printf("vdata[%2d]=%2d\n", j, vdata[j]);
#endif
  EG_free(xyzs);
  
  /* get the number of sample points per Edge */
  stat = EG_allocSeq(nstripe, &ncp);
  if (stat !=  EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Edges (EG_blend)!\n", nstripe);
    EG_free(vknot);
    return stat;
  }
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
      nface++;
    } else if (oclass == BODY) {
      loop = chldrn[0];
    } else {
      loop = secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (j = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTopo = %d (EG_blend)!\n",
                 i+1, j+1, stat);
        EG_free(vknot);
        EG_freeSeq(nstripe, ncp);
        return stat;
      }
      if (mtype == DEGENERATE) continue;
      
      mtype = -1;
      do {
        if ((mtype == BEZIER) || (mtype == BSPLINE)) EG_free(iinfo);
        geom = ref;
        stat = EG_getGeometry(geom, &oclass, &mtype, &ref, &iinfo, &rinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d getGeom = %d (EG_blend)!\n",
                   i+1, j+1, stat);
          EG_free(vknot);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (mtype != BSPLINE) EG_free(rinfo);
      } while ((mtype == TRIMMED) || (mtype == OFFSET));

      stat = EGADS_NOTFOUND;
      if (mtype == LINE) {
        stat = EG_setSeq(j, ncp, 3, NULL, NULL);
      } else if (mtype == CIRCLE) {
/*      k = 16*(data[1]-data[0])/3.1415926;
        if (k < 4) k = 4;
        stat = EG_setSeq(j, ncp,  k, NULL, NULL);  */
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == ELLIPSE) {
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == PARABOLA) {
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == HYPERBOLA) {
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == BEZIER) {
        k    = iinfo[2];
        if (k < 12) k = 12;
        stat = EG_setSeq(j, ncp, k, NULL, NULL);
        EG_free(iinfo);
      } else if (mtype == BSPLINE) {
        stat = EG_setSeq(j, ncp, 12, iinfo, rinfo);
        EG_free(iinfo);
        EG_free(rinfo);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d Edge %d setSeq = %d (EG_blend)!\n",
                 i+1, j+1, stat);
        EG_free(vknot);
        EG_freeSeq(nstripe, ncp);
        return stat;
      }
      j++;
    }
  }
  
  /* make a surface for each patch/stripe */
  surfs = (ego *) EG_alloc(nface*sizeof(ego));
  if (surfs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Surfaces (EG_blend)!\n", nface);
    EG_free(vknot);
    EG_freeSeq(nstripe, ncp);
    return EGADS_MALLOC;
  }
  sen = (int *) EG_alloc(nsec*sizeof(int));
  for (j = 0; j < nstripe; j++) {
    xyzs = (double *) EG_alloc(3*(ncp[j].ncp+2)*nsec*sizeof(double));
    if (xyzs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation for %dX%d points (EG_blend)!\n",
               ncp[j].ncp, nsec);
      for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
      if (sen != NULL) EG_free(sen);
      EG_free(surfs);
      EG_free(vknot);
      EG_freeSeq(nstripe, ncp);
      return EGADS_MALLOC;      
    }
    t1 = &xyzs[3* ncp[j].ncp   *nsec];
    tN = &xyzs[3*(ncp[j].ncp+1)*nsec];
    if (planar == 0) {
      sides[0] = t1;
      sides[2] = tN;
    }
    for (npt = i = 0; i < nsec; i++) {
      if (sen != NULL) sen[i] = 0;
      stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) continue;
      if (oclass == NODE) {
        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          xyzs[3*npt  ] = data[0];
          xyzs[3*npt+1] = data[1];
          xyzs[3*npt+2] = data[2];
        }
        if (i == 0) {
          t1[0] = t1[1] = t1[2] = 0.0;
          tN[0] = tN[1] = tN[2] = 0.0;
        } else {
          t1[3*i  ] = t1[3*i+1] = t1[3*i+2] = 0.0;
          tN[3*i  ] = tN[3*i+1] = tN[3*i+2] = 0.0;
        }
        continue;
      } else if (oclass == FACE) {
        loop = chldrn[0];
      } else if (oclass == BODY) {
        loop = chldrn[0];
      } else {
        loop = secs[i];
      }
      stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                            &senses);
      if (stat != EGADS_SUCCESS) continue;
      for (nn = jj = 0; jj < n; jj++) {
        stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                              &chldrn, &iinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Section %d EDGE %d getTOPO = %d (EG_blend)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          if (sen != NULL) EG_free(sen);
          for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
          EG_free(surfs);
          EG_free(vknot);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (mtype == DEGENERATE) continue;
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, v1, &k,
                              &dum, &iinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d EDGE %d getTOPOn1 = %d (EG_blend)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          if (sen != NULL) EG_free(sen);
          for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
          EG_free(surfs);
          EG_free(vknot);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (nchld == 1) {
          v2[0] = v1[0];
          v2[1] = v1[1];
          v2[2] = v1[2];
        } else {
          stat = EG_getTopology(chldrn[1], &ref, &oclass, &mtype, v2, &k,
                                &dum, &iinfo);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d EDGE %d getTOPOn2 = %d (EG_blend)!\n",
                     i+1, j+1, stat);
            EG_free(xyzs);
            if (sen != NULL) EG_free(sen);
            for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
            EG_free(surfs);
            EG_free(vknot);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
        }
        if (nn == j) {
          stat = EG_getRange(edges[jj], ts, &k);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d getRange = %d (EG_blend)!\n",
                     i+1, j+1, stat);
            EG_free(xyzs);
            if (sen != NULL) EG_free(sen);
            for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
            EG_free(surfs);
            EG_free(vknot);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
          if (sen != NULL) sen[i] = senses[jj];
          dt = ts[1] - ts[0];
          for (k = 0; k < ncp[j].ncp; k++, npt++) {
            if (senses[jj] == 1) {
              t = ts[0] + ncp[j].knots[k]*dt;
            } else {
              t = ts[1] - ncp[j].knots[k]*dt;
            }
            stat = EG_evaluate(edges[jj], &t, data);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Section %d Edge %d eval = %d (EG_blend)!\n",
                       i+1, j+1, stat);
              EG_free(xyzs);
              if (sen != NULL) EG_free(sen);
              for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
              EG_free(surfs);
              EG_free(vknot);
              EG_freeSeq(nstripe, ncp);
              return stat;
            }
            xyzs[3*npt  ] = data[0];
            xyzs[3*npt+1] = data[1];
            xyzs[3*npt+2] = data[2];
            if (k == 0) {
              if (senses[jj] == 1) {
                t1[3*i  ]     =  data[3]*dt;
                t1[3*i+1]     =  data[4]*dt;
                t1[3*i+2]     =  data[5]*dt;
                xyzs[3*npt  ] =  v1[0];
                xyzs[3*npt+1] =  v1[1];
                xyzs[3*npt+2] =  v1[2];
                EG_checkDirs(edges[jj], t+0.001*dt, data, &t1[3*i]);
              } else {
                t1[3*i  ]     = -data[3]*dt;
                t1[3*i+1]     = -data[4]*dt;
                t1[3*i+2]     = -data[5]*dt;
                xyzs[3*npt  ] =  v2[0];
                xyzs[3*npt+1] =  v2[1];
                xyzs[3*npt+2] =  v2[2];
                EG_checkDirs(edges[jj], t-0.001*dt, data, &t1[3*i]);
              }
            } else if (k == ncp[j].ncp-1) {
              if (senses[jj] == 1) {
                tN[3*i  ]     = -data[3]*dt;
                tN[3*i+1]     = -data[4]*dt;
                tN[3*i+2]     = -data[5]*dt;
                xyzs[3*npt  ] =  v2[0];
                xyzs[3*npt+1] =  v2[1];
                xyzs[3*npt+2] =  v2[2];
                EG_checkDirs(edges[jj], t-0.001*dt, data, &tN[3*i]);
              } else {
                tN[3*i  ]     =  data[3]*dt;
                tN[3*i+1]     =  data[4]*dt;
                tN[3*i+2]     =  data[5]*dt;
                xyzs[3*npt  ] =  v1[0];
                xyzs[3*npt+1] =  v1[1];
                xyzs[3*npt+2] =  v1[2];
                EG_checkDirs(edges[jj], t+0.001*dt, data, &tN[3*i]);
              }
            }
          }
          break;
        }
        nn++;
      }
    }
    stat = EG_loft2spline(context, vknot, sides, &ncp[j], nsec, xyzs, 1.e-8,
                          &surfs[j]);
    EG_free(xyzs);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Strip %d splined = %d (EG_blend)!\n", j+1, stat);
      if (sen != NULL) EG_free(sen);
      for (i = 0; i < j; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      EG_freeSeq(nstripe, ncp);
      return stat;
    }
    k    = j+1;
    stat = EG_attributeAdd(surfs[j], ".blendStrip", ATTRINT, 1,
                           &k, NULL, NULL);
    if (stat != EGADS_SUCCESS)
      if (outLevel > 0)
        printf(" EGADS Warning: Strip %d blendStrip = %d (EG_blend)!\n",
               k, stat);
    stat = EG_attributeAdd(surfs[j], ".blendSamples", ATTRREAL, ncp[j].ncp,
                           NULL, ncp[j].knots, NULL);
    if (stat != EGADS_SUCCESS)
      if (outLevel > 0)
        printf(" EGADS Warning: Strip %d blendSamples = %d (EG_blend)!\n",
               k, stat);
    if (sen != NULL) {
      stat = EG_attributeAdd(surfs[j], ".blendSenses",  ATTRINT, nsec,
                             sen, NULL, NULL);
      if (stat != EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: Strip %d blendSenses = %d (EG_blend)!\n",
                 k, stat);
    }
  }
  if (sen != NULL) EG_free(sen);
  EG_freeSeq(nstripe, ncp);
  
  /* the cap surfaces */
  plane[0] = plane[1] = 0;
  if (nstripe != nface) {
    k = nstripe;
    EG_getTopology(secs[0], &ref, &oclass, &mtype, data, &n, &chldrn,
                   &senses);
    if (oclass == FACE) {
      stat = EG_copyObject(ref, NULL, &surfs[k]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: copy CapSurf0 = %d (EG_blend)!\n", stat);
        for (i = 0; i < k; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      EG_getGeometry(surfs[k], &oclass, &mtype, &ref, &iinfo, &rinfo);
      if  (mtype == PLANE) plane[0] = 1;
      if ((mtype == BEZIER) || (mtype == BSPLINE)) EG_free(iinfo);
      EG_free(rinfo);
      k++;
    }
    EG_getTopology(secs[nsec-1], &ref, &oclass, &mtype, data, &n, &chldrn,
                   &senses);
    if (oclass == FACE) {
      stat = EG_copyObject(ref, NULL, &surfs[k]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: copy CapSurfN = %d (EG_blend)!\n", stat);
        for (i = 0; i < k; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      EG_getGeometry(surfs[k], &oclass, &mtype, &ref, &iinfo, &rinfo);
      if  (mtype == PLANE) plane[1] = 1;
      if ((mtype == BEZIER) || (mtype == BSPLINE)) EG_free(iinfo);
      EG_free(rinfo);
    }
  }
  
  /* take the patch-matched surfaces and build the topology */
  
  pcurve = (ego *) EG_alloc(4*sizeof(ego));
  if (pcurve == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for PCurves (EG_blend)!\n");
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return EGADS_MALLOC;
  }
  data[0] = 0.0;
  data[1] = 0.0;
  data[2] = 1.0;
  data[3] = 0.0;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurve[0]);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve0 = %d (EG_blend)!\n", stat);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return stat;
  }
  data[0] = 1.0;
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 1.0;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurve[1]);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve1 = %d (EG_blend)!\n", stat);
    EG_deleteObject(pcurve[0]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return stat;
  }
  data[0] = 0.0;
  data[1] = 1.0;
  data[2] = 1.0;
  data[3] = 0.0;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurve[2]);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve2 = %d (EG_blend)!\n", stat);
    for (i = 0; i < 2;     i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return stat;
  }
  data[0] = 0.0;
  data[1] = 0.0;
  data[2] = 0.0;
  data[3] = 1.0;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurve[3]);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve3 = %d (EG_blend)!\n", stat);
    for (i = 0; i < 3;     i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return stat;
  }

  /* make the nodes */
  if ((begRC >= 0) && (endRC >= 0)) {
    nnode = 2;
  } else if ((begRC >= 0) || (endRC >= 0)) {
    nnode =    nstripe + (1-closed) + 1;
  } else {
    nnode = 2*(nstripe + (1-closed));
  }
  nodes = (ego *) EG_alloc(nnode*sizeof(ego));
  if (nodes == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc Error on %d nodes (EG_blend)!\n", nnode);
    for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return EGADS_MALLOC;
  }

  if ((begRC >= 0) && (endRC >= 0)) {
    /* beginning and ends are degenerate */
    stat = EG_copyObject(secs[0], NULL, &nodes[0]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: copyObject = %d (EG_blend)!\n", stat);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    stat = EG_copyObject(secs[nsec-1], NULL, &nodes[1]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: copyObject = %d (EG_blend)!\n", stat);
      EG_deleteObject(nodes[0]);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
  } else if (begRC >= 0) {
    /* beginning is degenerate */
#ifdef DEBUG
    printf("  Node %d: input\n", 0);
#endif
    stat = EG_copyObject(secs[0], NULL, &nodes[0]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: copyObject = %d (EG_blend)!\n", stat);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    for (j = 0; j < nstripe; j++) {
      stat = EG_getRange(surfs[j], lims, &k);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: getRange = %d (EG_blend)!\n", stat);
        for (i = 0; i <= j;    i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      uv[0] = lims[0];
      uv[1] = lims[3];
      stat = EG_evaluate(surfs[j], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i <= j;    i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", j+1, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[j+1]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i <= j;    i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
    if (closed == 0) {
      uv[0] = lims[1];
      uv[1] = lims[3];
      stat = EG_evaluate(surfs[nstripe-1], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i <= nstripe; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;        i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;    i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", nstripe+1, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[nstripe+1]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i <= nstripe; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;        i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;    i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }      
    }
  } else if (endRC >= 0) {
    /* end is degenerate */
    for (j = 0; j < nstripe; j++) {
      stat = EG_getRange(surfs[j], lims, &k);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: getRange = %d (EG_blend)!\n", stat);
        for (i = 0; i < j;     i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      uv[0] = lims[0];
      uv[1] = lims[2];
      stat = EG_evaluate(surfs[j], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i < j;     i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", j, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[j]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i < j;     i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
    if (closed == 0) {
      uv[0] = lims[1];
      uv[1] = lims[2];
      stat = EG_evaluate(surfs[nstripe-1], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i < nstripe; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", nstripe, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[nstripe]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i < nstripe; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
#ifdef DEBUG
    printf("  Node %d: input\n", nnode-1);
#endif
    stat = EG_copyObject(secs[nsec-1], NULL, &nodes[nnode-1]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: copyObject = %d (EG_blend)!\n", stat);
      for (i = 0; i < nnode-1; i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    
  } else {
    /* both ends are "open" */
    for (j = 0; j < nstripe; j++) {
      stat = EG_getRange(surfs[j], lims, &k);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: getRange = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*j;   i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      uv[0] = lims[0];
      uv[1] = lims[2];
      stat = EG_evaluate(surfs[j], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*j;   i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", 2*j, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[2*j]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*j;   i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      uv[0] = lims[0];
      uv[1] = lims[3];
      stat = EG_evaluate(surfs[j], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*j+1; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", 2*j+1, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[2*j+1]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*j+1; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
    if (closed == 0) {
      uv[0] = lims[1];
      uv[1] = lims[2];
      stat = EG_evaluate(surfs[nstripe-1], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*nstripe; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;         i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;     i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", 2*nstripe, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[2*nstripe]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i < 2*nstripe; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;         i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;     i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      uv[0] = lims[1];
      uv[1] = lims[3];
      stat = EG_evaluate(surfs[nstripe-1], uv, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate = %d (EG_blend)!\n", stat);
        for (i = 0; i < nnode-1; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf("  Node %d: %lf %lf %lf\n", nnode-1, data[0], data[1], data[2]);
#endif
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[nnode-1]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology = %d (EG_blend)!\n", stat);
        for (i = 0; i < nnode-1; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }

  }

  /* make the edges */
  ncurv = 0;
  nedge = 3*nstripe + (1-closed);
  edges = (ego *) EG_alloc(nedge*sizeof(ego));
  curvs = (ego *) EG_alloc(nedge*sizeof(ego));
  if ((edges == NULL) || (curvs == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc Error on %d edges (EG_blend)!\n", nedge);
    if (curvs != NULL) EG_free(curvs);
    if (edges != NULL) EG_free(edges);
    for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
    EG_free(nodes);
    for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return EGADS_MALLOC;
  }
  for (j = 0; j < nstripe; j++) {
    stat = EG_getRange(surfs[j], lims, &k);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: getRange = %d (EG_blend)!\n", stat);
      for (i = 0; i < 3*j;   i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    stat = EG_isoCline(surfs[j], 0, lims[0], &curvs[ncurv]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: isoCline = %d (EG_blend)!\n", stat);
      for (i = 0; i < 3*j;   i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    EG_getRange(curvs[ncurv], ts, &k);
    ncurv++;
    if ((begRC >= 0) && (endRC >= 0)) {
      n0 = 0;
      n1 = 1;
    } else if (begRC >= 0) {
      n0 = 0;
      n1 = j+1;
    } else if (endRC >= 0) {
      n0 = j;
      n1 = nstripe + (1-closed);
    } else {
      n0 = 2*j;
      n1 = 2*j+1;
    }
    objs[0] = nodes[n0];
    objs[1] = nodes[n1];
#ifdef DEBUG
    printf(" making Edge %d from %d to %d\n", 3*j, n0, n1);
#endif
    stat = EG_makeTopology(context, curvs[ncurv-1], EDGE, TWONODE, ts, 2, objs,
                           NULL, &edges[3*j  ]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopology Edge0 = %d (EG_blend)!\n", stat);
      for (i = 0; i < 3*j;   i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    if (begRC >= 0) {
      ts[0] = 0.0;
      ts[1] = 1.0;
      stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1, &nodes[n0],
                             NULL, &edges[3*j+1]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology Edge1 = %d (EG_blend)!\n", stat);
        for (i = 0; i < 3*j+1; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf(" making Edge %d from %d to %d (degen)\n", 3*j+1, n0, n0);
#endif
    } else {
      stat = EG_isoCline(surfs[j], 1, lims[2], &curvs[ncurv]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: isoCline = %d (EG_blend)!\n", stat);
        for (i = 0; i < 3*j+1; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      EG_getRange(curvs[ncurv], ts, &k);
      ncurv++;
      if (endRC >= 0) {
        objs[0] = nodes[j];
        if ((j == nstripe-1) && (closed == 1)) {
          objs[1] = nodes[0];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+1, j, 0);
#endif
        } else {
          objs[1] = nodes[j+1];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+1, j, j+1);
#endif
        }
      } else {
        objs[0] = nodes[2*j];
        if ((j == nstripe-1) && (closed == 1)) {
          objs[1] = nodes[0];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+1, 2*j, 0);
#endif
        } else {
          objs[1] = nodes[2*j+2];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+1, 2*j, 2*j+2);
#endif
        }
      }
      stat = EG_makeTopology(context, curvs[ncurv-1], EDGE, TWONODE, ts, 2,
                             objs, NULL, &edges[3*j+1]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology Edge2 = %d (EG_blend)!\n", stat);
        for (i = 0; i < 3*j+1; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
    if (endRC >= 0) {
      ts[0] = 0.0;
      ts[1] = 1.0;
      stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1, &nodes[n1],
                             NULL, &edges[3*j+2]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology Edge3 = %d (EG_blend)!\n", stat);
        for (i = 0; i < 3*j+2; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
#ifdef DEBUG
      printf(" making Edge %d from %d to %d (degen)\n", 3*j+2, n1, n1);
#endif
    } else {
      stat = EG_isoCline(surfs[j], 1, lims[3], &curvs[ncurv]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: isoCline = %d (EG_blend)!\n", stat);
        for (i = 0; i < 3*j+2; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      EG_getRange(curvs[ncurv], ts, &k);
      ncurv++;
      if (begRC >= 0) {
        objs[0] = nodes[j+1];
        if ((j == nstripe-1) && (closed == 1)) {
          objs[1] = nodes[1];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+2, j+1, 1);
#endif
        } else {
          objs[1] = nodes[j+2];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+2, j+1, j+2);
#endif
        }
      } else {
        objs[0] = nodes[2*j+1];
        if ((j == nstripe-1) && (closed == 1)) {
          objs[1] = nodes[1];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+2, 2*j+1, 1);
#endif
        } else {
          objs[1] = nodes[2*j+3];
#ifdef DEBUG
          printf(" making Edge %d from %d to %d\n", 3*j+2, 2*j+1, 2*j+3);
#endif
        }
      }
      stat = EG_makeTopology(context, curvs[ncurv-1], EDGE, TWONODE, ts, 2,
                             objs, NULL, &edges[3*j+2]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology Edge4 = %d (EG_blend)!\n", stat);
        for (i = 0; i < 3*j+2; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
  }
  if (closed == 0) {
    stat = EG_isoCline(surfs[nstripe-1], 0, lims[1], &curvs[ncurv]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: isoCline = %d (EG_blend)!\n", stat);
      for (i = 0; i < nedge-1; i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv;   i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode;   i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    EG_getRange(curvs[ncurv], ts, &k);
    ncurv++;
    if ((begRC >= 0) && (endRC >= 0)) {
      n0 = 0;
      n1 = 1;
    } else if (begRC >= 0) {
      n0 = 0;
      n1 = nstripe+1;
    } else if (endRC >= 0) {
      n0 = nstripe-1 + (1-closed);
      n1 = nstripe   + (1-closed);
    } else {
      n0 = 2*nstripe;
      n1 = 2*nstripe+1;
    }
    objs[0] = nodes[n0];
    objs[1] = nodes[n1];
#ifdef DEBUG
    printf(" making EDGE %d from %d to %d\n", 3*nstripe, n0, n1);
#endif
    stat = EG_makeTopology(context, curvs[ncurv-1], EDGE, TWONODE, ts, 2,
                           objs, NULL, &edges[3*nstripe]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopology Edge5 = %d (EG_blend)!\n", stat);
      for (i = 0; i < nedge-1; i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv;   i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode;   i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;       i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface;   i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
  }
  
  /* make the loops */
  loops = (ego *) EG_alloc(nface*sizeof(ego));
  if (loops == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc Error on %d loops (EG_blend)!\n", nstripe);
    for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
    EG_free(edges);
    for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
    EG_free(curvs);
    for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
    EG_free(nodes);
    for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return EGADS_MALLOC;
  }
  for (j = 0; j < nstripe; j++) {
    if (j != nstripe-1) {
      k = 3*j+3;
    } else {
      if (closed == 1) {
        k = 0;
      } else {
        k = 3*nstripe;
      }
    }
    objs[0] = edges[3*j+1];
    objs[1] = edges[k];
    objs[2] = edges[3*j+2];
    objs[3] = edges[3*j  ];
    objs[4] = pcurve[0];
    objs[5] = pcurve[1];
    objs[6] = pcurve[2];
    objs[7] = pcurve[3];
#ifdef DEBUG
    printf(" making Loop %d from %d %d %d %d\n", j, 3*j+1, k, 3*j+2, 3*j);
#endif
    stat = EG_makeTopology(context, surfs[j], LOOP, CLOSED, NULL, 4, objs,
                           lsenses, &loops[j]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopology Loop = %d (EG_blend)!\n", stat);
      for (i = 0; i < j;     i++) EG_deleteObject(loops[i]);
      EG_free(loops);
      for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < 4;     i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
  }
  /* the cap loops */
  npcrv = 4;
  if (nstripe != nface) {
    k = nstripe;
    EG_getTopology(secs[0], &ref, &oclass, &mtype, data, &n, &chldrn,
                   &senses);
    if (oclass == FACE) {
      if (plane[0] == 0) {
        chldrn = (ego *) EG_reall(pcurve, (npcrv+nstripe)*sizeof(ego));
        if (chldrn == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: reAlloc0 PCurve (EG_blend)!\n");
          for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
          EG_free(loops);
          for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
          EG_free(edges);
          for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
          EG_free(curvs);
          for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
          EG_free(nodes);
          for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
          EG_free(pcurve);
          for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
          EG_free(surfs);
          EG_free(vknot);
          return EGADS_MALLOC;
        }
        pcurve = chldrn;
      }
      chldrn = (ego *) EG_alloc(2*nstripe*sizeof(ego));
      senses = (int *) EG_alloc(  nstripe*sizeof(int));
      if ((chldrn == NULL) || (senses == NULL)) {
        if (chldrn != NULL) EG_free(chldrn);
        if (senses != NULL) EG_free(senses);
        if (outLevel > 0)
          printf(" EGADS Error: Malloc0 PCurve (EG_blend)!\n");
        for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
        EG_free(loops);
        for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return EGADS_MALLOC;
      }
      for (j = 0; j < nstripe; j++) {
        chldrn[nstripe+j] = NULL;
        chldrn[j] = edges[3*j+1];
        senses[j] = 1;
        if (plane[0] == 1) continue;
        EG_getTopology(chldrn[j], &ref, &oclass, &mtype, data, &n, &dum,
                       &iinfo);
        stat = EG_otherCurve(surfs[k], ref, 1.e-7, &pcurve[npcrv]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: otherCurve0 Loop = %d (EG_blend)!\n", stat);
          EG_free(chldrn);
          EG_free(senses);
          for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
          EG_free(loops);
          for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
          EG_free(edges);
          for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
          EG_free(curvs);
          for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
          EG_free(nodes);
          for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
          EG_free(pcurve);
          for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
          EG_free(surfs);
          EG_free(vknot);
          return stat;
        }
        chldrn[nstripe+j] = pcurve[npcrv];
        npcrv++;
      }
      ref = surfs[k];
      if (plane[0] == 1) ref = NULL;
      j = OPEN;
      if (closed == 1) j = CLOSED;
      stat = EG_makeTopology(context, ref, LOOP, j, NULL, nstripe,
                             chldrn, senses, &loops[k]);
      EG_free(chldrn);
      EG_free(senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology0 Loop = %d (EG_blend)!\n", stat);
        for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
        EG_free(loops);
        for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      k++;
    }
    EG_getTopology(secs[nsec-1], &ref, &oclass, &mtype, data, &n, &chldrn,
                   &senses);
    if (oclass == FACE) {
      if (plane[1] == 0) {
        chldrn = (ego *) EG_reall(pcurve, (npcrv+nstripe)*sizeof(ego));
        if (chldrn == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: reAllocN PCurve (EG_blend)!\n");
          for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
          EG_free(loops);
          for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
          EG_free(edges);
          for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
          EG_free(curvs);
          for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
          EG_free(nodes);
          for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
          EG_free(pcurve);
          for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
          EG_free(surfs);
          EG_free(vknot);
          return EGADS_MALLOC;
        }
        pcurve = chldrn;
      }
      chldrn = (ego *) EG_alloc(2*nstripe*sizeof(ego));
      senses = (int *) EG_alloc(  nstripe*sizeof(int));
      if ((chldrn == NULL) || (senses == NULL)) {
        if (chldrn != NULL) EG_free(chldrn);
        if (senses != NULL) EG_free(senses);
        if (outLevel > 0)
          printf(" EGADS Error: MallocN PCurve (EG_blend)!\n");
        for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
        EG_free(loops);
        for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return EGADS_MALLOC;
      }
      for (j = 0; j < nstripe; j++) {
        chldrn[nstripe+j] = NULL;
        chldrn[j] = edges[3*j+2];
        senses[j] = 1;
        if (plane[1] == 1) continue;
        EG_getTopology(chldrn[j], &ref, &oclass, &mtype, data, &n, &dum,
                       &iinfo);
        if (ref == NULL) continue;
        stat = EG_otherCurve(surfs[k], ref, 1.e-7, &pcurve[npcrv]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: otherCurveN Loop = %d (EG_blend)!\n", stat);
          EG_free(chldrn);
          EG_free(senses);
          for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
          EG_free(loops);
          for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
          EG_free(edges);
          for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
          EG_free(curvs);
          for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
          EG_free(nodes);
          for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
          EG_free(pcurve);
          for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
          EG_free(surfs);
          EG_free(vknot);
          return stat;
        }
        chldrn[nstripe+j] = pcurve[npcrv];
        npcrv++;
      }
      ref = surfs[k];
      if (plane[1] == 1) ref = NULL;
      j = OPEN;
      if (closed == 1) j = CLOSED;
      stat = EG_makeTopology(context, ref, LOOP, j, NULL, nstripe,
                             chldrn, senses, &loops[k]);
      EG_free(chldrn);
      EG_free(senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopologyN Loop = %d (EG_blend)!\n", stat);
        for (i = 0; i < k;     i++) EG_deleteObject(loops[i]);
        EG_free(loops);
        for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
    }
  }

  /* make the faces */
  faces = (ego *) EG_alloc(nface*sizeof(ego));
  if (faces == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc Error on %d faces (EG_blend)!\n", nface);
    for (i = 0; i < nface; i++) EG_deleteObject(loops[i]);
    EG_free(loops);
    for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
    EG_free(edges);
    for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
    EG_free(curvs);
    for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
    EG_free(nodes);
    for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return EGADS_MALLOC;
  }
  for (j = 0; j < nstripe; j++) {
    k = SFORWARD;
    /* why is this necessary? */
    if ((nstripe == 1) && (closed == 1) && ((endRC >= 0) || (begRC >= 0)))
      k = SREVERSE;
    stat = EG_makeTopology(context, surfs[j], FACE, k, NULL, 1,
                           &loops[j], lsenses, &faces[j]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopology Face %d = %d (EG_blend)!\n",
               j, stat);
      for (i = 0; i < j;     i++) EG_deleteObject(faces[i]);
      EG_free(faces);
      for (i = 0; i < nface; i++) EG_deleteObject(loops[i]);
      EG_free(loops);
      for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
      EG_free(edges);
      for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
      EG_free(nodes);
      for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
      EG_free(pcurve);
      for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
      EG_free(surfs);
      EG_free(vknot);
      return stat;
    }
    /* copy the surface attributes */
    EG_attributeDup(surfs[j], faces[j]);
  }
  
  /* copy the attributes from the face caps */
  nchld = 0;
  if (nstripe != nface) {
    k = nstripe;
    EG_getTopology(secs[0], &ref, &oclass, &mtype, data, &n, &chldrn, &senses);
    if (oclass == FACE) {
      if ((nstripe == 1) && (closed == 1) && ((endRC >= 0) || (begRC >= 0)))
        mtype = -mtype;
      stat = EG_makeTopology(context, surfs[k], FACE, mtype, NULL, 1,
                             &loops[k], lsenses, &faces[k]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology Face %d = %d (EG_blend)!\n",
                 k, stat);
        for (i = 0; i < k;     i++) EG_deleteObject(faces[i]);
        EG_free(faces);
        for (i = 0; i < nface; i++) EG_deleteObject(loops[i]);
        EG_free(loops);
        for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      EG_attributeDup(secs[0], faces[k]);
      i    = -3;
      stat = EG_attributeAdd(faces[k], ".blendStrip", ATTRINT, 1,
                             &i, NULL, NULL);
      if (stat != EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: Strip %d blendStrip = %d (EG_blend)!\n",
                 i, stat);
      if (rc1 != NULL)
        if (rc1[0] == 0.0) {
          tips[nchld]  = k;
          ratio[nchld] = rc1[1];
          bstrp[nchld] = -1;
          nchld++;
        }
      k++;
    }
    EG_getTopology(secs[nsec-1], &ref, &oclass, &mtype, data, &n, &chldrn,
                   &senses);
    if (oclass == FACE) {
      if ((nstripe == 1) && (closed == 1) && ((endRC >= 0) || (begRC >= 0)))
        mtype = -mtype;
      stat = EG_makeTopology(context, surfs[k], FACE, mtype, NULL, 1,
                             &loops[k], lsenses, &faces[k]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopology Face %d = %d (EG_blend)!\n",
                 k, stat);
        for (i = 0; i < k;     i++) EG_deleteObject(faces[i]);
        EG_free(faces);
        for (i = 0; i < nface; i++) EG_deleteObject(loops[i]);
        EG_free(loops);
        for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
        EG_free(edges);
        for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
        EG_free(curvs);
        for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
        EG_free(nodes);
        for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
        EG_free(pcurve);
        for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
        EG_free(surfs);
        EG_free(vknot);
        return stat;
      }
      EG_attributeDup(secs[nsec-1], faces[k]);
      i    = -4;
      stat = EG_attributeAdd(faces[k], ".blendStrip", ATTRINT, 1,
                             &i, NULL, NULL);
      if (stat != EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: Strip %d blendStrip = %d (EG_blend)!\n",
                 i, stat);
      if (rcN != NULL)
        if (rcN[0] == 0.0) {
          tips[nchld]  = k;
          ratio[nchld] = rcN[1];
          bstrp[nchld] = -2;
          nchld++;
        }
    }
  }

  /* make the shell */
  j = OPEN;
  if (solid == 1) j = CLOSED;
  stat = EG_makeTopology(context, NULL, SHELL, j, NULL, nface, faces, NULL,
                         &shell);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeTopology Shell = %d (EG_blend)!\n", stat);
    for (i = 0; i < nface; i++) EG_deleteObject(faces[i]);
    EG_free(faces);
    for (i = 0; i < nface; i++) EG_deleteObject(loops[i]);
    EG_free(loops);
    for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
    EG_free(edges);
    for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
    EG_free(curvs);
    for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
    EG_free(nodes);
    for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
    EG_free(pcurve);
    for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
    EG_free(surfs);
    EG_free(vknot);
    return stat;
  }
  
  /* make the body */
  j = SHEETBODY;
  if (solid == 1) j = SOLIDBODY;
  stat = EG_makeTopology(context, NULL, BODY, j, NULL, 1, &shell, NULL,
                         result);

  /* free up all of our loose objects */
  EG_deleteObject(shell);
  for (i = 0; i < nface; i++) EG_deleteObject(faces[i]);
  EG_free(faces);
  for (i = 0; i < nface; i++) EG_deleteObject(loops[i]);
  EG_free(loops);
  for (i = 0; i < nedge; i++) EG_deleteObject(edges[i]);
  EG_free(edges);
  for (i = 0; i < ncurv; i++) EG_deleteObject(curvs[i]);
  EG_free(curvs);
  for (i = 0; i < nnode; i++) EG_deleteObject(nodes[i]);
  EG_free(nodes);
  for (i = 0; i < npcrv; i++) EG_deleteObject(pcurve[i]);
  EG_free(pcurve);
  for (i = 0; i < nface; i++) EG_deleteObject(surfs[i]);
  EG_free(surfs);
  
#ifdef DEBUG
  if (*result != NULL) {
    EG_getBodyTopos(*result, NULL, FACE, &nface, &faces);
    for (i = 0; i < nface; i++) {
      EG_getTolerance(faces[i], &t);
      printf(" Face %d tolerance = %le\n", i, t);
    }
    EG_free(faces);
    EG_getBodyTopos(*result, NULL, EDGE, &nedge, &edges);
    for (i = 0; i < nedge; i++) {
      stat = EG_getTolerance(edges[i], &t);
      printf(" Edge %d tolerance = %le\n", i, t);
    }
    EG_free(edges);
    EG_getBodyTopos(*result, NULL, NODE, &nnode, &nodes);
    for (i = 0; i < nnode; i++) {
      stat = EG_getTolerance(nodes[i], &t);
      printf(" Node %d tolerance = %le\n", i, t);
    }
    EG_free(nodes);
  }
#endif
  
  if ((stat != EGADS_SUCCESS) || (*result == NULL)) {
    if (stat == EGADS_SUCCESS) stat = EGADS_NULLOBJ;
    if (outLevel > 0)
      printf(" EGADS Error: makeTopology Body = %d (EG_blend)!\n", stat);
    EG_free(vknot);
    return stat;
  }
  
  /* modify cap Faces with wingTips */
  if (nchld != 0) {
    stat = EG_getBodyTopos(*result, NULL, FACE, &nface, &faces);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Warning: EG_getBodyToposWT = %d (EG_blend)!\n", stat);
      return EGADS_SUCCESS;
    }
    for (i = 0; i < nchld; i++) objs[i] = faces[tips[i]];
    stat = EG_wingTip(*result, nchld, ratio, objs, bstrp, &newbody);
    EG_free(faces);
    if (stat == EGADS_SUCCESS) {
      EG_deleteObject(*result);
      *result = newbody;
    } else {
      printf(" EGADS Warning: EG_wingTip = %d (EG_blend)!\n", stat);
    }
  }
  if (nsex < 0) {
    EG_free(vknot);
    return EGADS_SUCCESS;
  }

  /* look for and scribe C0 sections (multiplicity of 3 or more) */
  /*  for (j = 0; j < nsec; j++) printf(" knot %d: %lf\n", j, vknot[j]);  */
  for (j = 0; j < nsec-2; j++)
    if ((vknot[j] == vknot[j+1]) && (vknot[j+1] == vknot[j+2])) {
      if (j > 1)
        if ((vknot[j-1] == vknot[j]) && (vknot[j] == vknot[j+1])) continue;
/*    printf(" ** C0 at section %d (of %d), nstripe = %d! **\n",
             j+2, nsec, nstripe);  */
      stat = EG_getTopology(secs[j+1], &ref, &oclass, &mtype, data, &n, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Warning: EG_getTopology/Split = %d (EG_blend)!\n",
                 stat);
        continue;
      }
      if (oclass == NODE) {
        continue;
      } else if (oclass == FACE) {
        loop = chldrn[0];
      } else if (oclass == BODY) {
        loop = chldrn[0];
      } else {
        loop = secs[j+1];
      }
      stat = EG_getPlane(loop, &geom);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Warning: EG_getPlane = %d (EG_blend)!\n", stat);
        continue;
      }
      lims[0] = lims[2] = -1.e10;
      lims[1] = lims[3] =  1.e10;
      stat = EG_makeFace(geom, SFORWARD, lims, &ref);
      if ((stat != EGADS_SUCCESS) || (ref == NULL)) {
        if (outLevel > 0)
          printf(" EGADS Warning: EG_makeFace = %d (EG_blend)!\n", stat);
        EG_deleteObject(geom);
        continue;
      }
      stat = EG_intersection(*result, ref, &n, &dum, &shell);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Warning: EG_intersection = %d (EG_blend)!\n", stat);
        EG_deleteObject(ref);
        EG_deleteObject(geom);
        continue;
      }
      if (n == 0) {
        if (outLevel > 0)
          printf(" EGADS Warning: No intersected Faces (EG_blend)!\n");
        EG_deleteObject(ref);
        EG_deleteObject(geom);
        continue;
      }
//    printf(" number of intersected Edges = %d\n", n);
      stat = EG_imprintBody(*result, n, dum, &newbody);
      EG_deleteObject(shell);
      EG_free(dum);
      EG_deleteObject(ref);
      EG_deleteObject(geom);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Warning: EG_imprintBody = %d (EG_blend)!\n", stat);
      } else {
        EG_deleteObject(*result);
        *result = newbody;
      }
    }
  EG_free(vknot);
  
  return EGADS_SUCCESS;
}


int
EG_ruled(int nsec, const ego *secs, ego *result)
{
  int    i, j, k, n, jj, nn, outLevel, stat, nstripe, oclass, mtype, npt, nchld;
  int    closed, nface, ncap, npatch, cap, planar, *senses, *iinfo;
  double data[18], ts[2], t, dt, lims[4], v1[3], v2[3];
  double *rinfo, *xyzs, t1[6], tN[6];
  ego    context, loop, ref, geom, model, body, bLoop, eLoop;
  ego    *chldrn, *curvs, *surfs, *edges, *faces, *nodes, nds[2], *dum;
  egSequ *ncp;
  
  *result = NULL;
  bLoop   = eLoop = NULL;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL)               return EGADS_NOTCNTX;
  
  /* look at the input and check to see if OK */
  closed = 1;
  for (ncap = nstripe = i = 0; i < nsec; i++) {
    if (secs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Section Object %d (EG_ruled)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (secs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is not an EGO (EG_ruled)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (secs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d has no data (EG_ruled)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(secs[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Context Mismatch (EG_ruled)!\n", i+1);
      return EGADS_MIXCNTX;
    }
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d getTopology = %d (EG_ruled)!\n",
               i+1, stat);
      return stat;
    }
    if (oclass == NODE) {
      if ((i != 0) && (i != nsec-1)) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Node and not Bound (EG_ruled)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      loop = NULL;
    } else if (oclass == FACE) {
      if (n != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Face with %d Loops (EG_ruled)!\n",
                 i+1, n);
        return EGADS_TOPOERR;
      }
      loop = chldrn[0];
      if ((i == 0) || (i == nsec-1)) ncap++;
    } else if (oclass == BODY) {
      if (secs[i]->mtype != WIREBODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Not a WireBody (EG_ruled)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      loop = chldrn[0];
    } else if (oclass == LOOP) {
      loop = secs[i];
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is Not a Loop (EG_ruled)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    
    if (loop == NULL) continue;
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &jj, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Loop getTopology = %d (EG_ruled)!\n",
               i+1, stat);
      return stat;
    }
    if (mtype == OPEN) closed = 0;
    for (n = j = 0; j < jj; j++) {
      stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, data, &k, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Edge %d getTopo = %d (EG_ruled)!\n",
                 i+1, j+1, stat);
        return stat;
      }
      if (mtype != DEGENERATE) n++;
    }
    if (nstripe == 0) {
      nstripe = n;
    } else {
      if (n != nstripe) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d has %d Edges -- prev = %d (EG_ruled)!\n",
                 i+1, n, nstripe);
        return EGADS_TOPOERR;
      }
    }
  }
  if (nstripe == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Edges found (EG_ruled)!\n");
    return EGADS_TOPOERR;
  }
  if (closed == 0) ncap = 0;
  nface = nstripe*(nsec-1) + ncap;

  /* get the number of sample points per Edge */
  stat = EG_allocSeq(nstripe, &ncp);
  if (stat !=  EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Edges (EG_ruled)!\n", nstripe);
    return stat;
  }
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      loop = chldrn[0];
    } else {
      loop = secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (j = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTopo = %d (EG_ruled)!\n",
                 i+1, j+1, stat);
        EG_freeSeq(nstripe, ncp);
        return stat;
      }
      if (mtype == DEGENERATE) continue;
      
      mtype = -1;
      do {
        if ((mtype == BEZIER) || (mtype == BSPLINE)) EG_free(iinfo);
        geom = ref;
        stat = EG_getGeometry(geom, &oclass, &mtype, &ref, &iinfo, &rinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d getGeom = %d (EG_ruled)!\n",
                   i+1, j+1, stat);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (mtype != BSPLINE) EG_free(rinfo);
      } while ((mtype == TRIMMED) || (mtype == OFFSET));
      
      stat = EGADS_NOTFOUND;
      if (mtype == LINE) {
        stat = EG_setSeq(j, ncp, 3, NULL, NULL);
      } else if (mtype == CIRCLE) {
/*      k = 16*(data[1]-data[0])/3.1415926;
        if (k < 4) k = 4;
        stat = EG_setSeq(j, ncp,  k, NULL, NULL);  */
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == ELLIPSE) {
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == PARABOLA) {
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == HYPERBOLA) {
        stat = EG_setSeq(j, ncp, 12, NULL, NULL);
      } else if (mtype == BEZIER) {
        k    = iinfo[2];
        if (k < 12) k = 12;
        stat = EG_setSeq(j, ncp, k, NULL, NULL);
        EG_free(iinfo);
      } else if (mtype == BSPLINE) {
        stat = EG_setSeq(j, ncp, 12, iinfo, rinfo);
        EG_free(iinfo);
        EG_free(rinfo);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d Edge %d setSeq = %d (EG_ruled)!\n",
                 i+1, j+1, stat);
        EG_freeSeq(nstripe, ncp);
        return stat;
      }
      j++;
    }
  }
  
  /* make a surface for each patch */
  surfs = (ego *) EG_alloc(nface*sizeof(ego));
  faces = (ego *) EG_alloc(nface*sizeof(ego));
  if ((surfs == NULL) || (faces == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Surfaces (EG_ruled)!\n", nface);
    if (surfs != NULL) EG_free(surfs);
    if (faces != NULL) EG_free(faces);
    EG_freeSeq(nstripe, ncp);
    return EGADS_MALLOC;
  }
  curvs = NULL;
  if (ncap != 0) {
    curvs = (ego *) EG_alloc(2*nstripe*sizeof(ego));
    if (curvs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation for %d Curves (EG_ruled)!\n",
               ncap*nstripe);
      EG_free(surfs);
      EG_free(faces);
      EG_freeSeq(nstripe, ncp);
      return EGADS_MALLOC;
    }
    for (i = 0; i < 2*nstripe; i++) curvs[i] = NULL;
  }
  for (npatch = j = 0; j < nstripe; j++) {
    xyzs = (double *) EG_alloc(6*ncp[j].ncp*sizeof(double));
    if (xyzs == NULL)  {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation for %dX%d points (EG_ruled)!\n",
               ncp[j].ncp, nsec);
      if (xyzs  != NULL) EG_free(xyzs);
      if (curvs != NULL) {
        for (i = 0; i < 2*nstripe; i++)
          if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
        EG_free(curvs);
      }
      for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
      for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
      EG_free(faces);
      EG_free(surfs);
      EG_freeSeq(nstripe, ncp);
      return EGADS_MALLOC;
    }

    for (i = 0; i < nsec-1; i++) {
      npt  = cap = planar = 0;
      
      /* first section */
      stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) continue;
      if (oclass == NODE) {
        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          xyzs[3*npt  ] = data[0];
          xyzs[3*npt+1] = data[1];
          xyzs[3*npt+2] = data[2];
        }
        t1[0] = t1[1] = t1[2] = 0.0;
        tN[0] = tN[1] = tN[2] = 0.0;
        goto next;
      } else if (oclass == FACE) {
        loop = chldrn[0];
        if ((i == 0) && (ncap != 0)) cap = 1;
        if (ref->mtype != PLANE) planar = 1;
      } else if (oclass == BODY) {
        loop = chldrn[0];
        if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
      } else {
        loop = secs[i];
        if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
      }
      stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                            &senses);
      if (stat != EGADS_SUCCESS) continue;
      for (nn = jj = 0; jj < n; jj++) {
        stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                              &chldrn, &iinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Section %d EDGE %d getTOPO = %d (EG_ruled)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          if (curvs != NULL) {
            for (i = 0; i < 2*nstripe; i++)
              if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
            EG_free(curvs);
          }
          for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
          for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
          EG_free(faces);
          EG_free(surfs);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (mtype == DEGENERATE) continue;
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, v1, &k,
                              &dum, &iinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d EDGE %d getTOPOn1 = %d (EG_ruled)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          if (curvs != NULL) {
            for (i = 0; i < 2*nstripe; i++)
              if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
            EG_free(curvs);
          }
          for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
          for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
          EG_free(faces);
          EG_free(surfs);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (nchld == 1) {
          v2[0] = v1[0];
          v2[1] = v1[1];
          v2[2] = v1[2];
        } else {
          stat = EG_getTopology(chldrn[1], &ref, &oclass, &mtype, v2, &k,
                                &dum, &iinfo);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d EDGE %d getTOPOn2 = %d (EG_ruled)!\n",
                     i+1, j+1, stat);
            EG_free(xyzs);
            if (curvs != NULL) {
              for (i = 0; i < 2*nstripe; i++)
                if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
              EG_free(curvs);
            }
            for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
            for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
            EG_free(faces);
            EG_free(surfs);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
        }
        if (nn == j) {
          stat = EG_getRange(edges[jj], ts, &k);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d getRange = %d (EG_ruled)!\n",
                     i+1, j+1, stat);
            EG_free(xyzs);
            if (curvs != NULL) {
              for (i = 0; i < 2*nstripe; i++)
                if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
              EG_free(curvs);
            }
            for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
            for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
            EG_free(faces);
            EG_free(surfs);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
          dt = ts[1] - ts[0];
          for (k = 0; k < ncp[j].ncp; k++, npt++) {
            if (senses[jj] == 1) {
              t = ts[0] + ncp[j].knots[k]*dt;
            } else {
              t = ts[1] - ncp[j].knots[k]*dt;
            }
            stat = EG_evaluate(edges[jj], &t, data);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Section %d Edge %d eval = %d (EG_ruled)!\n",
                       i+1, j+1, stat);
              EG_free(xyzs);
              if (curvs != NULL) {
                for (i = 0; i < 2*nstripe; i++)
                  if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
                EG_free(curvs);
              }
              for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
              for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
              EG_free(faces);
              EG_free(surfs);
              EG_freeSeq(nstripe, ncp);
              return stat;
            }
            xyzs[3*npt  ] = data[0];
            xyzs[3*npt+1] = data[1];
            xyzs[3*npt+2] = data[2];
            if (k == 0) {
              if (senses[jj] == 1) {
                t1[0]         =  data[3]*dt;
                t1[1]         =  data[4]*dt;
                t1[2]         =  data[5]*dt;
                xyzs[3*npt  ] =  v1[0];
                xyzs[3*npt+1] =  v1[1];
                xyzs[3*npt+2] =  v1[2];
                EG_checkDirs(edges[jj], t+0.001*dt, data, &t1[0]);
              } else {
                t1[0]         = -data[3]*dt;
                t1[1]         = -data[4]*dt;
                t1[2]         = -data[5]*dt;
                xyzs[3*npt  ] =  v2[0];
                xyzs[3*npt+1] =  v2[1];
                xyzs[3*npt+2] =  v2[2];
                EG_checkDirs(edges[jj], t-0.001*dt, data, &t1[0]);
              }
            } else if (k == ncp[j].ncp-1) {
              if (senses[jj] == 1) {
                tN[0]         = -data[3]*dt;
                tN[1]         = -data[4]*dt;
                tN[2]         = -data[5]*dt;
                xyzs[3*npt  ] =  v2[0];
                xyzs[3*npt+1] =  v2[1];
                xyzs[3*npt+2] =  v2[2];
                EG_checkDirs(edges[jj], t-0.001*dt, data, &tN[0]);
              } else {
                tN[0]         =  data[3]*dt;
                tN[1]         =  data[4]*dt;
                tN[2]         =  data[5]*dt;
                xyzs[3*npt  ] =  v1[0];
                xyzs[3*npt+1] =  v1[1];
                xyzs[3*npt+2] =  v1[2];
                EG_checkDirs(edges[jj], t+0.001*dt, data, &tN[0]);
              }
            }
          }
          break;
        }
        nn++;
      }

next:
      /* second section */
      stat = EG_getTopology(secs[i+1], &ref, &oclass, &mtype, data, &n, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Info: Section+ %d/%d getTopology = %d (EG_ruled)!\n",
                 i+1, j+1, stat);
        continue;
      }
      if (oclass == NODE) {
        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          xyzs[3*npt  ] = data[0];
          xyzs[3*npt+1] = data[1];
          xyzs[3*npt+2] = data[2];
        }
        t1[3] = t1[4] = t1[5] = 0.0;
        tN[3] = tN[4] = tN[5] = 0.0;
        goto finis;
      } else if (oclass == FACE) {
        loop = chldrn[0];
        if ((i == nsec-2) && (ncap != 0)) cap |= 2;
        if (ref->mtype != PLANE) planar = 1;
      } else if (oclass == BODY) {
        loop = chldrn[0];
        if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
      } else {
        loop = secs[i+1];
        if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
      }
      stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Info: Section+ %d/%d getTopo Loop = %d (EG_ruled)!\n",
                 i+1, j+1, stat);
        continue;
      }
      for (nn = jj = 0; jj < n; jj++) {
        stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                              &chldrn, &iinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Section %d EDGE %d getTOPO = %d (EG_ruled)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          if (curvs != NULL) {
            for (i = 0; i < 2*nstripe; i++)
              if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
            EG_free(curvs);
          }
          for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
          for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
          EG_free(faces);
          EG_free(surfs);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (mtype == DEGENERATE) continue;
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, v1, &k,
                              &dum, &iinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d EDGE %d getTOPOn1 = %d (EG_ruled)!\n",
                   i+1, j+1, stat);
          EG_free(xyzs);
          if (curvs != NULL) {
            for (i = 0; i < 2*nstripe; i++)
              if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
            EG_free(curvs);
          }
          for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
          for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
          EG_free(faces);
          EG_free(surfs);
          EG_freeSeq(nstripe, ncp);
          return stat;
        }
        if (nchld == 1) {
          v2[0] = v1[0];
          v2[1] = v1[1];
          v2[2] = v1[2];
        } else {
          stat = EG_getTopology(chldrn[1], &ref, &oclass, &mtype, v2, &k,
                                &dum, &iinfo);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d EDGE %d getTOPOn2 = %d (EG_ruled)!\n",
                     i+1, j+1, stat);
            EG_free(xyzs);
            if (curvs != NULL) {
              for (i = 0; i < 2*nstripe; i++)
                if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
              EG_free(curvs);
            }
            for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
            for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
            EG_free(faces);
            EG_free(surfs);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
        }
        if (nn == j) {
          stat = EG_getRange(edges[jj], ts, &k);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d getRange = %d (EG_ruled)!\n",
                     i+1, j+1, stat);
            EG_free(xyzs);
            if (curvs != NULL) {
              for (i = 0; i < 2*nstripe; i++)
                if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
              EG_free(curvs);
            }
            for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
            for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
            EG_free(faces);
            EG_free(surfs);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
          dt = ts[1] - ts[0];
          for (k = 0; k < ncp[j].ncp; k++, npt++) {
            if (senses[jj] == 1) {
              t = ts[0] + ncp[j].knots[k]*dt;
            } else {
              t = ts[1] - ncp[j].knots[k]*dt;
            }
            stat = EG_evaluate(edges[jj], &t, data);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Section %d Edge %d eval = %d (EG_ruled)!\n",
                       i+1, j+1, stat);
              EG_free(xyzs);
              if (curvs != NULL) {
                for (i = 0; i < 2*nstripe; i++)
                  if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
                EG_free(curvs);
              }
              for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
              for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
              EG_free(faces);
              EG_free(surfs);
              EG_freeSeq(nstripe, ncp);
              return stat;
            }
            xyzs[3*npt  ] = data[0];
            xyzs[3*npt+1] = data[1];
            xyzs[3*npt+2] = data[2];
            if (k == 0) {
              if (senses[jj] == 1) {
                t1[3]         =  data[3]*dt;
                t1[4]         =  data[4]*dt;
                t1[5]         =  data[5]*dt;
                xyzs[3*npt  ] =  v1[0];
                xyzs[3*npt+1] =  v1[1];
                xyzs[3*npt+2] =  v1[2];
                EG_checkDirs(edges[jj], t+0.001*dt, data, &t1[3]);
              } else {
                t1[3]         = -data[3]*dt;
                t1[4]         = -data[4]*dt;
                t1[5]         = -data[5]*dt;
                xyzs[3*npt  ] =  v2[0];
                xyzs[3*npt+1] =  v2[1];
                xyzs[3*npt+2] =  v2[2];
                EG_checkDirs(edges[jj], t-0.001*dt, data, &t1[3]);
              }
            } else if (k == ncp[j].ncp-1) {
              if (senses[jj] == 1) {
                tN[3]         = -data[3]*dt;
                tN[4]         = -data[4]*dt;
                tN[5]         = -data[5]*dt;
                xyzs[3*npt  ] =  v2[0];
                xyzs[3*npt+1] =  v2[1];
                xyzs[3*npt+2] =  v2[2];
                EG_checkDirs(edges[jj], t-0.001*dt, data, &tN[3]);
              } else {
                tN[3]         =  data[3]*dt;
                tN[4]         =  data[4]*dt;
                tN[5]         =  data[5]*dt;
                xyzs[3*npt  ] =  v1[0];
                xyzs[3*npt+1] =  v1[1];
                xyzs[3*npt+2] =  v1[2];
                EG_checkDirs(edges[jj], t+0.001*dt, data, &tN[3]);
              }
            }
          }
          break;
        }
        nn++;
      }

finis:
      /* get the BSpline surface */
      if ((planar == 1) || (ncp[j].ncp == 3)) {
        stat = EG_spline2dAppx(context, 1, ncp[j].knots, NULL, NULL, NULL, NULL,
                               NULL, NULL, NULL, NULL, ncp[j].ncp, 2, xyzs,
                               1.e-8, &surfs[npatch]);
      } else {
        stat = EG_spline2dAppx(context, 1, ncp[j].knots, NULL, NULL, t1, tN,
                               NULL, NULL, NULL, NULL, ncp[j].ncp, 2, xyzs,
                               1.e-8, &surfs[npatch]);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d spline2d = %d (EG_ruled)!\n",
                 j+1, i+1, stat);
        EG_free(xyzs);
        if (curvs != NULL) {
          for (i = 0; i < 2*nstripe; i++)
            if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
          EG_free(curvs);
        }
        for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
        for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
        EG_free(faces);
        EG_free(surfs);
        EG_freeSeq(nstripe, ncp);
        return stat;
      }
      if ((cap != 0) && (curvs != NULL)) {
        if ((cap&1) != 0) {
          stat = EG_isoCline(surfs[npatch], 1, 0.0, &curvs[j]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCline- = %d (EG_ruled)!\n",
                     j+1, i+1, stat);
            EG_free(xyzs);
            if (curvs != NULL) {
              for (i = 0; i < 2*nstripe; i++)
                if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
              EG_free(curvs);
            }
            for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
            for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
            EG_free(faces);
            EG_free(surfs);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
        }
        if ((cap&2) != 0) {
          stat = EG_isoCline(surfs[npatch], 1, 1.0, &curvs[j+nstripe]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCline+ = %d (EG_ruled)!\n",
                     j+1, i+1, stat);
            EG_free(xyzs);
            if (curvs != NULL) {
              for (i = 0; i < 2*nstripe; i++)
                if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
              EG_free(curvs);
            }
            for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
            for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
            EG_free(faces);
            EG_free(surfs);
            EG_freeSeq(nstripe, ncp);
            return stat;
          }
        }
      }
      lims[0] = lims[2] = 0.0;
      lims[1] = lims[3] = 1.0;
      stat = EG_makeFace(surfs[npatch], SFORWARD, lims, &faces[npatch]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d makeFace = %d (EG_ruled)!\n",
                 j+1, i+1, stat);
        EG_free(xyzs);
        if (curvs != NULL) {
          for (i = 0; i < 2*nstripe; i++)
            if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
          EG_free(curvs);
        }
        for (i = 0; i <= npatch; i++) EG_deleteObject(faces[i]);
        for (i = 0; i <  npatch; i++) EG_deleteObject(surfs[i]);
        EG_free(faces);
        EG_free(surfs);
        EG_freeSeq(nstripe, ncp);
        return stat;
      }
      npatch++;
    }
    EG_free(xyzs);
  }
  EG_freeSeq(nstripe, ncp);
  
  /* add cap faces if necessary */
  if ((ncap != 0) && (curvs != NULL)) {
    nodes = (ego *) EG_alloc(2*nstripe*sizeof(ego));
    edges = (ego *) EG_alloc(4*nstripe*sizeof(ego));
    if ((nodes == NULL) || (edges == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge %d/%d makeFace = %d (EG_ruled)!\n",
               j+1, i+1, stat);
      if (nodes != NULL) EG_free(nodes);
      if (edges != NULL) EG_free(edges);
      for (i = 0; i < 2*nstripe; i++)
        if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
      EG_free(curvs);
      for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
      for (i = 0; i < npatch; i++) EG_deleteObject(surfs[i]);
      EG_free(faces);
      EG_free(surfs);
      return EGADS_MALLOC;
    }
    for (i = 0; i < 2*nstripe; i++) nodes[i] = NULL;
    for (i = 0; i < 4*nstripe; i++) edges[i] = NULL;

    stat = EG_getTopology(secs[0], &ref, &oclass, &mtype, data, &n,
                          &chldrn, &senses);
    if ((stat == EGADS_SUCCESS) && (oclass == FACE)) {
      t    = 0.0;
      stat = EG_evaluate(curvs[0], &t, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Evaluate- curve 0 = %d (EG_ruled)!\n", stat);
        goto caperr;
      }
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[0]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo- node 0 = %d (EG_ruled)!\n", stat);
        goto caperr;
      }
      for (i = 1; i < nstripe; i++) {
        t    = 1.0;
        stat = EG_evaluate(curvs[i-1], &t, data);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Evaluate- curve %d = %d (EG_ruled)!\n",
                   i, stat);
          goto caperr;
        }
        stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                               &nodes[i]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo- node %d = %d (EG_ruled)!\n",
                   i, stat);
          goto caperr;
        }
      }
      ts[0] = 0.0;
      ts[1] = 1.0;
      for (i = 0; i < nstripe; i++) {
        if (nstripe == 1) {
          n      = 1;
          nds[0] = nodes[0];
        } else {
          n      = 2;
          nds[0] = nodes[i];
          if (i == nstripe-1) {
            nds[1] = nodes[0];
          } else {
            nds[1] = nodes[i+1];
          }
        }
        stat = EG_makeTopology(context, curvs[i], EDGE, 0, ts, n, nds, NULL,
                               &edges[i]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo- edge %d = %d (EG_ruled)!\n",
                   i, stat);
          goto caperr;
        }
      }
      senses = (int *) EG_alloc(nstripe*sizeof(int));
      stat   = EGADS_MALLOC;
      if (senses == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: can not allocate %d senses- (EG_ruled)!\n",
                 nstripe);
        goto caperr;
      }

      for (i = 0; i < nstripe; i++) senses[i] = 1;
      if (ref->mtype == PLANE) {
        stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, nstripe,
                               edges, senses, &bLoop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Beg Loop = %d (EG_ruled)!\n", stat);
          EG_free(senses);
          goto caperr;
        }
        EG_free(senses);
/*@-nullpass@*/
        stat = EG_makeFace(bLoop, SFORWARD, NULL, &faces[npatch]);
/*@+nullpass@*/
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Beg Cap makeFace = %d (EG_ruled)!\n", stat);
          goto caperr;
        }
      } else {
/*      printf(" first section non-planar!\n");  */
        /* make PCurves */
        for (i = 0; i < nstripe; i++) {
          stat = EG_otherCurve(ref, curvs[i], 1.e-7, &edges[nstripe+i]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: %d otherCurve First = %d (EG_ruled)!\n",
                     i, stat);
            EG_free(senses);
            goto caperr;
          }
        }
        stat = EG_makeTopology(context, ref, LOOP, CLOSED, NULL, nstripe,
                               edges, senses, &bLoop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Beg Loop = %d (EG_ruled)!\n", stat);
            EG_free(senses);
          goto caperr;
        }
        stat = EG_makeTopology(context, ref, FACE, SFORWARD, NULL, 1,
                               &bLoop, senses, &faces[npatch]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Beg Cap makeTopo = %d (EG_ruled)!\n", stat);
          EG_free(senses);
          goto caperr;
        }
        EG_free(senses);
      }
      surfs[npatch] = NULL;
      EG_attributeDup(secs[0], faces[npatch]);
      npatch++;
    } else {
      bLoop = NULL;
    }

    stat = EG_getTopology(secs[nsec-1], &ref, &oclass, &mtype, data, &n,
                          &chldrn, &senses);
    if ((stat == EGADS_SUCCESS) && (oclass == FACE)) {
      t    = 0.0;
      stat = EG_evaluate(curvs[nstripe], &t, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Evaluate+ curve 0 = %d (EG_ruled)!\n", stat);
        goto caperr;
      }
      stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                             &nodes[nstripe]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo+ node 0 = %d (EG_ruled)!\n", stat);
        goto caperr;
      }
      for (i = 1; i < nstripe; i++) {
        t    = 1.0;
        stat = EG_evaluate(curvs[nstripe+i-1], &t, data);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Evaluate+ curve %d = %d (EG_ruled)!\n",
                   i, stat);
          goto caperr;
        }
        stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                               &nodes[nstripe+i]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo+ node %d = %d (EG_ruled)!\n",
                   i, stat);
          goto caperr;
        }
      }
      ts[0] = 0.0;
      ts[1] = 1.0;
      for (i = 0; i < nstripe; i++) {
        if (nstripe == 1) {
          n      = 1;
          nds[0] = nodes[nstripe];
        } else {
          n      = 2;
          nds[0] = nodes[nstripe+i];
          if (i == nstripe-1) {
            nds[1] = nodes[nstripe];
          } else {
            nds[1] = nodes[nstripe+i+1];
          }
        }
        stat = EG_makeTopology(context, curvs[nstripe+i], EDGE, 0, ts, n, nds,
                               NULL, &edges[2*nstripe+i]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo+ edge %d = %d (EG_ruled)!\n",
                   i, stat);
          goto caperr;
        }
      }
      senses = (int *) EG_alloc(nstripe*sizeof(int));
      stat   = EGADS_MALLOC;
      if (senses == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: can not allocate %d senses+ (EG_ruled)!\n",
                 nstripe);
        goto caperr;
      }
      for (i = 0; i < nstripe; i++) senses[i] = 1;
      if (ref->mtype == PLANE) {
        stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, nstripe,
                               &edges[2*nstripe], senses, &eLoop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo End Loop = %d (EG_ruled)!\n", stat);
          EG_free(senses);
          goto caperr;
        }
        EG_free(senses);
/*@-nullpass@*/
        stat = EG_makeFace(eLoop, SFORWARD, NULL, &faces[npatch]);
/*@+nullpass@*/
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: End Cap makeFace = %d (EG_ruled)!\n", stat);
          goto caperr;
        }
      } else {
/*      printf(" last section non-planar!\n");  */
        /* make PCurves */
        for (i = 0; i < nstripe; i++) {
          stat = EG_otherCurve(ref, curvs[nstripe+i], 1.e-7,
                               &edges[3*nstripe+i]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: %d otherCurve Last = %d (EG_ruled)!\n",
                     i, stat);
            EG_free(senses);
            goto caperr;
          }
        }
        stat = EG_makeTopology(context, ref, LOOP, CLOSED, NULL, nstripe,
                               &edges[2*nstripe], senses, &eLoop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo End Loop = %d (EG_ruled)!\n", stat);
          EG_free(senses);
          goto caperr;
        }
        stat = EG_makeTopology(context, ref, FACE, SFORWARD, NULL, 1,
                               &eLoop, senses, &faces[npatch]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: End Cap makeTopo = %d (EG_ruled)!\n", stat);
          EG_free(senses);
          goto caperr;
        }
        EG_free(senses);
      }
      surfs[npatch] = NULL;
      EG_attributeDup(secs[nsec-1], faces[npatch]);
      npatch++;
    } else {
      eLoop = NULL;
    }
  } else {
    bLoop = eLoop = NULL;
    nodes = edges = NULL;
  }
  
  body = NULL;
  if (npatch == 1) {
    stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1,
                           &faces[0], NULL, &body);
    if (stat == EGADS_SUCCESS)
      stat = EG_makeTopology(context, NULL, MODEL, 0, NULL, 1,
                             &body, NULL, &model);
  } else {
    /* sew all the faces together */
    stat = EG_sewFaces(npatch, faces, 0.0, 0, &model);
  }
  
  /* clean up all of our temps */
  for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
  for (i = 0; i < npatch; i++)
    if (surfs[i] != NULL) EG_deleteObject(surfs[i]);
  EG_free(faces);
  EG_free(surfs);
  if (curvs != NULL) {
    if (bLoop != NULL) EG_deleteObject(bLoop);
    if (eLoop != NULL) EG_deleteObject(eLoop);
    if (edges != NULL) {
      for (i = 0; i < 4*nstripe; i++)
        if (edges[i] != NULL) EG_deleteObject(edges[i]);
      EG_free(edges);
    }
    if (nodes != NULL) {
      for (i = 0; i < 2*nstripe; i++)
        if (nodes[i] != NULL) EG_deleteObject(nodes[i]);
      EG_free(nodes);
    }
    if (curvs != NULL) {
      for (i = 0; i < 2*nstripe; i++)
        if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
      EG_free(curvs);
    }
  }
  if (npatch == 1) {
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo Body/Model = %d (EG_ruled)!\n", stat);
      if (body != NULL) EG_deleteObject(body);
      return stat;
    }
  } else {
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: sewFaces = %d (EG_ruled)!\n", stat);
      return stat;
    }
  }
  
  /* get the body out of the model */
  stat = EG_getTopology(model, &ref, &oclass, &mtype, data, &n, &chldrn,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Model getTopology = %d (EG_ruled)!\n", stat);
    EG_deleteObject(model);
    return stat;
  }
  if (n != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: #Bodies = %d (EG_ruled)!\n", n);
    EG_deleteObject(model);
    return EGADS_TOPOERR;
  }
  stat = EG_copyObject(chldrn[0], NULL, result);
  EG_deleteObject(model);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Model copyObject = %d (EG_ruled)!\n", stat);
    return stat;
  }

  return EGADS_SUCCESS;
  
  /* cap construction error cleanup */
caperr:
  for (i = 0; i < npatch; i++) EG_deleteObject(faces[i]);
  for (i = 0; i < npatch; i++)
    if (surfs[i] != NULL) EG_deleteObject(surfs[i]);
  EG_free(faces);
  EG_free(surfs);
  if (bLoop != NULL) EG_deleteObject(bLoop);
  if (eLoop != NULL) EG_deleteObject(eLoop);
  if (edges != NULL) {
    for (i = 0; i < 4*nstripe; i++)
      if (edges[i] != NULL) EG_deleteObject(edges[i]);
    EG_free(edges);
  }
  if (nodes != NULL) {
    for (i = 0; i < 2*nstripe; i++)
      if (nodes[i] != NULL) EG_deleteObject(nodes[i]);
    EG_free(nodes);
  }
  if (curvs != NULL) {
    for (i = 0; i < 2*nstripe; i++)
      if (curvs[i] != NULL) EG_deleteObject(curvs[i]);
    EG_free(curvs);
  }
  
  return stat;
}
