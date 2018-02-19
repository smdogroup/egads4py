/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Quad Tessellation Functions
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"


#define MAXSIDE 2049

#define AREA2D(a,b,c)   ((a[0]-c[0])*(b[1]-c[1]) -  (a[1]-c[1])*(b[0]-c[0]))
#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define MAX(a,b)	(((a) > (b)) ? (a) : (b))
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))


typedef struct {
  int nodes[4];         /* quad indices into Node list */
} Quad;


typedef struct {
  double uv[2];         /* (u,v) for node */
  double duv[2];        /* delta for coordinate update */
  double area;          /* accumulated area; -1 is boundary node */
  double xyz[3];	/* xyz for the node */
} Node;


typedef struct {
  int    fillT;
  int    nvert;
  int    nquad;
  int    ntris;
  Node   *verts;
  Quad   *quads;
  int    *tris;
  int    sizes[8];
  int    last[MAXSIDE];
  int    outLevel;
  double flip;
  int    npatch;
  int    *vpatch;
  int    patch[17][2];
} qFill;


  extern int EG_evaluate( const egObject *geom, const double *param,
                          double *results );


/* Compute arclength basis functions for TFI use */

static void 
EG_arcBasis(qFill *q, int nx, int ny, int *sideptr[], double *abasis[2])
{
  int    i, j, k, i0, im, j0, jm, nny = ny+1;
  double xi, et;
  double anorm;

  for (j = k = 0; j <= ny; j += ny, k = 2) {		/* j const boundaries */
    abasis[0][j] = 0.0;					/* i == 0 */
    for (i = 1; i <= nx; i++) {
      i0 = sideptr[k][i  ];
      im = sideptr[k][i-1];
      abasis[0][nny*i+j] = abasis[0][nny*(i-1)+j] + 
                                sqrt((q->verts[i0].uv[0] - q->verts[im].uv[0]) *
                                     (q->verts[i0].uv[0] - q->verts[im].uv[0]) +
                                     (q->verts[i0].uv[1] - q->verts[im].uv[1]) *
                                     (q->verts[i0].uv[1] - q->verts[im].uv[1]));
    }
    if (abasis[0][nny*nx+j] > 1.0e-6) {			/* degenerate */
      anorm = 1.0/abasis[0][nny*nx+j];
      for (i = 0; i <= nx; i++) abasis[0][nny*i+j] *= anorm;
    } else {
      anorm = 1.0/(double)(nx);
      for (i = 0; i <= nx; i++) abasis[0][nny*i+j] = (double) i*anorm;
    }
  }
  for (i = 0; i <= nx; i++) {				/* boundaries */
    abasis[1][nny*i   ] = 0.0;
    abasis[1][nny*i+ny] = 1.0;
  }

  for (i = 0, k = 1; i <= nx; i+=nx, k = 3) {		/* i const boundaries */
    abasis[1][nny*i] = 0.0;
    for (j = 1; j <= ny; j++) {
      j0 = sideptr[k][j  ];
      jm = sideptr[k][j-1];
      abasis[1][nny*i+j] = abasis[1][nny*i+(j-1)] + 
                                sqrt((q->verts[j0].uv[0] - q->verts[jm].uv[0]) *
                                     (q->verts[j0].uv[0] - q->verts[jm].uv[0]) +
                                     (q->verts[j0].uv[1] - q->verts[jm].uv[1]) *
                                     (q->verts[j0].uv[1] - q->verts[jm].uv[1]));
    }
    if (abasis[1][nny*i+ny] > 1.0e-6) {			/* degenerate */
      anorm = 1.0/abasis[1][nny*i+ny];
      for (j = 0; j <= ny; j++) abasis[1][nny*i+j] *= anorm;
    } else {
      anorm = 1.0/(double)(ny);
      for (j = 0; j <= ny; j++) abasis[1][nny*i+j] = (double) j*anorm;
    }
  }
  for (j = 0; j <= ny; j++) {
    abasis[0][       j] = 0.0;
    abasis[0][nny*nx+j] = 1.0;
  }

  for (j = 1; j < ny; j++) {
    for (i = 1; i < nx; i++) {
      anorm = 1.0 - (abasis[0][nny*i +ny]-abasis[0][nny*i  ]) *
                    (abasis[1][nny*nx+ j]-abasis[1][      j]);

      xi = ( abasis[0][nny*i  ] -
             abasis[1][      j]*(abasis[0][nny*i +ny]-abasis[0][nny*i    ]) ) /
           anorm;

      et = ( abasis[1][      j] -
             abasis[0][nny*i  ]*(abasis[1][nny*nx+ j]-abasis[1][        j]) ) /
           anorm;

      abasis[0][nny*i+j] = xi;
      abasis[1][nny*i+j] = et;
    }
  }

}


/* remap into the actual UV space */

static void
EG_getside(int iuv, double t, int len, int *side, double *uvx,
           double *uv, double *uvi)
{
  int    i0, i1, j;
  double dis;

  for (j = 1; j < len; j++) {
    i0 = side[j-1];
    i1 = side[j  ];
    if (((t >= uvx[2*i0+iuv]) && (t <= uvx[2*i1+iuv])) ||
        ((t >= uvx[2*i1+iuv]) && (t <= uvx[2*i0+iuv]))) {
      dis = (t - uvx[2*i0+iuv])/(uvx[2*i1+iuv]-uvx[2*i0+iuv]);
      uvi[0] = uv[2*i0  ] - dis*(uv[2*i0  ] - uv[2*i1  ]);
      uvi[1] = uv[2*i0+1] - dis*(uv[2*i0+1] - uv[2*i1+1]);
      return;
    }
  }
}


static int
EG_dQuadTFI(int *elen, double *uv, int npts, double *uvx)
{
  int    i, j, k, len, ll, lr, ur, ul;
  int    cipt[4], *sideptr[4];
  double et, xi, uvi[2][2], uvj[2][2], smap[4];

  cipt[ 0] = 0;
  len      = elen[0];
  cipt[ 1] = len;
  len     += elen[1];
  cipt[ 2] = len;
  len     += elen[2];
  cipt[ 3] = len;

  /* set the exterior block sides */

  for (i = 0; i < 4; i++) sideptr[i] = NULL;
  for (j = i = 0; i < 4; i++) {
    len = elen[i] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      return -1;
    }
    if (i >= 2) {
      for (k = len-1; k > 0; k--, j++) {
        sideptr[i][k] = j;
      }
      sideptr[i][0] = j;
      if (i == 3) sideptr[i][0] = 0;
    } else {
      for (k = 0; k < len-1; k++, j++) {
        sideptr[i][k] = j;
      }
      sideptr[i][len-1] = j;
    }
  }

  /* create the quads and get coordinates via TFI */

  len = elen[0] + elen[1] + elen[2] + elen[3];
  ll  = 2*cipt[0];
  lr  = 2*cipt[1];
  ur  = 2*cipt[2];
  ul  = 2*cipt[3];
  for (k = len; k < npts; k++) {
    EG_getside(0, uvx[2*k  ], elen[3]+1, sideptr[3], uvx, uv, uvi[0]);
    EG_getside(0, uvx[2*k  ], elen[1]+1, sideptr[1], uvx, uv, uvi[1]);
    EG_getside(1, uvx[2*k+1], elen[0]+1, sideptr[0], uvx, uv, uvj[0]);
    EG_getside(1, uvx[2*k+1], elen[2]+1, sideptr[2], uvx, uv, uvj[1]);
    smap[3] = sqrt((uvi[0][0]-uv[ll  ])*(uvi[0][0]-uv[ll  ]) +
                   (uvi[0][1]-uv[ll+1])*(uvi[0][1]-uv[ll+1])) /
              sqrt((uv[ul  ] -uv[ll  ])*(uv[ul  ] -uv[ll  ]) +
                   (uv[ul+1] -uv[ll+1])*(uv[ul+1] -uv[ll+1]));
    smap[1] = sqrt((uvi[1][0]-uv[lr  ])*(uvi[1][0]-uv[lr  ]) +
                   (uvi[1][1]-uv[lr+1])*(uvi[1][1]-uv[lr+1])) /
              sqrt((uv[ur  ] -uv[lr  ])*(uv[ur  ] -uv[lr  ]) +
                   (uv[ur+1] -uv[lr+1])*(uv[ur+1] -uv[lr+1]));
    smap[0] = sqrt((uvj[0][0]-uv[ll  ])*(uvj[0][0]-uv[ll  ]) +
                   (uvj[0][1]-uv[ll+1])*(uvj[0][1]-uv[ll+1])) /
              sqrt((uv[lr  ] -uv[ll  ])*(uv[lr  ] -uv[ll  ]) +
                   (uv[lr+1] -uv[ll+1])*(uv[lr+1] -uv[ll+1]));
    smap[2] = sqrt((uvj[1][0]-uv[ul  ])*(uvj[1][0]-uv[ul  ]) +
                   (uvj[1][1]-uv[ul+1])*(uvj[1][1]-uv[ul+1])) /
              sqrt((uv[ur  ] -uv[ul  ])*(uv[ur  ] -uv[ul  ]) +
                   (uv[ur+1] -uv[ul+1])*(uv[ur+1] -uv[ul+1]));
    et = smap[3]*(1.0-uvx[2*k+1]) + smap[1]*uvx[2*k+1];
    xi = smap[0]*(1.0-uvx[2*k  ]) + smap[2]*uvx[2*k  ];

    uvx[2*k  ] = (1.0-xi)            * uvi[0][0] +
                 (    xi)            * uvi[1][0] +
                            (1.0-et) * uvj[0][0] +
                            (    et) * uvj[1][0] -
                 (1.0-xi) * (1.0-et) * uv[ll  ] -
                 (1.0-xi) * (    et) * uv[ul  ] -
                 (    xi) * (1.0-et) * uv[lr  ] -
                 (    xi) * (    et) * uv[ur  ];
    uvx[2*k+1] = (1.0-xi)            * uvi[0][1] +
                 (    xi)            * uvi[1][1] +
                            (1.0-et) * uvj[0][1] +
                            (    et) * uvj[1][1] -
                 (1.0-xi) * (1.0-et) * uv[ll+1] -
                 (1.0-xi) * (    et) * uv[ul+1] -
                 (    xi) * (1.0-et) * uv[lr+1] -
                 (    xi) * (    et) * uv[ur+1];
  }
  for (k = 0; k < 2*len; k++) uvx[k] = uv[k];

  /* free up our integrated sides */
  for (i = 0; i < 4; i++) EG_free(sideptr[i]);

  return 0;
}


/* get the vertex count for the suite of blocks */

static int
EG_getVertCnt(qFill *q, int len, int blocks[][6])
{
  int k, cnt;
  
  q->npatch = len;
  for (cnt = k = 0; k < len; k++) {
    q->patch[k][0] = q->sizes[blocks[k][0]] + 1;
    q->patch[k][1] = q->sizes[blocks[k][1]] + 1;
    cnt += q->patch[k][0]*q->patch[k][1];
  }

  return cnt;
}


/* sets the individual quads by looping through the blocks */

static void 
EG_setQuads(qFill *q, int len, int blocks[][6], int *sideptr[])
{
  int    i, j, k, i0, i1, i2, i3, ilast, nx, ny;
  int    ll, lr, ur, ul, j0, jm, ii, im, iv, sav;
  double et, xi;
  Node   *verts;
  Quad   *quads;

  verts = q->verts;
  quads = q->quads;
  q->nquad = iv = 0;
  for (k = 0; k < len; k++) {
    nx = q->sizes[blocks[k][0]];
    ny = q->sizes[blocks[k][1]];
    i0 = blocks[k][2];
    i1 = blocks[k][3];
    i2 = blocks[k][4];
    i3 = blocks[k][5];
    ll = sideptr[i2][0];
    lr = sideptr[i2][nx];
    ur = sideptr[i3][nx];
    ul = sideptr[i3][0];
    for (i = 0; i < nx+1; i++) q->last[i] = sideptr[i2][i];
    for (j = 0; j < ny; j++) {
      ii    = sideptr[i0][j+1];
      if (i1 > 0) {
        im  = sideptr[i1][j+1];
      } else {
        im = sideptr[-i1][ny-j-1];
      }
      et    = ((double) (j+1)) / ((double) ny);
      ilast = sideptr[i0][j+1];
      sav   = q->nquad;
      for (i = 0; i < nx; i++) {
        j0 = sideptr[i2][i+1];
        jm = sideptr[i3][i+1];
        xi = ((double) (i+1)) / ((double) nx);
        quads[q->nquad].nodes[0] = q->last[i  ];
        quads[q->nquad].nodes[1] = q->last[i+1];
        if (j == ny-1) {
          quads[q->nquad].nodes[2] = sideptr[i3][i+1];
          quads[q->nquad].nodes[3] = sideptr[i3][i  ];
        } else {
          if (i == nx-1) {
            if (i1 > 0) {
              quads[q->nquad].nodes[2] = sideptr[i1][j+1];
              q->last[i]  = ilast;
              q->last[nx] = sideptr[i1][j+1];
            } else {
              quads[q->nquad].nodes[2] = sideptr[-i1][ny-j-1];
              q->last[i]  = ilast;
              q->last[nx] = sideptr[-i1][ny-j-1];
            }
          } else {
            quads[q->nquad].nodes[2] = q->nvert;
            verts[q->nvert].uv[0]    = (1.0-xi)            * verts[ii].uv[0] +
                                       (    xi)            * verts[im].uv[0] +
                                                  (1.0-et) * verts[j0].uv[0] +
                                                  (    et) * verts[jm].uv[0] -
                                       (1.0-xi) * (1.0-et) * verts[ll].uv[0] -
                                       (1.0-xi) * (    et) * verts[ul].uv[0] -
                                       (    xi) * (1.0-et) * verts[lr].uv[0] -
                                       (    xi) * (    et) * verts[ur].uv[0];
            verts[q->nvert].uv[1]    = (1.0-xi)            * verts[ii].uv[1] +
                                       (    xi)            * verts[im].uv[1] +
                                                  (1.0-et) * verts[j0].uv[1] +
                                                  (    et) * verts[jm].uv[1] -
                                       (1.0-xi) * (1.0-et) * verts[ll].uv[1] -
                                       (1.0-xi) * (    et) * verts[ul].uv[1] -
                                       (    xi) * (1.0-et) * verts[lr].uv[1] -
                                       (    xi) * (    et) * verts[ur].uv[1];
            verts[q->nvert].area     = 0.0;
            q->nvert++;
          }
          quads[q->nquad].nodes[3] = ilast;
          q->last[i] = ilast;
          ilast      = q->nvert-1;
        }
        q->nquad++;
      }
      if (q->fillT == 0) {
        if (j == 0) {
          q->vpatch[iv] = quads[sav].nodes[0];
          iv++;
          for (i = 0; i < nx; i++, iv++)
            q->vpatch[iv] = quads[sav+i].nodes[1];
        }
        q->vpatch[iv] = quads[sav].nodes[3];
        iv++;
        for (i = 0; i < nx; i++, sav++, iv++)
          q->vpatch[iv] = quads[sav].nodes[2];
      }
    }
  }
}


/* perform the laplacian smoothing on the grid vertices */

static void
EG_smoothQuads(const egObject *face, qFill *q, int len, int npass)
{
  int           i, j, i0, i1, i2, i3, status, pass;
  double        qarea, sum, big, delta1, sums[2], x1[3], x2[3], xn[3];
  double        tAreaUV, tAreaXYZ, holdArea, results[18];
  Node          *verts;
  static double wXYZ = 0.75;

  verts = q->verts;

  /* outer iteration -- pass 1 (uv only) */

  for (i = 0; i < len; i++) {

    /* initialize deltas */
    for (j = 0; j < q->nvert; j++) {
      verts[j].duv[0] = verts[j].duv[1] = 0.0;
      if (verts[j].area > 0.0) verts[j].area = 0.0;
    }
    /* calculate and distribute change */
    for (j = 0; j < q->nquad; j++) {
      i0    = q->quads[j].nodes[0];
      i1    = q->quads[j].nodes[1];
      i2    = q->quads[j].nodes[2];
      i3    = q->quads[j].nodes[3];
      qarea = q->flip*(AREA2D(verts[i0].uv, verts[i1].uv, verts[i2].uv) +
                       AREA2D(verts[i0].uv, verts[i2].uv, verts[i3].uv));
      if (qarea <= 0.0) {
#ifdef DEBUG
        printf(" Quad %d: NonPositive Area = %le  %d %d %d %d\n",
               j, qarea, i0, i1, i2, i3);
        printf("          %lf %lf  %lf %lf  %lf %lf  %lf %lf\n",
               verts[i0].uv[0],verts[i0].uv[1], verts[i1].uv[0],verts[i1].uv[1],
               verts[i2].uv[0],verts[i2].uv[1], verts[i3].uv[0],verts[i3].uv[1]);
#endif
        qarea = -qarea;
      }

      sum = qarea*(verts[i0].uv[0] + verts[i1].uv[0] + verts[i2].uv[0] +
                   verts[i3].uv[0])/4.0;
      verts[i0].duv[0] += sum;
      verts[i1].duv[0] += sum;
      verts[i2].duv[0] += sum;
      verts[i3].duv[0] += sum;
      sum = qarea*(verts[i0].uv[1] + verts[i1].uv[1] + verts[i2].uv[1] +
                   verts[i3].uv[1])/4.0;
      verts[i0].duv[1] += sum;
      verts[i1].duv[1] += sum;
      verts[i2].duv[1] += sum;
      verts[i3].duv[1] += sum;

      if (verts[i0].area >= 0.0) verts[i0].area += qarea;
      if (verts[i1].area >= 0.0) verts[i1].area += qarea;
      if (verts[i2].area >= 0.0) verts[i2].area += qarea;
      if (verts[i3].area >= 0.0) verts[i3].area += qarea;
    }
    /* update distributions */
    big = 0.0;
    for (j = 0; j < q->nvert; j++) {
      if (verts[j].area <= 0.0) continue;
      sums[0] = verts[j].duv[0]/verts[j].area;
      sums[1] = verts[j].duv[1]/verts[j].area;
      sum     = fabs(sums[0] - verts[j].uv[0]);
      if (big < sum) big = sum;
      sum     = fabs(sums[1] - verts[j].uv[1]);
      if (big < sum) big = sum;
      verts[j].uv[0] = sums[0];
      verts[j].uv[1] = sums[1];
    }
    if (i == 0) {
      delta1 = big;
      if (delta1 == 0.0) break;
    } else {
      if (big/delta1 < 1.e-3) break;
    }
  }

  /* pseudo non-linear loop */

  for (pass = 0; pass < npass; pass++) {

    /* get xyz */
    for (j = 0; j < q->nvert; j++) {
      status = EG_evaluate(face, verts[j].uv, results);
      if (status != EGADS_SUCCESS) {
        if (q->outLevel > 0)
          printf(" EGADS Info: EG_evaluate = %d (EG_smoothQuad)!\n", 
                 status);
        return;
      }
      verts[j].xyz[0] = results[0];
      verts[j].xyz[1] = results[1];
      verts[j].xyz[2] = results[2];
    }

    tAreaUV = tAreaXYZ = 0.0;
    for (j  = 0; j < q->nquad; j++) {
      i0        = q->quads[j].nodes[0];
      i1        = q->quads[j].nodes[1];
      i2        = q->quads[j].nodes[2];
      i3        = q->quads[j].nodes[3];
      holdArea  = q->flip*(AREA2D(verts[i0].uv, verts[i1].uv, verts[i2].uv) +
                           AREA2D(verts[i0].uv, verts[i2].uv, verts[i3].uv));
      if (holdArea < 0.0) holdArea = -holdArea;
      tAreaUV  += holdArea;
      x1[0]     = verts[i1].xyz[0] - verts[i0].xyz[0];
      x2[0]     = verts[i2].xyz[0] - verts[i0].xyz[0];
      x1[1]     = verts[i1].xyz[1] - verts[i0].xyz[1];
      x2[1]     = verts[i2].xyz[1] - verts[i0].xyz[1];
      x1[2]     = verts[i1].xyz[2] - verts[i0].xyz[2];
      x2[2]     = verts[i2].xyz[2] - verts[i0].xyz[2];
      CROSS(xn, x1, x2);
      holdArea  = DOT(xn, xn);
      if (holdArea <  0.0) holdArea = -holdArea;
      tAreaXYZ += holdArea;
      x1[0]     = verts[i3].xyz[0] - verts[i0].xyz[0];
      x1[1]     = verts[i3].xyz[1] - verts[i0].xyz[1];
      x1[2]     = verts[i3].xyz[2] - verts[i0].xyz[2];
      CROSS(xn, x2, x1);
      holdArea  = DOT(xn, xn);
      if (holdArea <  0.0) holdArea = -holdArea;
      tAreaXYZ += holdArea;
    }
#ifdef DEBUG
    printf(" ** %d   Areas = %le  %le **\n", pass, tAreaUV, tAreaXYZ);
#endif

    /* outer iteration -- pass 2 (mix) */
    for (i = 0; i < len; i++) {

      /* initialize deltas */
      for (j = 0; j < q->nvert; j++) {
        verts[j].duv[0] = verts[j].duv[1] = 0.0;
        if (verts[j].area > 0.0) verts[j].area = 0.0;
      }
      /* calculate and distribute change */
      for (j = 0; j < q->nquad; j++) {
        i0    = q->quads[j].nodes[0];
        i1    = q->quads[j].nodes[1];
        i2    = q->quads[j].nodes[2];
        i3    = q->quads[j].nodes[3];
        qarea = q->flip*(AREA2D(verts[i0].uv, verts[i1].uv, verts[i2].uv) +
                         AREA2D(verts[i0].uv, verts[i2].uv, verts[i3].uv));
        if (qarea <= 0.0) {
          qarea = -qarea;
#ifdef DEBUG
          printf(" Quad %d: Neg Area = %le\n", j, qarea);
#endif
        }
        qarea *= (1.0-wXYZ)/tAreaUV;

        x1[0]     = verts[i1].xyz[0] - verts[i0].xyz[0];
        x2[0]     = verts[i2].xyz[0] - verts[i0].xyz[0];
        x1[1]     = verts[i1].xyz[1] - verts[i0].xyz[1];
        x2[1]     = verts[i2].xyz[1] - verts[i0].xyz[1];
        x1[2]     = verts[i1].xyz[2] - verts[i0].xyz[2];
        x2[2]     = verts[i2].xyz[2] - verts[i0].xyz[2];
        CROSS(xn, x1, x2);
        holdArea  = DOT(xn, xn);
        if (holdArea <  0.0) holdArea = -holdArea;
        qarea    += holdArea*wXYZ/tAreaXYZ;
        x1[0]     = verts[i3].xyz[0] - verts[i0].xyz[0];
        x1[1]     = verts[i3].xyz[1] - verts[i0].xyz[1];
        x1[2]     = verts[i3].xyz[2] - verts[i0].xyz[2];
        CROSS(xn, x2, x1);
        holdArea  = DOT(xn, xn);
        if (holdArea <  0.0) holdArea = -holdArea;
        qarea    += holdArea*wXYZ/tAreaXYZ;

        sum = qarea*(verts[i0].uv[0] + verts[i1].uv[0] + verts[i2].uv[0] +
                     verts[i3].uv[0])/4.0;
        verts[i0].duv[0] += sum;
        verts[i1].duv[0] += sum;
        verts[i2].duv[0] += sum;
        verts[i3].duv[0] += sum;
        sum = qarea*(verts[i0].uv[1] + verts[i1].uv[1] + verts[i2].uv[1] +
                     verts[i3].uv[1])/4.0;
        verts[i0].duv[1] += sum;
        verts[i1].duv[1] += sum;
        verts[i2].duv[1] += sum;
        verts[i3].duv[1] += sum;

        if (verts[i0].area >= 0.0) verts[i0].area += qarea;
        if (verts[i1].area >= 0.0) verts[i1].area += qarea;
        if (verts[i2].area >= 0.0) verts[i2].area += qarea;
        if (verts[i3].area >= 0.0) verts[i3].area += qarea;
      }
      /* update distributions */
      big = 0.0;
      for (j = 0; j < q->nvert; j++) {
        if (verts[j].area <= 0.0) continue;
        sums[0] = verts[j].duv[0]/verts[j].area;
        sums[1] = verts[j].duv[1]/verts[j].area;
        sum     = fabs(sums[0] - verts[j].uv[0]);
        if (big < sum) big = sum;
        sum     = fabs(sums[1] - verts[j].uv[1]);
        if (big < sum) big = sum;
        verts[j].uv[0] = sums[0];
        verts[j].uv[1] = sums[1];
      }
      if (i == 0) {
        delta1 = big;
        if (delta1 == 0.0) break;
      } else {
        if (big/delta1 < 1.e-3) break;
      }
    }

  }

}


/* triangle templating case */

static int
EG_triTemplate(/*@unused@*/ long tID, const egObject *face, qFill *q, int nsp,
               int *indices, int *elens, double *uv, int *npts, double **uvs)
{
  int    i, j, k, len, m, n, p;
  int    i0, i1, i2, i3, cipt[7], *sideptr[9];
  double d02, d13, cpts[7][2], *uvb;
  static int sides[9][3] = { 0,  0,  1,   2,  1,  2,   1,  2,  3,
                             0,  3,  4,   2,  4,  5,   1,  5,  0,
                             1,  1,  6,   0,  6,  5,   2,  3,  6 };

                      /*      size    sides         */
  static int blocks[3][6] = { 1, 0,   0, 7,  6, 5,
                              1, 2,   1, 8,  2, 6,
                              0, 2,   8, 4,  3, 7 };


  /* get and check sizes */

  m = (elens[0] + elens[1] - elens[2])/2;
  n =  elens[1] - m;
  p =  elens[0] - m;
  if (elens[2]-n != p) {
    m++;
    n = elens[1] - m;
    p = elens[0] - m;
  }
#ifdef DEBUG
  printf("%lX EGADS Info: triTemplate Case - %d %d %d\n",
         tID, elens[0], elens[1], elens[2]);
  printf("%lX             N = %d, M = %d, P = %d (%d)\n",
         tID, n, m, p, elens[2]-n);
#endif
  if ((m <= 1) || (n <= 1) || (p <= 1) || (elens[2]-n != p)) return -2;

  q->sizes[0] = m;
  q->sizes[1] = n;
  q->sizes[2] = p;
  for (i = 0; i < 3; i++)
    if (q->sizes[i] > MAXSIDE-1) return -3;

  /* set the 7 critical points -- 6 exterior */

  cpts[0][0] = uv[0];
  cpts[0][1] = uv[1];
  cipt[0]    = indices[0];
  len  = q->sizes[0];
  cpts[1][0] = uv[2*len  ];
  cpts[1][1] = uv[2*len+1];
  cipt[1]    = indices[len];
  len += q->sizes[2];
  cpts[2][0] = uv[2*len  ];
  cpts[2][1] = uv[2*len+1];
  cipt[2]    = indices[len];
  len += q->sizes[1];
  cpts[3][0] = uv[2*len  ];
  cpts[3][1] = uv[2*len+1];
  cipt[3]    = indices[len];
  len += q->sizes[0];
  cpts[4][0] = uv[2*len  ];
  cpts[4][1] = uv[2*len+1];
  cipt[4]    = indices[len];
  len += q->sizes[2];
  cpts[5][0] = uv[2*len  ];
  cpts[5][1] = uv[2*len+1];
  cipt[5]    = indices[len];

  cpts[6][0] = (cpts[1][0] + cpts[3][0] + cpts[5][0])/3.0;
  cpts[6][1] = (cpts[1][1] + cpts[3][1] + cpts[5][1])/3.0;
  
  if (q->fillT == 0) {
    len       = EG_getVertCnt(q, 3, blocks);
    q->vpatch = (int *) EG_alloc(len*sizeof(int));
    if (q->vpatch == NULL) return -1;
  }

  /* allocate our temporary storage */

  len      = MAX(elens[1], elens[2]);
  len      = MAX(elens[0], len);
  q->quads = (Quad *) EG_alloc(len*len*sizeof(Quad));
  if (q->quads == NULL) {
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }
  q->verts = (Node *) EG_alloc((len+1)*(len+1)*sizeof(Node));
  if (q->verts == NULL) {
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  /* initialize the vertices */

  q->nvert = elens[0] + elens[1] + elens[2];
  for (i = 0; i < q->nvert; i++) {
    j = indices[i];
    q->verts[j].uv[0] = uv[2*i  ];
    q->verts[j].uv[1] = uv[2*i+1];
    q->verts[j].area  = -1.0;
  }
  q->verts[q->nvert].uv[0] = cpts[6][0];
  q->verts[q->nvert].uv[1] = cpts[6][1];
  q->verts[q->nvert].area  = 0.0;
  cipt[6]                  = q->nvert;
  q->nvert++;

  /* set the exterior block sides */

  for (i = 0; i < 9; i++) sideptr[i] = NULL;
  for (j = i = 0; i < 6; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    if ((i <= 1) || (i == 5)) {
      for (k = len-1; k > 0; k--, j++) sideptr[i][k] = indices[j];
      sideptr[i][0] = indices[j];
    } else {
      for (k = 0; k < len-1; k++, j++) sideptr[i][k] = indices[j];
      sideptr[i][len-1] = indices[j];
    }
  }

  /* do the interior sides */

  for (i = 6; i < 9; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    i0 = sides[i][1];
    i1 = sides[i][2];
    sideptr[i][0] = cipt[i0];
    for (j = 1; j < len-1; j++) {
      q->verts[q->nvert].uv[0] = cpts[i0][0] +
                                   j*(cpts[i1][0]-cpts[i0][0])/(len-1);
      q->verts[q->nvert].uv[1] = cpts[i0][1] +
                                   j*(cpts[i1][1]-cpts[i0][1])/(len-1);
      q->verts[q->nvert].area  = 0.0;
      sideptr[i][j]            = q->nvert;
      q->nvert++;
    }
    sideptr[i][len-1] = cipt[i1];
  }
/*
  for (i = 0; i < 9; i++) {
    len = q->sizes[sides[i][0]] + 1;
    printf(" %d: ", i);
    for (j = 1; j < len-1; j++) printf(" %d", sideptr[i][j]);
    printf("\n");
  }
*/

  /* start filling the quads by specifying the 3 blocks */

  EG_setQuads(q, 3, blocks, sideptr);

  /* free up our integrated sides */
  for (i = 0; i < 9; i++) EG_free(sideptr[i]);

  /* get the actual storage that we return the data with */

  uvb = (double *) EG_alloc(2*q->nvert*sizeof(double));
  if (uvb == NULL) {
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;    
  }

  /* calculate the actual coordinates */

  len = MAX(elens[1], elens[2]);
  len = MAX(elens[0], len);
  EG_smoothQuads(face, q, len*len, nsp);

  /* fill the memory to be returned */

  for (j = 0; j < q->nvert; j++) {
    uvb[2*j  ] = q->verts[j].uv[0];
    uvb[2*j+1] = q->verts[j].uv[1];
  }
  
  if (q->fillT != 0) {
    q->ntris = 2*q->nquad;
    q->tris  = (int *) EG_alloc(3*q->ntris*sizeof(int));
    if (q->tris == NULL) {
      EG_free(uvb);
      EG_free(q->verts);
      EG_free(q->quads);
      return -1;
    }
    for (j = 0; j < q->nquad; j++) {
      i0  = q->quads[j].nodes[0];
      i1  = q->quads[j].nodes[1];
      i2  = q->quads[j].nodes[2];
      i3  = q->quads[j].nodes[3];
      d02 = ((uvb[2*i2  ]-uvb[2*i0  ])*(uvb[2*i2  ]-uvb[2*i0  ]) +
             (uvb[2*i2+1]-uvb[2*i0+1])*(uvb[2*i2+1]-uvb[2*i0+1]));
      d13 = ((uvb[2*i3  ]-uvb[2*i1  ])*(uvb[2*i3  ]-uvb[2*i1  ]) +
             (uvb[2*i3+1]-uvb[2*i1+1])*(uvb[2*i3+1]-uvb[2*i2+1]));
      if (d02 > d13) {
        q->tris[6*j  ] = i0;
        q->tris[6*j+1] = i1;
        q->tris[6*j+2] = i3;
        q->tris[6*j+3] = i1;
        q->tris[6*j+4] = i2;
        q->tris[6*j+5] = i3;
      } else {
        q->tris[6*j  ] = i0;
        q->tris[6*j+1] = i1;
        q->tris[6*j+2] = i2;
        q->tris[6*j+3] = i0;
        q->tris[6*j+4] = i2;
        q->tris[6*j+5] = i3;
      }
    }
  }

  /* cleanup and exit */

  *npts = q->nvert;
  *uvs  = uvb;
  EG_free(q->verts);

  return 0;
}


/* general blocking case */

static int
EG_quadFillG(const egObject *face, qFill *q, int nsp, int *indices, int *elens,
             double *uv, int *npts, double **uvs)
{
  int    i, j, k, len, N, M, P, Q, extra = -1;
  int    i0, i1, cipt[26], *sideptr[42];
  double sums[2], cpts[26][2], *uvb;
  static int sides[42][3] = {  0,  0,  4,   1,  4, 13,   2, 13, 19, 
                               3, 19, 20,   6, 20, 21,   4, 21, 22,
                               5, 22, 24,   2, 23, 24,   7, 18, 23,
                               1, 25, 18,   6, 11, 25,   7,  7, 11,
                               0,  3,  7,   5,  2,  3,   4,  1,  2,
                               3,  0,  1,   3,  4,  5,   0,  1,  5,
                               3, 13, 14,   1,  5, 14,   2, 14, 20,
                               4,  5,  6,   0,  2,  6,   6,  5,  8,
                               4,  8,  9,   6,  6,  9,   6, 14, 15,
                               1,  8, 15,   2, 15, 21,   4, 15, 16,
                               1,  9, 16,   2, 16, 22,   5, 16, 23,
                               7, 17, 16,   5, 17, 18,   7,  9, 12,
                               1, 12, 17,   7,  6, 10,   6, 10, 12,
                               5,  6,  7,   5, 10, 11,   5, 12, 25 };

                      /*       size    sides           */
  static int blocks[17][6] = { 0, 3,   15, 16,  0, 17,
                               1, 3,   16, 18,  1, 19,
                               2, 3,   18,  3,  2, 20,
                               0, 4,   14, 21, 17, 22,
                               6, 4,   21, 24, 23, 25,
                               1, 6,   23, 26, 19, 27,
                               2, 6,   26,  4, 20, 28,
                               1, 4,   24, 29, 27, 30,
                               2, 4,   29,  5, 28, 31,
                               2, 5,   32,  6, 31,  7,
                               7, 5,   34, 32, 33,  8,
                               1, 7,   35,-33, 30, 36,
                               6, 7,   37, 35, 25, 38,
                               0, 5,   13, 39, 22, 12,
                               7, 5,   39, 40, 37, 11,
                               6, 5,   40, 41, 38, 10,
                               1, 5,   41, 34, 36,  9 };

                      /*      center   stencil nodes     */     
  static int interior[10][6] = { 5,    1,  4, 14,  8,  6,
                                 6,    5,  7,  9, 10,  2,
                                 8,    5,  9, 15, -1, -1,
                                 9,    8, 12, 16,  6, -1,
                                10,   12,  6, 11, -1, -1,
                                12,    9, 25, 17, 10, -1,
                                14,   13, 15, 20,  5, -1,
                                15,   14, 16, 21,  8, -1,
                                16,   15, 17, 22,  9, 25,
                                17,   16, 18, 12, -1, -1 };
  
  /* get and check sizes */

  N =  elens[0];
  M =  elens[3];
  P =  elens[1] - M;
  Q = (elens[2] - N - P)/2;
  if (Q*2 != elens[2]-N-P)
    if (q->fillT == 0) {
      if (q->outLevel > 0) {
        printf(" EGADS Info: General case off by 1 - %d %d  %d %d\n",
               elens[0], elens[2], elens[1], elens[3]);
        printf("             N = %d, M = %d, P = %d, Q = %d\n",
               N, M, P, Q);
      }
      return -2;
    } else {
      extra = 0;
#ifdef DEBUG
      printf(" EGADS Info: General case - %d %d  %d %d\n",
             elens[0], elens[2], elens[1], elens[3]);
      printf("             N = %d, M = %d, P = %d, Q = %d\n",
             N, M, P, Q);
#endif
    }

  q->sizes[0] = q->sizes[1] = q->sizes[2] = N/3;
  if (q->sizes[0]+q->sizes[1]+q->sizes[2] != N) q->sizes[0]++;
  if (q->sizes[0]+q->sizes[1]+q->sizes[2] != N) q->sizes[2]++;
  q->sizes[3] = q->sizes[4] = q->sizes[5] = M/3;
  if (q->sizes[3]+q->sizes[4]+q->sizes[5] != M) q->sizes[3]++;
  if (q->sizes[3]+q->sizes[4]+q->sizes[5] != M) q->sizes[5]++;
  q->sizes[6] = P;
  q->sizes[7] = Q;
  if ((extra == 0) && (q->sizes[1] == 1)) return -2;
  for (i = 0; i < 8; i++) 
    if (q->sizes[i] > MAXSIDE-1) return -3;

  /* set the 26 critical points -- 16 exterior */

  cpts[ 0][0] = uv[0];
  cpts[ 0][1] = uv[1];
  cipt[ 0]    = indices[0];
  len  = q->sizes[0];
  cpts[ 4][0] = uv[2*len  ];
  cpts[ 4][1] = uv[2*len+1];
  cipt[ 4]    = indices[len];
  len += q->sizes[1];
  cpts[13][0] = uv[2*len  ];
  cpts[13][1] = uv[2*len+1];
  cipt[13]    = indices[len];
  len += q->sizes[2];
  cpts[19][0] = uv[2*len  ];
  cpts[19][1] = uv[2*len+1];
  cipt[19]    = indices[len];
  len += q->sizes[3];
  cpts[20][0] = uv[2*len  ];
  cpts[20][1] = uv[2*len+1];
  cipt[20]    = indices[len];
  len += q->sizes[6];
  cpts[21][0] = uv[2*len  ];
  cpts[21][1] = uv[2*len+1];
  cipt[21]    = indices[len];
  len += q->sizes[4];
  cpts[22][0] = uv[2*len  ];
  cpts[22][1] = uv[2*len+1];
  cipt[22]    = indices[len];
  len += q->sizes[5];
  cpts[24][0] = uv[2*len  ];
  cpts[24][1] = uv[2*len+1];
  cipt[24]    = indices[len];
  len += q->sizes[2];
  cpts[23][0] = uv[2*len  ];
  cpts[23][1] = uv[2*len+1];
  cipt[23]    = indices[len];
  len += q->sizes[7];
  cpts[18][0] = uv[2*len  ];
  cpts[18][1] = uv[2*len+1];
  cipt[18]    = indices[len];
  len += q->sizes[1];
  if (extra != -1) len++;
  cpts[25][0] = uv[2*len  ];
  cpts[25][1] = uv[2*len+1];
  cipt[25]    = indices[len];
  len += q->sizes[6];
  cpts[11][0] = uv[2*len  ];
  cpts[11][1] = uv[2*len+1];
  cipt[11]    = indices[len];
  len += q->sizes[7];
  cpts[ 7][0] = uv[2*len  ];
  cpts[ 7][1] = uv[2*len+1];
  cipt[ 7]    = indices[len];
  len += q->sizes[0];
  cpts[ 3][0] = uv[2*len  ];
  cpts[ 3][1] = uv[2*len+1];
  cipt[ 3]    = indices[len];
  len += q->sizes[5];
  cpts[ 2][0] = uv[2*len  ];
  cpts[ 2][1] = uv[2*len+1];
  cipt[ 2]    = indices[len];
  len += q->sizes[4];
  cpts[ 1][0] = uv[2*len  ];
  cpts[ 1][1] = uv[2*len+1];
  cipt[ 1]    = indices[len];

  /* guess the interior */
  for (j = 0; j < 10; j++) 
    cpts[interior[j][0]][0] = cpts[interior[j][0]][1] = 0.0;
  for (i = 0; i < 10; i++)
    for (j = 0; j < 10; j++) {
      sums[0] = sums[1] = 0.0;
      for (len = k = 0; k < 5; k++) {
        if (interior[j][k+1] < 0) continue;
        sums[0] += cpts[interior[j][k+1]][0];
        sums[1] += cpts[interior[j][k+1]][1];
        len++;
      }
      cpts[interior[j][0]][0] = sums[0]/len;
      cpts[interior[j][0]][1] = sums[1]/len;
    }

  if (q->fillT == 0) {
    len       = EG_getVertCnt(q, 17, blocks);
    q->vpatch = (int *) EG_alloc(len*sizeof(int));
    if (q->vpatch == NULL) return -1;
  }
  
  /* allocate our temporary storage */

  len      = MAX(elens[1], elens[2]);
  q->quads = (Quad *) EG_alloc(len*len*sizeof(Quad));
  if (q->quads == NULL) {
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }
  q->verts = (Node *) EG_alloc((len+1)*(len+1)*sizeof(Node));
  if (q->verts == NULL) {
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  /* initialize the vertices */

  q->nvert = elens[0] + elens[1] + elens[2] + elens[3];
  for (i = 0; i < q->nvert; i++) {
    j = indices[i];
    q->verts[j].uv[0] = uv[2*i  ];
    q->verts[j].uv[1] = uv[2*i+1];
    q->verts[j].area  = -1.0;
  }
  for (i = 0; i < 10; i++) {
    j = interior[i][0];
    q->verts[q->nvert].uv[0] = cpts[j][0];
    q->verts[q->nvert].uv[1] = cpts[j][1];
    q->verts[q->nvert].area  = 0.0;
    cipt[j]                  = q->nvert;
    q->nvert++;
  }

  /* set the exterior block sides */

  for (i = 0; i < 42; i++) sideptr[i] = NULL;
  for (j = i = 0; i < 16; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    if (i >= 7) {
      if ((i == 9) && (extra != -1)) {
        for (k = len-1; k > 0; k--, j++) {
          if ((k == len/2) && (extra == 0)) {
            extra = j;
            j++;
          }
          sideptr[i][k] = indices[j];
        }
        sideptr[i][0] = indices[j];
      } else {
        for (k = len-1; k > 0; k--, j++) sideptr[i][k] = indices[j];
        if (i == 15) {
          sideptr[i][0] = indices[0];
        } else {
          sideptr[i][0] = indices[j];
        }
      }
    } else {
      for (k = 0; k < len-1; k++, j++) sideptr[i][k] = indices[j];
      sideptr[i][len-1] = indices[j];
    }
  }

  /* do the interior sides */

  for (i = 16; i < 42; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    i0 = sides[i][1];
    i1 = sides[i][2];
    sideptr[i][0] = cipt[i0];
    for (j = 1; j < len-1; j++) {
      q->verts[q->nvert].uv[0] = cpts[i0][0] +
                                   j*(cpts[i1][0]-cpts[i0][0])/(len-1);
      q->verts[q->nvert].uv[1] = cpts[i0][1] +
                                   j*(cpts[i1][1]-cpts[i0][1])/(len-1);
      q->verts[q->nvert].area  = 0.0;
      sideptr[i][j]            = q->nvert;
      q->nvert++;
    }
    sideptr[i][len-1] = cipt[i1];
  }

  /* start filling the quads by specifying the 17 blocks */

  EG_setQuads(q, 17, blocks, sideptr);

  /* free up our integrated sides */
  for (i = 0; i < 42; i++) EG_free(sideptr[i]);

  /* get the actual storage that we return the data with */

  uvb = (double *) EG_alloc(2*q->nvert*sizeof(double));
  if (uvb == NULL) {
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;    
  }

  /* calculate the actual coordinates */

  len = elens[1]*elens[2];
  EG_smoothQuads(face, q, len, nsp);

  /* fill the memory to be returned */

  for (j = 0; j < q->nvert; j++) {
    uvb[2*j  ] = q->verts[j].uv[0];
    uvb[2*j+1] = q->verts[j].uv[1];
  }
  
  if (q->fillT != 0) {
    len = 6*q->nquad;
    if (extra != -1) len += 3;
    q->ntris = len/3;
    q->tris  = (int *) EG_alloc(len*sizeof(int));
    if (q->tris == NULL) {
      EG_free(uvb);
      EG_free(q->verts);
      EG_free(q->quads);
      return -1;
    }
    for (j = 0; j < q->nquad; j++) {
      q->tris[6*j  ] = q->quads[j].nodes[0];
      q->tris[6*j+1] = q->quads[j].nodes[1];
      q->tris[6*j+2] = q->quads[j].nodes[2];
      q->tris[6*j+3] = q->quads[j].nodes[0];
      q->tris[6*j+4] = q->quads[j].nodes[2];
      q->tris[6*j+5] = q->quads[j].nodes[3];
    }
    
    len = 2*q->nquad;
    if (extra != -1) {
      k = -1;
      for (j = 0; j < len; j++) {
        i = 0;
        if ((q->tris[3*j  ] == indices[extra-1]) ||
            (q->tris[3*j+1] == indices[extra-1]) ||
            (q->tris[3*j+2] == indices[extra-1])) i++;
        if ((q->tris[3*j  ] == indices[extra+1]) ||
            (q->tris[3*j+1] == indices[extra+1]) ||
            (q->tris[3*j+2] == indices[extra+1])) i++;
        if (i == 2) {
          k = j;
          break;
        }
      }
      if (k == -1) {
        EG_free(uvb);
        EG_free(q->tris);
        EG_free(q->verts);
        EG_free(q->quads);
#ifdef DEBUG
        printf(" Problem in Fixup!\n");
#endif
        return -4;
      }
      q->tris[3*len  ] = q->tris[3*k  ];
      q->tris[3*len+1] = q->tris[3*k+1];
      q->tris[3*len+2] = q->tris[3*k+2];
      if ((q->tris[3*k  ] != indices[extra-1]) &&
          (q->tris[3*k  ] != indices[extra+1])) {
        q->tris[3*k  +2] = indices[extra];
        q->tris[3*len+1] = indices[extra];
      } else if ((q->tris[3*k+1] != indices[extra-1]) &&
                 (q->tris[3*k+1] != indices[extra+1])) {
        q->tris[3*k  +2] = indices[extra];
        q->tris[3*len  ] = indices[extra];
      } else {
        q->tris[3*k  +1] = indices[extra];
        q->tris[3*len  ] = indices[extra];
      }
      len++;
    }
  }

  /* cleanup and exit */

  *npts = q->nvert;
  *uvs  = uvb;
  EG_free(q->verts);

  return 0;
}


/* No P case */

static int
EG_quadFillQ(const egObject *face, qFill *q, int nsp, int *indices, int *elens,
             double *uv, int *npts, double **uvs)
{
  int    i, j, k, len, N, M, P, Q, extra = -1;
  int    i0, i1, cipt[20], *sideptr[31];
  double sums[2], cpts[20][2], *uvb;
  static int sides[31][3] = {  0,  0,  1,   1,  1,  2,   2,  2,  3, 
                               3,  3,  4,   4,  4,  5,   5,  5,  6,
                               2,  7,  6,   7,  8,  7,   1,  9,  8,
                               7, 10,  9,   0, 11, 10,   5, 12, 11,
                               4, 13, 12,   3,  0, 13,   3,  1, 14,
                               0, 13, 14,   3,  2, 15,   1, 14, 15,
                               2, 15,  4,   4, 14, 16,   0, 12, 16,
                               4, 15, 17,   1, 16, 17,   2, 17,  5,
                               5, 16, 10,   7, 16, 18,   5, 18,  9,
                               7, 19, 17,   1, 18, 19,   5, 17,  7,
                               5, 19,  8 };

                      /*       size    sides           */
  static int blocks[12][6] = { 0, 3,   13, 14,  0, 15,
                               1, 3,   14, 16,  1, 17,
                               2, 3,   16,  3,  2, 18,
                               0, 4,   12, 19, 15, 20,
                               1, 4,   19, 21, 17, 22,
                               2, 4,   21,  4, 18, 23,
                               0, 5,   11, 24, 20, 10,
                               7, 5,   24, 26, 25,  9,
                               1, 7,   25,-27, 22, 28,
                               7, 5,   30, 29, 27,  7,
                               2, 5,   29,  5, 23,  6,
                               1, 5,   26, 30, 28,  8 };

                      /*      center   stencil nodes     */     
  static int interior[6][6] = { 14,    1, 13, 15, 16, -1,
                                15,    2, 14,  4, 17, -1,
                                16,   14, 12, 10, 17, 18,
                                17,   16, 19, 15,  7,  5,
                                18,   16,  9, 19, -1, -1,
                                19,   18, 17,  8, -1, -1 };

  /* get and check sizes */

  N =  elens[0];
  M =  elens[3];
  P =  elens[1] - M;
  Q = (elens[2] - N - P)/2;
  if (Q*2 != elens[2]-N-P)
    if (q->fillT == 0) {
      if (q->outLevel > 0) {
        printf(" EGADS Info: Q Case off by 1 - %d %d  %d %d\n",
               elens[0], elens[2], elens[1], elens[3]);
        printf("             N = %d, M = %d, P = %d, Q = %d\n",
               N, M, P, Q);
      }
      return -2;
    } else {
      extra = 0;
#ifdef DEBUG
      printf(" EGADS Info: Q Case - %d %d  %d %d\n",
             elens[0], elens[2], elens[1], elens[3]);
      printf("             N = %d, M = %d, P = %d, Q = %d\n",
             N, M, P, Q);
#endif
    }
  if (P != 0) return -2;

  q->sizes[0] = q->sizes[1] = q->sizes[2] = N/3;
  if (q->sizes[0]+q->sizes[1]+q->sizes[2] != N) q->sizes[0]++;
  if (q->sizes[0]+q->sizes[1]+q->sizes[2] != N) q->sizes[2]++;
  q->sizes[3] = q->sizes[4] = q->sizes[5] = M/3;
  if (q->sizes[3]+q->sizes[4]+q->sizes[5] != M) q->sizes[3]++;
  if (q->sizes[3]+q->sizes[4]+q->sizes[5] != M) q->sizes[5]++;
  q->sizes[6] = P;
  q->sizes[7] = Q;
  if ((extra == 0) && (q->sizes[1] == 1)) return -2;
  for (i = 0; i < 8; i++)
    if (q->sizes[i] > MAXSIDE-1) return -3;

  /* set the 20 critical points -- 14 exterior */

  cpts[ 0][0] = uv[0];
  cpts[ 0][1] = uv[1];
  cipt[ 0]    = indices[0];;
  len  = q->sizes[0];
  cpts[ 1][0] = uv[2*len  ];
  cpts[ 1][1] = uv[2*len+1];
  cipt[ 1]    = indices[len];
  len += q->sizes[1];
  cpts[ 2][0] = uv[2*len  ];
  cpts[ 2][1] = uv[2*len+1];
  cipt[ 2]    = indices[len];
  len += q->sizes[2];
  cpts[ 3][0] = uv[2*len  ];
  cpts[ 3][1] = uv[2*len+1];
  cipt[ 3]    = indices[len];
  len += q->sizes[3];
  cpts[ 4][0] = uv[2*len  ];
  cpts[ 4][1] = uv[2*len+1];
  cipt[ 4]    = indices[len];
  len += q->sizes[4];
  cpts[ 5][0] = uv[2*len  ];
  cpts[ 5][1] = uv[2*len+1];
  cipt[ 5]    = indices[len];
  len += q->sizes[5];
  cpts[ 6][0] = uv[2*len  ];
  cpts[ 6][1] = uv[2*len+1];
  cipt[ 6]    = indices[len];
  len += q->sizes[2];
  cpts[ 7][0] = uv[2*len  ];
  cpts[ 7][1] = uv[2*len+1];
  cipt[ 7]    = indices[len];
  len += q->sizes[7];
  cpts[ 8][0] = uv[2*len  ];
  cpts[ 8][1] = uv[2*len+1];
  cipt[ 8]    = indices[len];
  len += q->sizes[1];
  if (extra != -1) len++;
  cpts[ 9][0] = uv[2*len  ];
  cpts[ 9][1] = uv[2*len+1];
  cipt[ 9]    = indices[len];
  len += q->sizes[7];
  cpts[10][0] = uv[2*len  ];
  cpts[10][1] = uv[2*len+1];
  cipt[10]    = indices[len];
  len += q->sizes[0];
  cpts[11][0] = uv[2*len  ];
  cpts[11][1] = uv[2*len+1];
  cipt[11]    = indices[len];
  len += q->sizes[5];
  cpts[12][0] = uv[2*len  ];
  cpts[12][1] = uv[2*len+1];
  cipt[12]    = indices[len];
  len += q->sizes[4];
  cpts[13][0] = uv[2*len  ];
  cpts[13][1] = uv[2*len+1];
  cipt[13]    = indices[len];

  /* guess the interior */
  for (j = 0; j < 6; j++) 
    cpts[interior[j][0]][0] = cpts[interior[j][0]][1] = 0.0;
  for (i = 0; i < 10; i++)
    for (j = 0; j < 6; j++) {
      sums[0] = sums[1] = 0.0;
      for (len = k = 0; k < 5; k++) {
        if (interior[j][k+1] < 0) continue;
        sums[0] += cpts[interior[j][k+1]][0];
        sums[1] += cpts[interior[j][k+1]][1];
        len++;
      }
      cpts[interior[j][0]][0] = sums[0]/len;
      cpts[interior[j][0]][1] = sums[1]/len;
    }
 
  if (q->fillT == 0) {
    len       = EG_getVertCnt(q, 12, blocks);
    q->vpatch = (int *) EG_alloc(len*sizeof(int));
    if (q->vpatch == NULL) return -1;
  }

  /* allocate our temporary storage */

  len      = MAX(elens[1], elens[2]);
  q->quads = (Quad *) EG_alloc(len*len*sizeof(Quad));
  if (q->quads == NULL) {
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }
  q->verts = (Node *) EG_alloc((len+1)*(len+1)*sizeof(Node));
  if (q->verts == NULL) {
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  /* initialize the vertices */

  q->nvert = elens[0] + elens[1] + elens[2] + elens[3];
  for (i = 0; i < q->nvert; i++) {
    j = indices[i];
    q->verts[j].uv[0] = uv[2*i  ];
    q->verts[j].uv[1] = uv[2*i+1];
    q->verts[j].area  = -1.0;
  }
  for (i = 0; i < 6; i++) {
    j = interior[i][0];
    q->verts[q->nvert].uv[0] = cpts[j][0];
    q->verts[q->nvert].uv[1] = cpts[j][1];
    q->verts[q->nvert].area  = 0.0;
    cipt[j]                  = q->nvert;
    q->nvert++;
  }

  /* set the exterior block sides */

  for (i = 0; i < 31; i++) sideptr[i] = NULL;
  for (j = i = 0; i < 14; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    if (i >= 6) {
      if ((i == 8) && (extra != -1)) {
        for (k = len-1; k > 0; k--, j++) {
          if ((k == len/2) && (extra == 0)) {
            extra = j;
            j++;
          }
          sideptr[i][k] = indices[j];
        }
        sideptr[i][0] = indices[j];
      } else {
        for (k = len-1; k > 0; k--, j++) sideptr[i][k] = indices[j];
        if (i == 13) {
          sideptr[i][0] = indices[0];
        } else {
          sideptr[i][0] = indices[j];
        }
      }
    } else {
      for (k = 0; k < len-1; k++, j++) sideptr[i][k] = indices[j];
      sideptr[i][len-1] = indices[j];
    }
  }

  /* do the interior sides */

  for (i = 14; i < 31; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    i0 = sides[i][1];
    i1 = sides[i][2];
    sideptr[i][0] = cipt[i0];
    for (j = 1; j < len-1; j++) {
      q->verts[q->nvert].uv[0] = cpts[i0][0] +
                                   j*(cpts[i1][0]-cpts[i0][0])/(len-1);
      q->verts[q->nvert].uv[1] = cpts[i0][1] +
                                   j*(cpts[i1][1]-cpts[i0][1])/(len-1);
      q->verts[q->nvert].area  = 0.0;
      sideptr[i][j]            = q->nvert;
      q->nvert++;
    }
    sideptr[i][len-1] = cipt[i1];
  }

  /* start filling the quads by specifying the 12 blocks */

  EG_setQuads(q, 12, blocks, sideptr);

  /* free up our integrated sides */
  for (i = 0; i < 31; i++) EG_free(sideptr[i]);

  /* get the actual storage that we return the data with */

  uvb = (double *) EG_alloc(2*q->nvert*sizeof(double));
  if (uvb == NULL) {
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;    
  }

  /* calculate the actual coordinates */

  len = elens[1]*elens[2];
  EG_smoothQuads(face, q, len, nsp);

  /* fill the memory to be returned */

  for (j = 0; j < q->nvert; j++) {
    uvb[2*j  ] = q->verts[j].uv[0];
    uvb[2*j+1] = q->verts[j].uv[1];
  }
  
  if (q->fillT != 0) {
    len = 6*q->nquad;
    if (extra != -1) len += 3;
    q->ntris = len/3;
    q->tris  = (int *) EG_alloc(len*sizeof(int));
    if (q->tris == NULL) {
      EG_free(uvb);
      EG_free(q->verts);
      EG_free(q->quads);
      return -1;
    }

    for (j = 0; j < q->nquad; j++) {
      q->tris[6*j  ] = q->quads[j].nodes[0];
      q->tris[6*j+1] = q->quads[j].nodes[1];
      q->tris[6*j+2] = q->quads[j].nodes[2];
      q->tris[6*j+3] = q->quads[j].nodes[0];
      q->tris[6*j+4] = q->quads[j].nodes[2];
      q->tris[6*j+5] = q->quads[j].nodes[3];
    }
    
    /* fixup for an extra point */
    
    len = 2*q->nquad;
    if (extra != -1) {
      k = -1;
      for (j = 0; j < len; j++) {
        i = 0;
        if ((q->tris[3*j  ] == indices[extra-1]) ||
            (q->tris[3*j+1] == indices[extra-1]) ||
            (q->tris[3*j+2] == indices[extra-1])) i++;
        if ((q->tris[3*j  ] == indices[extra+1]) ||
            (q->tris[3*j+1] == indices[extra+1]) ||
            (q->tris[3*j+2] == indices[extra+1])) i++;
        if (i == 2) {
          k = j;
          break;
        }
      }
      if (k == -1) {
        EG_free(uvb);
        EG_free(q->tris);
        EG_free(q->verts);
        EG_free(q->quads);
        return -4;
      }
      q->tris[3*len  ] = q->tris[3*k  ];
      q->tris[3*len+1] = q->tris[3*k+1];
      q->tris[3*len+2] = q->tris[3*k+2];
      if ((q->tris[3*k  ] != indices[extra-1]) &&
          (q->tris[3*k  ] != indices[extra+1])) {
        q->tris[3*k  +2] = indices[extra];
        q->tris[3*len+1] = indices[extra];
      } else if ((q->tris[3*k+1] != indices[extra-1]) &&
                 (q->tris[3*k+1] != indices[extra+1])) {
        q->tris[3*k  +2] = indices[extra];
        q->tris[3*len  ] = indices[extra];
      } else {
        q->tris[3*k  +1] = indices[extra];
        q->tris[3*len  ] = indices[extra];
      }
      len++;
    }
  }

  /* cleanup and exit */

  *npts = q->nvert;
  *uvs  = uvb;
  EG_free(q->verts);
  
  return 0;
}


/* No Q case */

static int
EG_quadFillP(const egObject *face, qFill *q, int nsp, int *indices, int *elens,
             double *uv, int *npts, double **uvs)
{
  int    i, j, k, len, N, M, P, Q, extra = -1;
  int    i0, i1, cipt[21], *sideptr[33];
  double sums[2], cpts[21][2], *uvb;
  static int sides[33][3] = {  0,  0,  1,   1,  1,  2,   2,  2,  3, 
                               3,  3,  4,   6,  4,  5,   4,  5,  6,
                               5,  6,  7,   2,  8,  7,   1,  9,  8,
                               6, 10,  9,   0, 11, 10,   5, 12, 11,
                               4, 13, 12,   3,  0, 13,   3,  1, 14,
                               0, 13, 14,   3,  2, 15,   1, 14, 15,
                               2, 15,  4,   4, 14, 18,   0, 12, 18,
                               6, 14, 16,   4, 16, 19,   6, 15, 17,
                               1, 16, 17,   2, 17,  5,   1, 19, 20,
                               2, 20,  6,   5, 18, 10,   6, 18, 19,
                               5, 19,  9,   4, 17, 20,   5, 20,  8 };

                      /*       size    sides           */
  static int blocks[13][6] = { 0, 3,   13, 14,  0, 15,
                               1, 3,   14, 16,  1, 17,
                               2, 3,   16,  3,  2, 18,
                               0, 4,   12, 19, 15, 20,
                               6, 4,   19, 22, 21, 29,
                               1, 6,   21, 23, 17, 24,
                               2, 6,   23,  4, 18, 25,
                               1, 4,   22, 31, 24, 26,
                               2, 4,   31,  5, 25, 27,
                               0, 5,   11, 28, 20, 10,
                               6, 5,   28, 30, 29,  9,
                               1, 5,   30, 32, 26,  8,
                               2, 5,   32,  6, 27,  7 };

                      /*      center   stencil nodes     */     
  static int interior[7][6] = { 14,    1, 13, 18, 16, 15,
                                15,    2,  4, 14, 17, -1,
                                16,   14, 17, 19, -1, -1,
                                17,    5, 15, 16, 20, -1,
                                18,   10, 12, 14, 19, -1,
                                19,    9, 16, 18, 20, -1,
                                20,    6,  8, 17, 19, -1 };

  /* get and check sizes */

  N =  elens[0];
  M =  elens[3];
  P =  elens[1] - M;
  Q = (elens[2] - N - P)/2;
  if (Q*2 != elens[2]-N-P)
    if (q->fillT == 0) {
      if (q->outLevel > 0) {
        printf(" EGADS Info: P Case off by 1 - %d %d  %d %d\n",
               elens[0], elens[2], elens[1], elens[3]);
        printf("             N = %d, M = %d, P = %d, Q = %d\n",
               N, M, P, Q);
      }
      return -2;
    } else {
      extra = 0;
#ifdef DEBUG
      printf(" EGADS Info: P Case - %d %d  %d %d\n",
             elens[0], elens[2], elens[1], elens[3]);
      printf("             N = %d, M = %d, P = %d, Q = %d\n",
             N, M, P, Q);
#endif
    }
  if (Q != 0) return -2;

  q->sizes[0] = q->sizes[1] = q->sizes[2] = N/3;
  if (q->sizes[0]+q->sizes[1]+q->sizes[2] != N) q->sizes[0]++;
  if (q->sizes[0]+q->sizes[1]+q->sizes[2] != N) q->sizes[2]++;
  q->sizes[3] = q->sizes[4] = q->sizes[5] = M/3;
  if (q->sizes[3]+q->sizes[4]+q->sizes[5] != M) q->sizes[3]++;
  if (q->sizes[3]+q->sizes[4]+q->sizes[5] != M) q->sizes[5]++;
  q->sizes[6] = P;
  q->sizes[7] = Q;
  if ((extra == 0) && (q->sizes[1] == 1)) return -2;
  for (i = 0; i < 8; i++)
    if (q->sizes[i] > MAXSIDE-1) return -3;

  /* set the 21 critical points -- 14 exterior */

  cpts[ 0][0] = uv[0];
  cpts[ 0][1] = uv[1];
  cipt[ 0]    = indices[0];
  len  = q->sizes[0];
  cpts[ 1][0] = uv[2*len  ];
  cpts[ 1][1] = uv[2*len+1];
  cipt[ 1]    = indices[len];
  len += q->sizes[1];
  cpts[ 2][0] = uv[2*len  ];
  cpts[ 2][1] = uv[2*len+1];
  cipt[ 2]    = indices[len];
  len += q->sizes[2];
  cpts[ 3][0] = uv[2*len  ];
  cpts[ 3][1] = uv[2*len+1];
  cipt[ 3]    = indices[len];
  len += q->sizes[3];
  cpts[ 4][0] = uv[2*len  ];
  cpts[ 4][1] = uv[2*len+1];
  cipt[ 4]    = indices[len];
  len += q->sizes[6];
  cpts[ 5][0] = uv[2*len  ];
  cpts[ 5][1] = uv[2*len+1];
  cipt[ 5]    = indices[len];
  len += q->sizes[4];
  cpts[ 6][0] = uv[2*len  ];
  cpts[ 6][1] = uv[2*len+1];
  cipt[ 6]    = indices[len];
  len += q->sizes[5];
  cpts[ 7][0] = uv[2*len  ];
  cpts[ 7][1] = uv[2*len+1];
  cipt[ 7]    = indices[len];
  len += q->sizes[2];
  cpts[ 8][0] = uv[2*len  ];
  cpts[ 8][1] = uv[2*len+1];
  cipt[ 8]    = indices[len];
  len += q->sizes[1];
  if (extra != -1) len++;
  cpts[ 9][0] = uv[2*len  ];
  cpts[ 9][1] = uv[2*len+1];
  cipt[ 9]    = indices[len];
  len += q->sizes[6];
  cpts[10][0] = uv[2*len  ];
  cpts[10][1] = uv[2*len+1];
  cipt[10]    = indices[len];
  len += q->sizes[0];
  cpts[11][0] = uv[2*len  ];
  cpts[11][1] = uv[2*len+1];
  cipt[11]    = indices[len];
  len += q->sizes[5];
  cpts[12][0] = uv[2*len  ];
  cpts[12][1] = uv[2*len+1];
  cipt[12]    = indices[len];
  len += q->sizes[4];
  cpts[13][0] = uv[2*len  ];
  cpts[13][1] = uv[2*len+1];
  cipt[13]    = indices[len];

  /* guess the interior */
  for (j = 0; j < 7; j++) 
    cpts[interior[j][0]][0] = cpts[interior[j][0]][1] = 0.0;
  for (i = 0; i < 10; i++)
    for (j = 0; j < 7; j++) {
      sums[0] = sums[1] = 0.0;
      for (len = k = 0; k < 5; k++) {
        if (interior[j][k+1] < 0) continue;
        sums[0] += cpts[interior[j][k+1]][0];
        sums[1] += cpts[interior[j][k+1]][1];
        len++;
      }
      cpts[interior[j][0]][0] = sums[0]/len;
      cpts[interior[j][0]][1] = sums[1]/len;
    }
  
  if (q->fillT == 0) {
    len       = EG_getVertCnt(q, 13, blocks);
    q->vpatch = (int *) EG_alloc(len*sizeof(int));
    if (q->vpatch == NULL) return -1;
  }

  /* allocate our temporary storage */

  len      = MAX(elens[1], elens[2]);
  q->quads = (Quad *) EG_alloc(len*len*sizeof(Quad));
  if (q->quads == NULL) {
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }
  q->verts = (Node *) EG_alloc((len+1)*(len+1)*sizeof(Node));
  if (q->verts == NULL) {
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  /* initialize the vertices */

  q->nvert = elens[0] + elens[1] + elens[2] + elens[3];
  for (i = 0; i < q->nvert; i++) {
    j = indices[i];
    q->verts[j].uv[0] = uv[2*i  ];
    q->verts[j].uv[1] = uv[2*i+1];
    q->verts[j].area  = -1.0;
  }
  for (i = 0; i < 7; i++) {
    j = interior[i][0];
    q->verts[q->nvert].uv[0] = cpts[j][0];
    q->verts[q->nvert].uv[1] = cpts[j][1];
    q->verts[q->nvert].area  = 0.0;
    cipt[j]                  = q->nvert;
    q->nvert++;
  }

  /* set the exterior block sides */

  for (i = 0; i < 33; i++) sideptr[i] = NULL;
  for (j = i = 0; i < 14; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    if (i >= 7) {
      if ((i == 8) && (extra != -1)) {
        for (k = len-1; k > 0; k--, j++) {
          if ((k == len/2) && (extra == 0)) {
            extra = j;
            j++;
          }
          sideptr[i][k] = indices[j];
        }
        sideptr[i][0] = indices[j];
      } else {
        for (k = len-1; k > 0; k--, j++) sideptr[i][k] = indices[j];
        if (i == 13) {
          sideptr[i][0] = indices[0];
        } else {
          sideptr[i][0] = indices[j];
        }
      }
    } else {
      for (k = 0; k < len-1; k++, j++) sideptr[i][k] = indices[j];
      sideptr[i][len-1] = indices[j];
    }
  }

  /* do the interior sides */

  for (i = 14; i < 33; i++) {
    len = q->sizes[sides[i][0]] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    i0 = sides[i][1];
    i1 = sides[i][2];
    sideptr[i][0] = cipt[i0];
    for (j = 1; j < len-1; j++) {
      q->verts[q->nvert].uv[0] = cpts[i0][0] +
                                   j*(cpts[i1][0]-cpts[i0][0])/(len-1);
      q->verts[q->nvert].uv[1] = cpts[i0][1] +
                                   j*(cpts[i1][1]-cpts[i0][1])/(len-1);
      q->verts[q->nvert].area  = 0.0;
      sideptr[i][j]            = q->nvert;
      q->nvert++;
    }
    sideptr[i][len-1] = cipt[i1];
  }

  /* start filling the quads by specifying the 13 blocks */

  EG_setQuads(q, 13, blocks, sideptr);

  /* free up our integrated sides */
  for (i = 0; i < 33; i++) EG_free(sideptr[i]);

  /* get the actual storage that we return the data with */

  uvb = (double *) EG_alloc(2*q->nvert*sizeof(double));
  if (uvb == NULL) {
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;    
  }

  /* calculate the actual coordinates */

  len = elens[1]*elens[2];
  EG_smoothQuads(face, q, len, nsp);

  /* fill the memory to be returned */

  for (j = 0; j < q->nvert; j++) {
    uvb[2*j  ] = q->verts[j].uv[0];
    uvb[2*j+1] = q->verts[j].uv[1];
  }
  
  if (q->fillT != 0) {
    len = 6*q->nquad;
    if (extra != -1) len += 3;
    q->ntris = len/3;
    q->tris  = (int *) EG_alloc(len*sizeof(int));
    if (q->tris == NULL) {
      EG_free(uvb);
      EG_free(q->verts);
      EG_free(q->quads);
      return -1;
    }
    for (j = 0; j < q->nquad; j++) {
      q->tris[6*j  ] = q->quads[j].nodes[0];
      q->tris[6*j+1] = q->quads[j].nodes[1];
      q->tris[6*j+2] = q->quads[j].nodes[2];
      q->tris[6*j+3] = q->quads[j].nodes[0];
      q->tris[6*j+4] = q->quads[j].nodes[2];
      q->tris[6*j+5] = q->quads[j].nodes[3];
    }
    
    /* fixup for an extra point */
    
    len = 2*q->nquad;
    if (extra != -1) {
      k = -1;
      for (j = 0; j < len; j++) {
        i = 0;
        if ((q->tris[3*j  ] == indices[extra-1]) ||
            (q->tris[3*j+1] == indices[extra-1]) ||
            (q->tris[3*j+2] == indices[extra-1])) i++;
        if ((q->tris[3*j  ] == indices[extra+1]) ||
            (q->tris[3*j+1] == indices[extra+1]) ||
            (q->tris[3*j+2] == indices[extra+1])) i++;
        if (i == 2) {
          k = j;
          break;
        }
      }
      if (k == -1) {
        EG_free(uvb);
        EG_free(q->tris);
        EG_free(q->verts);
        EG_free(q->quads);
        return -4;
      }
      q->tris[3*len  ] = q->tris[3*k  ];
      q->tris[3*len+1] = q->tris[3*k+1];
      q->tris[3*len+2] = q->tris[3*k+2];
      if ((q->tris[3*k  ] != indices[extra-1]) &&
          (q->tris[3*k  ] != indices[extra+1])) {
        q->tris[3*k  +2] = indices[extra];
        q->tris[3*len+1] = indices[extra];
      } else if ((q->tris[3*k+1] != indices[extra-1]) &&
                 (q->tris[3*k+1] != indices[extra+1])) {
        q->tris[3*k  +2] = indices[extra];
        q->tris[3*len  ] = indices[extra];
      } else {
        q->tris[3*k  +1] = indices[extra];
        q->tris[3*len  ] = indices[extra];
      }
      len++;
    }
  }

  /* cleanup and exit */

  *npts = q->nvert;
  *uvs  = uvb;
  EG_free(q->verts);

  return 0;
}


/* TFI case */

static int
EG_quadFillT(qFill *q, int *elens, double *uv, int *npts, double **uvs)
{
  int    i, j, k, m, nx, ny, len, ilast, extra = -1, eside = -1;
  int    cipt[4], *sideptr[4];
  int    ll, lr, ur, ul, j0, jm, i0, im, iv, sav;
  double et, xi, *uvb, *uv0, *uv1, *uv2, *abasis[2];
  Node   *verts;

  nx = elens[0];
  ny = elens[1];
  if (q->fillT != 0) {
    if (elens[2] < nx) {
      nx    = elens[2];
      extra = elens[0]/2;
      eside = 0;
    } else if (elens[2] > nx) {
      extra = elens[0] + elens[1] + elens[2]/2;
      eside = 2;
    }

    if (elens[3] < ny) {
      ny    = elens[3];
      extra = elens[0] + elens[1]/2;
      eside = 1;
    } else if (elens[3] > ny) {
      extra = elens[0] + elens[1] + elens[2] + elens[3]/2;
      eside = 3;
    }
    if ((nx == 1) && (ny == 1)) return -2;
  }

#ifdef DEBUG
  printf(" TFI Case -- Sizes: nx = %d,  ny = %d\n", nx, ny);
#endif
  if (nx >= MAXSIDE) return -3;
  if (ny >= MAXSIDE) return -3;

  cipt[ 0] = 0;
  len      = nx;
  if (eside == 0) len++;
  cipt[ 1] = len;
  len     += ny;
  if (eside == 1) len++;
  cipt[ 2] = len;
  len     += nx;
  if (eside == 2) len++;
  cipt[ 3] = len;
  q->sizes[0] = q->sizes[2] = nx;
  q->sizes[1] = q->sizes[3] = ny;
  
  if (q->fillT == 0) {
    len       = (nx+1)*(ny+1);
    q->vpatch = (int *) EG_alloc(len*sizeof(int));
    if (q->vpatch == NULL) return -1;
    q->npatch      = 1;
    q->patch[0][0] = nx + 1;
    q->patch[0][1] = ny + 1;
  }

  /* allocate our temporary storage */

  len      = nx*ny;
  q->quads = (Quad *) EG_alloc(len*sizeof(Quad));
  if (q->quads == NULL) {
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }
  len      = (nx+1)*(ny+1) + 1;
  q->verts = verts = (Node *) EG_alloc(len*sizeof(Node));
  if (q->verts == NULL) {
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  /* initialize the vertices */

  q->nvert = elens[0] + elens[1] + elens[2] + elens[3];
  for (i = 0; i < q->nvert; i++) {
    verts[i].uv[0] = uv[2*i  ];
    verts[i].uv[1] = uv[2*i+1];
    verts[i].area  = -1.0;
  }

  /* set the exterior block sides */

  for (i = 0; i < 4; i++) sideptr[i] = NULL;
  for (j = i = 0; i < 4; i++) {
    len = q->sizes[i] + 1;
    sideptr[i] = (int *) EG_alloc(len*sizeof(int));
    if (sideptr[i] == NULL) {
      for (k = 0; k < i; k++) EG_free(sideptr[k]);
      EG_free(q->verts);
      EG_free(q->quads);
      if (q->fillT == 0) EG_free(q->vpatch);
      return -1;
    }
    if (i >= 2) {
      for (k = len-1; k > 0; k--, j++) {
        if (j == extra) j++;
        sideptr[i][k] = j;
      }
      sideptr[i][0] = j;
      if (i == 3) sideptr[i][0] = 0;
    } else {
      for (k = 0; k < len-1; k++, j++) {
        if (j == extra) j++;
        sideptr[i][k] = j;
      }
      sideptr[i][len-1] = j;
    }
  }

  /* compute arclength basis functions */

  abasis[0] = (double *) EG_alloc((nx+1)*(ny+1)*sizeof(double));
  if (abasis[0] == NULL) {
    for (i = 0; i < 4; i++) EG_free(sideptr[i]);
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  abasis[1] = (double *) EG_alloc((nx+1)*(ny+1)*sizeof(double));
  if (abasis[1] == NULL) {
    EG_free(abasis[0]);
    for (i = 0; i < 4; i++) EG_free(sideptr[i]);
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;
  }

  /* Finally compute the basis functions */
  EG_arcBasis(q, nx, ny, sideptr, abasis);

  /* create the quads and get coordinates via TFI */

  q->nquad = iv = 0;
  for (i = 0; i < nx+1; i++) q->last[i] = sideptr[0][i];
  ll = cipt[0];
  lr = cipt[1];
  ur = cipt[2];
  ul = cipt[3];
  for (j = 0; j < ny; j++) {
    i0    = sideptr[3][j+1];
    im    = sideptr[1][j+1];
    ilast = i0;
    sav   = q->nquad;
    for (i = 0; i < nx; i++) {
      q->quads[q->nquad].nodes[0] = q->last[i  ];
      q->quads[q->nquad].nodes[1] = q->last[i+1];
      if (j == ny-1) {
        q->quads[q->nquad].nodes[2] = sideptr[2][i+1];
        q->quads[q->nquad].nodes[3] = sideptr[2][i  ];
      } else {
        if (i == nx-1) {
          q->quads[q->nquad].nodes[2] = sideptr[1][j+1];
          q->last[i]  = ilast;
          q->last[nx] = sideptr[1][j+1];
        } else {
          j0 = sideptr[0][i+1];
          jm = sideptr[2][i+1];
          k  = (i+1)*(ny+1) + (j+1);
          xi = abasis[0][k];
          et = abasis[1][k];
          q->quads[q->nquad].nodes[2] = q->nvert;
          verts[q->nvert].uv[0]       = (1.0-xi)            * verts[i0].uv[0] +
                                        (    xi)            * verts[im].uv[0] +
                                                   (1.0-et) * verts[j0].uv[0] +
                                                   (    et) * verts[jm].uv[0] -
                                        (1.0-xi) * (1.0-et) * verts[ll].uv[0] -
                                        (1.0-xi) * (    et) * verts[ul].uv[0] -
                                        (    xi) * (1.0-et) * verts[lr].uv[0] -
                                        (    xi) * (    et) * verts[ur].uv[0];
          verts[q->nvert].uv[1]       = (1.0-xi)            * verts[i0].uv[1] +
                                        (    xi)            * verts[im].uv[1] +
                                                   (1.0-et) * verts[j0].uv[1] +
                                                   (    et) * verts[jm].uv[1] -
                                        (1.0-xi) * (1.0-et) * verts[ll].uv[1] -
                                        (1.0-xi) * (    et) * verts[ul].uv[1] -
                                        (    xi) * (1.0-et) * verts[lr].uv[1] -
                                        (    xi) * (    et) * verts[ur].uv[1];
          verts[q->nvert].area        = 0.0;
          q->nvert++;
        }
        q->quads[q->nquad].nodes[3]   = ilast;
        q->last[i] = ilast;
        ilast      = q->nvert-1;
      }
      q->nquad++;
    }
    if (q->fillT == 0) {
      if (j == 0) {
        q->vpatch[iv] = q->quads[sav].nodes[0];
        iv++;
        for (i = 0; i < nx; i++, iv++)
          q->vpatch[iv] = q->quads[sav+i].nodes[1];
      }
      q->vpatch[iv] = q->quads[sav].nodes[3];
      iv++;
      for (i = 0; i < nx; i++, sav++, iv++)
        q->vpatch[iv] = q->quads[sav].nodes[2];
    }
  }

  /* free up our basis functions */
  EG_free(abasis[1]);
  EG_free(abasis[0]);

  /* free up our integrated sides */
  for (i = 0; i < 4; i++) EG_free(sideptr[i]);

  /* get the actual storage that we return the data with */

  uvb = (double *) EG_alloc(2*q->nvert*sizeof(double));
  if (uvb == NULL) {
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    return -1;    
  }

  /* fill the memory to be returned */

  for (j = 0; j < q->nvert; j++) {
    uvb[2*j  ] = verts[j].uv[0];
    uvb[2*j+1] = verts[j].uv[1];
  }
  for (k = i = 0; i < q->nquad; i++) {
    m   = k;
    uv0 = &uvb[2*q->quads[i].nodes[0]];
    uv1 = &uvb[2*q->quads[i].nodes[1]];
    uv2 = &uvb[2*q->quads[i].nodes[2]];
    if (AREA2D(uv0, uv1, uv2) <= 0.0) k++;
    if (m != k) continue;
    uv1 = &uvb[2*q->quads[i].nodes[2]];
    uv2 = &uvb[2*q->quads[i].nodes[3]];
    if (AREA2D(uv0, uv1, uv2) <= 0.0) k++;
  }
  if (k != 0) {
    if (q->outLevel > 1)
      printf(" Bad mapping: %d non-positive of %d quads\n", k, q->nquad);
    EG_free(uvb);
    EG_free(q->verts);
    EG_free(q->quads);
    if (q->fillT == 0) EG_free(q->vpatch);
    *npts = 0;
    *uvs  = NULL;
    return -6;
  }
  
  if (q->fillT != 0) {
    len = 6*q->nquad;
    if (extra != -1) len += 3;
    q->ntris = len/3;
    q->tris  = (int *) EG_alloc(len*sizeof(int));
    if (q->tris == NULL) {
      EG_free(uvb);
      EG_free(q->verts);
      EG_free(q->quads);
      return -1;
    }
    for (j = 0; j < q->nquad; j++) {
      q->tris[6*j  ] = q->quads[j].nodes[0];
      q->tris[6*j+1] = q->quads[j].nodes[1];
      q->tris[6*j+2] = q->quads[j].nodes[2];
      q->tris[6*j+3] = q->quads[j].nodes[0];
      q->tris[6*j+4] = q->quads[j].nodes[2];
      q->tris[6*j+5] = q->quads[j].nodes[3];
    }
    
    /* fixup for an extra point */
    
    len = 2*q->nquad;
    if (extra != -1) {
      k = -1;
      for (j = 0; j < len; j++) {
        i = 0;
        if ((q->tris[3*j  ] == extra-1) || (q->tris[3*j+1] == extra-1) ||
            (q->tris[3*j+2] == extra-1)) i++;
        if ((q->tris[3*j  ] == extra+1) || (q->tris[3*j+1] == extra+1) ||
            (q->tris[3*j+2] == extra+1)) i++;
        if (i == 2) {
          k = j;
          break;
        }
      }
      if (k == -1) {
        EG_free(uvb);
        EG_free(q->tris);
        EG_free(q->verts);
        EG_free(q->quads);
        return -4;
      }
      q->tris[3*len  ] = q->tris[3*k  ];
      q->tris[3*len+1] = q->tris[3*k+1];
      q->tris[3*len+2] = q->tris[3*k+2];
      if ((q->tris[3*k  ] != extra-1) && (q->tris[3*k  ] != extra+1)) {
        q->tris[3*k  +2] = extra;
        q->tris[3*len+1] = extra;
      } else if ((q->tris[3*k+1] != extra-1) && (q->tris[3*k+1] != extra+1)) {
        q->tris[3*k  +2] = extra;
        q->tris[3*len  ] = extra;
      } else {
        q->tris[3*k  +1] = extra;
        q->tris[3*len  ] = extra;
      }
      len++;
    }
  }

  /* cleanup and exit */

  *npts = q->nvert;
  *uvs  = uvb;
  EG_free(q->verts);
  EG_free(q->quads);

  return 0;
}


/* takes a simple quad loop and fills it with quads based on a sub-blocking
 *       scheme that supports differing sizes per side
 *
 * where: parms[0] - Edge Tol
 *        parms[1] - Side Ratio
 *        parms[2] - # smoothing passes
 *
 *        elens[0] - number of segments on the left side
 *        elens[1] - number of segments on the bottom
 *        elens[2] - number of segments on the right
 *        elens[3] - number of segments on the top side
 *
 *        uv[]     - input (u,v) pairs going around the loop CCW with
 *                   no duplicates at corners starting at UL corner (the
 *                   start of side 0)
 *                   len = 2*(sum of elens)
 *
 *        npts     - total number of resultant points (returned)
 *        uvs      - pointer to (u,v) pairs of all of the points (returned,
 *                   should be free'd);  note -> the first set matches uv
 *        npat     - resultant number of quad patchs (max = 17)
 *        pat      - filled patch sizes (at least 2*17 in len)
 *        vpat     - pointer to patch indices of uvs that support the 
 *                   quad patch(s) (returned, should be free'd)
 *
 * return codes:  0 - success
 *               -1 - malloc error
 *               -2 - elen error
 *               -3 - block side too big
 *               -4 - extra edge not found
 *               -6 - neg area tris
 *               -7 - mismatched sides
 */

int
EG_quadFill(const egObject *face, double *parms, int *elens, double *uv, 
            int *npts, double **uvs, int *npat, int *pats, int **vpats)
{
  int    i, j, k, m, N, M, P, Q, len, ret, lens[4], *indices;
  int    outLevel, iv, nx, save, align = 0, unmap = 0, nsp = 0;
  double *uvx, *uv0, *uv1, *uv2, dist, sav[2], xylim[2][2], slim[4][2][2];
  double edgeTOL = 0.05, sideRAT = 3.0;
  qFill  q;

  *npts    = *npat = 0;
  *uvs     = NULL;
  *vpats   = NULL;
  outLevel = q.outLevel = EG_outLevel(face);
  q.flip   = 1.0;
  q.fillT  = 0;
  q.ntris  = 0;

  /* note: all zeros gives the default values */
  if ((parms[0] >= 0.001) && (parms[0] <= 0.5)) edgeTOL = parms[0];
  if ((parms[1] > 0.0) && (parms[1] <= 1000.0)) sideRAT = parms[1];
  if ((parms[2] > 0.5) && (parms[2] <= 100.0))  nsp     = parms[2]+0.1;

  /* can we use a simple TFI scheme? */

  if ((elens[0] == elens[2]) && (elens[1] == elens[3])) {
    ret = EG_quadFillT(&q, elens, uv, npts, uvs);
    if (ret != EGADS_SUCCESS) return ret;
    *npat   = q.npatch;
    pats[0] = q.patch[0][0];
    pats[1] = q.patch[0][1];
    *vpats  = q.vpatch;
    return EGADS_SUCCESS;
  } else if ((elens[0] == elens[2]) && (abs(elens[1]-elens[3]) == 1)) {
    if (outLevel > 0)
      printf(" EGADS Info: TFI off by 1 on top/bottom!\n");
  } else if ((elens[1] == elens[3]) && (abs(elens[0]-elens[2]) == 1)) {
    if (outLevel > 0)
      printf(" EGADS Info: TFI off by 1 on left/right!\n");
  }

  sav[0] = elens[0];
  sav[1] = elens[2];
  if (MAX(sav[0],sav[1])/MIN(sav[0],sav[1]) > sideRAT) {
    if (outLevel > 0)
      printf(" EGADS Info: Edge ratio0 %lf exceeded: %g %g\n", 
             sideRAT, sav[0], sav[1]);
    return -7;
  }
  sav[0] = elens[1];
  sav[1] = elens[3];
  if (MAX(sav[0],sav[1])/MIN(sav[0],sav[1]) > sideRAT) {
    if (outLevel > 0)
      printf(" EGADS Info: Edge ratio1 %lf exceeded: %g %g\n",
             sideRAT, sav[0], sav[1]);
    return -7;
  }

  /* no -- use our 3 templates */

  len     = elens[0] + elens[1] + elens[2] + elens[3] + 1;
  indices = (int *) EG_alloc(len*sizeof(int));
  if (indices == NULL) return -1;
  uvx = (double *) EG_alloc(4*len*sizeof(double));
  if (uvx == NULL) {
    EG_free(indices);
    return -1;
  }
  for (i = 0; i < len-1; i++) {
    j          = len+i;
    uvx[2*j  ] = uv[2*i  ];
    uvx[2*j+1] = uv[2*i+1];
  }

  /* determine if our quad sides align with U & V */

  xylim[0][0] = xylim[1][0] = uv[0];
  xylim[0][1] = xylim[1][1] = uv[1];
  for (i = 1; i < len-1; i++) {
    if (xylim[0][0] > uv[2*i  ]) xylim[0][0] = uv[2*i  ];
    if (xylim[1][0] < uv[2*i  ]) xylim[1][0] = uv[2*i  ];
    if (xylim[0][1] > uv[2*i+1]) xylim[0][1] = uv[2*i+1];
    if (xylim[1][1] < uv[2*i+1]) xylim[1][1] = uv[2*i+1];
  }
  slim[0][0][0] = slim[1][0][0] = slim[2][0][0] = slim[3][0][0] = xylim[1][0];
  slim[0][1][0] = slim[1][1][0] = slim[2][1][0] = slim[3][1][0] = xylim[0][0];
  slim[0][0][1] = slim[1][0][1] = slim[2][0][1] = slim[3][0][1] = xylim[1][1];
  slim[0][1][1] = slim[1][1][1] = slim[2][1][1] = slim[3][1][1] = xylim[0][1];
  for (i = 0; i <= elens[0]; i++) {
    if (slim[0][0][0] > uv[2*i  ]) slim[0][0][0] = uv[2*i  ];
    if (slim[0][1][0] < uv[2*i  ]) slim[0][1][0] = uv[2*i  ];
    if (slim[0][0][1] > uv[2*i+1]) slim[0][0][1] = uv[2*i+1];
    if (slim[0][1][1] < uv[2*i+1]) slim[0][1][1] = uv[2*i+1];
  }
  for (i = elens[0]; i <= elens[0]+elens[1]; i++) {
    if (slim[1][0][0] > uv[2*i  ]) slim[1][0][0] = uv[2*i  ];
    if (slim[1][1][0] < uv[2*i  ]) slim[1][1][0] = uv[2*i  ];
    if (slim[1][0][1] > uv[2*i+1]) slim[1][0][1] = uv[2*i+1];
    if (slim[1][1][1] < uv[2*i+1]) slim[1][1][1] = uv[2*i+1];
  }
  for (i = elens[0]+elens[1]; i <= elens[0]+elens[1]+elens[2]; i++) {
    if (slim[2][0][0] > uv[2*i  ]) slim[2][0][0] = uv[2*i  ];
    if (slim[2][1][0] < uv[2*i  ]) slim[2][1][0] = uv[2*i  ];
    if (slim[2][0][1] > uv[2*i+1]) slim[2][0][1] = uv[2*i+1];
    if (slim[2][1][1] < uv[2*i+1]) slim[2][1][1] = uv[2*i+1];
  }
  for (i = elens[0]+elens[1]+elens[2]; i < len-1; i++) {
    if (slim[3][0][0] > uv[2*i  ]) slim[3][0][0] = uv[2*i  ];
    if (slim[3][1][0] < uv[2*i  ]) slim[3][1][0] = uv[2*i  ];
    if (slim[3][0][1] > uv[2*i+1]) slim[3][0][1] = uv[2*i+1];
    if (slim[3][1][1] < uv[2*i+1]) slim[3][1][1] = uv[2*i+1];
  }
  if (slim[3][0][0] > uv[0]) slim[3][0][0] = uv[0];
  if (slim[3][1][0] < uv[0]) slim[3][1][0] = uv[0];
  if (slim[3][0][1] > uv[1]) slim[3][0][1] = uv[0];
  if (slim[3][1][1] < uv[1]) slim[3][1][1] = uv[0];
  /* check side range vs face range */
  for (i = 0; i < 4; i++) {
    if (((slim[i][1][0]-slim[i][0][0])/(xylim[1][0]-xylim[0][0]) >= edgeTOL) &&
        ((slim[i][1][1]-slim[i][0][1])/(xylim[1][1]-xylim[0][1]) >= edgeTOL)) {
      align = 1;
      break;
    }
  }
  
  /* fill up our 0 to 1 mapping in UV */

  for (i = 0; i < len; i++) indices[i] = i;
  if ((unmap == 0) && (align == 0)) {
    for (j = i = 0; i < elens[0]; i++) {
      uv[2*j  ] = 0.0;
      uv[2*j+1] = 1.0 - i/((double) elens[0]);
      j++;
    }
    for (i = 0; i < elens[1]; i++) {
      uv[2*j  ] = 0.0 + i/((double) elens[1]);
      uv[2*j+1] = 0.0;
      j++;
    }
    for (i = 0; i < elens[2]; i++) {
      uv[2*j  ] = 1.0;
      uv[2*j+1] = 0.0 + i/((double) elens[2]);
      j++;
    }
    for (i = 0; i < elens[3]; i++) {
      uv[2*j  ] = 1.0 - i/((double) elens[3]);
      uv[2*j+1] = 1.0;
      j++;
    }
  }
  for (i = 0; i < 2*len-2; i++) uvx[i] = uv[i];
#ifdef DEBUG
  printf("quadFill unmap = %d,  align = %d!\n", unmap, align);
#endif

  /* rotate sides to get the biggest delta on side 2 */

  lens[0] = elens[0];
  lens[1] = elens[1];
  lens[2] = elens[2];
  lens[3] = elens[3];
  if (abs(elens[0]-elens[2]) >= abs(elens[1]-elens[3])) {
    if (elens[2] < elens[0]) {
#ifdef DEBUG
      printf("  Side 0  -> Side 2!\n");
#endif
      len = lens[0]+lens[1];
      for (i = 0; i < lens[2]+lens[3]; i++) {
        indices[i] = len+i;
        uvx[2*i  ] = uv[2*(len+i)  ];
        uvx[2*i+1] = uv[2*(len+i)+1];
      }
      len = lens[2]+lens[3];
      for (i = len; i < lens[0]+lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = i-len;
        uvx[2*i  ] = uv[2*(i-len)  ];
        uvx[2*i+1] = uv[2*(i-len)+1];
      }
      lens[2] = elens[0];
      lens[3] = elens[1];
      lens[0] = elens[2];
      lens[1] = elens[3];
    }
  } else {
    if (elens[3] > elens[1]) {
#ifdef DEBUG
      printf("  Side 3  -> Side 2!\n");
#endif
      len = lens[0];
      for (i = 0; i < lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = len+i;
	uvx[2*i  ] = uv[2*(len+i)  ];
	uvx[2*i+1] = uv[2*(len+i)+1];
      }
      len = lens[1]+lens[2]+lens[3];
      for (i = len; i < lens[0]+lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = i-len;
        uvx[2*i  ] = uv[2*(i-len)  ];
        uvx[2*i+1] = uv[2*(i-len)+1];
      }
      lens[3] = elens[0];
      lens[0] = elens[1];
      lens[1] = elens[2];
      lens[2] = elens[3];
    } else {
#ifdef DEBUG
      printf("  Side 1  -> Side 2!\n");
#endif
      len = lens[0]+lens[1]+lens[2];
      for (i = 0; i < lens[3]; i++) {
        indices[i] = len+i;
        uvx[2*i  ] = uv[2*(len+i)  ];
        uvx[2*i+1] = uv[2*(len+i)+1];
      }
      len = lens[3];
      for (i = len; i < lens[0]+lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = i-len;
        uvx[2*i  ] = uv[2*(i-len)  ];
        uvx[2*i+1] = uv[2*(i-len)+1];
      }
      lens[1] = elens[0];
      lens[2] = elens[1];
      lens[3] = elens[2];
      lens[0] = elens[3];
    }
  }

  /* make side 1 bigger than 3 */

  if (lens[1] < lens[3]) {
#ifdef DEBUG
    printf("  Side 1 <-> Side 3!\n");
#endif
    len = lens[0] + lens[1] + lens[2] + lens[3] + 1;
    uvx[2*len-2]   = uvx[0];
    uvx[2*len-1]   = uvx[1];
    indices[len-1] = indices[0];
    for (i = 0; i < len/2; i++) {
      j = len - i - 1;
      sav[0]     = uvx[2*i  ];
      sav[1]     = uvx[2*i+1];
      ret        = indices[i];
      uvx[2*i  ] = uvx[2*j  ];
      uvx[2*i+1] = uvx[2*j+1];
      indices[i] = indices[j];
      uvx[2*j  ] = sav[0];
      uvx[2*j+1] = sav[1];
      indices[j] = ret;
    }
    len--;
    for (j = 0; j < lens[0]; j++) {
      sav[0] = uvx[2*len-2];
      sav[1] = uvx[2*len-1];
      ret    = indices[len-1];
      for(i = len-1; i > 0; i--) {
        uvx[2*i  ] = uvx[2*i-2];
        uvx[2*i+1] = uvx[2*i-1];
        indices[i] = indices[i-1];
      }
      uvx[0]     = sav[0];
      uvx[1]     = sav[1];
      indices[0] = ret;
    }
    i       = lens[3];
    lens[3] = lens[1];
    lens[1] = i;
    q.flip  = -1.0;
  }

  /* get the template & go */

  len =  lens[0]+lens[1]+lens[2]+lens[3];
  N   =  lens[0];
  M   =  lens[3];
  P   =  lens[1] - M;
  Q   = (lens[2] - N - P)/2;
  if ((N < 3) || (M < 3) || (P < 0) || (Q < 0)) {
    for (i = 0; i < len; i++) {
      j         = len+i+1;
      uv[2*i  ] = uvx[2*j  ];
      uv[2*i+1] = uvx[2*j+1];
    }
    EG_free(uvx);
    EG_free(indices);
    if (outLevel > 0)
      printf(" EGADS Info: Too small ->  %d %d (>3)   %d %d\n", 
             N, M, P, Q);
    return -2;
  }

  if (P == 0) {
    ret = EG_quadFillQ(face, &q, nsp, indices, lens, uvx, npts, uvs);
  } else if (Q == 0) {
    ret = EG_quadFillP(face, &q, nsp, indices, lens, uvx, npts, uvs);
  } else {
    ret = EG_quadFillG(face, &q, nsp, indices, lens, uvx, npts, uvs);
  }
  for (i = 0; i < len; i++) {
    j         = len+i+1;
    uv[2*i  ] = uvx[2*j  ];
    uv[2*i+1] = uvx[2*j+1];
  }
  EG_free(uvx);
  EG_free(indices);

  /* remap back to our UV */

  if ((ret == 0) && (unmap == 0) && (align == 0) && (*uvs != NULL)) {
    ret = EG_dQuadTFI(elens, uv, *npts, *uvs);
    if (ret != 0) {
      EG_free(*uvs);
      EG_free(q.quads);
      EG_free(q.vpatch);
      *npts = 0;
      *uvs  = NULL;
    }
  }

  /* fix orientation if flipped direction */

  uvx = *uvs;
  if (uvx == NULL) return -99;
  if ((ret == 0) && (q.flip < 0.0))
    for (iv = k = 0; k < q.npatch; k++) {
      nx = q.patch[k][0];
      for (j = 0; j < q.patch[k][1]; j++) {
        for (i = 0; i < nx/2; i++) {
          m              = nx - i - 1;
          save           = q.vpatch[iv+i];
          q.vpatch[iv+i] = q.vpatch[iv+m];
          q.vpatch[iv+m] = save;
        }
        iv += nx;
      }
    }

  /* make sure we are OK */

  if (ret == 0) {
    for (k = i = 0; i < q.nquad; i++) {
      m    = k;
      uv0  = &uvx[2*q.quads[i].nodes[0]];
      uv1  = &uvx[2*q.quads[i].nodes[1]];
      uv2  = &uvx[2*q.quads[i].nodes[2]];
      dist = AREA2D(uv0, uv1, uv2)*q.flip;
      if (dist*0.0 != 0.0) k++;		/* special Nan, Ind, ... checker */
      if (dist     <= 0.0) k++;
      if (m != k) continue;
      uv1 = &uvx[2*q.quads[i].nodes[2]];
      uv2 = &uvx[2*q.quads[i].nodes[3]];
      dist = AREA2D(uv0, uv1, uv2)*q.flip;
      if (dist*0.0 != 0.0) k++;
      if (dist     <= 0.0) k++;
    }
    EG_free(q.quads);
    if (k != 0) {
      if (outLevel > 0)
        printf(" EGADS Info: Bad mapping - %d non-positive of %d quads\n",
               k, q.nquad);
      EG_free(*uvs);
      EG_free(q.vpatch);
      *npts = 0;
      *uvs  = NULL;
      return -6;
    }
  }

  if (ret == 0) {
    for (k = 0; k < q.npatch; k++) {
      pats[2*k  ] = q.patch[k][0];
      pats[2*k+1] = q.patch[k][1];
    }
    *npat  = q.npatch;
    *vpats = q.vpatch;
  }

  return ret;
}


/* takes a simple quad loop and fills it with tris based on a sub-blocking
 *       scheme that supports differing sizes per side
 *
 * where: parms[0] - Edge Tol
 *        parms[1] - Side Ratio
 *        parms[2] - # smoothing passes
 *
 *        elens[0] - number of segments on the left side
 *        elens[1] - number of segments on the bottom
 *        elens[2] - number of segments on the right
 *        elens[3] - number of segments on the top side
 *
 *        uv[]     - input (u,v) pairs going around the loop CCW with
 *                   no duplicates at corners starting at UL corner (the
 *                   start of side 0)
 *                   len = 2*(sum of elens)
 *
 *        npts     - total number of resultant points (returned)
 *                   on input -- unmap flag
 *        uvs      - pointer to (u,v) pairs of all of the points (returned,
 *                   should be free'd);  note -> the first set matches uv
 *        ntris    - resultant number of tris
 *        tris     - pointer to index triads of uvs that support the
 *                   split quad patch(s) (returned, should be free'd)
 *
 * return codes:  0 - success
 *               -1 - malloc error
 *               -2 - elen error
 *               -3 - block side too big
 *               -4 - extra edge not found
 *               -6 - neg area tris
 *               -7 - mismatched sides
 */

int
EG_quad2tris(long tID, const egObject *face, double *parms, int *elens,
             double *uv, int *npts, double **uvs, int *ntris, int **tris)
{
  int    i, j, k, N, M, P, Q, len, ret, lens[4], *indices;
  int    outLevel, unmap, align = 0, nsp = 0;
  double *uvx, *uv0, *uv1, *uv2, dist, sav[2], xylim[2][2], slim[4][2][2];
  double edgeTOL = 0.05, sideRAT = 3.0;
  qFill  q;
  
  unmap    = *npts;

  *npts    = *ntris = 0;
  *uvs     = NULL;
  *tris    = NULL;
  outLevel = q.outLevel = EG_outLevel(face);
  q.flip   = 1.0;
  q.fillT  = 1;
  q.ntris  = 0;
  q.quads  = NULL;
  
  /* note: all zeros gives the default values */
  if ((parms[0] >= 0.001) && (parms[0] <= 0.5)) edgeTOL = parms[0];
  if ((parms[1] > 0.0) && (parms[1] <= 1000.0)) sideRAT = parms[1];
  if ((parms[2] > 0.5) && (parms[2] <= 100.0))  nsp     = parms[2]+0.1;
  
  /* can we use a simple TFI scheme? */
  
  ret = 9999;
  if ((elens[0] == elens[2]) && (elens[1] == elens[3])) {
    ret = EG_quadFillT(&q, elens, uv, npts, uvs);
  } else if ((elens[0] == elens[2]) && (abs(elens[1]-elens[3]) == 1)) {
    ret = EG_quadFillT(&q, elens, uv, npts, uvs);
  } else if ((elens[1] == elens[3]) && (abs(elens[0]-elens[2]) == 1)) {
    ret = EG_quadFillT(&q, elens, uv, npts, uvs);
  }
  if (ret != 9999) {
    if (ret != EGADS_SUCCESS) {
      if (outLevel > 1)
        printf("%lX EGADS Info: TFI return = %d\n", tID, ret);
      return ret;
    }
    *ntris = q.ntris;
    *tris  = q.tris;
    return ret;
  }
  
  sav[0] = elens[0];
  sav[1] = elens[2];
  if (MAX(sav[0],sav[1])/MIN(sav[0],sav[1]) > sideRAT) {
    if (outLevel > 1)
      printf("%lX EGADS Info: Edge ratio0 %lf exceeded: %g %g\n",
             tID, sideRAT, sav[0], sav[1]);
    return -7;
  }
  sav[0] = elens[1];
  sav[1] = elens[3];
  if (MAX(sav[0],sav[1])/MIN(sav[0],sav[1]) > sideRAT) {
    if (outLevel > 1)
      printf("%lX EGADS Info: Edge ratio1 %lf exceeded: %g %g\n",
             tID, sideRAT, sav[0], sav[1]);
    return -7;
  }
  
  /* no -- use our 3 templates */
  
  len     = elens[0] + elens[1] + elens[2] + elens[3] + 1;
  indices = (int *) EG_alloc(len*sizeof(int));
  if (indices == NULL) return -1;
  uvx = (double *) EG_alloc(4*len*sizeof(double));
  if (uvx == NULL) {
    EG_free(indices);
    return -1;
  }
  for (i = 0; i < len-1; i++) {
    j          = len+i;
    uvx[2*j  ] = uv[2*i  ];
    uvx[2*j+1] = uv[2*i+1];
  }
  
  /* determine if our quad sides align with U & V */
  
  xylim[0][0] = xylim[1][0] = uv[0];
  xylim[0][1] = xylim[1][1] = uv[1];
  for (i = 1; i < len-1; i++) {
    if (xylim[0][0] > uv[2*i  ]) xylim[0][0] = uv[2*i  ];
    if (xylim[1][0] < uv[2*i  ]) xylim[1][0] = uv[2*i  ];
    if (xylim[0][1] > uv[2*i+1]) xylim[0][1] = uv[2*i+1];
    if (xylim[1][1] < uv[2*i+1]) xylim[1][1] = uv[2*i+1];
  }
  slim[0][0][0] = slim[1][0][0] = slim[2][0][0] = slim[3][0][0] = xylim[1][0];
  slim[0][1][0] = slim[1][1][0] = slim[2][1][0] = slim[3][1][0] = xylim[0][0];
  slim[0][0][1] = slim[1][0][1] = slim[2][0][1] = slim[3][0][1] = xylim[1][1];
  slim[0][1][1] = slim[1][1][1] = slim[2][1][1] = slim[3][1][1] = xylim[0][1];
  for (i = 0; i <= elens[0]; i++) {
    if (slim[0][0][0] > uv[2*i  ]) slim[0][0][0] = uv[2*i  ];
    if (slim[0][1][0] < uv[2*i  ]) slim[0][1][0] = uv[2*i  ];
    if (slim[0][0][1] > uv[2*i+1]) slim[0][0][1] = uv[2*i+1];
    if (slim[0][1][1] < uv[2*i+1]) slim[0][1][1] = uv[2*i+1];
  }
  for (i = elens[0]; i <= elens[0]+elens[1]; i++) {
    if (slim[1][0][0] > uv[2*i  ]) slim[1][0][0] = uv[2*i  ];
    if (slim[1][1][0] < uv[2*i  ]) slim[1][1][0] = uv[2*i  ];
    if (slim[1][0][1] > uv[2*i+1]) slim[1][0][1] = uv[2*i+1];
    if (slim[1][1][1] < uv[2*i+1]) slim[1][1][1] = uv[2*i+1];
  }
  for (i = elens[0]+elens[1]; i <= elens[0]+elens[1]+elens[2]; i++) {
    if (slim[2][0][0] > uv[2*i  ]) slim[2][0][0] = uv[2*i  ];
    if (slim[2][1][0] < uv[2*i  ]) slim[2][1][0] = uv[2*i  ];
    if (slim[2][0][1] > uv[2*i+1]) slim[2][0][1] = uv[2*i+1];
    if (slim[2][1][1] < uv[2*i+1]) slim[2][1][1] = uv[2*i+1];
  }
  for (i = elens[0]+elens[1]+elens[2]; i < len-1; i++) {
    if (slim[3][0][0] > uv[2*i  ]) slim[3][0][0] = uv[2*i  ];
    if (slim[3][1][0] < uv[2*i  ]) slim[3][1][0] = uv[2*i  ];
    if (slim[3][0][1] > uv[2*i+1]) slim[3][0][1] = uv[2*i+1];
    if (slim[3][1][1] < uv[2*i+1]) slim[3][1][1] = uv[2*i+1];
  }
  if (slim[3][0][0] > uv[0]) slim[3][0][0] = uv[0];
  if (slim[3][1][0] < uv[0]) slim[3][1][0] = uv[0];
  if (slim[3][0][1] > uv[1]) slim[3][0][1] = uv[0];
  if (slim[3][1][1] < uv[1]) slim[3][1][1] = uv[0];
  /* check side range vs face range */
  for (i = 0; i < 4; i++) {
    if (((slim[i][1][0]-slim[i][0][0])/(xylim[1][0]-xylim[0][0]) >= edgeTOL) &&
        ((slim[i][1][1]-slim[i][0][1])/(xylim[1][1]-xylim[0][1]) >= edgeTOL)) {
      align = 1;
      break;
    }
  }
  
  /* fill up our 0 to 1 mapping in UV */
  
  for (i = 0; i < len; i++) indices[i] = i;
  if ((unmap == 0) && (align == 0)) {
    for (j = i = 0; i < elens[0]; i++) {
      uv[2*j  ] = 0.0;
      uv[2*j+1] = 1.0 - i/((double) elens[0]);
      j++;
    }
    for (i = 0; i < elens[1]; i++) {
      uv[2*j  ] = 0.0 + i/((double) elens[1]);
      uv[2*j+1] = 0.0;
      j++;
    }
    for (i = 0; i < elens[2]; i++) {
      uv[2*j  ] = 1.0;
      uv[2*j+1] = 0.0 + i/((double) elens[2]);
      j++;
    }
    for (i = 0; i < elens[3]; i++) {
      uv[2*j  ] = 1.0 - i/((double) elens[3]);
      uv[2*j+1] = 1.0;
      j++;
    }
  }
  for (i = 0; i < 2*len-2; i++) uvx[i] = uv[i];
#ifdef DEBUG
  printf("quad2tris unmap = %d,  align = %d!\n", unmap, align);
#endif
  
  /* rotate sides to get the biggest delta on side 2 */
  
  lens[0] = elens[0];
  lens[1] = elens[1];
  lens[2] = elens[2];
  lens[3] = elens[3];
  if (abs(elens[0]-elens[2]) >= abs(elens[1]-elens[3])) {
    if (elens[2] < elens[0]) {
#ifdef DEBUG
      printf("  Side 0  -> Side 2!\n");
#endif
      len = lens[0]+lens[1];
      for (i = 0; i < lens[2]+lens[3]; i++) {
        indices[i] = len+i;
        uvx[2*i  ] = uv[2*(len+i)  ];
        uvx[2*i+1] = uv[2*(len+i)+1];
      }
      len = lens[2]+lens[3];
      for (i = len; i < lens[0]+lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = i-len;
        uvx[2*i  ] = uv[2*(i-len)  ];
        uvx[2*i+1] = uv[2*(i-len)+1];
      }
      lens[2] = elens[0];
      lens[3] = elens[1];
      lens[0] = elens[2];
      lens[1] = elens[3];
    }
  } else {
    if (elens[3] > elens[1]) {
#ifdef DEBUG
      printf("  Side 3  -> Side 2!\n");
#endif
      len = lens[0];
      for (i = 0; i < lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = len+i;
	uvx[2*i  ] = uv[2*(len+i)  ];
	uvx[2*i+1] = uv[2*(len+i)+1];
      }
      len = lens[1]+lens[2]+lens[3];
      for (i = len; i < lens[0]+lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = i-len;
        uvx[2*i  ] = uv[2*(i-len)  ];
        uvx[2*i+1] = uv[2*(i-len)+1];
      }
      lens[3] = elens[0];
      lens[0] = elens[1];
      lens[1] = elens[2];
      lens[2] = elens[3];
    } else {
#ifdef DEBUG
      printf("  Side 1  -> Side 2!\n");
#endif
      len = lens[0]+lens[1]+lens[2];
      for (i = 0; i < lens[3]; i++) {
        indices[i] = len+i;
        uvx[2*i  ] = uv[2*(len+i)  ];
        uvx[2*i+1] = uv[2*(len+i)+1];
      }
      len = lens[3];
      for (i = len; i < lens[0]+lens[1]+lens[2]+lens[3]; i++) {
        indices[i] = i-len;
        uvx[2*i  ] = uv[2*(i-len)  ];
        uvx[2*i+1] = uv[2*(i-len)+1];
      }
      lens[1] = elens[0];
      lens[2] = elens[1];
      lens[3] = elens[2];
      lens[0] = elens[3];
    }
  }
  
  /* make side 1 bigger than 3 */
  
  if (lens[1] < lens[3]) {
#ifdef DEBUG
    printf("  Side 1 <-> Side 3!\n");
#endif
    len = lens[0] + lens[1] + lens[2] + lens[3] + 1;
    uvx[2*len-2]   = uvx[0];
    uvx[2*len-1]   = uvx[1];
    indices[len-1] = indices[0];
    for (i = 0; i < len/2; i++) {
      j = len - i - 1;
      sav[0]     = uvx[2*i  ];
      sav[1]     = uvx[2*i+1];
      ret        = indices[i];
      uvx[2*i  ] = uvx[2*j  ];
      uvx[2*i+1] = uvx[2*j+1];
      indices[i] = indices[j];
      uvx[2*j  ] = sav[0];
      uvx[2*j+1] = sav[1];
      indices[j] = ret;
    }
    len--;
    for (j = 0; j < lens[0]; j++) {
      sav[0] = uvx[2*len-2];
      sav[1] = uvx[2*len-1];
      ret    = indices[len-1];
      for(i = len-1; i > 0; i--) {
        uvx[2*i  ] = uvx[2*i-2];
        uvx[2*i+1] = uvx[2*i-1];
        indices[i] = indices[i-1];
      }
      uvx[0]     = sav[0];
      uvx[1]     = sav[1];
      indices[0] = ret;
    }
    i       = lens[3];
    lens[3] = lens[1];
    lens[1] = i;
    q.flip  = -1.0;
  }
  
  /* get the template & go */
  
  len =  lens[0]+lens[1]+lens[2]+lens[3];
  N   =  lens[0];
  M   =  lens[3];
  P   =  lens[1] - M;
  Q   = (lens[2] - N - P)/2;
  if ((N < 3) || (M < 3) || (P < 0) || (Q < 0)) {
    for (i = 0; i < len; i++) {
      j         = len+i+1;
      uv[2*i  ] = uvx[2*j  ];
      uv[2*i+1] = uvx[2*j+1];
    }
    EG_free(uvx);
    EG_free(indices);
    if (outLevel > 1)
      printf("%lX EGADS Info: Too small ->  %d %d (>3)   %d %d\n",
             tID, N, M, P, Q);
    return -2;
  }
  
  if (P == 0) {
    ret = EG_quadFillQ(face, &q, nsp, indices, lens, uvx, npts, uvs);
  } else if (Q == 0) {
    ret = EG_quadFillP(face, &q, nsp, indices, lens, uvx, npts, uvs);
  } else {
    ret = EG_quadFillG(face, &q, nsp, indices, lens, uvx, npts, uvs);
  }
  for (i = 0; i < len; i++) {
    j         = len+i+1;
    uv[2*i  ] = uvx[2*j  ];
    uv[2*i+1] = uvx[2*j+1];
  }
  EG_free(uvx);
  EG_free(indices);
  if (q.quads != NULL) EG_free(q.quads);
  if (ret != 0) {
    if (outLevel > 1)
      printf("%lX EGADS Info: quadFill = %d    N,M,P,Q = %d %d %d %d\n",
             tID, ret, N, M, P, Q);
    return ret;
  }
  if ((q.tris == NULL) || (*uvs == NULL)) {
    if (*uvs   != NULL) EG_free(*uvs);
    if (q.tris != NULL) EG_free(q.tris);
    return -99;
  }
  
  /* remap back to our UV */
  
  if ((unmap == 0) && (align == 0)) {
    ret = EG_dQuadTFI(elens, uv, *npts, *uvs);
    if (ret != 0) {
      EG_free(*uvs);
      EG_free(q.tris);
      *npts = 0;
      *uvs  = NULL;
      return ret;
    }
  }
  
  /* fix orientation if flipped direction */
  
  indices = q.tris;
  uvx     = *uvs;
  if (q.flip < 0.0)
    for (i = 0; i < q.ntris; i++) {
      len            = indices[3*i+2];
      indices[3*i+2] = indices[3*i+1];
      indices[3*i+1] = len;
    }
  
  /* make sure we are OK */
  
  for (k = i = 0; i < q.ntris; i++) {
    uv0 = &uvx[2*indices[3*i  ]];
    uv1 = &uvx[2*indices[3*i+1]];
    uv2 = &uvx[2*indices[3*i+2]];
    dist = AREA2D(uv0, uv1, uv2);
    if (dist*0.0 != 0.0) k++;		/* special Nan, Ind, ... checker */
    if (dist     <= 0.0) k++;
  }
  if (k != 0) {
    if (outLevel > 1)
      printf("%lX EGADS Info: Bad mapping - %d non-positive of %d quads\n",
             tID, k, q.nquad);
    EG_free(*uvs);
    EG_free(q.tris);
    *npts = 0;
    *uvs  = NULL;
    return -6;
  }
      
  *ntris = q.ntris;
  *tris  = q.tris;
  return ret;
}


static int
EG_quadfil3(long tID, const egObject *face, double *parms, int *elens,
            double *uv, int *npts, double **uvs, int *ntris, int **tris)
{
  int    i, j, k, len, ret, lens[4], degen, collap[2], *indices;
  double *uvx, *uv0, *uv1, *uv2, xylim[2][2], slim[3][2][2];
  double edgeTOL = 0.05;
  
  *ntris = 0;
  *uvs   = NULL;
  *tris  = NULL;
  if (elens[0] < 3) return -2;
  if (elens[1] < 3) return -2;
  if (elens[2] < 3) return -2;
  
  if ((parms[0] >= 0.001) && (parms[0] <= 0.5)) edgeTOL = parms[0];
  
  /* check side limits that have 2 min/maxs touching */
  
  if ((elens[0] < 2) || (elens[1] < 2) || (elens[2] < 2)) return -2;
  len = elens[0] + elens[1] + elens[2];
  xylim[0][0] = xylim[1][0] = uv[0];
  xylim[0][1] = xylim[1][1] = uv[1];
  for (i = 1; i < len; i++) {
    if (xylim[0][0] > uv[2*i  ]) xylim[0][0] = uv[2*i  ];
    if (xylim[1][0] < uv[2*i  ]) xylim[1][0] = uv[2*i  ];
    if (xylim[0][1] > uv[2*i+1]) xylim[0][1] = uv[2*i+1];
    if (xylim[1][1] < uv[2*i+1]) xylim[1][1] = uv[2*i+1];
  }
  slim[0][0][0] = slim[1][0][0] = slim[2][0][0] = xylim[1][0];
  slim[0][1][0] = slim[1][1][0] = slim[2][1][0] = xylim[0][0];
  slim[0][0][1] = slim[1][0][1] = slim[2][0][1] = xylim[1][1];
  slim[0][1][1] = slim[1][1][1] = slim[2][1][1] = xylim[0][1];
  for (i = 1; i < elens[0]; i++) {
    if (slim[0][0][0] > uv[2*i  ]) slim[0][0][0] = uv[2*i  ];
    if (slim[0][1][0] < uv[2*i  ]) slim[0][1][0] = uv[2*i  ];
    if (slim[0][0][1] > uv[2*i+1]) slim[0][0][1] = uv[2*i+1];
    if (slim[0][1][1] < uv[2*i+1]) slim[0][1][1] = uv[2*i+1];
  }
  for (i = elens[0]+1; i < elens[0]+elens[1]; i++) {
    if (slim[1][0][0] > uv[2*i  ]) slim[1][0][0] = uv[2*i  ];
    if (slim[1][1][0] < uv[2*i  ]) slim[1][1][0] = uv[2*i  ];
    if (slim[1][0][1] > uv[2*i+1]) slim[1][0][1] = uv[2*i+1];
    if (slim[1][1][1] < uv[2*i+1]) slim[1][1][1] = uv[2*i+1];
  }
  for (i = elens[0]+elens[1]+1; i < len; i++) {
    if (slim[2][0][0] > uv[2*i  ]) slim[2][0][0] = uv[2*i  ];
    if (slim[2][1][0] < uv[2*i  ]) slim[2][1][0] = uv[2*i  ];
    if (slim[2][0][1] > uv[2*i+1]) slim[2][0][1] = uv[2*i+1];
    if (slim[2][1][1] < uv[2*i+1]) slim[2][1][1] = uv[2*i+1];
  }
  for (i = 0; i < 3; i++) {
    lens[i] = 0;
    if ((slim[i][1][0]-slim[i][0][0])/(xylim[1][0]-xylim[0][0]) < edgeTOL)
      if ((slim[i][1][0]-xylim[0][0])/(xylim[1][0]-xylim[0][0]) < 0.5) {
        lens[i] |= 1;
      } else {
        lens[i] |= 2;
      }
    if ((slim[i][1][1]-slim[i][0][1])/(xylim[1][1]-xylim[0][1]) < edgeTOL)
      if ((slim[i][1][1]-xylim[0][1])/(xylim[1][1]-xylim[0][1]) < 0.5) {
        lens[i] |= 4;
      } else {
        lens[i] |= 8;
      }
  }
  for (i = 0; i < 3; i++)
    if ((lens[i] != 0) && (lens[i] != 1) && (lens[i] != 2) &&
        (lens[i] != 4) && (lens[i] != 8)) {
      printf(" %lX Odd Degenerate test for side %d: %d!\n", tID, i, lens[i]);
      return -2;
    }
#ifdef DEBUG
  printf(" side markers: %d %d    %d %d    %d %d\n", elens[0], lens[0],
         elens[1], lens[1], elens[2], lens[2]);
#endif
  
  /* do we need to test the end points? */
  degen = j = -1;
  if (((lens[0]==1) && (lens[1]==2)) || ((lens[0]==2) && (lens[1]==1)) ||
      ((lens[0]==4) && (lens[1]==8)) || ((lens[0]==8) && (lens[1]==4))) {
    ret   = 2;
    degen = 1;
    if ((lens[0]==4) || (lens[1]==4)) j = 1;
  }
  if (((lens[0]==1) && (lens[2]==2)) || ((lens[0]==2) && (lens[2]==1)) ||
      ((lens[0]==4) && (lens[2]==8)) || ((lens[0]==8) && (lens[2]==4))) {
    ret   = 1;
    degen = 0;
    if ((lens[0]==4) || (lens[2]==4)) j = 1;
  }
  if (((lens[1]==1) && (lens[2]==2)) || ((lens[1]==2) && (lens[2]==1)) ||
      ((lens[1]==4) && (lens[2]==8)) || ((lens[1]==8) && (lens[2]==4))) {
    ret   = 0;
    degen = 2;
    if ((lens[1]==4) || (lens[2]==4)) j = 1;
  }
  if (degen == -1) return -5;
  
  /* expand the degenerate point */
  
  len += elens[ret];
  uvx = (double *) EG_alloc(2*len*sizeof(double));
  if (uvx == NULL) return -1;
  
  if (degen == 0) {
    for (k = 0, i = elens[0]+elens[1]; i >= elens[0]; i--) {
      if (j == 1) {
        uvx[2*k  ] = uv[0];
        uvx[2*k+1] = uv[2*i+1];
      } else {
        uvx[2*k  ] = uv[2*i];
        uvx[2*k+1] = uv[1];
      }
      k++;
    }
    for (i = 1; i < elens[0]+elens[1]+elens[2]; i++) {
      uvx[2*k  ] = uv[2*i  ];
      uvx[2*k+1] = uv[2*i+1];
      k++;
    }
    collap[0] = 0;
    collap[1] = elens[1];
    lens[0]   = elens[1];
    lens[1]   = elens[0];
    lens[2]   = elens[1];
    lens[3]   = elens[2];
#ifdef DEBUG
    printf(" degen 0: lens = %d %d\n", k, lens[0]+lens[1]+lens[2]+lens[3]);
#endif
  } else if (degen == 1) {
    for (k = i = 0; i < elens[0]; i++) {
      uvx[2*k  ] = uv[2*i  ];
      uvx[2*k+1] = uv[2*i+1];
      k++;
    }
    if (j == 1) {
      uvx[2*k  ] = uv[2*elens[0]  ];
      uvx[2*k+1] = uv[1];
    } else {
      uvx[2*k  ] = uv[0];
      uvx[2*k+1] = uv[2*elens[0]+1];
    }
    k++;
    for (i = elens[0]+elens[1]+elens[2]-1; i >= elens[0]+elens[1]; i--) {
      if (j == 1) {
        uvx[2*k  ] = uv[2*elens[0]  ];
        uvx[2*k+1] = uv[2*i+1];
      } else {
        uvx[2*k  ] = uv[2*i];
        uvx[2*k+1] = uv[2*elens[0]+1];
      }
      k++;
    }
    for (i = elens[0]+1; i < elens[0]+elens[1]+elens[2]; i++) {
      uvx[2*k  ] = uv[2*i  ];
      uvx[2*k+1] = uv[2*i+1];
      k++;
    }
    collap[0] = elens[0];
    collap[1] = elens[0]+elens[2];
    lens[0]   = elens[0];
    lens[1]   = elens[2];
    lens[2]   = elens[1];
    lens[3]   = elens[2];
#ifdef DEBUG
    printf(" degen 1: lens = %d %d\n", k, lens[0]+lens[1]+lens[2]+lens[3]);
#endif
  } else {
    for (k = i = 0; i < elens[0]+elens[1]; i++) {
      uvx[2*k  ] = uv[2*i  ];
      uvx[2*k+1] = uv[2*i+1];
      k++;
    }
    for (i = elens[0]; i >= 0; i--) {
      if (j == 1) {
        uvx[2*k  ] = uv[2*(elens[0]+elens[1])  ];
        uvx[2*k+1] = uv[2*i+1];
      } else {
        uvx[2*k  ] = uv[2*i];
        uvx[2*k+1] = uv[2*(elens[0]+elens[1])+1];
      }
      k++;
    }
    for (i = elens[0]+elens[1]+1; i < elens[0]+elens[1]+elens[2]; i++) {
      uvx[2*k  ] = uv[2*i  ];
      uvx[2*k+1] = uv[2*i+1];
      k++;
    }
    collap[0] = elens[0]+elens[1];
    collap[1] = elens[0]+elens[1]+elens[0];
    lens[0]   = elens[0];
    lens[1]   = elens[1];
    lens[2]   = elens[0];
    lens[3]   = elens[2];
#ifdef DEBUG
    printf(" degen 2: lens = %d %d\n", k, lens[0]+lens[1]+lens[2]+lens[3]);
#endif
  }
  
  /* fill the area */
  
  ret = EG_quad2tris(tID, face, parms, lens, uvx, npts, uvs, ntris, tris);
  EG_free(uvx);
  if (ret != 0) return ret;
  
  /* re-collapse the degenerate side */
  
  uvx = *uvs;
  if (uvx == NULL) return -99;
  for (k = i = 0; i < *npts; i++)
    if ((i <= collap[0]) || (i > collap[1])) {
      uvx[2*k  ] = uvx[2*i  ];
      uvx[2*k+1] = uvx[2*i+1];
      k++;
    }
  len = collap[0];
  uvx[2*len  ] = uv[2*len  ];
  uvx[2*len+1] = uv[2*len+1];
  *npts = k;
  
  indices = *tris;
  if (indices == NULL) return -99;
  for (k = i = 0; i < *ntris; i++) {
    if ((indices[3*i  ] > collap[0]) && (indices[3*i  ] <= collap[1])) {
      indices[3*i  ] = collap[0];
    } else if (indices[3*i  ] > collap[1]) {
      indices[3*i  ] = indices[3*i  ] - collap[1] + collap[0];
    }
    if ((indices[3*i+1] > collap[0]) && (indices[3*i+1] <= collap[1])) {
      indices[3*i+1] = collap[0];
    } else if (indices[3*i+1] > collap[1]) {
      indices[3*i+1] = indices[3*i+1] - collap[1] + collap[0];
    }
    if ((indices[3*i+2] > collap[0]) && (indices[3*i+2] <= collap[1])) {
      indices[3*i+2] = collap[0];
    } else if (indices[3*i+2] > collap[1]) {
      indices[3*i+2] = indices[3*i+2] - collap[1] + collap[0];
    }
    if ((indices[3*i  ] != indices[3*i+1]) &&
        (indices[3*i  ] != indices[3*i+2]) &&
        (indices[3*i+1] != indices[3*i+2])) {
      indices[3*k  ] = indices[3*i  ];
      indices[3*k+1] = indices[3*i+1];
      indices[3*k+2] = indices[3*i+2];
      k++;
    }
    
  }
  *ntris = k;
  
  /* make sure we are OK */
  
  for (k = i = 0; i < *ntris; i++) {
    uv0 = &uvx[2*indices[3*i  ]];
    uv1 = &uvx[2*indices[3*i+1]];
    uv2 = &uvx[2*indices[3*i+2]];
    if (AREA2D(uv0, uv1, uv2) <= 0.0) k++;
  }
  if (k == 0) return 0;
  
#ifdef DEBUG
  printf(" Bad mapping3: %d of %d non-positive tris\n", k, *ntris);
#endif
  EG_free(*uvs);
  EG_free(*tris);
  *npts = *ntris = 0;
  *uvs  = NULL;
  *tris = NULL;
  return -8;
}


static int
EG_quadfil3as4(long tID, const egObject *face, double *parms, int *elens,
               double *uv, int *npts, double **uvs, int *ntris, int **tris)
{
  int    i, k, ret, lens[4], *indices;
  double *uvx, *uv0, *uv1, *uv2;
  
#ifdef DEBUG
  printf("EG_quadfil3as4: Filling 3-sided Face as if it were 4-sided\n");
  printf("EG_quadfil3as4: %d %d %d\n", elens[0], elens[1], elens[2]);
#endif
  
  *ntris = 0;
  *uvs  = NULL;
  *tris = NULL;
  if (elens[0] < 3) return -2;
  if (elens[1] < 3) return -2;
  if (elens[2] < 3) return -2;
  
  if (abs(elens[1]-elens[2]) < 2) {	/* Similar edge balance */
    ret = 0;
  } else if (abs(elens[0]-elens[2]) < 2) {
    ret = 1;
  } else if (abs(elens[0]-elens[1]) < 2) {
    ret = 2;
    
    /* Split side with max points */
    
  } else if ((elens[0] >= elens[1]) && (elens[0] >= elens[2])) {
    ret = 0;
  } else if  (elens[1] >= elens[2]) {
    ret = 1;
  } else {
    ret = 2;
  }
  
  if (ret == 0) {
    lens[0] = elens[0]/2;
    lens[1] = elens[0]-lens[0];
    lens[2] = elens[1];
    lens[3] = elens[2];
  } else if (ret == 1) {
    lens[0] = elens[0];
    lens[1] = elens[1]/2;
    lens[2] = elens[1]-lens[1];
    lens[3] = elens[2];
  } else if (ret == 2) {
    lens[0] = elens[0];
    lens[1] = elens[1];
    lens[2] = elens[2]/2;
    lens[3] = elens[2]-lens[2];
  }
  
#ifdef DEBUG
  printf("EG_quadfil3as4: Quad with %d %d %d %d\n",
         lens[0], lens[1], lens[2], lens[3]);
#endif
  
  /* fill the area */
  
  ret = EG_quad2tris(tID, face, parms, lens, uv, npts, uvs, ntris, tris);
  if (ret != 0) return ret;
  
  uvx     = *uvs;
  indices = *tris;
  if ((indices == NULL) || (uvx == NULL)) return -99;
  for (k = i = 0; i < *ntris; i++) {
    uv0 = &uvx[2*indices[3*i  ]];
    uv1 = &uvx[2*indices[3*i+1]];
    uv2 = &uvx[2*indices[3*i+2]];
    if (AREA2D(uv0, uv1, uv2) <= 0.0) k++;
  }
  if (k == 0) return 0;
  
#ifdef DEBUG
  printf(" Bad mapping3: %d of %d non-positive tris\n", k, *ntris);
#endif
  EG_free(*uvs);
  EG_free(*tris);
  *npts = *ntris = 0;
  *uvs  = NULL;
  *tris = NULL;
  return -8;
}


/* takes a simple 3 sided loop determines if a node is degenerate and if so
 *       fills it with tris based on a quad sub-blocking scheme that
 *       supports differing sizes per side
 *
 * where: parms[0] - Edge Tol
 *        parms[1] - Side Ratio
 *        parms[2] - # smoothing passes
 *
 *        elens[0] - number of segments on the left side
 *        elens[1] - number of segments on the bottom
 *        elens[2] - number of segments on the right
 *
 *        uv[]     - input (u,v) pairs going around the loop CCW with
 *                   no duplicates at corners starting at UL corner (the
 *                   start of side 0)
 *                   len = 2*(sum of elens)
 *
 *        npts     - total number of resultant points (returned)
 *        uvs      - pointer to (u,v) pairs of all of the points (returned,
 *                   should be free'd);  note -> the first set matches uv
 *        ntris    - resultant number of tris
 *        tris     - pointer to index triads of uvs that support the
 *                   split quad patch(s) (returned, should be free'd)
 *
 * return codes:  0 - success
 *               -1 - malloc error
 *               -2 - elen error
 *               -3 - block side too big
 *               -4 - extra edge not found
 *               -5 - not degenerate
 *               -6 - neg area tris
 */

int
EG_quad2tris3(long tID, const egObject *face, double *parms, int *elens,
              double *uv, int *npts, double **uvs, int *ntris, int **tris,
              int *flag)
{
  int    i, j, k, n, ret, nsp, *indices;
  double *uvx, *uv0, *uv1, *uv2, dist;
  qFill  q;
  
  *npts = 0;
  *flag = 1;
  /* try and make a conical mapping */
  ret   = EG_quadfil3(tID, face, parms, elens, uv, npts, uvs, ntris, tris);
  if (ret == -8) {
    *npts = 1;
    ret   = EG_quadfil3(tID, face, parms, elens, uv, npts, uvs, ntris, tris);
    if (ret == -8) ret = -6;
  }
  
  /* try the triangle template */
  if (ret != 0) {
    *flag      = 0;
    *npts      = 0;
    n          = elens[0]+elens[1]+elens[2];
    indices    = (int *) EG_alloc((n+1)*sizeof(int));
    if (indices == NULL) return -1;
    for (i = 0; i < n; i++) indices[i] = i;
    indices[n] = 0;
    q.outLevel = EG_outLevel(face);
    q.flip     = 1.0;
    q.fillT    = 1;
    q.ntris    = 0;
    q.quads    = NULL;
    nsp        = parms[2]+0.1;
    ret = EG_triTemplate(tID, face, &q, nsp, indices, elens, uv, npts, uvs);
    /* make tris & free up q */
    EG_free(indices);
    if (q.quads != NULL) EG_free(q.quads);
    if (ret != 0) {
      if (q.outLevel > 1)
        printf("%lX EGADS Info: EG_triTemplate = %d!\n", tID, ret);
    } else {
      if ((q.tris == NULL) || (*uvs == NULL)) {
        if (*uvs   != NULL) EG_free(*uvs);
        if (q.tris != NULL) EG_free(q.tris);
        ret = -1;
      } else {
        /* fix orientation if flipped direction */
        indices = q.tris;
        uvx     = *uvs;
        if (q.flip < 0.0)
          for (i = 0; i < q.ntris; i++) {
            j              = indices[3*i+2];
            indices[3*i+2] = indices[3*i+1];
            indices[3*i+1] = j;
          }
        /* make sure we are OK */
        for (k = i = 0; i < q.ntris; i++) {
          uv0 = &uvx[2*indices[3*i  ]];
          uv1 = &uvx[2*indices[3*i+1]];
          uv2 = &uvx[2*indices[3*i+2]];
          dist = AREA2D(uv0, uv1, uv2);
          if (dist*0.0 != 0.0) k++;         /* special Nan, Ind, ... checker */
          if (dist     <= 0.0) k++;
        }
        if (k != 0) {
          if (q.outLevel > 1)
            printf("%lX EGADS Info: Bad mapping - %d non-positive of %d quads\n",
                   tID, k, q.nquad);
          EG_free(*uvs);
          EG_free(q.tris);
          *npts = 0;
          *uvs  = NULL;
          ret   = -6;
        } else {
          *ntris = q.ntris;
          *tris  = q.tris;
        }
      }
    }
  }
  
  /* finally see if splitting an Edge gives us a good mapping */
  if (ret != 0) {
    *npts = 0;
    ret   = EG_quadfil3as4(tID, face, parms, elens, uv, npts, uvs, ntris, tris);
  }
  
  return ret;
}
