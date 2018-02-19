/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Tessellation Functions
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>	/* Needed in some systems for DBL_MAX definition */
#include <limits.h>	/* Needed in other systems for DBL_EPSILON */

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsTris.h"
#include "emp.h"


#define INSERTKNOTS

#define NOTFILLED	-1
#define TOL		 1.e-7
#define UVTOL            1.e-4


#define AREA2D(a,b,c)   ((a[0]-c[0])*(b[1]-c[1]) - (a[1]-c[1])*(b[0]-c[0]))
#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define DIST2(a,b)      ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
#define DOT2(a,b)       (((a)[0]*(b)[0]) + ((a)[1]*(b)[1]))
#define VSUB2(a,b,c)    (c)[0] = (a)[0] - (b)[0]; (c)[1] = (a)[1] - (b)[1];
#define MAX(a,b)        (((a) > (b)) ?  (a) : (b))
#define MIN(a,b)        (((a) < (b)) ?  (a) : (b))
#define ABS(a)          (((a) <   0) ? -(a) : (a))


  extern int  EG_getTolerance( const egObject *topo, double *tol );
  extern int  EG_getTopology( const egObject *topo, egObject **geom, int *oclas,
                              int *type, /*@null@*/ double *limits, int *nChild,
                              egObject ***children, int **senses );
  extern int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                               int oclass, int *ntopo, egObject ***topos );
  extern int  EG_indexBodyTopo( const egObject *body, const egObject *src );
  extern int  EG_sameBodyTopo( const egObject *bod1, const egObject *bod2 );
  extern int  EG_getEdgeUV( const egObject *face, const egObject *edge,
                            int sense, double t, double *result );
  extern int  EG_getRange( const egObject *geom, double *range, int *pflag );
  extern int  EG_getGeometry( const egObject *geom, int *oclass, int *type,
                              egObject **rGeom, int **ivec, double **rvec );
  extern int  EG_evaluate( const egObject *geom, const double *param,
                           double *result );
  extern int  EG_invEvaluate( const egObject *geom, double *xyz, double *param,
                              double *result );
  extern int  EG_curvature( const egObject *geom, const double *param,
                            double *result );
  extern int  EG_isSame( const egObject *geom1, const egObject *geom2 );
  extern int  EG_attributeRet( const ego obj, const char *name, int *atype,
                               int *len, /*@null@*/ const int    **ints,
                                         /*@null@*/ const double **reals,
                                         /*@null@*/ const char   **str );

  extern int  EG_tessellate( int outLevel, triStruct *ts, long tID );
  extern int  EG_quadFill( const egObject *face, double *parms, int *elens,
                           double *uv, int *npts, double **uvs, int *npat,
                           int *pats, int **vpats );
  extern int  EG_baryFrame( egTess2D *tess2d );
  extern int  EG_baryTess( egTess2D tess2d, const double *uv, double *w );
  extern void EG_mapTessTs( egTess1D src, egTess1D dst );
  extern int  EG_relPosTs( egObject *geom, int n, const double *rel, double *ts,
                            double *xyzs );
#ifdef INSERTKNOTS
  extern int  EG_mapSequen(egObject *src, egObject *dst, egObject **result);
#endif



static int
EG_faceConnIndex(egFconn conn, int face)
{
  int i;

  if (conn.nface == 1) {
    if (conn.index == face) return 1;
  } else {
    for (i = 0; i < conn.nface; i++)
      if (conn.faces[i] == face) return i+1;
  }
    
  return 0;
}


#ifdef CHECK
static void
EG_checkTriangulation(egTessel *btess)
{
  int i, j, k, n, n1, n2, nf, iface, itri, ie, iv, side;
  static int sides[3][2] = {1,2, 2,0, 0,1};

  for (iface = 1; iface <= btess->nFace; iface++) {
    for (itri = 1; itri <= btess->tess2d[iface-1].ntris; itri++) {
      for (j = 0; j < 3; j++) {
        if ((btess->tess2d[iface-1].tris[3*itri+j-3] >
             btess->tess2d[iface-1].npts) ||
            (btess->tess2d[iface-1].tris[3*itri+j-3] <= 0))
          printf(" checkTriangulation: Face %d, Tri %d[%d] = %d!\n",
                 iface, itri, j, btess->tess2d[iface-1].tris[3*itri+j-3]);
        n = btess->tess2d[iface-1].tric[3*itri+j-3];
        if (n > btess->tess2d[iface-1].ntris) {
          printf(" checkTriangulation: Face %d, Nei %d[%d] = %d (%d)!\n",
                 iface, itri, j, n, btess->tess2d[iface-1].ntris);
        } else if (n == 0) {
          printf(" checkTriangulation: Face %d, No Neighbor %d[%d]\n",
                 iface, itri, j);          
        } else if (n > 0) {
          side = -1;
          if (btess->tess2d[iface-1].tric[3*n-3] == itri) side = 0;
          if (btess->tess2d[iface-1].tric[3*n-2] == itri) side = 1;
          if (btess->tess2d[iface-1].tric[3*n-1] == itri) side = 2;
          if (side == -1) {
            printf(" checkTriangulation: Face %d, Tri Nei %d[%d] = %d!\n",
                   iface, itri, j, n);
            printf("                             Tri Nei %d[0] = %d\n",
                    n, btess->tess2d[iface-1].tric[3*n-3]);
            printf("                             Tri Nei %d[1] = %d\n",
                    n, btess->tess2d[iface-1].tric[3*n-2]);
            printf("                             Tri Nei %d[2] = %d\n",
                    n, btess->tess2d[iface-1].tric[3*n-1]);
          } else {
            n1 = btess->tess2d[iface-1].tris[3*itri+sides[j][0]-3];
            n2 = btess->tess2d[iface-1].tris[3*itri+sides[j][1]-3];
            if (((n1 != btess->tess2d[iface-1].tris[3*n+sides[side][0]-3]) ||
                 (n2 != btess->tess2d[iface-1].tris[3*n+sides[side][1]-3])) &&
                ((n1 != btess->tess2d[iface-1].tris[3*n+sides[side][1]-3]) ||
                 (n2 != btess->tess2d[iface-1].tris[3*n+sides[side][0]-3]))) {
              printf(" checkTriangulation: Face %d, Tri Nei %d[%d] = %d!\n",
                     iface, itri, j, n);
              printf("                             verts = %d %d, %d %d\n", 
                     n1, n2, btess->tess2d[iface-1].tris[3*n+sides[side][0]-3],
                             btess->tess2d[iface-1].tris[3*n+sides[side][1]-3]);
            }
          }
        } else {
          n1 = btess->tess2d[iface-1].tris[3*itri+sides[j][0]-3];
          n2 = btess->tess2d[iface-1].tris[3*itri+sides[j][1]-3];
          ie = -n;
          iv =  0;
          if (btess->tess2d[iface-1].ptype[n1-1] == -1) {
            printf(" checkTriangulation: Face %d, Tri Nei1 %d[%d] Interior Vert!\n",
                   iface, itri, j);            
          } else if (btess->tess2d[iface-1].ptype[n1-1] > 0) {
            if (btess->tess2d[iface-1].pindex[n1-1] != ie) {
              printf(" checkTriangulation: Face %d, Tri Nei1 %d[%d] Edge %d %d!\n",
                     iface, itri, j, ie, btess->tess2d[iface-1].pindex[n1-1]);
            } else {
              iv = btess->tess2d[iface-1].ptype[n1-1];
            }
          }
          if (btess->tess2d[iface-1].ptype[n2-1] == -1) {
            printf(" checkTriangulation: Face %d, Tri Nei2 %d[%d] Interior Vert!\n",
                   iface, itri, j);
            iv = 0;
          } else if (btess->tess2d[iface-1].ptype[n2-1] > 0) {
            if (btess->tess2d[iface-1].pindex[n2-1] != ie) {
              printf(" checkTriangulation: Face %d, Tri Nei2 %d[%d] Edge %d %d!\n",
                     iface, itri, j, ie, btess->tess2d[iface-1].pindex[n2-1]);
              iv = 0;
            } else {
              if ((iv != 0) && (iv > btess->tess2d[iface-1].ptype[n2-1]))
                iv = btess->tess2d[iface-1].ptype[n2-1];
            }
          } else {
            iv = 0;
          }
          if ((ie < 1) || (ie > btess->nEdge)) {
            printf(" checkTriangulation: Face %d, Tri Nei %d[%d] = %d (%d)!\n",
                   iface, itri, j, ie, btess->nEdge);
          } else {
            if (iv == 0) {
              for (i = 0; i <  btess->tess1d[ie-1].npts-1; 
                          i += btess->tess1d[ie-1].npts-2) {
                nf = btess->tess1d[ie-1].faces[0].nface;
                if (nf > 0) {
                  k = EG_faceConnIndex(btess->tess1d[ie-1].faces[0], iface);
                  if (k != 0)
                    if (btess->tess1d[ie-1].faces[0].tric[i*nf+k-1] == itri) break;
                }
                nf = btess->tess1d[ie-1].faces[1].nface;
                if (nf > 0) {
                  k = EG_faceConnIndex(btess->tess1d[ie-1].faces[1], iface);
                  if (k != 0) 
                    if (btess->tess1d[ie-1].faces[1].tric[i*nf+k-1] == itri) break;
                }
              }
              if (i > btess->tess1d[ie-1].npts-1)
                printf(" checkTriangulation: Face %d, Tri Nei %d[%d] Not Found in %d!\n",
                       iface, itri, j, ie);
            } else {
              i  = 0;
              nf = btess->tess1d[ie-1].faces[0].nface;
              if (nf > 0) {
                k = EG_faceConnIndex(btess->tess1d[ie-1].faces[0], iface);
                if (k != 0)
                  if (btess->tess1d[ie-1].faces[0].tric[(iv-1)*nf+k-1] == itri) i++;
              }
              nf = btess->tess1d[ie-1].faces[1].nface;
              if (nf > 0) {
                k = EG_faceConnIndex(btess->tess1d[ie-1].faces[1], iface);
                if (k != 0)
                  if (btess->tess1d[ie-1].faces[1].tric[(iv-1)*nf+k-1] == itri) i++;
              }
              if (i == 0) {
                printf(" checkTriangulation: Face %d, Tri Nei %d[%d] Edge %d =",
                       iface, itri, j, ie);
                nf = btess->tess1d[ie-1].faces[0].nface;
                k  = EG_faceConnIndex(btess->tess1d[ie-1].faces[0], iface);
                if (k != 0)
                  printf(" %d", btess->tess1d[ie-1].faces[0].tric[(iv-1)*nf+k-1]);
                nf = btess->tess1d[ie-1].faces[1].nface;
                k  = EG_faceConnIndex(btess->tess1d[ie-1].faces[1], iface);
                if (k != 0)                
                  printf(" %d", btess->tess1d[ie-1].faces[1].tric[(iv-1)*nf+k-1]);
                printf("!\n");
              }
            }
          }
        }
      }
    }
  }
}
#endif


/* 
 * determine if this line segment crosses any active segments 
 * pass:      0 - first pass; conservative algorithm
 * 	      1 - second pass; use dirty tricks
 */

static int
EG_crossSeg(int index, const double *mid, int i2, const double *vertices, 
            int pass, fillArea *fa)
{
  int    i, i0, i1, iF0, iF1;
  double angle, cosan, sinan, dist2, distF, eps, ty0, ty1, frac;
  double uv0[2], uv1[2], uv2[2], x[2];
  double uvF0[2], uvF1[2];

  uv2[0] = vertices[2*i2  ];
  uv2[1] = vertices[2*i2+1];

  /*  Store away coordinates of front */
  iF0 	  = fa->front[index].i0;
  iF1 	  = fa->front[index].i1;
  uvF0[0] = vertices[2*iF0  ];
  uvF0[1] = vertices[2*iF0+1];
  uvF1[0] = vertices[2*iF1  ];
  uvF1[1] = vertices[2*iF1+1];

  dist2  = DIST2( mid,  uv2);
  distF  = DIST2(uvF0, uvF1);
  eps    = (dist2 + distF) * DBL_EPSILON;

  /* transform so that we are in mid-uv2 coordinate frame */
  angle = atan2(uv2[1]-mid[1], uv2[0]-mid[0]);
  cosan = cos(angle);
  sinan = sin(angle);

  /* look at the current front */

  for (i = 0; i < fa->nfront; i++) {
    if ((i == index) || (fa->front[i].sright == NOTFILLED)) continue;
    if (fa->front[i].snew == 0) continue;
    i0 = fa->front[i].i0;
    i1 = fa->front[i].i1;
    if ((i0 == i2) || (i1 == i2)) continue;
    uv0[0] = vertices[2*i0  ];
    uv0[1] = vertices[2*i0+1];
    uv1[0] = vertices[2*i1  ];
    uv1[1] = vertices[2*i1+1];

    /* look to see if the transformed y's from uv2-mid cross 0.0 */
    ty0 = (uv0[1]-mid[1])*cosan - (uv0[0]-mid[0])*sinan;
    ty1 = (uv1[1]-mid[1])*cosan - (uv1[0]-mid[0])*sinan;
    if ((ty0 == 0.0) && (ty1 == 0.0)) return 1;
    if  (ty0*ty1 >= 0.0) continue;

    /* get fraction of line for crossing */
    frac = -ty0/(ty1-ty0);
    if ((frac < 0.0) || (frac > 1.0)) continue;

    /* get the actual coordinates */
    x[0] = uv0[0] + frac*(uv1[0]-uv0[0]);
    x[1] = uv0[1] + frac*(uv1[1]-uv0[1]);

    /* are we in the range for the line seg? */
    frac = (x[0]-mid[0])*cosan + (x[1]-mid[1])*sinan;
    if ((frac > 0.0) && (frac*frac < dist2*(1.0+TOL))) return 2;
  }

  /* look at our original loops */

  for (i = 0; i < fa->nsegs; i++) {
    double area10, area01, area11, area00;

    i0 = fa->segs[2*i  ];
    i1 = fa->segs[2*i+1];

    if (((i0 == fa->front[index].i0) && (i1 == fa->front[index].i1)) ||
        ((i0 == fa->front[index].i1) && (i1 == fa->front[index].i0))) 
      continue;

    uv0[0] = vertices[2*i0  ];
    uv0[1] = vertices[2*i0+1];
    uv1[0] = vertices[2*i1  ];
    uv1[1] = vertices[2*i1+1];

    if (pass != 0) {
      area10 = AREA2D(uv2, uv1, uvF0);
      area00 = AREA2D(uv2, uv0, uvF0);
      area10 = ABS(area10);
      area00 = ABS(area00);
    }
    if ((pass != 0) && area10 < eps && area00 < eps) {
      /* I2 and Boundary Segment are collinear with IF0 (Front.I0) */
      double del0[2], del1[2], del2[2];

      VSUB2(uv2, uvF0, del2);
      VSUB2(uv1, uvF0, del1);
      VSUB2(uv0, uvF0, del0);
      /*  See if I1 is between IF0 and I2 */
      if (i1 != iF0 && DOT2(del2, del1) > 0 &&
	  DOT2(del2, del2) > DOT2(del1, del1)) return 5;
      /*  See if I0 is between IF0 and I2 */
      if (i0 != iF0 && DOT2(del2, del0) > 0 && 
	  DOT2(del2, del2) > DOT2(del0, del0)) return 6;
    }
    if (pass != 0) {
      area11 = AREA2D(uv2, uv1, uvF1);
      area01 = AREA2D(uv2, uv0, uvF1);
      area11 = ABS(area11);
      area01 = ABS(area01);
    }
    if ((pass != 0) && area11 < eps && area01 < eps) {
      /* I2 and Boundary Segment are collinear with IF1 (Front.I1) */
      double del0[2], del1[2], del2[2];

      VSUB2(uv2, uvF1, del2);
      VSUB2(uv1, uvF1, del1);
      VSUB2(uv0, uvF1, del0);
      /*  See if I1 is between IF1 and I2 */
      if (i1 != iF1 && DOT2(del2, del1) > 0 &&
	  DOT2(del2, del2) > DOT2(del1, del1)) return 7;
      /*  See if I0 is between IF1 and I2 */
      if (i0 != iF1 && DOT2(del2, del0) > 0 && 
	  DOT2(del2, del2) > DOT2(del0, del0)) return 8;
    }

    if ((i1 == i2) || (i0 == i2)) continue;
    /* look to see if the transformed y's from uv2-mid cross 0.0 */
    ty0 = (uv0[1]-mid[1])*cosan - (uv0[0]-mid[0])*sinan;
    ty1 = (uv1[1]-mid[1])*cosan - (uv1[0]-mid[0])*sinan;
    if ((ty0 == 0.0) && (ty1 == 0.0)) return 3;
    if  (ty0*ty1 >= 0.0) continue;

    /* get fraction of line for crossing */
    frac = -ty0/(ty1-ty0);
    if ((frac < 0.0) || (frac > 1.0)) continue;

    /* get the actual coordinates */
    x[0] = uv0[0] + frac*(uv1[0]-uv0[0]);
    x[1] = uv0[1] + frac*(uv1[1]-uv0[1]);

    /* are we in the range for the line seg? */
    frac = (x[0]-mid[0])*cosan + (x[1]-mid[1])*sinan;
    if ((frac > 0.0) && (frac*frac < dist2*(1.0+TOL))) return 4;
  }

  return 0;
}


/* Input specified as contours.
 * Outer contour must be counterclockwise.
 * All inner contours must be clockwise.
 *  
 * Every contour is specified by giving all its points in order. No
 * point shoud be repeated. i.e. if the outer contour is a square,
 * only the four distinct endpoints should be specified in order.
 *  
 * ncontours: #contours
 * cntr: An array describing the number of points in each
 *	 contour. Thus, cntr[i] = #points in the i'th contour.
 * vertices: Input array of vertices. Vertices for each contour
 *           immediately follow those for previous one. Array location
 *           vertices[0] must NOT be used (i.e. i/p starts from
 *           vertices[1] instead. The output triangles are
 *	     specified  WRT the indices of these vertices.
 * triangles: Output array to hold triangles (allocated before the call)
 * pass:      0 - first pass; conservative algorithm
 * 	      1 - second pass; use dirty tricks
 *  
 * The number of output triangles produced for a polygon with n points is:
 *    (n - 2) + 2*(#holes)
 *
 * returns: -1 degenerate contour (zero length segment)
 *           0 allocation error
 *           + number of triangles
 */
int
EG_fillArea(int ncontours, const int *cntr, const double *vertices,
            int *triangles, int *n_fig8, int pass, fillArea *fa)
{
  int    i, j, i0, i1, i2, index, indx2, k, l, npts, neg;
  int    start, next, left, right, ntri, mtri;
  double side2, dist, d, area, uv0[2], uv1[2], uv2[2], mid[2];
  Front  *tmp;
  int    *itmp;

  *n_fig8 = 0;
  for (i = 0; i < ncontours; i++) if (cntr[i] < 3) return -1;
  for (fa->nfront = i = 0; i < ncontours; i++) fa->nfront += cntr[i];
  if (fa->nfront == 0) return -1;
  fa->npts = fa->nsegs = fa->nfront;

  mtri = fa->nfront - 2 + 2*(ncontours-1);
  ntri = 0;

  /* allocate the memory for the front */

  if (fa->front == NULL) {
    fa->mfront = CHUNK;
    while (fa->mfront < fa->nfront) fa->mfront += CHUNK;
    fa->front = (Front *) EG_alloc(fa->mfront*sizeof(Front));
    if (fa->front == NULL) return 0;
    fa->segs  = (int *) EG_alloc(2*fa->mfront*sizeof(int));
    if (fa->segs  == NULL) {
      EG_free(fa->front);
      fa->front = NULL;
      return 0;
    }
  } else {
    if (fa->mfront < fa->nfront) {
      i = fa->mfront;
      while (i < fa->nfront) i += CHUNK;
      tmp = (Front *) EG_reall(fa->front, i*sizeof(Front));
      if (tmp == NULL) return 0;
      itmp = (int *) EG_reall(fa->segs, 2*i*sizeof(int));
      if (itmp == NULL) {
        EG_free(tmp);
        return 0;
      }
      fa->mfront = i;
      fa->front  =  tmp;
      fa->segs   = itmp;
    }
  }

  /* allocate the memory for our point markers */
  npts = fa->nfront+1;
  if (fa->pts == NULL) {
    fa->mpts = CHUNK;
    while (fa->mpts < npts) fa->mpts += CHUNK;
    fa->pts = (int *) EG_alloc(fa->mpts*sizeof(int));
    if (fa->pts == NULL) return 0;
  } else {
    if (fa->mpts < npts) {
      i = fa->mpts;
      while (i < npts) i += CHUNK;
      itmp = (int *) EG_reall(fa->pts, i*sizeof(int));
      if (itmp == NULL) return 0;
      fa->mpts = i;
      fa->pts  = itmp;
    }
  }

  /* initialize the front */
  for (start = index = i = 0; i < ncontours; i++) {
    left = start + cntr[i] - 1;
    for (j = 0; j < cntr[i]; j++, index++) {
      fa->segs[2*index  ]     = left      + 1;
      fa->segs[2*index+1]     = start + j + 1;
      fa->front[index].sleft  = left;
      fa->front[index].i0     = left      + 1;
      fa->front[index].i1     = start + j + 1;
      fa->front[index].sright = start + j + 1;
      fa->front[index].snew   = 0;
      left = start + j;
    }
    fa->front[index-1].sright = start;

    /* look for fig 8 nodes in the contour */
    for (j = 0; j < cntr[i]-1; j++) {
      i0 = start + j + 1;
      for (k = j+1; k < cntr[i]; k++) {
        i1 = start + k + 1;
        if ((vertices[2*i0  ] == vertices[2*i1  ]) &&
            (vertices[2*i0+1] == vertices[2*i1+1])) {
          if (i0+1 == i1) {
            printf(" EGADS Internal: Null in loop %d -> %d %d\n", i, i0, i1);
            continue;
          }
          printf(" EGADS Internal: Fig 8 in loop %d (%d) -> %d %d (removed)\n",
                 i, ncontours, i0, i1);
	  /* figure 8's in the external loop decrease the triangle count */
          if (i == 0) (*n_fig8)++;  /* . . . . sometimes                 */
          for (l = 0; l < index; l++) {
            if (fa->front[l].i0 == i1) fa->front[l].i0 = i0;
            if (fa->front[l].i1 == i1) fa->front[l].i1 = i0;
          }
        }
      }
    }
    start += cntr[i];
  }

  /* collapse the front while building the triangle list*/

  neg = 0;
  do {

    /* count the number of vertex hits (right-hand links) */

    for (i = 0; i < npts; i++) fa->pts[i] = 0;
    for (i = 0; i < fa->nfront; i++) 
      if (fa->front[i].sright != NOTFILLED) fa->pts[fa->front[i].i1]++;

    /* remove any simple isolated triangles */

    for (j = i = 0; i < fa->nfront; i++) {
      if (fa->front[i].sright == NOTFILLED) continue;
      i0    = fa->front[i].i0;
      i1    = fa->front[i].i1;
      right = fa->front[i].sright;
      left  = fa->front[right].sright;
      if (fa->front[left].i1 == i0) {
        i2     = fa->front[right].i1;
        uv0[0] = vertices[2*i0  ];
        uv0[1] = vertices[2*i0+1];
        uv1[0] = vertices[2*i1  ];
        uv1[1] = vertices[2*i1+1];
        uv2[0] = vertices[2*i2  ];
        uv2[1] = vertices[2*i2+1];
        area   = AREA2D(uv0, uv1, uv2);
        if ((neg == 0) && (area <= 0.0)) continue;
        if (fa->front[left].sright != i) {
          start = fa->front[left].sright;
          fa->front[start].sleft  = fa->front[i].sleft;
          start = fa->front[i].sleft;
          fa->front[start].sright = fa->front[left].sright;
        }
        triangles[3*ntri  ]    = i0;
        triangles[3*ntri+1]    = i1;
        triangles[3*ntri+2]    = i2;
        fa->front[i].sleft     = fa->front[i].sright     = NOTFILLED;
        fa->front[right].sleft = fa->front[right].sright = NOTFILLED;
        fa->front[left].sleft  = fa->front[left].sright  = NOTFILLED;
        ntri++;
        j++;
        if (ntri >= mtri) break;
        neg = 0;
      }
    }
    if (j != 0) continue;

    /* look for triangles hidden by "figure 8" vetrices */

    for (j = i = 0; i < fa->nfront; i++) {
      if (fa->front[i].sright == NOTFILLED) continue;
      i0 = fa->front[i].i0;
      i1 = fa->front[i].i1;
      if (fa->pts[i1] == 1) continue;
      for (k = 0; k < fa->nfront; k++) {
        if (fa->front[k].sright == NOTFILLED) continue;
        if (k == fa->front[i].sright) continue;
        if (fa->front[k].i0 != i1) continue;
        i2 = fa->front[k].i1;
        uv0[0] = vertices[2*i0  ];
        uv0[1] = vertices[2*i0+1];
        uv1[0] = vertices[2*i1  ];
        uv1[1] = vertices[2*i1+1];
        uv2[0] = vertices[2*i2  ];
        uv2[1] = vertices[2*i2+1];
        area   = AREA2D(uv0, uv1, uv2);
        if ((neg == 0) && (area <= 0.0)) continue;
        for (l = 0; l < fa->nfront; l++) {
          if (fa->front[l].sright == NOTFILLED) continue;
          if (fa->front[l].sleft  == NOTFILLED) continue;
          if ((fa->front[l].i0 == i2) && (fa->front[l].i1 == i0)) {
            if (fa->front[i].sleft != l) {
              index = fa->front[i].sleft;
              indx2 = fa->front[l].sright;
              fa->front[i].sleft      = l;
              fa->front[l].sright     = i;
              fa->front[index].sright = indx2;
              fa->front[indx2].sleft  = index;
            }
            if (fa->front[i].sright != k) {
              index = fa->front[i].sright;
              indx2 = fa->front[k].sleft;
              fa->front[i].sright     = k;
              fa->front[k].sleft      = i;
              fa->front[index].sleft  = indx2;
              fa->front[indx2].sright = index;
            }
            if (fa->front[k].sright != l) {
              index = fa->front[k].sright;
              indx2 = fa->front[l].sleft;
              fa->front[k].sright     = l;
              fa->front[l].sleft      = k;
              fa->front[index].sleft  = indx2;
              fa->front[indx2].sright = index;
            }

            left   = fa->front[i].sleft;
            right  = fa->front[i].sright;
            triangles[3*ntri  ]    = i0;
            triangles[3*ntri+1]    = i1;
            triangles[3*ntri+2]    = i2;
            fa->front[i].sleft     = fa->front[i].sright     = NOTFILLED;
            fa->front[right].sleft = fa->front[right].sright = NOTFILLED;
            fa->front[left].sleft  = fa->front[left].sright  = NOTFILLED;
            ntri++;
            j++;
            if (ntri >= mtri) break;
            neg = 0;
          }
        }
        if (ntri >= mtri) break;
      }
      if (ntri >= mtri) break;
    }
    if (j != 0) continue;

    /* get smallest segment left */

    for (i = 0; i < fa->nfront; i++) fa->front[i].mark = 0;
small:
    index = -1;
    side2 = DBL_MAX;
    for (i = 0; i < fa->nfront; i++) {
      if (fa->front[i].sright == NOTFILLED) continue;
      if (fa->front[i].mark == 1) continue;
      i0     = fa->front[i].i0;
      i1     = fa->front[i].i1;
      uv0[0] = vertices[2*i0  ];
      uv0[1] = vertices[2*i0+1];
      uv1[0] = vertices[2*i1  ];
      uv1[1] = vertices[2*i1+1];
      d      = DIST2(uv0, uv1);
      if (d < side2) {
        side2 = d;
        index = i;
      }
    }
    if (index == -1) {
      for (k = 0; k < *n_fig8; k++) 
        if (ntri + 2*k == mtri) break;
      if (neg == 0) {
        neg = 1;
        continue;
      }
      printf(" EGADS Internal: can't find segment!\n");
      goto error;
    }

    /* find the best candidate -- closest to midpoint and correct area */

    i0     = fa->front[index].i0;
    i1     = fa->front[index].i1;
    uv0[0] = vertices[2*i0  ];
    uv0[1] = vertices[2*i0+1];
    uv1[0] = vertices[2*i1  ];
    uv1[1] = vertices[2*i1+1];
    mid[0] = 0.5*(uv0[0] + uv1[0]);
    mid[1] = 0.5*(uv0[1] + uv1[1]);

    indx2 = -1;
    dist  = DBL_MAX;
    for (i = 0; i < fa->nfront; i++) {
      if ((i == index) || (fa->front[i].sright == NOTFILLED)) continue;
      i2 = fa->front[i].i1;
      if ((i2 == i0) || (i2 == i1)) continue;
      uv2[0] = vertices[2*i2  ];
      uv2[1] = vertices[2*i2+1];
      area   = AREA2D(uv0, uv1, uv2);
      if (area > 0.0) {
        d = DIST2(mid, uv2)/area;
        if (d < dist) {
          if (EG_crossSeg(index, mid, i2, vertices, pass, fa)) continue;
          dist  = d;
          indx2 = i;
        }
      }
    }
    /* may not find a candidate for segments that are too small
               retry with next largest (and hope for closure later) */
    if (indx2 == -1) {
      fa->front[index].mark = 1;
      goto small;
    }

    /* construct the triangle */

    i2 = fa->front[indx2].i1;
    triangles[3*ntri  ] = i0;
    triangles[3*ntri+1] = i1;
    triangles[3*ntri+2] = i2;
    ntri++;
    neg = 0;

    /* patch up the front */

    left  = fa->front[index].sleft;
    right = fa->front[index].sright;

    if (i2 == fa->front[left].i0) {
      /* 1) candate is in the left segment */

      fa->front[left].sright = right;
      fa->front[left].i1     = i1;
      fa->front[left].snew   = 1;
      fa->front[right].sleft = left;
      fa->front[index].sleft = fa->front[index].sright = NOTFILLED;

    } else if (i2 == fa->front[right].i1) {
      /* 2) candate is in the right segment */

      fa->front[left].sright = right;
      fa->front[right].sleft = left;
      fa->front[right].i0    = i0;
      fa->front[right].snew  = 1;
      fa->front[index].sleft = fa->front[index].sright = NOTFILLED;

    } else {
      /* 3) some other situation */

      start = 0;

      /* "figure 8" vertices? */

      if (fa->pts[i0] != 1) 
        for (i = 0; i < fa->nfront; i++) {
          if (fa->front[i].sright == NOTFILLED) continue;
          if (fa->front[i].i0 != i2) continue;
          if (fa->front[i].i1 != i0) continue;
          j = fa->front[i].sright;
          fa->front[left].sright = j;
          fa->front[j].sleft     = left;
          fa->front[index].sleft = i;
          fa->front[i].sright    = index;
          left = i;
          fa->front[left].sright = right;
          fa->front[left].i1     = i1;
          fa->front[left].snew   = 1;
          fa->front[right].sleft = left;
          fa->front[index].sleft = fa->front[index].sright = NOTFILLED;
          start = 1;
          break;
        }

      if ((fa->pts[i1] != 1) && (start == 0))
        for (i = 0; i < fa->nfront; i++) {
          if (fa->front[i].sright == NOTFILLED) continue;
          if (fa->front[i].i0 != i1) continue;
          if (fa->front[i].i1 != i2) continue;
          j = fa->front[i].sleft;
          fa->front[right].sleft  = j;
          fa->front[j].sright     = right;
          fa->front[index].sright = i;
          fa->front[i].sleft      = index;
          right = i;
          fa->front[left].sright = right;
          fa->front[right].sleft = left;
          fa->front[right].i0    = i0;
          fa->front[right].snew  = 1;
          fa->front[index].sleft = fa->front[index].sright = NOTFILLED;
          start = 1;
          break;
        }

      /* no, add a segment */

      if (start == 0) {

        next = -1;
        for (i = 0; i < fa->nfront; i++)
          if (fa->front[i].sright == NOTFILLED) {
            next = i;
            break;
          }

        if (next == -1) {
          if (fa->nfront >= fa->mfront) {
            i = fa->mfront + CHUNK;
            tmp = (Front *) EG_reall(fa->front, i*sizeof(Front));
            if (tmp == NULL) return 0;
            itmp = (int *) EG_reall(fa->segs, 2*i*sizeof(int));
            if (itmp == NULL) {
              EG_free(tmp);
              return 0;
            }
            fa->mfront = i;
            fa->front  =  tmp;
            fa->segs   = itmp;
          }
          next = fa->nfront;
          fa->nfront++;
        }

        start = fa->front[indx2].sright;
        fa->front[index].i1     = i2;
        fa->front[index].sright = start;
        fa->front[index].snew   = 1;
        fa->front[start].sleft  = index;
        fa->front[indx2].sright = next;
        fa->front[right].sleft  = next;
        fa->front[next].sleft   = indx2;
        fa->front[next].i0      = i2;
        fa->front[next].i1      = i1;
        fa->front[next].sright  = right;
        fa->front[next].snew    = 1;
      }
    }

  } while (ntri < mtri);

error:
  for (j = i = 0; i < fa->nfront; i++) 
    if (fa->front[i].sright != NOTFILLED) j++;

  if (j != 0) {
#ifdef DEBUG
    printf(" EGADS Internal: # unused segments = %d\n", j);
#endif
    ntri = 0;
  }

  return ntri;
}


void
EG_makeConnect(int k1, int k2, int *tri, int *kedge, int *ntable,
               connect *etable, int face)
{
  int kn1, kn2, iface, oface, look;

  if (k1 > k2) {
    kn1 = k2-1;
    kn2 = k1-1;
  } else {
    kn1 = k1-1;
    kn2 = k2-1;
  }

  /* add to edge table */

  if (ntable[kn1] == NOTFILLED) {

    /* virgin node */

    *kedge               += 1;
    ntable[kn1]           = *kedge;
    etable[*kedge].node1  = kn1;
    etable[*kedge].node2  = kn2;
    etable[*kedge].tri    = tri;
    etable[*kedge].thread = NOTFILLED;
    return;

  } else {

    /* old node */
    iface = ntable[kn1];

again:
    if (etable[iface].node2 == kn2) {
      if (etable[iface].tri != NULL) {
        look = *etable[iface].tri;
        *etable[iface].tri = *tri;
        *tri = look;
        etable[iface].tri = NULL;
      } else {
        printf("EGADS Internal: Face %d", face);
        printf(", Side %d %d complete [but %d] (EG_makeConnect)!\n",
                k1+1, k2+1, *tri);
      }
      return;
    } else {
      oface = iface;
      iface = etable[oface].thread;
    }

    /* try next position in thread */
    if (iface == NOTFILLED) {
      *kedge               += 1;
      etable[oface].thread  = *kedge;
      etable[*kedge].node1  = kn1;
      etable[*kedge].node2  = kn2;
      etable[*kedge].tri    = tri;
      etable[*kedge].thread = NOTFILLED;
      return;
    }

    goto again;
  }
}


int
EG_makeNeighbors(triStruct *ts, int f)
{
  int     *ntab, nside, j;
  connect *etab;
  
  ntab = (int *) EG_alloc(ts->nverts*sizeof(int));
  if (ntab == NULL) {
    printf(" EGADS Error: Vert Table Malloc (EG_makeTessBody)!\n");
    return EGADS_MALLOC;    
  }
  etab = (connect *) EG_alloc(ts->ntris*3*sizeof(connect));
  if (etab == NULL) {
    printf(" EGADS Error: Edge Table Malloc (EG_makeTessBody)!\n");
    EG_free(ntab);
    return EGADS_MALLOC;    
  }

  nside = -1;
  for (j = 0; j < ts->nverts; j++) ntab[j] = NOTFILLED;
  for (j = 0; j < ts->ntris;  j++) {
    EG_makeConnect( ts->tris[j].indices[1], ts->tris[j].indices[2],
                   &ts->tris[j].neighbors[0], &nside, ntab, etab, f);
    EG_makeConnect( ts->tris[j].indices[0], ts->tris[j].indices[2], 
                   &ts->tris[j].neighbors[1], &nside, ntab, etab, f);
    EG_makeConnect( ts->tris[j].indices[0], ts->tris[j].indices[1],
                   &ts->tris[j].neighbors[2], &nside, ntab, etab, f);
  }

  for (j = 0; j < ts->nsegs; j++)
    EG_makeConnect( ts->segs[j].indices[0], ts->segs[j].indices[1],
                   &ts->segs[j].neighbor, &nside, ntab, etab, f);
                   
  /* report any unconnected triangle sides */

  for (j = 0; j <= nside; j++) {
    if (etab[j].tri == NULL) continue;
    printf(" EGADS Info: Face %d, Unconnected Side %d %d = %d\n",
           f, etab[j].node1+1, etab[j].node2+1, *etab[j].tri);
    *etab[j].tri = 0;
  }

  EG_free(etab);
  EG_free(ntab);
  return EGADS_SUCCESS;
}


static void
EG_updateTris(triStruct *ts, egTessel *btess, int fIndex)
{
  int    i, j, k, m, n, nf, edge, *ptype, *pindex, *frame, *frlps, *tris, *tric;
  double *xyz, *uv;
  
  xyz    = (double *) EG_alloc(3*ts->nverts*sizeof(double));
  uv     = (double *) EG_alloc(2*ts->nverts*sizeof(double));
  ptype  = (int *)    EG_alloc(  ts->nverts*sizeof(int));
  pindex = (int *)    EG_alloc(  ts->nverts*sizeof(int));
  frame  = (int *)    EG_alloc(3*ts->nframe*sizeof(int));
  frlps  = (int *)    EG_alloc(  ts->nloop* sizeof(int));
  tris   = (int *)    EG_alloc(3*ts->ntris* sizeof(int));
  tric   = (int *)    EG_alloc(3*ts->ntris* sizeof(int));
  if ((xyz    == NULL) || (uv    == NULL) || (ptype == NULL) ||
      (pindex == NULL) || (tris  == NULL) || (tric  == NULL) ||
      (frame  == NULL) || (frlps == NULL)) {
    printf(" EGADS Error: Cannot Allocate Tessellation Memory for %d!\n",
           fIndex);
    if (tric   != NULL) EG_free(tric);
    if (tris   != NULL) EG_free(tris);
    if (frlps  != NULL) EG_free(frlps);
    if (frame  != NULL) EG_free(frame);
    if (pindex != NULL) EG_free(pindex);
    if (ptype  != NULL) EG_free(ptype);
    if (uv     != NULL) EG_free(uv);
    if (xyz    != NULL) EG_free(xyz);
    return;
  }

  /* fix up the vertices */

  for (i = 0; i < ts->nverts; i++) {
    xyz[3*i  ] = ts->verts[i].xyz[0];
    xyz[3*i+1] = ts->verts[i].xyz[1];
    xyz[3*i+2] = ts->verts[i].xyz[2];
    uv[2*i  ]  = ts->verts[i].uv[0];
    uv[2*i+1]  = ts->verts[i].uv[1];
    ptype[i]   = pindex[i] = -1;
    if (ts->verts[i].type == NODE) {
      ptype[i]  = 0;
      pindex[i] = ts->verts[i].index;
    } else if (ts->verts[i].type == EDGE) {
      ptype[i]  = ts->verts[i].index;
      pindex[i] = ts->verts[i].edge;
    }
  }
  btess->tess2d[fIndex-1].xyz    = xyz;
  btess->tess2d[fIndex-1].uv     = uv;
  btess->tess2d[fIndex-1].ptype  = ptype;
  btess->tess2d[fIndex-1].pindex = pindex;
  btess->tess2d[fIndex-1].npts   = ts->nverts;
  
  /* fix up the initial triangulation frame */
  for (i = 0; i < 3*ts->nframe; i++) frame[i] = ts->frame[i];
  for (i = 0; i <   ts->nloop;  i++) {
    frlps[i] = ts->loop[i];
    if (i == 0) continue;
    frlps[i] += frlps[i-1];
  }
  
  /* fix up the triangles */
  
  for (i = 0; i < ts->ntris; i++) {
    tris[3*i  ] = ts->tris[i].indices[0];
    tris[3*i+1] = ts->tris[i].indices[1];
    tris[3*i+2] = ts->tris[i].indices[2];
    tric[3*i  ] = ts->tris[i].neighbors[0];
    tric[3*i+1] = ts->tris[i].neighbors[1];
    tric[3*i+2] = ts->tris[i].neighbors[2];
  }
  for (i = 0; i < ts->ntris; i++)
    for (j = 0; j < 3; j++) 
      if (tric[3*i+j] < 0) {
        n    = -tric[3*i+j];
        edge = abs(ts->segs[n-1].edge);
        k    = ts->segs[n-1].index-1;
        if (ts->segs[n-1].edge > 0) {
          nf = btess->tess1d[edge-1].faces[1].nface;
          m  = EG_faceConnIndex(btess->tess1d[edge-1].faces[1], fIndex);
          if (m == 0) {
            printf(" EGADS Internal: Face %d not found in Edge (+) %d!\n",
                   fIndex, edge);
          } else {
            btess->tess1d[edge-1].faces[1].tric[k*nf+m-1] = i+1;
          }
        } else {
          nf = btess->tess1d[edge-1].faces[0].nface;
          m  = EG_faceConnIndex(btess->tess1d[edge-1].faces[0], fIndex);
          if (m == 0) {
            printf(" EGADS Internal: Face %d not found in Edge (-) %d!\n",
                   fIndex, edge);
          } else {
            btess->tess1d[edge-1].faces[0].tric[k*nf+m-1] = i+1;
          }

        }
        tric[3*i+j] = -edge;
      }
  btess->tess2d[fIndex-1].frame  = frame;
  btess->tess2d[fIndex-1].frlps  = frlps;
  btess->tess2d[fIndex-1].tris   = tris;
  btess->tess2d[fIndex-1].tric   = tric;
  btess->tess2d[fIndex-1].nframe = ts->nframe;
  btess->tess2d[fIndex-1].nfrlps = ts->nloop;
  btess->tess2d[fIndex-1].ntris  = ts->ntris;
}


int
EG_getTessEdge(const egObject *tess, int indx, int *len,
               const double **xyz, const double **t)
{
  int          outLevel, index, stat, aType, alen;
  egTessel     *btess;
  egObject     *obj;
  const int    *ints;
  const double *reals;
  const char   *str;

  *len = 0;
  *xyz = *t = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  index = abs(indx);
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_getTessEdge)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_getTessEdge)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_getTessEdge)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_getTessEdge)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_getTessEdge)!\n");
    return EGADS_NODATA;  
  }
  if ((index < 1) || (index > btess->nEdge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_getTessEdge)!\n",
             index, btess->nEdge);
    return EGADS_INDEXERR;
  }
  
  /* mapping? */
  if (indx < 0) {
    stat = EG_attributeRet(obj, ".eMap", &aType, &alen, &ints, &reals, &str);
    if (stat == EGADS_SUCCESS) {
      if (aType != ATTRINT) {
        if (outLevel > 0)
          printf(" EGADS Error: Map Attribute Not an Int (EG_getTessEdge)!\n");
        return EGADS_ATTRERR;
      }
      if (alen != btess->nEdge) {
        if (outLevel > 0)
          printf(" EGADS Error: Len Mismatch %d %d (EG_getTessEdge)!\n",
                 alen, btess->nFace);
        return EGADS_TOPOERR;
      }
      index = ints[index-1];
      if ((index < 1) || (index > btess->nEdge)) {
        if (outLevel > 0)
          printf(" EGADS Error: Mapped Index = %d [1-%d] (EG_getTessEdge)!\n",
                 index, btess->nEdge);
        return EGADS_INDEXERR;
      }
    }
  }

  *len = btess->tess1d[index-1].npts;
  *xyz = btess->tess1d[index-1].xyz;
  *t   = btess->tess1d[index-1].t;
  return EGADS_SUCCESS;
}


int
EG_getTessFace(const egObject *tess, int indx, int *len, const double **xyz,
               const double **uv, const int **ptype, const int **pindex, 
               int *ntri, const int **tris, const int **tric)
{
  int          outLevel, index, stat, aType, alen;
  egTessel     *btess;
  egObject     *obj;
  const int    *ints;
  const double *reals;
  const char   *str;

  *len   = *ntri   = 0;
  *xyz   = *uv     = NULL;
  *ptype = *pindex = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  index = abs(indx);
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_getTessFace)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_getTessFace)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_getTessFace)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_getTessFace)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_getTessFace)!\n");
    return EGADS_NODATA;  
  }
  if ((index < 1) || (index > btess->nFace)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_getTessFace)!\n",
             index, btess->nFace);
    return EGADS_INDEXERR;
  }
  
  /* mapping? */
  if (indx < 0) {
    stat = EG_attributeRet(obj, ".fMap", &aType, &alen, &ints, &reals, &str);
    if (stat == EGADS_SUCCESS) {
      if (aType != ATTRINT) {
        if (outLevel > 0)
          printf(" EGADS Error: Map Attribute Not an Int (EG_getTessFace)!\n");
        return EGADS_ATTRERR;
      }
      if (alen != btess->nFace) {
        if (outLevel > 0)
          printf(" EGADS Error: Len Mismatch %d %d (EG_getTessFace)!\n",
                 alen, btess->nFace);
        return EGADS_TOPOERR;
      }
      index = ints[index-1];
      if ((index < 1) || (index > btess->nFace)) {
        if (outLevel > 0)
          printf(" EGADS Error: Mapped Index = %d [1-%d] (EG_getTessFace)!\n",
                 index, btess->nFace);
        return EGADS_INDEXERR;
      }
    }
  }

  *len    = btess->tess2d[index-1].npts;
  *xyz    = btess->tess2d[index-1].xyz;
  *uv     = btess->tess2d[index-1].uv;
  *ptype  = btess->tess2d[index-1].ptype;
  *pindex = btess->tess2d[index-1].pindex;
  *ntri   = btess->tess2d[index-1].ntris;
  *tris   = btess->tess2d[index-1].tris;
  *tric   = btess->tess2d[index-1].tric;
/*  test out frame data
  {
    int i, n;
    for (n = i = 0; i < 3*btess->tess2d[index-1].nframe; i++)
      if (btess->tess2d[index-1].frame[i] > n)
        n = btess->tess2d[index-1].frame[i];
    for (i = 0; i < n; i++)
      printf(" %d: %d %d\n", i+1, btess->tess2d[index-1].ptype[i],
             btess->tess2d[index-1].pindex[i]);
  }
*/
  return EGADS_SUCCESS;
}


int
EG_getTessLoops(const egObject *tess, int fIndex, int *nloop,
                const int **lIndices)
{
  egTessel *btess;
  egObject *obj;
  
  *nloop    = 0;
  *lIndices = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  btess = (egTessel *) tess->blind;
  if (btess == NULL)                return EGADS_NOTFOUND;
  obj = btess->src;
  if (obj == NULL)                  return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (obj->oclass != BODY)          return EGADS_NOTBODY;
  if (btess->tess2d == NULL)        return EGADS_NODATA;
  if ((fIndex < 1) ||
      (fIndex > btess->nFace))      return EGADS_INDEXERR;

  *nloop    = btess->tess2d[fIndex-1].nfrlps;
  *lIndices = btess->tess2d[fIndex-1].frlps;
  return EGADS_SUCCESS;
}


static int
EG_loopQuad(double edgeTol, int *lens, double *uvs, triStruct *ts)
{
  int    i, j, k, n, len, *types, edge, ne;
  double uvmm[4], uratio, vratio, edgeTOL = 0.05;

  if ((edgeTol >= 0.001) && (edgeTol <= 0.5)) edgeTOL = edgeTol;

  len  = 1;
  edge = ts->segs[0].edge;
  for (i = 1; i < ts->nsegs; i++)
    if (ts->segs[i].edge != edge) {
      len++;
      edge = ts->segs[i].edge;
    }
  
  if (len <= 4) {
    for (k = i = 0; i < len; i++) {
      lens[i] = 0;
      edge    = ts->segs[k].edge;
      while (edge == ts->segs[k].edge) {
        lens[i]++;
        k++;
        if (k == ts->nsegs) break;
      }
    }
    return len;
  }
  
  types = (int *) EG_alloc(2*len*sizeof(int));
  if (types == NULL) return 0;

  for (k = i = 0; i < len; i++) {
    n    = k;
    ne   = 1;
    edge = ts->segs[k].edge;
    while (edge == ts->segs[n].edge) {
      ne++;
      n++;
      if (n == ts->nsegs) break;
    }
    
    types[2*i  ] =   -1;
    types[2*i+1] = ne-1;
    
    uvmm[0] = uvs[2*k  ];
    uvmm[1] = uvs[2*k+1];
    uvmm[2] = uvs[2*k  ];
    uvmm[3] = uvs[2*k+1];
    k++;
    for (j = 1; j < ne; j++) {
      if ((i == len-1) && (j == ne-1)) k = 0;
      if (uvmm[0] > uvs[2*k  ]) uvmm[0] = uvs[2*k  ];
      if (uvmm[1] > uvs[2*k+1]) uvmm[1] = uvs[2*k+1];
      if (uvmm[2] < uvs[2*k  ]) uvmm[2] = uvs[2*k  ];
      if (uvmm[3] < uvs[2*k+1]) uvmm[3] = uvs[2*k+1];
      if (j != ne-1) k++;
    }
    uratio = (uvmm[2]-uvmm[0])/(ts->range[1]-ts->range[0]);
    vratio = (uvmm[3]-uvmm[1])/(ts->range[3]-ts->range[2]);
    if ((uratio < edgeTOL) && (vratio < edgeTOL))
      if (uratio > vratio) uratio = edgeTOL;
    
    if (uratio < edgeTOL) {
      types[2*i] = 0;
      if (uvmm[2] > 0.5*(ts->range[1]+ts->range[0])) types[2*i] = 2;
    } else if (vratio < edgeTOL) {
      types[2*i] = 1;
      if (uvmm[3] > 0.5*(ts->range[3]+ts->range[2])) types[2*i] = 3;
    }
    if (types[2*i] == -1) {
      /* not aligned */
      EG_free(types);
      return 0;
    }
  }
  
  k = 0;
  lens[k] = types[1];
  for (i = 1; i < len; i++) {
    if (types[2*i-2] == types[2*i]) {
      lens[k] += types[2*i+1];
    } else {
      k++;
      if (k == 5) {
        EG_free(types);
        return 0;
      }
      lens[k] = types[2*i+1];
    }
  }
  
  /* check for quad w/ wrap-around */
  if (k == 4)
    if (types[0] != types[2*len-2]) k = 0;
  
  /* currently we can not handle degenerates with long loops */
  if (k == 3)
    if (types[0] == types[2*len-2]) k = 0;
  
  EG_free(types);
  if ((k != 3) && (k != 4)) return 0;
  return k+1;
}


int
EG_fillTris(egObject *body, int iFace, egObject *face, egObject *tess, 
            triStruct *ts, fillArea *fa, long tID)
{
  int      i, j, k, m, n, stat, nedge, nloop, oclass, mtype, or, np, degen;
  int      *senses, *sns, *tris, *tmp, npts, ntot, ntri, sen, n_fig8, nd, st;
  int      outLevel, mm, mp, *lsenses, lor;
  double   range[4], trange[2], uvm[2], uvp[2], smallu, smallv, *uvs, *intEdg;
  egObject *geom, **edges, **loops, **nds;
  egTessel *btess;
  triVert  *tv;
  triTri   *tt;
  triSeg   *tsg;
  static double scl[3][2] = {1.0, 1.0,  10.0, 1.0,  0.1, 10.0};
  const  double *xyzs, *tps;

  outLevel = EG_outLevel(body);
  btess    = (egTessel *) tess->blind;

  /* get the Loops */

  stat = EG_getTopology(face, &geom, &oclass, &or, range, &nloop, &loops,
                        &lsenses);
  if (stat != EGADS_SUCCESS) return stat;
  smallu = 0.00005*(range[1] - range[0]);
  smallv = 0.00005*(range[3] - range[2]);
#ifdef DEBUG
  printf("%lX Face %d: nLoop = %d   Range = %lf %lf  %lf %lf\n",
         tID, iFace, nloop, range[0], range[1], range[2], range[3]);
#endif
  ts->fIndex   = iFace;
  ts->face     = face;
  ts->orUV     = or;
  ts->planar   = 0;
  ts->orCnt    = 0;
  ts->lens[0]  = ts->lens[1] = ts->lens[2] = ts->lens[3] = ts->lens[4] = 0;
  ts->range[0] = range[0];
  ts->range[1] = range[1];
  ts->range[2] = range[2];
  ts->range[3] = range[3];
  ts->uvs      = NULL;
  if (geom->mtype == PLANE) ts->planar = 1;

  /* get the point count */

  for (ntot = i = 0; i < nloop; i++) {
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, NULL, &nedge,
                          &edges, &senses);
    if (stat != EGADS_SUCCESS) return stat;
    for (j = 0; j < nedge; j++) {
      k = EG_indexBodyTopo(body, edges[j]);
      if (k <= EGADS_SUCCESS) {
        printf("%lX EGADS Error: Face %d -> Can not find Edge = %d!\n",
               tID, iFace, k);
        return EGADS_NOTFOUND;
      }
      stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, trange, &nd,
                            &nds, &sns);
      if (stat != EGADS_SUCCESS) return stat;
      if (mtype == DEGENERATE) continue;
      stat = EG_getTessEdge(tess, k, &npts, &xyzs, &tps);
      if (stat != EGADS_SUCCESS) return stat;
      ntot += npts-1;
    }
  }
  
  ntri = ntot-2 + 2*(nloop-1);
#ifdef DEBUG
  printf("%lX:       total points = %d,  total tris = %d\n", tID, ntot, ntri);
#endif

  /* get enough storage for the verts & boundary segs */
  
  n  = ntot/CHUNK + 1;
  n *= CHUNK;
  if (ts->verts == NULL) {
    ts->verts = (triVert *) EG_alloc(n*sizeof(triVert));
    if (ts->verts == NULL) return EGADS_MALLOC;
    ts->mverts = n;
  } else {
    if (n > ts->mverts) {
      tv = (triVert *) EG_reall(ts->verts, n*sizeof(triVert));
      if (tv == NULL) return EGADS_MALLOC;
      ts->verts  = tv;
      ts->mverts = n;
    }
  }
  ts->nverts = ntot;

  n  = ntot/CHUNK + 1;
  n *= CHUNK;
  if (ts->segs == NULL) {
    ts->segs = (triSeg *) EG_alloc(n*sizeof(triSeg));
    if (ts->segs == NULL) return EGADS_MALLOC;
    ts->msegs = n;
  } else {
    if (n > ts->msegs) {
      tsg = (triSeg *) EG_reall(ts->segs, n*sizeof(triSeg));
      if (tsg == NULL) return EGADS_MALLOC;
      ts->segs  = tsg;
      ts->msegs = n;
    }
  }
  ts->nsegs = ntot;
  
  /* get memory for the loops */

  if (nloop > ts->mloop) {
    if (ts->loop == NULL) {
      ts->loop = (int *) EG_alloc(nloop*sizeof(int));
      if (ts->loop == NULL) return EGADS_MALLOC;
#ifdef DEBUG
      printf(" Alloc Loop: with %d\n", nloop);
#endif
    } else {
      tmp = (int *) EG_reall(ts->loop, nloop*sizeof(int));
      if (tmp == NULL) return EGADS_MALLOC;
      ts->loop = tmp;
#ifdef DEBUG
      printf(" Realloc Loop: now %d (%d)\n", ts->mloop, nloop);
#endif
    }
    ts->mloop = nloop;
  }
  ts->nloop = nloop;
  uvs = (double *) EG_alloc((ntot*2+2)*sizeof(double));
  if (uvs == NULL) return EGADS_MALLOC;
 
  /* fill in the loops & mark the boundary segments */
  
  np     = 1;
  uvs[0] = uvs[1] = 0.0;
  for (i = 0; i < nloop; i++) {
    st   = np;
    stat = EG_getTopology(loops[i], &geom, &oclass, &mtype, NULL, &nedge,
                          &edges, &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(uvs);
      return stat;
    }
/*
    if ((nedge > 4) && (nloop == 1)) {
      printf("%lX EGADS Info: Face %d -- nedge = %d\n", tID, iFace, nedge);
      for (m = 0; m < nedge-1; m++)
        if (EG_isSame(edges[m], edges[m+1]) == EGADS_SUCCESS)
          printf("%lX EGADS Info: Face %d -> Edge %d is Same as %d!\n",
                 tID, iFace, i+1, i+2);
      if (EG_isSame(edges[nedge-1], edges[0]) == EGADS_SUCCESS)
          printf("%lX EGADS Info: Face %d -> Edge %d is Same as %d!\n",
                 tID, iFace, nedge, 1);
    }
*/
    lor = 1;
    if ((lsenses[i] == 2) || (lsenses[i] == -2)) lor = -1;
    n = 0;
    if (or*lor == SREVERSE) n = nedge-1;
    for (degen = ntot = j = 0; j < nedge; j++, n += or*lor) {
      k = EG_indexBodyTopo(body, edges[n]);
      if (k <= EGADS_SUCCESS) {
        printf("%lX EGADS Error: Face %d -> Can not find Edge = %d!\n",
               tID, iFace, k);
        EG_free(uvs);
        return EGADS_NOTFOUND;
      }
      stat = EG_getTopology(edges[n], &geom, &oclass, &mtype, trange, &nd,
                            &nds, &sns);
      if (stat != EGADS_SUCCESS) {
        EG_free(uvs);
        return stat;
      }
      if (mtype == DEGENERATE) {
        degen = 1;
        continue;
      }
      stat = EG_getTessEdge(tess, k, &npts, &xyzs, &tps);
      if (stat != EGADS_SUCCESS) {
        EG_free(uvs);
        return stat;
      }
      sen = senses[n]*or*lor;
    
      /* internal Edge? */
      intEdg = NULL;
      for (m = 0; m < nedge; m++) {
        if (m == n) continue;
        if (k == EG_indexBodyTopo(body, edges[m])) {
          uvm[0] = uvm[1] = -1.0;
          uvp[0] = uvp[1] =  1.0;
          EG_getEdgeUV(face, edges[n], -1, 0.5*(trange[0]+trange[1]), uvm);
          EG_getEdgeUV(face, edges[n],  1, 0.5*(trange[0]+trange[1]), uvp);
          if ((uvm[0] == uvp[0]) && (uvm[1] == uvp[1]) && (intEdg == NULL)) {
            printf(" %lX EGADS Info: ", tID);
            printf("Face #%d -> Edge #%d (%d) Internally in Loop %d %d, sen = %d!\n",
                   iFace, k, nedge, n+1, m+1, sen);
            intEdg = (double *) EG_alloc(4*npts*sizeof(double));
            if (intEdg == NULL) continue;
            for (m = 0; m < npts; m++) {
              stat = EG_getEdgeUV(face, edges[n], senses[n]*lor, tps[m],
                                  &intEdg[4*m]);
              if (stat != EGADS_SUCCESS) {
                printf("%lX EGADS Error: getEdgeUV! = %d  for Face %d/%d, Edge = %d\n",
                       tID, stat, iFace, i+1, n+1);
                EG_free(intEdg);
                EG_free(uvs);
                return stat;
              }
            }
            for (m = 0; m < npts; m++) {
              mm = m - 1;
              mp = m + 1;
              if (mm <  0)    mm = 0;
              if (mp >= npts) mp = npts-1;
              uvm[0] = intEdg[4*mp  ] - intEdg[4*mm  ];
              uvm[1] = intEdg[4*mp+1] - intEdg[4*mm+1];
              uvp[0] = atan2(-uvm[1], uvm[0]);
              intEdg[4*m+2] = sen*smallu*sin(uvp[0]);
              intEdg[4*m+3] = sen*smallv*cos(uvp[0]);
            }
          }
        }
      }
  
      /* fill frame for this Edge */
      if (sen == 1) {
        for (m = 0; m < npts-1; m++, np++) {
          ts->verts[np-1].type   = EDGE;
          ts->verts[np-1].edge   = k;
          ts->verts[np-1].index  = m+1;
          ts->verts[np-1].xyz[0] = xyzs[3*m  ];
          ts->verts[np-1].xyz[1] = xyzs[3*m+1];
          ts->verts[np-1].xyz[2] = xyzs[3*m+2];
          if (intEdg == NULL) {
            stat = EG_getEdgeUV(face, edges[n], senses[n]*lor, tps[m], &uvs[2*np]);
            if (stat != EGADS_SUCCESS) {
              printf("%lX EGADS Error: getEdgeUV+ = %d  for Face %d/%d, Edge = %d\n",
                     tID, stat, iFace, i+1, n+1);
              EG_free(uvs);
              return stat;
            }
            ts->verts[np-1].uv[0] = uvs[2*np  ];
            ts->verts[np-1].uv[1] = uvs[2*np+1];
          } else {
            ts->verts[np-1].uv[0] = intEdg[4*m  ];
            ts->verts[np-1].uv[1] = intEdg[4*m+1];
            uvs[2*np  ] = intEdg[4*m  ] + intEdg[4*m+2];
            uvs[2*np+1] = intEdg[4*m+1] + intEdg[4*m+3];
          }
          if (m == 0) {
            ts->verts[np-1].type   = NODE;
            ts->verts[np-1].edge   = 0;
            ts->verts[np-1].index  = EG_indexBodyTopo(body, nds[0]);
            if (degen == 1) {
#ifdef DEBUG
              printf("%lX  Face %d, Vertex %d: Node = %d is Degen!\n",
                     tID, iFace, np, ts->verts[np-1].index);
#endif
              ts->verts[np-1].edge = -1;
              degen = 0;
            }
          }
          ts->segs[np-1].indices[0] =  np;
          ts->segs[np-1].indices[1] =  np+1;
          ts->segs[np-1].neighbor   = -np;
          ts->segs[np-1].edge       =  senses[n]*lor*k;
          ts->segs[np-1].index      =  m+1;
#ifdef DEBUG
          printf("%lX:         %lf %lf\n", tID, uvs[2*np  ], uvs[2*np+1]);
#endif
        }
      } else {
        for (m = npts-1; m > 0; m--, np++) {
          ts->verts[np-1].type   = EDGE;
          ts->verts[np-1].edge   = k;
          ts->verts[np-1].index  = m+1;
          ts->verts[np-1].xyz[0] = xyzs[3*m  ];
          ts->verts[np-1].xyz[1] = xyzs[3*m+1];
          ts->verts[np-1].xyz[2] = xyzs[3*m+2];
          if (intEdg == NULL) {
            stat = EG_getEdgeUV(face, edges[n], senses[n]*lor, tps[m], &uvs[2*np]);
            if (stat != EGADS_SUCCESS) {
              printf("%lX EGADS Error: getEdgeUV- = %d  for Face %d/%d, Edge = %d\n",
                     tID, stat, iFace, i+1, n+1);
              EG_free(uvs);
              return stat;
            }
            ts->verts[np-1].uv[0] = uvs[2*np  ];
            ts->verts[np-1].uv[1] = uvs[2*np+1];
          } else {
            ts->verts[np-1].uv[0] = intEdg[4*m  ];
            ts->verts[np-1].uv[1] = intEdg[4*m+1];
            uvs[2*np  ] = intEdg[4*m  ] + intEdg[4*m+2];
            uvs[2*np+1] = intEdg[4*m+1] + intEdg[4*m+3];
          }
          if (m == npts-1) {
            ts->verts[np-1].type   = NODE;
            ts->verts[np-1].edge   = 0;
            if (mtype == TWONODE) {
              ts->verts[np-1].index = EG_indexBodyTopo(body, nds[1]);
            } else {
              ts->verts[np-1].index = EG_indexBodyTopo(body, nds[0]);
            }
            if (degen == 1) {
#ifdef DEBUG
              printf("%lX: Face %d, Vertex %d: Node = %d is Degen!\n",
                     tID, iFace, np, ts->verts[np-1].index);
#endif
              ts->verts[np-1].edge = -1;
              degen = 0;
            }
          }
          ts->segs[np-1].indices[0] =  np;
          ts->segs[np-1].indices[1] =  np+1;
          ts->segs[np-1].neighbor   = -np;
          ts->segs[np-1].edge       =  senses[n]*lor*k;
          ts->segs[np-1].index      =  m;
#ifdef DEBUG
          printf("%lX:          %lf %lf\n", tID, uvs[2*np  ], uvs[2*np+1]);
#endif
        }
      }
#ifdef DEBUG
      printf("%lX **** End Edge %d sen = %d ****\n", tID, k+1, sen);
#endif
      if (intEdg != NULL) EG_free(intEdg);
      ntot += npts-1;
    }
    if (np > 1) ts->segs[np-2].indices[1] = st;
    if (degen == 1) {
      if (ts->verts[st-1].edge != 0) {
        printf("%lX EGADS Error: Degen setting w/ Face %d  Marker = %d %d %d\n",
               tID, iFace, ts->verts[st-1].type, ts->verts[st-1].edge,
                           ts->verts[st-1].index);
      } else {
#ifdef DEBUG
        printf("%lX Face %d, Vertex %d: Node = %d is Degen!\n",
               tID, iFace, st, ts->verts[st-1].index);
#endif
        ts->verts[st-1].edge = -1;
      }
    }
#ifdef DEBUG
    printf("%lX **** End Loop %d: nedge = %d  %d ****\n",
           tID, i+1, nedge, ntot);
#endif
    ts->loop[i] = ntot;
  }
  
  /* handle quadding */
  if ((nloop == 1) && (ts->qparm[0] >= 0.0)) {
    /* avoid all together? */
    stat = EG_loopQuad(ts->qparm[0], ts->lens, &uvs[2], ts);
    if ((stat >= 3) && (stat <= 5) && (outLevel > 1))
      printf("%lX Face %d: quad = %d -- %d %d %d %d %d = %d\n", tID, iFace,
             stat, ts->lens[0], ts->lens[1], ts->lens[2], ts->lens[3],
             ts->lens[4], ntot);
    if ((stat != 3) && (stat != 4)) {
      ts->lens[0] = ts->lens[1] = ts->lens[2] = ts->lens[3] = ts->lens[4] = 0;
    } else {
      ts->uvs = uvs;
    }
  }
  
  /* fill in the interior with triangles */
  
  tris = (int *) EG_alloc(3*ntri*sizeof(int));
  if (tris == NULL) {
    /*@-kepttrans@*/
    EG_free(uvs);
    /*@+kepttrans@*/
    return EGADS_MALLOC;
  }
  
  n = EG_fillArea(nloop, ts->loop, uvs, tris, &n_fig8, 0, fa);
  
  /* adjust for figure 8 configurations */
  if (n_fig8 != 0) {
    printf("%lX EG_fillArea Warning: Face %d -> Found %d figure 8's!\n",
           tID, iFace, n_fig8);
    for (i = 0; i < n_fig8; i++) if (n+2*i == ntri) ntri = n;
  }
#ifdef DEBUG  
  printf("%lX:   EG_fillArea = %d (%d),  #loops = %d, or = %d,  #fig8 = %d\n", 
         tID, n, ntri, nloop, or, n_fig8);
#endif
         
  if (n != ntri) {
    range[0] = range[2] = uvs[2];
    range[1] = range[3] = uvs[3];
    for (i = 2; i < np; i++) {
      if (uvs[2*i  ] < range[0]) range[0] = uvs[2*i  ];
      if (uvs[2*i+1] < range[1]) range[1] = uvs[2*i+1];
      if (uvs[2*i  ] > range[2]) range[2] = uvs[2*i  ];
      if (uvs[2*i+1] > range[3]) range[3] = uvs[2*i+1];
    }
    for (i = 1; i < np; i++) {
      uvs[2*i  ] = (uvs[2*i  ]-range[0])/(range[2]-range[0]);
      uvs[2*i+1] = (uvs[2*i+1]-range[1])/(range[3]-range[1]);
    }
    for (j = 0; j < 3; j++) {
      for (i = 1; i < np; i++) {
        uvs[2*i  ] *= scl[j][0];
        uvs[2*i+1] *= scl[j][1];
      }
      n = EG_fillArea(nloop, ts->loop, uvs, tris, &n_fig8, 1, fa);
      printf("%lX EGADS Internal: Face %d -> Renormalizing %d, ntris = %d (%d)!\n",
             tID, iFace, j, ntri, n);
      if (n == ntri) break;
    }
  }
  /*@-kepttrans@*/
  if (ts->uvs == NULL) EG_free(uvs);
  /*@+kepttrans@*/
  if (n != ntri) {
    EG_free(tris);
    if (ts->uvs != NULL) EG_free(ts->uvs);
    return EGADS_DEGEN;
  }
  
  /* fill up the triangles */
  
  n  = ntri/CHUNK + 1;
  n *= CHUNK;
  if (ts->tris == NULL) {
    ts->tris = (triTri *) EG_alloc(n*sizeof(triTri));
    if (ts->tris == NULL) {
      EG_free(tris);
      if (ts->uvs != NULL) EG_free(ts->uvs);
      return EGADS_MALLOC;
    }
    ts->mtris = n;
  } else {
    if (n > ts->mtris) {
      tt = (triTri *) EG_reall(ts->tris, n*sizeof(triTri));
      if (tt == NULL) {
        EG_free(tris);
        if (ts->uvs != NULL) EG_free(ts->uvs);
        return EGADS_MALLOC;
      }
      ts->tris  = tt;
      ts->mtris = n;
    }
  } 
  
  for (i = 0; i < ntri; i++) {
    ts->tris[i].mark         = 0;
    ts->tris[i].indices[0]   = tris[3*i  ];
    ts->tris[i].indices[1]   = tris[3*i+1];
    ts->tris[i].indices[2]   = tris[3*i+2];
    ts->tris[i].neighbors[0] = i+1;
    ts->tris[i].neighbors[1] = i+1;
    ts->tris[i].neighbors[2] = i+1;
  }
  ts->ntris = ntri;
  EG_free(tris);
  
  /* flip tri orientation if face is reversed */
  if (or == SREVERSE)
    for (i = 0; i < ts->ntris; i++) {
      j                      = ts->tris[i].indices[1];
      ts->tris[i].indices[1] = ts->tris[i].indices[2];
      ts->tris[i].indices[2] = j;
    }
    
  /* connect the triangles and make the neighbor info */
  
  stat = EG_makeNeighbors(ts, iFace);
  if (stat != EGADS_SUCCESS) {
    if (ts->uvs != NULL) EG_free(ts->uvs);
    return stat;
  }

  /* enhance the tessellation */

  stat = EG_tessellate(outLevel, ts, tID);
  if (stat == EGADS_SUCCESS) {
    /* set it in the tessellation structure */
    EG_updateTris(ts, btess, iFace);
  }

  return stat;
}


void
EG_cleanupTess(egTessel *btess)
{
  int i;
  
  if (btess->xyzs != NULL) EG_free(btess->xyzs);
  
  if (btess->tess1d != NULL) {
    for (i = 0; i < btess->nEdge; i++) {
      if (btess->tess1d[i].faces[0].faces != NULL)
        EG_free(btess->tess1d[i].faces[0].faces);
      if (btess->tess1d[i].faces[1].faces != NULL)
        EG_free(btess->tess1d[i].faces[1].faces);
      if (btess->tess1d[i].faces[0].tric  != NULL)
        EG_free(btess->tess1d[i].faces[0].tric);
      if (btess->tess1d[i].faces[1].tric  != NULL)
        EG_free(btess->tess1d[i].faces[1].tric);
      if (btess->tess1d[i].xyz    != NULL)
        EG_free(btess->tess1d[i].xyz);
      if (btess->tess1d[i].t      != NULL)
        EG_free(btess->tess1d[i].t);
      if (btess->tess1d[i].global != NULL)
        EG_free(btess->tess1d[i].global);
    }
    EG_free(btess->tess1d);
  }
  
  if (btess->tess2d != NULL) {
    for (i = 0; i < 2*btess->nFace; i++) {
      if (btess->tess2d[i].mKnots != NULL)
        EG_deleteObject(btess->tess2d[i].mKnots);
      if (btess->tess2d[i].xyz    != NULL) 
        EG_free(btess->tess2d[i].xyz);
      if (btess->tess2d[i].uv     != NULL) 
        EG_free(btess->tess2d[i].uv);
      if (btess->tess2d[i].global != NULL)
        EG_free(btess->tess2d[i].global);
      if (btess->tess2d[i].ptype  != NULL) 
        EG_free(btess->tess2d[i].ptype);
      if (btess->tess2d[i].pindex != NULL) 
        EG_free(btess->tess2d[i].pindex);
      if (btess->tess2d[i].bary   != NULL)
        EG_free(btess->tess2d[i].bary);
      if (btess->tess2d[i].frame  != NULL)
        EG_free(btess->tess2d[i].frame);
      if (btess->tess2d[i].frlps  != NULL)
        EG_free(btess->tess2d[i].frlps);
      if (btess->tess2d[i].tris   != NULL) 
        EG_free(btess->tess2d[i].tris);
      if (btess->tess2d[i].tric   != NULL) 
        EG_free(btess->tess2d[i].tric);
    }
    EG_free(btess->tess2d);
  }
  if (btess->globals != NULL) EG_free(btess->globals);

}


static double
EG_curvNorm(egObject *face, double *uv, int sense, double d, double *dx)
{
#ifndef USEFACPROJ
  int    stat;
  double rad, c2, s2, curv, x1[3], x2[3], nrme[3], ds[3], result[18];

  /* get normal at mid-point in Edge segment based on Face UVs */
  stat  = EG_evaluate(face, uv, result);
  if (stat != EGADS_SUCCESS) return -2.0;
  x1[0] = result[3];
  x1[1] = result[4];
  x1[2] = result[5];
  x2[0] = result[6];
  x2[1] = result[7];
  x2[2] = result[8];
  CROSS(nrme, x1, x2);
  rad   = DOT(nrme, nrme);
  if (rad == 0.0) return -2.0;
  rad      = 1.0/sqrt(rad);
  nrme[0] *= rad;
  nrme[1] *= rad;
  nrme[2] *= rad;
  
  /* get segment-distance (rad) and orthogonal-direction to Edge segment */
  rad    = sqrt(d);
  dx[0] /= rad;
  dx[1] /= rad;
  dx[2] /= rad;
  CROSS(ds, dx, nrme);
  if (sense == 1) {
    ds[0] = -ds[0];
    ds[1] = -ds[1];
    ds[2] = -ds[2];
  }

  /* use Face curvature at Edge to estimate interior dot of normals */
  stat  = EG_curvature(face, uv, result);
  if (stat != EGADS_SUCCESS) return -2.0;
  x1[0] = result[1];
  x1[1] = result[2];
  x1[2] = result[3];
  c2    = DOT(ds,x1)*DOT(ds,x1);
  s2    = 1.0 - c2;
  curv  = result[0]*c2 + result[4]*s2;
  /* 1/4 segment-distance is used as orthogonal length
     and remember curvature = 1/RadOfCurv */
  rad  *= 0.25*fabs(curv);
  if (rad > PI) rad = PI;
  return cos(rad);
  
#else
  
  int      stat, count, oclass, mtype, nc, *sen;
  double   area, a00, a10, a11, b0, b1, det, dist, ldist, lims[4];
  double   result[18], x1[3], x2[3], nrme[3], nrmi[3], ds[3], dX[3];
  egObject *ref, **chld;
  
  /* get face limits */
  stat = EG_getTopology(face, &ref, &oclass, &mtype, lims, &nc, &chld, &sen);
  if (stat != EGADS_SUCCESS) return -2.0;
  
  /* get normal at mid-point in UV */
  stat  = EG_evaluate(face, uv, result);
  if (stat != EGADS_SUCCESS) return -2.0;
  x1[0] = result[3];
  x1[1] = result[4];
  x1[2] = result[5];
  x2[0] = result[6];
  x2[1] = result[7];
  x2[2] = result[8];
  CROSS(nrme, x1, x2);
  area = DOT(nrme, nrme);
  if (area == 0.0) return -2.0;
  area     = 1.0/sqrt(area);
  nrme[0] *= area;
  nrme[1] *= area;
  nrme[2] *= area;
            
  /* get interior Face normal */
  area   = sqrt(d);
  dx[0] /= area;
  dx[1] /= area;
  dx[2] /= area;
  CROSS(ds, dx, nrme);
  if (sense == 1) {
    ds[0] = -ds[0];
    ds[1] = -ds[1];
    ds[2] = -ds[2];
  }

  /* target interior position */
  area /= 4.0;
  x1[0] = result[0] + area*ds[0];
  x1[1] = result[1] + area*ds[1];
  x1[2] = result[2] + area*ds[2];

  /* newton raphson nearest loop -- avoid invEvaluate if possible */
  for (count = 0; count < 20; count++) {
    stat  = EGADS_SUCCESS;
    dX[0] = result[0] - x1[0];
    dX[1] = result[1] - x1[1];
    dX[2] = result[2] - x1[2];
    dist  = sqrt(dX[0]*dX[0] + dX[1]*dX[1] + dX[2]*dX[2]);
    if (count != 0) {
      if (fabs(ldist-dist) < 1.e-7) break;
      /* diverging */
      if ((ldist < dist) && (count > 1)) {
        stat = EGADS_RANGERR;
        break;
      }
    }

    b0  = -dX[0]*result[3] - dX[1]*result[4] - dX[2]*result[5];
    b1  = -dX[0]*result[6] - dX[1]*result[7] - dX[2]*result[8];
    a00 = result[ 3]*result[ 3] + result[ 4]*result[ 4] + result[ 5]*result[ 5] +
               dX[0]*result[ 9] +      dX[1]*result[10] +      dX[2]*result[11];
    a10 = result[ 3]*result[ 6] + result[ 4]*result[ 7] + result[ 5]*result[ 8] +
               dX[0]*result[12] +      dX[1]*result[13] +      dX[2]*result[14];
    a11 = result[ 6]*result[ 6] + result[ 7]*result[ 7] + result[ 8]*result[ 8] +
               dX[0]*result[15] +      dX[1]*result[16] +      dX[2]*result[17];
    
    det = a00*a11 - a10*a10;
    if (det == 0.0) {
      stat = EGADS_DEGEN;
      break;
    }
    det    = 1.0/det;
    uv[0] += det*(b0*a11 - b1*a10);
    uv[1] += det*(b1*a00 - b0*a10);
    ldist  = dist;
    stat = EGADS_OUTSIDE;
    if (uv[0] < lims[0]) break;
    if (uv[0] > lims[1]) break;
    if (uv[1] < lims[2]) break;
    if (uv[1] > lims[3]) break;
    stat = EG_evaluate(face, uv, result);
    if (stat != EGADS_SUCCESS) break;
  }

  if ((stat != EGADS_SUCCESS) || (count == 20)) {
    stat = EG_invEvaluate(face, x1, uv, x2);
    if (stat != EGADS_SUCCESS) return -2.0;
    stat = EG_evaluate(face, uv, result);
    if (stat != EGADS_SUCCESS) return -2.0;
  }

  x1[0] = result[3];
  x1[1] = result[4];
  x1[2] = result[5];
  x2[0] = result[6];
  x2[1] = result[7];
  x2[2] = result[8];
  CROSS(nrmi, x1, x2);
  area = DOT(nrmi, nrmi);
  if (area == 0.0) return -2.0;
  area     = 1.0/sqrt(area);
  nrmi[0] *= area;
  nrmi[1] *= area;
  nrmi[2] *= area;
            
  /* dot the normals */
  return DOT(nrme, nrmi);
#endif
}


static int
EG_otherFaces(egTess1D tess1d, ego *faces, int facex, ego edge,
              double t0, double th, double t1)
{
  int    n, nf, face, sense, stat;
  double uv0[2], uvh[2], uv1[2], dirw[2], dot;
  
  for (n = 0; n < 2; n++) {
                sense =  1;
    if (n == 0) sense = -1;
    for (nf = 0; nf < tess1d.faces[n].nface; nf++) {
      face = tess1d.faces[n].index;
      if (tess1d.faces[n].nface > 1) face = tess1d.faces[n].faces[nf];
      if ((face <= 0) || (face == facex)) continue;
      stat = EG_getEdgeUV(faces[face-1], edge, sense, t0, uv0);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_getEdgeUV(faces[face-1], edge, sense, th, uvh);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_getEdgeUV(faces[face-1], edge, sense, t1, uv1);
      if (stat != EGADS_SUCCESS) return stat;
      dirw[0] = uv1[0] - uv0[0];
      dirw[1] = uv1[1] - uv0[1];
      dot = (uvh[0] - uv0[0])*dirw[0] + (uvh[1] - uv0[1])*dirw[1];
      if (dot < 0.0) return EGADS_DEGEN;
      dot = (uv1[0] - uvh[0])*dirw[0] + (uv1[1] - uvh[1])*dirw[1];
      if (dot < 0.0) return EGADS_DEGEN;
    }
  }
  
  return EGADS_SUCCESS;
}


static int
EG_tessEdge(egTessel *btess, egObject **faces, int j, egObject *edge, long tID)
{
  int      i, k, n, npts, stat, outLevel, oclass, mtype, nnode, btype, *info;
  int      nf, ntype, ndum, face, sense, *senses, aStat, aType, aLen;
  double   xyz[MAXELEN][3], t[MAXELEN], aux[MAXELEN][3], mindist, tol, dist;
  double   d, dotnrm, dot, limits[2], mid[3], range[4], dx[3], result[18];
  double   params[3], uv[2], *prv;
  egObject *body, *geom, *ref, *rref, **nodes, **dum;
#ifdef SPLITDEGEN
  int      be[2];
#endif
  const int    *aInts;
  const double *aReals;
  const char   *aStr;
  
  body      = btess->src;
  outLevel  = EG_outLevel(body);
  params[0] = btess->params[0];
  params[1] = btess->params[1];
  params[2] = btess->params[2];
#ifdef PROGRESS
  if (outLevel > 0) {
    printf("    tessellating Edge %3d of %3d\r", j+1, btess->nEdge);
    fflush(stdout);
  }
#endif
  
  /* adjust the parameters? */
  stat = EG_getBodyTopos(body, edge, FACE, &n, &dum);
  if (stat == EGADS_SUCCESS) {
    for (i = 0; i < n; i++) {
      aStat = EG_attributeRet(dum[i], ".tParams", &aType, &aLen, &aInts,
                              &aReals, &aStr);
      if (aStat == EGADS_SUCCESS)
        if (aType != ATTRREAL) {
          printf("%lX EGADS Warning: tParams NONReal Attribute (EG_tessEdges)!\n",
                 tID);
        } else {
          for (k = 0; k < aLen; k++) {
            if (k == 3) break;
            if ((aReals[k] < params[k]) && (aReals[k] > 0.0))
              params[k] = aReals[k];
          }
        }
      aStat = EG_attributeRet(dum[i], ".tParam", &aType, &aLen, &aInts,
                              &aReals, &aStr);
      if (aStat == EGADS_SUCCESS)
        if (aType != ATTRREAL) {
          printf("%lX EGADS Warning: tParam NONReal Attribute (EG_tessEdges)!\n",
                 tID);
        } else {
          for (k = 0; k < aLen; k++) {
            if (k == 3) break;
            if (aReals[k] > 0.0) params[k] = aReals[k];
          }
        }
    }
    EG_free(dum);
  }
  aStat = EG_attributeRet(edge, ".tParams", &aType, &aLen, &aInts,
                          &aReals, &aStr);
  if (aStat == EGADS_SUCCESS)
    if (aType != ATTRREAL) {
      printf("%lX EGADS Warning: tParams NonReal Attribute (EG_tessEdges)!\n",
             tID);
    } else {
      for (i = 0; i < aLen; i++) {
        if (i == 3) break;
        if ((aReals[i] < params[i]) && (aReals[i] > 0.0))
          params[i] = aReals[i];
      }
    }
  aStat = EG_attributeRet(edge, ".tParam", &aType, &aLen, &aInts,
                          &aReals, &aStr);
  if (aStat == EGADS_SUCCESS)
    if (aType != ATTRREAL) {
      printf("%lX EGADS Warning: tParam NonReal Attribute (EG_tessEdges)!\n",
             tID);
    } else {
      for (i = 0; i < aLen; i++) {
        if (i == 3) break;
        if (aReals[i] > 0.0) params[i] = aReals[i];
      }
    }
  
  dist = fabs(params[2]);
  if (dist > 30.0) dist = 30.0;
  if (dist <  0.5) dist =  0.5;
  dotnrm = cos(PI*dist/180.0);
  
  stat = EG_getTopology(edge, &geom, &oclass, &mtype, limits,
                        &nnode, &nodes, &senses);
  if (stat != EGADS_SUCCESS) return stat;
#ifdef DEBUG
  printf("%lX Edge %d: type = %d, geom type = %d, limits = %lf %lf, nnode = %d\n",
         tID, j+1, mtype, geom->mtype, limits[0], limits[1], nnode);
#endif
  
  /* set end points */
  stat = EG_getTopology(nodes[0], &ref, &oclass, &ntype, xyz[0],
                        &ndum, &dum, &senses);
  if (stat != EGADS_SUCCESS) return stat;
  npts      = 2;
  t[0]      = limits[0];
  xyz[1][0] = xyz[0][0];
  xyz[1][1] = xyz[0][1];
  xyz[1][2] = xyz[0][2];
  t[1]      = limits[1];
  btess->tess1d[j].nodes[0] = EG_indexBodyTopo(body, nodes[0]);
  btess->tess1d[j].nodes[1] = btess->tess1d[j].nodes[0];
  if (mtype == TWONODE) {
    stat = EG_getTopology(nodes[1], &ref, &oclass, &ntype, xyz[1],
                          &ndum, &dum, &senses);
    if (stat != EGADS_SUCCESS) return stat;
    btess->tess1d[j].nodes[1] = EG_indexBodyTopo(body, nodes[1]);
  }
  
  /* degenerate -- finish up */
  if (mtype == DEGENERATE) {
    btess->tess1d[j].xyz = (double *) EG_alloc(3*npts*sizeof(double));
    if (btess->tess1d[j].xyz == NULL) {
      if (outLevel > 0)
        printf("%lX EGADS Error: Alloc %d Pts Edge %d (EG_tessEdges)!\n",
               tID, npts, j+1);
      return EGADS_MALLOC;
    }
    btess->tess1d[j].t = (double *) EG_alloc(npts*sizeof(double));
    if (btess->tess1d[j].t == NULL) {
      EG_free(btess->tess1d[j].xyz);
      btess->tess1d[j].xyz = NULL;
      if (outLevel > 0)
        printf("%lX EGADS Error: Alloc %d Ts Edge %d (EG_tessEdges)!\n",
               tID, npts, j+1);
      return EGADS_MALLOC;
    }
    for (i = 0; i < npts; i++) {
      btess->tess1d[j].xyz[3*i  ] = xyz[i][0];
      btess->tess1d[j].xyz[3*i+1] = xyz[i][1];
      btess->tess1d[j].xyz[3*i+2] = xyz[i][2];
      btess->tess1d[j].t[i]       = t[i];
    }
    btess->tess1d[j].npts = npts;
    return EGADS_SUCCESS;
  }
  
  /* should we use the point distribution given us? */
  aStat = EG_attributeRet(edge, ".tPos", &aType, &aLen, &aInts, &aReals, &aStr);
  if (aStat == EGADS_SUCCESS)
    if ((aType == ATTRSTRING) || (aType == ATTRCSYS)) {
      printf("%lX EGADS Warning: tPos NonReal Attribute (EG_tessEdges)!\n",
             tID);
    } else if (aType == ATTRINT) {
      if ((aLen == 1) && (aInts[0] == 0)) goto fill1D;
      printf("%lX EGADS Warning: tPos NonReal Attribute (EG_tessEdges)!\n",
             tID);
    } else {
      for (k = i = 0; i < aLen; i++)
        if ((aReals[i] <= t[0]) || (aReals[i] >= t[1])) {
          printf("%lX EGADS Error: %d tPos[%d] = %lf [%lf-%lf] (EG_tessEdges)!\n",
                 tID, j+1, i, aReals[i], t[0], t[1]);
          k++;
        }
      if (k == 0)
        for (i = 0; i < aLen-1; i++)
          if (aReals[i] >= aReals[i+1]) {
            printf("%lX EGADS Error: %d tPos[%d] = %lf %lf NotMono (EG_tessEdges)!\n",
                   tID, j+1, i, aReals[i], aReals[i+1]);
            k++;
          }
      if (aLen+2 > MAXELEN) {
        printf("%lX EGADS Error: Edge %d Attr len %d > %d (EG_tessEdges)!\n",
               tID, j+1, aLen, MAXELEN-2);
        k++;
      }
      if (k == 0) {
        xyz[aLen+1][0] = xyz[1][0];
        xyz[aLen+1][1] = xyz[1][1];
        xyz[aLen+1][2] = xyz[1][2];
        t[aLen+1]      = t[1];
        for (i = 0; i < aLen; i++) {
          t[i+1] = aReals[i];
          stat   = EG_evaluate(edge, &t[i+1], result);
          if (stat != EGADS_SUCCESS) return stat;
          xyz[i+1][0] = result[0];
          xyz[i+1][1] = result[1];
          xyz[i+1][2] = result[2];
        }
        npts = aLen+2;
        goto fill1D;
      }
    }
  aStat = EG_attributeRet(edge, ".rPos", &aType, &aLen, &aInts, &aReals, &aStr);
  if (aStat == EGADS_SUCCESS)
    if ((aType == ATTRSTRING) || (aType == ATTRCSYS)) {
      printf("%lX EGADS Warning: rPos NonReal Attribute (EG_tessEdges)!\n",
             tID);
    } else if (aType == ATTRINT) {
      if ((aLen == 1) && (aInts[0] == 0)) goto fill1D;
      printf("%lX EGADS Warning: rPos NonReal Attribute (EG_tessEdges)!\n",
             tID);
    } else {
      for (k = i = 0; i < aLen; i++)
        if ((aReals[i] <= 0.0) || (aReals[i] >= 1.0)) {
          printf("%lX EGADS Error: %d rPos[%d] = %lf [0 - 1] (EG_tessEdges)!\n",
                 tID, j+1, i, aReals[i]);
          k++;
        }
      if (k == 0)
        for (i = 0; i < aLen-1; i++)
          if (aReals[i] >= aReals[i+1]) {
            printf("%lX EGADS Error: %d rPos[%d] = %lf %lf NotMono (EG_tessEdges)!\n",
                   tID, j+1, i, aReals[i], aReals[i+1]);
            k++;
          }
      if (aLen+2 > MAXELEN) {
        printf("%lX EGADS Error: Edge %d Attr Len %d > %d (EG_tessEdges)!\n",
               tID, j+1, aLen, MAXELEN-2);
        k++;
      }
      if (k == 0) {
        xyz[aLen+1][0] = xyz[1][0];
        xyz[aLen+1][1] = xyz[1][1];
        xyz[aLen+1][2] = xyz[1][2];
        t[aLen+1]      = t[1];
/*
        for (i = 0; i < aLen; i++) {
          t[i+1] = t[0] + aReals[i]*(t[aLen+1] - t[0]);
          stat   = EG_evaluate(edge, &t[i+1], result);
          if (stat != EGADS_SUCCESS) return stat;
          xyz[i+1][0] = result[0];
          xyz[i+1][1] = result[1];
          xyz[i+1][2] = result[2];
        }
 */
        stat = EG_relPosTs(edge, aLen+2, aReals, t, (double *) xyz);
        if (stat == EGADS_SUCCESS) {
          npts = aLen+2;
          goto fill1D;
        }
        printf("%lX EGADS Error: EG_relPosTs = %d (EG_tessEdges)!\n",
               tID, stat);
      }
    }
  
  /* get minimum distance */
  stat = EG_evaluate(edge, &t[0], result);
  if (stat != EGADS_SUCCESS) return stat;
  mindist = (xyz[0][0]-result[0])*(xyz[0][0]-result[0]) +
            (xyz[0][1]-result[1])*(xyz[0][1]-result[1]) +
            (xyz[0][2]-result[2])*(xyz[0][2]-result[2]);
  dist = sqrt(result[3]*result[3] + result[4]*result[4] +
              result[5]*result[5]);
  if (dist == 0.0) dist = 1.0;
  aux[0][0] = result[3]/dist;
  aux[0][1] = result[4]/dist;
  aux[0][2] = result[5]/dist;
  stat = EG_evaluate(edge, &t[1], result);
  if (stat != EGADS_SUCCESS) return stat;
  dist = sqrt(result[3]*result[3] + result[4]*result[4] +
              result[5]*result[5]);
  if (dist == 0.0) dist = 1.0;
  aux[1][0] = result[3]/dist;
  aux[1][1] = result[4]/dist;
  aux[1][2] = result[5]/dist;
  dist = (xyz[1][0]-result[0])*(xyz[1][0]-result[0]) +
         (xyz[1][1]-result[1])*(xyz[1][1]-result[1]) +
         (xyz[1][2]-result[2])*(xyz[1][2]-result[2]);
  if (dist > mindist) mindist = dist;
  mindist = sqrt(mindist);
  if (0.1*params[1] > mindist) mindist = 0.1*params[1];
#ifdef DEBUG
  printf("%lX     minDist = %le\n", tID, mindist);
#endif
  
  /* periodic -- add a vertex */
  if (mtype == ONENODE) {
    xyz[2][0] = xyz[1][0];
    xyz[2][1] = xyz[1][1];
    xyz[2][2] = xyz[1][2];
    aux[2][0] = aux[1][0];
    aux[2][1] = aux[1][1];
    aux[2][2] = aux[1][2];
    t[2]      = t[1];
    t[1]      = 0.5*(t[0]+t[2]);
    stat      = EG_evaluate(edge, &t[1], result);
    if (stat != EGADS_SUCCESS) return stat;
    dist      = sqrt(result[3]*result[3] + result[4]*result[4] +
                     result[5]*result[5]);
    if (dist == 0.0) dist = 1.0;
    xyz[1][0] = result[0];
    xyz[1][1] = result[1];
    xyz[1][2] = result[2];
    aux[1][0] = result[3]/dist;
    aux[1][1] = result[4]/dist;
    aux[1][2] = result[5]/dist;
    npts      = 3;
  }
  
  /* non-linear curve types */
  if (geom->mtype != LINE) {
    
    /* angle criteria - aux is normalized tangent */
    
    if (params[2] != 0.0) {
      stat = EG_evaluate(edge, &t[0], result);
      if (stat != EGADS_SUCCESS) return stat;
      dist = sqrt(result[3]*result[3] + result[4]*result[4] +
                  result[5]*result[5]);
      if (dist == 0.0) dist = 1.0;
      aux[0][0] = result[3]/dist;
      aux[0][1] = result[4]/dist;
      aux[0][2] = result[5]/dist;
      stat = EG_evaluate(edge, &t[npts-1], result);
      if (stat != EGADS_SUCCESS) return stat;
      dist = sqrt(result[3]*result[3] + result[4]*result[4] +
                  result[5]*result[5]);
      if (dist == 0) dist = 1.0;
      aux[npts-1][0] = result[3]/dist;
      aux[npts-1][1] = result[4]/dist;
      aux[npts-1][2] = result[5]/dist;
      
      while (npts < MAXELEN) {
        /* find the segment with the largest angle */
        k   = -1;
        dot =  1.0;
        for (i = 0; i < npts-1; i++) {
          dist = (xyz[i][0]-xyz[i+1][0])*(xyz[i][0]-xyz[i+1][0]) +
          (xyz[i][1]-xyz[i+1][1])*(xyz[i][1]-xyz[i+1][1]) +
          (xyz[i][2]-xyz[i+1][2])*(xyz[i][2]-xyz[i+1][2]);
          if (dist < mindist*mindist) continue;
          d = aux[i][0]*aux[i+1][0] + aux[i][1]*aux[i+1][1] +
          aux[i][2]*aux[i+1][2];
          if (d < dot) {
            dot = d;
            k   = i;
          }
        }
        if ((dot > dotnrm) || (k == -1)) break;
        /* insert */
        for (i = npts-1; i > k; i--) {
          xyz[i+1][0] = xyz[i][0];
          xyz[i+1][1] = xyz[i][1];
          xyz[i+1][2] = xyz[i][2];
          aux[i+1][0] = aux[i][0];
          aux[i+1][1] = aux[i][1];
          aux[i+1][2] = aux[i][2];
          t[i+1]      = t[i];
        }
        t[k+1] = 0.5*(t[k]+t[k+2]);
        stat   = EG_evaluate(edge, &t[k+1], result);
        if (stat != EGADS_SUCCESS) return stat;
        dist   = sqrt(result[3]*result[3] + result[4]*result[4] +
                      result[5]*result[5]);
        if (dist == 0.0) dist = 1.0;
        xyz[k+1][0] = result[0];
        xyz[k+1][1] = result[1];
        xyz[k+1][2] = result[2];
        aux[k+1][0] = result[3]/dist;
        aux[k+1][1] = result[4]/dist;
        aux[k+1][2] = result[5]/dist;
        npts++;
      }
#ifdef DEBUG
      printf("%lX     Angle  Phase npts = %d @ %lf (%lf)\n",
             tID, npts, dot, dotnrm);
#endif
    }
    
    /* sag - aux is midpoint value */
    if (params[1] > 0.0) {
      for (i = 0; i < npts-1; i++) {
        d    = 0.5*(t[i]+t[i+1]);
        stat = EG_evaluate(edge, &d, result);
        if (stat != EGADS_SUCCESS) return stat;
        aux[i][0] = result[0];
        aux[i][1] = result[1];
        aux[i][2] = result[2];
      }
      while (npts < MAXELEN) {
        /* find the biggest deviation */
        k    = -1;
        dist = 0.0;
        for (i = 0; i < npts-1; i++) {
          dot = (xyz[i][0]-xyz[i+1][0])*(xyz[i][0]-xyz[i+1][0]) +
          (xyz[i][1]-xyz[i+1][1])*(xyz[i][1]-xyz[i+1][1]) +
          (xyz[i][2]-xyz[i+1][2])*(xyz[i][2]-xyz[i+1][2]);
          if (dot < mindist*mindist) continue;
          mid[0] = 0.5*(xyz[i][0] + xyz[i+1][0]);
          mid[1] = 0.5*(xyz[i][1] + xyz[i+1][1]);
          mid[2] = 0.5*(xyz[i][2] + xyz[i+1][2]);
          d      = (aux[i][0]-mid[0])*(aux[i][0]-mid[0]) +
          (aux[i][1]-mid[1])*(aux[i][1]-mid[1]) +
          (aux[i][2]-mid[2])*(aux[i][2]-mid[2]);
          if (d > dist) {
            dist = d;
            k    = i;
          }
        }
        if ((dist < params[1]*params[1]) || (k == -1)) break;
        /* insert */
        for (i = npts-1; i > k; i--) {
          xyz[i+1][0] = xyz[i][0];
          xyz[i+1][1] = xyz[i][1];
          xyz[i+1][2] = xyz[i][2];
          aux[i+1][0] = aux[i][0];
          aux[i+1][1] = aux[i][1];
          aux[i+1][2] = aux[i][2];
          t[i+1]      = t[i];
        }
        t[k+1]      = 0.5*(t[k]+t[k+2]);
        xyz[k+1][0] = aux[k][0];
        xyz[k+1][1] = aux[k][1];
        xyz[k+1][2] = aux[k][2];
        d    = 0.5*(t[k+1]+t[k+2]);
        stat = EG_evaluate(edge, &d, result);
        if (stat != EGADS_SUCCESS) return stat;
        aux[k+1][0] = result[0];
        aux[k+1][1] = result[1];
        aux[k+1][2] = result[2];
        d    = 0.5*(t[k]+t[k+1]);
        stat = EG_evaluate(edge, &d, result);
        if (stat != EGADS_SUCCESS) return stat;
        aux[k][0] = result[0];
        aux[k][1] = result[1];
        aux[k][2] = result[2];
        npts++;
      }
#ifdef DEBUG
      printf("%lX     Sag    Phase npts = %d @ %lf (%lf)\n",
             tID, npts, sqrt(dist), params[1]);
#endif
    }
  }
  
  /* look at non-planar faces for curvature -- aux is uv */
  
  if (params[2] > 0.0)
    for (n = 0; n < 2; n++) {
      sense =  1;
      if (n == 0) sense = -1;
      for (nf = 0; nf < btess->tess1d[j].faces[n].nface; nf++) {
        face = btess->tess1d[j].faces[n].index;
        if (btess->tess1d[j].faces[n].nface > 1)
          face = btess->tess1d[j].faces[n].faces[nf];
        if (face <= 0) continue;
        stat = EG_getTopology(faces[face-1], &ref, &oclass, &ntype, range,
                              &ndum, &dum, &senses);
        if (stat != EGADS_SUCCESS) continue;
        if (ref == NULL) continue;
        if (ref->mtype == PLANE) continue;
        stat = EG_getTolerance(faces[face-1], &tol);
        if (stat != EGADS_SUCCESS) continue;
        if (params[1] > tol) tol = params[1];
        result[0] = result[2] =  1.e10;
        result[1] = result[3] = -1.e10;
        for (i = 0; i < npts; i++) {
          aux[i][2] = 1.0;
          stat = EG_getEdgeUV(faces[face-1], edge, sense, t[i],
                              aux[i]);
          if (stat != EGADS_SUCCESS) {
            aux[i][2] = 0.0;
          } else {
            if (aux[i][0] < result[0]) result[0] = aux[i][0];
            if (aux[i][0] > result[1]) result[1] = aux[i][0];
            if (aux[i][1] < result[2]) result[2] = aux[i][1];
            if (aux[i][1] > result[3]) result[3] = aux[i][1];
          }
        }
        /* an iso-cline on a BSPline surface with multiplicity of Knots? */
        if ((fabs(result[1]-result[0]) < UVTOL) ||
            (fabs(result[3]-result[2]) < UVTOL)) {
          stat = EG_getGeometry(ref, &oclass, &btype, &rref, &info, &prv);
          if (stat == EGADS_SUCCESS) {
            if (btype == BSPLINE) {
              double *uKnots, *vKnots;
              uKnots =  prv;
              vKnots = &prv[info[3]];
              if (fabs(result[1]-result[0]) < UVTOL) {
                if ((fabs(result[0]-uKnots[0])         > UVTOL) &&
                    (fabs(result[0]-uKnots[info[3]-1]) > UVTOL)) {
                  for (k = i = 0; i < info[3]; i++)
                    if (fabs(result[0]-uKnots[i]) < UVTOL) k++;
                  if ((outLevel > 1) && (k > 0))
                    printf("%lX Edge %d/Face %d: UisoCline %lf  multip = %d!\n",
                           tID, j+1, face, result[0], k);
                  if (k > 2) {
                    EG_free(info);
                    EG_free(prv);
                    continue;
                  }
                }
              }
              if (fabs(result[3]-result[2]) < UVTOL) {
                if ((fabs(result[2]-vKnots[0])         > UVTOL) &&
                    (fabs(result[2]-vKnots[info[6]-1]) > UVTOL)) {
                  for (k = i = 0; i < info[6]; i++)
                    if (fabs(result[2]-vKnots[i]) < UVTOL) k++;
                  if ((outLevel > 1) && (k > 0))
                    printf("%lX Edge %d/Face %d: VisoCline %lf  multip = %d!\n",
                           tID, j+1, face, result[2], k);
                  if (k > 2) {
                    EG_free(info);
                    EG_free(prv);
                    continue;
                  }
                }
              }
            }
            EG_free(info);
            EG_free(prv);
          }
        }
      
        for (i = 0; i < npts-1; i++) {
          if (aux[i][2]   <= 0.0) continue;
          if (aux[i+1][2] == 0.0) continue;
          dx[0] = xyz[i+1][0] - xyz[i][0];
          dx[1] = xyz[i+1][1] - xyz[i][1];
          dx[2] = xyz[i+1][2] - xyz[i][2];
          d     = DOT(dx, dx);
          if (d < tol*tol) {
            aux[i][2] = -1.0;
            continue;
          }
          /* get normal at mid-point in UV */
          stat = EG_getEdgeUV(faces[face-1], edge, sense,
                              0.5*(t[i]+t[i+1]), uv);
          if (stat != EGADS_SUCCESS) {
            aux[i][2] = 0.0;
          } else {
            dot = EG_curvNorm(faces[face-1], uv, sense*ntype, d, dx);
            if ((dot > dotnrm) || (dot < -1.1)) aux[i][2] = -1.0;
            stat = EG_otherFaces(btess->tess1d[j], faces, face, edge,
                                 t[i], 0.5*(t[i]+t[i+1]), t[i+1]);
            if (stat != EGADS_SUCCESS) aux[i][2] = -1.0;
          }
        }
        
        while (npts < MAXELEN) {
          /* find the largest segment with Face curvature too big */
          k    = -1;
          dist =  tol*tol;
          for (i = 0; i < npts-1; i++) {
            if (aux[i][2]   <= 0.0) continue;
            if (aux[i+1][2] == 0.0) continue;
            dx[0] = xyz[i+1][0] - xyz[i][0];
            dx[1] = xyz[i+1][1] - xyz[i][1];
            dx[2] = xyz[i+1][2] - xyz[i][2];
            d     = DOT(dx, dx);
            if (d < tol*tol) {
              aux[i][2] = -1.0;
              continue;
            }
            if (d < dist) continue;
            dist = d;
            k    = i;
          }
          if (k == -1) break;
          
          /* insert */
          for (i = npts-1; i > k; i--) {
            xyz[i+1][0] = xyz[i][0];
            xyz[i+1][1] = xyz[i][1];
            xyz[i+1][2] = xyz[i][2];
            aux[i+1][0] = aux[i][0];
            aux[i+1][1] = aux[i][1];
            aux[i+1][2] = aux[i][2];
            t[i+1]      = t[i];
          }
          t[k+1] = 0.5*(t[k]+t[k+2]);
          stat   = EG_evaluate(edge, &t[k+1], result);
          if (stat != EGADS_SUCCESS) return stat;
          xyz[k+1][0] = result[0];
          xyz[k+1][1] = result[1];
          xyz[k+1][2] = result[2];
          aux[k+1][2] = 1.0;
          stat = EG_getEdgeUV(faces[face-1], edge, sense, t[k+1],
                              aux[k+1]);
          if (stat != EGADS_SUCCESS) aux[k+1][2] = 0.0;
          stat = EG_getEdgeUV(faces[face-1], edge, sense,
                              0.5*(t[k]+t[k+1]), uv);
          if (stat != EGADS_SUCCESS) {
            aux[k][2] = 0.0;
          } else {
            dx[0] = xyz[k+1][0] - xyz[k][0];
            dx[1] = xyz[k+1][1] - xyz[k][1];
            dx[2] = xyz[k+1][2] - xyz[k][2];
            d     = DOT(dx, dx);
            dot   = EG_curvNorm(faces[face-1], uv,sense*ntype, d, dx);
            if ((dot > dotnrm) || (dot < -1.1)) aux[k][2] = -1.0;
            stat  = EG_otherFaces(btess->tess1d[j], faces, face, edge,
                                  t[k], 0.5*(t[k]+t[k+1]), t[k+1]);
            if (stat != EGADS_SUCCESS) aux[k][2] = -1.0;
          }
          stat = EG_getEdgeUV(faces[face-1], edge, sense,
                              0.5*(t[k+1]+t[k+2]), uv);
          if (stat != EGADS_SUCCESS) {
            aux[k+1][2] = 0.0;
          } else {
            dx[0] = xyz[k+2][0] - xyz[k+1][0];
            dx[1] = xyz[k+2][1] - xyz[k+1][1];
            dx[2] = xyz[k+2][2] - xyz[k+1][2];
            d     = DOT(dx, dx);
            dot   = EG_curvNorm(faces[face-1], uv, sense*ntype, d, dx);
            if ((dot > dotnrm) || (dot < -1.1)) aux[k+1][2] = -1.0;
            stat  = EG_otherFaces(btess->tess1d[j], faces, face, edge,
                                  t[k+1], 0.5*(t[k+1]+t[k+2]), t[k+2]);
            if (stat != EGADS_SUCCESS) aux[k+1][2] = -1.0;
          }
          npts++;
        }
#ifdef DEBUG
        printf("%lX     FacNrm Phase npts = %d @ %lf  Face = %d\n",
               tID, npts, dotnrm, face);
#endif
      }
    }
  
#ifdef SPLITDEGEN
  /* split edge segments that approach degenerate Nodes */
  
  be[0] = be[1] = 0;
  for (n = 0; n < 2; n++) {
    sense =  1;
    if (n == 0) sense = -1;
    for (nf = 0; nf < btess->tess1d[j].faces[n].nface; nf++) {
      face = btess->tess1d[j].faces[n].index;
      if (btess->tess1d[j].faces[n].nface > 1)
        face = btess->tess1d[j].faces[n].faces[nf];
      if (face <= 0) continue;
      stat = EG_getTopology(faces[face-1], &ref, &oclass, &ntype, range,
                            &ndum, &dum, &senses);
      if (stat != EGADS_SUCCESS) continue;
      if (ref == NULL) continue;
      if (ref->mtype == PLANE) continue;
      
      /* look at beginning of Edge */
      if (be[0] < 2) {
        stat = EG_getEdgeUV(faces[face-1], edge, sense, t[0], aux[0]);
        if (stat == EGADS_SUCCESS) {
          stat = EG_evaluate(faces[face-1], aux[0], result);
          if (stat == EGADS_SUCCESS) {
            if ((sqrt(result[3]*result[3] + result[4]*result[4] +
                      result[5]*result[5]) < DEGENUV) ||
                (sqrt(result[6]*result[6] + result[7]*result[7] +
                      result[8]*result[8]) < DEGENUV)) {
/*            printf("%lX Face #%d: Edge %d, sense %d beg is degenerate!\n",
                     tID, face, j+1, sense);  */
              be[0]++;
              if (npts < MAXELEN) {
                d    = 0.5*(t[0]+t[1]);
                stat = EG_evaluate(edge, &d, result);
                if (stat == EGADS_SUCCESS) {
                  /* insert */
                  for (i = npts-1; i > 0; i--) {
                    xyz[i+1][0] = xyz[i][0];
                    xyz[i+1][1] = xyz[i][1];
                    xyz[i+1][2] = xyz[i][2];
                    t[i+1]      = t[i];
                  }
                  xyz[1][0] = result[0];
                  xyz[1][1] = result[1];
                  xyz[1][2] = result[2];
                  t[1]      = d;
                  npts++;
                }
              }
            }
          }
        }
      }
      
      /* end of Edge */
      if (be[1] >= 2) continue;
      stat = EG_getEdgeUV(faces[face-1], edge, sense, t[npts-1], aux[1]);
      if (stat != EGADS_SUCCESS) continue;
      stat = EG_evaluate(faces[face-1], aux[1], result);
      if (stat != EGADS_SUCCESS) continue;
      if ((sqrt(result[3]*result[3] + result[4]*result[4] +
                result[5]*result[5]) < DEGENUV) ||
          (sqrt(result[6]*result[6] + result[7]*result[7] +
                result[8]*result[8]) < DEGENUV)) {
/*      printf("%lX Face #%d: Edge %d, sense %d end is degenerate!\n",
               tID, face, j+1, sense);  */
        be[1]++;
        if (npts >= MAXELEN) continue;
        /* insert */
        d    = 0.5*(t[npts-2]+t[npts-1]);
        stat = EG_evaluate(edge, &d, result);
        if (stat != EGADS_SUCCESS) continue;
        for (i = npts-1; i > npts-2; i--) {
          xyz[i+1][0] = xyz[i][0];
          xyz[i+1][1] = xyz[i][1];
          xyz[i+1][2] = xyz[i][2];
          t[i+1]      = t[i];
        }
        xyz[npts-1][0] = result[0];
        xyz[npts-1][1] = result[1];
        xyz[npts-1][2] = result[2];
        t[npts-1]      = d;
        npts++;
      }
    }
    if ((be[0] >= 2) && (be[1] >= 2)) break;
  }
#endif
  
  /* max side -- for all curve types */
  
  if (params[0] > 0.0) {
    for (i = 0; i < npts-1; i++)
      aux[i][0] = (xyz[i][0]-xyz[i+1][0])*(xyz[i][0]-xyz[i+1][0]) +
      (xyz[i][1]-xyz[i+1][1])*(xyz[i][1]-xyz[i+1][1]) +
      (xyz[i][2]-xyz[i+1][2])*(xyz[i][2]-xyz[i+1][2]);
    aux[npts-1][0] = 0.0;
    /* check face uvs
    for (i = 0; i < npts-2; i++) {
      stat = EG_otherFaces(btess->tess1d[j], faces, 0, edge,
                           t[i], t[i+1], t[i+2]);
      if (stat == EGADS_SUCCESS) continue;
      printf("%lX *** fold in Edge %d   %d/%d ***\n", tID, j+1, i+1, npts);
    }  */
    while (npts < MAXELEN) {
      /* find the biggest segment */
      k    = 0;
      dist = aux[0][0];
      for (i = 1; i < npts-1; i++) {
        d = aux[i][0];
        if (d > dist) {
          dist = d;
          k    = i;
        }
      }
      if (dist < params[0]*params[0]) break;
/*    stat = EG_otherFaces(btess->tess1d[j], faces, 0, edge,
                           t[k], 0.5*(t[k]+t[k+1]), t[k+1]);
      if (stat != EGADS_SUCCESS) {
        printf("%lX *** avoid fold in Edge %d  %d/%d ***\n", tID, j+1, k+1, npts);
        aux[k][0] = 0.0;
        continue;
      }  */
      /* insert */
      for (i = npts-1; i > k; i--) {
        xyz[i+1][0] = xyz[i][0];
        xyz[i+1][1] = xyz[i][1];
        xyz[i+1][2] = xyz[i][2];
        aux[i+1][0] = aux[i][0];
        t[i+1]      = t[i];
      }
      t[k+1] = 0.5*(t[k]+t[k+2]);
      stat   = EG_evaluate(edge, &t[k+1], result);
      if (stat != EGADS_SUCCESS) return stat;
      xyz[k+1][0] = result[0];
      xyz[k+1][1] = result[1];
      xyz[k+1][2] = result[2];
      npts++;
      d = (xyz[k][0]-xyz[k+1][0])*(xyz[k][0]-xyz[k+1][0]) +
      (xyz[k][1]-xyz[k+1][1])*(xyz[k][1]-xyz[k+1][1]) +
      (xyz[k][2]-xyz[k+1][2])*(xyz[k][2]-xyz[k+1][2]);
      aux[k][0] = d;
      if (d < 0.0625*params[0]*params[0]) break;
      d = (xyz[k+2][0]-xyz[k+1][0])*(xyz[k+2][0]-xyz[k+1][0]) +
      (xyz[k+2][1]-xyz[k+1][1])*(xyz[k+2][1]-xyz[k+1][1]) +
      (xyz[k+2][2]-xyz[k+1][2])*(xyz[k+2][2]-xyz[k+1][2]);
      aux[k+1][0] = d;
      if (d < 0.0625*params[0]*params[0]) break;
    }
  }
#ifdef DEBUG
  if (params[0] > 0.0)
    printf("%lX     MxSide Phase npts = %d @ %lf (%lf)\n",
           tID, npts, sqrt(dist), params[0]);
#endif
  
  /* fill in the 1D structure */
fill1D:
  btess->tess1d[j].xyz = (double *) EG_alloc(3*npts*sizeof(double));
  if (btess->tess1d[j].xyz == NULL) {
    if (outLevel > 0)
      printf("%lX EGADS Error: Alloc %d Pts Edge %d (EG_tessEdges)!\n",
             tID, npts, j+1);
    return EGADS_MALLOC;
  }
  btess->tess1d[j].t = (double *) EG_alloc(npts*sizeof(double));
  if (btess->tess1d[j].t == NULL) {
    EG_free(btess->tess1d[j].xyz);
    btess->tess1d[j].xyz = NULL;
    if (outLevel > 0)
      printf("%lX EGADS Error: Alloc %d Ts Edge %d (EG_tessEdges)!\n",
             tID, npts, j+1);
    return EGADS_MALLOC;
  }
  nf = btess->tess1d[j].faces[0].nface;
  if (nf > 0) {
    btess->tess1d[j].faces[0].tric = (int *) EG_alloc((nf*(npts-1))*sizeof(int));
    if (btess->tess1d[j].faces[0].tric == NULL) {
      EG_free(btess->tess1d[j].t);
      btess->tess1d[j].t   = NULL;
      EG_free(btess->tess1d[j].xyz);
      btess->tess1d[j].xyz = NULL;
      if (outLevel > 0)
        printf("%lX EGADS Error: Alloc %d Tric- Edge %d (EG_tessEdges)!\n",
               tID, npts, j+1);
      return EGADS_MALLOC;
    }
  }
  nf = btess->tess1d[j].faces[1].nface;
  if (nf > 0) {
    btess->tess1d[j].faces[1].tric = (int *) EG_alloc((nf*(npts-1))*sizeof(int));
    if (btess->tess1d[j].faces[1].tric == NULL) {
      if (btess->tess1d[j].faces[0].tric != NULL)
        EG_free(btess->tess1d[j].faces[0].tric);
      btess->tess1d[j].faces[0].tric = NULL;
      EG_free(btess->tess1d[j].t);
      btess->tess1d[j].t   = NULL;
      EG_free(btess->tess1d[j].xyz);
      btess->tess1d[j].xyz = NULL;
      if (outLevel > 0)
        printf("%lX EGADS Error: Alloc %d Tric+ Edge %d (EG_tessEdges)!\n",
               tID, npts, j+1);
      return EGADS_MALLOC;
    }
  }
  for (i = 0; i < npts; i++) {
    btess->tess1d[j].xyz[3*i  ] = xyz[i][0];
    btess->tess1d[j].xyz[3*i+1] = xyz[i][1];
    btess->tess1d[j].xyz[3*i+2] = xyz[i][2];
    btess->tess1d[j].t[i]       = t[i];
  }
  for (i = 0; i < npts-1; i++) {
    nf = btess->tess1d[j].faces[0].nface;
    for (k = 0; k < nf; k++)
      btess->tess1d[j].faces[0].tric[i*nf+k] = 0;
    nf = btess->tess1d[j].faces[1].nface;
    for (k = 0; k < nf; k++)
      btess->tess1d[j].faces[1].tric[i*nf+k] = 0;
  }
  btess->tess1d[j].npts = npts;
  
  return EGADS_SUCCESS;
}


static void
EG_edgeThread(void *struc)
{
  int     index, stat;
#ifdef PROGRESS
  int     outLevel;
#endif
  long    ID;
  EMPtess *tthread;
  
  tthread = (EMPtess *) struc;
#ifdef PROGRESS
  outLevel = EG_outLevel(tthread->body);
#endif
  
  /* get our identifier */
  ID = EMP_ThreadID();
  
  /* look for work */
  for (;;) {
    
    /* only one thread at a time here -- controlled by a mutex! */
    if (tthread->mutex != NULL) EMP_LockSet(tthread->mutex);
    if (tthread->mark == NULL) {
      index = tthread->index;
    } else {
      for (index = tthread->index; index < tthread->end; index++) {
        if (tthread->mark[index] == 0) continue;
        break;
      }
    }
    tthread->index = index+1;
    if (tthread->mutex != NULL) EMP_LockRelease(tthread->mutex);
    if (index >= tthread->end) break;
#ifdef PROGRESS
    if (outLevel > 0) {
      printf("    tessellating Edge %3d of %3d\r", index+1, tthread->end);
      fflush(stdout);
    }
#endif

    /* do the work */
    stat = EG_tessEdge(tthread->btess, tthread->faces, index,
                       tthread->edges[index], ID);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: Edge %d -> EG_tessEdge = %d (EG_edgeThread)!\n",
             index+1, stat);
  }
  
  /* exhausted all work -- exit */
  if (ID != tthread->master) EMP_ThreadExit();
}


static int
EG_tessEdges(egTessel *btess, /*@null@*/ int *retess)
{
  int      i, j, k, n, stat, outLevel, nedge, oclass, mtype, np;
  int      nface, nloop, ndum, *senses, *finds, *lsense, lor;
  double   limits[4];
  long     start;
  void     **threads = NULL;
  egObject *body, *geom, **faces, **loops, **edges, **dum;
  EMPtess  tthread;

  body     = btess->src;
  outLevel = EG_outLevel(body);
  
  stat = EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getBodyTopos(body, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) return stat;
  
  if (retess == NULL) {
    btess->tess1d = (egTess1D *) EG_alloc(nedge*sizeof(egTess1D));
    if (btess->tess1d == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Alloc %d Edges (EG_tessEdges)!\n", nedge);
      EG_free(faces);
      EG_free(edges);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nedge; j++) {
      btess->tess1d[j].obj            = edges[j];
      btess->tess1d[j].faces[0].index = 0;
      btess->tess1d[j].faces[0].nface = 0;
      btess->tess1d[j].faces[0].faces = NULL;
      btess->tess1d[j].faces[0].tric  = NULL;
      btess->tess1d[j].faces[1].index = 0;
      btess->tess1d[j].faces[1].nface = 0;
      btess->tess1d[j].faces[1].faces = NULL;
      btess->tess1d[j].faces[1].tric  = NULL;
      btess->tess1d[j].nodes[0]       = 0;
      btess->tess1d[j].nodes[1]       = 0;
      btess->tess1d[j].xyz            = NULL;
      btess->tess1d[j].t              = NULL;
      btess->tess1d[j].global         = NULL;
      btess->tess1d[j].npts           = 0;
    }
    btess->nEdge = nedge;
  
    /* get the face indices (if any) */
    for (i = 0; i < nface; i++) {
      stat = EG_getTopology(faces[i], &geom, &oclass, &mtype, limits,
                            &nloop, &loops, &lsense);
      if (stat != EGADS_SUCCESS) continue;
      for (j = 0; j < nloop; j++) {
        lor = 1;
        if ((lsense[j] == 2) || (lsense[j] == -2)) lor = -1;
        stat = EG_getTopology(loops[j], &geom, &oclass, &mtype, limits,
                              &ndum, &dum, &senses);
        if (stat != EGADS_SUCCESS) continue;
        for (k = 0; k < ndum; k++) {
          n = EG_indexBodyTopo(body, dum[k]);
          if (n <= EGADS_SUCCESS) continue;
          if (senses[k]*lor < 0) {
            if (btess->tess1d[n-1].faces[0].nface != 0) {
              if (btess->tess1d[n-1].faces[0].nface == 1) {
                btess->tess1d[n-1].faces[0].faces = (int *) EG_alloc(2*sizeof(int));
                if (btess->tess1d[n-1].faces[0].faces == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: Alloc (-) Edge %d (EG_tessEdges)!\n", n);
                  EG_free(faces);
                  EG_free(edges);
                  return EGADS_MALLOC;
                }
                btess->tess1d[n-1].faces[0].faces[0] = btess->tess1d[n-1].faces[0].index;
                btess->tess1d[n-1].faces[0].faces[1] = i+1;
              } else {
                finds = (int *) EG_reall( btess->tess1d[n-1].faces[0].faces,
                                         (btess->tess1d[n-1].faces[0].nface+1)*sizeof(int));
                if (finds == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: ReAlloc (-) Edge %d (EG_tessEdges)!\n", n);
                  EG_free(faces);
                  EG_free(edges);
                  return EGADS_MALLOC;
                } 
                finds[btess->tess1d[n-1].faces[0].nface] = i+1;
                btess->tess1d[n-1].faces[0].faces = finds;
              }
            }
            btess->tess1d[n-1].faces[0].index = i+1;
            btess->tess1d[n-1].faces[0].nface++;
          } else {
            if (btess->tess1d[n-1].faces[1].nface != 0) {
              if (btess->tess1d[n-1].faces[1].nface == 1) {
                btess->tess1d[n-1].faces[1].faces = (int *) EG_alloc(2*sizeof(int));
                if (btess->tess1d[n-1].faces[1].faces == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: Alloc (+) Edge %d (EG_tessEdges)!\n", n);
                  EG_free(faces);
                  EG_free(edges);
                  return EGADS_MALLOC;
                }
                btess->tess1d[n-1].faces[1].faces[0] = btess->tess1d[n-1].faces[1].index;
                btess->tess1d[n-1].faces[1].faces[1] = i+1;
              } else {
                finds = (int *) EG_reall( btess->tess1d[n-1].faces[1].faces,
                                         (btess->tess1d[n-1].faces[1].nface+1)*sizeof(int));
                if (finds == NULL) {
                  if (outLevel > 0)
                    printf(" EGADS Error: ReAlloc (+) Edge %d (EG_tessEdges)!\n", n);
                  EG_free(faces);
                  EG_free(edges);
                  return EGADS_MALLOC;
                }
                finds[btess->tess1d[n-1].faces[1].nface] = i+1;
                btess->tess1d[n-1].faces[1].faces = finds;
              }
            }
            btess->tess1d[n-1].faces[1].index = i+1;
            btess->tess1d[n-1].faces[1].nface++;
          }
        }
      }
    }
    /* report any non-manifold Edges */
    if (outLevel > 1)
      for (j = 0; j < nedge; j++) {
        if (btess->tess1d[j].faces[0].nface > 1) {
          printf(" EGADS Internal: Non-manifold Edge %d (-) with Faces", j+1);
          if (btess->tess1d[j].faces[0].faces != NULL)
            for (k = 0; k < btess->tess1d[j].faces[0].nface; k++)
              printf(" %d", btess->tess1d[j].faces[0].faces[k]);
          printf("!\n");
        }
        if (btess->tess1d[j].faces[1].nface > 1) {
          printf(" EGADS Internal: Non-manifold Edge %d (+) with Faces", j+1);
          if (btess->tess1d[j].faces[1].faces != NULL)
            for (k = 0; k < btess->tess1d[j].faces[1].nface; k++)
              printf(" %d", btess->tess1d[j].faces[1].faces[k]);
          printf("!\n");
        }
      }
  }

  /* set up for explicit multithreading */
  tthread.mutex     = NULL;
  tthread.master    = EMP_ThreadID();
  tthread.index     = 0;
  tthread.end       = nedge;
  tthread.mark      = retess;
  tthread.tess      = NULL;
  tthread.btess     = btess;
  tthread.body      = body;
  tthread.faces     = faces;
  tthread.edges     = edges;
  tthread.params    = NULL;
  tthread.tparam    = NULL;
  tthread.qparam[0] = tthread.qparam[1] = tthread.qparam[2] = 0.0;
  
  np = EMP_Init(&start);
  if (outLevel > 1) printf(" EMP NumProcs = %d!\n", np);
  
  if (np > 1) {
    /* create the mutex to handle list synchronization */
    tthread.mutex = EMP_LockCreate();
    if (tthread.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(tthread.mutex);
        np = 1;
      }
    }
  }
  
  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EG_edgeThread, &tthread);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  EG_edgeThread(&tthread);
  
  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

#ifdef PROGRESS
  if (outLevel > 0) printf("\n");
#endif
  
  /* cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  if (tthread.mutex != NULL) EMP_LockDestroy(tthread.mutex);
  if (threads != NULL) free(threads);
  if (outLevel > 1)
    printf(" EMP Number of Seconds on Edge Thread Block = %ld\n",
             EMP_Done(&start));

  EG_free(faces);
  EG_free(edges);
  return EGADS_SUCCESS;
}


int
EG_makeTessGeom(egObject *obj, double *params, int *sizes, egObject **tess)
{
  int      i, j, k, stat, outLevel, np, nu, nv = 0;
  double   *dtess, uv[2], result[18];
  egTessel *btess;
  egObject *gtess, *context;

  *tess = NULL;
  if  (obj == NULL)               return EGADS_NULLOBJ;
  if  (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((obj->oclass != SURFACE) && (obj->oclass != CURVE)) return EGADS_NOTGEOM;
  outLevel = EG_outLevel(obj);
  context  = EG_context(obj);

  nu = abs(sizes[0]);
  np = nu;
  if (obj->oclass == SURFACE) {
    nv = abs(sizes[1]);
    if ((nu < 2) || (nv < 2)) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface size = %d %d (EG_makeTessGeom)!\n",
               nu, nv);
      return EGADS_INDEXERR;
    }
    np *= nv;
  } else {
    if (nu < 2) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve len = %d (EG_makeTessGeom)!\n", nu);
      return EGADS_INDEXERR;
    }
  }
  
  btess = (egTessel *) EG_alloc(sizeof(egTessel));
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Blind Malloc (EG_makeTessGeom)!\n");
    return EGADS_MALLOC;
  }
  btess->src     = obj;
  btess->xyzs    = NULL;
  btess->tess1d  = NULL;
  btess->tess2d  = NULL;
  btess->globals = NULL;
  btess->nGlobal = 0;
  btess->nEdge   = 0;
  btess->nFace   = 0;
  btess->nu      = nu;
  btess->nv      = nv;
  btess->done    = 1;
  
  /* get the storage for the tessellation */
  dtess = (double *) EG_alloc(3*np*sizeof(double));
  if (dtess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Data Malloc (EG_makeTessGeom)!\n");
    EG_free(btess);
    return EGADS_MALLOC;
  }
  
  stat = EG_makeObject(context, &gtess);
  if (stat != EGADS_SUCCESS) {
    EG_free(dtess);
    EG_free(btess);
    return stat;
  }
  gtess->oclass = TESSELLATION;
  gtess->mtype  = obj->oclass;
  gtess->blind  = btess;
  EG_referenceObject(gtess, context);
  EG_referenceTopObj(obj,   gtess);
  
  /* fill the data */
  
  if (obj->oclass == SURFACE) {
  
    for (k = j = 0; j < nv; j++) {
      if (sizes[1] < 0) {
        uv[1] = params[2] + (nv-j-1)*(params[3]-params[2])/(nv-1);
      } else {
        uv[1] = params[2] +        j*(params[3]-params[2])/(nv-1);
      }
      for (i = 0; i < nu; i++, k++) {
        if (sizes[0] < 0) {
          uv[0] = params[0] + (nu-i-1)*(params[1]-params[0])/(nu-1);
        } else {
          uv[0] = params[0] +        i*(params[1]-params[0])/(nu-1);
        }
        stat  = EG_evaluate(obj, uv, result);
        dtess[3*k  ] = result[0];
        dtess[3*k+1] = result[1];
        dtess[3*k+2] = result[2];
        if (stat == EGADS_SUCCESS) continue;
        if (outLevel > 0)
          printf(" EGADS Warning: %d/%d, %d/%d eval ret = %d  (EG_makeTessGeom)!\n",
                 i+1, nv, j+1, nv, stat);
      }
    }

  } else {
  
    for (i = 0; i < nu; i++) {
      if (sizes[0] < 0) {
        uv[0] = params[0] + (nu-i-1)*(params[1]-params[0])/(nu-1);
      } else {
        uv[0] = params[0] +        i*(params[1]-params[0])/(nu-1);
      }
      stat  = EG_evaluate(obj, uv, result);
      dtess[3*i  ] = result[0];
      dtess[3*i+1] = result[1];
      dtess[3*i+2] = result[2];
      if (stat == EGADS_SUCCESS) continue;
      if (outLevel > 0)
        printf(" EGADS Warning: %d/%d evaluate ret = %d  (EG_makeTessGeom)!\n",
               i+1, nv, stat);
    }

  }
  
  btess->xyzs      = dtess;
  btess->params[0] = params[0];
  btess->params[1] = params[1];
  btess->params[2] = nu;
  if (nv == 0) {
    btess->params[3] = btess->params[4] = btess->params[5] = 0.0;
  } else {
    btess->params[3] = params[2];
    btess->params[4] = params[3];
    btess->params[5] = nv;
  }
  for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = 0.0;

  *tess = gtess;
  return EGADS_SUCCESS;
}


int
EG_getTessGeom(const egObject *tess, int *sizes, double **xyz)
{
  int      outLevel;
  egTessel *btess;
  egObject *obj;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_getTessGeom)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_getTessGeom)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_getTessGeom)!\n");
    return EGADS_NOTOBJ;
  }
  if ((obj->oclass != SURFACE) && (obj->oclass != CURVE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not a Curve/Surface (EG_getTessGeom)!\n");
    return EGADS_NOTGEOM;
  }

  sizes[0] = btess->nu;
  sizes[1] = btess->nv;
  *xyz     = btess->xyzs;
  return EGADS_SUCCESS;
}


static void
EG_deleteQuads(egTessel *btess, int iface)
{
  int i, j;
  
  i = btess->nFace + iface - 1;
  if (btess->tess2d[i].xyz    != NULL) EG_free(btess->tess2d[i].xyz);
  if (btess->tess2d[i].uv     != NULL) EG_free(btess->tess2d[i].uv);
  if (btess->tess2d[i].ptype  != NULL) EG_free(btess->tess2d[i].ptype);
  if (btess->tess2d[i].pindex != NULL) EG_free(btess->tess2d[i].pindex);
  for (j = 0; j < btess->tess2d[i].npatch; j++) {
    if (btess->tess2d[i].patch[j].ipts != NULL) 
      EG_free(btess->tess2d[i].patch[j].ipts);
    if (btess->tess2d[i].patch[j].bounds != NULL) 
      EG_free(btess->tess2d[i].patch[j].bounds);
  }
  EG_free(btess->tess2d[i].patch);
  btess->tess2d[i].xyz    = NULL;
  btess->tess2d[i].uv     = NULL;
  btess->tess2d[i].ptype  = NULL;
  btess->tess2d[i].pindex = NULL;
  btess->tess2d[i].npts   = 0;
  btess->tess2d[i].patch  = NULL;
  btess->tess2d[i].npatch = 0;
}


int
EG_moveEdgeVert(egObject *tess, int eIndex, int vIndex, double t)
{
  int      i, j, m, nf, stat, outLevel, nedge, nface, iface, itri, ivrt;
  int      sense;
  double   result[9], uv[2];
  egTessel *btess;
  egObject *obj, **edges, **faces;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_moveEdgeVert)!\n");  
    return EGADS_NOTFOUND;
  }
  if (btess->done != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad State (EG_moveEdgeVerts)!\n");
    return EGADS_TESSTATE;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_moveEdgeVert)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_moveEdgeVert)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_moveEdgeVert)!\n");
    return EGADS_NOTBODY;
  }
  if (obj->mtype == WIREBODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source is WireBody (EG_moveEdgeVert)!\n");
    return EGADS_TOPOERR;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_moveEdgeVert)!\n");
    return EGADS_NODATA;  
  }
  if ((eIndex < 1) || (eIndex > btess->nEdge)) {
    if (outLevel > 0)
      printf(" EGADS Error: eIndex = %d [1-%d] (EG_moveEdgeVert)!\n",
             eIndex, btess->nEdge);
    return EGADS_INDEXERR;
  }
  if ((vIndex < 2) || (eIndex >= btess->tess1d[eIndex-1].npts)) {
    if (outLevel > 0)
      printf(" EGADS Error: vIndex = %d [2-%d] (EG_moveEdgeVert)!\n",
             vIndex, btess->tess1d[eIndex-1].npts-1);
    return EGADS_INDEXERR;
  }
  if ((t <= btess->tess1d[eIndex-1].t[vIndex-2]) ||
      (t >= btess->tess1d[eIndex-1].t[vIndex])) {
    if (outLevel > 0)
      printf(" EGADS Error: t = %lf [%lf-%lf] (EG_moveEdgeVert)!\n",
             t, btess->tess1d[eIndex-1].t[vIndex-2], 
                btess->tess1d[eIndex-1].t[vIndex]);
    return EGADS_RANGERR;
  }
  stat = EG_getBodyTopos(btess->src, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getBodyTopos(btess->src, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    EG_free(edges);
    return stat;
  }
  
  stat = EG_evaluate(edges[eIndex-1], &t, result);
  if (stat != EGADS_SUCCESS) {
    EG_free(faces);
    EG_free(edges);
    return stat;
  }
  /* make sure we can get UVs */
  for (m = 0; m < 2; m++) {
    nf = btess->tess1d[eIndex-1].faces[m].nface;
    for (j = 0; j < nf; j++) {
      iface = btess->tess1d[eIndex-1].faces[m].index;
      if (nf > 1) iface = btess->tess1d[eIndex-1].faces[m].faces[j];
      if (iface != 0) {
        sense = faces[iface-1]->mtype;
        if (EG_faceConnIndex(btess->tess1d[eIndex-1].faces[1-m], iface) == 0)
          sense = 0;
        if (m == 0) sense = -sense;
        stat = EG_getEdgeUV(faces[iface-1], edges[eIndex-1], sense, t, uv);
        if (stat != EGADS_SUCCESS) {
          EG_free(faces);
          EG_free(edges);
          return stat;
        }
      }
    }
  }
  
  /* got everything -- update the tessellation */
  btess->tess1d[eIndex-1].xyz[3*vIndex-3] = result[0];
  btess->tess1d[eIndex-1].xyz[3*vIndex-2] = result[1];
  btess->tess1d[eIndex-1].xyz[3*vIndex-1] = result[2];
  btess->tess1d[eIndex-1].t[vIndex-1]     = t;
  for (m = 0; m < 2; m++) {
    nf = btess->tess1d[eIndex-1].faces[m].nface;
    for (j = 0; j < nf; j++) {
      iface = btess->tess1d[eIndex-1].faces[m].index;
      if (nf > 1) iface = btess->tess1d[eIndex-1].faces[m].faces[j];
      if (iface == 0) continue;
      sense = faces[iface-1]->mtype;
      if (EG_faceConnIndex(btess->tess1d[eIndex-1].faces[1-m], iface) == 0)
        sense = 0;
      if (m == 0) sense = -sense;
      EG_getEdgeUV(faces[iface-1], edges[eIndex-1], sense, t, uv);
      itri = btess->tess1d[eIndex-1].faces[m].tric[(vIndex-1)*nf+j] - 1;
      for (i = 0; i < 3; i++) {
        ivrt = btess->tess2d[iface-1].tris[3*itri+i] - 1;
        if ((btess->tess2d[iface-1].pindex[ivrt] == eIndex) &&
            (btess->tess2d[iface-1].ptype[ivrt]  == vIndex)) {
          btess->tess2d[iface-1].xyz[3*ivrt  ] = result[0];
          btess->tess2d[iface-1].xyz[3*ivrt+1] = result[1];
          btess->tess2d[iface-1].xyz[3*ivrt+2] = result[2];
          btess->tess2d[iface-1].uv[2*ivrt  ]  = uv[0];
          btess->tess2d[iface-1].uv[2*ivrt+1]  = uv[1];
          break;
        }
      }
      /* delete any quads & invalidate the frame */
      if (btess->tess2d[iface-1].bary  != NULL)
        EG_free(btess->tess2d[iface-1].bary);
      btess->tess2d[iface-1].bary = NULL;
      if (btess->tess2d[iface-1].frame != NULL)
        EG_free(btess->tess2d[iface-1].frame);
      btess->tess2d[iface-1].frame = NULL;
      if (btess->tess2d[iface-1].frlps != NULL)
        EG_free(btess->tess2d[iface-1].frlps);
      btess->tess2d[iface-1].frlps = NULL;
      EG_deleteQuads(btess, iface);
    }
  }
  EG_free(faces);
  EG_free(edges);
  
  return EGADS_SUCCESS;
}


void
EG_cleanupTessMaps(egTessel *btess)
{
  int i;
  
  if (btess->xyzs != NULL) {
    EG_free(btess->xyzs);
    btess->xyzs = NULL;
  }
  if (btess->tess1d != NULL)
    for (i = 0; i < btess->nEdge; i++) {
      if (btess->tess1d[i].global == NULL) continue;
      EG_free(btess->tess1d[i].global);
      btess->tess1d[i].global = NULL;
    }
  if (btess->tess2d != NULL)
    for (i = 0; i < 2*btess->nFace; i++) {
      if (btess->tess2d[i].global == NULL) continue;
      EG_free(btess->tess2d[i].global);
      btess->tess2d[i].global = NULL;
    }
  if (btess->globals != NULL) {
    EG_free(btess->globals);
    btess->globals = NULL;
    btess->nGlobal = 0;
  }
}


int
EG_deleteEdgeVert(egObject *tess, int eIndex, int vIndex, int dir)
{
  int      i, k, m, n, nf, outLevel, iface, iv[2], it, ivert;
  int      n1, n2, ie, i1, i2, i3, pt1, pi1, pt2, pi2, ref, nfr;
  egTessel *btess;
  egObject *obj;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  if ((dir != -1) && (dir != 1)) {
    if (outLevel > 0)
      printf(" EGADS Error: Collapse Dir = %d (EG_deleteEdgeVert)!\n",
             dir);  
    return EGADS_RANGERR;
  }
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_deleteEdgeVert)!\n");  
    return EGADS_NOTFOUND;
  }
  if (btess->done != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad State (EG_deleteEdgeVerts)!\n");
    return EGADS_TESSTATE;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_deleteEdgeVert)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_deleteEdgeVert)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_deleteEdgeVert)!\n");
    return EGADS_NOTBODY;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_deleteEdgeVert)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_deleteEdgeVert)!\n");
    return EGADS_NODATA;  
  }
  if ((eIndex < 1) || (eIndex > btess->nEdge)) {
    if (outLevel > 0)
      printf(" EGADS Error: eIndex = %d [1-%d] (EG_deleteEdgeVert)!\n",
             eIndex, btess->nEdge);
    return EGADS_INDEXERR;
  }
  if ((vIndex < 2) || (eIndex >= btess->tess1d[eIndex-1].npts)) {
    if (outLevel > 0)
      printf(" EGADS Error: vIndex = %d [2-%d] (EG_deleteEdgeVert)!\n",
             vIndex, btess->tess1d[eIndex-1].npts-1);
    return EGADS_INDEXERR;
  }
  
  /* cleanup mappings */
  EG_cleanupTessMaps(btess);
 
  /* fix up each face */
  for (m = 0; m < 2; m++) {
    nf = btess->tess1d[eIndex-1].faces[m].nface;
    for (n = 0; n < nf; n++) {
      iface = btess->tess1d[eIndex-1].faces[m].index;
      if (nf > 1) iface = btess->tess1d[eIndex-1].faces[m].faces[n];
      if (iface == 0) continue;
      if (dir == -1) {
        it = btess->tess1d[eIndex-1].faces[m].tric[nf*(vIndex-2)+n];
      } else {
        it = btess->tess1d[eIndex-1].faces[m].tric[nf*(vIndex-1)+n];
      }
      iv[0] = iv[1] = -1;
      /* find the vert to remove */
      for (i = 0; i < 3; i++) {
        ivert = btess->tess2d[iface-1].tris[3*it+i-3];
        if ((btess->tess2d[iface-1].pindex[ivert-1] == eIndex) &&
            (btess->tess2d[iface-1].ptype[ivert-1]  == vIndex)) {
          iv[0] = ivert;
          break;
        }
      }
      /* find the vert to collpse to */
      for (i = 0; i < 3; i++) {
        ivert = btess->tess2d[iface-1].tris[3*it+i-3];
        if ((btess->tess2d[iface-1].pindex[ivert-1] == eIndex) &&
            (btess->tess2d[iface-1].ptype[ivert-1]  == vIndex+dir)) {
          iv[1] = ivert;
          break;
        }
      }
      if ((iv[0] == -1) || (iv[1] == -1)) {
        printf(" EGADS Internal: EG_deleteEdgeVert Verts = %d %d!\n",
               iv[0], iv[1]);
        return EGADS_GEOMERR;
      }

      pt1 =       vIndex;
      pi1 = pi2 = eIndex;
      pt2 = vIndex + dir;
      if (pt2 == 1) {
        pt2 = 0;
        pi2 = btess->tess1d[eIndex-1].nodes[0];
      }
      if (pt2 == btess->tess1d[eIndex-1].npts) {
        pt2 = 0;
        pi2 = btess->tess1d[eIndex-1].nodes[1];
      }

      /* patch up the neighbors for the removed triangle */
      i1 = btess->tess2d[iface-1].tris[3*it-3]-1;
      i2 = btess->tess2d[iface-1].tris[3*it-2]-1;
      i3 = btess->tess2d[iface-1].tris[3*it-1]-1;        
      if (((btess->tess2d[iface-1].pindex[i2] == pi1) &&
           (btess->tess2d[iface-1].ptype[i2]  == pt1) &&
           (btess->tess2d[iface-1].pindex[i3] == pi2) &&
           (btess->tess2d[iface-1].ptype[i3]  == pt2)) || 
          ((btess->tess2d[iface-1].pindex[i2] == pi2) &&
           (btess->tess2d[iface-1].ptype[i2]  == pt2) &&
           (btess->tess2d[iface-1].pindex[i3] == pi1) &&
           (btess->tess2d[iface-1].ptype[i3]  == pt1))) {
        n1 = btess->tess2d[iface-1].tric[3*it-2];
        n2 = btess->tess2d[iface-1].tric[3*it-1];
      } else if (((btess->tess2d[iface-1].pindex[i1] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt1) &&
                  (btess->tess2d[iface-1].pindex[i3] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i3]  == pt2)) || 
                 ((btess->tess2d[iface-1].pindex[i1] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt2) &&
                  (btess->tess2d[iface-1].pindex[i3] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i3]  == pt1))) {
        n1 = btess->tess2d[iface-1].tric[3*it-3];
        n2 = btess->tess2d[iface-1].tric[3*it-1];
      } else if (((btess->tess2d[iface-1].pindex[i1] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt1) &&
                  (btess->tess2d[iface-1].pindex[i2] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i2]  == pt2)) || 
                 ((btess->tess2d[iface-1].pindex[i1] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt2) &&
                  (btess->tess2d[iface-1].pindex[i2] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i2]  == pt1))) {
        n1 = btess->tess2d[iface-1].tric[3*it-3];
        n2 = btess->tess2d[iface-1].tric[3*it-2];
      } else {
        printf(" EGADS Internal: Can not find segment for %d %d  %d %d - %d!\n",
               pt1, pi1, pt2, pi2, btess->tess1d[eIndex-1].npts);
        return EGADS_GEOMERR;
      }
#ifdef DEBUG
      printf(" verts = %d (%d %d) -> %d (%d %d)   tri = %d   ns = %d %d,  dir = %d\n", 
             iv[0], pt1, pi1, iv[1], pt2, pi2, it, n1, n2, dir);
      printf("      tri = %d   %d %d %d   %d %d %d\n", it,
             btess->tess2d[iface-1].tris[3*it-3], btess->tess2d[iface-1].tris[3*it-2],
             btess->tess2d[iface-1].tris[3*it-1], btess->tess2d[iface-1].tric[3*it-3],
             btess->tess2d[iface-1].tris[3*it-2], btess->tess2d[iface-1].tric[3*it-1]);
      printf("      tri = %d   %d %d %d   %d %d %d\n", n1,
             btess->tess2d[iface-1].tris[3*n1-3], btess->tess2d[iface-1].tris[3*n1-2],
             btess->tess2d[iface-1].tris[3*n1-1], btess->tess2d[iface-1].tric[3*n1-3],
             btess->tess2d[iface-1].tris[3*n1-2], btess->tess2d[iface-1].tric[3*n1-1]);
      printf("      tri = %d   %d %d %d   %d %d %d\n", n2,
             btess->tess2d[iface-1].tris[3*n2-3], btess->tess2d[iface-1].tris[3*n2-2],
             btess->tess2d[iface-1].tris[3*n2-1], btess->tess2d[iface-1].tric[3*n2-3],
             btess->tess2d[iface-1].tris[3*n2-2], btess->tess2d[iface-1].tric[3*n2-1]);
#endif
      if (n1 > 0) {
        for (i = 0; i < 3; i++)
          if (btess->tess2d[iface-1].tric[3*n1+i-3] == it) {
            btess->tess2d[iface-1].tric[3*n1+i-3] = n2;
            break;
          }
      } else if (n1 < 0) {
        ie  = -n1;
        ref = EG_faceConnIndex(btess->tess1d[ie-1].faces[0], iface);
        nfr = btess->tess1d[ie-1].faces[0].nface;
        if (ref != 0)
          for (k = 0; k < btess->tess1d[ie-1].npts-1; k++)
            if (btess->tess1d[ie-1].faces[0].tric[nfr*k+ref] == it)
              btess->tess1d[ie-1].faces[0].tric[nfr*k+ref] = n2;
        ref = EG_faceConnIndex(btess->tess1d[ie-1].faces[1], iface);
        nfr = btess->tess1d[ie-1].faces[1].nface;
        if (ref != 0)
          for (k = 0; k < btess->tess1d[ie-1].npts-1; k++)
            if (btess->tess1d[ie-1].faces[1].tric[nfr*k+ref] == it)
              btess->tess1d[ie-1].faces[1].tric[nfr*k+ref] = n2;
      }
      if (n2 > 0) {
        for (i = 0; i < 3; i++)
          if (btess->tess2d[iface-1].tric[3*n2+i-3] == it) {
            btess->tess2d[iface-1].tric[3*n2+i-3] = n1;
            break;
          }
      } else if (n2 < 0) {
        ie  = -n2;
        ref = EG_faceConnIndex(btess->tess1d[ie-1].faces[0], iface);
        nfr = btess->tess1d[ie-1].faces[0].nface;
        if (ref != 0)
          for (k = 0; k < btess->tess1d[ie-1].npts-1; k++)
            if (btess->tess1d[ie-1].faces[0].tric[nfr*k+ref] == it)
              btess->tess1d[ie-1].faces[0].tric[nfr*k+ref] = n1;
        ref = EG_faceConnIndex(btess->tess1d[ie-1].faces[1], iface);
        nfr = btess->tess1d[ie-1].faces[1].nface;
        if (ref != 0)
          for (k = 0; k < btess->tess1d[ie-1].npts-1; k++)
            if (btess->tess1d[ie-1].faces[1].tric[nfr*k+ref] == it)
              btess->tess1d[ie-1].faces[1].tric[nfr*k+ref] = n1;
      }
      /* collapse the vert from the triangulation by subsitution*/
      for (i = 0; i < btess->tess2d[iface-1].ntris; i++) {
        if (btess->tess2d[iface-1].tris[3*i  ] == iv[0])
          btess->tess2d[iface-1].tris[3*i  ] = iv[1];
        if (btess->tess2d[iface-1].tris[3*i+1] == iv[0])
          btess->tess2d[iface-1].tris[3*i+1] = iv[1];
        if (btess->tess2d[iface-1].tris[3*i+2] == iv[0])
          btess->tess2d[iface-1].tris[3*i+2] = iv[1];
      }
  
      /* compress the face */
      for (i = 0; i < btess->tess2d[iface-1].npts; i++)
        if ((btess->tess2d[iface-1].pindex[i] == eIndex) &&
            (btess->tess2d[iface-1].ptype[i]  >= vIndex))
          btess->tess2d[iface-1].ptype[i]--;

      for (i = 0; i < btess->tess2d[iface-1].ntris; i++) {
        if (btess->tess2d[iface-1].tris[3*i  ] > iv[0])
          btess->tess2d[iface-1].tris[3*i  ]--;
        if (btess->tess2d[iface-1].tris[3*i+1] > iv[0])
          btess->tess2d[iface-1].tris[3*i+1]--;
        if (btess->tess2d[iface-1].tris[3*i+2] > iv[0])
          btess->tess2d[iface-1].tris[3*i+2]--;
        if (btess->tess2d[iface-1].tric[3*i  ] > it)
          btess->tess2d[iface-1].tric[3*i  ]--;
        if (btess->tess2d[iface-1].tric[3*i+1] > it)
          btess->tess2d[iface-1].tric[3*i+1]--;
        if (btess->tess2d[iface-1].tric[3*i+2] > it)
          btess->tess2d[iface-1].tric[3*i+2]--;
      }
      for (ie = 0; ie < btess->nEdge; ie++) {
        nfr = btess->tess1d[ie-1].faces[0].nface;
        for (i = 0; i < nfr; i++) {
          k = btess->tess1d[ie-1].faces[0].index;
          if (nfr > 1) k = btess->tess1d[ie-1].faces[0].faces[i];
          if (iface != k) continue;
          if (btess->tess1d[ie-1].faces[0].tric == NULL) continue;
          for (k = 0; k < btess->tess1d[ie-1].npts-1; k++)
            if (btess->tess1d[ie-1].faces[0].tric[nfr*k+i] > it)
              btess->tess1d[ie-1].faces[0].tric[nfr*k+i]--;
        }
        nfr = btess->tess1d[ie-1].faces[1].nface;
        for (i = 0; i < nfr; i++) {
          k = btess->tess1d[ie-1].faces[1].index;
          if (nfr > 1) k = btess->tess1d[ie-1].faces[1].faces[i];
          if (iface != k) continue;
          if (btess->tess1d[ie-1].faces[1].tric == NULL) continue;
          for (k = 0; k < btess->tess1d[ie-1].npts-1; k++)
            if (btess->tess1d[ie-1].faces[1].tric[nfr*k+i] > it)
              btess->tess1d[ie-1].faces[1].tric[nfr*k+i]--;
        }
      }
      btess->tess2d[iface-1].npts--;
      for (i = iv[0]-1; i < btess->tess2d[iface-1].npts; i++) {
        btess->tess2d[iface-1].xyz[3*i  ] = btess->tess2d[iface-1].xyz[3*i+3];
        btess->tess2d[iface-1].xyz[3*i+1] = btess->tess2d[iface-1].xyz[3*i+4];
        btess->tess2d[iface-1].xyz[3*i+2] = btess->tess2d[iface-1].xyz[3*i+5];
        btess->tess2d[iface-1].uv[2*i  ]  = btess->tess2d[iface-1].uv[2*i+2];
        btess->tess2d[iface-1].uv[2*i+1]  = btess->tess2d[iface-1].uv[2*i+3];
        btess->tess2d[iface-1].ptype[i]   = btess->tess2d[iface-1].ptype[i+1];
        btess->tess2d[iface-1].pindex[i]  = btess->tess2d[iface-1].pindex[i+1];
      }
      btess->tess2d[iface-1].ntris--;
      for (i = it-1; i < btess->tess2d[iface-1].ntris; i++) {
        btess->tess2d[iface-1].tris[3*i  ] = btess->tess2d[iface-1].tris[3*i+3];
        btess->tess2d[iface-1].tris[3*i+1] = btess->tess2d[iface-1].tris[3*i+4];
        btess->tess2d[iface-1].tris[3*i+2] = btess->tess2d[iface-1].tris[3*i+5];
        btess->tess2d[iface-1].tric[3*i  ] = btess->tess2d[iface-1].tric[3*i+3];
        btess->tess2d[iface-1].tric[3*i+1] = btess->tess2d[iface-1].tric[3*i+4];
        btess->tess2d[iface-1].tric[3*i+2] = btess->tess2d[iface-1].tric[3*i+5];
      }

      /* remove any quads & and mark frame as destroyed */
      if (btess->tess2d[iface-1].bary  != NULL)
        EG_free(btess->tess2d[iface-1].bary);
      btess->tess2d[iface-1].bary = NULL;
      if (btess->tess2d[iface-1].frame != NULL)
        EG_free(btess->tess2d[iface-1].frame);
      btess->tess2d[iface-1].frame = NULL;
      if (btess->tess2d[iface-1].frlps != NULL)
        EG_free(btess->tess2d[iface-1].frlps);
      btess->tess2d[iface-1].frlps = NULL;
      EG_deleteQuads(btess, iface);
    }
  }
  
  /* compress the Edge storage */
                 k = vIndex-1;
  if (dir == -1) k = vIndex-2;
  btess->tess1d[eIndex-1].npts--;
  for (i = k; i < btess->tess1d[eIndex-1].npts; i++) {
    if (i != btess->tess1d[eIndex-1].npts-1)
      for (m = 0; m < 2; m++) {
        nf = btess->tess1d[eIndex-1].faces[m].nface;
        for (n = 0; n < nf; n++)
          btess->tess1d[eIndex-1].faces[m].tric[nf*i+n] = 
            btess->tess1d[eIndex-1].faces[m].tric[nf*(i+1)+n];
    }
    btess->tess1d[eIndex-1].xyz[3*i  ] = btess->tess1d[eIndex-1].xyz[3*i+3];
    btess->tess1d[eIndex-1].xyz[3*i+1] = btess->tess1d[eIndex-1].xyz[3*i+4];
    btess->tess1d[eIndex-1].xyz[3*i+2] = btess->tess1d[eIndex-1].xyz[3*i+5];
    btess->tess1d[eIndex-1].t[i]       = btess->tess1d[eIndex-1].t[i+1];
  }
  
#ifdef CHECK
  EG_checkTriangulation(btess);
#endif

  return EGADS_SUCCESS;
}


int
EG_insertEdgeVerts(egObject *tess, int eIndex, int vIndex, int npts,
                   double *t)
{
  int      i, j, k, m, nf, nx, stat, outLevel, nedge, nface, iface, itri;
  int      n0, n1, v0, v1, vert, vn, nl, nn, sense, pt1, pi1, pt2, pi2;
  int      i1, i2, i3, cnt, stripe;
  int      *etric[2], *pindex, *ptype, *tris, *tric;
  double   result[9], *vals, *xyzs, *ts, *xyz, *uv;
  egTessel *btess;
  egObject *obj, **edges, **faces;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  if (npts <= 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Inserts (EG_insertEdgeVerts)!\n");
    return EGADS_RANGERR;
  }
  for (i = 0; i < npts-1; i++)
    if (t[i+1] <= t[i]) {
      if (outLevel > 0)
        printf(" EGADS Error: Ts are NOT monitonic (EG_insertEdgeVerts)!\n");  
      return EGADS_RANGERR;
    }

  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_insertEdgeVerts)!\n");  
    return EGADS_NOTFOUND;
  }
  if (btess->done != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad State (EG_insertEdgeVerts)!\n");
    return EGADS_TESSTATE;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_insertEdgeVerts)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_insertEdgeVerts)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_insertEdgeVerts)!\n");
    return EGADS_NOTBODY;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_insertEdgeVert)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Edge Tessellations (EG_insertEdgeVerts)!\n");
    return EGADS_NODATA;  
  }
  if ((eIndex < 1) || (eIndex > btess->nEdge)) {
    if (outLevel > 0)
      printf(" EGADS Error: eIndex = %d [1-%d] (EG_insertEdgeVerts)!\n",
             eIndex, btess->nEdge);
    return EGADS_INDEXERR;
  }
  if ((vIndex < 1) || (vIndex >= btess->tess1d[eIndex-1].npts)) {
    if (outLevel > 0)
      printf(" EGADS Error: vIndex = %d [1-%d] (EG_insertEdgeVerts)!\n",
             vIndex, btess->tess1d[eIndex-1].npts-1);
    return EGADS_INDEXERR;
  }
  if ((t[0]      <= btess->tess1d[eIndex-1].t[vIndex-1]) ||
      (t[npts-1] >= btess->tess1d[eIndex-1].t[vIndex])) {
    if (outLevel > 0)
      printf(" EGADS Error: t = %lf %lf [%lf-%lf] (EG_insertEdgeVerts)!\n",
             t[0], t[npts-1], btess->tess1d[eIndex-1].t[vIndex-1], 
                              btess->tess1d[eIndex-1].t[vIndex]);
    return EGADS_RANGERR;
  }
  
  /* make sure we are not inserting along a DEGEN Edge */
  for (cnt = m = 0; m < 2; m++) {
    nf = btess->tess1d[eIndex-1].faces[m].nface;
    for (nx = 0; nx < nf; nx++) {
      iface = btess->tess1d[eIndex-1].faces[m].index;
      if (nf > 1) iface = btess->tess1d[eIndex-1].faces[m].faces[nx];
      if (iface == 0) continue;
      itri = btess->tess1d[eIndex-1].faces[m].tric[(vIndex-1)*nf+nx];
      i1   = btess->tess2d[iface-1].tris[3*itri-3]-1;
      i2   = btess->tess2d[iface-1].tris[3*itri-2]-1;
      i3   = btess->tess2d[iface-1].tris[3*itri-1]-1;
      if (btess->tess2d[iface-1].pindex[i1] == btess->tess2d[iface-1].pindex[i2])
        if ((btess->tess2d[iface-1].ptype[i1] == 0) &&
            (btess->tess2d[iface-1].ptype[i2] == 0)) {
          if (outLevel > 0) {
            printf(" EGADS Error: Degen EDGE (EG_insertEdgeVerts)!\n");
            printf("        Face %d: tri = %d, %d/%d  %d/%d  %d/%d\n", iface, itri,
                   btess->tess2d[iface-1].ptype[i1], 
                   btess->tess2d[iface-1].pindex[i1],
                   btess->tess2d[iface-1].ptype[i2], 
                   btess->tess2d[iface-1].pindex[i2],
                   btess->tess2d[iface-1].ptype[i3], 
                   btess->tess2d[iface-1].pindex[i3]);
            return EGADS_TOPOERR;
          }
        }
      if (btess->tess2d[iface-1].pindex[i2] == btess->tess2d[iface-1].pindex[i3])
        if ((btess->tess2d[iface-1].ptype[i2] == 0) &&
            (btess->tess2d[iface-1].ptype[i3] == 0)) {
          if (outLevel > 0) {
            printf(" EGADS Error: Degen EDGE (EG_insertEdgeVerts)!\n");
            printf("        Face %d: tri = %d, %d/%d  %d/%d  %d/%d\n", iface, itri,
                   btess->tess2d[iface-1].ptype[i1], 
                   btess->tess2d[iface-1].pindex[i1],
                   btess->tess2d[iface-1].ptype[i2], 
                   btess->tess2d[iface-1].pindex[i2],
                   btess->tess2d[iface-1].ptype[i3], 
                   btess->tess2d[iface-1].pindex[i3]);
            return EGADS_TOPOERR;
          }
        }
      if (btess->tess2d[iface-1].pindex[i1] == btess->tess2d[iface-1].pindex[i3])
        if ((btess->tess2d[iface-1].ptype[i1] == 0) &&
            (btess->tess2d[iface-1].ptype[i3] == 0)) {
          if (outLevel > 0) {
            printf(" EGADS Error: Degen EDGE (EG_insertEdgeVerts)!\n");
            printf("        Face %d: tri = %d, %d/%d  %d/%d  %d/%d\n", iface, itri,
                   btess->tess2d[iface-1].ptype[i1], 
                   btess->tess2d[iface-1].pindex[i1],
                   btess->tess2d[iface-1].ptype[i2], 
                   btess->tess2d[iface-1].pindex[i2],
                   btess->tess2d[iface-1].ptype[i3], 
                   btess->tess2d[iface-1].pindex[i3]);
            return EGADS_TOPOERR; 
          }
        }
      cnt++;
    }
  }

  stripe = 3 + 2*cnt;
  vals   = (double *) EG_alloc(stripe*npts*sizeof(double));
  if (vals == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on Tmp %d %d (EG_insertEdgeVerts)!\n",
             npts, stripe);
    return EGADS_MALLOC;
  }
  stat = EG_getBodyTopos(btess->src, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) {
    EG_free(vals);
    return stat;
  }
  stat = EG_getBodyTopos(btess->src, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    EG_free(edges);
    EG_free(vals);
    return stat;
  }

  /* get the new data on the Edge and Faces */
  for (i = 0; i < npts; i++) {
    stat = EG_evaluate(edges[eIndex-1], &t[i], result);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      EG_free(edges);
      EG_free(vals);
      return stat;
    }
    vals[stripe*i  ] = result[0];
    vals[stripe*i+1] = result[1];
    vals[stripe*i+2] = result[2];
    for (cnt = m = 0; m < 2; m++) {
      nf = btess->tess1d[eIndex-1].faces[m].nface;
      for (nx = 0; nx < nf; nx++) {
        iface = btess->tess1d[eIndex-1].faces[m].index;
        if (nf > 1) iface = btess->tess1d[eIndex-1].faces[m].faces[nx];
        if (iface == 0) continue;
        sense = faces[iface-1]->mtype;
        if (EG_faceConnIndex(btess->tess1d[eIndex-1].faces[1-m], iface) == 0)
          sense = 0;
        if (m == 0) sense = -sense;
        stat = EG_getEdgeUV(faces[iface-1], edges[eIndex-1], sense, 
                            t[i], &vals[stripe*i+3+2*cnt]);
        if (stat != EGADS_SUCCESS) {
          EG_free(faces);
          EG_free(edges);
          EG_free(vals);
          return stat;
        }
        cnt++;
      }
    }
  }
  EG_free(faces);
  EG_free(edges);
  
  /* get all of the Edge memory we will need */
  xyzs = (double *) EG_alloc(3*(npts+btess->tess1d[eIndex-1].npts)*
                             sizeof(double));
  ts   = (double *) EG_alloc(  (npts+btess->tess1d[eIndex-1].npts)*
                             sizeof(double));
  if ((xyzs == NULL) || (ts == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on Edge %d %d (EG_insertEdgeVerts)!\n",
             npts, btess->tess1d[eIndex-1].npts);
    if (ts   != NULL) EG_free(ts);
    if (xyzs != NULL) EG_free(xyzs);
    EG_free(vals);
    return EGADS_MALLOC;
  }
  etric[0] = etric[1] = NULL;
  nf = btess->tess1d[eIndex-1].faces[0].nface;
  if (nf > 0) {
    etric[0] = (int *) EG_alloc(nf*(npts+btess->tess1d[eIndex-1].npts-1)*
                                sizeof(int));
    if (etric[0] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc on Edge- %d %d (EG_insertEdgeVerts)!\n",
               npts, btess->tess1d[eIndex-1].npts-1);
      EG_free(ts);
      EG_free(xyzs);
      EG_free(vals);
      return EGADS_MALLOC;
    }
  }
  nf = btess->tess1d[eIndex-1].faces[1].nface;
  if (nf > 0) {
    etric[1] = (int *) EG_alloc(nf*(npts+btess->tess1d[eIndex-1].npts-1)*
                                sizeof(int));
    if (etric[1] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc on Edge+ %d %d (EG_insertEdgeVerts)!\n",
               npts, btess->tess1d[eIndex-1].npts-1);
      if (etric[0] != NULL) EG_free(etric[0]);
      EG_free(ts);
      EG_free(xyzs);
      EG_free(vals);
      return EGADS_MALLOC;
    }
  }
  
  /* cleanup mappings */
  EG_cleanupTessMaps(btess);
  
  /* set the new Edge tessellation information */
  for (j = i = 0; i < btess->tess1d[eIndex-1].npts; i++, j++) {
    xyzs[3*j  ] = btess->tess1d[eIndex-1].xyz[3*i  ];
    xyzs[3*j+1] = btess->tess1d[eIndex-1].xyz[3*i+1];
    xyzs[3*j+2] = btess->tess1d[eIndex-1].xyz[3*i+2];
    ts[j]       = btess->tess1d[eIndex-1].t[i];
    if (i != btess->tess1d[eIndex-1].npts-1)
      for (m = 0; m < 2; m++) {
        nf = btess->tess1d[eIndex-1].faces[m].nface;
        for (nx = 0; nx < nf; nx++)
          etric[m][j*nf+nx] = btess->tess1d[eIndex-1].faces[m].tric[i*nf+nx];
      }
    if (i != vIndex-1) continue;
    for (k = 0; k < npts; k++) {
      j++;
      xyzs[3*j  ] = vals[stripe*k  ];
      xyzs[3*j+1] = vals[stripe*k+1];
      xyzs[3*j+2] = vals[stripe*k+2];
      ts[j]       = t[k];
      for (m = 0; m < 2; m++) {
        nf = btess->tess1d[eIndex-1].faces[m].nface;
        for (nx = 0; nx < nf; nx++) etric[m][j*nf+nx] = 0;
      }
    }
  }

  /* do each Face touched by the Edge */
  for (cnt = m = 0; m < 2; m++) {
    nf = btess->tess1d[eIndex-1].faces[m].nface;
    for (nx = 0; nx < nf; nx++) {
      iface = btess->tess1d[eIndex-1].faces[m].index;
      if (nf > 1) iface = btess->tess1d[eIndex-1].faces[m].faces[nx];
      if (iface == 0) continue;
      xyz    = (double *) EG_alloc(3*(npts+btess->tess2d[iface-1].npts)*
                                   sizeof(double));
      uv     = (double *) EG_alloc(2*(npts+btess->tess2d[iface-1].npts)*
                                   sizeof(double));
      ptype  = (int *)    EG_alloc(  (npts+btess->tess2d[iface-1].npts)*
                                   sizeof(int));
      pindex = (int *)    EG_alloc(  (npts+btess->tess2d[iface-1].npts)*
                                   sizeof(int));
      tris   = (int *)    EG_alloc(3*(npts+btess->tess2d[iface-1].ntris)*
                                   sizeof(int));
      tric   = (int *)    EG_alloc(3*(npts+btess->tess2d[iface-1].ntris)*
                                   sizeof(int));
      if ((xyz == NULL)    || (uv == NULL)   || (ptype == NULL) ||
          (pindex == NULL) || (tris == NULL) || (tric == NULL)) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc on Edge %d %d (EG_insertEdgeVerts)!\n",
                 npts, btess->tess1d[eIndex-1].npts);
        if (tric     != NULL) EG_free(tric);
        if (tris     != NULL) EG_free(tris);
        if (pindex   != NULL) EG_free(pindex);
        if (ptype    != NULL) EG_free(ptype);
        if (uv       != NULL) EG_free(uv);
        if (xyz      != NULL) EG_free(xyz);
        if (etric[0] != NULL) EG_free(etric[0]);
        if (etric[1] != NULL) EG_free(etric[1]);
        EG_free(ts);
        EG_free(xyzs);
        EG_free(vals);
        if (cnt != 0) EG_deleteObject(tess);
        return EGADS_MALLOC;
      }
      for (i = 0; i < btess->tess2d[iface-1].npts; i++) {
        xyz[3*i  ] = btess->tess2d[iface-1].xyz[3*i  ];
        xyz[3*i+1] = btess->tess2d[iface-1].xyz[3*i+1];
        xyz[3*i+2] = btess->tess2d[iface-1].xyz[3*i+2];
        uv[2*i  ]  = btess->tess2d[iface-1].uv[2*i  ];
        uv[2*i+1]  = btess->tess2d[iface-1].uv[2*i+1];
        ptype[i]   = btess->tess2d[iface-1].ptype[i];
        pindex[i]  = btess->tess2d[iface-1].pindex[i];
        if (pindex[i] == eIndex) 
          if (ptype[i] > vIndex) ptype[i] += npts;
      }
      j = btess->tess2d[iface-1].npts;
      for (i = 0; i < npts; i++) {
        xyz[3*(j+i)  ] = vals[stripe*i  ];
        xyz[3*(j+i)+1] = vals[stripe*i+1];
        xyz[3*(j+i)+2] = vals[stripe*i+2];
        uv[2*(j+i)  ]  = vals[stripe*i+3+2*cnt  ];
        uv[2*(j+i)+1]  = vals[stripe*i+3+2*cnt+1];
        ptype[j+i]     = vIndex + i+1;
        pindex[j+i]    = eIndex;
      }
      for (i = 0; i < btess->tess2d[iface-1].ntris; i++) {
        tris[3*i  ] = btess->tess2d[iface-1].tris[3*i  ];
        tris[3*i+1] = btess->tess2d[iface-1].tris[3*i+1];
        tris[3*i+2] = btess->tess2d[iface-1].tris[3*i+2];
        tric[3*i  ] = btess->tess2d[iface-1].tric[3*i  ];
        tric[3*i+1] = btess->tess2d[iface-1].tric[3*i+1];
        tric[3*i+2] = btess->tess2d[iface-1].tric[3*i+2];
      }
      
      /* adjust the Face tessellation */
      sense = 1;
      if (etric[m] == NULL) {
        printf(" EGADS Internal: Can not find Triangle %d!\n", m);
        continue;
      }
      itri  = etric[m][(vIndex-1)*nf+nx];
      pt1   =       vIndex;
      pi1   = pi2 = eIndex;
      pt2   = vIndex + 1;
      if (vIndex == 1) {
        pt1 = 0;
        pi1 = btess->tess1d[eIndex-1].nodes[0];
      }
      if (pt2 == btess->tess1d[eIndex-1].npts) {
        pt2 = 0;
        pi2 = btess->tess1d[eIndex-1].nodes[1];
      }
      i1 = tris[3*itri-3]-1;
      i2 = tris[3*itri-2]-1;
      i3 = tris[3*itri-1]-1;        
      if (((btess->tess2d[iface-1].pindex[i2] == pi1) &&
           (btess->tess2d[iface-1].ptype[i2]  == pt1) &&
           (btess->tess2d[iface-1].pindex[i3] == pi2) &&
           (btess->tess2d[iface-1].ptype[i3]  == pt2)) || 
          ((btess->tess2d[iface-1].pindex[i2] == pi2) &&
           (btess->tess2d[iface-1].ptype[i2]  == pt2) &&
           (btess->tess2d[iface-1].pindex[i3] == pi1) &&
           (btess->tess2d[iface-1].ptype[i3]  == pt1))) {
        vert = i1 + 1;
        v0   = i2 + 1;
        v1   = i3 + 1;
        n0   = tric[3*itri-2];
        n1   = tric[3*itri-1];
/*      printf("    0: neighbor = %d\n", tric[m][3*itri-3]);  */
      } else if (((btess->tess2d[iface-1].pindex[i1] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt1) &&
                  (btess->tess2d[iface-1].pindex[i3] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i3]  == pt2)) || 
                 ((btess->tess2d[iface-1].pindex[i1] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt2) &&
                  (btess->tess2d[iface-1].pindex[i3] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i3]  == pt1))) {
        v1   = i1 + 1;
        vert = i2 + 1;
        v0   = i3 + 1;
        n1   = tric[3*itri-3];
        n0   = tric[3*itri-1];
/*      printf("    1: neighbor = %d\n", tric[m][3*itri-2]);  */
      } else if (((btess->tess2d[iface-1].pindex[i1] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt1) &&
                  (btess->tess2d[iface-1].pindex[i2] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i2]  == pt2)) || 
                 ((btess->tess2d[iface-1].pindex[i1] == pi2) &&
                  (btess->tess2d[iface-1].ptype[i1]  == pt2) &&
                  (btess->tess2d[iface-1].pindex[i2] == pi1) &&
                  (btess->tess2d[iface-1].ptype[i2]  == pt1))) {
        v0   = i1 + 1;
        v1   = i2 + 1;
        vert = i3 + 1;
        n0   = tric[3*itri-3];
        n1   = tric[3*itri-2];
/*      printf("    2: neighbor = %d\n", tric[m][3*itri-1]);  */
      } else {
        printf(" EGADS Internal: Can not find segment for %d %d  %d %d - %d!\n",
               pt1, pi1, pt2, pi2, btess->tess1d[eIndex-1].npts);
        /* fill in some values -- we are dropping through */
        v0   = i1 + 1;
        v1   = i2 + 1;
        vert = i3 + 1;
        n0   = tric[3*itri-3];
        n1   = tric[3*itri-2];
      }
      if ((btess->tess2d[iface-1].ptype[v1-1]  == pt1) &&
          (btess->tess2d[iface-1].pindex[v1-1] == pi1)) {
        i     =  v0;
        v0    =  v1;
        v1    =  i;
#ifndef __clang_analyzer__
        i     =  n0;
#endif
        n0    =  n1;
#ifndef __clang_analyzer__
        n1    =  i;
#endif
        sense = -1;
      }
/*    printf("       %d: %d", iface, itri); 
      printf(" vert = %d (%d/%d)  v0 = %d (%d/%d)  v1 = %d (%d/%d), sense = %d\n", 
             vert, btess->tess2d[iface-1].ptype[vert-1], 
                   btess->tess2d[iface-1].pindex[vert-1], 
             v0,   btess->tess2d[iface-1].ptype[v0-1], 
                   btess->tess2d[iface-1].pindex[v0-1], 
             v1,   btess->tess2d[iface-1].ptype[v1-1], 
                   btess->tess2d[iface-1].pindex[v1-1], sense);  */
      for (i = 0; i < 3; i++) {
        if (btess->tess2d[iface-1].tris[3*itri+i-3] == v1) 
          tris[3*itri+i-3] = btess->tess2d[iface-1].npts  + 1;
        if (btess->tess2d[iface-1].tris[3*itri+i-3] == v0) 
          tric[3*itri+i-3] = btess->tess2d[iface-1].ntris + 1;
      } 
      nl = itri;
      for (i = 0; i < npts; i++) {
        j  = btess->tess2d[iface-1].ntris + i;
        v0 = btess->tess2d[iface-1].npts  + i + 1;
        vn = btess->tess2d[iface-1].npts  + i + 2;
        nn = j+2;
        if (i == npts-1) {
          vn = v1;
          nn = n0;
        }
        tris[3*j  ] =  vert;
        tric[3*j  ] = -eIndex;
        if (sense == 1) {
          tris[3*j+1] = v0;
          tris[3*j+2] = vn;
          tric[3*j+1] = nn;
          tric[3*j+2] = nl;
        } else {
          tris[3*j+1] = vn;
          tris[3*j+2] = v0;
          tric[3*j+1] = nl;
          tric[3*j+2] = nn;
        }
        etric[m][nf*(vIndex+i)+nx] = j + 1;
        nl = j+1;
      }
      if (n0 > 0) {
        for (i = 0; i < 3; i++)
          if (btess->tess2d[iface-1].tric[3*n0+i-3] == itri) 
            tric[3*n0+i-3] = btess->tess2d[iface-1].ntris+npts;
      } else if (n0 < 0) {
        j = -n0 - 1;
        if (btess->tess1d[j].faces[0].index == iface) {
          for (k = i = 0; i < btess->tess1d[j].npts-1; i++)
            if (etric[m][nf*i+nx] == itri) k++;
          for (i = 0; i < btess->tess1d[j].npts-1; i++) {
            if ((k > 1) && (i >= vIndex-1) && (i < vIndex+npts-1)) continue;
            if (etric[m][nf*i+nx] == itri)
              etric[m][nf*i+nx] = btess->tess2d[iface-1].ntris+npts;
          }
        }
      }

      /* update the Face pointers */
      if (btess->tess2d[iface-1].xyz    != NULL) 
        EG_free(btess->tess2d[iface-1].xyz);
      if (btess->tess2d[iface-1].uv     != NULL) 
        EG_free(btess->tess2d[iface-1].uv);
      if (btess->tess2d[iface-1].ptype  != NULL) 
        EG_free(btess->tess2d[iface-1].ptype);
      if (btess->tess2d[iface-1].pindex != NULL) 
        EG_free(btess->tess2d[iface-1].pindex);
      if (btess->tess2d[iface-1].bary   != NULL)
        EG_free(btess->tess2d[iface-1].bary);
      if (btess->tess2d[iface-1].frame  != NULL)
        EG_free(btess->tess2d[iface-1].frame);
      if (btess->tess2d[iface-1].frlps != NULL)
        EG_free(btess->tess2d[iface-1].frlps);
      if (btess->tess2d[iface-1].tris   != NULL)
        EG_free(btess->tess2d[iface-1].tris);
      if (btess->tess2d[iface-1].tric   != NULL) 
        EG_free(btess->tess2d[iface-1].tric);
      btess->tess2d[iface-1].xyz    = xyz;
      btess->tess2d[iface-1].uv     = uv;
      btess->tess2d[iface-1].ptype  = ptype;
      btess->tess2d[iface-1].pindex = pindex;
      btess->tess2d[iface-1].bary   = NULL;
      btess->tess2d[iface-1].frame  = NULL;
      btess->tess2d[iface-1].frlps  = NULL;
      btess->tess2d[iface-1].tris   = tris;
      btess->tess2d[iface-1].tric   = tric;
      btess->tess2d[iface-1].ntris += npts;
      btess->tess2d[iface-1].npts  += npts;

      /* delete any quads and mark the frame as invalid */
      EG_deleteQuads(btess, iface);

      cnt++;
    }
  }
  EG_free(vals);
  
  /* set the updated Edge tessellation */
  if (btess->tess1d[eIndex-1].faces[0].tric != NULL) 
    EG_free(btess->tess1d[eIndex-1].faces[0].tric);
  if (btess->tess1d[eIndex-1].faces[1].tric != NULL) 
    EG_free(btess->tess1d[eIndex-1].faces[1].tric);
  btess->tess1d[eIndex-1].faces[0].tric = etric[0];
  btess->tess1d[eIndex-1].faces[1].tric = etric[1];
  if (btess->tess1d[eIndex-1].xyz  != NULL) 
    EG_free(btess->tess1d[eIndex-1].xyz);
  if (btess->tess1d[eIndex-1].t    != NULL) 
    EG_free(btess->tess1d[eIndex-1].t);
  btess->tess1d[eIndex-1].xyz   = xyzs;
  btess->tess1d[eIndex-1].t     = ts;
  btess->tess1d[eIndex-1].npts += npts;
  
#ifdef CHECK
  EG_checkTriangulation(btess);
#endif

  return EGADS_SUCCESS;
}


static void
EG_tessThread(void *struc)
{
  int          i, index, stat, aStat, aType, aLen;
#ifdef PROGRESS
  int          outLevel;
#endif
  long         ID;
  double       dist, params[3];
  const int    *aInts;
  const double *aReals;
  const char   *aStr;
  triStruct    tst;
  fillArea     fast;
  EMPtess      *tthread;
  
  tthread = (EMPtess *) struc;
#ifdef PROGRESS
  outLevel = EG_outLevel(tthread->body);
#endif

  /* get our identifier */
  ID = EMP_ThreadID();
  
  dist = fabs(tthread->params[2]);
  if (dist > 30.0) dist = 30.0;
  if (dist <  0.5) dist =  0.5;
  tst.maxlen   = tthread->params[0];
  tst.chord    = tthread->params[1];
  tst.dotnrm   = cos(PI*dist/180.0);
  tst.minlen   = tthread->tparam[0];
  tst.maxPts   = tthread->tparam[1];
  tst.qparm[0] = tthread->qparam[0];
  tst.qparm[1] = tthread->qparam[1];
  tst.qparm[2] = tthread->qparam[2];
  tst.mverts   = tst.nverts = 0;
  tst.verts    = NULL;
  tst.mtris    = tst.ntris  = 0;
  tst.tris     = NULL;
  tst.msegs    = tst.nsegs  = 0;
  tst.segs     = NULL;
  tst.mframe   = tst.nframe = 0;
  tst.frame    = NULL;
  tst.mloop    = tst.nloop  = 0;
  tst.loop     = NULL;
  tst.numElem  = -1;
  tst.hashTab  = NULL;
  
  fast.pts     = NULL;
  fast.segs    = NULL;
  fast.front   = NULL;
 
  /* look for work */
  for (;;) {
    
    /* only one thread at a time here -- controlled by a mutex! */
    if (tthread->mutex != NULL) EMP_LockSet(tthread->mutex);
    if (tthread->mark == NULL) {
      index = tthread->index;
    } else {
      for (index = tthread->index; index < tthread->end; index++) {
        if (tthread->mark[index] == 0) continue;
        break;
      }
    }
    tthread->index = index+1;
    if (tthread->mutex != NULL) EMP_LockRelease(tthread->mutex);
    if (index >= tthread->end) break;
#ifdef PROGRESS
    if (outLevel > 0) {
      printf("    tessellating Face %3d of %3d\r", index+1, tthread->end);
      fflush(stdout);
    }
#endif
    
    /* adjust the parameters? */
    aStat = EG_attributeRet(tthread->faces[index], ".tParams", &aType, &aLen,
                            &aInts, &aReals, &aStr);
    if (aStat == EGADS_SUCCESS)
      if (aType != ATTRREAL) {
        printf(" EGADS Warning: tParams NonReal Attribute (EG_tessThread)!\n");
      } else {
        params[0] = tthread->params[0];
        params[1] = tthread->params[1];
        params[2] = tthread->params[2];
        for (i = 0; i < aLen; i++) {
          if (i == 3) break;
          if ((aReals[i] < params[i]) && (aReals[i] > 0.0))
            params[i] = aReals[i];
        }
        dist = fabs(params[2]);
        if (dist > 30.0) dist = 30.0;
        if (dist <  0.5) dist =  0.5;
        tst.maxlen = params[0];
        tst.chord  = params[1];
        tst.dotnrm = cos(PI*dist/180.0);
      }
    aStat = EG_attributeRet(tthread->faces[index], ".tParam", &aType, &aLen,
                            &aInts, &aReals, &aStr);
    if (aStat == EGADS_SUCCESS)
      if (aType != ATTRREAL) {
        printf(" EGADS Warning: tParam NonReal Attribute (EG_tessThread)!\n");
      } else {
        params[0] = tthread->params[0];
        params[1] = tthread->params[1];
        params[2] = tthread->params[2];
        for (i = 0; i < aLen; i++) {
          if (i == 3) break;
          if (aReals[i] > 0.0) params[i] = aReals[i];
        }
        dist = fabs(params[2]);
        if (dist > 30.0) dist = 30.0;
        if (dist <  0.5) dist =  0.5;
        tst.maxlen = params[0];
        tst.chord  = params[1];
        tst.dotnrm = cos(PI*dist/180.0);
      }
    tst.qparm[0] = tthread->qparam[0];
    tst.qparm[1] = tthread->qparam[1];
    tst.qparm[2] = tthread->qparam[2];
    aStat = EG_attributeRet(tthread->faces[index], ".qParams", &aType, &aLen,
                            &aInts, &aReals, &aStr);
    if (aStat == EGADS_SUCCESS)
      if ((aType == ATTRINT) || (aType == ATTRCSYS)) {
        printf(" EGADS Warning: qParam Integer/CSys Attribute (EG_tessThread)!\n");
      } else if (aType == ATTRREAL) {
        for (i = 0; i < aLen; i++) {
          if (i == 3) break;
          tst.qparm[i] = aReals[i];
        }
      } else {
        tst.qparm[0] = -1.0;
      }

    /* do the work */
    stat = EG_fillTris(tthread->body, index+1, tthread->faces[index],
                       tthread->tess, &tst, &fast, ID);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: Face %d -> EG_fillTris = %d (EG_tessThread)!\n",
             index+1, stat);
    
    /* restore the tParams parameters */
    dist = fabs(tthread->params[2]);
    if (dist > 30.0) dist = 30.0;
    if (dist <  0.5) dist =  0.5;
    tst.maxlen = tthread->params[0];
    tst.chord  = tthread->params[1];
    tst.dotnrm = cos(PI*dist/180.0);
  }
  
  /* exhausted all work -- cleanup & exit */
  if (tst.verts  != NULL) EG_free(tst.verts);
  if (tst.tris   != NULL) EG_free(tst.tris);
  if (tst.segs   != NULL) EG_free(tst.segs);
  if (tst.frame  != NULL) EG_free(tst.frame);
  if (tst.loop   != NULL) EG_free(tst.loop);
  
  if (fast.segs  != NULL) EG_free(fast.segs);
  if (fast.pts   != NULL) EG_free(fast.pts);
  if (fast.front != NULL) EG_free(fast.front);
  
  if (ID != tthread->master) EMP_ThreadExit();
}


int
EG_makeTessBody(egObject *object, double *paramx, egObject **tess)
{
  int      i, j, stat, outLevel, nface, np, aStat, aType, aLen;
  double   params[3];
  void     **threads = NULL;
  long     start;
  egTessel *btess;
  egObject *ttess, *context, **faces;
  egCntxt  *cntx;
  EMPtess  tthread;
  const int    *aInts;
  const double *aReals;
  const char   *aStr;

  *tess = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != BODY)       return EGADS_NOTBODY;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  if (context == NULL)              return EGADS_NULLOBJ;
  cntx     = (egCntxt *) context->blind;
  
  /* get global settings */
  params[0] = fabs(paramx[0]);
  params[1] =      paramx[1];
  params[2] =      paramx[2];
  aStat     = EG_attributeRet(object, ".tParams", &aType, &aLen, &aInts,
                              &aReals, &aStr);
  if (aStat == EGADS_SUCCESS)
    if (aType != ATTRREAL) {
      printf(" EGADS Warning: tParams NonReal Attribute (EG_makeTessBody)!\n");
    } else {
      for (i = 0; i < aLen; i++) {
        if (i == 3) break;
        if ((aReals[i] < params[i]) && (aReals[i] > 0.0))
          params[i] = aReals[i];
      }
    }
  
  /* get quadding parameters, if any */
  aStat = EG_attributeRet(object, ".qParams", &aType, &aLen, &aInts,
                          &aReals, &aStr);
  
  btess = (egTessel *) EG_alloc(sizeof(egTessel));
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Blind Malloc (EG_makeTessBody)!\n");
    return EGADS_MALLOC;
  }
  btess->src       = object;
  btess->xyzs      = NULL;
  btess->tess1d    = NULL;
  btess->tess2d    = NULL;
  btess->globals   = NULL;
  btess->nGlobal   = 0;
  btess->nEdge     = 0;
  btess->nFace     = 0;
  btess->nu        = 0;
  btess->nv        = 0;
  btess->done      = 1;
  btess->params[0] = params[0];
  btess->params[1] = params[1];
  btess->params[2] = params[2];
  for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = 0.0;
  if (cntx != NULL)
    for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = cntx->tess[i];
  
  /* do the Edges & make the Tessellation Object */
  
  stat = EG_tessEdges(btess, NULL);
  if (stat != EGADS_SUCCESS) {
    EG_cleanupTess(btess);
    EG_free(btess);
    return stat;  
  }
  stat = EG_makeObject(context, &ttess);
  if (stat != EGADS_SUCCESS) {
    EG_cleanupTess(btess);
    EG_free(btess);
    return stat;
  }
  ttess->oclass = TESSELLATION;
  ttess->blind  = btess;
  EG_referenceObject(ttess,  context);
  EG_referenceTopObj(object, ttess);
  *tess = ttess;
  
  /* Wire Body or Edges Only */
  if ((object->mtype == WIREBODY) || (paramx[0] < 0.0)) return EGADS_SUCCESS;

  /* Need Face triangulations */
  stat = EG_getBodyTopos(object, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos = %d (EG_makeTessBody)!\n",
           stat);
    EG_deleteObject(ttess);
    *tess = NULL;
    return stat;
  }
  btess->tess2d = (egTess2D *) EG_alloc(2*nface*sizeof(egTess2D));
  if (btess->tess2d == NULL) {
    printf(" EGADS Error: Alloc %d Faces (EG_makeTessBody)!\n", nface);  
    EG_deleteObject(ttess);
    *tess = NULL;
    return EGADS_MALLOC;
  }
  for (j = 0; j < 2*nface; j++) {
    btess->tess2d[j].mKnots = NULL;
    btess->tess2d[j].xyz    = NULL;
    btess->tess2d[j].uv     = NULL;
    btess->tess2d[j].global = NULL;
    btess->tess2d[j].ptype  = NULL;
    btess->tess2d[j].pindex = NULL;
    btess->tess2d[j].bary   = NULL;
    btess->tess2d[j].frame  = NULL;
    btess->tess2d[j].frlps  = NULL;
    btess->tess2d[j].tris   = NULL;
    btess->tess2d[j].tric   = NULL;
    btess->tess2d[j].patch  = NULL;
    btess->tess2d[j].npts   = 0;
    btess->tess2d[j].nframe = 0;
    btess->tess2d[j].nfrlps = 0;
    btess->tess2d[j].ntris  = 0;
    btess->tess2d[j].npatch = 0;
  }
  btess->nFace = nface;
  
  /* set up for explicit multithreading */
  tthread.mutex     = NULL;
  tthread.master    = EMP_ThreadID();
  tthread.index     = 0;
  tthread.end       = nface;
  tthread.mark      = NULL;
  tthread.tess      = ttess;
  tthread.btess     = NULL;
  tthread.body      = object;
  tthread.faces     = faces;
  tthread.edges     = NULL;
  tthread.params    = params;
  tthread.tparam    = btess->tparam;
  tthread.qparam[0] = tthread.qparam[1] = tthread.qparam[2] = 0.0;
  if (aStat == EGADS_SUCCESS) {
    if ((aType == ATTRINT) || (aType == ATTRCSYS)) {
      printf(" EGADS Warning: qParam Integer/CSys Attribute (EG_makeTessBody)!\n");
    } else if (aType == ATTRREAL) {
      for (i = 0; i < aLen; i++) {
        if (i == 3) break;
        tthread.qparam[i] = aReals[i];
      }
    } else {
      tthread.qparam[0] = -1.0;
    }
  }
  
  np = EMP_Init(&start);
  if (outLevel > 1) printf(" EMP NumProcs = %d!\n", np);
  
  if (np > 1) {
    /* create the mutex to handle list synchronization */
    tthread.mutex = EMP_LockCreate();
    if (tthread.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(tthread.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EG_tessThread, &tthread);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  EG_tessThread(&tthread);
  
  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);
#ifdef PROGRESS
  if (outLevel > 0) printf("\n");
#endif
#ifdef CHECK
  EG_checkTriangulation(btess);
#endif

  /* cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  if (tthread.mutex != NULL) EMP_LockDestroy(tthread.mutex);
  if (threads != NULL) free(threads);
  EG_free(faces);
  if (outLevel > 1)
    printf(" EMP Number of Seconds on Face Thread Block = %ld\n",
           EMP_Done(&start));

  return EGADS_SUCCESS;
}


int
EG_remakeTess(egObject *tess, int nobj, egObject **objs, double *paramx)
{
  int      i, j, mx, stat, outLevel, iface, nface, hit, np, aStat, aType, aLen;
  int      *ed, *marker = NULL;
  double   params[3];
  void     **threads = NULL;
  long     start;
  double   save[3];
  egObject *context, *object, **faces;
  egTessel *btess;
  egCntxt  *cntx;
  EMPtess  tthread;
  const int    *aInts;
  const double *aReals;
  const char   *aStr;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;
  btess  = tess->blind;
  if (btess->done != 1)             return EGADS_TESSTATE;
  object = btess->src;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != BODY)       return EGADS_NOTBODY;
  if (nobj <= 0)                    return EGADS_NODATA;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  if (context == NULL)              return EGADS_NULLOBJ;
  cntx     = (egCntxt *) context->blind;

  /* get global settings */
  params[0] = paramx[0];
  params[1] = paramx[1];
  params[2] = paramx[2];
  aStat     = EG_attributeRet(object, ".tParams", &aType, &aLen, &aInts,
                              &aReals, &aStr);
  if (aStat == EGADS_SUCCESS)
    if (aType != ATTRREAL) {
      printf(" EGADS Warning: tParams NonReal Attribute (EG_remakeTess)!\n");
    } else {
      for (i = 0; i < aLen; i++) {
        if (i == 3) break;
        if ((aReals[i] < params[i]) && (aReals[i] > 0.0))
          params[i] = aReals[i];
      }
    }

  /* get quadding parameters, if any */
  aStat = EG_attributeRet(object, ".qParams", &aType, &aLen, &aInts,
                          &aReals, &aStr);

  for (hit = j = 0; j < nobj; j++) {
    if (objs[j] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Object[%d] (EG_remakeTess)!\n",
               j+1);
      return EGADS_NULLOBJ;
    }
    if (objs[j]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Not an Object[%d] (EG_remakeTess)!\n",
               j+1);
      return EGADS_NOTOBJ;
    }
    if ((objs[j]->oclass != EDGE) && (objs[j]->oclass != FACE)) {
      if (outLevel > 0)
        printf(" EGADS Error: Not Edge/Face[%d] (EG_remakeTess)!\n",
               j+1);
      return EGADS_NOTOBJ;
    }
    stat = EG_indexBodyTopo(object, objs[j]);
    if (stat == EGADS_NOTFOUND) {
      if (outLevel > 0)
        printf(" EGADS Error: Object[%d] Not in Body (EG_remakeTess)!\n",
               j+1);
      return stat;
    }
    if (objs[j]->oclass == FACE) continue;
    if (objs[j]->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge[%d] is DEGENERATE (EG_remakeTess)!\n",
               j+1);
      return EGADS_DEGEN;
    }
    hit++;
  }
  
  /* cleanup local/global mappings */
  EG_cleanupTessMaps(btess);

  /* mark faces */
  
  if (btess->nFace != 0) {
    marker = (int *) EG_alloc(btess->nFace*sizeof(int));
    if (marker == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: MALLOC on %d Faces (EG_remakeTess)!\n",
               btess->nFace);
      return EGADS_MALLOC;
    }
    for (j = 0; j < btess->nFace; j++) marker[j] = 0;
    for (j = 0; j < nobj; j++) {
      i = EG_indexBodyTopo(object, objs[j]);
      if (objs[j]->oclass == EDGE) {
        for (mx = 0; mx < 2; mx++) {
          iface = btess->tess1d[i-1].faces[mx].index;
          if (iface == 0) continue;
          marker[iface-1] = 1;
          EG_deleteQuads(btess, iface);
        }
      } else {
        marker[i-1] = 1;
      }
    }
  }
  
  /* do egdes */
  
  if (hit != 0) {
    ed = (int *) EG_alloc(btess->nEdge*sizeof(int));
    if (ed == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: MALLOC on %d Faces (EG_remakeTess)!\n",
               btess->nFace);
      EG_free(marker);
      return EGADS_MALLOC;
    }
    for (j = 0; j < btess->nEdge; j++) ed[j] = 0;
    for (j = 0; j < nobj; j++) {
      if (objs[j]->oclass != EDGE) continue;
      i = EG_indexBodyTopo(object, objs[j]);
      if (btess->tess1d[i-1].xyz != NULL) EG_free(btess->tess1d[i-1].xyz);
      if (btess->tess1d[i-1].t   != NULL) EG_free(btess->tess1d[i-1].t);
      if (btess->tess1d[i-1].faces[0].tric  != NULL)
        EG_free(btess->tess1d[i-1].faces[0].tric);
      if (btess->tess1d[i-1].faces[1].tric  != NULL)
        EG_free(btess->tess1d[i-1].faces[1].tric);
      btess->tess1d[i-1].faces[0].tric = NULL;
      btess->tess1d[i-1].faces[1].tric = NULL;
      btess->tess1d[i-1].xyz           = NULL;
      btess->tess1d[i-1].t             = NULL;
      btess->tess1d[i-1].npts          = 0;
      ed[i-1] = 1;
    }
    save[0] = btess->params[0];
    save[1] = btess->params[1];
    save[2] = btess->params[2];
    btess->params[0] = params[0];
    btess->params[1] = params[1];
    btess->params[2] = params[2];
    stat = EG_tessEdges(btess, ed);
    btess->params[0] = save[0];
    btess->params[1] = save[1];
    btess->params[2] = save[2];
    EG_free(ed);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_tessEdges =  %d (EG_remakeTess)!\n",
               stat);
      EG_free(marker);
      return stat;
    }
  }
  if (marker == NULL) return EGADS_SUCCESS;
 
  stat = EG_getBodyTopos(object, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_getBodyTopos = %d (EG_remakeTess)!\n",
           stat);
    EG_free(marker);
    return stat;
  }
  /* cleanup old Face tessellations */
  for (j = 0; j < btess->nFace; j++) {
    if (marker[j] == 0) continue;
    
    if (btess->tess2d[j].xyz    != NULL) EG_free(btess->tess2d[j].xyz);
    if (btess->tess2d[j].uv     != NULL) EG_free(btess->tess2d[j].uv);
    if (btess->tess2d[j].ptype  != NULL) EG_free(btess->tess2d[j].ptype);
    if (btess->tess2d[j].pindex != NULL) EG_free(btess->tess2d[j].pindex);
    if (btess->tess2d[j].bary   != NULL) EG_free(btess->tess2d[j].bary);
    if (btess->tess2d[j].frlps  != NULL) EG_free(btess->tess2d[j].frlps);
    if (btess->tess2d[j].frame  != NULL) EG_free(btess->tess2d[j].frame);
    if (btess->tess2d[j].tris   != NULL) EG_free(btess->tess2d[j].tris);
    if (btess->tess2d[j].tric   != NULL) EG_free(btess->tess2d[j].tric);
    btess->tess2d[j].xyz    = NULL;
    btess->tess2d[j].uv     = NULL;
    btess->tess2d[j].ptype  = NULL;
    btess->tess2d[j].pindex = NULL;
    btess->tess2d[j].bary   = NULL;
    btess->tess2d[j].frlps  = NULL;
    btess->tess2d[j].frame  = NULL;
    btess->tess2d[j].tris   = NULL;
    btess->tess2d[j].tric   = NULL;
    btess->tess2d[j].npts   = 0;
    btess->tess2d[j].nfrlps = 0;
    btess->tess2d[j].nframe = 0;
    btess->tess2d[j].ntris  = 0;
  }
  for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = 0.0;
  if (cntx != NULL)
    for (i = 0; i < MTESSPARAM; i++) btess->tparam[i] = cntx->tess[i];

  /* set up for explicit multithreading */
  tthread.mutex     = NULL;
  tthread.master    = EMP_ThreadID();
  tthread.index     = 0;
  tthread.end       = btess->nFace;
  tthread.mark      = marker;
  tthread.tess      = tess;
  tthread.btess     = NULL;
  tthread.body      = object;
  tthread.faces     = faces;
  tthread.edges     = NULL;
  tthread.params    = params;
  tthread.tparam    = btess->tparam;
  tthread.qparam[0] = tthread.qparam[1] = tthread.qparam[2] = 0.0;
  if (aStat == EGADS_SUCCESS) {
    if ((aType == ATTRINT) || (aType == ATTRCSYS)) {
      printf(" EGADS Warning: qParam Integer/Csys Attribute (EG_remakeTess)!\n");
    } else if (aType == ATTRREAL) {
      for (i = 0; i < aLen; i++) {
        if (i == 3) break;
        tthread.qparam[i] = aReals[i];
      }
    } else {
      tthread.qparam[0] = -1.0;
    }
  }
  
  np = EMP_Init(&start);
  if (outLevel > 1) printf(" EMP NumProcs = %d!\n", np);
  
  if (np > 1) {
    /* create the mutex to handle list synchronization */
    tthread.mutex = EMP_LockCreate();
    if (tthread.mutex == NULL) {
      printf(" EMP Error: mutex creation = NULL!\n");
      np = 1;
    } else {
      /* get storage for our extra threads */
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(tthread.mutex);
        np = 1;
      }
    }
  }
  
  /* create the threads and get going! */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EG_tessThread, &tthread);
      if (threads[i] == NULL)
        printf(" EMP Error Creating Thread #%d!\n", i+1);
    }
  /* now run the thread block from the original thread */
  EG_tessThread(&tthread);
  
  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);
#ifdef PROGRESS
  if (outLevel > 0) printf("\n");
#endif
#ifdef CHECK
  EG_checkTriangulation(btess);
#endif
  
  /* cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);
  if (tthread.mutex != NULL) EMP_LockDestroy(tthread.mutex);
  if (threads != NULL) free(threads);
  EG_free(faces);
  EG_free(marker);
  if (outLevel > 1)
    printf(" EMP Number of Seconds on Face Thread Block = %ld\n",
           EMP_Done(&start));
  
  return EGADS_SUCCESS;
}


int
EG_mapTessBody(egObject *tess, egObject *body, egObject **mapTess)
{
  int       i, j, k, n, i0, i1, i2, ifrm, stat, outLevel, nnode, nedge, nface;
  int       oclass, mtype, alen, aType, nc, sen, iloop, last, nxt, lst, OK;
  int       *senses, *ptype, *pindex, *frame, *frlps, *tris, *tric;
  const int *nMap, *eMap, *fMap;
  double    param, tr, range[4], result[18], *xyz, *uv, *s;
  egObject  *context, *tessb, *mapObj, *geom, *obj2D;
  egObject  **nodes, **edges, **faces, **child;
  egTessel  *btess, *mtess;
  egBary    *bary;
  
  *mapTess = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;
  btess = tess->blind;
  tessb = btess->src;
  if (tessb == NULL)                return EGADS_NULLOBJ;
  if (tessb->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if (tessb->oclass != BODY)        return EGADS_NOTBODY;
  if (body == NULL)                 return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (body->oclass != BODY)         return EGADS_NOTBODY;
  context  = EG_context(body);
  outLevel = EG_outLevel(body);

  /* get mappings if the attributes exist */
  stat = EG_attributeRet(body, ".fMap", &aType, &alen, &fMap, NULL, NULL);
  if (stat == EGADS_SUCCESS) {
    if (aType != ATTRINT) {
      printf(" EGADS Error: fMap Attribute Not an Int (EG_mapTessBody)!\n");
      return EGADS_ATTRERR;
    }
    if (alen != btess->nFace) {
      printf(" EGADS Error: fLen Mismatch %d %d (EG_mapTessBody)!\n",
             alen, btess->nFace);
      return EGADS_TOPOERR;
    }
    for (i = 0; i < btess->nFace; i++)
      if ((fMap[i] < 1) || (fMap[i] > btess->nFace)) {
        printf(" EGADS Error: Face Map %d = %d [1-%d] (EG_mapTessBody)!\n",
               i+1, fMap[i], btess->nFace);
        return EGADS_TOPOERR;
      }
  }
  stat = EG_attributeRet(body, ".eMap", &aType, &alen, &eMap, NULL, NULL);
  if (stat == EGADS_SUCCESS) {
    if (aType != ATTRINT) {
      printf(" EGADS Error: eMap Attribute Not an Int (EG_mapTessBody)!\n");
      return EGADS_ATTRERR;
    }
    if (alen != btess->nEdge) {
      printf(" EGADS Error: eLen Mismatch %d %d (EG_mapTessBody)!\n",
             alen, btess->nEdge);
      return EGADS_TOPOERR;
    }
    for (i = 0; i < btess->nEdge; i++)
      if ((eMap[i] < 1) || (eMap[i] > btess->nEdge)) {
        printf(" EGADS Error: Edge Map %d = %d [1-%d] (EG_mapTessBody)!\n",
               i+1, eMap[i], btess->nEdge);
        return EGADS_TOPOERR;
      }
  }
  stat = EG_attributeRet(body, ".nMap", &aType, &alen, &nMap, NULL, NULL);
  if (stat == EGADS_SUCCESS) {
    if (aType != ATTRINT) {
      printf(" EGADS Error: nMap Attribute Not an Int (EG_mapTessBody)!\n");
      return EGADS_ATTRERR;
    }
    stat = EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: getBodyTopo NODES = %d (EG_mapTessBody)!\n",
               stat);
      return stat;
    }
    if (alen != nnode) {
      printf(" EGADS Error: nLen Mismatch %d %d (EG_mapTessBody)!\n",
             alen, nnode);
      return EGADS_TOPOERR;
    }
    for (i = 0; i < nnode; i++)
      if ((nMap[i] < 1) || (nMap[i] > nnode)) {
        printf(" EGADS Error: Node Map %d = %d [1-%d] (EG_mapTessBody)!\n",
               i+1, nMap[i], nnode);
        return EGADS_TOPOERR;
      }
    EG_free(nodes);
  }

  /* if we have one we should have them all */
  if (fMap == NULL) {
    if ((nMap != NULL) || (eMap != NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Node and/or Egde Map exists (EG_mapTessBody)!\n");
      return EGADS_NODATA;
    }
    stat = EG_sameBodyTopo(tessb, body);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: sameBodyTopo = %d (EG_mapTessBody)!\n", stat);
      return stat;
    }
  } else {
    if ((nMap == NULL) || (eMap == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Node and/or Egde Map missing (EG_mapTessBody)!\n");
      return EGADS_NODATA;
    }
  }
  
  /* make sure we have all of the Face frames */
  for (j = 0; j < btess->nFace; j++)
    if (btess->tess2d[j].frame == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Frame %d is NULL (EG_mapTessBody)!\n", j+1);
      return EGADS_NOTFOUND;
    }
  
  /* make the object data */
  mtess = (egTessel *) EG_alloc(sizeof(egTessel));
  if (mtess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Blind Malloc (EG_mapTessBody)!\n");
    return EGADS_MALLOC;
  }
  mtess->src       = body;
  mtess->xyzs      = NULL;
  mtess->tess1d    = NULL;
  mtess->tess2d    = NULL;
  mtess->globals   = NULL;
  mtess->nGlobal   = 0;
  mtess->nEdge     = btess->nEdge;
  mtess->nFace     = btess->nFace;
  mtess->nu        = 0;
  mtess->nv        = 0;
  mtess->done      = 1;
  mtess->params[0] = btess->params[0];
  mtess->params[1] = btess->params[1];
  mtess->params[2] = btess->params[2];

  /* set up for the Edges */
  stat = EG_getBodyTopos(body, NULL, EDGE, &nedge, &edges);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getBodyTopo Edges = %d (EG_mapTessBody)!\n",
             stat);
    EG_cleanupTess(mtess);
    EG_free(mtess);
    return stat;
  }

  mtess->tess1d = (egTess1D *) EG_alloc(mtess->nEdge*sizeof(egTess1D));
  if (mtess->tess1d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Alloc %d Edges (EG_mapTessBody)!\n", mtess->nEdge);
    EG_cleanupTess(mtess);
    EG_free(mtess);
    return EGADS_MALLOC;
  }
  for (j = 0; j < mtess->nEdge; j++) {
    k = j;
    if (eMap != NULL) k = eMap[j] - 1;
    mtess->tess1d[k].obj            = edges[k];
    mtess->tess1d[k].faces[0].index = btess->tess1d[j].faces[0].index;
    mtess->tess1d[k].faces[0].nface = btess->tess1d[j].faces[0].nface;
    mtess->tess1d[k].faces[0].faces = NULL;
    mtess->tess1d[k].faces[0].tric  = NULL;
    mtess->tess1d[k].faces[1].index = btess->tess1d[j].faces[1].index;
    mtess->tess1d[k].faces[1].nface = btess->tess1d[j].faces[1].nface;
    mtess->tess1d[k].faces[1].faces = NULL;
    mtess->tess1d[k].faces[1].tric  = NULL;
    mtess->tess1d[k].nodes[0]       = btess->tess1d[j].nodes[0];
    mtess->tess1d[k].nodes[1]       = btess->tess1d[j].nodes[1];
    mtess->tess1d[k].xyz            = NULL;
    mtess->tess1d[k].t              = NULL;
    mtess->tess1d[k].global         = NULL;
    mtess->tess1d[k].npts           = btess->tess1d[j].npts;
    if (fMap != NULL) {
      i0 = mtess->tess1d[k].faces[0].index;
      if (i0 != 0) mtess->tess1d[k].faces[0].index = fMap[i0-1];
      i0 = mtess->tess1d[k].faces[1].index;
      if (i0 != 0) mtess->tess1d[k].faces[1].index = fMap[i0-1];
    }
    if (nMap != NULL) {
      i0 = mtess->tess1d[k].nodes[0];
      if (i0 != 0) mtess->tess1d[k].nodes[0] = nMap[i0-1];
      i0 = mtess->tess1d[k].nodes[1];
      if (i0 != 0) mtess->tess1d[k].nodes[1] = nMap[i0-1];
    }
  }
  EG_free(edges);

  /* set up for the Faces */
  if (mtess->nFace > 0) {
    mtess->tess2d = (egTess2D *) EG_alloc(2*mtess->nFace*sizeof(egTess2D));
    if (mtess->tess2d == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Alloc %d Faces (EG_mapTessBody)!\n", mtess->nFace);
      EG_cleanupTess(mtess);
      EG_free(mtess);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 2*mtess->nFace; j++) {
      k = j;
      if (fMap != NULL)
        if (j < mtess->nFace) {
          k = fMap[j] - 1;
        } else {
          k = fMap[j-mtess->nFace] + mtess->nFace - 1;
        }
      mtess->tess2d[k].mKnots = NULL;
      mtess->tess2d[k].xyz    = NULL;
      mtess->tess2d[k].uv     = NULL;
      mtess->tess2d[k].global = NULL;
      mtess->tess2d[k].ptype  = NULL;
      mtess->tess2d[k].pindex = NULL;
      mtess->tess2d[k].bary   = NULL;
      mtess->tess2d[k].frame  = NULL;
      mtess->tess2d[k].frlps  = NULL;
      mtess->tess2d[k].tris   = NULL;
      mtess->tess2d[k].tric   = NULL;
      mtess->tess2d[k].patch  = NULL;
      mtess->tess2d[k].npts   = btess->tess2d[j].npts;
      mtess->tess2d[k].nframe = btess->tess2d[j].nframe;
      mtess->tess2d[k].nfrlps = btess->tess2d[j].nfrlps;
      mtess->tess2d[k].ntris  = btess->tess2d[j].ntris;
      mtess->tess2d[k].npatch = 0;
    }
  }
  
  stat = EG_getBodyTopos(body, NULL, NODE, &nnode, &nodes);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getBodyTopo Nodes = %d (EG_mapTessBody)!\n",
             stat);
    EG_cleanupTess(mtess);
    EG_free(mtess);
    return stat;
  }
  stat = EG_getBodyTopos(body, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getBodyTopo Faces = %d (EG_mapTessBody)!\n",
             stat);
    EG_free(nodes);
    EG_cleanupTess(mtess);
    EG_free(mtess);
    return stat;
  }

  /* fill the Edges */

  for (j = 0; j < mtess->nEdge; j++) {
    k = j;
    if (eMap != NULL) k = eMap[j] - 1;
    n = mtess->tess1d[k].faces[0].nface;
    if (n > 1) {
      mtess->tess1d[k].faces[0].faces = (int *) EG_alloc(n*sizeof(int));
      if (mtess->tess1d[k].faces[0].faces == NULL) {
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return EGADS_MALLOC;
      }
      for (i = 0; i < n; i++) {
        i0 = btess->tess1d[j].faces[0].faces[i];
        if ((i0 != 0) && (fMap != NULL)) i0 = fMap[i0-1];
        mtess->tess1d[k].faces[0].faces[i] = i0;
      }
    }
    n *= mtess->tess1d[k].npts-1;
    if (n > 1) {
      mtess->tess1d[k].faces[0].tric = (int *) EG_alloc(n*sizeof(int));
      if (mtess->tess1d[k].faces[0].tric == NULL) {
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return EGADS_MALLOC;
      }
      for (i = 0; i < n; i++)
        mtess->tess1d[k].faces[0].tric[i] = btess->tess1d[j].faces[0].tric[i];
    }

    n = mtess->tess1d[k].faces[1].nface;
    if (n > 1) {
      mtess->tess1d[k].faces[1].faces = (int *) EG_alloc(n*sizeof(int));
      if (mtess->tess1d[k].faces[1].faces == NULL) {
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return EGADS_MALLOC;
      }
      for (i = 0; i < n; i++) {
        i0 = btess->tess1d[j].faces[1].faces[i];
        if ((i0 != 0) && (fMap != NULL)) i0 = fMap[i0-1];
        mtess->tess1d[k].faces[1].faces[i] = i0;
      }
    }
    n *= mtess->tess1d[k].npts-1;
    if (n > 1) {
      mtess->tess1d[k].faces[1].tric = (int *) EG_alloc(n*sizeof(int));
      if (mtess->tess1d[k].faces[1].tric == NULL) {
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return EGADS_MALLOC;
      }
      for (i = 0; i < n; i++)
        mtess->tess1d[k].faces[1].tric[i] = btess->tess1d[j].faces[1].tric[i];
    }
    n = mtess->tess1d[k].npts;
    mtess->tess1d[k].xyz = (double *) EG_alloc(3*n*sizeof(double));
    mtess->tess1d[k].t   = (double *) EG_alloc(  n*sizeof(double));
    if ((mtess->tess1d[k].xyz == NULL) || (mtess->tess1d[k].xyz == NULL)) {
      EG_free(faces);
      EG_free(nodes);
      EG_cleanupTess(mtess);
      EG_free(mtess);
      return EGADS_MALLOC;
    }
    
    stat = EG_getTopology(mtess->tess1d[k].obj, &geom, &oclass, &mtype, range,
                          &nc, &child, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: getTopology Edge %d = %d (EG_mapTessBody)!\n",
               k+1, stat);
      EG_free(faces);
      EG_free(nodes);
      EG_cleanupTess(mtess);
      EG_free(mtess);
      return stat;
    }
    tr = btess->tess1d[j].t[n-1] - btess->tess1d[j].t[0];
    for (i = 0; i < n; i++) {
      param                 = (btess->tess1d[j].t[i]-btess->tess1d[j].t[0])/tr;
      mtess->tess1d[k].t[i] = range[0] + param*(range[1]-range[0]);
    }
    stat = EG_getTopology(child[0], &geom, &oclass, &mtype,
                          &mtess->tess1d[k].xyz[0], &i, &edges, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: getTopology N0 Edge %d = %d (EG_mapTessBody)!\n",
               j+1, stat);
      EG_free(faces);
      EG_free(nodes);
      EG_cleanupTess(mtess);
      EG_free(mtess);
      return stat;
    }
    for (i = 1; i < n-1; i++) {
      stat = EG_evaluate(mtess->tess1d[k].obj, &mtess->tess1d[k].t[i], result);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: evaluate %d Edge %d/%d = %d (EG_mapTessBody)!\n",
                 i+1, j+1, k+1, stat);
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return stat;
      }
      mtess->tess1d[k].xyz[3*i  ] = result[0];
      mtess->tess1d[k].xyz[3*i+1] = result[1];
      mtess->tess1d[k].xyz[3*i+2] = result[2];
    }
    stat = EG_getTopology(child[nc-1], &geom, &oclass, &mtype,
                          &mtess->tess1d[k].xyz[3*n-3], &i, &edges, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: getTopology N1 Edge %d = %d (EG_mapTessBody)!\n",
               j+1, stat);
      EG_free(faces);
      EG_free(nodes);
      EG_cleanupTess(mtess);
      EG_free(mtess);
      return stat;
    }
    /* check for matching relative arcLengths */
    if (n > 2) {
      s = (double *) malloc(2*n*sizeof(double));
      if (s != NULL) {
        range[0] = range[1] = 0.0;
        for (i = 0; i < n-1; i++) {
          s[2*i  ] = sqrt((mtess->tess1d[k].xyz[3*i+3]-
                           mtess->tess1d[k].xyz[3*i  ])*
                          (mtess->tess1d[k].xyz[3*i+3]-
                           mtess->tess1d[k].xyz[3*i  ]) +
                          (mtess->tess1d[k].xyz[3*i+4]-
                           mtess->tess1d[k].xyz[3*i+1])*
                          (mtess->tess1d[k].xyz[3*i+4]-
                           mtess->tess1d[k].xyz[3*i+1]) +
                          (mtess->tess1d[k].xyz[3*i+5]-
                           mtess->tess1d[k].xyz[3*i+2])*
                          (mtess->tess1d[k].xyz[3*i+5]-
                           mtess->tess1d[k].xyz[3*i+2]));
          s[2*i+1] = sqrt((btess->tess1d[j].xyz[3*i+3]-
                           btess->tess1d[j].xyz[3*i  ])*
                          (btess->tess1d[j].xyz[3*i+3]-
                           btess->tess1d[j].xyz[3*i  ]) +
                          (btess->tess1d[j].xyz[3*i+4]-
                           btess->tess1d[j].xyz[3*i+1])*
                          (btess->tess1d[j].xyz[3*i+4]-
                           btess->tess1d[j].xyz[3*i+1]) +
                          (btess->tess1d[j].xyz[3*i+5]-
                           btess->tess1d[j].xyz[3*i+2])*
                          (btess->tess1d[j].xyz[3*i+5]-
                           btess->tess1d[j].xyz[3*i+2]));
          range[0] += s[2*i  ];
          range[1] += s[2*i+1];
        }
        if (range[0] > 0.0) {
          for (i0 = i = 0; i < n-1; i++)
            if (fabs(s[2*i  ]/range[0] - s[2*i+1]/range[1]) > 1.e-7) i0++;
          /* do we need to remap the discretization? */
          if (i0 != 0) {
/*          printf(" Edge %d: nseg = %d/%d\n", j+1, i0, n-1);  */
            EG_mapTessTs(btess->tess1d[j], mtess->tess1d[k]);
          }
        }
        free(s);
      }
    }

  }
  
  /* fill the Faces */

  if (mtess->tess2d != NULL)
    for (j = 0; j < mtess->nFace; j++) {
      
      k = j;
      if (fMap != NULL) k = fMap[j] - 1;
      
      /* check BSpline surfaces */
      stat = EG_getTopology(faces[k], &geom, &oclass, &mtype, range, &i, &child,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTopology Face %d = %d (EG_mapTessBody)!\n",
                 k+1, stat);
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return stat;
      }
      if (geom->mtype == BSPLINE) {
        int      *ivec = NULL, *ivecb = NULL;
        double   *rvec = NULL, *rvecb = NULL, sknot, sknotb;
        egObject *rGeom, *geomm, **faceb;
        
        OK = -1;
        EG_getGeometry(geom, &oclass, &mtype, &rGeom, &ivec, &rvec);
        stat = EG_getBodyTopos(btess->src, NULL, FACE, &n, &faceb);
        if (stat == EGADS_SUCCESS) {
          tr   = 0.0;
          stat = EG_getTopology(faceb[j], &geomm, &oclass, &mtype, range, &i,
                                &child, &senses);
          EG_free(faceb);
          if (stat == EGADS_SUCCESS) {
            if (geomm->mtype == BSPLINE) {
              if (outLevel > 1)
                printf(" EGADS Info: Face %d both BSpline (EG_mapTessBody)!\n",
                       j+1);
              EG_getGeometry(geomm, &oclass, &mtype, &rGeom, &ivecb, &rvecb);
              if ((ivec != NULL) && (ivecb != NULL) &&
                  (rvec != NULL) && (rvecb != NULL))
                if ((ivec[3] == ivecb[3]) && (ivec[6] == ivecb[6])) {
                  for (i = 0; i < ivec[3]; i++) {
                    sknot  = (rvec[i] -rvec[0]) /(rvec[ivec[3]-1] -rvec[0]);
                    sknotb = (rvecb[i]-rvecb[0])/(rvecb[ivec[3]-1]-rvecb[0]);
                    if (fabs(sknot-sknotb) > KNDIFF) {
                      if (fabs(sknot-sknotb) > tr) tr = fabs(sknot-sknotb);
                      OK = 1;
                    }
                  }
                  for (i = ivec[3]; i < ivec[3]+ivec[6]; i++) {
                    sknot  = (rvec[i] -rvec[ivec[3]]) /
                             (rvec[ivec[3]+ivec[6]-1] -rvec[ivec[3]]);
                    sknotb = (rvecb[i]-rvecb[ivec[3]])/
                             (rvecb[ivec[3]+ivec[6]-1]-rvecb[ivec[3]]);
                    if (fabs(sknot-sknotb) > KNDIFF) {
                      if (fabs(sknot-sknotb) > tr) tr = fabs(sknot-sknotb);
                      if (OK == -1) {
                        OK = 2;
                      } else {
                        if (OK == 1) OK = 3;
                      }
                    }
                  }
                  if (OK == -1) OK = 0;
                } else {
                  printf(" EGADS Warning: Face %d BSpline #knots (EG_mapTessBody)!\n",
                         j+1);
                }
              EG_free(ivecb);
              EG_free(rvecb);
            } else {
              printf(" EGADS Warning: Face %d BSpline Mismatch (EG_mapTessBody)!\n",
                     j+1);
            }
          }
          if (OK == -1) {
            OK = 4;
          } else {
            if (OK != 0)
              printf(" EGADS Info: Face %d BSp knots OK=%d  %le (EG_mapTessBody)!\n",
                     j+1, OK, tr);
          }
#ifdef INSERTKNOTS
          if ((OK > 0) && (OK < 4)) {
            stat = EG_mapSequen(geomm, geom, &btess->tess2d[k].mKnots);
            if (stat != EGADS_SUCCESS)
              printf(" EGADS Error: mapSequen Face %d = %d (EG_mapTessBody)!\n",
                     j+1, stat);
          }
#endif
        }
        EG_free(ivec);
        EG_free(rvec);
      }
      
      /* compute the barycentric coordinates if not there */
      if (btess->tess2d[j].bary == NULL) {
        stat = EG_baryFrame(&btess->tess2d[j]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: baryFrame Face %d = %d (EG_mapTessBody)!\n",
                   j+1, stat);
          EG_free(faces);
          EG_free(nodes);
          EG_cleanupTess(mtess);
          EG_free(mtess);
          return stat;
        }
      }
      
      /* allocate our space */
      xyz    = (double *) EG_alloc(3*mtess->tess2d[k].npts*sizeof(double));
      uv     = (double *) EG_alloc(2*mtess->tess2d[k].npts*sizeof(double));
      ptype  = (int *)    EG_alloc(  mtess->tess2d[k].npts*sizeof(int));
      pindex = (int *)    EG_alloc(  mtess->tess2d[k].npts*sizeof(int));
      bary   = (egBary *) EG_alloc(  mtess->tess2d[k].npts*sizeof(egBary));
      frame  = (int *)    EG_alloc(3*mtess->tess2d[k].nframe*sizeof(int));
      frlps  = (int *)    EG_alloc(  mtess->tess2d[k].nfrlps*sizeof(int));
      tris   = (int *)    EG_alloc(3*mtess->tess2d[k].ntris* sizeof(int));
      tric   = (int *)    EG_alloc(3*mtess->tess2d[k].ntris* sizeof(int));
      if ((xyz    == NULL) || (uv   == NULL) || (ptype == NULL) ||
          (pindex == NULL) || (tris == NULL) || (tric  == NULL) ||
          (frame  == NULL) || (bary == NULL) || (frlps == NULL)) {
        printf(" EGADS Error: Can't Alloc Tess Memory for %d (EG_mapTessBody)!\n",
               k+1);
        if (bary   != NULL) EG_free(bary);
        if (tric   != NULL) EG_free(tric);
        if (tris   != NULL) EG_free(tris);
        if (frlps  != NULL) EG_free(frlps);
        if (frame  != NULL) EG_free(frame);
        if (pindex != NULL) EG_free(pindex);
        if (ptype  != NULL) EG_free(ptype);
        if (uv     != NULL) EG_free(uv);
        if (xyz    != NULL) EG_free(xyz);
        EG_free(faces);
        EG_free(nodes);
        EG_cleanupTess(mtess);
        EG_free(mtess);
        return EGADS_MALLOC;
      }
      
      /* copy the things that don't change */
      for (i = 0; i < mtess->tess2d[k].npts; i++) {
        ptype[i]    = btess->tess2d[j].ptype[i];
        pindex[i]   = btess->tess2d[j].pindex[i];
        bary[i]     = btess->tess2d[j].bary[i];
        if (ptype[i] == 0) {
          if (nMap != NULL) pindex[i] = nMap[pindex[i]-1];
        } else if (ptype[i] > 0) {
          if (eMap != NULL) pindex[i] = eMap[pindex[i]-1];
        }
      }
      for (i = 0; i < 3*mtess->tess2d[k].nframe; i++)
        frame[i]    = btess->tess2d[j].frame[i];
      for (i = 0; i < mtess->tess2d[k].nfrlps; i++)
        frlps[i]    = btess->tess2d[j].frlps[i];
      for (i = 0; i < mtess->tess2d[k].ntris; i++) {
        tris[3*i  ] = btess->tess2d[j].tris[3*i  ];
        tris[3*i+1] = btess->tess2d[j].tris[3*i+1];
        tris[3*i+2] = btess->tess2d[j].tris[3*i+2];
        tric[3*i  ] = btess->tess2d[j].tric[3*i  ];
        tric[3*i+1] = btess->tess2d[j].tric[3*i+1];
        tric[3*i+2] = btess->tess2d[j].tric[3*i+2];
        if ((tric[3*i  ] < 0) && (eMap != 0)) tric[3*i  ] = -eMap[-tric[3*i  ]-1];
        if ((tric[3*i+1] < 0) && (eMap != 0)) tric[3*i+1] = -eMap[-tric[3*i+1]-1];
        if ((tric[3*i+2] < 0) && (eMap != 0)) tric[3*i+2] = -eMap[-tric[3*i+2]-1];
      }
      mtess->tess2d[k].ptype  = ptype;
      mtess->tess2d[k].pindex = pindex;
      mtess->tess2d[k].bary   = bary;
      mtess->tess2d[k].frame  = frame;
      mtess->tess2d[k].frlps  = frlps;
      mtess->tess2d[k].tris   = tris;
      mtess->tess2d[k].tric   = tric;
      
      /* fill in the frame data */
      
      for (last = iloop = i = 0;
           i < mtess->tess2d[k].frlps[mtess->tess2d[k].nfrlps-1]; i++) {
        /* next loop? */
        if (i == mtess->tess2d[k].frlps[iloop]) {
          last = i;
          iloop++;
        }
        if (ptype[i] == 0) {
          /* a node */
          stat = EG_getTopology(nodes[pindex[i]-1], &geom, &oclass, &mtype,
                                &xyz[3*i], &nc, &edges, &senses);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: getTopo N %d Face %d = %d (EG_mapTessBody)!\n",
                     pindex[i], k+1, stat);
            EG_free(uv);
            EG_free(xyz);
            EG_free(faces);
            EG_free(nodes);
            EG_cleanupTess(mtess);
            EG_free(mtess);
            return stat;
          }
          /* set UV -- if Edge in twice, don't use & look for other */
          stat = EG_getBodyTopos(body, nodes[pindex[i]-1], EDGE, &nc, &edges);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: getBodyTopo E Face %d = %d (EG_mapTessBody)!\n",
                     k+1, stat);
            EG_free(uv);
            EG_free(xyz);
            EG_free(faces);
            EG_free(nodes);
            EG_cleanupTess(mtess);
            EG_free(mtess);
            return stat;
          }
          for (sen = n = 0; n < nc; n++) {
            i2 = EG_indexBodyTopo(body, edges[n]);
            i0 = EG_faceConnIndex(mtess->tess1d[i2-1].faces[0], k+1);
            i1 = EG_faceConnIndex(mtess->tess1d[i2-1].faces[1], k+1);
            if ((i0 == 0) && (i1 == 0)) continue;
            if ((i0 != 0) && (i1 != 0)) continue;
            ifrm = 0;
            if (mtess->tess1d[i2-1].nodes[0] == mtess->tess1d[i2-1].nodes[1]) {
              lst = i-1;
              if (lst < last) lst = mtess->tess2d[k].frlps[iloop]-1;
              nxt = i+1;
              if (nxt == mtess->tess2d[k].frlps[iloop]) nxt = last;
              if ((ptype[lst] > 0) && (pindex[lst] == i2)) {
                if (ptype[lst] != 2) ifrm = mtess->tess1d[i2-1].npts-1;
              } else if ((ptype[nxt] > 0) && (pindex[nxt] == i2)) {
                if (ptype[nxt] != 2) ifrm = mtess->tess1d[i2-1].npts-1;
              } else {
                if (edges[n]->mtype != DEGENERATE)
                  printf(" EGADS Info: Can't find E=%d dir F=%d (EG_mapTessBody)!\n",
                         i2, k+1);
/*               printf("   pointer = %ld   npts = %d\n",
                        mtess->tess1d[i2-1].t, mtess->tess1d[i2-1].npts);  */
              }
            } else if (mtess->tess1d[i2-1].nodes[1] == pindex[i]) {
              ifrm = mtess->tess1d[i2-1].npts-1;
            }
            stat = EG_getEdgeUV(faces[k], edges[n], 0,
                                mtess->tess1d[i2-1].t[ifrm], &uv[2*i]);
/*          if (edges[n]->mtype == DEGENERATE)
              printf(" EGADS Info: Degen UV = %lf %lf\n", uv[2*i], uv[2*i+1]); */
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: getEdgeUV Node %d = %d (EG_mapTessBody)!\n",
                       pindex[i], stat);
              EG_free(edges);
              EG_free(uv);
              EG_free(xyz);
              EG_free(faces);
              EG_free(nodes);
              EG_cleanupTess(mtess);
              EG_free(mtess);
              return stat;
            }
            
            /* are we at a degenerate node? */
            result[3] = result[4] = result[5] = 1.0;
            result[6] = result[7] = result[8] = 1.0;
            EG_evaluate(faces[k], &uv[2*i], result);
            if (sqrt(result[3]*result[3] + result[4]*result[4] +
                     result[5]*result[5]) < DEGENUV)
              uv[2*i  ] = btess->tess2d[j].uv[2*i  ];
            if (sqrt(result[6]*result[6] + result[7]*result[7] +
                     result[8]*result[8]) < DEGENUV)
              uv[2*i+1] = btess->tess2d[j].uv[2*i+1];
/*          printf(" Node %d: du = %le   dv = %le\n", pindex[i],
                   sqrt(result[3]*result[3]+result[4]*result[4]+ result[5]*result[5]),
                   sqrt(result[6]*result[6]+result[7]*result[7]+ result[8]*result[8]));
 */
            sen++;
            break;
          }
          EG_free(edges);
          if (sen == 0)
            printf(" EGADS Info: Can't Find Node %d dir F=%d (EG_mapTessBody)!\n",
                   pindex[i], k+1);
          
        } else if (ptype[i] > 0) {
          
          /* an edge vertex */
          stat = EG_evaluate( mtess->tess1d[pindex[i]-1].obj,
                             &mtess->tess1d[pindex[i]-1].t[ptype[i]-1], result);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: evaluate %d EDGE %d = %d (EG_mapTessBody)!\n",
                     ptype[i], pindex[i], stat);
            EG_free(uv);
            EG_free(xyz);
            EG_free(faces);
            EG_free(nodes);
            EG_cleanupTess(mtess);
            EG_free(mtess);
            return stat;
          }
          xyz[3*i  ] = result[0];
          xyz[3*i+1] = result[1];
          xyz[3*i+2] = result[2];
          
          i0  = EG_faceConnIndex(mtess->tess1d[pindex[i]-1].faces[0], k+1);
          i1  = EG_faceConnIndex(mtess->tess1d[pindex[i]-1].faces[1], k+1);
          sen = 0;
          if ((i0 != 0) && (i1 != 0)) {
            nxt = i+1;
            if (nxt == mtess->tess2d[k].frlps[iloop]) nxt = last;
            if (pindex[i] == pindex[nxt]) {
              sen = ptype[nxt] - ptype[i];
            } else {
              if (mtess->tess1d[pindex[i]-1].nodes[0] ==
                  mtess->tess1d[pindex[i]-1].nodes[1]) {
                /* potential problem with periodic Edge and 3 pts! */
                if (ptype[i] == 2) {
                  sen = -1;
                } else {
                  sen =  1;
                }
              } else if (pindex[nxt] == mtess->tess1d[pindex[i]-1].nodes[0]) {
                sen = -1;
              } else if (pindex[nxt] == mtess->tess1d[pindex[i]-1].nodes[1]) {
                sen =  1;
              } else {
                printf(" EGADS Info: Can't Edge %d dir F=%d (EG_mapTessBody)!\n",
                       pindex[i], k+1);
              }
            }
/*          printf("  face = %d, sense = %d\n", k+1, faces[k]->mtype);  */
            sen *= faces[k]->mtype;
          }
          stat = EG_getEdgeUV(faces[k], mtess->tess1d[pindex[i]-1].obj, sen,
                              mtess->tess1d[pindex[i]-1].t[ptype[i]-1], &uv[2*i]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: getEdgeUV %d EDGE %d = %d (EG_mapTessBody)!\n",
                     ptype[i], pindex[i], stat);
            EG_free(uv);
            EG_free(xyz);
            EG_free(faces);
            EG_free(nodes);
            EG_cleanupTess(mtess);
            EG_free(mtess);
            return stat;
          }
        } else {
          printf(" EGADS Info: Bad Frame Vert %d (%d %d) F=%d (EG_mapTessBody)!\n",
                 i+1, ptype[i], pindex[i], k+1);
        }
      }
      
      /* fill in the interior uvs & xyzs */
      obj2D = faces[k];
#ifdef INSERTKNOTS
      if (mtess->tess2d[k].mKnots != NULL) obj2D = mtess->tess2d[k].mKnots;
#endif
      for (i = mtess->tess2d[k].frlps[mtess->tess2d[k].nfrlps-1];
           i < mtess->tess2d[k].npts; i++) {
        ifrm = bary[i].tri - 1;
        i0   = frame[3*ifrm  ] - 1;
        i1   = frame[3*ifrm+1] - 1;
        i2   = frame[3*ifrm+2] - 1;
        if ((ptype[i0] == -1) || (ptype[i1] == -1) || (ptype[i2] == -1)) {
          if (outLevel > 0)
            printf(" EGADS Error: %d Frame types = %d (%d)  %d (%d)  %d (%d)!\n",
                   k+1, ptype[i0], i0+1, ptype[i1], i1+1, ptype[i2], i2+1);
          EG_free(uv);
          EG_free(xyz);
          EG_free(faces);
          EG_free(nodes);
          EG_cleanupTess(mtess);
          EG_free(mtess);
          return EGADS_INDEXERR;
        }
        uv[2*i  ] = uv[2*i0  ]*bary[i].w[0] + uv[2*i1  ]*bary[i].w[1] +
                    uv[2*i2  ]*(1.0-bary[i].w[0]-bary[i].w[1]);
        uv[2*i+1] = uv[2*i0+1]*bary[i].w[0] + uv[2*i1+1]*bary[i].w[1] +
                    uv[2*i2+1]*(1.0-bary[i].w[0]-bary[i].w[1]);
        stat = EG_evaluate(obj2D, &uv[2*i], result);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: evaluate %d Face %d = %d (EG_mapTessBody)!\n",
                   i+1, k+1, stat);
          EG_free(uv);
          EG_free(xyz);
          EG_free(faces);
          EG_free(nodes);
          EG_cleanupTess(mtess);
          EG_free(mtess);
          return stat;
        }
        xyz[3*i  ] = result[0];
        xyz[3*i+1] = result[1];
        xyz[3*i+2] = result[2];
      }
      mtess->tess2d[k].xyz = xyz;
      mtess->tess2d[k].uv  = uv;
    }
  EG_free(faces);
  EG_free(nodes);
  
  /* create the object */

  stat = EG_makeObject(context, &mapObj);
  if (stat != EGADS_SUCCESS) {
    EG_cleanupTess(mtess);
    EG_free(mtess);
    return stat;
  }
  mapObj->oclass = TESSELLATION;
  mapObj->blind  = mtess;
  EG_referenceObject(mapObj, context);
  EG_referenceTopObj(body,   mapObj);
  *mapTess = mapObj;
 
  return EGADS_SUCCESS;
}


int
EG_locateTessBody(const egObject *tess, int npts, const int *ifaces,
                  const double *uvs, /*@null@*/ int *itris, double *results)
{
  int          i, iface, stat, nface, aType, alen;
  double       data[18];
  egObject     *tessb, *obj2D, **faces;
  egTessel     *btess;
  const int    *ints = NULL;
  const double *reals;
  const char   *str;
  
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;
  btess = tess->blind;
  tessb = btess->src;
  if (tessb == NULL)                return EGADS_NULLOBJ;
  if (tessb->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if (tessb->oclass != BODY)        return EGADS_NOTBODY;
  if (btess->tess2d == NULL)        return EGADS_NODATA;

  stat = EG_attributeRet(tessb, ".fMap", &aType, &alen, &ints, &reals, &str);
  if (stat == EGADS_SUCCESS) {
    if (aType != ATTRINT) {
      printf(" EGADS Error: Map Attribute Not an Int (EG_locateTessBody)!\n");
      return EGADS_ATTRERR;
    }
    if (alen != btess->nFace) {
      printf(" EGADS Error: Len Mismatch %d %d (EG_locateTessBody)!\n",
             alen, btess->nFace);
      return EGADS_TOPOERR;
    }
  }

  if (itris == NULL) {
    stat = EG_getBodyTopos(tessb, NULL, FACE, &nface, &faces);
    if (stat != EGADS_SUCCESS) return stat;

    for (i = 0; i < npts; i++) {
      iface = abs(ifaces[i]);
      if ((iface < 1) || (iface > btess->nFace)) {
        printf(" EGADS Error: %d = %d [1-%d] (EG_locateTessBody)!\n",
               i+1, iface, btess->nFace);
        EG_free(faces);
        return EGADS_INDEXERR;
      }
      if ((stat == EGADS_SUCCESS) && (ifaces[i] < 0)) {
        if (ints == NULL) {
          printf(" EGADS Error: %d Mapped w/ No Mapping (EG_locateTessBody)!\n",
                 i+1);
          EG_free(faces);
          return EGADS_INDEXERR;
        }
        iface = ints[iface-1];
        if ((iface < 1) || (iface > btess->nFace)) {
          printf(" EGADS Error: Mapped %d = %d [1-%d] (EG_locateTessBody)!\n",
                 i+1, iface, btess->nFace);
          EG_free(faces);
          return EGADS_INDEXERR;
        }
      }
      obj2D = faces[iface-1];
      if (btess->tess2d[iface-1].mKnots != NULL)
        obj2D = btess->tess2d[iface-1].mKnots;
      stat = EG_evaluate(obj2D, &uvs[2*i], data);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        return stat;
      }
      results[3*i  ] = data[0];
      results[3*i+1] = data[1];
      results[3*i+2] = data[2];
    }
    EG_free(faces);
    
    return EGADS_SUCCESS;
  }

  for (i = 0; i < npts; i++) {
    iface = abs(ifaces[i]);
    if ((iface < 1) || (iface > btess->nFace)) {
      printf(" EGADS Error: %d = %d [1-%d] (EG_locateTessBody)!\n",
             i+1, iface, btess->nFace);
      return EGADS_INDEXERR;
    }
    if ((stat == EGADS_SUCCESS) && (ifaces[i] < 0)) {
      if (ints == NULL) {
        printf(" EGADS Error: %d Mapped w/ No Mapping (EG_locateTessBody)!\n",
               i+1);
        return EGADS_INDEXERR;
      }
      iface = ints[iface-1];
      if ((iface < 1) || (iface > btess->nFace)) {
        printf(" EGADS Error: Mapped %d = %d [1-%d] (EG_locateTessBody)!\n",
               i+1, iface, btess->nFace);
        return EGADS_INDEXERR;
      }
    }
    itris[i] = EG_baryTess(btess->tess2d[iface-1], &uvs[2*i], &results[3*i]);
  }

  return EGADS_SUCCESS;
}


int
EG_getTessQuads(const egObject *tess, int *nquad, int **fIndices)
{
  int      i, n, outLevel, *ivec;
  egTessel *btess;
  egObject *obj;
  
  *nquad    = 0;
  *fIndices = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_getTessQuads)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_getTessQuads)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_getTessQuads)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_getTessQuads)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_getTessQuads)!\n");
    return EGADS_NODATA;  
  }

  for (n = i = 0; i < btess->nFace; i++)
    if (btess->tess2d[i+btess->nFace].xyz != NULL) n++;
  if (n == 0) return EGADS_SUCCESS;
  
  ivec = (int *) EG_alloc(n*sizeof(int));
  if (ivec == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d Faces (EG_getTessQuads)!\n", n);
    return EGADS_MALLOC;
  }
  
  for (n = i = 0; i < btess->nFace; i++)
    if (btess->tess2d[i+btess->nFace].xyz != NULL) {
      ivec[n] = i+1;
      n++;
    }
  *nquad    = n;
  *fIndices = ivec;
  
  return EGADS_SUCCESS;
}


static int
EG_quadLoop(egTessel *btess, int outLevel, int nedge, int *eindex, 
            int *senses, double *parms, int *lim)
{
  int    i, j, ie0, ie1, imax, nside;
  double t0[3], t1[3], dist, dmax, edgeTOL = 0.05;
  
  if ((parms[0] >= 0.001) && (parms[0] <= 0.5)) edgeTOL = parms[0];
  nside = nedge;

  while (nside > 4) {
  
    /* merge the 2 Edges with the smallest delta in tangent */
    dmax = -1.0;
    imax = -1; 
    for (i = 0; i < 4; i++) {
      ie0 = eindex[lim[i]  ]-1;
      ie1 = eindex[lim[i]+1]-1;
      if (senses[lim[i]  ] == 1) {
        j     = btess->tess1d[ie0].npts-2;
        t0[0] = btess->tess1d[ie0].xyz[3*j+3]-btess->tess1d[ie0].xyz[3*j  ];
        t0[1] = btess->tess1d[ie0].xyz[3*j+4]-btess->tess1d[ie0].xyz[3*j+1];
        t0[2] = btess->tess1d[ie0].xyz[3*j+5]-btess->tess1d[ie0].xyz[3*j+2];
      } else {
        t0[0] = btess->tess1d[ie0].xyz[0]-btess->tess1d[ie0].xyz[3];
        t0[1] = btess->tess1d[ie0].xyz[1]-btess->tess1d[ie0].xyz[4];
        t0[2] = btess->tess1d[ie0].xyz[2]-btess->tess1d[ie0].xyz[5];
      }
      dist = sqrt(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]);
      if (dist != 0.0) {
        t0[0] /= dist;
        t0[1] /= dist;
        t0[2] /= dist;
      }
      if (senses[lim[i]+1] == 1) {
        t1[0] = btess->tess1d[ie1].xyz[3]-btess->tess1d[ie1].xyz[0];
        t1[1] = btess->tess1d[ie1].xyz[4]-btess->tess1d[ie1].xyz[1];
        t1[2] = btess->tess1d[ie1].xyz[5]-btess->tess1d[ie1].xyz[2];
      } else {
        j     = btess->tess1d[ie1].npts-2;
        t1[0] = btess->tess1d[ie1].xyz[3*j  ]-btess->tess1d[ie1].xyz[3*j+3];
        t1[1] = btess->tess1d[ie1].xyz[3*j+1]-btess->tess1d[ie1].xyz[3*j+4];
        t1[2] = btess->tess1d[ie1].xyz[3*j+2]-btess->tess1d[ie1].xyz[3*j+5];
      }
      dist = sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2]);
      if (dist != 0.0) {
        t1[0] /= dist;
        t1[1] /= dist;
        t1[2] /= dist;
      }
      dist = t0[0]*t1[0] + t0[1]*t1[1] + t0[2]*t1[2];
      if (outLevel > 1)
        printf("  Dot between %d %d = %lf\n", ie0+1, ie1+1, dist);
      if (dist > dmax) {
        dmax = dist;
        imax = i;
      }
    }
    if (imax == -1) return EGADS_INDEXERR;
    if (dmax < 1.0-edgeTOL) return EGADS_INDEXERR;

    for (i = imax; i < 3; i++) lim[i] = lim[i+1];
    lim[3]++;
    nside--;
    if (outLevel > 1)
      printf("  endIndex = %d %d %d %d,  nSide = %d\n", 
             lim[0], lim[1], lim[2], lim[3], nside);
  }

  return EGADS_SUCCESS;
}


int
EG_makeQuads(egObject *tess, double *parms, int index)
{
  int      i, j, k, l, m, n, outLevel, stat, oclass, mtype, ftype;
  int      nface, nloop, nedge, sens, npt, nx, *eindex, lim[4];
  int      npts, npat, save, iv, iv1, nside, pats[34], lens[4];
  int      *ptype, *pindex, *pin, *vpats, *senses, *ntable;
  double   *uvs, *quv, *xyz, *xyzs, limits[4], res[18], area;
  connect  *etable;
  egTessel *btess;
  egPatch  *patch;
  egObject *obj, *geom, **faces, **loops, **edges;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_makeQuads)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_makeQuads)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_makeQuads)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_makeQuads)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_makeQuads)!\n");
    return EGADS_NODATA;  
  }
  if ((index < 1) || (index > btess->nFace)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_makeQuads)!\n",
             index, btess->nFace);
    return EGADS_INDEXERR;
  }

  /* quad patch based on current Edge tessellations */
  
  stat = EG_getBodyTopos(obj, NULL, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getTopology(faces[index-1], &geom, &oclass, &ftype, limits,
                        &nloop, &loops, &senses);
  if (stat != EGADS_SUCCESS) {
    EG_free(faces);
    return stat;
  }
  if (nloop != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Face %d has %d loops (EG_makeQuads)!\n",
             index, nloop);
    EG_free(faces);
    return EGADS_TOPOERR;
  }
  stat = EG_getTopology(loops[0], &geom, &oclass, &mtype, limits,
                        &nedge, &edges, &senses);
  if (stat != EGADS_SUCCESS) {
    EG_free(faces);
    return stat;
  }
  if (nedge < 4) {
    if (outLevel > 0)
      printf(" EGADS Error: %d Edges in Face %d (EG_makeQuads)!\n", 
             nedge, index);
    EG_free(faces);
    return EGADS_INDEXERR;
  }

  /* get Edge Indices */
  eindex = (int *) EG_alloc(nedge*sizeof(int));
  if (eindex == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d Edges (EG_makeQuads)!\n",
             nedge);
    EG_free(faces);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nedge; i++) {
    if (edges[i]->mtype == DEGENERATE) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge in Face %d is Degenerate (EG_makeQuads)!\n", 
               index);
      EG_free(eindex);
      EG_free(faces);
      return EGADS_INDEXERR;
    }
    eindex[i] = 0;
    for (j = 0; j < btess->nEdge; j++)
      if (edges[i] == btess->tess1d[j].obj) {
        eindex[i] = j+1;
        break;
      }
    if (eindex[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Not Found in Tess (EG_makeQuads)!\n");
      EG_free(eindex);
      EG_free(faces);
      return EGADS_NOTFOUND;
    }
  }
  
  /* block off the 4 sides if available */
  lim[0] = lens[0] = lens[1] = lens[2] = lens[3] = 0;
  lim[1] = 1;
  lim[2] = 2;
  lim[3] = 3;
  if (nedge > 4) {
    stat = EG_quadLoop(btess, outLevel, nedge, eindex, senses, 
                       parms, lim);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: %d Edges in Face %d (EG_makeQuads)!\n", 
               nedge, index);
      EG_free(eindex);
      EG_free(faces);
      return stat;
    }
  }
  for (npts = l = i = 0; i < nedge; i++) {
    j     = eindex[i]-1;
    npts += btess->tess1d[j].npts-1;
    if (ftype == SFORWARD) {
      lens[l]   += btess->tess1d[j].npts-1;
    } else {
      lens[3-l] += btess->tess1d[j].npts-1;
    }
    if (lim[l] == i) l++;
  }

  /* allocate the info for the frame of the blocking */
  xyzs = (double *) EG_alloc(3*npts*sizeof(double));
  if (xyzs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d XYZs (EG_makeQuads)!\n",
             npts);
    EG_free(eindex);
    EG_free(faces);
    return EGADS_MALLOC;
  }
  uvs = (double *) EG_alloc(2*npts*sizeof(double));
  if (uvs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d Points (EG_makeQuads)!\n",
             npts);
    EG_free(xyzs);
    EG_free(eindex);
    EG_free(faces);
    return EGADS_MALLOC;
  }
  pin = (int *) EG_alloc(3*npts*sizeof(int));
  if (pin == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc on %d Pindex (EG_makeQuads)!\n",
             npts);
    EG_free(uvs);
    EG_free(xyzs);
    EG_free(eindex);
    EG_free(faces);
    return EGADS_MALLOC;
  }
 
  /* fill in uvs around the loop */

  for (npts = i = 0; i < nedge; i++) {
    if (ftype == SFORWARD) {
      j    =  eindex[i] - 1;
      sens =  senses[i];
      m    =  senses[i];
    } else {
      j    =  eindex[nedge-i-1] - 1;
      sens = -senses[nedge-i-1];
      m    =  senses[nedge-i-1];
    }
    if (sens == 1) {
      for (k = 0; k < btess->tess1d[j].npts-1; k++, npts++) {
        stat = EG_getEdgeUV(faces[index-1], btess->tess1d[j].obj, 
                            m, btess->tess1d[j].t[k], &uvs[2*npts]);
        if (stat != EGADS_SUCCESS) {
          EG_free(uvs);
          EG_free(pin);
          EG_free(xyzs);
          EG_free(eindex);
          EG_free(faces);
          return stat;        
        }
        xyzs[3*npts  ] = btess->tess1d[j].xyz[3*k  ];
        xyzs[3*npts+1] = btess->tess1d[j].xyz[3*k+1];
        xyzs[3*npts+2] = btess->tess1d[j].xyz[3*k+2];
        pin[3*npts  ]  =  j+1;
        pin[3*npts+1]  =  k+1;
        pin[3*npts+2]  = -j-1;
        if (k == 0) {
          pin[3*npts  ] = btess->tess1d[j].nodes[0];
          pin[3*npts+1] = 0;
        }
      }
    } else {
      for (k = btess->tess1d[j].npts-1; k >= 1; k--, npts++) {
        stat = EG_getEdgeUV(faces[index-1], btess->tess1d[j].obj, 
                            m, btess->tess1d[j].t[k], &uvs[2*npts]);
        if (stat != EGADS_SUCCESS) {
          EG_free(uvs);
          EG_free(pin);
          EG_free(xyzs);
          EG_free(eindex);
          EG_free(faces);
          return stat;        
        }
        xyzs[3*npts  ] = btess->tess1d[j].xyz[3*k  ];
        xyzs[3*npts+1] = btess->tess1d[j].xyz[3*k+1];
        xyzs[3*npts+2] = btess->tess1d[j].xyz[3*k+2];
        pin[3*npts  ]  =  j+1;
        pin[3*npts+1]  =  k+1;
        pin[3*npts+2]  = -j-1;
        if (k == btess->tess1d[j].npts-1) {
          pin[3*npts  ] = btess->tess1d[j].nodes[1];
          pin[3*npts+1] = 0;
        }
      }
    }
  }
  EG_free(eindex);
  
  i = npts-1;
  area = (uvs[0]+uvs[2*i])*(uvs[1]-uvs[2*i+1]);
  for (i = 0; i < npts-1; i++)
    area += (uvs[2*i+2]+uvs[2*i])*(uvs[2*i+3]-uvs[2*i+1]);
  area /= 2.0;
  if (outLevel > 1) 
    printf(" makeQuads: loop area = %lf,  ori = %d\n", area, ftype);

  stat = EG_quadFill(faces[index-1], parms, lens, uvs, 
                     &npt, &quv, &npat, pats, &vpats);
  EG_free(uvs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: quadFill = %d (EG_makeQuads)!\n",
             stat);
    EG_free(pin);
    EG_free(xyzs);
    EG_free(faces);
    return EGADS_CONSTERR;
  }

  xyz    = (double *) EG_alloc(3*npt*sizeof(double));
  ptype  = (int *)    EG_alloc(  npt*sizeof(int));
  pindex = (int *)    EG_alloc(  npt*sizeof(int));
  if ((xyz == NULL) || (ptype == NULL) || (pindex == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc npts = %d (EG_makeQuads)!\n",
             npt);
    if (pindex != NULL) EG_free(pindex);
    if (ptype  != NULL) EG_free(ptype);
    if (xyz    != NULL) EG_free(xyz);
    EG_free(vpats);
    EG_free(quv);
    EG_free(pin);
    EG_free(xyzs);
    EG_free(faces);
    return EGADS_MALLOC;
  }
  for (i = 0; i < npts; i++) {
    pindex[i]  = pin[3*i  ];
    ptype[i]   = pin[3*i+1];
    xyz[3*i  ] = xyzs[3*i  ];
    xyz[3*i+1] = xyzs[3*i+1];
    xyz[3*i+2] = xyzs[3*i+2];
  }
  EG_free(xyzs);  
  for (i = npts; i < npt; i++) {
    pindex[i] = ptype[i] = -1;
    EG_evaluate(faces[index-1], &quv[2*i], res);
    xyz[3*i  ] = res[0];
    xyz[3*i+1] = res[1];
    xyz[3*i+2] = res[2];
  }
  EG_free(faces);
  patch = (egPatch *) EG_alloc(npat*sizeof(egPatch));
  if (patch == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Malloc npatchs = %d (EG_makeQuads)!\n",
             npat);
    EG_free(pindex);
    EG_free(ptype);
    EG_free(xyz);
    EG_free(vpats);
    EG_free(quv);
    EG_free(pin);
    return EGADS_MALLOC;
  }

  /* put back in face orientation */
  if (ftype != SFORWARD)
    for (iv = k = 0; k < npat; k++) {
      nx = pats[2*k  ];
      for (j = 0; j < pats[2*k+1]; j++) {
        for (i = 0; i < nx/2; i++) {
          m           = nx - i - 1;
          save        = vpats[iv+i];
          vpats[iv+i] = vpats[iv+m];
          vpats[iv+m] = save;
        }
        iv += nx;
      }
    }
    
  for (nx = m = 0; m < npat; m++) 
    nx += 2*pats[2*m] + 2*pats[2*m+1] - 4;
  ntable = (int *) EG_alloc(npt*sizeof(int));
  if (ntable == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Vert Table Malloc (EG_makeQuads)!\n");
    EG_free(patch);
    EG_free(pindex);
    EG_free(ptype);
    EG_free(xyz);
    EG_free(vpats);
    EG_free(quv);
    EG_free(pin);
    return EGADS_MALLOC;    
  }
  etable = (connect *) EG_alloc((nx+1)*sizeof(connect));
  if (etable == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge Table Malloc (EG_makeQuads)!\n");
    EG_free(ntable);
    EG_free(patch);
    EG_free(pindex);
    EG_free(ptype);
    EG_free(xyz);
    EG_free(vpats);
    EG_free(quv);
    EG_free(pin);
    return EGADS_MALLOC;    
  }
  for (j = 0; j < npt; j++) ntable[j] = NOTFILLED;

  /* fill in the patch */
  for (k = m = 0; m < npat; m++) {
    if (outLevel > 1)
      printf("  Patch %d: size = %d %d\n", m+1, pats[2*m  ], pats[2*m+1]);
    patch[m].nu     = pats[2*m  ];
    patch[m].nv     = pats[2*m+1];
    patch[m].ipts   = (int *) EG_alloc(pats[2*m]*pats[2*m+1]*sizeof(int));
    patch[m].bounds = (int *) EG_alloc((2*(pats[2*m  ]-1) + 
                                        2*(pats[2*m+1]-1))*sizeof(int));
    for (n = j = 0; j < pats[2*m+1]; j++)
      for (i = 0; i < pats[2*m]; i++, n++, k++) 
        if (patch[m].ipts != NULL) patch[m].ipts[n] = vpats[k] + 1;
  }

  /* connect the patches */
  nside = -1;
  for (j = 0; j < npts-1; j++)
    EG_makeConnect(j+1, j+2, &pin[3*j+2], 
                   &nside, ntable, etable, index);
  EG_makeConnect(npts, 1, &pin[3*npts-1], 
                 &nside, ntable, etable, index);
  for (l = n = m = 0; m < npat; m++) {
    if (patch[m].bounds == NULL) continue;
    for (k = i = 0; i < patch[m].nu-1; i++, k++) {
      iv  = vpats[l+i  ] + 1;
      iv1 = vpats[l+i+1] + 1;
      patch[m].bounds[k] = n+i+1;
      EG_makeConnect(iv, iv1, &patch[m].bounds[k], 
                     &nside, ntable, etable, index);
    }
    for (i = 0; i < patch[m].nv-1; i++, k++) {
      iv  = vpats[l+(i+1)*patch[m].nu-1] + 1;
      iv1 = vpats[l+(i+2)*patch[m].nu-1] + 1;
      patch[m].bounds[k] = n+(i+1)*(patch[m].nu-1);
      EG_makeConnect(iv, iv1, &patch[m].bounds[k], 
                     &nside, ntable, etable, index);
    }
    for (i = 0; i < patch[m].nu-1; i++, k++) {
      iv  = vpats[l+patch[m].nu*patch[m].nv-i-1] + 1;
      iv1 = vpats[l+patch[m].nu*patch[m].nv-i-2] + 1;
      patch[m].bounds[k] = n+(patch[m].nu-1)*(patch[m].nv-1)-i;
      EG_makeConnect(iv, iv1, &patch[m].bounds[k], 
                     &nside, ntable, etable, index);
    }
    for (i = 0; i < patch[m].nv-1; i++, k++) {
      iv  = vpats[l+(patch[m].nv-i-1)*patch[m].nu] + 1;
      iv1 = vpats[l+(patch[m].nv-i-2)*patch[m].nu] + 1;
      patch[m].bounds[k] = n+(patch[m].nv-i-2)*(patch[m].nu-1);
      EG_makeConnect(iv, iv1, &patch[m].bounds[k], 
                     &nside, ntable, etable, index);
    }
    n += (patch[m].nu-1)*(patch[m].nv-1);
    l +=  patch[m].nu   * patch[m].nv;
  }

  /* report any unconnected boundary sides */
  for (j = 0; j <= nside; j++) {
    if (etable[j].tri == NULL) continue;
    printf(" EGADS Info: Face %d, Unconnected Quad Side %d %d = %d\n",
           index, etable[j].node1, etable[j].node2, *etable[j].tri);
    *etable[j].tri = 0;
  }

  EG_free(etable);
  EG_free(ntable);
  EG_free(vpats);
  EG_free(pin);
  
  /* delete any existing quads */
  EG_deleteQuads(btess, index);

  /* save away the patches */

  i = btess->nFace + index - 1;
  btess->tess2d[i].xyz    = xyz;
  btess->tess2d[i].uv     = quv;
  btess->tess2d[i].ptype  = ptype;
  btess->tess2d[i].pindex = pindex;
  btess->tess2d[i].npts   = npt;
  btess->tess2d[i].patch  = patch;
  btess->tess2d[i].npatch = npat;
  
  return EGADS_SUCCESS;
}


int
EG_getQuads(const egObject *tess, int index, int *len, const double **xyz, 
            const double **uv, const int **ptype, const int **pindex, 
            int *npatch)
{
  int      i, outLevel;
  egTessel *btess;
  egObject *obj;

  *len   = *npatch = 0;
  *xyz   = *uv     = NULL;
  *ptype = *pindex = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_getQuads)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_getQuads)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_getQuads)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_getQuads)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_getQuads)!\n");
    return EGADS_NODATA;  
  }
  if ((index < 1) || (index > btess->nFace)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_getQuads)!\n",
             index, btess->nFace);
    return EGADS_INDEXERR;
  }

  i       = btess->nFace + index - 1;
  *len    = btess->tess2d[i].npts;
  *xyz    = btess->tess2d[i].xyz;
  *uv     = btess->tess2d[i].uv;
  *ptype  = btess->tess2d[i].ptype;
  *pindex = btess->tess2d[i].pindex;
  *npatch = btess->tess2d[i].npatch;

  return EGADS_SUCCESS;
}


int
EG_getPatch(const egObject *tess, int index, int patch, int *nu, int *nv, 
            const int **ipts, const int **bounds)
{
  int      i, outLevel;
  egTessel *btess;
  egObject *obj;

  *nu   = *nv     = 0;
  *ipts = *bounds = NULL;
  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);
  
  btess = (egTessel *) tess->blind;
  if (btess == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Blind Object (EG_getPatch)!\n");  
    return EGADS_NOTFOUND;
  }
  obj = btess->src;
  if (obj == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Source Object (EG_getPatch)!\n");
    return EGADS_NULLOBJ;
  }
  if (obj->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not an Object (EG_getPatch)!\n");
    return EGADS_NOTOBJ;
  }
  if (obj->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Source Not Body (EG_getPatch)!\n");
    return EGADS_NOTBODY;
  }
  if (btess->tess2d == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Face Tessellations (EG_getPatch)!\n");
    return EGADS_NODATA;  
  }
  if ((index < 1) || (index > btess->nFace)) {
    if (outLevel > 0)
      printf(" EGADS Error: Index = %d [1-%d] (EG_getPatch)!\n",
             index, btess->nFace);
    return EGADS_INDEXERR;
  }
  i = btess->nFace + index - 1;
  if (btess->tess2d[i].patch == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: No Patch Data (EG_getPatch)!\n");
    return EGADS_NODATA;
  }
  if ((patch < 1) || (patch > btess->tess2d[i].npatch)) {
    if (outLevel > 0)
      printf(" EGADS Error: Patch index = %d [1-%d] (EG_getPatch)!\n",
             patch, btess->tess2d[i].npatch);
    return EGADS_INDEXERR;
  }

  *nu     = btess->tess2d[i].patch[patch-1].nu;
  *nv     = btess->tess2d[i].patch[patch-1].nv;
  *ipts   = btess->tess2d[i].patch[patch-1].ipts;
  *bounds = btess->tess2d[i].patch[patch-1].bounds;
  
  return EGADS_SUCCESS;
}
