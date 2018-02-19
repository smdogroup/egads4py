/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Blend Derivative Functions
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

//#define SHARPER

#define PI               3.1415926535897931159979635

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])


template<int N>
int EG_spline2dEval(int *ivec, SurrealS<N> *data, const double *uv,
                    SurrealS<N> *point);

template<int N>
int EG_spline2dDeriv(int *ivec, SurrealS<N> *data, int der, const double *uv,
                     SurrealS<N> *deriv);

template<int N>
int EG_spline2dAprx(int endc, int imax, int jmax, const SurrealS<N> *xyz,
                    const SurrealS<N> *uknot,     const SurrealS<N> *vknot,
                    const int *vdata,
                    const SurrealS<N> *wesT,      const SurrealS<N> *easT,
                    const SurrealS<N> *south,           SurrealS<N> *snor,
                    const SurrealS<N> *north,           SurrealS<N> *nnor,
                    double tol, int *header,            SurrealS<N> **rdata);


// also defined in egadsSpline.c
typedef struct {
  int    ncp;                   /* number of control points */
  int    nave;                  /* number used to average positions */
  double *knots;                /* the knot positions */
} egSequ;


typedef struct {
  ego         *objs;
  SurrealS<1> *tbeg;
  SurrealS<1> *tend;
  SurrealS<1> *xyzs;
  int         imax;
  double      *ts;
  int         header[7];
  SurrealS<1> *dataS;
} egStrip;


typedef struct {
  int         magic;
  int         outLevel;
  int         nSec;
  int         planar;
  ego         *secs;
  int         begRC;
  int         endRC;
  SurrealS<1> *rc1;
  SurrealS<1> *rcN;
  SurrealS<1> *vKnots;
  int         *vData;
  int         nStrips;
  egSequ      *sequ;
  egStrip     *strips;
} egCache;



extern "C" int  EG_outLevel( const ego object );
extern "C" int  EG_isPlanar( const ego object );
extern "C" int  EG_allocSeq(int nstripe, egSequ **sequ);
extern "C" void EG_freeSeq(int nstripe, /*@only@*/ egSequ *sequ);
extern "C" int  EG_setSeq(int stripe, egSequ *seq, int num,
                          /*@null@*/ int *iinfo, /*@null@*/ double *rinfo);


extern "C" void
EG_sens_free(/*@null@*/ /*@only@*/ void *cache)
{
  int     i;
  egCache *eCache;

  if (cache == NULL) return;
  eCache = (egCache *) cache;
  if (eCache->magic  != MAGIC) {
    printf(" EGADS ERROR: EG_sens_free called with BAD Cache!\n");
    return;
  }
  
  if (eCache->secs   != NULL) EG_free(eCache->secs);
  if (eCache->rc1    != NULL) EG_free(eCache->rc1);
  if (eCache->rcN    != NULL) EG_free(eCache->rcN);
  if (eCache->vKnots != NULL) EG_free(eCache->vKnots);
  if (eCache->vData  != NULL) EG_free(eCache->vData);
  if (eCache->sequ   != NULL) EG_freeSeq(eCache->nStrips, eCache->sequ);

  if (eCache->strips != NULL) {
    for (i = 0; i < eCache->nStrips+2; i++) {
      if (eCache->strips[i].objs  != NULL) EG_free(eCache->strips[i].objs);
      if (eCache->strips[i].tbeg  != NULL) EG_free(eCache->strips[i].tbeg);
      if (eCache->strips[i].tend  != NULL) EG_free(eCache->strips[i].tend);
      if (eCache->strips[i].xyzs  != NULL) EG_free(eCache->strips[i].xyzs);
      if (eCache->strips[i].ts    != NULL) EG_free(eCache->strips[i].ts);
      if (eCache->strips[i].dataS != NULL) EG_free(eCache->strips[i].dataS);
    }
    EG_free(eCache->strips);
  }

  eCache->magic = 0;
  EG_free(cache);
}


extern "C" int
EG_blend_init(int nsec, const ego *secs,        /*@null@*/ const double *rc1,
              /*@null@*/ const double *rc1_dot, /*@null@*/ const double *rcN,
              /*@null@*/ const double *rcN_dot, int *nStrips, void **cache)
{
  int     i, j, k, n, jj, stat, oclass, mtype, outLevel, nstripe, nchld;
  int     begRC, endRC, planar, *senses, *iinfo;
  double  uvlim[4], *rinfo;
  ego     ref, geom, loop, *chldrn, *edges;
  egCache *eCache;

  *cache   = NULL;
  *nStrips = 0;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  outLevel = EG_outLevel(secs[0]);
  
  /* look at the input and check to see if OK */
  begRC  = endRC = -1;
  planar = 0;
  for (nstripe = i = 0; i < nsec; i++) {
    if (secs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Section Object %d (EG_blend_init)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (secs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is not an EGO (EG_blend_init)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (secs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d has no data (EG_blend_init)!\n", i+1);
      return EGADS_NODATA;
    }
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, uvlim, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d getTopology = %d (EG_blend_init)!\n",
               i+1, stat);
      return stat;
    }
    if (oclass == NODE) {
      if ((i != 0) && (i != nsec-1)) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Node and not Bound (EG_blend_init)!\n",
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
          printf(" EGADS Error: Section %d is Face and not Bound (EG_blend_init)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      if (n != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Face with %d Loops (EG_blend_init)!\n",
                 i+1, n);
        return EGADS_TOPOERR;
      }
      loop = chldrn[0];
      if (ref->mtype != PLANE) planar = 1;
    } else if (oclass == BODY) {
      if (secs[i]->mtype != WIREBODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Not a WireBody (EG_blend_init)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      loop = chldrn[0];
      if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
    } else if (oclass == LOOP) {
      loop = secs[i];
      if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is Not a Loop (EG_blend_init)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
    if (loop == NULL) continue;
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, uvlim, &jj, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Loop getTopology = %d (EG_blend_init)!\n",
               i+1, stat);
      return stat;
    }
    for (n = j = 0; j < jj; j++) {
      stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, uvlim, &k, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Edge %d getTopo = %d (EG_blend_init)!\n",
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
          printf(" EGADS Error: Sec %d has %d Edges -- prev = %d (EG_blend_init)!\n",
                 i+1, n, nstripe);
        return EGADS_TOPOERR;
      }
    }
  }
  if (nstripe == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Edges found (EG_blend_init)!\n");
    return EGADS_TOPOERR;    
  }
  if (((begRC == 1) || (endRC == 1)) && (nsec <= 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 2 for Nose Treatment (EG_blend_init)!\n");
    return EGADS_GEOMERR;
  }
  if ((begRC == 1) && (endRC == 1) && (nsec <= 3)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 3 for 2Nose Treatment (EG_blend_init)!\n");
    return EGADS_GEOMERR;
  }
  
  /* get our storage */
  eCache = (egCache *) EG_alloc(sizeof(egCache));
  if (eCache == NULL) {
    printf(" EGADS Error: Cannot allocate cache (EG_blend_init)!\n");
    return EGADS_MALLOC;
  }
  eCache->magic  = MAGIC;
  eCache->secs   = NULL;
  eCache->rc1    = NULL;
  eCache->rcN    = NULL;
  eCache->vKnots = NULL;
  eCache->vData  = NULL;
  eCache->sequ   = NULL;
  eCache->strips = (egStrip *) EG_alloc((nstripe+2)*sizeof(egStrip));
  if (eCache->strips == NULL) {
    EG_sens_free(eCache);
    printf(" EGADS Error: Cannot allocate %d Strips (EG_blend_init)!\n",
           nstripe);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nstripe+2; i++) {
    eCache->strips[i].imax  = 0;
    eCache->strips[i].objs  = NULL;
    eCache->strips[i].tbeg  = NULL;
    eCache->strips[i].tend  = NULL;
    eCache->strips[i].xyzs  = NULL;
    eCache->strips[i].ts    = NULL;
    eCache->strips[i].dataS = NULL;
  }
  eCache->outLevel = outLevel;
  eCache->nStrips  = nstripe;
  eCache->nSec     = nsec;
  eCache->planar   = planar;
  eCache->secs     = (ego *) EG_alloc(nsec*sizeof(ego));
  if (eCache->secs == NULL) {
    EG_sens_free(eCache);
    printf(" EGADS Error: Cannot allocate %d Sections (EG_blend_init)!\n",
           nsec);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nsec; i++) eCache->secs[i] = secs[i];
  if (rc1 != NULL) {
    jj = 8;
    if (begRC != 1) {
      jj = 4;
      if (rc1[0] == 0.0) {
        jj = 0;
        if (nstripe <= 3) {
          jj    = 2;
          begRC = 2;
        }
      }
    }
    if (jj > 0) {
      eCache->rc1 = (SurrealS<1> *) EG_alloc(jj*sizeof(SurrealS<1>));
      if (eCache->rc1 == NULL) {
        EG_sens_free(eCache);
        printf(" EGADS Error: Cannot allocate Beg Data (EG_blend_init)!\n");
        return EGADS_MALLOC;
      }
      for (i = 0; i < jj; i++) {
        eCache->rc1[i]         = rc1[i];
        eCache->rc1[i].deriv() = rc1_dot[i];
      }
    }
  }
  if (rcN != NULL) {
    jj = 8;
    if (endRC != 1) {
      jj = 4;
      if (rcN[0] == 0.0) {
        jj = 0;
        if (nstripe <= 3) {
          jj    = 2;
          endRC = 2;
        }
      }
    }
    if (jj > 0) {
      eCache->rcN = (SurrealS<1> *) EG_alloc(jj*sizeof(SurrealS<1>));
      if (eCache->rcN == NULL) {
        EG_sens_free(eCache);
        printf(" EGADS Error: Cannot allocate End Data (EG_blend_init)!\n");
        return EGADS_MALLOC;
      }
      for (i = 0; i < jj; i++) {
        eCache->rcN[i]         = rcN[i];
        eCache->rcN[i].deriv() = rcN_dot[i];
      }
    }
  }
  eCache->begRC  = begRC;
  eCache->endRC  = endRC;
  eCache->vKnots = (SurrealS<1> *) EG_alloc(nsec*sizeof(SurrealS<1>));
  if (eCache->vKnots == NULL) {
    EG_sens_free(eCache);
    printf(" EGADS Error: Cannot allocate vKnots %d (EG_blend_init)!\n", nsec);
    return EGADS_MALLOC;
  }
  eCache->vData = (int *) EG_alloc(nsec*sizeof(int));
  if (eCache->vData == NULL) {
    EG_sens_free(eCache);
    printf(" EGADS Error: Cannot allocate vData %d (EG_blend_init)!\n", nsec);
    return EGADS_MALLOC;
  }
  for (i = 0; i < nstripe; i++) {
    eCache->strips[i].objs = (ego *) EG_alloc(nsec*sizeof(ego));
    if (eCache->strips[i].objs == NULL) {
      EG_sens_free(eCache);
      printf(" EGADS Error: %d Cannot allocate Objs %d (EG_blend_init)!\n",
             i, nsec);
      return EGADS_MALLOC;
    }
  }
  
  /* get the number of sample points per Edge */
  stat = EG_allocSeq(nstripe, &eCache->sequ);
  if (stat !=  EGADS_SUCCESS) {
    EG_sens_free(eCache);
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Edges (EG_blend_init)!\n", nstripe);
    return stat;
  }
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, uvlim, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      for (j = 0; j < nstripe; j++) eCache->strips[j].objs[i] = secs[i];
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      loop = chldrn[0];
    } else {
      loop = secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, uvlim, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (j = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, uvlim, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTopo = %d (EG_blend_init)!\n",
                 i+1, j+1, stat);
        EG_sens_free(eCache);
        return stat;
      }
      if (mtype == DEGENERATE) continue;
      eCache->strips[j].objs[i] = edges[jj];
      
      mtype = -1;
      do {
        if ((mtype == BEZIER) || (mtype == BSPLINE)) EG_free(iinfo);
        geom = ref;
        stat = EG_getGeometry(geom, &oclass, &mtype, &ref, &iinfo, &rinfo);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d getGeom = %d (EG_blend_init)!\n",
                   i+1, j+1, stat);
          EG_sens_free(eCache);
          return stat;
        }
        if (mtype != BSPLINE) EG_free(rinfo);
      } while ((mtype == TRIMMED) || (mtype == OFFSET));

      stat = EGADS_NOTFOUND;
      if (mtype == LINE) {
        stat = EG_setSeq(j, eCache->sequ, 3, NULL, NULL);
      } else if (mtype == CIRCLE) {
/*      k = 16*(uvlim[1]-uvlim[0])/3.1415926;
        if (k < 4) k = 4;
        stat = EG_setSeq(j, eCache->sequ,  k, NULL, NULL);  */
        stat = EG_setSeq(j, eCache->sequ, 12, NULL, NULL);
      } else if (mtype == ELLIPSE) {
        stat = EG_setSeq(j, eCache->sequ, 12, NULL, NULL);
      } else if (mtype == PARABOLA) {
        stat = EG_setSeq(j, eCache->sequ, 12, NULL, NULL);
      } else if (mtype == HYPERBOLA) {
        stat = EG_setSeq(j, eCache->sequ, 12, NULL, NULL);
      } else if (mtype == BEZIER) {
        k    = iinfo[2];
        if (k < 12) k = 12;
        stat = EG_setSeq(j, eCache->sequ, k, NULL, NULL);
        EG_free(iinfo);
      } else if (mtype == BSPLINE) {
        stat = EG_setSeq(j, eCache->sequ, 12, iinfo, rinfo);
        EG_free(iinfo);
        EG_free(rinfo);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d Edge %d setSeq = %d (EG_blend_init)!\n",
                 i+1, j+1, stat);
        EG_sens_free(eCache);
        return stat;
      }
      j++;
    }
  }
  
  *nStrips = nstripe;
  *cache   = eCache;
  return EGADS_SUCCESS;
}


extern "C" int
EG_blend_pos(void *cache, int sIndex, ego **objs, int *imax, double **ts)
{
  int     i, j, k, n, jj, nn, npt, stat, oclass, mtype, nchld, *senses, *sen;
  double  dt, *t, uvlim[4];
  ego     ref, loop, *chldrn, *edges;
  egCache *eCache;
  
  *objs = NULL;
  *imax = 0;
  *ts   = NULL;
  if (cache == NULL)          return EGADS_NULLOBJ;
  eCache = (egCache *) cache;
  if (eCache->magic != MAGIC) return EGADS_NOTOBJ;
  if ((sIndex < 1) || (sIndex > eCache->nStrips)) return EGADS_INDEXERR;
  
  j     = sIndex-1;
  *objs = eCache->strips[j].objs;

  /* make the t map */
  t = (double *) EG_alloc(eCache->sequ[j].ncp*eCache->nSec*sizeof(double));
  if (t == NULL) {
    if (eCache->outLevel > 0)
      printf(" EGADS Error: Malloc on %d x %d (EG_blend_pos)!\n",
             eCache->sequ[j].ncp, eCache->nSec);
    return EGADS_MALLOC;
  }
  for (npt = i = 0; i < eCache->nSec; i++) {
    stat = EG_getTopology(eCache->secs[i], &ref, &oclass, &mtype, uvlim, &n,
                          &chldrn, &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      for (k = 0; k < eCache->sequ[j].ncp; k++, npt++) t[npt] = 0.0;
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      loop = chldrn[0];
    } else {
      loop = eCache->secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, uvlim, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (nn = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, uvlim, &nchld,
                            &chldrn, &sen);
      if (stat != EGADS_SUCCESS) {
        if (eCache->outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTOPO = %d (EG_blend_pos)!\n",
                 i+1, j+1, stat);
        EG_free(t);
        return stat;
      }
      if (mtype == DEGENERATE) continue;
      if (nn == j) {
        dt = uvlim[1] - uvlim[0];
        for (k = 0; k < eCache->sequ[j].ncp; k++, npt++)
          if (senses[jj] == 1) {
            t[npt] = uvlim[0] + eCache->sequ[j].knots[k]*dt;
          } else {
            t[npt] = uvlim[1] - eCache->sequ[j].knots[k]*dt;
          }
        break;
      }
      nn++;
    }
  }
  
  *imax = eCache->sequ[j].ncp;
  *ts   = t;
  eCache->strips[j].ts = t;
  return EGADS_SUCCESS;
}


static int
EG_blend_tip(egCache *eCache, int cap)
{
  int         i, j, k, m, n, stat, nstrip, nKnot, sense, open, fill;
  int         lsen[3], te = -1;
  double      dist, tol, uv[2], norm[3], alen[3], last[3], frac = 0.1;
  SurrealS<1> ratio, len, t, mlen;
  SurrealS<1> x0[3], x1[3], lnrm[3], derivS[18], *uKnots, *data, *east, *west;

  nstrip = eCache->nStrips;
  uv[1]  = 0.0;
  ratio  = 1.0;
  fill   = nstrip;
  if (eCache->rc1 != NULL) ratio = eCache->rc1[1];
  if (cap == 1) {
    uv[1] = 1.0;
    fill++;
    if (eCache->rcN != NULL) ratio = eCache->rcN[1];
  }
  
  /* get cap normal */
  norm[0] = norm[1] = norm[2] = 0.0;
  alen[0] = alen[1] = alen[2] = 0.0;
  last[0] = last[1] = last[2] = 0.0;
  lsen[0] = lsen[1] = lsen[2] = 0;
  for (i = 0; i < nstrip; i++) {
    for (j = 0; j < 4; j++) {
      uv[0] = j/3.0;
      stat  = EG_spline2dDeriv(eCache->strips[i].header,
                               eCache->strips[i].dataS, 2, uv, derivS);
      if (stat != EGADS_SUCCESS) continue;
      norm[0] += derivS[ 9].value();
      norm[1] += derivS[10].value();
      norm[2] += derivS[11].value();
      if (j != 0)
        alen[i] += sqrt((derivS[0].value()-last[0])*(derivS[0].value()-last[0])+
                        (derivS[1].value()-last[1])*(derivS[1].value()-last[1])+
                        (derivS[2].value()-last[2])*(derivS[2].value()-last[2]));
      last[0]  = derivS[ 0].value();
      last[1]  = derivS[ 1].value();
      last[2]  = derivS[ 2].value();
    }
    if (cap == 0) {
      lsen[i] = -1;
      if (eCache->sequ[i].knots[0] > 0.5) lsen[i] =  1;
    } else {
      lsen[i] =  1;
      if (eCache->sequ[i].knots[0] > 0.5) lsen[i] = -1;
    }
  }
  dist = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
  if (cap == 0) dist = -dist;
  norm[0] /= dist;
  norm[1] /= dist;
  norm[2] /= dist;
#ifdef DEBUG
  printf("  Cap %d: norm = %lf %lf %lf\n", cap, norm[0], norm[1], norm[2]);
  printf("          lsenses = %d %d %d\n", lsen[0], lsen[1], lsen[2]);
#endif
  
  /* find the TE strip when 3 exist */
  tol = 1.e100;
  if (nstrip == 3) {
    te   = 0;
    dist = alen[0];
    if (alen[1] < dist) {
      te   = 1;
      dist = alen[1];
    }
    if (alen[2] < dist) te = 2;
    tol = 0.01*alen[te];
  }
#ifdef DEBUG
  printf("     te = %d, alen = %lf %lf %lf\n", te, alen[0], alen[1], alen[2]);
#endif
  
  /* compute the uKnots */
  uKnots = NULL;
  for (nKnot = i = 0; i < nstrip; i++) {
    if (i == te) continue;
    if (nKnot == 0) {
      nKnot = eCache->sequ[i].ncp;
      uKnots = (SurrealS<1> *) EG_alloc(nKnot*sizeof(SurrealS<1>));
      if (uKnots == NULL) {
        printf(" EGADS Error: Allocting %d Knots for Strip %d (EG_blend_tip)!\n",
               nKnot, i+1);
        return EGADS_MALLOC;
      }
      for (j = 0; j < nKnot; j++) uKnots[j] = eCache->sequ[i].knots[j];
      sense = lsen[i];
    } else {
      if (nKnot != eCache->sequ[i].ncp) {
        EG_free(uKnots);
        printf(" EGADS Error: nKnot = %d %d for Strip %d (EG_blend_tip)!\n",
               nKnot, eCache->sequ[i].ncp, i+1);
        return EGADS_GEOMERR;
      }
/*@-nullderef@*/
      if (lsen[i] == sense) {
        for (j = 0; j < nKnot; j++)
          uKnots[nKnot-j-1] += 1.0 - eCache->sequ[i].knots[j];
      } else {
        for (j = 0; j < nKnot; j++) uKnots[j] += eCache->sequ[i].knots[j];
      }
      for (j = 0; j < nKnot; j++) uKnots[j] /= 2.0;
/*@+nullderef@*/
    }
  }

  /* make the grid of points - nKnots by 3/5 */
#ifdef SHARPER
  data = (SurrealS<1> *) EG_alloc(5*3*nKnot*sizeof(SurrealS<1>));
#else
  data = (SurrealS<1> *) EG_alloc(7*3*nKnot*sizeof(SurrealS<1>));
#endif
  if (data == NULL) {
    EG_free(uKnots);
    printf(" EGADS Error: %d x 3 Allocate for Face %d (EG_blend_tip)!\n",
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

  /* fill in the bounds */
  for (n = i = 0; i < nstrip; i++) {
    if (i == te) continue;
    for (j = 0; j < nKnot; j++) {
      uv[0] = uKnots[j].value();
      if ((n != 0) && (lsen[i] == sense)) uv[0] = 1.0 - uKnots[j].value();
      stat = EG_spline2dDeriv(eCache->strips[i].header,
                              eCache->strips[i].dataS, 2, uv, derivS);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: %d/%d Eval = %d for Strip %d (EG_blend_tip)!\n",
               j, nKnot, stat, i+1);
        EG_free(data);
        EG_free(uKnots);
        return stat;
      }
      data[(n*nKnot+j)*3  ] = derivS[0];
      data[(n*nKnot+j)*3+1] = derivS[1];
      data[(n*nKnot+j)*3+2] = derivS[2];
      x1[0] = derivS[ 9];
      x1[1] = derivS[10];
      x1[2] = derivS[11];
      len   = sqrt(DOT(x1,x1));
      if (len != 0.0) {
        x1[0] /= len;
        x1[1] /= len;
        x1[2] /= len;
      }
      len  = 1.0;
      if (n == 0) {
        if (DOT(x1, norm) < 0.0) len = -1.0;
        west[3*j  ] = len*x1[0];
        west[3*j+1] = len*x1[1];
        west[3*j+2] = len*x1[2];
      } else {
        if (DOT(x1, norm) > 0.0) len = -1.0;
        east[3*j  ] = len*x1[0];
        east[3*j+1] = len*x1[1];
        east[3*j+2] = len*x1[2];
      }
    }
#ifdef SHARPER
    n = 2;
#else
    n = 4;
#endif
  }

  /* fill in the middle */
  for (open = j = 0; j < nKnot; j++) {
    x0[0] =      data[3*j  ] - data[(n*nKnot+j)*3  ];
    x0[1] =      data[3*j+1] - data[(n*nKnot+j)*3+1];
    x0[2] =      data[3*j+2] - data[(n*nKnot+j)*3+2];
    x1[0] = 0.5*(data[3*j  ] + data[(n*nKnot+j)*3  ]);
    x1[1] = 0.5*(data[3*j+1] + data[(n*nKnot+j)*3+1]);
    x1[2] = 0.5*(data[3*j+2] + data[(n*nKnot+j)*3+2]);
    len   = sqrt(DOT(x0,x0));
    if (len != 0.0) {
      x0[0] /= len;
      x0[1] /= len;
      x0[2] /= len;
    }
    lnrm[0] = west[3*j  ] - east[3*j  ];
    lnrm[1] = west[3*j+1] - east[3*j+1];
    lnrm[2] = west[3*j+2] - east[3*j+2];
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
    if ((j == 0)       && (len > tol)) open = -1;
    if ((j == nKnot-1) && (len > tol)) open =  1;
    len *= 0.5*ratio;
    if (n == 2) {
      if ((j == 0) || (j == nKnot-1)) len = 0.0;
      data[(  nKnot+j)*3  ] = x1[0] + len*lnrm[0];
      data[(  nKnot+j)*3+1] = x1[1] + len*lnrm[1];
      data[(  nKnot+j)*3+2] = x1[2] + len*lnrm[2];
    } else {
      if ((j == 0)       && (open != -1)) len = 0.0;
      if ((j == nKnot-1) && (open !=  1)) len = 0.0;
      mlen = len;
      if (len != 0.0) {
        mlen = sqrt(west[3*j  ]*west[3*j  ] + west[3*j+1]*west[3*j+1] +
                    west[3*j+2]*west[3*j+2]);
        mlen = frac*len/mlen;
      }
      data[(  nKnot+j)*3  ] = data[3*j  ] + mlen*west[3*j  ];
      data[(  nKnot+j)*3+1] = data[3*j+1] + mlen*west[3*j+1];
      data[(  nKnot+j)*3+2] = data[3*j+2] + mlen*west[3*j+2];
      data[(2*nKnot+j)*3  ] = x1[0] + len*lnrm[0];
      data[(2*nKnot+j)*3+1] = x1[1] + len*lnrm[1];
      data[(2*nKnot+j)*3+2] = x1[2] + len*lnrm[2];
      if (len != 0.0) {
        mlen = sqrt(east[3*j  ]*east[3*j  ] + east[3*j+1]*east[3*j+1] +
                    east[3*j+2]*east[3*j+2]);
        mlen = -frac*len/mlen;
      }
      data[(3*nKnot+j)*3  ] = data[(n*nKnot+j)*3  ] + mlen*east[3*j  ];
      data[(3*nKnot+j)*3+1] = data[(n*nKnot+j)*3+1] + mlen*east[3*j+1];
      data[(3*nKnot+j)*3+2] = data[(n*nKnot+j)*3+2] + mlen*east[3*j+2];
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
      len   = sin(t/j);
      if (n == 2) {
        data[(  nKnot+k)*3  ] = (1.0-len)*x1[0] + len*data[(  nKnot+k)*3  ];
        data[(  nKnot+k)*3+1] = (1.0-len)*x1[1] + len*data[(  nKnot+k)*3+1];
        data[(  nKnot+k)*3+2] = (1.0-len)*x1[2] + len*data[(  nKnot+k)*3+2];
      } else {
        x0[0] = (1.0-frac)*data[3*k  ] + frac*data[(n*nKnot+k)*3  ];
        x0[1] = (1.0-frac)*data[3*k+1] + frac*data[(n*nKnot+k)*3+1];
        x0[2] = (1.0-frac)*data[3*k+2] + frac*data[(n*nKnot+k)*3+2];
        data[(  nKnot+k)*3  ] = (1.0-len)*x0[0] + len*data[(  nKnot+k)*3  ];
        data[(  nKnot+k)*3+1] = (1.0-len)*x0[1] + len*data[(  nKnot+k)*3+1];
        data[(  nKnot+k)*3+2] = (1.0-len)*x0[2] + len*data[(  nKnot+k)*3+2];
        data[(2*nKnot+k)*3  ] = (1.0-len)*x1[0] + len*data[(2*nKnot+k)*3  ];
        data[(2*nKnot+k)*3+1] = (1.0-len)*x1[1] + len*data[(2*nKnot+k)*3+1];
        data[(2*nKnot+k)*3+2] = (1.0-len)*x1[2] + len*data[(2*nKnot+k)*3+2];
        x0[0] = frac*data[3*k  ] + (1.0-frac)*data[(n*nKnot+k)*3  ];
        x0[1] = frac*data[3*k+1] + (1.0-frac)*data[(n*nKnot+k)*3+1];
        x0[2] = frac*data[3*k+2] + (1.0-frac)*data[(n*nKnot+k)*3+2];
        data[(3*nKnot+k)*3  ] = (1.0-len)*x0[0] + len*data[(3*nKnot+k)*3  ];
        data[(3*nKnot+k)*3+1] = (1.0-len)*x0[1] + len*data[(3*nKnot+k)*3+1];
        data[(3*nKnot+k)*3+2] = (1.0-len)*x0[2] + len*data[(3*nKnot+k)*3+2];
      }
    }
  }

  stat = EG_spline2dAprx(0, nKnot, n+1, data, uKnots, (SurrealS<1> *) NULL,
                         NULL, west, east,     (SurrealS<1> *) NULL,
                         (SurrealS<1> *) NULL, (SurrealS<1> *) NULL,
                         (SurrealS<1> *) NULL, 1.e-8,
                          eCache->strips[fill].header,
                         &eCache->strips[fill].dataS);
  if (stat != EGADS_SUCCESS)
    printf(" EGADS Warning: Cap %d EG_spline2dAprx = %d (EG_blend_sens)!\n",
           cap+1, stat);
  
  EG_free(data);
  EG_free(uKnots);
  return EGADS_SUCCESS;
}


extern "C" int
EG_blend_sens(void *cache, int sIndex, const double *xyz, const double *xyz_dot,
              /*@null@*/ const double *tbeg, /*@null@*/ const double *tbeg_dot,
              /*@null@*/ const double *tend, /*@null@*/ const double *tend_dot)
{
  int         i, j, k, n, ii, jj, stat, imax, ix, nsec, n0, left, rite, *vdata;
  egCache     *eCache;
  SurrealS<1> d2xdt2L, d2ydt2L, d2zdt2L, d2sdt2L, nnor[3], x0[3], dist;
  SurrealS<1> d2xdt2R, d2ydt2R, d2zdt2R, d2sdt2R, snor[3], x1[3], *uknots;
  SurrealS<1> dx, dy, *vknot, *xyzs, *south, *souT, *north, *norT, *easT, *wesT;

  if  (cache == NULL)                             return EGADS_NULLOBJ;
  eCache = (egCache *) cache;
  if  (eCache->magic != MAGIC)                    return EGADS_NOTOBJ;
  if ((sIndex < 1) || (sIndex > eCache->nStrips)) return EGADS_INDEXERR;
  if (eCache->strips[sIndex-1].imax != 0)         return EGADS_EXISTS;
  nsec = eCache->nSec;
  
  /* fill in this strip */
  imax = eCache->sequ[sIndex-1].ncp;
  if (tbeg != NULL) {
    eCache->strips[sIndex-1].tbeg = (SurrealS<1> *)
                                    EG_alloc(3*nsec*sizeof(SurrealS<1>));
    if (eCache->strips[sIndex-1].tbeg == NULL) {
      if (eCache->outLevel > 0)
        printf(" EGADS Error: Malloc on tangent beg %d (EG_blend_sens)!\n",
               nsec);
      return EGADS_MALLOC;
    }
  }
  if (tend != NULL) {
    eCache->strips[sIndex-1].tend = (SurrealS<1> *)
                                    EG_alloc(3*nsec*sizeof(SurrealS<1>));
    if (eCache->strips[sIndex-1].tend == NULL) {
      if (eCache->outLevel > 0)
        printf(" EGADS Error: Malloc on tangent end %d (EG_blend_sens)!\n", nsec);
      if (tbeg != NULL) EG_free(eCache->strips[sIndex-1].tbeg);
      eCache->strips[sIndex-1].tbeg = NULL;
      return EGADS_MALLOC;
    }
  }
  eCache->strips[sIndex-1].xyzs = (SurrealS<1> *)
                                  EG_alloc(3*nsec*imax*sizeof(SurrealS<1>));
  if (eCache->strips[sIndex-1].xyzs == NULL) {
    if (eCache->outLevel > 0)
      printf(" EGADS Error: Malloc on xyzs %d x %d (EG_blend_sens)!\n",
             nsec, imax);
    if (tend != NULL) EG_free(eCache->strips[sIndex-1].tend);
    eCache->strips[sIndex-1].tend = NULL;
    if (tbeg != NULL) EG_free(eCache->strips[sIndex-1].tbeg);
    eCache->strips[sIndex-1].tbeg = NULL;
    return EGADS_MALLOC;
  }
  for (i = 0; i < nsec; i++) {
    if (tbeg != NULL) {
      eCache->strips[sIndex-1].tbeg[3*i  ] = tbeg[3*i  ];
      eCache->strips[sIndex-1].tbeg[3*i+1] = tbeg[3*i+1];
      eCache->strips[sIndex-1].tbeg[3*i+2] = tbeg[3*i+2];
      if (tbeg_dot != NULL) {
        eCache->strips[sIndex-1].tbeg[3*i  ].deriv() = tbeg_dot[3*i  ];
        eCache->strips[sIndex-1].tbeg[3*i+1].deriv() = tbeg_dot[3*i+1];
        eCache->strips[sIndex-1].tbeg[3*i+2].deriv() = tbeg_dot[3*i+2];
      } else {
        eCache->strips[sIndex-1].tbeg[3*i  ].deriv() = 0.0;
        eCache->strips[sIndex-1].tbeg[3*i+1].deriv() = 0.0;
        eCache->strips[sIndex-1].tbeg[3*i+2].deriv() = 0.0;
      }
    }
    if (tend != NULL) {
      eCache->strips[sIndex-1].tend[3*i  ] = -tend[3*i  ];
      eCache->strips[sIndex-1].tend[3*i+1] = -tend[3*i+1];
      eCache->strips[sIndex-1].tend[3*i+2] = -tend[3*i+2];
      if (tend_dot != NULL) {
        eCache->strips[sIndex-1].tend[3*i  ].deriv() = -tend_dot[3*i  ];
        eCache->strips[sIndex-1].tend[3*i+1].deriv() = -tend_dot[3*i+1];
        eCache->strips[sIndex-1].tend[3*i+2].deriv() = -tend_dot[3*i+2];
      } else {
        eCache->strips[sIndex-1].tend[3*i  ].deriv() = 0.0;
        eCache->strips[sIndex-1].tend[3*i+1].deriv() = 0.0;
        eCache->strips[sIndex-1].tend[3*i+2].deriv() = 0.0;
      }
    }
    for (j = 0; j < imax; j++) {
      k = i*imax + j;
      eCache->strips[sIndex-1].xyzs[3*k  ]         = xyz[3*k  ];
      eCache->strips[sIndex-1].xyzs[3*k  ].deriv() = xyz_dot[3*k  ];
      eCache->strips[sIndex-1].xyzs[3*k+1]         = xyz[3*k+1];
      eCache->strips[sIndex-1].xyzs[3*k+1].deriv() = xyz_dot[3*k+1];
      eCache->strips[sIndex-1].xyzs[3*k+2]         = xyz[3*k+2];
      eCache->strips[sIndex-1].xyzs[3*k+2].deriv() = xyz_dot[3*k+2];
    }
  }
  eCache->strips[sIndex-1].imax = imax;
  
  /* check if complete */
  for (i = 0; i < eCache->nStrips; i++)
    if (eCache->strips[i].imax == 0) return EGADS_SUCCESS;
  
  /* complete the blend -- fill in vKnots & vData */
  vdata = eCache->vData;
  vknot = eCache->vKnots;
  for (j = 0; j < nsec; j++) vknot[j] = vknot[j].deriv() = 0.0;
  /* inconsistent between this and EG_blend! */
  for (ix = n0 = i = 0; i < eCache->nStrips; i++) {
    xyzs = eCache->strips[i].xyzs;
    imax = eCache->strips[i].imax;
    if (imax > ix) ix = imax;
    for (ii = 0; ii < imax; ii++) {
      dy = 0.0;
      for (j = 1; j < nsec; j++)
        dy += sqrt((xyzs[3*(ii+j*imax)  ]-xyzs[3*(ii+(j-1)*imax)  ])*
                   (xyzs[3*(ii+j*imax)  ]-xyzs[3*(ii+(j-1)*imax)  ]) +
                   (xyzs[3*(ii+j*imax)+1]-xyzs[3*(ii+(j-1)*imax)+1])*
                   (xyzs[3*(ii+j*imax)+1]-xyzs[3*(ii+(j-1)*imax)+1]) +
                   (xyzs[3*(ii+j*imax)+2]-xyzs[3*(ii+(j-1)*imax)+2])*
                   (xyzs[3*(ii+j*imax)+2]-xyzs[3*(ii+(j-1)*imax)+2]));
      if (dy == 0.0) continue;
      n0++;
      dx = 0.0;
      for (j = 1; j < nsec; j++) {
        dx += sqrt((xyzs[3*(ii+j*imax)  ]-xyzs[3*(ii+(j-1)*imax)  ])*
                   (xyzs[3*(ii+j*imax)  ]-xyzs[3*(ii+(j-1)*imax)  ]) +
                   (xyzs[3*(ii+j*imax)+1]-xyzs[3*(ii+(j-1)*imax)+1])*
                   (xyzs[3*(ii+j*imax)+1]-xyzs[3*(ii+(j-1)*imax)+1]) +
                   (xyzs[3*(ii+j*imax)+2]-xyzs[3*(ii+(j-1)*imax)+2])*
                   (xyzs[3*(ii+j*imax)+2]-xyzs[3*(ii+(j-1)*imax)+2]))/dy;
        vknot[j] += dx;
      }
    }
  }
  for (n = j = 1; j < nsec; j++)
    if (vknot[j].value() < vknot[j-1].value()) n++;
  if (n == 1) {
    /* arc-length spaced */
    for (j = 0; j < nsec; j++) vknot[j] /= n0;
  } else {
    /* equally spaced */
    for (jj = j = 0; j < nsec-1; j++)
      if (eCache->secs[j] != eCache->secs[j+1]) jj++;
    for (k = j = 0; j < nsec; j++) {
      dy               = k;
      vknot[j]         = dy/jj;
      vknot[j].deriv() = 0.0;
      if (eCache->secs[j] != eCache->secs[j+1]) k++;
    }
  }

  vdata[0] = vdata[nsec-1] = 1;
  for (j = 1; j < nsec-1; j++)
    if ((vknot[j].value() == vknot[j+1].value()) &&
        (vknot[j].value() == vknot[j+2].value())) {
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
    printf("vdata[%2d]=%2d   %lf\n", j, vdata[j], vknot[j].value());
#endif
  
  /* for all multiplicity 2 knots, determine where the flat spot is */
  for (j = 1; j < nsec-1; j++)
    if ((vdata[j] == -2) && (vdata[j+1] == +2)) {
      if (j < 2) {
        /* flat on left */
        vdata[j+1] = 1;
      } else if (j > nsec-4) {
        /* flat on right */
        vdata[j  ] = 1;
      } else if (vknot[j-1].value() == vknot[j-2].value()) {
        /* flat on left because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on left at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j+1] = 1;
      } else if (vknot[j+3].value() == vknot[j+2].value()) {
        /* flat on right because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on right at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j  ] = 1;
      } else {
        /* find second derivatives on both sides */
        for (left = rite = i = 0; i < eCache->nStrips; i++) {
          xyzs = eCache->strips[i].xyzs;
          imax = eCache->strips[i].imax;
          for (ii = 0; ii < imax; ii++) {
            d2xdt2L = ((xyzs[3*(ii+(j  )*imax)  ]-
                        xyzs[3*(ii+(j-1)*imax)  ])/(vknot[j  ]-vknot[j-1])
                    -  (xyzs[3*(ii+(j-1)*imax)  ]-
                        xyzs[3*(ii+(j-2)*imax)  ])/(vknot[j-1]-vknot[j-2]))
                    / (vknot[j  ]-vknot[j-2]);
            d2ydt2L = ((xyzs[3*(ii+(j  )*imax)+1]-
                        xyzs[3*(ii+(j-1)*imax)+1])/(vknot[j  ]-vknot[j-1])
                    -  (xyzs[3*(ii+(j-1)*imax)+1]-
                        xyzs[3*(ii+(j-2)*imax)+1])/(vknot[j-1]-vknot[j-2]))
                    / (vknot[j  ]-vknot[j-2]);
            d2zdt2L = ((xyzs[3*(ii+(j  )*imax)+2]-
                        xyzs[3*(ii+(j-1)*imax)+2])/(vknot[j  ]-vknot[j-1])
                    -  (xyzs[3*(ii+(j-1)*imax)+2]-
                        xyzs[3*(ii+(j-2)*imax)+2])/(vknot[j-1]-vknot[j-2]))
                    / (vknot[j  ]-vknot[j-2]);
          
            d2xdt2R = ((xyzs[3*(ii+(j+3)*imax)  ]-
                        xyzs[3*(ii+(j+2)*imax)  ])/(vknot[j+3]-vknot[j+2])
                    -  (xyzs[3*(ii+(j+2)*imax)  ]-
                        xyzs[3*(ii+(j+1)*imax)  ])/(vknot[j+2]-vknot[j+1]))
                    / (vknot[j+3]-vknot[j+1]);
            d2ydt2R = ((xyzs[3*(ii+(j+3)*imax)+1]-
                        xyzs[3*(ii+(j+2)*imax)+1])/(vknot[j+3]-vknot[j+2])
                    -  (xyzs[3*(ii+(j+2)*imax)+1]-
                        xyzs[3*(ii+(j+1)*imax)+1])/(vknot[j+2]-vknot[j+1]))
                    / (vknot[j+3]-vknot[j+1]);
            d2zdt2R = ((xyzs[3*(ii+(j+3)*imax)+2]-
                        xyzs[3*(ii+(j+2)*imax)+2])/(vknot[j+3]-vknot[j+2])
                    -  (xyzs[3*(ii+(j+2)*imax)+2]-
                        xyzs[3*(ii+(j+1)*imax)+2])/(vknot[j+2]-vknot[j+1]))
                    / (vknot[j+3]-vknot[j+1]);
          
            d2sdt2L = d2xdt2L*d2xdt2L + d2ydt2L*d2ydt2L + d2zdt2L*d2zdt2L;
            d2sdt2R = d2xdt2R*d2xdt2R + d2ydt2R*d2ydt2R + d2zdt2R*d2zdt2R;
            if (d2sdt2L.value() > d2sdt2R.value()) {
              rite++;
            } else {
              left++;
            }
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
  
  uknots = new SurrealS<1>[ix];
  if (uknots == NULL) {
    printf(" EGADS Fatal Error: Malloc on uKnots %d (EG_blend_sens)!\n", ix);
    EG_sens_free(eCache);
    return EGADS_MALLOC;
  }
  
  /* get the splines */
  for (i = 0; i < eCache->nStrips; i++) {
    xyzs = eCache->strips[i].xyzs;
    imax = eCache->strips[i].imax;
    wesT = eCache->strips[i].tbeg;
    easT = eCache->strips[i].tend;
    if ((imax == 3) || (eCache->planar != 0)) {
      wesT = NULL;
      easT = NULL;
    }
    south = souT = north = norT = NULL;
    if (eCache->begRC == 1) {
      south = eCache->rc1;
      x0[0] = eCache->rc1[1];
      x0[1] = eCache->rc1[2];
      x0[2] = eCache->rc1[3];
      x1[0] = eCache->rc1[5];
      x1[1] = eCache->rc1[6];
      x1[2] = eCache->rc1[7];
      CROSS(snor, x0, x1);
      dist     = 1.0/sqrt(DOT(snor, snor));
      snor[0] *= dist;
      snor[1] *= dist;
      snor[2] *= dist;
      souT     = snor;
    } else if ((eCache->rc1 != NULL) && (eCache->begRC != 2)) {
      snor[0] = eCache->rc1[1]*eCache->rc1[0];
      snor[1] = eCache->rc1[2]*eCache->rc1[0];
      snor[2] = eCache->rc1[3]*eCache->rc1[0];
      souT    = snor;
    }
    if (eCache->endRC == 1) {
      north = eCache->rcN;
      x0[0] = eCache->rcN[1];
      x0[1] = eCache->rcN[2];
      x0[2] = eCache->rcN[3];
      x1[0] = eCache->rcN[5];
      x1[1] = eCache->rcN[6];
      x1[2] = eCache->rcN[7];
      CROSS(nnor, x0, x1);
      dist     = 1.0/sqrt(DOT(nnor, nnor));
      nnor[0] *= dist;
      nnor[1] *= dist;
      nnor[2] *= dist;
      norT     = nnor;
    } else if ((eCache->rcN != NULL) && (eCache->endRC != 2)) {
      nnor[0] = eCache->rcN[1]*eCache->rcN[0];
      nnor[1] = eCache->rcN[2]*eCache->rcN[0];
      nnor[2] = eCache->rcN[3]*eCache->rcN[0];
      norT    = nnor;
    }
    for (j = 0; j < imax; j++) {
      uknots[j]         = eCache->sequ[i].knots[j];
      uknots[j].deriv() = 0.0;
    }
    stat = EG_spline2dAprx(1, imax, nsec, xyzs, uknots, eCache->vKnots,
                           eCache->vData, wesT, easT, south, souT, north, norT,
                           1.e-8,  eCache->strips[i].header,
                                  &eCache->strips[i].dataS);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: Strip %d EG_spline2dAprx = %d (EG_blend_sens)!\n",
             i+1, stat);
  }
  delete [] uknots;
  
  /* get the wing tip splines (if any) */
  if (eCache->begRC == 2) EG_blend_tip(eCache, 0);
  if (eCache->endRC == 2) EG_blend_tip(eCache, 1);
  
  return EGADS_OUTSIDE;
}


extern "C" int
EG_blend_seval(void *cache, int sIndex, const double *uv, double *xyz, double *vel)
{
  int         i, j, stat;
  SurrealS<1> pointS[3];
  egCache     *eCache;
  
  if (cache == NULL) return EGADS_NULLOBJ;
  eCache = (egCache *) cache;
  if (eCache->magic != MAGIC) return EGADS_NOTOBJ;
  if ((sIndex == -1) || (sIndex == -2)) {
    j = eCache->nStrips - sIndex - 1;
  } else {
    if ((sIndex < 1) || (sIndex > eCache->nStrips)) return EGADS_INDEXERR;
    j = sIndex - 1;
  }
  if (eCache->strips[j].dataS == NULL) return EGADS_EMPTY;
  
  stat = EG_spline2dEval(eCache->strips[j].header, eCache->strips[j].dataS, uv,
                         pointS);
  if (stat == EGADS_SUCCESS) {
    for (i = 0; i < 3; i++) {
      xyz[i] = pointS[i].value();
      vel[i] = pointS[i].deriv();
    }
  }
  
  return stat;
}
