/*
 ******************************************************************************
 * EGADS: Electronic Geometry Aircraft Design System                          *
 *                                                                            *
 * This module creates and manages parameterizations                          *
 * Written by John Dannenhoffer @ Syracuse University &                       *
 *            Bob Haimes @ MIT                                                *
 *                                                                            *
 * Copyright 2011-2016, Massachusetts Institute of Technology                 *
 * Licensed under The GNU Lesser General Public License, version 2.1          *
 * See http://www.opensource.org/licenses/lgpl-2.1.php                        *
 ******************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsTris.h"
#include "prm.h"

extern int EG_fillArea(int nconts, const int *cntr, const double *vertices,
                       int *tris, int *nfig8, int pass, fillArea *fa);

#ifdef GRAFIC
#include "grafic.h"
#endif


/*
 * used by the sparse matrix linked list rountines
 */
typedef struct sparseM {
  double         value;
  int            j;
  struct sparseM *next;
} sparseM;


typedef struct {
  int     ni;
  sparseM *is;
} sparseMat;


/*
 * internal routines defined below
 */
static void   surfaceNormal           (int, int, int, prmXYZ[], double*, double[]);
static int    smoothUVcoords          (int, prmTri[], int, prmUV[], double*);
static void   layDownTriangle         (prmXYZ[], prmUV[], int, int, int);
static void   layDownTriangle2        (prmUV[], double, double, int, double, double,
                                       int, double, double, int);
static int    axialProjection         (int, prmTri[], int, prmXYZ[], prmUV[]);
static int    floaterParameterization (int, prmTri[], int, prmXYZ[], prmUV[]);
static int    jonesTransformation     (int, prmTri[], prmUVF[], int, prmUV[]);
static int    planarProjection        (int, prmTri[], int, prmXYZ[], prmUV[]);
static int    simpleProjection        (int, prmTri[], int, prmXYZ[], prmUV[]);
static int    unroll                  (int, prmTri[], int, prmXYZ[], prmUV[]);

static int    fillHoles               (int, prmTri[], int, prmUV[], int*, prmTri*[]);
static int    insideLoop              (double, double, int, prmUV[], int[]);
static double interp1D                (double, int, double[], double[]);
#ifdef GRAFIC
static void   plotTris                (int, prmTri[], int, prmUV[], char[]);
static void   plotTrisImage           (int*, void*, void*, void*, void*, void*,
                                       void*, void*, void*, void*, void*, float*,
                                       char*, int);
#endif

/*
 * external routines used by other prm files
 */
extern int    printSMF                (FILE*, double[], int[], double[], double[]);
extern int    sparseBCG               (double[], int[], double[], double[]);

/*
 * global constants
 */
#define  HUGEQ           1.0e+40
#define  EPS04           1.0e-04         /* tolerance on U and V */
#define  EPS12           1.0e-12         /* tolerance on area */
#define  EPS20           1.0e-20
#define  TWOPI           2.0*PI
#define  MAXLOOPS        500

/*
* useful macros
*/
#define  SQR(A)          ((A) * (A))
#define  MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define  MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define  MAX3(A,B,C)     MAX((A), MAX((B), (C)))
#define  MIN3(A,B,C)     MIN((A), MIN((B), (C)))
#define  MINMAX(A,B,C)   MIN(MAX((A), (B)), (C))
#define  ACOS(A)         acos(MINMAX(-1, (A), +1))
#define  NORM012(A)      (sqrt( SQR(A[0]) + SQR(A[1]) + SQR(A[2])))
#define  DSWAP(A,B)      dswap_temp = A; A = B; B = dswap_temp;
static   double          dswap_temp;

/*
 * macros for progress printing
 */
#ifdef PROGRESS
   #define PPRINT0(FORMAT) \
           printf(FORMAT); printf("\n");
   #define PPRINT1(FORMAT,A) \
           printf(FORMAT, A); printf("\n");
   #define PPRINT2(FORMAT,A,B) \
           printf(FORMAT, A, B); printf("\n");
   #define PPRINT3(FORMAT,A,B,C) \
           printf(FORMAT, A, B, C); printf("\n");
   #define PPRINT4(FORMAT,A,B,C,D) \
           printf(FORMAT, A, B, C, D); printf("\n");
   #define PPRINT5(FORMAT,A,B,C,D,E) \
           printf(FORMAT, A, B, C, D, E); printf("\n");
#else
   #define PPRINT0(FORMAT)
   #define PPRINT1(FORMAT,A)
   #define PPRINT2(FORMAT,A,B)
   #define PPRINT3(FORMAT,A,B,C)
   #define PPRINT4(FORMAT,A,B,C,D)
   #define PPRINT5(FORMAT,A,B,C,D,E)
#endif

/*
 * macros for debug printing
 */
#ifdef   DEBUG
   #define DOPEN \
           {if (dbg_fp == NULL) dbg_fp = fopen("prmUV.dbg", "w");}
   #define DFLUSH \
           {fprintf(dbg_fp, "\n"); fflush(dbg_fp);}
   #define DPRINT0(FORMAT) \
           {DOPEN; fprintf(dbg_fp, FORMAT); DFLUSH;}
   #define DPRINT0x(FORMAT) \
           {DOPEN; fprintf(dbg_fp, FORMAT);}
   #define DPRINT1(FORMAT,A) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A); DFLUSH;}
   #define DPRINT1x(FORMAT,A) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A);}
   #define DPRINT2(FORMAT,A,B) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B); DFLUSH;}
   #define DPRINT2x(FORMAT,A,B) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B);}
   #define DPRINT3(FORMAT,A,B,C) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C); DFLUSH;}
   #define DPRINT4(FORMAT,A,B,C,D) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D); DFLUSH;}
   #define DPRINT5(FORMAT,A,B,C,D,E) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E); DFLUSH;}
   #define DPRINT6(FORMAT,A,B,C,D,E,F) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F); DFLUSH;}
   #define DPRINT7(FORMAT,A,B,C,D,E,F,G) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F, G); DFLUSH;}
   #define DPRINT8(FORMAT,A,B,C,D,E,F,G,H) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F, G, H); DFLUSH;}
   #define DPRINT9(FORMAT,A,B,C,D,E,F,G,H,I) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F, G, H, I); DFLUSH;}
   #define DPRINT10(FORMAT,A,B,C,D,E,F,G,H,I,J) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F, G, H, I, J); DFLUSH;}
   #define DPRINT11(FORMAT,A,B,C,D,E,F,G,H,I,J,K) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F, G, H, I, J, K); DFLUSH;}
   #define DPRINT12(FORMAT,A,B,C,D,E,F,G,H,I,J,K,L)                         \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E, F, G, H, I, J, K, L); DFLUSH;}
   #define DCLOSE \
           {if (dbg_fp != NULL) fclose(dbg_fp); dbg_fp = NULL;}
   #define ROUTINE(NAME) char routine[] = #NAME ;
#else
   #define DPRINT0(FORMAT)
   #define DPRINT0x(FORMAT)
   #define DPRINT1(FORMAT,A)
   #define DPRINT1x(FORMAT,A)
   #define DPRINT2(FORMAT,A,B)
   #define DPRINT2x(FORMAT,A,B)
   #define DPRINT3(FORMAT,A,B,C)
   #define DPRINT4(FORMAT,A,B,C,D)
   #define DPRINT5(FORMAT,A,B,C,D,E)
   #define DPRINT6(FORMAT,A,B,C,D,E,F)
   #define DPRINT7(FORMAT,A,B,C,D,E,F,G)
   #define DPRINT8(FORMAT,A,B,C,D,E,F,G,H)
   #define DPRINT9(FORMAT,A,B,C,D,E,F,G,H,I)
   #define DPRINT10(FORMAT,A,B,C,D,E,F,G,H,I,J)
   #define DPRINT11(FORMAT,A,B,C,D,E,F,G,H,I,J,K)
   #define DPRINT12(FORMAT,A,B,C,D,E,F,G,H,I,J,K,L)
   #define DCLOSE
   #define ROUTINE(NAME)
#endif

#ifdef DEBUG
static FILE     *dbg_fp = NULL;         /* file pointer for debug information */
#endif

/*
 * macros for verbose output
 */
#ifdef VERBOSE
   #define DPRINT0(A)  \
           printf( A); \
           printf( "\n")
   #define DPRINT1(A,B)  \
           printf( A,B); \
           printf( "\n")
   #define DPRINT2(A,B,C)  \
           printf( A,B,C); \
           printf( "\n")
#endif

/*
 * macros for error checking
 */
#define CHECK_STATUS \
        if (status < EGADS_SUCCESS) {\
            DPRINT2("ERROR:: BAD STATUS = %d in routine %s", status, routine);\
            goto cleanup;\
        }
#define MALLOC(PTR,TYPE,SIZE) \
        DPRINT3("mallocing %s in routine %s (size=%d)", #PTR, routine, SIZE);\
        PTR = (TYPE *) EG_alloc((SIZE) * sizeof(TYPE)); \
        if (PTR == NULL) {\
            DPRINT2("ERROR:: MALLOC PROBLEM for %s in routine %s", #PTR, routine);\
            status = EGADS_MALLOC;\
            goto cleanup;\
        }
#define RALLOC(PTR,TYPE,SIZE) \
        DPRINT3("rallocing %s in routine %s (size=%d)", #PTR, routine, SIZE);\
        realloc_temp = EG_reall(PTR, (SIZE) * sizeof(TYPE)); \
        if (PTR == NULL) {\
            DPRINT2("ERROR:: RALLOC PROBLEM for %s in routine %s", #PTR, routine);\
            status = EGADS_MALLOC;\
            goto cleanup;\
        } else {\
            PTR = (TYPE *)realloc_temp;\
        }
#define FREE(PTR) \
        DPRINT2("freeing %s in routine %s", #PTR, routine);\
        EG_free(PTR);\
        PTR = NULL;
//static void *realloc_temp = NULL;            /* used by RALLOC macro */


/*
 ********************************************************************************
 *                                                                              *
 * sparse matrix functions -- note no real removal of matrix entries            *
 *                                                                              *
 ********************************************************************************
 */

static int
initSmat(int ni, sparseMat *mat)
{
  int i;
  
  mat->ni = ni;
  mat->is = NULL;
  if (ni <= 0) return EGADS_INDEXERR;
  
  mat->is = (sparseM *) EG_alloc(ni*sizeof(sparseM));
  if (mat->is == NULL) return EGADS_MALLOC;
  for (i = 0; i < ni; i++) {
    mat->is[i].value = 0.0;
    mat->is[i].j     = i;
    mat->is[i].next  = NULL;
  }
  return EGADS_SUCCESS;
}


static void
freeSmat(sparseMat *mat)
{
  int     i;
  sparseM *next, *save;
  
  if (mat->ni == 0) return;
  for (i = 0; i < mat->ni; i++) {
    next = mat->is[i].next;
    while (next != NULL) {
      save = next->next;
      EG_free(next);
      next = save;
    }
  }
  
  EG_free(mat->is);
  mat->is = NULL;
  mat->ni = 0;
}


static int
countSmat(sparseMat *mat)
{
  int     i, n;
  sparseM *next;

  for (n = i = 0; i < mat->ni; i++) {
    n++;
    next = mat->is[i].next;
    while (next != NULL) {
      n++;
      next = next->next;
    }
  }
  
  return n;
}


static int
diagSmat(sparseMat *mat, int i)
{
  sparseM *next, *save;
  
  if ((i < 0) || (i >= mat->ni)) return EGADS_INDEXERR;

  next = mat->is[i].next;
  while (next != NULL) {
    save = next->next;
    EG_free(next);
    next = save;
  }
  mat->is[i].next = NULL;
  
  return EGADS_SUCCESS;
}


static int
setSmat(sparseMat *mat, int i, int j, double value)
{
  sparseM *next, *save;
  
  if ((i < 0) || (i >= mat->ni)) return EGADS_INDEXERR;
  if ((j < 0) || (j >= mat->ni)) return EGADS_INDEXERR;
  if (i == j) {
    mat->is[i].value = value;
    return EGADS_SUCCESS;
  }
  
  /* look for existing entry */
  save = &mat->is[i];
  next =  mat->is[i].next;
  while (next != NULL) {
    if (next->j == j) {
      next->value = value;
      return EGADS_SUCCESS;
    }
    save = next;
    next = next->next;
  }
  
  /* none found -- make one */
  next = (sparseM *) EG_alloc(sizeof(sparseM));
  if (next == NULL) return EGADS_MALLOC;
  next->value = value;
  next->j     = j;
  next->next  = NULL;
  save->next  = next;
  
  return EGADS_SUCCESS;
}


static int
getSmat(sparseMat *mat, int i, int j, double *value)
{
  sparseM *next;
  
  *value = 0.0;
  if ((i < 0) || (i >= mat->ni)) return EGADS_INDEXERR;
  if ((j < 0) || (j >= mat->ni)) return EGADS_INDEXERR;
  if (i == j) {
    *value = mat->is[i].value;
    return EGADS_SUCCESS;
  }
  
  /* look for existing entry */
  next = mat->is[i].next;
  while (next != NULL) {
    if (next->j == j) {
      *value = next->value;
      return EGADS_SUCCESS;
    }
    next = next->next;
  }
  
  /* none found */
  return EGADS_NOTFOUND;
}


static int
sumSmat(sparseMat *mat, int i, double *sum)
{
  sparseM *next;
  
  *sum = 0.0;
  if ((i < 0) || (i >= mat->ni)) return EGADS_INDEXERR;

  next = mat->is[i].next;
  while (next != NULL) {
    *sum = next->value + *sum;
    next = next->next;
  }
  
  return EGADS_SUCCESS;
}


static int
divSmat(sparseMat *mat, int i)
{
  double  sum = 0.0;
  sparseM *next;

  if ((i < 0) || (i >= mat->ni)) return EGADS_INDEXERR;
  
  sum  = mat->is[i].value;
  next = mat->is[i].next;
  while (next != NULL) {
    sum += next->value;
    next = next->next;
  }
  
  sum = -sum;
  mat->is[i].value /= sum;
  next = mat->is[i].next;
  while (next != NULL) {
    next->value /= sum;
    next = next->next;
  }
  
  return EGADS_SUCCESS;
}


static int
fillSmat(sparseMat *mat, int i, int *kk, double *asmf, int *ismf)
{
  int     k;
  sparseM *next;
  
  if ((i < 0) || (i >= mat->ni)) return EGADS_INDEXERR;
  
  k = *kk;
  
  next = mat->is[i].next;
  while (next != NULL) {
    k++;
    asmf[k] = next->value;
    ismf[k] = next->j;
    next    = next->next;
  }
  
  *kk = k;
  return EGADS_SUCCESS;
}

/*
 ********************************************************************************
 *                                                                              *
 * surfaceNormal -- compute surface normal for a Triangle                       *
 *                                                                              *
 ********************************************************************************
 */
static void
surfaceNormal(int      iv0,                  /* (in)   first  Vertex (bias-0) */
              int      iv1,                  /* (in)   second Vertex (bias-0) */
              int      iv2,                  /* (in)   third  Vertex (bias-0) */
              prmXYZ   xyz[],                /* (in)   array of Coordinates */
              double   *area,                /* (out)  area of Triangle */
              double   norm[])               /* (out)  components of the unit normal */
{

//    ROUTINE(surfaceNormal);

    /* ----------------------------------------------------------------------- */

    if ((iv0 < 0) || (iv1 < 0) || (iv2 < 0)) {
        *area   = 0.0;
        norm[0] = 0.0;
        norm[1] = 0.0;
        norm[2] = 0.0;

    } else {
        norm[0] = (xyz[iv1].y - xyz[iv0].y) * (xyz[iv2].z - xyz[iv0].z)
                - (xyz[iv2].y - xyz[iv0].y) * (xyz[iv1].z - xyz[iv0].z);
        norm[1] = (xyz[iv1].z - xyz[iv0].z) * (xyz[iv2].x - xyz[iv0].x)
                - (xyz[iv2].z - xyz[iv0].z) * (xyz[iv1].x - xyz[iv0].x);
        norm[2] = (xyz[iv1].x - xyz[iv0].x) * (xyz[iv2].y - xyz[iv0].y)
                - (xyz[iv2].x - xyz[iv0].x) * (xyz[iv1].y - xyz[iv0].y);

        *area = NORM012(norm);

        if ((*area) > 0) {
            norm[0] /= (*area);
            norm[1] /= (*area);
            norm[2] /= (*area);
        }
    }

// cleanup:
}



/*
 ********************************************************************************
 *                                                                              *
 * smoothUVcoords -- smooth UV coordinates                                      *
 *                                                                              *
 ********************************************************************************
 */
static int
smoothUVcoords(int      ntri,                /* (in)   number of Triangles */
               prmTri   tri[],               /* (in)   array  of Triangles */
               int      nvrt,                /* (in)   number of Vertices */
               prmUV    uv[],                /* (both) array of Parameters */
               double   *amin)               /* (out)  minimum area */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    double      *du = NULL;
    double      *dv = NULL;
    double      *nn = NULL;

    int         ivrt, iv0, iv1, iv2;
    int         itri;
    double      usum, vsum, area2;

    double      omega = 0.50;

    ROUTINE(smoothUVcoords);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * get work arrays
     */
    MALLOC(du, double, nvrt);
    MALLOC(dv, double, nvrt);
    MALLOC(nn, double, nvrt);

    /*
     * initialize the sums (du, dv, and nn)
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        du[ivrt] = 0;
        dv[ivrt] = 0;
        nn[ivrt] = 0;
    }

    /*
     * accumulate sums (du, dv, and nn) at each of the Vertices
     */
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        usum = uv[iv0].u + uv[iv1].u + uv[iv2].u;
        vsum = uv[iv0].v + uv[iv1].v + uv[iv2].v;

        du[iv0] += usum - 3 * uv[iv0].u;
        du[iv1] += usum - 3 * uv[iv1].u;
        du[iv2] += usum - 3 * uv[iv2].u;
        dv[iv0] += vsum - 3 * uv[iv0].v;
        dv[iv1] += vsum - 3 * uv[iv1].v;
        dv[iv2] += vsum - 3 * uv[iv2].v;
        nn[iv0] += 2;
        nn[iv1] += 2;
        nn[iv2] += 2;
    }

    /*
     * zero out du and dv on the boundary
     */
    for (itri = 0; itri < ntri; itri++) {
        if (tri[itri].neigh[1] <= 0 || tri[itri].neigh[2] <= 0) {
            du[tri[itri].indices[0]-1] = 0.0;
            dv[tri[itri].indices[0]-1] = 0.0;
        }

        if (tri[itri].neigh[2] <= 0 || tri[itri].neigh[0] <= 0) {
            du[tri[itri].indices[1]-1] = 0.0;
            dv[tri[itri].indices[1]-1] = 0.0;
        }

        if (tri[itri].neigh[0] <= 0 || tri[itri].neigh[1] <= 0) {
            du[tri[itri].indices[2]-1] = 0.0;
            dv[tri[itri].indices[2]-1] = 0.0;
        }
    }

    /*
     * update the Vertex locations with under-relaxation
     */

    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uv[ivrt].u += omega * du[ivrt] / nn[ivrt];
        uv[ivrt].v += omega * dv[ivrt] / nn[ivrt];
    }

    /*
     * keep track of the smallest area
     */
    *amin = +HUGEQ;
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        area2 = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
              - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
        *amin = MIN(*amin, area2);
    }

 cleanup:
    FREE(du);
    FREE(dv);
    FREE(nn);

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * layDownTriangle -- lay down a Triangle to unroll                             *
 *                                                                              *
 ********************************************************************************
 */
static void
layDownTriangle(prmXYZ    xyz[],             /* (both) array of Coordinates */
                prmUV     uv[],              /* (both) array of Parameters */
                int       iv0,               /* (in)   base Vertex 1 */
                int       iv1,               /* (in)   base Vertex 2 */
                int       iv2)               /* (in)   Vertex to be defined */
{
    double len01, len12, len20, phi, tht, area;

//    ROUTINE(layDownTriangle);

    /* ----------------------------------------------------------------------- */

    /*
     * compute lengths
     */
    len01 = sqrt(  SQR(xyz[iv0].x - xyz[iv1].x)
                 + SQR(xyz[iv0].y - xyz[iv1].y)
                 + SQR(xyz[iv0].z - xyz[iv1].z));
    len12 = sqrt(  SQR(xyz[iv1].x - xyz[iv2].x)
                 + SQR(xyz[iv1].y - xyz[iv2].y)
                 + SQR(xyz[iv1].z - xyz[iv2].z));
    len20 = sqrt(  SQR(xyz[iv2].x - xyz[iv0].x)
                 + SQR(xyz[iv2].y - xyz[iv0].y)
                 + SQR(xyz[iv2].z - xyz[iv0].z));

    /*
     * place Vertex iv2 (using law of cosines)
     */
    phi = atan2(uv[iv1].v - uv[iv0].v, uv[iv1].u - uv[iv0].u);
    tht = ACOS( (SQR(len01) - SQR(len12) + SQR(len20))
              / (    len01  *   2.0      *     len20 ));

    uv[iv2].u = uv[iv0].u + len20 * cos(phi - tht);
    uv[iv2].v = uv[iv0].v + len20 * sin(phi - tht);

    /*
     * make Triangle on other side of 0-1 if a negative area
     */
    area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
         - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);

    if (area < 0.0) {
        uv[iv2].u = uv[iv0].u + len20 * cos(phi + tht);
        uv[iv2].v = uv[iv0].v + len20 * sin(phi + tht);
    }
}



/*
 ********************************************************************************
 *                                                                              *
 * layDownTriangle2 -- lay down a Triangle for jones' transformation            *
 *                                                                              *
 ********************************************************************************
 */
static void
layDownTriangle2(prmUV     uv[],
                 double    u0,
                 double    v0,
                 int       iv0,
                 double    u1,
                 double    v1,
                 int       iv1,
                 double    u2,
                 double    v2,
                 int       iv2)
{
    double len01, len12, len20, phi, tht;

//    ROUTINE(layDownTriangle2);

    /* ----------------------------------------------------------------------- */

    /*
     * compute lengths
     */
    len01 = sqrt(SQR(u0 - u1) + SQR(v0 - v1));
    len12 = sqrt(SQR(u1 - u2) + SQR(v1 - v2));
    len20 = sqrt(SQR(u2 - u0) + SQR(v2 - v0));

    /*
     * place Vertex iv2 (using law of cosines)
     */
    phi = atan2(uv[iv1].v - uv[iv0].v,
                uv[iv1].u - uv[iv0].u);
    tht = ACOS((SQR(len20) - SQR(len12) + SQR(len01)) / (2.0 * len20 * len01));

    uv[iv2].u = uv[iv0].u + len20 * cos(phi + tht);
    uv[iv2].v = uv[iv0].v + len20 * sin(phi + tht);
}



/*
 ********************************************************************************
 *                                                                              *
 * axialProjection -- use axial projection (if possible)                        *
 *                                                                              *
 ********************************************************************************
 */
static int
axialProjection(int          ntri,           /* (in)   number of Triangles */
                prmTri       tri[],          /* (in)   array  of Triangles */
                int          nvrt,           /* (in)   number of Vertices */
                prmXYZ       xyz[],          /* (in)   array  of Coordinates */
                prmUV        uv[])           /* (out)  array  of Parameters */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_AXIAL */

    double      area;
    int         change;
    double      dx, dy, dz;
    int         itri;
    int         ivrt, iv0, iv1, iv2;
    int         nneg;
    int         pass;
    double      theta0, theta1, theta2;
    double      xcent, ycent, zcent;
    double      xmin,  xmax,  ymin,  ymax,  zmin,  zmax;
    double      xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
    double      xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;

    ROUTINE(axialProjection);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * find the minimum and maximum extents of the Vertices in the tessellation
     */
    xmin = xyz[0].x;
    xmax = xyz[0].x;
    ymin = xyz[0].y;
    ymax = xyz[0].y;
    zmin = xyz[0].z;
    zmax = xyz[0].z;

    for (ivrt = 1; ivrt < nvrt; ivrt++) {
        xmin = MIN(xmin, xyz[ivrt].x);
        xmax = MAX(xmax, xyz[ivrt].x);
        ymin = MIN(ymin, xyz[ivrt].y);
        ymax = MAX(ymax, xyz[ivrt].y);
        zmin = MIN(zmin, xyz[ivrt].z);
        zmax = MAX(zmax, xyz[ivrt].z);
    }
    DPRINT3("xmin=%10.5f   ymin=%10.5f   zmin=%10.5f", xmin, ymin, zmin);
    DPRINT3("xmax=%10.5f   ymax=%10.5f   zmax=%10.5f", xmax, ymax, zmax);

    /*
     * the x-extent is larger than the other two extents, so try an
     *    axis that is aligned with the x-axis
     */
    if ((xmax-xmin) > (ymax-ymin) && (xmax-xmin) > (zmax-zmin)) {
        DPRINT0("potential axis in x-direction");

        /*
         * find the center of the Vertices within the last 1 percent of the x extent
         */
        dx   = 0.01 * (xmax - xmin);

        ymin1 = +HUGEQ;
        ymax1 = -HUGEQ;
        ymin2 = +HUGEQ;
        ymax2 = -HUGEQ;
        zmin1 = +HUGEQ;
        zmax1 = -HUGEQ;
        zmin2 = +HUGEQ;
        zmax2 = -HUGEQ;

        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (xyz[ivrt].x < xmin+dx) {
                ymin1 = MIN(ymin1, xyz[ivrt].y);
                ymax1 = MAX(ymax1, xyz[ivrt].y);
                zmin1 = MIN(zmin1, xyz[ivrt].z);
                zmax1 = MAX(zmax1, xyz[ivrt].z);
            } else if (xyz[ivrt].x > xmax-dx) {
                ymin2 = MIN(ymin2, xyz[ivrt].y);
                ymax2 = MAX(ymax2, xyz[ivrt].y);
                zmin2 = MIN(zmin2, xyz[ivrt].z);
                zmax2 = MAX(zmax2, xyz[ivrt].z);
            }
        }

        ymin = MAX(ymin1, ymin2);
        ymax = MIN(ymax1, ymax2);
        zmin = MAX(zmin1, zmin2);
        zmax = MIN(zmax1, zmax2);

        /*
         * try placing the axis at the center and at each of the extents
         */
        for (pass = 0; pass < 5; pass++) {
            if (pass == 0) {
                xcent = (xmin + xmax) / 2;
                ycent = (ymin + ymax) / 2;
                zcent = (zmin + zmax) / 2;
            } else if (pass == 1) {
                xcent = (xmin + xmax) / 2;
                ycent =  ymin            ;
                zcent = (zmin + zmax) / 2;
            } else if (pass == 2) {
                xcent = (xmin + xmax) / 2;
                ycent =         ymax     ;
                zcent = (zmin + zmax) / 2;
            } else if (pass == 3) {
                xcent = (xmin + xmax) / 2;
                ycent = (ymin + ymax) / 2;
                zcent =  zmin            ;
            } else {
                xcent = (xmin + xmax) / 2;
                ycent = (ymin + ymax) / 2;
                zcent =         zmax     ;
            }
            DPRINT3("xcent=%10.5f  ycent=%10.5f  zcent=%10.5f", xcent, ycent, zcent);

            /*
             * assign u=x and v=theta for every Vertex
             */
            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                uv[ivrt].u = xmin - 1;
            }

            uv[0].u = xyz[0].x;
            uv[0].v = atan2(xyz[0].y - ycent, xyz[0].z - zcent);

            change = 1;
            while (change > 0) {
                change = 0;

                for (itri = 0; itri < ntri; itri++) {
                    iv0 = tri[itri].indices[0] - 1;
                    iv1 = tri[itri].indices[1] - 1;
                    iv2 = tri[itri].indices[2] - 1;


                    if (uv[iv0].u >= xmin) {
                        theta0 = uv[iv0].v;

                        if (uv[iv1].u < xmin) {
                            theta1 = atan2(xyz[iv1].y - ycent, xyz[iv1].z - zcent);

                            if (theta1 < theta0 - PI) theta1 += TWOPI;
                            if (theta1 > theta0 + PI) theta1 -= TWOPI;

                            uv[iv1].u = xyz[iv1].x;
                            uv[iv1].v = theta1;
                            change++;
                        }

                        if (uv[iv2].u < xmin) {
                            theta2 = atan2(xyz[iv2].y - ycent, xyz[iv2].z - zcent);

                            if (theta2 < theta0 - PI) theta2 += TWOPI;
                            if (theta2 > theta0 + PI) theta2 -= TWOPI;

                            uv[iv2].u = xyz[iv2].x;
                            uv[iv2].v = theta2;
                            change++;
                        }
                    }

                    if (uv[iv1].u >= xmin) {
                        theta1 = uv[iv1].v;

                        if (uv[iv2].u < xmin) {
                            theta2 = atan2(xyz[iv2].y - ycent, xyz[iv2].z - zcent);

                            if (theta2 < theta1 - PI) theta2 += TWOPI;
                            if (theta2 > theta1 + PI) theta2 -= TWOPI;

                            uv[iv2].u = xyz[iv2].x;
                            uv[iv2].v = theta2;
                            change++;
                        }

                        if (uv[iv0].u < xmin) {
                            theta0 = atan2(xyz[iv0].y - ycent, xyz[iv0].z - zcent);

                            if (theta0 < theta1 - PI) theta0 += TWOPI;
                            if (theta0 > theta1 + PI) theta0 -= TWOPI;

                            uv[iv0].u = xyz[iv0].x;
                            uv[iv0].v = theta0;
                            change++;
                        }
                    }

                    if (uv[iv2].u >= xmin) {
                        theta2 = uv[iv2].v;

                        if (uv[iv0].u < xmin) {
                            theta0 = atan2(xyz[iv0].y - ycent, xyz[iv0].z - zcent);

                            if (theta0 < theta2 - PI) theta0 += TWOPI;
                            if (theta0 > theta2 + PI) theta0 -= TWOPI;

                            uv[iv0].u = xyz[iv0].x;
                            uv[iv0].v = theta0;
                            change++;
                        }

                        if (uv[iv1].u < xmin) {
                            theta1 = atan2(xyz[iv1].y - ycent, xyz[iv1].z - zcent);

                            if (theta1 < theta2 - PI) theta1 += TWOPI;
                            if (theta1 > theta2 + PI) theta1 -= TWOPI;

                            uv[iv1].u = xyz[iv1].x;
                            uv[iv1].v = theta1;
                            change++;
                        }
                    }
                }
                DPRINT1("   change=%d", change);
            }

            /*
             * make sure that all Triangles have outward pointing normals
             */
            nneg = 0;
            for (itri = 0; itri < ntri; itri++) {
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                     - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
                if (area <= 0) {
                    nneg++;
                }
            }

            /*
             * if we got here, then the x-axial projection was successful
             */
            if (nneg == 0) {
                DPRINT0("   x-axial projection was successful");

                status = PRM_OK_AXIAL;
                goto cleanup;
            } else {
                DPRINT2("   there were %d non-positive areas (ntri=%d)", nneg, ntri);
            }
        }
    }

    /*
     * the y-extent is larger than the other two extents, so try an
     *    axis that is aligned with the y-axis
     */
    if ((ymax-ymin) > (zmax-zmin) && (ymax-ymin) > (xmax-xmin)) {
        DPRINT0("potential axis in y-direction");

        /*
         * find the center of the Vertices within the last 1 percent of the x extent
         */
        dy   = 0.01 * (ymax - ymin);

        zmin1 = +HUGEQ;
        zmax1 = -HUGEQ;
        zmin2 = +HUGEQ;
        zmax2 = -HUGEQ;
        xmin1 = +HUGEQ;
        xmax1 = -HUGEQ;
        xmin2 = +HUGEQ;
        xmax2 = -HUGEQ;

        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (xyz[ivrt].y < ymin+dy) {
                zmin1 = MIN(zmin1, xyz[ivrt].z);
                zmax1 = MAX(zmax1, xyz[ivrt].z);
                xmin1 = MIN(xmin1, xyz[ivrt].x);
                xmax1 = MAX(xmax1, xyz[ivrt].x);
            } else if (xyz[ivrt].y > ymax-dy) {
                zmin2 = MIN(zmin2, xyz[ivrt].z);
                zmax2 = MAX(zmax2, xyz[ivrt].z);
                xmin2 = MIN(xmin2, xyz[ivrt].x);
                xmax2 = MAX(xmax2, xyz[ivrt].x);
            }
        }

        zmin = MAX(zmin1, zmin2);
        zmax = MIN(zmax1, zmax2);
        xmin = MAX(xmin1, xmin2);
        xmax = MIN(xmax1, xmax2);

        /*
         * try placing the axis at the center and at each of the extents
         */
        for (pass = 0; pass < 5; pass++) {
            if (pass == 0) {
                ycent = (ymin + ymax) / 2;
                zcent = (zmin + zmax) / 2;
                xcent = (xmin + xmax) / 2;
            } else if (pass == 1) {
                ycent = (ymin + ymax) / 2;
                zcent =  zmin            ;
                xcent = (xmin + xmax) / 2;
            } else if (pass == 2) {
                ycent = (ymin + ymax) / 2;
                zcent =         zmax     ;
                xcent = (xmin + xmax) / 2;
            } else if (pass == 3) {
                ycent = (ymin + ymax) / 2;
                zcent = (zmin + zmax) / 2;
                xcent =  xmin            ;
            } else {
                ycent = (ymin + ymax) / 2;
                zcent = (zmin + zmax) / 2;
                xcent =         xmax     ;
            }
            DPRINT3("xcent=%10.5f  ycent=%10.5f  zcent=%10.5f", xcent, ycent, zcent);

            /*
             * assign u=y and v=theta for every Veretx
             */
            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                uv[ivrt].u = ymin - 1;
            }

            uv[0].u = xyz[0].y;
            uv[0].v = atan2(xyz[0].z - zcent, xyz[0].x - xcent);

            change = 1;
            while (change > 0) {
                change = 0;

                for (itri = 0; itri < ntri; itri++) {
                    iv0 = tri[itri].indices[0] - 1;
                    iv1 = tri[itri].indices[1] - 1;
                    iv2 = tri[itri].indices[2] - 1;


                    if (uv[iv0].u >= ymin) {
                        theta0 = uv[iv0].v;

                        if (uv[iv1].u < ymin) {
                            theta1 = atan2(xyz[iv1].z - zcent, xyz[iv1].x - xcent);

                            if (theta1 < theta0 - PI) theta1 += TWOPI;
                            if (theta1 > theta0 + PI) theta1 -= TWOPI;

                            uv[iv1].u = xyz[iv1].y;
                            uv[iv1].v = theta1;
                            change++;
                        }

                        if (uv[iv2].u < ymin) {
                            theta2 = atan2(xyz[iv2].z - zcent, xyz[iv2].x - xcent);

                            if (theta2 < theta0 - PI) theta2 += TWOPI;
                            if (theta2 > theta0 + PI) theta2 -= TWOPI;

                            uv[iv2].u = xyz[iv2].y;
                            uv[iv2].v = theta2;
                            change++;
                        }
                    }

                    if (uv[iv1].u >= ymin) {
                        theta1 = uv[iv1].v;

                        if (uv[iv2].u < ymin) {
                            theta2 = atan2(xyz[iv2].z - zcent, xyz[iv2].x - xcent);

                            if (theta2 < theta1 - PI) theta2 += TWOPI;
                            if (theta2 > theta1 + PI) theta2 -= TWOPI;

                            uv[iv2].u = xyz[iv2].y;
                            uv[iv2].v = theta2;
                            change++;
                        }

                        if (uv[iv0].u < ymin) {
                            theta0 = atan2(xyz[iv0].z - zcent, xyz[iv0].x - xcent);

                            if (theta0 < theta1 - PI) theta0 += TWOPI;
                            if (theta0 > theta1 + PI) theta0 -= TWOPI;

                            uv[iv0].u = xyz[iv0].y;
                            uv[iv0].v = theta0;
                            change++;
                        }
                    }

                    if (uv[iv2].u >= ymin) {
                        theta2 = uv[iv2].v;

                        if (uv[iv0].u < ymin) {
                            theta0 = atan2(xyz[iv0].z - zcent, xyz[iv0].x - xcent);

                            if (theta0 < theta2 - PI) theta0 += TWOPI;
                            if (theta0 > theta2 + PI) theta0 -= TWOPI;

                            uv[iv0].u = xyz[iv0].y;
                            uv[iv0].v = theta0;
                            change++;
                        }

                        if (uv[iv1].u < ymin) {
                            theta1 = atan2(xyz[iv1].z - zcent, xyz[iv1].x - xcent);

                            if (theta1 < theta2 - PI) theta1 += TWOPI;
                            if (theta1 > theta2 + PI) theta1 -= TWOPI;

                            uv[iv1].u = xyz[iv1].y;
                            uv[iv1].v = theta1;
                            change++;
                        }
                    }
                }
                DPRINT1("   change=%d", change);
            }

            /*
             * make sure that all Triangles have outward pointing normals
             */
            nneg = 0;
            for (itri = 0; itri < ntri; itri++) {
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                     - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
                if (area <= 0) {
                    nneg++;
                }
            }

            /*
             * if we got here, then the y-axial projection was successful
             */
            if (nneg == 0) {
                DPRINT0("   y-axial projection was successful");

                status = PRM_OK_AXIAL;
                goto cleanup;
            } else {
                DPRINT2("   there were %d non-positive areas (ntri=%d)", nneg, ntri);
            }
        }
    }

    /*
     * the z-extent is larger than the other two extents, so try an
     *    axis that is aligned with the z-axis
     */
    if ((zmax-zmin) > (xmax-xmin) && (zmax-zmin) > (ymax-ymin)) {
        DPRINT0("potential axis in z-direction");

        /*
         * find the center of the Vertices within the last 1 percent of the x extent
         */
        dz   = 0.01 * (zmax - zmin);

        xmin1 = +HUGEQ;
        xmax1 = -HUGEQ;
        xmin2 = +HUGEQ;
        xmax2 = -HUGEQ;
        ymin1 = +HUGEQ;
        ymax1 = -HUGEQ;
        ymin2 = +HUGEQ;
        ymax2 = -HUGEQ;

        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (xyz[ivrt].z < zmin+dz) {
                xmin1 = MIN(xmin1, xyz[ivrt].x);
                xmax1 = MAX(xmax1, xyz[ivrt].x);
                ymin1 = MIN(ymin1, xyz[ivrt].y);
                ymax1 = MAX(ymax1, xyz[ivrt].y);
            } else if (xyz[ivrt].z > zmax-dz) {
                xmin2 = MIN(xmin2, xyz[ivrt].x);
                xmax2 = MAX(xmax2, xyz[ivrt].x);
                ymin2 = MIN(ymin2, xyz[ivrt].y);
                ymax2 = MAX(ymax2, xyz[ivrt].y);
            }
        }

        xmin = MAX(xmin1, xmin2);
        xmax = MIN(xmax1, xmax2);
        ymin = MAX(ymin1, ymin2);
        ymax = MIN(ymax1, ymax2);

        /*
         * try placing the axis at the center and at each of the extents
         */
        for (pass = 0; pass < 5; pass++) {
            if (pass == 0) {
                zcent = (zmin + zmax) / 2;
                xcent = (xmin + xmax) / 2;
                ycent = (ymin + ymax) / 2;
            } else if (pass == 1) {
                zcent = (zmin + zmax) / 2;
                xcent =  xmin            ;
                ycent = (ymin + ymax) / 2;
            } else if (pass == 2) {
                zcent = (zmin + zmax) / 2;
                xcent =         xmax     ;
                ycent = (ymin + ymax) / 2;
            } else if (pass == 3) {
                zcent = (zmin + zmax) / 2;
                xcent = (xmin + xmax) / 2;
                ycent =  ymin            ;
            } else {
                zcent = (zmin + zmax) / 2;
                xcent = (xmin + xmax) / 2;
                ycent =         ymax     ;
            }
            DPRINT3("xcent=%10.5f  ycent=%10.5f  zcent=%10.5f", xcent, ycent, zcent);

            /*
             * assign u=x and v=theta for every Vertex
             */
            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                uv[ivrt].u = zmin - 1;
            }

            uv[0].u = xyz[0].z;
            uv[0].v = atan2(xyz[0].x - xcent, xyz[0].y - ycent);

            change = 1;
            while (change > 0) {
                change = 0;

                for (itri = 0; itri < ntri; itri++) {
                    iv0 = tri[itri].indices[0] - 1;
                    iv1 = tri[itri].indices[1] - 1;
                    iv2 = tri[itri].indices[2] - 1;


                    if (uv[iv0].u >= zmin) {
                        theta0 = uv[iv0].v;

                        if (uv[iv1].u < zmin) {
                            theta1 = atan2(xyz[iv1].x - xcent, xyz[iv1].y - ycent);

                            if (theta1 < theta0 - PI) theta1 += TWOPI;
                            if (theta1 > theta0 + PI) theta1 -= TWOPI;

                            uv[iv1].u = xyz[iv1].z;
                            uv[iv1].v = theta1;
                            change++;
                        }

                        if (uv[iv2].u < zmin) {
                            theta2 = atan2(xyz[iv2].x - xcent, xyz[iv2].y - ycent);

                            if (theta2 < theta0 - PI) theta2 += TWOPI;
                            if (theta2 > theta0 + PI) theta2 -= TWOPI;

                            uv[iv2].u = xyz[iv2].z;
                            uv[iv2].v = theta2;
                            change++;
                        }
                    }

                    if (uv[iv1].u >= zmin) {
                        theta1 = uv[iv1].v;

                        if (uv[iv2].u < zmin) {
                            theta2 = atan2(xyz[iv2].x - xcent, xyz[iv2].y - ycent);

                            if (theta2 < theta1 - PI) theta2 += TWOPI;
                            if (theta2 > theta1 + PI) theta2 -= TWOPI;

                            uv[iv2].u = xyz[iv2].z;
                            uv[iv2].v = theta2;
                            change++;
                        }

                        if (uv[iv0].u < zmin) {
                            theta0 = atan2(xyz[iv0].x - xcent, xyz[iv0].y - ycent);

                            if (theta0 < theta1 - PI) theta0 += TWOPI;
                            if (theta0 > theta1 + PI) theta0 -= TWOPI;

                            uv[iv0].u = xyz[iv0].z;
                            uv[iv0].v = theta0;
                            change++;
                        }
                    }

                    if (uv[iv2].u >= zmin) {
                        theta2 = uv[iv2].v;

                        if (uv[iv0].u < zmin) {
                            theta0 = atan2(xyz[iv0].x - xcent, xyz[iv0].y - ycent);

                            if (theta0 < theta2 - PI) theta0 += TWOPI;
                            if (theta0 > theta2 + PI) theta0 -= TWOPI;

                            uv[iv0].u = xyz[iv0].z;
                            uv[iv0].v = theta0;
                            change++;
                        }

                        if (uv[iv1].u < zmin) {
                            theta1 = atan2(xyz[iv1].x - xcent, xyz[iv1].y - ycent);

                            if (theta1 < theta2 - PI) theta1 += TWOPI;
                            if (theta1 > theta2 + PI) theta1 -= TWOPI;

                            uv[iv1].u = xyz[iv1].z;
                            uv[iv1].v = theta1;
                            change++;
                        }
                    }
                }
                DPRINT1("   change=%d", change);
            }

            /*
             * make sure that all Triangles have outward pointing normals
             */
            nneg = 0;
            for (itri = 0; itri < ntri; itri++) {
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                     - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
                if (area <= 0) {
                    nneg++;
                }
            }

            /*
             * if we got here, then the z-axial projection was successful
             */
            if (nneg == 0) {
                DPRINT0("   z-axial projection was successful");

                status = PRM_OK_AXIAL;
                goto cleanup;
            } else {
                DPRINT2("   there were %d non-positive areas (ntri=%d)", nneg, ntri);
            }
        }
    }

    status = PRM_NOGLOBALUV;

 cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * floaterParameterization -- use floater's transformation (if possible)        *
 *                                                                              *
 ********************************************************************************
 */
static int
floaterParameterization(int        ntri,     /* (in)   number of Triangles */
                        prmTri     tri[],    /* (in)   array  of Triangles */
                        int        nvrt,     /* (in)   number of Vertices */
                        prmXYZ     xyz[],    /* (in)   array  of Coordinates */
                        prmUV      uv[])     /* (out)  array of  Parameters */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_FLOATER */

    typedef struct {
        int     vrt;                         /* first Vertex in Loop */
        int     nvrt;                        /* number of Vertex in Loop */
        double  alen;                        /* length (in XYZ) of Loop */
    } lup_t;

    lup_t       *lups = NULL;                /* array  of Loops */
    int         nlup;                        /* number of Loops */

    typedef struct {
        int     next;                        /* next boundary Vertex */
        int     lup;                         /* Loop number */
        double  ang;                         /* angle at Loop Vertices */
    } vrt_t;
    vrt_t       *vrts = NULL;                /* array of boundary Vertices */

    sparseMat   amat;                        /* full matrix */
    double      *asmf   = NULL;
    int         *ismf   = NULL;
    double      *urhs   = NULL;
    double      *vrhs   = NULL;
    prmTri      *newTri = NULL;
    int         *ihole  = NULL;
    double      *ahole  = NULL;

    double      ang0, ang1, ang2, ang3, ang4, ang5;
    double      area, amin;
    double      d1sq, d2sq, d3sq;
    double      d01sq, d12sq, d20sq, d23sq;
    double      d30sq, d04sq, d41sq, d15sq, d52sq;
    double      du, dv;
    double      dist;
    double      err=1, erru, errv;
    int         found;
    int         i, j, k, ii, nn, im1, ip1;
    int         imin;
    int         iter, itmax=1000;
    int         itri, oldNtri, tmp;
    int         ivrt, oldIvrt, iv0, iv1, iv2, iv3, iv4, iv5;
    int         lup;
    int         newNtri, count;
    int         nneg;
    int         outer;
    double      xold, yold, zold, told, gold, hold;
    double      xnew, ynew, znew, tnew, unew, vnew;
    int         nhole;

    double      errtol = 0.000001;
    double      frac = 0.25;
    double      omega = 1.6;

    ROUTINE(floaterParameterization);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */
  
    amat.ni = 0;
    amat.is = NULL;

    /*
     * allocate storage that will be used to keep track of the Loop
     *    number and the next Vertex in the Loop
     */
    MALLOC(vrts, vrt_t,  nvrt    );
    MALLOC(lups, lup_t,  MAXLOOPS);
    MALLOC(urhs, double, nvrt    );
    MALLOC(vrhs, double, nvrt    );

    /*
     * for each Vertex that is along a boundary of the Quilt, keep track of
     *    the next Vertex on the boundary
     */
    DPRINT0("finding next Vertex");
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        vrts[ivrt].next = -1;
        vrts[ivrt].lup  = -1;
    }

    nn = 0;
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        if (tri[itri].neigh[0] <= 0) {
            nn++;
            if (vrts[iv1].next == -1) {
                vrts[iv1].next = iv2;
            }
        }
        if (tri[itri].neigh[1] <= 0) {
            nn++;
            if (vrts[iv2].next == -1) {
                vrts[iv2].next = iv0;
            }
        }
        if (tri[itri].neigh[2] <= 0) {
            nn++;
            if (vrts[iv0].next == -1) {
                vrts[iv0].next = iv1;
            }
        }
    }

    /*
     * if there are no boundary Triangles, then return without setting UVtype
     */
    if (nn == 0) {
        DPRINT0("no boundary Tiangles were found");
        status = PRM_NOGLOBALUV;
        goto cleanup;
    }

    /*
     * associate each boundary Vertex with a Loop
     */
    DPRINT0("associating boundary Vertices with Loops");
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        DPRINT2("vrts[%5d].next = %5d", ivrt, vrts[ivrt].next);
    }
    nlup = 0;
    while (1) {

        /*
         * find a Vertex that has a next but not a Loop
         */
        found = 0;
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if ((vrts[ivrt].next >= 0) && (vrts[ivrt].lup < 0)) {
                vrts[ivrt].lup = nlup;
                vrts[ivrt].ang = 0;

                lups[nlup].vrt  = ivrt;
                lups[nlup].nvrt = 1;
                lups[nlup].alen = 0;
                nlup++;

                found++;
                break;
            }
        }

        /*
         * if no unassigned boundary Vertex is found, then exit
         */
        if (found == 0) {
            break;
        }
        DPRINT4("lups[%d].vrt=%x, .nvrt=%d, alen=%f",
                nlup-1, lups[nlup-1].vrt,
                        lups[nlup-1].nvrt,
                        lups[nlup-1].alen);

        /*
         * assign all the Vertices in the Loop
         */
        count = 0;
        oldIvrt = ivrt;
        ivrt    = vrts[ivrt].next;
        while (ivrt != lups[nlup-1].vrt) {
            vrts[ivrt].lup = nlup - 1;
            vrts[ivrt].ang = 0;

            lups[nlup-1].nvrt++;
            lups[nlup-1].alen += sqrt(  SQR(xyz[ivrt].x - xyz[oldIvrt].x)
                                      + SQR(xyz[ivrt].y - xyz[oldIvrt].y)
                                      + SQR(xyz[ivrt].z - xyz[oldIvrt].z));

            oldIvrt = ivrt;
            ivrt    = vrts[ivrt].next;
            count++;
            if (count > nvrt+1) {
                status = PRM_CANNOTFORMLOOP;
                goto cleanup;
            }
        }

        DPRINT4("lup[%d].vrt=%x   .nvrt=%d   .alen=%f",
                nlup-1, lups[nlup-1].vrt, lups[nlup-1].nvrt, lups[nlup-1].alen);
    }

    /*
     * the outer Loop is the Loop with the greatest length (in XYZ)
     */
    DPRINT0("selecting outer Loop");
    outer = 0;
    for (lup = 1; lup < nlup; lup++) {
        if (lups[lup].alen > lups[outer].alen) {
            outer = lup;
        }
    }
    DPRINT1("outer Loop is %d", outer);

    /*
     * count the number of new Triangles we will have when we fill the
     *    inner Loops
     */
    newNtri = ntri;
    for (lup = 0; lup < nlup; lup++) {
        if (lup != outer) {
            newNtri += lups[lup].nvrt - 2;
        }
    }
    DPRINT1("newNtri=%d", newNtri);

    /*
     * get new tri array that is big enough for the new Triangles
     */
    MALLOC(newTri, prmTri, newNtri);
    memcpy(newTri, tri, ntri*sizeof(prmTri));

    /*
     * fill the inner Loops with Triangles
     */
    DPRINT0("filling inner Loops with Triangles");
    oldNtri = ntri;
    if (nlup > 1) {

        /*
         * allocate arrays for closing holes
         */
        MALLOC(ihole, int,    nvrt);
        MALLOC(ahole, double, nvrt);

        /*
         * fill the inner Loops
         */
        for (lup = 0; lup < nlup; lup++) {
            if (lup == outer) continue;
            DPRINT1("filling Loop %d", lup);

            /*
             * make a list of the hole Vertices (ie, the Vertices that make
             *    up the Loop
             */

            nhole    = lups[lup].nvrt;
            ihole[0] = lups[lup].vrt;

            for (i = 1; i < nhole; i++) {
                ihole[i] = vrts[ihole[i-1]].next;
            }
            DPRINT1("nhole=%d", nhole);

            /*
             * one-by-one, cut off the hole Vertex with the smallest angle
             */
            while (nhole > 3) {

                /*
                 * compute the angle at each hole Vertex
                 */
                for (i = 0; i < nhole; i++) {
                    im1 = i - 1;
                    ip1 = i + 1;
                    if (im1 <  0    ) im1 += nhole;
                    if (ip1 == nhole) ip1 -= nhole;

                    d1sq = SQR(xyz[ihole[i  ]].x - xyz[ihole[im1]].x)
                         + SQR(xyz[ihole[i  ]].y - xyz[ihole[im1]].y)
                         + SQR(xyz[ihole[i  ]].z - xyz[ihole[im1]].z);
                    d2sq = SQR(xyz[ihole[ip1]].x - xyz[ihole[i  ]].x)
                         + SQR(xyz[ihole[ip1]].y - xyz[ihole[i  ]].y)
                         + SQR(xyz[ihole[ip1]].z - xyz[ihole[i  ]].z);
                    d3sq = SQR(xyz[ihole[im1]].x - xyz[ihole[ip1]].x)
                         + SQR(xyz[ihole[im1]].y - xyz[ihole[ip1]].y)
                         + SQR(xyz[ihole[im1]].z - xyz[ihole[ip1]].z);

                    ahole[i] = ACOS((d1sq + d2sq - d3sq) / 2 / sqrt(d1sq * d2sq));
                    DPRINT5("i=%5d %5d %5d %5d %10.5f", nhole, i, im1, ip1, ahole[i]);
                }

                /*
                 * minimum angle (which will be cut off)
                 */
                amin = ahole[0];
                imin = 0;

                for (i = 2; i < nhole; i++) {
                    if (ahole[i] < amin) {
                        amin = ahole[i];
                        imin = i;
                    }
                }
                DPRINT2("imin=%d   amin=%f", imin, amin);

                /*
                 * make new Triangle
                 */
                im1 = imin - 1;
                ip1 = imin + 1;
                if (im1 <  0    ) im1 += nhole;
                if (ip1 == nhole) ip1 -= nhole;

                newTri[ntri].indices[0] = ihole[imin] + 1;
                newTri[ntri].indices[1] = ihole[im1 ] + 1;
                newTri[ntri].indices[2] = ihole[ip1 ] + 1;
                newTri[ntri].neigh[0] = 0;
                newTri[ntri].neigh[1] = 0;
                newTri[ntri].neigh[2] = 0;
                ntri++;

                /*
                 * remove imin from hole list
                 */
                nhole--;
                for (j = imin; j < nhole; j++) {
                    ihole[j] = ihole[j+1];
                    ahole[j] = ahole[j+1];
                }
            }

            /*
             * make the final Triangle with the last 3 entries in the hole list
             */
            newTri[ntri].indices[0] = ihole[3] + 1;
            newTri[ntri].indices[1] = ihole[2] + 1;
            newTri[ntri].indices[2] = ihole[1] + 1;
            newTri[ntri].neigh[0] = 0;
            newTri[ntri].neigh[1] = 0;
            newTri[ntri].neigh[2] = 0;
            ntri++;
        }

        /*
         * set up the Triangle neighbor information for the new Triangles
         *    (note the backwards loop since it is most likely that the neighbor
         *    will be near the end of the list)
         */
        for (i = oldNtri; i < ntri; i++) {
            DPRINT1("neighbors for Triangle %d", i);
            iv0 = newTri[i].indices[0];
            iv1 = newTri[i].indices[1];
            iv2 = newTri[i].indices[2];

            for (j = ntri-1; j >= 0; j--) {
                if        (newTri[j].indices[2] == iv1 && newTri[j].indices[1] == iv2) {
                    newTri[i].neigh[0] = j + 1;
                    newTri[j].neigh[0] = i + 1;
                } else if (newTri[j].indices[0] == iv1 && newTri[j].indices[2] == iv2) {
                    newTri[i].neigh[0] = j + 1;
                    newTri[j].neigh[1] = i + 1;
                } else if (newTri[j].indices[1] == iv1 && newTri[j].indices[0] == iv2) {
                    newTri[i].neigh[0] = j + 1;
                    newTri[j].neigh[2] = i + 1;
                }
            }

            for (j = ntri-1; j >= 0; j--) {
                if        (newTri[j].indices[2] == iv2 && newTri[j].indices[1] == iv0) {
                    newTri[i].neigh[1] = j + 1;
                    newTri[j].neigh[0] = i + 1;
                } else if (newTri[j].indices[0] == iv2 && newTri[j].indices[2] == iv0) {
                    newTri[i].neigh[1] = j + 1;
                    newTri[j].neigh[1] = i + 1;
                } else if (newTri[j].indices[1] == iv2 && newTri[j].indices[0] == iv0) {
                    newTri[i].neigh[1] = j + 1;
                    newTri[j].neigh[2] = i + 1;
                }
            }

            for (j = ntri-1; j >= 0; j--) {
                if        (newTri[j].indices[2] == iv0 && newTri[j].indices[1] == iv1) {
                    newTri[i].neigh[2] = j + 1;
                    newTri[j].neigh[0] = i + 1;
                } else if (newTri[j].indices[0] == iv0 && newTri[j].indices[2] == iv1) {
                    newTri[i].neigh[2] = j + 1;
                    newTri[j].neigh[1] = i + 1;
                } else if (newTri[j].indices[1] == iv0 && newTri[j].indices[0] == iv1) {
                    newTri[i].neigh[2] = j + 1;
                    newTri[j].neigh[2] = i + 1;
                }
            }
        }
    }

    /*
     * compute the angle at each boundary Vertex associated with the outer Loop
     */
    DPRINT0("computing boundary angles");
    for (itri = 0; itri < ntri; itri++) {
        iv0 = newTri[itri].indices[0] - 1;
        iv1 = newTri[itri].indices[1] - 1;
        iv2 = newTri[itri].indices[2] - 1;

        if (vrts[iv0].lup >= 0 || vrts[iv1].lup >= 0 || vrts[iv2].lup >= 0) {
            d01sq = SQR(xyz[iv0].x - xyz[iv1].x)
                  + SQR(xyz[iv0].y - xyz[iv1].y)
                  + SQR(xyz[iv0].z - xyz[iv1].z);
            d12sq = SQR(xyz[iv1].x - xyz[iv2].x)
                  + SQR(xyz[iv1].y - xyz[iv2].y)
                  + SQR(xyz[iv1].z - xyz[iv2].z);
            d20sq = SQR(xyz[iv2].x - xyz[iv0].x)
                  + SQR(xyz[iv2].y - xyz[iv0].y)
                  + SQR(xyz[iv2].z - xyz[iv0].z);

            ang0 = ACOS((d20sq + d01sq - d12sq) / 2 / sqrt(d20sq * d01sq));
            ang1 = ACOS((d01sq + d12sq - d20sq) / 2 / sqrt(d01sq * d12sq));
            ang2 = ACOS((d12sq + d20sq - d01sq) / 2 / sqrt(d12sq * d20sq));

            if (vrts[iv0].lup >= 0) {vrts[iv0].ang += ang0;}
            if (vrts[iv1].lup >= 0) {vrts[iv1].ang += ang1;}
            if (vrts[iv2].lup >= 0) {vrts[iv2].ang += ang2;}
        }
    }

    /*
     * lay out the outer Loop (in UV) by marching along the Loop
     */
    DPRINT0("laying out outer Loop");
    ivrt       = lups[outer].vrt;
    uv[ivrt].u = 0;
    uv[ivrt].v = 0;

    xold = xyz[ivrt].x;
    yold = xyz[ivrt].y;
    zold = xyz[ivrt].z;
    gold = 0;
    hold = 0;
    told = 0;

    do {
        ivrt = vrts[ivrt].next;

        xnew = xyz[ivrt].x;
        ynew = xyz[ivrt].y;
        znew = xyz[ivrt].z;

        dist = sqrt(SQR(xnew - xold) + SQR(ynew - yold) + SQR(znew - zold));

        unew = gold + dist * cos(told);
        vnew = hold + dist * sin(told);
        tnew = told + PI - vrts[ivrt].ang;

        uv[ivrt].u = unew;
        uv[ivrt].v = vnew;

        xold = xnew;
        yold = ynew;
        zold = znew;
        gold = unew;
        hold = vnew;
        told = tnew;

    } while (ivrt != lups[outer].vrt);

    /*
     * determine how close we came to closing the Loop (in UV)
     */
    ivrt = lups[outer].vrt;
    dist = sqrt(SQR(uv[ivrt].u) + SQR(uv[ivrt].v));
    DPRINT1("closing dist=%f", dist);

    /*
     * if the Loop is almost closed, close it by linearly adjusting the UV
     */
    DPRINT0("closing outer Loop");
    if (dist < frac*lups[outer].alen) {
        du   = uv[ivrt].u / lups[outer].nvrt;
        dv   = uv[ivrt].v / lups[outer].nvrt;
        DPRINT2("adjusting by du=%f,  dv=%f", du, dv);

        ii = 0;
        do {
            ivrt = vrts[ivrt].next;
            ii++;

            uv[ivrt].u -= (ii * du);
            uv[ivrt].v -= (ii * dv);
        } while (ivrt != lups[outer].vrt);

    /*
     * otherwise, remap the outer Loop to the unit circle
     */
    } else {
        DPRINT0("mapping to circle");
        ivrt = lups[outer].vrt;

        xold = xyz[ivrt].x;
        yold = xyz[ivrt].y;
        zold = xyz[ivrt].z;
        told = 0.0;

        uv[ivrt].u = cos(told);
        uv[ivrt].v = sin(told);

        ivrt = vrts[ivrt].next;
        while (ivrt != lups[outer].vrt) {
            xnew = xyz[ivrt].x;
            ynew = xyz[ivrt].y;
            znew = xyz[ivrt].z;

            dist = sqrt(SQR(xnew - xold) + SQR(ynew - yold)
                                         + SQR(znew - zold));
            told = told + TWOPI * dist / lups[outer].alen;

            uv[ivrt].u = cos(told);
            uv[ivrt].v = sin(told);

            ivrt = vrts[ivrt].next;
            xold = xnew;
            yold = ynew;
            zold = znew;
        }
    }

    /*
     * initialize the amat matrix (which will be used to solve for the UV
     *    at all interior Vertices)
     */
    status = initSmat(nvrt, &amat);
    CHECK_STATUS;

    /*
     * find the mean value weights at each Vertex due to its neighbors.  the
     *    weigt function is the "mean value coordinates" function as described
     *    by M. Floater
     *
     *                   iv3--------iv2---------iv5
     *                     \         /\a5       /
     *                      \       /a2\       /
     *                       \     /    \     /
     *                        \a3 / itri \   /
     *                         \ /a0    a1\ /
     *                         iv0--------iv1
     *                           \     a4 /
     *                            \      /
     *                             \    /
     *                              \  /
     *                               \/
     *                              iv4
     */
    DPRINT0("finding mean value weights");
    for (itri = 0; itri < ntri; itri++) {
        iv0 = newTri[itri].indices[0] - 1;
        iv1 = newTri[itri].indices[1] - 1;
        iv2 = newTri[itri].indices[2] - 1;

        d01sq = SQR(xyz[iv0].x - xyz[iv1].x)
              + SQR(xyz[iv0].y - xyz[iv1].y)
              + SQR(xyz[iv0].z - xyz[iv1].z);
        d12sq = SQR(xyz[iv1].x - xyz[iv2].x)
              + SQR(xyz[iv1].y - xyz[iv2].y)
              + SQR(xyz[iv1].z - xyz[iv2].z);
        d20sq = SQR(xyz[iv2].x - xyz[iv0].x)
              + SQR(xyz[iv2].y - xyz[iv0].y)
              + SQR(xyz[iv2].z - xyz[iv0].z);

        ang0 = ACOS((d20sq + d01sq - d12sq) / 2 / sqrt(d20sq * d01sq));
        ang1 = ACOS((d01sq + d12sq - d20sq) / 2 / sqrt(d01sq * d12sq));
        ang2 = ACOS((d12sq + d20sq - d01sq) / 2 / sqrt(d12sq * d20sq));

        if (newTri[itri].neigh[1] > 0) {
            tmp = newTri[itri].neigh[1] -  1;
            if (       newTri[tmp].neigh[0]-1 == itri) {
                iv3 =  newTri[tmp].indices[0]-1;
            } else if (newTri[tmp].neigh[1]-1 == itri) {
                iv3 =  newTri[tmp].indices[1]-1;
            } else {
                iv3 =  newTri[tmp].indices[2]-1;
            }

            d23sq = SQR(xyz[iv2].x - xyz[iv3].x)
                  + SQR(xyz[iv2].y - xyz[iv3].y)
                  + SQR(xyz[iv2].z - xyz[iv3].z);
            d30sq = SQR(xyz[iv3].x - xyz[iv0].x)
                  + SQR(xyz[iv3].y - xyz[iv0].y)
                  + SQR(xyz[iv3].z - xyz[iv0].z);

            ang3 = ACOS((d30sq + d20sq - d23sq) / 2 / sqrt(d30sq * d20sq));

            status = setSmat(&amat, iv0, iv2,
                             (tan(ang0/2) + tan(ang3/2)) / sqrt(d20sq));
            CHECK_STATUS;
        }

        if (newTri[itri].neigh[2] > 0) {
            tmp = newTri[itri].neigh[2] -  1;
            if (       newTri[tmp].neigh[0]-1 == itri) {
                iv4 =  newTri[tmp].indices[0]-1;
            } else if (newTri[tmp].neigh[1]-1 == itri) {
                iv4 =  newTri[tmp].indices[1]-1;
            } else {
                iv4 =  newTri[tmp].indices[2]-1;
            }

            d04sq = SQR(xyz[iv0].x - xyz[iv4].x)
                  + SQR(xyz[iv0].y - xyz[iv4].y)
                  + SQR(xyz[iv0].z - xyz[iv4].z);
            d41sq = SQR(xyz[iv4].x - xyz[iv1].x)
                  + SQR(xyz[iv4].y - xyz[iv1].y)
                  + SQR(xyz[iv4].z - xyz[iv1].z);

            ang4 = ACOS((d41sq + d01sq - d04sq) / 2 / sqrt(d41sq * d01sq));

            status = setSmat(&amat, iv1, iv0,
                             (tan(ang1/2) + tan(ang4/2)) / sqrt(d01sq));
            CHECK_STATUS;
        }

        if (newTri[itri].neigh[0] > 0) {
            tmp = newTri[itri].neigh[0] -  1;
            if (       newTri[tmp].neigh[0]-1 == itri) {
                iv5 =  newTri[tmp].indices[0]-1;
            } else if (newTri[tmp].neigh[1]-1 == itri) {
                iv5 =  newTri[tmp].indices[1]-1;
            } else {
                iv5 =  newTri[tmp].indices[2]-1;
            }

            d15sq = SQR(xyz[iv1].x - xyz[iv5].x)
                  + SQR(xyz[iv1].y - xyz[iv5].y)
                  + SQR(xyz[iv1].z - xyz[iv5].z);
            d52sq = SQR(xyz[iv5].x - xyz[iv2].x)
                  + SQR(xyz[iv5].y - xyz[iv2].y)
                  + SQR(xyz[iv5].z - xyz[iv2].z);

            ang5 = ACOS((d52sq + d12sq - d15sq) / 2 / sqrt(d52sq * d12sq));

            status = setSmat(&amat, iv2, iv1,
                             (tan(ang2/2) + tan(ang5/2)) / sqrt(d12sq));
            CHECK_STATUS;
        }
    }

    /*
     * set up the final amat and the right-hand sides
     */
    for (i = 0; i < nvrt; i++) {

        /*
         * for interior Vertices, normalize the weights, set the diagonal to 1,
         *    and zero-out the right-hand sides
         */
        if (vrts[i].lup != outer) {
            status = divSmat(&amat, i);
            CHECK_STATUS;
            status = setSmat(&amat, i, i, 1.0);
            CHECK_STATUS;

            urhs[i] = 0;
            vrhs[i] = 0;

            uv[i].u = 0;
            uv[i].v = 0;

        /*
         * for boundary Vertices, zero-out all amat elements, set the diagonal
         *    to 1, and store the boundary values in the RHS
         */
        } else {
            status = diagSmat(&amat, i);
            CHECK_STATUS;
            status = setSmat(&amat, i, i, 1.0);
            CHECK_STATUS;

            urhs[i] = uv[i].u;
            vrhs[i] = uv[i].v;
        }
    }

    /*
     * count the number of non-zero entries in amat
     */
    nn = countSmat(&amat);

    /*
     * allocate matrices for sparse-matrix form.  this form is as
     *    described in 'Numerical Recipes'
     */
    DPRINT0("setting up sparse matrix");
    MALLOC(asmf, double, nn+2);
    MALLOC(ismf, int,    nn+2);

    /*
     * store the matrix in sparse-matrix form (which will make
     *    the iterations below much more efficient)
     */
    for (j = 0; j < nvrt; j++) {
      status = getSmat(&amat, j, j, &asmf[j]);
      CHECK_STATUS;
    }
  
    ismf[0] = nvrt + 2;
  
    k = nvrt + 1;
    for (i = 0; i < nvrt; i++) {
      status = fillSmat(&amat, i, &k, asmf, ismf);
      CHECK_STATUS;
      ismf[i+1] = k + 1;
    }

    /*
     * solve for the Gs and Hs using successive-over-relaxation
     */
    DPRINT0("solving sparse matrix");
    for (iter = 0; iter < itmax; iter++) {

        /*
         * apply successive-over-relaxation
         */
        for (i = 0; i < nvrt; i++) {
            du = urhs[i] - asmf[i] * uv[i].u;
            dv = vrhs[i] - asmf[i] * uv[i].v;

            for (k = ismf[i]; k < ismf[i+1]; k++) {
                j    = ismf[k];
                du  -= asmf[k] * uv[j].u;
                dv  -= asmf[k] * uv[j].v;
            }

            uv[i].u += omega * du / asmf[i];
            uv[i].v += omega * dv / asmf[i];
        }

        /*
         * compute the norm of the residual
         */
        err = 0;
        for (i = 0; i < nvrt; i++) {
            erru = urhs[i] - asmf[i] * uv[i].u;
            errv = vrhs[i] - asmf[i] * uv[i].v;

            for (k = ismf[i]; k < ismf[i+1]; k++) {
                j     = ismf[k];
                erru -= asmf[k] * uv[j].u;
                errv -= asmf[k] * uv[j].v;
            }

            err += SQR(erru) + SQR(errv);
        }
        err = sqrt(err);

        /*
         * exit if converged
         */
        DPRINT2("iter=%5d   err=%15.8e", iter, err);
        if (err < errtol) break;
    }

    /*
     * find the number of Triangles with negative areas (in UV)
     */
    DPRINT0("checking final Triangles");
    nneg = 0;
    for (itri = 0; itri < oldNtri; itri++) {
        iv0 = newTri[itri].indices[0] - 1;
        iv1 = newTri[itri].indices[1] - 1;
        iv2 = newTri[itri].indices[2] - 1;

        area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
             - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);

        if (area < 0) nneg++;
    }
    DPRINT2("floater parameterization has nneg=%d of %d", nneg, ntri);

    /*
     * if we converged and none of the Triangles have negative areas,
     *    then set the UVtype
     */
    if (err < errtol && nneg == 0) {
        status = PRM_OK_FLOATER;
    } else {
        status = PRM_NOGLOBALUV;
    }

 cleanup:
    freeSmat(&amat);
    FREE(vrts);
    FREE(lups);
    FREE(urhs);
    FREE(vrhs);
    FREE(asmf);
    FREE(ismf);
    FREE(newTri);
    FREE(ihole);
    FREE(ahole);

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * jonesTransformation -- use jones' transformation (if possible)               *
 *                                                                              *
 ********************************************************************************
 */
static int
jonesTransformation(int      ntri,           /* (in)   number of Triangles */
                    prmTri   tri[],          /* (in)   array  of Triangles */
                    prmUVF   uvf[],          /* (in)   array  of UVF for Triangles */
                    int      nvrt,           /* (in)   number of Vertices */
                    prmUV    uv[])           /* (out)  array  of Parameters */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_JONES */

    int         *itmap = NULL;               /* Triangle mapping status */
                                             /*          (0=not mapped, 1=scheduled, 2=mapped) */
    int         *ivmap = NULL;               /* Vertex   mapping status */
                                             /*          (0=not mapped, 1=mapped) */
    int         *list  = NULL;               /* FIFO list of Triangles to process */

    double      area;
    int         itri, jtri;
    int         ivrt, iv0, iv1, iv2;
    int         jlist, nlist, nneg, side;

    ROUTINE(jonesTransformation);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

   /*
     * get the arrays that will keep track of the mapping status of
     *    the Vertices and Triangles
     */
    MALLOC(ivmap, int, nvrt);
    MALLOC(itmap, int, ntri);
    MALLOC(list,  int, ntri);

    /*
     * mark all Vertices and Triangles as not mapped
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        ivmap[ivrt] = 0;
    }
    for (itri = 0; itri < ntri; itri++) {
        itmap[itri] = 0;
    }

    /*
     * set up the coordinates for the first Triangle by copying UV
     *    from the first Triangle
     */
    itri = 0;

    iv0 = tri[itri].indices[0] - 1;
    iv1 = tri[itri].indices[1] - 1;
    iv2 = tri[itri].indices[2] - 1;

    uv[iv0].u = uvf[itri].u0;
    uv[iv0].v = uvf[itri].v0;

    uv[iv1].u = uvf[itri].u1;
    uv[iv1].v = uvf[itri].v1;

    uv[iv2].u = uvf[itri].u2;
    uv[iv2].v = uvf[itri].v2;

    ivmap[iv0] = 1;
    ivmap[iv1] = 1;
    ivmap[iv2] = 1;

    itmap[itri] = 2;

    /*
     * flip this Triangle over if its area was negative
     */
    area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
         - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);

    if (area < 0) {
        uv[iv0].u *= -1;
        uv[iv1].u *= -1;
        uv[iv2].u *= -1;

        DPRINT1("triangle flipped because area = %f", area);
    }

    /*
     * initialize the first-in-first-out list of Triangles to process
     */
    nlist = 0;                  /* number of entries in list */
    jlist = 0;                  /* next list item to process */

    for (side = 0; side < 3; side++) {
        jtri = tri[itri].neigh[side] - 1;

        if (jtri >= 0 && itmap[jtri] == 0) {
            itmap[jtri]   = 1;
            list[nlist++] = jtri;
        }
    }

    /*
     * process each list item (first-in-first-out)
     */
    while (jlist < nlist) {
        itri = list[jlist++];
        iv0  = tri[itri].indices[0] - 1;
        iv1  = tri[itri].indices[1] - 1;
        iv2  = tri[itri].indices[2] - 1;

        if        ((ivmap[iv0] == 1) && (ivmap[iv1] == 1) && (ivmap[iv2] == 1)) {
            itmap[itri] = 2;
        } else if ((ivmap[iv0] == 1) && (ivmap[iv1] == 1)) {
            layDownTriangle2(uv, uvf[itri].u0, uvf[itri].v0, iv0,
                                 uvf[itri].u1, uvf[itri].v1, iv1,
                                 uvf[itri].u2, uvf[itri].v2, iv2);

            ivmap[iv2 ] = 1;
            itmap[itri] = 2;
        } else if ((ivmap[iv1] == 1) && (ivmap[iv2] == 1)) {
            layDownTriangle2(uv, uvf[itri].u1, uvf[itri].v1, iv1,
                                 uvf[itri].u2, uvf[itri].v2, iv2,
                                 uvf[itri].u0, uvf[itri].v0, iv0);

            ivmap[iv0 ] = 1;
            itmap[itri] = 2;
        } else if ((ivmap[iv2] == 1) && (ivmap[iv0] == 1)) {
            layDownTriangle2(uv, uvf[itri].u2, uvf[itri].v2, iv2,
                                 uvf[itri].u0, uvf[itri].v0, iv0,
                                 uvf[itri].u1, uvf[itri].v1, iv1);

            ivmap[iv1 ] = 1;
            itmap[itri] = 2;
        }

        /*
         * add the neighbors to the list to process
         */
        for (side = 0; side < 3; side++) {
            jtri = tri[itri].neigh[side] - 1;

            if (jtri >= 0 && itmap[jtri] == 0) {
                itmap[jtri]   = 1;
                list[nlist++] = jtri;
            }
        }
    }

    /*
     * count the number of Triangles with a positive area
     */
    nneg = 0;
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
             - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
        if (area < EPS12) nneg++;
    }
    DPRINT2("jones transformation has nneg=%d of %d", nneg, ntri);

    /*
     * if all the areas are positive, set UVtype
     *    and use the newly-slitted tessellation
     */
    if (nneg == 0) {
        status = PRM_OK_JONES;
    } else {

#ifdef GRAFIC
        plotTris(ntri, tri, nvrt, uv,
                 "~U~V~Jones transformation with negative areas");
#endif

        status = PRM_NOGLOBALUV;
    }

 cleanup:
    FREE(ivmap);
    FREE(itmap);
    FREE(list );

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * planarProjection -- use planar projection (if possible)                      *
 *                                                                              *
 ********************************************************************************
 */
static int
planarProjection(int      ntri,              /* (in)   number of Triangles */
                 prmTri   tri[],             /* (in)   array  of Triangles */
                 int      nvrt,              /* (in)   number of Vertices */
                 prmXYZ   xyz[],             /* (in)   array  of Coordinates */
                 prmUV    uv[])              /* (out)  array  of Parameters */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_PLANAR */

    double      area, areax, areay, areaz, areamax;
    int         itri;
    int         ivrt, iv0, iv1, iv2;
    double      norm[3];
    int         nneg, npos;

    ROUTINE(planarProjection);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * find the total area vector
     */
    areax   = 0;
    areay   = 0;
    areaz   = 0;
    areamax = 0;

    for (itri = 0; itri < ntri; itri++) {
        surfaceNormal(tri[itri].indices[0]-1,
                      tri[itri].indices[1]-1,
                      tri[itri].indices[2]-1, xyz, &area, norm);

        areax += area * norm[0];
        areay += area * norm[1];
        areaz += area * norm[2];

        areamax = MAX(areamax, area);
    }

    DPRINT3("areax=%10.5f   areay=%10.5f   areaz=%10.5f",
            areax, areay, areaz);

    /*
     * exit if the area is very small (which probably means that there is
     *    no projection direction
     */
    area  = sqrt(SQR(areax) + SQR(areay) + SQR(areaz));
    if (area < areamax) {
        goto cleanup;
    } else {
        areax /= area;
        areay /= area;
        areaz /= area;
    }

    DPRINT3("areax=%10.5f, areay=%10.5f, areaz=%10.5f", areax, areay, areaz);

    /*
     * convert all xyz to uv
     */
    if        (fabs(areax) <= fabs(areay) && fabs(areax) <= fabs(areaz)) {
        DPRINT0("typeX::  y,z->U  x->V");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = areay * xyz[ivrt].z - areaz * xyz[ivrt].y;
            uv[ivrt].v =         xyz[ivrt].x;
        }
    } else if (fabs(areay) <= fabs(areaz) && fabs(areay) <= fabs(areax)) {
        DPRINT0("typeY::  z,x->U  y->V");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = areaz * xyz[ivrt].x - areax * xyz[ivrt].z;
            uv[ivrt].v =         xyz[ivrt].y;
        }
    } else {
        DPRINT0("typeZ::  x,y->U  z->V");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = areax * xyz[ivrt].y - areay * xyz[ivrt].x;
            uv[ivrt].v =         xyz[ivrt].z;
        }
    }

    /*
     * make sure that all Triangles have outward pointing normals
     */
    nneg = 0;
    npos = 0;
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
             - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
        if (area < 0) {
            nneg++;
        } else if (area > 0) {
            npos++;
        }
    }
    DPRINT3("ntri=%d   nneg=%d   npos=%d", ntri, nneg, npos);

    /*
     * if we got here, then the projection to the rotated z-plane was successful
     */
    if (npos == ntri) {
        DPRINT0("   positive planar projection was successful");
        status = PRM_OK_PLANAR;
    } else if (nneg == ntri) {
        DPRINT0("   negative planar projection was successful");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].v = -(uv[ivrt].v);
        }
        status = PRM_OK_PLANAR;
    } else {
        DPRINT0("   planar projection failed");
    }

 cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * simpleProjection -- use simple projection (if possible)                      *
 *                                                                              *
 ********************************************************************************
 */
static int
simpleProjection(int      ntri,              /* (in)   number of Triangles */
                 prmTri   tri[],             /* (in)   array  of Triangles */
                 int      nvrt,              /* (in)   number of Vertices */
                 prmXYZ   xyz[],             /* (in)   array  of Coordinates */
                 prmUV    uv[])              /* (out)  array  of Parameters */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_SIMPLE */

    double      area, norm[3];
    int         itri, ivrt, iposx, inegx, iposy, inegy, iposz, inegz;
    int         iv0, iv1, iv2, npos, nneg;

    double      sin10 = 0.17365;

    ROUTINE(simpleProjection);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * initially all projections are assumed to be possible
     */
    iposx = 1;
    inegx = 1;
    iposy = 1;
    inegy = 1;
    iposz = 1;
    inegz = 1;

    /*
     * loop through Triangles and eliminate possible projections
     */
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        surfaceNormal(iv0, iv1, iv2, xyz, &area, norm);

        if (norm[0] < +sin10) iposx = 0;
        if (norm[0] > -sin10) inegx = 0;
        if (norm[1] < +sin10) iposy = 0;
        if (norm[1] > -sin10) inegy = 0;
        if (norm[2] < +sin10) iposz = 0;
        if (norm[2] > -sin10) inegz = 0;
    }
    DPRINT6("iposx=%d inegx=%d  iposy=%d inegy=%d  iposz=%d inegz=%d",
            iposx, inegx, iposy, inegy, iposz, inegz);

    /*
     * perform the first avalable transformation
     */
    if        (iposx == 1) {
        DPRINT0("   simple +x projection");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = xyz[ivrt].y;
            uv[ivrt].v = xyz[ivrt].z;
        }
        status = PRM_OK_SIMPLE;
    } else if (inegx == 1) {
        DPRINT0("   simple -x projection");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = xyz[ivrt].z;
            uv[ivrt].v = xyz[ivrt].y;
        }
        status = PRM_OK_SIMPLE;
    } else if (iposy == 1) {
        DPRINT0("   simple +y projection");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = xyz[ivrt].z;
            uv[ivrt].v = xyz[ivrt].x;
        }
        status = PRM_OK_SIMPLE;
    } else if (inegy == 1) {
        DPRINT0("   simple -y projection");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = xyz[ivrt].x;
            uv[ivrt].v = xyz[ivrt].z;
        }
        status = PRM_OK_SIMPLE;
    } else if (iposz == 1) {
        DPRINT0("   simple +z projection");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = xyz[ivrt].x;
            uv[ivrt].v = xyz[ivrt].y;
        }
        status = PRM_OK_SIMPLE;
    } else if (inegz == 1) {
        DPRINT0("   simple -z projection");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            uv[ivrt].u = xyz[ivrt].y;
            uv[ivrt].v = xyz[ivrt].x;
        }
        status = PRM_OK_SIMPLE;
    } else {
        DPRINT0("   simple projection failed");
    }

    /*
     * make sure that all Triangles have outward pointing normals
     */
    if (status == PRM_OK_SIMPLE) {
        nneg = 0;
        npos = 0;
        for (itri = 0; itri < ntri; itri++) {
            iv0 = tri[itri].indices[0] - 1;
            iv1 = tri[itri].indices[1] - 1;
            iv2 = tri[itri].indices[2] - 1;

            area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                 - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
            if (area < 0) {
                nneg++;
            } else if (area > 0) {
                npos++;
            }
        }
        DPRINT3("ntri=%d   nneg=%d   npos=%d", ntri, nneg, npos);

        if (nneg != 0) {
            status = PRM_NOGLOBALUV;
            goto cleanup;
        }
    }

cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * unroll -- unroll (if possible)                                               *
 *                                                                              *
 ********************************************************************************
 */
static int
unroll(int      ntri,                        /* (in)   number of Triangles */
       prmTri   tri[],                       /* (in)   array  of Triangles */
       int      nvrt,                        /* (in)   number of Vertices */
       prmXYZ   xyz[],                       /* (in)   array  of Coordinates */
       prmUV    uv[])                        /* (out)  array  of Parameters */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_UNROLLING  */

    int         *itmap = NULL;               /* Triangle mapping status */
                                             /*          (0=not mapped, 1=scheduled, 2=mapped) */
    int         *ivmap = NULL;               /* Vertex   mapping status */
                                             /*          (0=not mapped, 1=mapped) */
    int         *list  = NULL;               /* FIFO list of Triangles to process */

    double      area;
    int         done[5][5][5];
    double      error, error_min;
    double      ubar, vbar, unew, vnew;
    int         itheta;
    int         i, j, k;
    int         ibeg;
    int         itri, jtri;
    int         ivrt, iv0, iv1, iv2;
    int         jlist;
    int         nlist;
    int         nneg;
    int         side;
    double      theta, theta_min, costht, sintht;
    double      xmin, ymin, zmin;
    double      xmax, ymax, zmax;

    ROUTINE(unroll);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * get the arrays that will keep track of the mapping status of
     *    the Vertices and Triangles
     */
    MALLOC(ivmap, int, nvrt);
    MALLOC(itmap, int, ntri);
    MALLOC(list,  int, ntri);

    /*
     * start unrolling at one Triangle in each "region" of the Quilt's
     *    bounding box
     */
    for (k = 0; k < 5; k++) {
        for (j = 0; j < 5; j++) {
            for (i = 0; i < 5; i++) {
                done[i][j][k] = 0;
            }
        }
    }

    xmin = +HUGEQ;
    xmax = -HUGEQ;
    ymin = +HUGEQ;
    ymax = -HUGEQ;
    zmin = +HUGEQ;
    zmax = -HUGEQ;

    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        xmin = MIN(xmin, xyz[ivrt].x);
        xmax = MAX(xmax, xyz[ivrt].x);
        ymin = MIN(ymin, xyz[ivrt].y);
        ymax = MAX(ymax, xyz[ivrt].y);
        zmin = MIN(zmin, xyz[ivrt].z);
        zmax = MAX(zmax, xyz[ivrt].z);
    }

    for (ibeg = 0; ibeg < ntri; ibeg++) {
        iv0 = tri[ibeg].indices[0] - 1;

        i = MINMAX(0, (int)(4.99 * (xyz[iv0].x - xmin) / (xmax - xmin)), 4);
        j = MINMAX(0, (int)(4.99 * (xyz[iv0].y - ymin) / (ymax - ymin)), 4);
        k = MINMAX(0, (int)(4.99 * (xyz[iv0].z - zmin) / (zmax - zmin)), 4);

        if (done[i][j][k] == 0) {
            done[i][j][k] = 1;

            /*
             * mark all Vertices and Triangles as not mapped
             */
            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                ivmap[ivrt] = 0;
            }
            for (itri = 0; itri < ntri; itri++) {
                itmap[itri] = 0;
            }

            /*
             * set up the coordinates for the first Triangle
             */
            itri = ibeg;
            iv0 = tri[itri].indices[0] - 1;
            iv1 = tri[itri].indices[1] - 1;
            iv2 = tri[itri].indices[2] - 1;

            uv[   iv0].u = 0.0;
            uv[   iv0].v = 0.0;
            ivmap[iv0]   = 1;

            uv[   iv1].u = sqrt(  SQR(xyz[iv0].x - xyz[iv1].x)
                                + SQR(xyz[iv0].y - xyz[iv1].y)
                                + SQR(xyz[iv0].z - xyz[iv1].z));
            uv[   iv1].v = 0.0;
            ivmap[iv1]   = 1;

            layDownTriangle(xyz, uv, iv0, iv1, iv2);
            ivmap[iv2]   = 1;
            itmap[itri]   = 2;

            /*
             * initialize the first-in-first-out list of Triangles to process
             */
            nlist = 0;                  /* number of entries in list */
            jlist = 0;                  /* next list item to process */

            for (side = 0; side < 3; side++) {
                jtri = tri[itri].neigh[side] - 1;

                if (jtri >= 0 && itmap[jtri] == 0) {
                    itmap[jtri]   = 1;
                    list[nlist++] = jtri;
                }
            }

            /*
             * process each list item (first-in-first-out)
             */
            while (jlist < nlist) {
                itri = list[jlist++];
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                if        ((ivmap[iv0] == 1) && (ivmap[iv1] == 1) && (ivmap[iv2] == 1)) {
                    itmap[itri] = 2;
                } else if ((ivmap[iv0] == 1) && (ivmap[iv1] == 1)) {
                    layDownTriangle(xyz, uv, iv0, iv1, iv2);
                    ivmap[iv2] = 1;
                    itmap[itri] = 2;
                } else if ((ivmap[iv1] == 1) && (ivmap[iv2] == 1)) {
                    layDownTriangle(xyz, uv, iv1, iv2, iv0);
                    ivmap[iv0] = 1;
                    itmap[itri] = 2;
                } else if ((ivmap[iv2] == 1) && (ivmap[iv0] == 1)) {
                    layDownTriangle(xyz, uv, iv2, iv0, iv1);
                    ivmap[iv1] = 1;
                    itmap[itri] = 2;
                }

                /*
                 * add the neighbors to the list to process
                 */
                for (side = 0; side < 3; side++) {
                    jtri = tri[itri].neigh[side] - 1;

                    if (jtri >= 0 && itmap[jtri] == 0) {
                        itmap[jtri]   = 1;
                        list[nlist++] = jtri;
                    }
                }
            }

            /*
             * rotate and translate u and v so that their centroid is at
             *    the origin and so that a "best-fit line" is in the u direction
             */
            ubar = 0.0;
            vbar = 0.0;

            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                ubar += uv[ivrt].u;
                vbar += uv[ivrt].v;
            }

            ubar /= nvrt;
            vbar /= nvrt;

            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                uv[ivrt].u -= ubar;
                uv[ivrt].v -= vbar;
            }

            error_min = +HUGEQ;
            theta_min = 0.0;

            for (itheta = -19; itheta <= 20; itheta++) {
                theta = PI * itheta / 40.0;
                error = 0.0;

                costht = cos(theta);
                sintht = sin(theta);

                for (ivrt = 0; ivrt < nvrt; ivrt++) {
                    error += SQR(uv[ivrt].u * sintht + uv[ivrt].v * costht);
                }

                if (error < error_min) {
                    error_min = error;
                    theta_min = theta;
                }
            }

            costht = cos(theta_min);
            sintht = sin(theta_min);

            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                unew = uv[ivrt].u * costht - uv[ivrt].v * sintht;
                vnew = uv[ivrt].u * sintht + uv[ivrt].v * costht;

                uv[ivrt].u = unew;
                uv[ivrt].v = vnew;
            }

            /*
             * count the number of Triangles with a positive area
             */
            nneg = 0;
            for (itri = 0; itri < ntri; itri++) {
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                     - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);
                if (area < EPS12) nneg++;
            }
            DPRINT2("unrolling from Triangle %d has nneg = %d", ibeg, nneg);

            /*
             * if all the areas are positive, set UVtype
             *    and copy the temporary uv values back into myQuilt->uv
             */
            if (nneg == 0) {
                status = PRM_OK_UNROLLING;
                goto cleanup;
            }
        }
    }

    /*
     * if we got here, there is no good unrolling, so simply return
     *    without updating UVtype
     */

    DPRINT0("no valid unrolling was found");
    status = PRM_NOGLOBALUV;

 cleanup:
    FREE(ivmap);
    FREE(itmap);
    FREE(list );

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * fillHoles -- fill holes in a Triangulation                                   *
 *                                                                              *
 ********************************************************************************
 */
static int
fillHoles(int      ntri,                     /* (in)   number of Triangles */
          prmTri   tri[],                    /* (in)   array  of Triangles */
          int      nvrt,                     /* (in)   number of Vertices */
          prmUV    uv[],                     /* (in)   array  of parameters */
          int      *newNtri,                 /* (out)  new number of Tri */
          prmTri   *newTri[])                /* (out)  new array  of Tri */
{
    int      status = EGADS_SUCCESS;         /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    typedef struct {
        int     ibeg;                        /* first Vertex in Loop */
        int     ilvl;                        /* level (starts at 0) */
        int     ipar;                        /* parent Loop (or -1 for outer Loop */
        int     nvrt;                        /* number of Vertices in Loop */
    } lup_t;

    int         nlup;                        /* number of Loops */
    lup_t       *lup     = NULL;             /* array  of Loops */

    int         *prev    = NULL;
    int         *lupid   = NULL;
    int         *ihole   = NULL;
    int         *tempTri = NULL;
    int         *cont    = NULL;
    double      *uvcont  = NULL;

    double      uu, vv;
    int         found, ilup, ivrt, jlup, i, j, k, nn, ncont, ipass, nfig8;
    int         itri, oldNtri, iv0, iv1, iv2;
    fillArea    fast;


    ROUTINE(fillHoles);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * allocate storage that will be used to keep track of the Loop
     *    number and the previous Vertex in the Loop
     */
    MALLOC(prev,    int,    nvrt    );
    MALLOC(lupid,   int,    nvrt    );
    MALLOC(ihole,   int,    nvrt    );
    MALLOC(lup,     lup_t,  MAXLOOPS);
    MALLOC(tempTri, int,    3*nvrt  );
    MALLOC(cont,    int,    MAXLOOPS);
    MALLOC(uvcont,  double, 2*nvrt  );

    /*
     * for each Vertex that is along a boundary, keep track of
     *    the prev Vertex on the boundary
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        prev[ ivrt] = -1;
        lupid[ivrt] = -1;
    }

    nn = 0;
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        if (tri[itri].neigh[0] <= 0) {
            nn++;
            if (prev[iv2] == -1) {
                prev[iv2] = iv1;
            }
        }
        if (tri[itri].neigh[1] <= 0) {
            nn++;
            if (prev[iv0] == -1) {
                prev[iv0] = iv2;
            }
        }
        if (tri[itri].neigh[2] <= 0) {
            nn++;
            if (prev[iv1] == -1) {
                prev[iv1] = iv0;
            }
        }
    }

    /*
     * if there are no boundary Triangles, then make a copy
     *    of the inputs and return
     */
    if (nn == 0) {
        DPRINT0("no boundary Triangles");
        MALLOC(*newTri, prmTri, ntri);
        memcpy(*newTri, tri,    ntri*sizeof(prmTri));
        *newNtri = ntri;
        goto cleanup;
    }

    /*
     * associate each boundary Vertex with a Loop
     */
    nlup = 0;
    while (1) {

        /*
         * find a Vertex that has a prev but not a Loop
         */
        found = 0;
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if ((prev[ivrt] >= 0) && (lupid[ivrt] < 0)) {
                lupid[ivrt] = nlup;

                lup[nlup].ibeg =  ivrt;
                lup[nlup].ilvl =  0;
                lup[nlup].ipar = -1;
                lup[nlup].nvrt =  1;
                nlup++;

                found++;
                break;
            }
        }

        /*
         * if no unassigned boundary Vertex is found, then exit
         */
        if (found == 0) {
            break;
        }

        /*
         * assign all the Vertices in the Loop
         */
        ivrt = prev[ivrt];
        while (ivrt != lup[nlup-1].ibeg) {
            lupid[ivrt] = nlup - 1;

            lup[nlup-1].nvrt++;

            ivrt = prev[ivrt];
        }
    }

    /*
     * if there is only one Loop, then just make a copy
     *    of the inputs and return
     */
    if (nlup <= 1) {
        DPRINT0("only 1 (outer) Loop");
        MALLOC(*newTri, prmTri, ntri);
        memcpy(*newTri, tri,    ntri*sizeof(prmTri));
        *newNtri =ntri;
        goto cleanup;
    }

    /*
     * determine the level number (starting at 0 for outer Loop) for each Loop
     */
    for (ilup = 0; ilup < nlup; ilup++) {
        ivrt = lup[ilup].ibeg;
        uu   = uv[ivrt].u;
        vv   = uv[ivrt].v;

        for (jlup = 0; jlup < nlup; jlup++) {
            if (jlup == ilup) continue;

            if (insideLoop(uu, vv, lup[jlup].ibeg, uv, prev) > 0) {
                lup[ilup].ilvl ++;
            }
        }
    }

    /*
     * find the parent of each Loop
     */
    for (ilup = 0; ilup < nlup; ilup++) {
        ivrt = lup[ilup].ibeg;
        uu   = uv[ivrt].u;
        vv   = uv[ivrt].v;

        for (jlup = 0; jlup < nlup; jlup++) {
            if (jlup == ilup) continue;

            if (lup[ilup].ilvl == lup[jlup].ilvl+1) {
                if (insideLoop(uu, vv, lup[jlup].ibeg, uv, prev) > 0) {
                    lup[ilup].ipar = jlup;
                    break;
                }
            }
        }

        DPRINT5("ilup=%2d, .ibeg=%4d, .nvrt=%4d, .ilvl=%2d, .ipar=%2d", ilup,
                lup[ilup].ibeg,
                lup[ilup].nvrt,
                lup[ilup].ilvl,
                lup[ilup].ipar);
        ivrt = lup[ilup].ibeg;
        for (i = 0; i < lup[ilup].nvrt; i++) {
            DPRINT3("   %5d %10.5f %10.5f", ivrt, uv[ivrt].u, uv[ivrt].v);
            ivrt = prev[ivrt];
        }
    }

    /*
     * get new tri array that is big enough for the new Triangles
     */
    *newNtri = ntri + 2 * nvrt;
    MALLOC(*newTri, prmTri, *newNtri);
    memcpy(*newTri, tri,        ntri*sizeof(prmTri));

    /*
     * fill the Loops that have "odd" level numbers
     */
    oldNtri  = ntri;
    *newNtri = ntri;

    for (ilup = 0; ilup < nlup; ilup++) {
        if (lup[ilup].ilvl%2  == 0) continue;

        DPRINT1("filling Loop %d", ilup);
        ncont     = 1;
        cont[0]   = lup[ilup].nvrt;
        uvcont[0] = 0;
        uvcont[1] = 0;
        k         = 1;

        ihole[k] = lup[ilup].ibeg;
        for (i = 1; i <= cont[0]; i++) {
            uvcont[2*k  ] = uv[  ihole[k]].u;
            uvcont[2*k+1] = uv[  ihole[k]].v;
            ihole[   k+1] = prev[ihole[k]];
            k++;
        }

        /*
         * splice in each of the children Loops
         */
        for (jlup = 0; jlup < nlup; jlup++) {
            if (lup[jlup].ipar == ilup) {

                DPRINT1("   contains Loop %d", jlup);
                ncont++;
                cont[ncont-1] = lup[jlup].nvrt;

                ihole[k] = lup[jlup].ibeg;
                for (i = 0; i < cont[ncont-1]; i++) {
                    uvcont[2*k  ] = uv[  ihole[k]].u;
                    uvcont[2*k+1] = uv[  ihole[k]].v;
                    ihole[   k+1] = prev[ihole[k]];
                    k++;
                }
            }
        }

        /*
         * fill the Loop (with possible holes)
         */
        fast.pts   = NULL;
        fast.segs  = NULL;
        fast.front = NULL;
        for (ipass = 0; ipass <= 1; ipass++) {
            status = EG_fillArea(ncont, cont, uvcont, tempTri, &nfig8, ipass,
                                  &fast);
            DPRINT1("      EG_fillArea: status=%d", status);
            if (status > 0) break;
        }
        if (fast.segs  != NULL) EG_free(fast.segs);
        if (fast.pts   != NULL) EG_free(fast.pts);
        if (fast.front != NULL) EG_free(fast.front);

        /*
         * add these Triangle to the end of the list
         */
        for (itri = 0; itri < status; itri++) {
            (*newTri)[*newNtri].indices[0] = ihole[tempTri[3*itri  ]] + 1;
            (*newTri)[*newNtri].indices[1] = ihole[tempTri[3*itri+1]] + 1;
            (*newTri)[*newNtri].indices[2] = ihole[tempTri[3*itri+2]] + 1;
            (*newTri)[*newNtri].neigh[  0] = 0;
            (*newTri)[*newNtri].neigh[  1] = 0;
            (*newTri)[*newNtri].neigh[  2] = 0;
            (*newNtri)++;
        }
    }

    /*
     * set up the Triangle neighbor information for the new Triangles
     *    (note the backwards loop since it is most likely that the neighnor
     *    will be near the end of the list)
     */
    for (i = oldNtri; i < *newNtri; i++) {
        iv0 = (*newTri)[i].indices[0];
        iv1 = (*newTri)[i].indices[1];
        iv2 = (*newTri)[i].indices[2];

        for (j = *newNtri-1; j >= 0; j--) {
            if        ((*newTri)[j].indices[2] == iv1 &&
                       (*newTri)[j].indices[1] == iv2   ) {
                (*newTri)[i].neigh[0] = j + 1;
                (*newTri)[j].neigh[0] = i + 1;
                break;
            } else if ((*newTri)[j].indices[0] == iv1 &&
                       (*newTri)[j].indices[2] == iv2   ) {
                (*newTri)[i].neigh[0] = j + 1;
                (*newTri)[j].neigh[1] = i + 1;
                break;
            } else if ((*newTri)[j].indices[1] == iv1 &&
                       (*newTri)[j].indices[0] == iv2   ) {
                (*newTri)[i].neigh[0] = j + 1;
                (*newTri)[j].neigh[2] = i + 1;
                break;
            }
        }

        for (j = *newNtri-1; j >= 0; j--) {
            if        ((*newTri)[j].indices[2] == iv2 &&
                       (*newTri)[j].indices[1] == iv0   ) {
                (*newTri)[i].neigh[1] = j + 1;
                (*newTri)[j].neigh[0] = i + 1;
                break;
            } else if ((*newTri)[j].indices[0] == iv2 &&
                       (*newTri)[j].indices[2] == iv0   ) {
                (*newTri)[i].neigh[1] = j + 1;
                (*newTri)[j].neigh[1] = i + 1;
                break;
            } else if ((*newTri)[j].indices[1] == iv2 &&
                       (*newTri)[j].indices[0] == iv0   ) {
                (*newTri)[i].neigh[1] = j + 1;
                (*newTri)[j].neigh[2] = i + 1;
                break;
            }
        }

        for (j = *newNtri-1; j >= 0; j--) {
            if        ((*newTri)[j].indices[2] == iv0 &&
                       (*newTri)[j].indices[1] == iv1   ) {
                (*newTri)[i].neigh[2] = j + 1;
                (*newTri)[j].neigh[0] = i + 1;
                break;
            } else if ((*newTri)[j].indices[0] == iv0 &&
                       (*newTri)[j].indices[2] == iv1   ) {
                (*newTri)[i].neigh[2] = j + 1;
                (*newTri)[j].neigh[1] = i + 1;
                break;
            } else if ((*newTri)[j].indices[1] == iv0 &&
                       (*newTri)[j].indices[0] == iv1   ) {
                (*newTri)[i].neigh[2] = j + 1;
                (*newTri)[j].neigh[2] = i + 1;
                break;
            }
        }
        DPRINT4("neighbors for Triangle %4d (%4d,%4d,%4d)",
                i, (*newTri)[i].neigh[0], (*newTri)[i].neigh[1], (*newTri)[i].neigh[2]);
    }

 cleanup:
    FREE(uvcont );
    FREE(cont   );
    FREE(tempTri);
    FREE(lup    );
    FREE(ihole  );
    FREE(lupid  );
    FREE(prev   );

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * insideLoop -- check if point is inside the Loop                              *
 *                                                                              *
 ********************************************************************************
 */
static int
insideLoop(double   uu,                      /* (in)   U-coord of point to check */
           double   vv,                      /* (in)   V-coord of point to check */
           int      ibeg,                    /* (in)   index of first Vertex in loop */
           prmUV    uv[],                    /* (in)   UV-coords of Vertices */
           int      next[])                  /* (in)   next Vertex in Loop (or -1) */
{
    int         inout = 0;                   /* (out)  returned value: */
                                             /*        1 if inside Loop */
                                             /*        0 otherwise */

    int         ivrt, im1, ncross;
    double      utmp;

//    ROUTINE(insideLoop);
//    DPRINT4("%s(uu=%f, vv=%f, ibeg=%d) {",
//            routine, uu, vv, ibeg);

    /* ----------------------------------------------------------------------- */

    /*
     * the following algorithm is from O"Rourke"s 'Computational Geometry in C'
     *
     * count the number of segments in the Loop that cross a horizontal
     *    (constant V) ray from (uu,vv) to the left
     */
    ncross = 0;
    im1    = ibeg;
    ivrt   = next[im1];

    while (1) {
        if ((uv[ivrt].v > vv && uv[im1 ].v <= vv) ||
            (uv[im1 ].v > vv && uv[ivrt].v <= vv)   ) {
            utmp = ( uv[ivrt].u * (uv[im1 ].v - vv)
                   - uv[im1 ].u * (uv[ivrt].v - vv))
                 / (              (uv[im1 ].v - vv)
                   -              (uv[ivrt].v - vv));
            if (utmp > uu) {
                ncross++;
            }
        }

        if (ivrt == ibeg) break;

        im1  = ivrt;
        ivrt = next[im1];
    }

    /*
     * if there are an odd number of crosses, then we are inside
     */
    inout = ncross % 2;

// cleanup:
//    DPRINT2("%s --> inout=%d}", routine, inout);
    return inout;
}



/*
 ********************************************************************************
 *                                                                              *
 * interp1D -- 1-D interpolation into a table                                   *
 *                                                                              *
 ********************************************************************************
 */
static double
interp1D(double   xx,                        /* (in)   independent var for interp */
         int      ntab,                      /* (in)   size of table */
         double   xtab[],                    /* (in)   table of independent variables */
         double   ytab[])                    /* (in)   table of   dependent variables */
{
    double      yy;                          /* (out)  interpolated value */

    int         ileft, irite, imid;
    double      frac;

//    ROUTINE(interp1D);
//    DPRINT3("%s(xx=%f, ntab=%d) {",
//            routine, xx, ntab);

    /* ----------------------------------------------------------------------- */

    /*
     * binary search for interval
     */
    ileft = 0;
    irite = ntab - 1;

    while (irite > ileft+1) {
        imid = (ileft + irite) / 2;
        if (xx < xtab[imid]) {
            irite = imid;
        } else {
            ileft = imid;
        }
    }

    /*
     * linear interpolation
     */
    frac = (xx - xtab[ileft]) / (xtab[irite] - xtab[ileft]);
    yy   = (1 - frac) * ytab[ileft] + frac * ytab[irite];

// cleanup:
//    DPRINT2("%s --> yy=%f}", routine, yy);
    return yy;
}



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * plotTris -- plot the UV coordinates (using GRAFIC)                           *
 *                                                                              *
 ********************************************************************************
 */
extern void
plotTris(int      ntri,                      /* (in)   number of Triangles */
         prmTri   tri[],                     /* (in)   array  of Triangles */
         int      nvrt,                      /* (in)   number of Vertices */
         prmUV    uv[],                      /* (in)   array  of Parameters */
         char     pltitl[])                  /* (in)   plot title */
{
    int         io_kbd = 5;
    int         io_scr = 6;
    int         indgr   = 1 + 2 + 4 + 16 + 64;

    ROUTINE(plotTris);
    DPRINT3("%s(ntri=%d, nvrt=%d) {",
            routine, ntri, nvrt);

    /* ----------------------------------------------------------------------- */

    grinit_(&io_kbd, &io_scr, "plotTris",
                       strlen("plotTris"));
    grctrl_(plotTrisImage, &indgr, pltitl,
            (void*)(&ntri), (void*)(tri), (void*)(&nvrt),
            (void*)(uv),    NULL, NULL,
            NULL, NULL, NULL, NULL, strlen(pltitl));

// cleanup:
    DPRINT1("%s}", routine);
}
#endif



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * plotTrisImage -- used by plotTris                                            *
 *                                                                              *
 ********************************************************************************
 */
static void
plotTrisImage(int   *ifunct,                 /* (in)   GRAFIC function indicator */
              void  *ntriP,                  /* (in)   number of Triangles */
              void  *triP,                   /* (in)   array  of Triangles */
              void  *nvrtP,                  /* (in)   number of Vertices */
              void  *uvP,                    /* (in)   array  of Parametric Coords */
              void  *a4,                     /* (in)   dummy GRAFIC argument */
              void  *a5,                     /* (in)   dummy GRAFIC argument */
              void  *a6,                     /* (in)   dummy GRAFIC argument */
              void  *a7,                     /* (in)   dummy GRAFIC argument */
              void  *a8,                     /* (in)   dummy GRAFIC argument */
              void  *a9,                     /* (in)   dummy GRAFIC argument */
              float *scale,                  /* (out)  array of scales */
              char  *text,                   /* (out)  help text */
              int   textlen)                 /* (in)   length of text */
{
    int    *ntri   = (int    *) ntriP;
    int    *nvrt   = (int    *) nvrtP;
    prmTri *tri    = (prmTri *) triP;
    prmUV  *uv     = (prmUV  *) uvP;

    int    itri, ivrt, iv0, iv1, iv2, nneg;
    double umin, umax, vmin, vmax;
    float  u[3], v[3], area;

    int    ithree  = 3;
    int    iblack  = GR_BLACK;
    int    ired    = GR_RED;
    int    igreen  = GR_GREEN;

//    ROUTINE(plotTrisImage);

    /* ----------------------------------------------------------------------- */

    /* ---------- return scales ----------*/
    if (*ifunct == 0) {
        umin = uv[0].u;
        umax = uv[0].u;
        vmin = uv[0].v;
        vmax = uv[0].v;

        for (ivrt = 1; ivrt < *nvrt; ivrt++) {
            umin = MIN(umin, uv[ivrt].u);
            umax = MAX(umax, uv[ivrt].u);
            vmin = MIN(vmin, uv[ivrt].v);
            vmax = MAX(vmax, uv[ivrt].v);
        }

        scale[0] = umin;
        scale[1] = umax;
        scale[2] = vmin;
        scale[3] = vmax;

        strncpy(text, "X", textlen-1);

    /* ---------- plot image ---------- */
    } else if (*ifunct == 1) {

        /*
         * outlines in green and black
         */
        for (itri = 0; itri < *ntri; itri++) {
            iv0 = tri[itri].indices[0] - 1;
            iv1 = tri[itri].indices[1] - 1;
            iv2 = tri[itri].indices[2] - 1;

            u[0] = uv[iv0].u;
            u[1] = uv[iv1].u;
            u[2] = uv[iv2].u;
            v[0] = uv[iv0].v;
            v[1] = uv[iv1].v;
            v[2] = uv[iv2].v;

            grmov2_(&u[0], &v[0]);

            if (tri[itri].neigh[2] > 0) {
                grcolr_(&igreen);
            } else {
                grcolr_(&iblack);
            }
            grdrw2_(&u[1], &v[1]);

            if (tri[itri].neigh[0] > 0) {
                grcolr_(&igreen);
            } else {
                grcolr_(&iblack);
            }
            grdrw2_(&u[2], &v[2]);

            if (tri[itri].neigh[1] > 0) {
                grcolr_(&igreen);
            } else {
                grcolr_(&iblack);
            }
            grdrw2_(&u[0], &v[0]);
        }

        /*
         * cells with negative areas in red
         */
        nneg = 0;
        for (itri = 0; itri < *ntri; itri++) {
            iv0 = tri[itri].indices[0] - 1;
            iv1 = tri[itri].indices[1] - 1;
            iv2 = tri[itri].indices[2] - 1;

            u[0] = uv[iv0].u;
            u[1] = uv[iv1].u;
            u[2] = uv[iv2].u;
            v[0] = uv[iv0].v;
            v[1] = uv[iv1].v;
            v[2] = uv[iv2].v;

            area = (u[1] - u[0]) * (v[2] - v[0])
                 - (v[1] - v[0]) * (u[2] - u[0]);

            if (area <= 0) {
                grfil2_(u, v, &ithree, &ired);
                nneg++;

                DPRINT2("itri=%d   area=%f", itri, area);
                DPRINT3("%5d %15.8f %15.8f", iv0, u[0], v[0]);
                DPRINT3("%5d %15.8f %15.8f", iv1, u[1], v[1]);
                DPRINT3("%5d %15.8f %15.8f", iv2, u[2], v[2]);
            }
        }

        if (nneg > 0) {
            printf("%d triangles with area <= 0\n", nneg);
        }

        grcolr_(&iblack);
    }
}
#endif



/*
 ********************************************************************************
 *                                                                              *
 * printSMF -- print a sparse matrix                                            *
 *                                                                              *
 ********************************************************************************
 */
extern int
printSMF(FILE     *fp,                       /* (in)   file pointer */
         double   asmf[],                    /* (in)   sparse-matrix data */
         int      ismf[],                    /* (in)   sparse-matrix indices */
         double   x[],                       /* (in)   initial  guess */
         double   rhs[])                     /* (in)   right-hand side */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        >0  number of iterations taken */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_ZEROPIVOT */

    int         i, k, n;

//    ROUTINE(printSMF);
//    DPRINT2("%s(ismf[0]=%d) {",
//            routine, ismf[0]);

    /* ----------------------------------------------------------------------- */

    /*
     * extract the matrix size from ismf
     */
    n = ismf[0] - 1;

    /*
     * print the matrix
     */
    for (i = 0; i < n; i++) {
        fprintf(fp, "%5d    %5d %10.5f %10.5f %10.5f\n",
                i, ismf[i], asmf[i], x[i], rhs[i]);
        for (k = ismf[i]; k < ismf[i+1]; k++) {
//$$$            if (asmf[k] != 0) {
                fprintf(fp, "   %5d %5d %10.5f\n", k, ismf[k], asmf[k]);
//$$$            }
        }
    }

// cleanup:
//    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * sparseBCG -- solve sparse matrix with biconjugate gradient                   *
 *                                                                              *
 ********************************************************************************
 */
extern int
sparseBCG(double   asmf[],                   /* (in)   sparse-matrix data */
          int      ismf[],                   /* (in)   sparse-matrix indices */
          double   x[],                      /* (in)   initial  guess */
                                             /* (out)  solution to A * x = rhs */
          double   rhs[])                    /* (in)   right-hand side */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        >0  number of iterations taken */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_ZEROPIVOT */

    double      *p  = NULL;
    double      *pp = NULL;
    double      *r  = NULL;
    double      *rr = NULL;
    double      *z  = NULL;
    double      *zz = NULL;

    double      err, rmin, rmax, check;
    double      ak,        akden;
    double      bk, bknum, bkden;
    int         i, j, k, n, ipass, iter;

    double      tol     = 1e-8;              /* convergence tolerance */
    int         maxiter = 10000;             /* maximum number iterations */
    int         maxpass = 10;                /* maximum number of passes */

    ROUTINE(sparseBCG);
    DPRINT2("%s(ismf[0]=%d) {",
            routine, ismf[0]);

    /* ----------------------------------------------------------------------- */

    /*
     * extract the matrix size from ismf
     */
    n = ismf[0] - 1;

    /*
     * if all rhs are zero, just return the trivial result
     */
    rmin = +HUGEQ;
    rmax = -HUGEQ;
    for (i = 0; i < n; i++) {
        rmin = MIN(rmin, rhs[i]);
        rmax = MAX(rmax, rhs[i]);
    }
    if (fabs(rmin) < tol && fabs(rmax) < tol) {
        for (i = 0; i < n; i++) {
            x[i] = 0;
        }
        goto cleanup;
    }

    /*
     * make sure that no diagonal elements are zero
     */
    for (i = 0; i < n; i++) {
        if (fabs(asmf[i]) < EPS20) {
            status = PRM_ZEROPIVOT;
            goto cleanup;
        }
    }

    /*
     * allocate the needed matrixes
     */
    MALLOC(p,  double, n);
    MALLOC(pp, double, n);
    MALLOC(r,  double, n);
    MALLOC(rr, double, n);
    MALLOC(z,  double, n);
    MALLOC(zz, double, n);

    /*
     * solve the matrix using biconjugate gradient technique.
     *    code below is adapted from "Numerical Recipes"
     *
     * multiple passes might be needed if we ever get a zero denominator
     */
    ipass = 0;
 new_pass:
    ipass++;

    if (ipass > maxpass) {
        status = PRM_NOTCONVERGED;
        goto cleanup;
    }

    /*
     * calculate the initial residual
     */
    for (i = 0; i < n; i++) {                /* A * x -> r */
        r[i] = asmf[i] * x[i];

        for (k = ismf[i]; k < ismf[i+1]; k++) {
            j     = ismf[k];
            r[i] += asmf[k] * x[j];
        }
    }

    for (j = 0; j < n; j++) {
        r[ j] = rhs[j] - r[j];
        rr[j] = r[j];
    }

    for (j = 0; j < n; j++) {
        z[j] = r[j] / asmf[j];
    }

    /*
     * jump out if we already have a solution
     */
    err = 0;
    for (j = 0; j < n; j++) {
        err += SQR(r[j]);
    }
    err = sqrt(err);
    if (err < tol) {
        goto cleanup;
    }

    /*
     * main iteration loop
     */
    bkden = 0;

    for (iter = 0; iter < maxiter; iter++) {
        DPRINT2("iter=%5d  err=%15.8e", iter, err);

        for (j = 0; j < n; j++) {
            zz[j] = rr[j] / asmf[j];
        }

        /*
         * calculate coefficient bk and direction vectors p and pp
         */
        bknum = 0;
        for (j = 0; j < n; j++) {
            bknum += z[j] * rr[j];
        }

        if (iter == 0) {
            for (j = 0; j < n; j++) {
                p[ j] = z[ j];
                pp[j] = zz[j];
            }
        } else {
            bk = bknum / bkden;
            for (j = 0; j < n; j++) {
                p[ j] = bk * p[ j] +z[ j];
                pp[j] = bk * pp[j] +zz[j];
            }
        }

        /*
         * calculate coefficient ak, new iterate x, and new residuals r and rr
         */
        bkden = bknum;

        if (fabs(bkden) < EPS20) {
            DPRINT2("restarting because bknum = %f (ipass=%d)", bknum, ipass);
            goto new_pass;
        }

        for (i = 0; i < n; i++) {       /* A * p -> z */
            z[i] = asmf[i] * p[i];

            for (k = ismf[i]; k < ismf[i+1]; k++) {
                j     = ismf[k];
                z[i] += asmf[k] * p[j];
            }
        }

        akden = 0;
        for (j = 0; j < n; j++) {
            akden += z[j] * pp[j];
        }

        if (fabs(akden) < EPS20) {
            DPRINT2("restarting because akden = %f (ipass=%d)", akden, ipass);
            goto new_pass;
        }

        for (i = 0; i < n; i++) {        /* AT * pp -> zz */
            zz[i] = asmf[i] * pp[i];
        }

        for (i = 0; i < n; i++) {
            for (k = ismf[i]; k < ismf[i+1]; k++) {
                j      = ismf[k];
                zz[j] += asmf[k] * pp[i];
            }
        }

        ak = bknum / akden;
        for (j = 0; j < n; j++) {
            x[ j] += ak * p[ j];
            r[ j] -= ak * z[ j];
            rr[j] -= ak * zz[j];
        }

        for (j = 0; j < n; j++) {
            z[j] = r[j] / asmf[j];
        }

        /*
         * solve a! * z = r and check stopping criterion
         */
        err = 0;
        for (j = 0; j < n; j++) {
            err += SQR(r[j]);
        }
        err = sqrt(err);

        if (err < tol) {
            goto cleanup;
        }
    }

    DPRINT1("exceeded maxiter=%d", maxiter);
    status = PRM_NOTCONVERGED;

 cleanup:
    if (status == PRM_ZEROPIVOT || status == PRM_NOTCONVERGED) {
        DPRINT1("Sparse matrix that caused status=%d", status);

        for (i = 0; i < ismf[0]-1; i++) {
            check = -rhs[i] + asmf[i] * x[i];
            for (j = ismf[i]; j < ismf[i+1]; j++) {
                check += asmf[j] * x[ismf[j]];
            }

            DPRINT5("%5d      %12.6f %12.6f %12.6f %12.4e",
                    i, asmf[i], x[i], rhs[i], check);

            for (j = ismf[i]; j < ismf[i+1]; j++) {
                DPRINT2("     %5d %12.6f", ismf[j], asmf[j]);
            }
        }
    }

    FREE(zz);
    FREE(z );
    FREE(rr);
    FREE(r );
    FREE(pp);
    FREE(p );

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_CreateU -- create a parameterization from a set of Vertices              *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_CreateU(int      nvrt,                   /* (in)   number of Vertices */
            double   u[],                    /* (out)  array  of Parameters */
            prmXYZ   xyz[],                  /* (in)   array  of Coordinates */
            double   tol,                    /* (in)   tolerance for periodic check */
            int      *periodic)              /* (out)  = 0  no preiodicity */
                                             /*        = 1  periodic in U */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    int         ivrt;

    ROUTINE(prm_CreateU);
    DPRINT2("%s(nvrt=%d) {",
            routine, nvrt);

    /* ----------------------------------------------------------------------- */

    u[0] = 0;

    for (ivrt = 1; ivrt < nvrt; ivrt++) {
        u[ivrt] = u[ivrt-1] + sqrt(  SQR(xyz[ivrt].x - xyz[ivrt-1].x)
                                   + SQR(xyz[ivrt].y - xyz[ivrt-1].y)
                                   + SQR(xyz[ivrt].z - xyz[ivrt-1].z));
    }

    if (fabs(xyz[0].x - xyz[nvrt-1].x) < tol &&
        fabs(xyz[0].y - xyz[nvrt-1].y) < tol &&
        fabs(xyz[0].z - xyz[nvrt-1].z) < tol   ) {
        *periodic = 1;
        printf("   INFO:: U-periodic point found (nvrt=%d)\n", nvrt);
    } else {
        *periodic = 0;
    }

// cleanup:
    DPRINT1("periodic=%d", *periodic);
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_CreateUV -- create a parameterization from a set of Triangles            *
 *                 and Vertices                                                 *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_CreateUV(int      UVtype,                /* (in)   type of transformation (see status) */
             int      ntri,                  /* (in)   number of Triangles */
             prmTri   tri[],                 /* (in)   array  of Triangles */
  /*@null@*/ prmUVF   uvf[],                 /* (in)   array  of UVF for Tris */
             int      nvrt,                  /* (in)   number of Vertices */
  /*@null@*/ int      ptype[],               /* (in)   array  of point types */
  /*@null@*/ int      pindx[],               /* (in)   array  of point indices */
             prmUV    uv[],                  /* (out)  array  of Parameters */
             prmXYZ   xyz[],                 /* (in)   array  of Coordinates */
             int      *periodic,             /* (out)  = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        = 2  periodic in V */
                                             /*        =-1  periodicity problem */
             int      *ppnts[])              /* (out)  pointer to indices of periodic points
                                                       NOTE: user must EG_free after use */
{
    int         status = PRM_NOGLOBALUV;     /* (out)  return status */
                                             /*        PRM_NOGLOBALUV */
                                             /*        PRM_OK_ONEFACE */
                                             /*        PRM_OK_SIMPLE */
                                             /*        PRM_OK_JONES */
                                             /*        PRM_OK_AXIAL */
                                             /*        PRM_OK_PLANAR */
                                             /*        PRM_OK_FLOATER */
                                             /*        PRM_OK_UNROLLING */
                                             /*        <0  error condition */

    double      amin, area, dU, dV, dUper, dVper;
//$$$    double      dot, dx1, dy1, dz1, dx2, dy2, dz2;
    int         i, found, ipass, itri, ivrt, jvrt, iv0, iv1, iv2;
    int         faces[5], nfaces, nneg, nper, n, swap;

    int         *points = NULL;

    int         maxpass = 10000;

    ROUTINE(prm_CreateUV);
    DPRINT4("%s(UVtype=%d, ntri=%d, nvrt=%d) {",
            routine, UVtype, ntri, nvrt);

    for (itri = 0; itri < ntri; itri++) {
        DPRINT8("tri %5d: %5d %5d %5d  %5d %5d %5d  %5d", itri,
                tri[itri].indices[0],
                tri[itri].indices[1],
                tri[itri].indices[2],
                tri[itri].neigh[0],
                tri[itri].neigh[1],
                tri[itri].neigh[2],
                tri[itri].own);
    }
    if (uvf != NULL) {
        for (itri = 0; itri < ntri; itri++) {
            DPRINT7("tri %5d: %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e", itri,
                    uvf[itri].u0,
                    uvf[itri].v0,
                    uvf[itri].u1,
                    uvf[itri].v1,
                    uvf[itri].u2,
                    uvf[itri].v2);
        }
    }
    if (ptype != NULL && pindx != NULL) {
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            DPRINT6("vrt %5d: %5d %5d  %15.8e %15.8e %15.8e", ivrt,
                    ptype[ivrt],
                    pindx[ivrt],
                    xyz[ivrt].x,
                    xyz[ivrt].y,
                    xyz[ivrt].z);
        }
    }

    /* ----------------------------------------------------------------------- */

    *ppnts = NULL;

    /*
     * count the number of Faces (up to 5) associated with the Triangles
     */
    if (uvf != NULL && (UVtype == 0 || UVtype == PRM_OK_ONEFACE)) {
        nfaces = 0;
        for (itri = 0; itri < ntri; itri++) {
            found = 0;
            for (i = 0; i < nfaces; i++) {
                if (tri[itri].own == faces[i]) {
                    found = 1;
                    break;
                }
            }

            if (found == 0 && nfaces < 5) {
                faces[nfaces++] = tri[itri].own;
            }
        }
        for (i = 0; i < nfaces; i++) {
            DPRINT2("faces[%d] = %d", i, faces[i]);
        }

        /*
         * Surface with one Face (but only if all positive areas)
         */
        if (nfaces == 1) {
            status = PRM_OK_ONEFACE;
            DPRINT1("Setting UVtype = %d", status);

            for (itri = 0; itri < ntri; itri++) {
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                uv[iv0].u = uvf[itri].u0;
                uv[iv0].v = uvf[itri].v0;
                uv[iv1].u = uvf[itri].u1;
                uv[iv1].v = uvf[itri].v1;
                uv[iv2].u = uvf[itri].u2;
                uv[iv2].v = uvf[itri].v2;
            }

            /*
             * find the Triangle with the minimum area (in UV)
             */
            amin = HUGEQ;
            nneg = 0;

            for (itri = 0; itri < ntri; itri++) {
                iv0 = tri[itri].indices[0] - 1;
                iv1 = tri[itri].indices[1] - 1;
                iv2 = tri[itri].indices[2] - 1;

                area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                     - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);

                if (area < 0) nneg++;
                amin = MIN(amin, area);
            }

            /*
             * if more than half the Triangles have negative areas, try swapping
             *    u and v and then recomputing the smallest area
             */
            if (nneg > ntri/2) {
                for (ivrt = 0; ivrt < nvrt; ivrt++) {
                    DSWAP(uv[ivrt].u, uv[ivrt].v);
                }

                amin = HUGEQ;

                for (itri = 0; itri < ntri; itri++) {
                    iv0 = tri[itri].indices[0] - 1;
                    iv1 = tri[itri].indices[1] - 1;
                    iv2 = tri[itri].indices[2] - 1;

                    area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
                         - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);

                    amin = MIN(amin, area);
                }
            }

            /*
             * smooth if we have non-positive areas (indicating a bad tessellation)
             */
            for (ipass = 0; ipass <= maxpass; ipass++) {
                if (amin > EPS12) break;

                smoothUVcoords(ntri, tri, nvrt, uv, &amin);
            }

            /*
             * if we still have small or negative areas, reset UVtype so
             *    hopefully one of the techniques below will succeed
             */
            if (amin <= EPS12) {
                status = PRM_NOGLOBALUV;
                DPRINT1("Resetting UVtype = %d", status);
            }
        }
    } else {
        DPRINT0("uvf not given");
/*      nfaces = 9999;     */
    }

    /*
     * simple projection
     */
    if (status == PRM_NOGLOBALUV && (UVtype == 0 || UVtype == PRM_OK_SIMPLE)) {
        status = simpleProjection(ntri, tri, nvrt, xyz, uv);
        if (status != PRM_NOGLOBALUV) {
            CHECK_STATUS;
        }
    }

retry:
    /*
     * jones transformation
     */
    if (status == PRM_NOGLOBALUV && uvf != NULL && (UVtype == 0 || UVtype == PRM_OK_JONES)) {
        status = jonesTransformation(ntri, tri, uvf, nvrt, uv);
        if (status != PRM_NOGLOBALUV) {
            CHECK_STATUS;
        }
    }

    /*
     * axial projection
     */
    if (status == PRM_NOGLOBALUV && (UVtype == 0 || UVtype == PRM_OK_AXIAL)) {
        status = axialProjection(ntri, tri, nvrt, xyz, uv);
        if (status != PRM_NOGLOBALUV) {
            CHECK_STATUS;
        }
    }

    /*
     * planar projection
     */
    if (status == PRM_NOGLOBALUV && (UVtype == 0 || UVtype == PRM_OK_PLANAR)) {
        status = planarProjection(ntri, tri, nvrt, xyz, uv);
        if (status != PRM_NOGLOBALUV) {
            CHECK_STATUS;
        }
    }

    /*
     * apply floater's parameterization
     */
    if (status == PRM_NOGLOBALUV && (UVtype == 0 || UVtype == PRM_OK_FLOATER)) {
        status = floaterParameterization(ntri, tri, nvrt, xyz, uv);
        if (status != PRM_NOGLOBALUV) {
            CHECK_STATUS;
        }
    }

    /*
     * unrolling
     */
    if (status == PRM_NOGLOBALUV && (UVtype == 0 || UVtype == PRM_OK_UNROLLING)) {
        status = unroll(ntri, tri, nvrt, xyz, uv);
        if (status != PRM_NOGLOBALUV) {
            CHECK_STATUS;
        }
    }

    /*
     * check for periodicity (and correct if needed)
     */
    dUper = -1;
    dVper = -1;
    nper  = 0;
    if (ptype != NULL && pindx != NULL) {
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (ptype[ivrt] >= 0) {
                for (jvrt = ivrt+1; jvrt < nvrt; jvrt++) {
                    if (ptype[ivrt] == ptype[jvrt] &&
                        pindx[ivrt] == pindx[jvrt]    ) {
                        dU = fabs(uv[ivrt].u - uv[jvrt].u);
                        dV = fabs(uv[ivrt].v - uv[jvrt].v);

                        if (dUper < 0) {
                            dUper = dU;
                            dVper = dV;
                            nper++;
                        } else if (fabs(dU - dUper) < EPS04 &&
                                   fabs(dV - dVper) < EPS04   ) {
                            nper++;
                        } else {
                            printf("WARNING:: correcting periodic problem for Vertices %d and %d\n",
                                   ivrt, jvrt);
                            printf("          ptype=%5d  pindx=%5d\n", ptype[ivrt], pindx[ivrt]);
                            printf("          dUold=%10.5f  dVold=%10.5f\n", dU,    dV   );
                            printf("          dUper=%10.5f  dVper=%10.5f\n", dUper, dVper);

                            if (uv[ivrt].u < uv[jvrt].u) {
                                uv[ivrt].u -= (dUper - dU) / 2;
                                uv[jvrt].u += (dUper - dU) / 2;
                            } else {
                                uv[ivrt].u += (dUper - dU) / 2;
                                uv[jvrt].u -= (dUper - dU) / 2;
                            }

                            if (uv[ivrt].v < uv[jvrt].v) {
                                uv[ivrt].v -= (dVper - dV) / 2;
                                uv[jvrt].v += (dVper - dV) / 2;
                            } else {
                                uv[ivrt].v += (dVper - dV) / 2;
                                uv[jvrt].v -= (dVper - dV) / 2;
                            }

                            printf("          dUnew=%10.5f  dVnew=%10.5f\n",
                                   fabs(uv[ivrt].u - uv[jvrt].u),
                                   fabs(uv[ivrt].v - uv[jvrt].v));

                            nper++;
                        }
                    }
                }
            }
        }
    }

    if (nper > 0) {
        MALLOC(points, int, nper+1);
        *ppnts = points;

        if        (dUper >= EPS04 && dVper >= EPS04) {
            printf("   INFO:: %d UV-periodic points found (dUper=%10.5f dVper=%10.5f)\n",
                   nper, dUper, dVper);

            *periodic = 0;
            points[0] = 0;
        } else if (dUper >= EPS04) {
            printf("   INFO:: %d U-periodic points found (dUper=%10.5f dVper=%10.5f)\n",
                   nper, dUper, dVper);

            *periodic = 1;
            points[0] = nper;
        } else if (dVper >= EPS04) {
            printf("   INFO:: %d V-periodic points found (dUper=%10.5f dVper=%10.5f)\n",
                   nper, dUper, dVper);

            *periodic = 2;
            points[0] = nper;
        } else if (status != PRM_OK_PLANAR) {
            printf("   WARNING:: %d BAD periodic points found (dUper=%10.5f dVper=%10.5f)\n",
                   nper, dUper, dVper);

            *periodic = -1;
            points[0] = 0;
        } else {
            status = PRM_NOGLOBALUV;
            goto retry;
        }

/*      n = 0;       */
        if (*periodic == 1 && ptype != NULL && pindx != NULL) {
            n = 1;
            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                if (ptype[ivrt] >= 0) {
                    for (jvrt = ivrt+1; jvrt < nvrt; jvrt++) {
                        if (ptype[ivrt] == ptype[jvrt] &&
                            pindx[ivrt] == pindx[jvrt]   ) {
                            if (uv[ivrt].u < uv[jvrt].u) {
                                points[n] = ivrt;
                            } else {
                                points[n] = jvrt;
                            }
                            n++;
                        }
                    }
                }
            }
            if (n != nper+1) {
                DPRINT2("n=%d, nper+1=%d", n, nper+1);
                printf("ERROR:: prm_CreateUV: n=%d, nper+1=%d\n", n, nper+1);
                status = PRM_INTERNAL;
                goto cleanup;
            }

            for (ivrt = 1; ivrt < n-1; ivrt++) {
                for (jvrt = ivrt+1; jvrt < n; jvrt++) {
                    if (uv[points[ivrt]].v > uv[points[jvrt]].v) {
                        swap         = points[ivrt];
                        points[ivrt] = points[jvrt];
                        points[jvrt] = swap;
                    }
                }
            }
        } else if (*periodic == 2 && ptype != NULL && pindx != NULL) {
            n = 1;
            for (ivrt = 0; ivrt < nvrt; ivrt++) {
                if (ptype[ivrt] >= 0) {
                    for (jvrt = ivrt+1; jvrt < nvrt; jvrt++) {
                        if (ptype[ivrt] == ptype[jvrt] &&
                            pindx[ivrt] == pindx[jvrt]   ) {
                            if (uv[ivrt].v < uv[jvrt].v) {
                                points[n] = ivrt;
                            } else {
                                points[n] = jvrt;
                            }
                            n++;
                        }
                    }
                }
            }
            if (n != nper+1) {
                DPRINT2("n=%d, nper+1=%d", n, nper+1);
                printf("ERROR:: prm_CreateUV: n=%d, nper+1=%d\n", n, nper+1);
                status = PRM_INTERNAL;
                goto cleanup;
            }

            for (ivrt = 1; ivrt < n-1; ivrt++) {
                for (jvrt = ivrt+1; jvrt < n; jvrt++) {
                    if (uv[points[ivrt]].u > uv[points[jvrt]].u) {
                        swap         = points[ivrt];
                        points[ivrt] = points[jvrt];
                        points[jvrt] = swap;
                    }
                }
            }
        }
    } else {
        *periodic = 0;
    }

//$$$
//$$$    /*
//$$$     * if periodic line is not straight (in XYZ), then the fitter will
//$$$     *    fail to converge (since linear interpolation is used along the
//$$$     *    periodic seam).  therefore for now, just suppress the periodicity
//$$$     */
//$$$    if (*periodic > 0 && points != NULL) {
//$$$        dx1 = xyz[points[2]].x - xyz[points[1]].x;
//$$$        dy1 = xyz[points[2]].y - xyz[points[1]].y;
//$$$        dz1 = xyz[points[2]].z - xyz[points[1]].z;
//$$$
//$$$        for (i = 3; i <= nper; i++) {
//$$$            dx2 = xyz[points[i]].x - xyz[points[1]].x;
//$$$            dy2 = xyz[points[i]].y - xyz[points[1]].y;
//$$$            dz2 = xyz[points[i]].z - xyz[points[1]].z;
//$$$
//$$$            dot =     (dx1 * dx2 + dy1 * dy2 + dz1 * dz2)
//$$$                / sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1)
//$$$                / sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
//$$$            if (dot < 0.999) {
//$$$                printf("   WARNING:: periodicity suppressed because periodic line is not straight\n");
//$$$
//$$$                FREE(*ppnts);
//$$$                *periodic = 0;
//$$$                *ppnts    = NULL;
//$$$
//$$$                break;
//$$$            }
//$$$        }
//$$$    }
//$$$

 cleanup:
    if        (status == PRM_OK_ONEFACE) {
        PPRINT0("   UV via single Face");
    } else if (status == PRM_OK_SIMPLE) {
        PPRINT0("   UV via simple projection");
    } else if (status == PRM_OK_JONES) {
        PPRINT0("   UV via jones transformation");
    } else if (status == PRM_OK_AXIAL) {
        PPRINT0("   UV via axial projection");
    } else if (status == PRM_OK_FLOATER) {
        PPRINT0("   UV via floater parameterization");
    } else if (status == PRM_OK_PLANAR) {
        PPRINT0("   UV via planar projection");
    } else if (status == PRM_OK_UNROLLING) {
        PPRINT0("   UV via unrolling");
    }

    if ((*periodic <= 0) || (points == NULL)) {
        DPRINT1("periodic=%d", *periodic);
    } else {
        DPRINT2("periodic=%d,  points[0]=%d",
                *periodic, points[0]);
    }
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_SmoothU -- reparameterize based upon arc-lengths                         *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_SmoothU(int      periodic,               /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
            int      nvrt,                   /* (in)   number of Vertices */
            int      nvar,                   /* (in)   number of dependent vars */
            double   u[],                    /* (both) array  of parameters */
            double   xyz[])                  /* (in)   array  of Coordinates */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */

    ROUTINE(prm_SmoothU);
    DPRINT4("%s(periodic=%d, nvrt=%d, nvar=%d) {",
            routine, periodic, nvrt, nvar);

    /* ----------------------------------------------------------------------- */

    printf("prm_SmoothU is not implemented yet\n");

    printf("periodic=%d, nvrt=%d, nvar=%d, u[0]=%f, xyz[0]=%f\n",
           periodic, nvrt, nvar, u[0], xyz[0]);

//cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_SmoothUV -- reparameterize based upon arc-lengths and angles             *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_SmoothUV(int      type,                  /* (in)   smoothing type =1 for boundary only
                                                                      =2 for interior only
                                                                      =3 for both */
             int      periodic,              /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        = 2  periodic in V */
  /*@null@*/ int      ppnts[],               /* (in)   indices of periodic points */
             int      ntri,                  /* (in)   number of Triangles */
             prmTri   tri[],                 /* (in)   array  of Triangles */
             int      nvrt,                  /* (in)   number of Vertices */
             int      nvar,                  /* (in)   number of dependent vars */
             prmUV    uv[],                  /* (both) array  of Parameters */
             double   xyz[])                 /* (in)   array  of Coordinates */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_NEGATIVEAREAS */

    int         *prev   = NULL;              /* prev Vertex if on a loop, -1 otherwise */
    int         *next   = NULL;              /* next Vertex if on a loop, -1 otherwise */
    int         *corn   = NULL;              /* =1 if corner Vertex */

    int         *itab   = NULL;              /* table of Vertices    around original boundary */
    double      *ttab   = NULL;              /* table of arc-lengths around original boundary */
    double      *stab   = NULL;              /* table of uv-lengths  around original boundary */
    double      *utab   = NULL;              /* table of U           around original boundary */
    double      *vtab   = NULL;              /* table of U           around original boundary */

    sparseMat   amat;                        /* full matrix */
    double      *unew   = NULL;              /* temp storage for U-parameters */
    double      *vnew   = NULL;              /* temp storage for V-parameters */
    double      *urhs   = NULL;              /* RHS for U-parameters */
    double      *vrhs   = NULL;              /* RHS for V-parameters */
    double      *asmf   = NULL;              /* sparse matrix data */
    int         *ismf   = NULL;              /* sparse matrix indices */
    prmTri      *newTri = NULL;              /* array of Tris with holes filled */
    int         newNtri = 0;                 /* number of Tris with holes filled */

    int         ivrt, itri, im1, ip1, itmp;
    int         i, j, k, nn, nneg, ibeg, iend, imin, ntab;
    int         iv0, iv1, iv2, iv3, iv4, iv5, ivar;
//$$$    double      umin, umax, vmin, vmax;
    double      elxyzm, elxyzp, eluvm, cosxyz, cosmin, frac, snew, du, dv, duvmax;
    double      d01sq, d12sq, d20sq, d23sq, d30sq, d04sq, d41sq, d15sq, d52sq;
    double      ang0, ang1, ang2, ang3, ang4, ang5, sum, area;

    double      cos20 = 0.93969;             /* cos(20 deg) */

    ROUTINE(prm_SmoothUV);
    DPRINT6("%s(type=%d, periodic=%d, ntri=%d, nvrt=%d, nvar=%d) {",
            routine, type, periodic, ntri, nvrt, nvar);

    /* ----------------------------------------------------------------------- */

    amat.ni = 0;
    amat.is = NULL;
    if (periodic > 0 && ppnts == NULL) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

//$$$    /*
//$$$     * find the extrema of the UV
//$$$     */
//$$$    umin = uv[0].u;
//$$$    umax = uv[0].u;
//$$$    vmin = uv[0].v;
//$$$    vmax = uv[0].v;
//$$$
//$$$    for (ivrt = 0; ivrt < nvrt; ivrt++) {
//$$$        umin = MIN(umin, uv[ivrt].u);
//$$$        umax = MAX(umax, uv[ivrt].u);
//$$$        vmin = MIN(vmin, uv[ivrt].v);
//$$$        vmax = MAX(vmax, uv[ivrt].v);
//$$$    }
//$$$
//$$$    /*
//$$$     * normalize UV
//$$$     */
//$$$    for (ivrt = 0; ivrt < nvrt; ivrt++) {
//$$$        uv[ivrt].u = (uv[ivrt].u - umin) / (umax - umin);
//$$$        uv[ivrt].v = (uv[ivrt].v - vmin) / (vmax - vmin);
//$$$    }

#ifdef GRAFIC
    plotTris(ntri, tri, nvrt, uv,
             "~U~V~original triangulation (normalized)");
#endif

    /*
     * fill in any holes in the Triangulation
     */
    status = fillHoles(ntri, tri, nvrt, uv, &newNtri, &newTri);
    CHECK_STATUS;
    if (newTri == NULL) {
        status = EGADS_MALLOC;
        goto cleanup;
    }

#ifdef GRAFIC
    plotTris(newNtri, newTri, nvrt, uv,
             "~U~V~filled-in triangulation");
#endif

    status = initSmat(nvrt, &amat);
    CHECK_STATUS;
  
    /*
     * malloc the necessary arrays
     */
    MALLOC(prev,  int,    nvrt);
    MALLOC(next,  int,    nvrt);
    MALLOC(corn,  int,    nvrt);
    MALLOC(itab,  int,    nvrt+1);
    MALLOC(ttab,  double, nvrt+1);
    MALLOC(stab,  double, nvrt+1);
    MALLOC(utab,  double, nvrt+1);
    MALLOC(vtab,  double, nvrt+1);
    MALLOC(unew,  double, nvrt);
    MALLOC(vnew,  double, nvrt);
    MALLOC(urhs,  double, nvrt);
    MALLOC(vrhs,  double, nvrt);

    /*
     * determine the next and previous Vertex along a boundary (or -1 if interior)
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        unew[ivrt] = uv[ivrt].u;
        vnew[ivrt] = uv[ivrt].v;
        next[ivrt] = -1;
        prev[ivrt] = -1;
        corn[ivrt] = -1;
    }

    for (itri = 0; itri < newNtri; itri++) {
        iv0 = newTri[itri].indices[0] - 1;
        iv1 = newTri[itri].indices[1] - 1;
        iv2 = newTri[itri].indices[2] - 1;

        if (newTri[itri].neigh[0] <= 0) {
            prev[ iv2] = iv1;
            next[ iv1] = iv2;
            corn[ iv2] = 0;
            corn[ iv1] = 0;
        }
        if (newTri[itri].neigh[1] <= 0) {
            prev[ iv0] = iv2;
            next[ iv2] = iv0;
            corn[ iv0] = 0;
            corn[ iv2] = 0;
        }
        if (newTri[itri].neigh[2] <= 0) {
            prev[ iv1] = iv0;
            next[ iv0] = iv1;
            corn[ iv1] = 0;
            corn[ iv0] = 0;
        }
    }

    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        if (next[ivrt] >= 0) {
            DPRINT4("prev[%5d]=%5d   next[%5d]=%5d", ivrt, prev[ivrt], ivrt, next[ivrt]);
        }
    }

    /*
     *   ***** boundary smoothing *****
     */
    if ((type == 1 || type == 3) && (periodic == 0)) {
        DPRINT0("Boundary smoothing");

        /*
         * starting at the first boundary Vertex, set up splines for the
         *    initial UV, arclength, and cosang
         */
        ibeg = 0;
        while (next[ibeg] < 0) {
            ibeg++;
        }
        DPRINT1("ibeg=%d", ibeg);

        /*
         * march around the boundary and find the Vertex with the smallest cosxyz
         */
        cosmin = HUGEQ;
        imin   = 99999;
        ivrt   = ibeg;
        do {
            im1 = prev[ivrt];
            ip1 = next[ivrt];

            elxyzm = 0;
            elxyzp = 0;
            for (ivar = 0; ivar < nvar; ivar++) {
                elxyzm += SQR(xyz[nvar*ivrt+ivar] - xyz[nvar*im1+ivar]);
                elxyzp += SQR(xyz[nvar*ivrt+ivar] - xyz[nvar*ip1+ivar]);
            }
            elxyzm = sqrt(elxyzm);
            elxyzp = sqrt(elxyzp);

            cosxyz = 0;
            for (ivar = 0; ivar < nvar; ivar++) {
                cosxyz += (xyz[nvar*ivrt+ivar] - xyz[nvar*im1 +ivar])
                        * (xyz[nvar*ip1 +ivar] - xyz[nvar*ivrt+ivar]);
            }
            cosxyz /= (elxyzm * elxyzp);

            if (cosxyz < cos20) {
                corn[ivrt] = 1;
            }

            if (cosxyz < cosmin) {
                cosmin = cosxyz;
                imin   = ivrt;
            }

            ivrt = ip1;
        } while (ivrt != ibeg);

        ibeg = imin;
        DPRINT1("ibeg=%d", ibeg);

        /*
         * starting at the new ibeg, build tables around the boundary
         */
        ntab = 0;

        ivrt = ibeg;
        do {
            im1 = prev[ivrt];

            elxyzm = 0;
            for (ivar = 0; ivar < nvar; ivar++) {
                elxyzm += SQR(xyz[nvar*ivrt+ivar] - xyz[nvar*im1+ivar]);
            }
            elxyzm = sqrt(elxyzm);

            eluvm  = sqrt(  SQR(uv[ ivrt].u - uv[ im1].u)
                          + SQR(uv[ ivrt].v - uv[ im1].v));

            itab[ntab] = ivrt;
            if (ntab == 0) {
                ttab[ntab] = 0;
                stab[ntab] = 0;
            } else {
                ttab[ntab] = ttab[ntab-1] + elxyzm;
                stab[ntab] = stab[ntab-1] + eluvm;
            }
            utab[ntab] = uv[ivrt].u;
            vtab[ntab] = uv[ivrt].v;
            ntab++;

            ivrt = next[ivrt];
        } while (ivrt != ibeg);

        ivrt   = ibeg;
        im1    = prev[ivrt];

        elxyzm = 0;
        for (ivar = 0; ivar < nvar; ivar++) {
            elxyzm += SQR(xyz[nvar*ivrt+ivar] - xyz[nvar*im1+ivar]);
        }
        elxyzm = sqrt(elxyzm);

        eluvm  = sqrt(  SQR(uv[ ivrt].u - uv[ im1].u)
                      + SQR(uv[ ivrt].v - uv[ im1].v));

        itab[ntab] = itab[0];
        ttab[ntab] = ttab[ntab-1] + elxyzm;
        stab[ntab] = stab[ntab-1] + eluvm;
        utab[ntab] = uv[ivrt].u;
        vtab[ntab] = uv[ivrt].v;
        ntab++;

        corn[itab[0     ]] = 1;
        corn[itab[ntab-1]] = 1;

        for (i = 0; i < ntab; i++) {
            DPRINT7("tab %5d: %5d %5d %10.5f %10.5f %10.5f %10.5f",
                    i, itab[i], corn[itab[i]], ttab[i], stab[i], utab[i], vtab[i]);
        }

#ifdef GRAFIC
        {
            int    io_kbd = 5;
            int    io_scr = 6;
            int    indgr  = 1 + 4 + 16 + 64;
            int    ilin[2], isym[2], nper[2], nline;
            float  xplot[1000], yplot[1000];

            grinit_(&io_kbd, &io_scr, "prm_SmoothUV",
                               strlen("prm_SmoothUV"));

            ilin[0] = 1; isym[0] = 1; nper[0] = ntab; nline = 1;
            for (i = 0; i < ntab; i++) {
                xplot[i] = utab[i];
                yplot[i] = vtab[i];
            }
            grline_(ilin, isym, &nline,                "~u~v~boundary spline",
                    &indgr, xplot, yplot, nper, strlen("~u~v~boundary spline"));

            ilin[0] = 1; isym[0] = -1; nper[0] = ntab; nline = 2;
            ilin[1] = 2; isym[1] = -2; nper[1] = ntab;
            for (i = 0; i < ntab; i++) {
                xplot[i     ] = ttab[i];
                yplot[i     ] = utab[i];
                xplot[i+ntab] = ttab[i];
                yplot[i+ntab] = vtab[i];
            }
            grline_(ilin, isym, &nline,                "~t~u,v~boundary spline",
                    &indgr, xplot, yplot, nper, strlen("~t~u,v~boundary spline"));

            ilin[0] = -1; isym[0] = +1; nper[0] = ntab; nline = 2;
            ilin[1] = +2; isym[1] = -2; nper[1] = 2;
            for (i = 0; i < ntab; i++) {
                xplot[i] = ttab[i];
                yplot[i] = stab[i];
            }
            xplot[ntab  ] = ttab[0     ];
            yplot[ntab  ] = stab[0     ];
            xplot[ntab+1] = ttab[ntab-1];
            yplot[ntab+1] = stab[ntab-1];

            grline_(ilin, isym, &nline,                "~t~s~boundary spline",
                    &indgr, xplot, yplot, nper, strlen("~t~s~boundary spline"));
        }
#endif

        /*
         * change boundary Vertex locations so that spacings in UV mimic those in XYZ
         */
        ibeg = 0;
        while (ibeg < ntab-1) {
            iend = ibeg + 1;

            while(corn[itab[iend]] != 1) {
                iend++;
            }
            DPRINT2("ibeg=%d  iend=%d", ibeg, iend);

            for (i = ibeg+1; i < iend; i++) {
                ivrt = itab[i];
                frac = (ttab[i] - ttab[ibeg]) / (ttab[iend] - ttab[ibeg]);
                snew = (1 - frac) * stab[ibeg] + (frac) * stab[iend];

                uv[ivrt].u = interp1D(snew, ntab, stab, utab);
                uv[ivrt].v = interp1D(snew, ntab, stab, vtab);
            }

            ibeg = iend;
        }

#ifdef GRAFIC
        plotTris(ntri, tri, nvrt, uv,
                 "~U~V~after boundary smoothing");
#endif
    }

    /*
     *   ***** interior smoothing *****
     */
    if (type == 2 || type == 3) {
        DPRINT0("Interior smoothing");

        /*
         * initialize the amat matrix (which will be used to store influence info)
         */
        for (j = 0; j < nvrt; j++) {
            status = setSmat(&amat, j, j, 1.0);
            CHECK_STATUS;
            unew[j] = uv[j].u;
            vnew[j] = uv[j].v;
            urhs[j] = uv[j].u;
            vrhs[j] = uv[j].v;
        }

        /*
         * find the mean value weights at each Vertex due to its neighbors.  the
         *    weigt function is the "mean value coordinates" function as described
         *    by M. Floater
         *
         *                   iv3--------iv2---------iv5
         *                     \         /\a5       /
         *                      \       /a2\       /
         *                       \     /    \     /
         *                        \a3 / itri \   /
         *                         \ /a0    a1\ /
         *                         iv0--------iv1
         *                           \     a4 /
         *                            \      /
         *                             \    /
         *                              \  /
         *                               \/
         *                              iv4
         */
        for (itri = 0; itri < newNtri; itri++) {
            iv0 = newTri[itri].indices[0] - 1;
            iv1 = newTri[itri].indices[1] - 1;
            iv2 = newTri[itri].indices[2] - 1;

            d01sq = 0;
            d12sq = 0;
            d20sq = 0;
            for (ivar = 0; ivar < nvar; ivar++) {
                d01sq += SQR(xyz[nvar*iv0+ivar] - xyz[nvar*iv1+ivar]);
                d12sq += SQR(xyz[nvar*iv1+ivar] - xyz[nvar*iv2+ivar]);
                d20sq += SQR(xyz[nvar*iv2+ivar] - xyz[nvar*iv0+ivar]);
            }

            ang0 = ACOS((d20sq + d01sq - d12sq) / 2 / sqrt(d20sq * d01sq));
            ang1 = ACOS((d01sq + d12sq - d20sq) / 2 / sqrt(d01sq * d12sq));
            ang2 = ACOS((d12sq + d20sq - d01sq) / 2 / sqrt(d12sq * d20sq));

            if (newTri[itri].neigh[1] > 0) {
                itmp = newTri[itri].neigh[1] -  1;
                if (       newTri[itmp].neigh[0]-1 == itri) {
                    iv3 =  newTri[itmp].indices[0]-1;
                } else if (newTri[itmp].neigh[1]-1 == itri) {
                    iv3 =  newTri[itmp].indices[1]-1;
                } else {
                    iv3 =  newTri[itmp].indices[2]-1;
                }

                d23sq = 0;
                d30sq = 0;
                for (ivar = 0; ivar < nvar; ivar++) {
                    d23sq += SQR(xyz[nvar*iv2+ivar] - xyz[nvar*iv3+ivar]);
                    d30sq += SQR(xyz[nvar*iv3+ivar] - xyz[nvar*iv0+ivar]);
                }

                ang3 = ACOS((d30sq + d20sq - d23sq) / 2 / sqrt(d30sq * d20sq));

                status = setSmat(&amat, iv0, iv2,
                                 (tan(ang0/2) + tan(ang3/2)) / sqrt(d20sq));
                CHECK_STATUS;
            }

            if (newTri[itri].neigh[2] > 0) {
                itmp = newTri[itri].neigh[2] -  1;
                if (       newTri[itmp].neigh[0]-1 == itri) {
                    iv4 =  newTri[itmp].indices[0]-1;
                } else if (newTri[itmp].neigh[1]-1 == itri) {
                    iv4 =  newTri[itmp].indices[1]-1;
                } else {
                    iv4 =  newTri[itmp].indices[2]-1;
                }

                d04sq = 0;
                d41sq = 0;
                for (ivar = 0; ivar < nvar; ivar++) {
                    d04sq += SQR(xyz[nvar*iv0+ivar] - xyz[nvar*iv4+ivar]);
                    d41sq += SQR(xyz[nvar*iv4+ivar] - xyz[nvar*iv1+ivar]);
                }

                ang4 = ACOS((d41sq + d01sq - d04sq) / 2 / sqrt(d41sq * d01sq));

                status = setSmat(&amat, iv1, iv0,
                                 (tan(ang1/2) + tan(ang4/2)) / sqrt(d01sq));
                CHECK_STATUS;
            }

            if (newTri[itri].neigh[0] > 0) {
                itmp = newTri[itri].neigh[0] -  1;
                if (       newTri[itmp].neigh[0]-1 == itri) {
                    iv5 =  newTri[itmp].indices[0]-1;
                } else if (newTri[itmp].neigh[1]-1 == itri) {
                    iv5 =  newTri[itmp].indices[1]-1;
                } else {
                    iv5 =  newTri[itmp].indices[2]-1;
                }

                d15sq = 0;
                d52sq = 0;
                for (ivar = 0; ivar < nvar; ivar++) {
                    d15sq += SQR(xyz[nvar*iv1+ivar] - xyz[nvar*iv5+ivar]);
                    d52sq += SQR(xyz[nvar*iv5+ivar] - xyz[nvar*iv2+ivar]);
                }

                ang5 = ACOS((d52sq + d12sq - d15sq) / 2 / sqrt(d52sq * d12sq));

                status = setSmat(&amat, iv2, iv1,
                                 (tan(ang2/2) + tan(ang5/2)) / sqrt(d12sq));
                CHECK_STATUS;
            }
        }

        /*
         * set up the final amat and the right-hand sides
         */
        for (i = 0; i < nvrt; i++) {
            if (next[i] < 0) {
                status = sumSmat(&amat, i,    &sum);
                CHECK_STATUS;
                status = setSmat(&amat, i, i, -sum);
                CHECK_STATUS;
                urhs[i] = 0;
                vrhs[i] = 0;
            } else {
                status = diagSmat(&amat, i);
                CHECK_STATUS;
                status = setSmat(&amat, i, i, 1.0);
                CHECK_STATUS;
                urhs[i] = uv[i].u;
                vrhs[i] = uv[i].v;
            }
        }

        /*
         * count the number of non-zero entries in amat (so that we know how
         *    big the sparse-matrix forms should be)
         */
        nn = countSmat(&amat);

        /*
         * allocate matrices for sparse-matrix form.  this form is as
         *    described in 'Numerical Recipes'
         */
        MALLOC(asmf, double, nn+1);
        MALLOC(ismf, int,    nn+1);

        /*
         * store the matrix in sparse-matrix form
         */
      
        for (j = 0; j < nvrt; j++) {
            status = getSmat(&amat, j, j, &asmf[j]);
            CHECK_STATUS;
        }

        ismf[0] = nvrt + 1;

        k = nvrt;
        for (i = 0; i < nvrt; i++) {
            status = fillSmat(&amat, i, &k, asmf, ismf);
            CHECK_STATUS;
            ismf[i+1] = k + 1;
        }

        /*
         * solve the matrix equations using biconjugate-gradient technique
         */
        DPRINT0("U smoothing");
        status = sparseBCG(asmf, ismf, unew, urhs);
        CHECK_STATUS;

        DPRINT0("V smoothing");
        status = sparseBCG(asmf, ismf, vnew, vrhs);
        CHECK_STATUS;

        DPRINT0("Old, new, and change in  Vertex locations");
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            DPRINT8("%5d%5d  %10.5f%10.5f%10.5f  %10.5f%10.5f%10.5f",
                    ivrt, next[ivrt],
                    uv[ivrt].u, unew[ivrt], uv[ivrt].u-unew[ivrt],
                    uv[ivrt].v, vnew[ivrt], uv[ivrt].v-vnew[ivrt]);
        }

        /*
         * store the solution into uv for all the interior Vertices
         */
        duvmax = 0;
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (next[ivrt] < 0) {
                du = fabs(uv[ivrt].u - unew[ivrt]);
                dv = fabs(uv[ivrt].v - vnew[ivrt]);

                uv[ivrt].u = unew[ivrt];
                uv[ivrt].v = vnew[ivrt];

                duvmax = MAX3(duvmax, du, dv);
            }
        }
        DPRINT1("duvmax = %10.5f", duvmax);

        FREE(ismf);
        FREE(asmf);

#ifdef GRAFIC
        plotTris(ntri, tri, nvrt, uv,
                 "~U~V~after interior smoothing");
#endif
    }

//$$$    /*
//$$$     * un-normalize UV
//$$$     */
//$$$    for (ivrt = 0; ivrt < nvrt; ivrt++) {
//$$$        uv[ivrt].u = umin + uv[ivrt].u * (umax - umin);
//$$$        uv[ivrt].v = vmin + uv[ivrt].v * (vmax - vmin);
//$$$    }

    /*
     * find the number of Triangles with negative areas (in UV)
     */
    nneg = 0;
    for (itri = 0; itri < ntri; itri++) {
        iv0 = tri[itri].indices[0] - 1;
        iv1 = tri[itri].indices[1] - 1;
        iv2 = tri[itri].indices[2] - 1;

        area = (uv[iv1].u - uv[iv0].u) * (uv[iv2].v - uv[iv0].v)
             - (uv[iv1].v - uv[iv0].v) * (uv[iv2].u - uv[iv0].u);

        if (area < 0) nneg++;
    }
    DPRINT2("smoothing has nneg=%d of %d", nneg, ntri);

    if (nneg > 0) {
        status = PRM_NEGATIVEAREAS;
    } else {
        status = EGADS_SUCCESS;
    }

cleanup:
    freeSmat(&amat);
    FREE(newTri);
    FREE(ismf  );
    FREE(asmf  );
    FREE(vrhs  );
    FREE(urhs  );
    FREE(vnew  );
    FREE(unew  );
    FREE(vtab  );
    FREE(utab  );
    FREE(stab  );
    FREE(ttab  );
    FREE(itab  );
    FREE(corn  );
    FREE(next  );
    FREE(prev  );

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_NormalizeU -- normalize a parameterization                               *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_NormalizeU(double   halo,                /* (in)   halo size */
               int      periodic,            /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
               int      nvrt,                /* (in)   number of Vertices */
               double   u[])                 /* (in)   array  of Parameters */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      umin = +HUGEQ, umax = -HUGEQ;
    double      haloU, au, bu;
    int         i;

    ROUTINE(prm_NormalizeU);
    DPRINT4("%s(halo=%f, periodic=%d, nvrt=%d) {",
            routine, halo, periodic, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * find the extrema
     */
    for (i = 0; i < nvrt; i++) {
        umin = MIN(umin, u[i]);
        umax = MAX(umax, u[i]);
    }

    if (periodic == 1) {
        haloU = 0;
    } else {
        haloU = halo;
    }

    /*
     * set up slope and intercept of transformation
     */

    au = (haloU * (umax + umin) - umin) / (umax - umin);
    bu = (1 - 2 * haloU)                / (umax - umin);

    /*
     * apply the linear transformation
     */
    for (i = 0; i < nvrt; i++) {
        u[i] = au + bu * u[i];
    }

// cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_NormalizeUV -- normalize a parameterization                              *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_NormalizeUV(double   halo,               /* (in)   halo size */
                int      periodic,           /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        = 2  periodic in V */
                int      nvrt,               /* (in)   number of Vertices */
                prmUV    uv[])               /* (both) array  of Parameters */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      umin = +HUGEQ, umax = -HUGEQ, vmin = +HUGEQ, vmax = -HUGEQ;
    double      haloU, haloV, au, av, bu, bv;
    int         i;

    ROUTINE(prm_NormalizeUV);
    DPRINT4("%s(halo=%f, periodic=%d, nvrt=%d) {",
            routine, halo, periodic, nvrt);

    /* ----------------------------------------------------------------------- */

    /*
     * find the extrema
     */
    for (i = 0; i < nvrt; i++) {
        umin = MIN(umin, uv[i].u);
        umax = MAX(umax, uv[i].u);
        vmin = MIN(vmin, uv[i].v);
        vmax = MAX(vmax, uv[i].v);
    }

    if        (periodic == 1) {
        haloU = 0;
        haloV = halo;
    } else if (periodic == 2) {
        haloU = halo;
        haloV = 0;
    } else {
        haloU = halo;
        haloV = halo;
    }

    DPRINT2("haloU=%10.5f  haloV=%10.5f", haloU, haloV);

    /*
     * set up slope and intercept of transformation
     */

    au = (haloU * (umax + umin) - umin) / (umax - umin);
    av = (haloV * (vmax + vmin) - vmin) / (vmax - vmin);
    bu = (1 - 2 * haloU)                / (umax - umin);
    bv = (1 - 2 * haloV)                / (vmax - vmin);

    /*
     * apply the linear transformation
     */
    umin = +HUGEQ;
    umax = -HUGEQ;
    vmin = +HUGEQ;
    vmax = -HUGEQ;

    for (i = 0; i < nvrt; i++) {
        uv[i].u = au + bu * uv[i].u;
        uv[i].v = av + bv * uv[i].v;

        umin = MIN(umin, uv[i].u);
        umax = MAX(umax, uv[i].u);
        vmin = MIN(vmin, uv[i].v);
        vmax = MAX(vmax, uv[i].v);
    }

    DPRINT4("umin=%10.5f  umax=%10.5f  vmin=%10.5f  vmax=%10.5f",
            umin, umax, vmin, vmax);

// cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}
