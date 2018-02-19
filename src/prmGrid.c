/*
 ******************************************************************************
 * EGADS: Electronic Geometry Aircraft Design System                          *
 *                                                                            *
 * This module creates and manages grid fits                                  *
 * Written by John Dannenhoffer @ Syracuse University                         *
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

#include "prm.h"
#include "egadsTypes.h"
#include "egadsInternals.h"

#ifdef GRAFIC
#include "grafic.h"
#endif


/*
 * internal routines defined below
 */
static double computeDotmin           (gridTree*, int, prmUV[], double[], int, prmTri[]);
static int    copyGrid                (gridTree*, gridTree*);
static int    divideCell              (gridTree*, int, int);
static void   evalBicubic             (double, double, int, double[], double[],
                                       double[], double[], double[],
                                       /*@null@*/double[], /*@null@*/double[],
                                       /*@null@*/double[], /*@null@*/double[],
                                       /*@null@*/double[]);
#ifdef GRAFIC
static double evalToLevel             (gridTree*, int, int, prmUV);
#endif
static int    finestCell              (gridTree*, double, double);
static int    globalRefine2d          (int, gridTree*, int, int, /*@null@*/int[],
                                       prmUV[], double[], double*, double*);
static int    initialTree2d           (gridTree*, int, int, int, /*@null@*/int[],
                                       prmUV[], double[], double*, double*);
#ifdef GRAFIC
static void   plotGrid                (int, int, int, int, prmUV[], double[], double[], double[]);
static void   plotGridImage           (int*, void*, void*, void*, void*, void*,
                                       void*, void*, void*, void*, void*,
                                       float*, char*, int);
#endif
#ifdef GRAFIC
static void   plotLevels              (gridTree*);
static void   plotLevelsImage         (int*, void*, void*, void*, void*, void*,
                                       void*, void*, void*, void*, void*,
                                       float*, char*, int);
#endif
#ifdef GRAFIC2
static void   plotVrts                (int, int, prmUV[], double[]);
static void   plotVrtsImage           (int*, void*, void*, void*, void*, void*,
                                       void*, void*, void*, void*, void*,
                                       float*, char*, int);
#endif
#ifdef DEBUG
static void   printTree               (FILE*, gridTree*);
#endif
static int    splineWithGaps2d        (int, prmUV[], double[], int, /*@null@*/int[], int, int,
                                       double[], double[], double[], double[]);

/*
 * external routines defined in prmUV
 */
extern int    printSMF                (FILE*, double[], int[], double[], double[]);
extern int    sparseBCG               (double[], int[], double[], double[]);

/*
 * global constants
 */
static int    globalSizeLimit = 257;

#define  HUGEQ           1.0e+40
#define  EPS12           1.0e-12
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
//#undef   DEBUG
#ifdef   DEBUG
   #define DOPEN \
           {if (dbg_fp == NULL) dbg_fp = fopen("prmGrid.dbg", "w");}
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
   #define DPRINT3x(FORMAT,A,B,C) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C);}
   #define DPRINT4(FORMAT,A,B,C,D) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D); DFLUSH;}
   #define DPRINT4x(FORMAT,A,B,C,D) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D);}
   #define DPRINT5(FORMAT,A,B,C,D,E) \
           {DOPEN; fprintf(dbg_fp, FORMAT, A, B, C, D, E);  DFLUSH;}
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
   #define DPRINT3x(FORMAT,A,B,C)
   #define DPRINT4(FORMAT,A,B,C,D)
   #define DPRINT4x(FORMAT,A,B,C,D)
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
   #define GI_OUT0(A) \
           DPRINT0(A) \
           printf( A) \
           printf( "\n")
   #define GI_OUT1(A,B) \
           DPRINT1(A,B) \
           printf( A,B) \
           printf( "\n")
   #define GI_OUT2(A,B,C) \
           DPRINT2(A,B,C) \
           printf( A,B,C) \
           printf( "\n")
#else
   #define GI_OUT0(A)
   #define GI_OUT1(A,B)
   #define GI_OUT2(A,B,C)
#endif

/*
 * macros for error checking
 */
#define CHECK_STATUS \
        if (status < EGADS_SUCCESS) {\
            GI_OUT2("ERROR:: BAD STATUS=%d in routine %s", status, routine);\
            goto cleanup;\
        }
#define MALLOC(PTR,TYPE,SIZE) \
        DPRINT3("mallocing %s in routine %s (size=%d)", #PTR, routine, SIZE);\
        PTR = (TYPE *) EG_alloc((SIZE) * sizeof(TYPE)); \
        if (PTR == NULL) {\
            GI_OUT2("ERROR:: MALLOC PROBLEM for %s in routine %s", #PTR, routine);\
            status = EGADS_MALLOC;\
            goto cleanup;\
        }
#define RALLOC(PTR,TYPE,SIZE) \
        DPRINT3("rallocing %s in routine %s (size=%d)", #PTR, routine, SIZE);\
        realloc_temp = EG_reall(PTR, (SIZE) * sizeof(TYPE)); \
        if (PTR == NULL) {\
            GI_OUT2("ERROR:: RALLOC PROBLEM for %s in routine %s", #PTR, routine);\
            status = EGADS_MALLOC;\
            goto cleanup;\
        } else {\
            PTR = (TYPE *)realloc_temp;\
        }
#define FREE(PTR) \
        DPRINT2("freeing %s in routine %s", #PTR, routine);\
        EG_free(PTR);\
        PTR = NULL;
static void *realloc_temp = NULL;            /* used by RALLOC macro */
#define MEMCPY(DEST,SRC,SIZE) \
        if ((SIZE) > 0) {\
            memcpy(DEST,SRC,SIZE);\
        }



/*
 ********************************************************************************
 *                                                                              *
 * computeDotmin -- compute minimum dot product between surface normals         *
 *                                                                              *
 ********************************************************************************
 */
static double
computeDotmin(gridTree *tree,                /* (in)   pointer to Grid Tree */
              int      nvar,                 /* (in)   number of dependent vars */
              prmUV    uv[],                 /* (in)   array  of Parameters */
              double   var[],                /* (in)   array  of dependent vars */
              int      ntri,                 /* (in)   number of Triangles (optional) */
              prmTri   tri[])                /* (in)   array  of Triangles (optional) */
{
    double      dotmin = 1;                  /* (out)  minimum dot product */

    int         itri, ip0, ip1, ip2;
    double      dot, norm1, norm1x, norm1y, norm1z, norm2, norm2x, norm2y, norm2z;
    double      xyz[3], dxyzdu[3], dxyzdv[3];
    prmUV       uuvv;

    ROUTINE(computeDotmin);
    DPRINT3("%s(nvar=%d, ntri=%d) {",
            routine, nvar, ntri);

    /* ----------------------------------------------------------------------- */

    /*
     * loop through the Triangles and compare its surface normal to the
     *    normal obtained by evaluating the Grid
     */
    for (itri = 0; itri < ntri; itri++) {

        /*
         * discrete normal of Triangle
         */
        ip0 = tri[itri].indices[0] - 1;
        ip1 = tri[itri].indices[1] - 1;
        ip2 = tri[itri].indices[2] - 1;

        norm1x = (var[nvar*ip1+1] - var[nvar*ip0+1]) * (var[nvar*ip2+2] - var[nvar*ip0+2])
               - (var[nvar*ip1+2] - var[nvar*ip0+2]) * (var[nvar*ip2+1] - var[nvar*ip0+1]);
        norm1y = (var[nvar*ip1+2] - var[nvar*ip0+2]) * (var[nvar*ip2+0] - var[nvar*ip0+0])
               - (var[nvar*ip1+0] - var[nvar*ip0+0]) * (var[nvar*ip2+2] - var[nvar*ip0+2]);
        norm1z = (var[nvar*ip1+0] - var[nvar*ip0+0]) * (var[nvar*ip2+1] - var[nvar*ip0+1])
               - (var[nvar*ip1+1] - var[nvar*ip0+1]) * (var[nvar*ip2+0] - var[nvar*ip0+0]);

        norm1  = sqrt(norm1x * norm1x + norm1y * norm1y + norm1z * norm1z);

        norm1x /= norm1;
        norm1y /= norm1;
        norm1z /= norm1;

        /*
         * fitted normal at centroid of Triangle
         */
        uuvv.u = (uv[ip0].u + uv[ip1].u + uv[ip2].u) / 3;
        uuvv.v = (uv[ip0].v + uv[ip1].v + uv[ip2].v) / 3;

        prm_EvalGrid(*tree, uuvv, xyz, dxyzdu, dxyzdv, NULL, NULL, NULL);

        norm2x = dxyzdu[1] * dxyzdv[2] - dxyzdu[2] * dxyzdv[1];
        norm2y = dxyzdu[2] * dxyzdv[0] - dxyzdu[0] * dxyzdv[2];
        norm2z = dxyzdu[0] * dxyzdv[1] - dxyzdu[1] * dxyzdv[0];

        norm2  = sqrt(norm2x * norm2x + norm2y * norm2y + norm2z * norm2z);

        norm2x /= norm2;
        norm2y /= norm2;
        norm2z /= norm2;

        dot = (norm1x * norm2x + norm1y * norm2y + norm1z * norm2z);
        DPRINT8("itri=%5d  norm1=(%6.3f,%6.3f,%6.3f)  norm2=(%6.3f,%6.3f,%6.3f)  dot=%10.5f",
                itri, norm1x, norm1y, norm1z, norm2x, norm2y, norm2z, dot);

        if (dot < dotmin) {
            dotmin = dot;
        }
    }

// cleanup:
    DPRINT2("%s --> dotmin=%f}", routine, dotmin);
    return dotmin;
}



/*
 ********************************************************************************
 *                                                                              *
 * copyGrid -- copy one Grid Tree to another                                    *
 *                                                                              *
 ********************************************************************************
 */
static int
copyGrid(gridTree *dest,                     /* (out)  pointer to destination Tree */
         gridTree *src)                      /* (in)   pointer to source Tree */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */

    ROUTINE(copyGrid);
    DPRINT1("%s() {",
            routine);

    /* ----------------------------------------------------------------------- */

    /*
     * free any arrays in the desitnation
     */
    FREE(dest->cell);
    FREE(dest->knot);

    /*
     * copy the scalar data
     */
#ifndef __clang_analyzer__
    dest->nvar = src->nvar;
    dest->nu   = src->nu;
    dest->nv   = src->nv;
    dest->ncel = src->ncel;
    dest->nknt = src->nknt;
#else
    dest->nvar = dest->nu = dest->nv = dest->ncel = dest->nknt = 1;
#endif

    /*
     * allocate new arrays
     */
    MALLOC(dest->cell, gridCell,   dest->ncel                );
    MALLOC(dest->knot, double,    (dest->nknt)*4*(dest->nvar));

    if (src->cell != NULL)
    MEMCPY(dest->cell, src->cell, (dest->ncel)               *sizeof(gridCell));
    if (src->knot != NULL)
    MEMCPY(dest->knot, src->knot, (dest->nknt)*4*(dest->nvar)*sizeof(double  ));

 cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * divideCell -- divide the given Cell                                          *
 *                                                                              *
 ********************************************************************************
 */
static int
divideCell(gridTree *tree,                   /* (in)   pointer to Grid Tree */
           int      icel,                    /* (in)   Cell to divide */
           int      dtype)                   /* (in)   division type */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_BADPARAM */
                                             /*        PRM_BADDIVISION */

    int         iknt, nknt, ncel, ivar, jcel;
    int         isw, ise, inw, ine;
    int         knot02, knot12, knot22;
    int         knot01, knot11, knot21;
    int         knot00, knot10, knot20;
    int         nborws, nborwn, nbores, nboren;
    int         nborsw, nborse, nbornw, nborne;
    double      umin, umax, vmin, vmax;

    int         nvar4 = 4 * tree->nvar;

    ROUTINE(divideCell);
    DPRINT3("%s(icel=%d, dtype=%d) {",
            routine, icel, dtype);

    /* ----------------------------------------------------------------------- */

    /*
     * make sure that Cell is divisible
     */
    if (icel < 0 || icel >= tree->ncel) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (tree->cell[icel].dtype > 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (dtype < 1 || dtype > 3) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

    ncel = tree->ncel;
    nknt = tree->nknt;

    /*
     * make maps of Knot numbers and neighbor Cells
     *
     *              nbornw     nborne
     *          0,2 ------ 1,2 ------ 2,2
     *           |          |          |
     *    nborwn |          |          | nboren
     *           |          |          |
     *          0,1 ------ 1,1 ------ 2,1
     *           |          |          |
     *    nborws |          |          | nbores
     *           |          |          |
     *          0,0 ------ 1,0 ------ 2,0
     *              nborsw     nborse
     */

    knot00 = -1;
    knot10 = -1;
    knot20 = -1;
    knot01 = -1;
    knot11 = -1;
    knot21 = -1;
    knot02 = -1;
    knot12 = -1;
    knot22 = -1;

    jcel = tree->cell[icel].wCell;
    if (jcel < 0) {
        nborws = jcel;
        nborwn = jcel;
    } else if (tree->cell[jcel].dtype == 0) {
        nborws = 0;
        nborwn = 0;
    } else if (tree->cell[jcel].dtype == 1) {
        nborws = tree->cell[jcel].child + 1;
        nborwn = tree->cell[jcel].child + 1;

        if (nborws > 0) {
            knot00 = tree->cell[nborws].seKnot;
            knot02 = tree->cell[nborwn].neKnot;
        }
    } else if (tree->cell[jcel].dtype == 2) {
        nborws = tree->cell[jcel].child;
        nborwn = tree->cell[jcel].child + 1;

        if (nborws > 0) {
            knot00 = tree->cell[nborws].seKnot;
            knot01 = tree->cell[nborws].neKnot;
            knot02 = tree->cell[nborwn].neKnot;
        }
    } else {
        nborws = tree->cell[jcel].child + 1;
        nborwn = tree->cell[jcel].child + 3;

        if (nborws > 0) {
            knot00 = tree->cell[nborws].seKnot;
            knot01 = tree->cell[nborws].neKnot;
            knot02 = tree->cell[nborwn].neKnot;
        }
    }

    jcel = tree->cell[icel].eCell;
    if (jcel < 0) {
        nbores = jcel;
        nboren = jcel;
    } else if (tree->cell[jcel].dtype == 0) {
        nbores = 0;
        nboren = 0;
    } else if (tree->cell[jcel].dtype == 1) {
        nbores = tree->cell[jcel].child;
        nboren = tree->cell[jcel].child;

        if (nbores > 0) {
            knot20 = tree->cell[nbores].swKnot;
            knot22 = tree->cell[nboren].nwKnot;
        }
    } else if (tree->cell[jcel].dtype == 2) {
        nbores = tree->cell[jcel].child;
        nboren = tree->cell[jcel].child + 1;

        if (nbores > 0) {
            knot20 = tree->cell[nbores].swKnot;
            knot21 = tree->cell[nbores].nwKnot;
            knot22 = tree->cell[nboren].nwKnot;
        }
    } else {
        nbores = tree->cell[jcel].child;
        nboren = tree->cell[jcel].child + 2;

        if (nbores > 0) {
            knot20 = tree->cell[nbores].swKnot;
            knot21 = tree->cell[nbores].nwKnot;
            knot22 = tree->cell[nboren].nwKnot;
        }
    }

    jcel = tree->cell[icel].sCell;
    if (jcel < 0) {
        nborsw = jcel;
        nborse = jcel;
    } else if (tree->cell[jcel].dtype == 0) {
        nborsw = 0;
        nborse = 0;
    } else if (tree->cell[jcel].dtype == 1) {
        nborsw = tree->cell[jcel].child;
        nborse = tree->cell[jcel].child + 1;

        if (nborsw > 0) {
            knot00 = tree->cell[nborsw].nwKnot;
            knot10 = tree->cell[nborsw].neKnot;
            knot20 = tree->cell[nborse].neKnot;
        }
    } else if (tree->cell[jcel].dtype == 2) {
        nborsw = tree->cell[jcel].child + 1;
        nborse = tree->cell[jcel].child + 1;

        if (nborsw > 0) {
            knot00 = tree->cell[nborsw].nwKnot;
            knot20 = tree->cell[nborse].neKnot;
        }
    } else {
        nborsw = tree->cell[jcel].child + 2;
        nborse = tree->cell[jcel].child + 3;

        if (nborsw > 0) {
            knot00 = tree->cell[nborsw].nwKnot;
            knot10 = tree->cell[nborsw].neKnot;
            knot20 = tree->cell[nborse].neKnot;
        }
    }

    jcel = tree->cell[icel].nCell;
    if (jcel < 0) {
        nbornw = jcel;
        nborne = jcel;
    } else if (tree->cell[jcel].dtype == 0) {
        nbornw = 0;
        nborne = 0;
    } else if (tree->cell[jcel].dtype == 1) {
        nbornw = tree->cell[jcel].child;
        nborne = tree->cell[jcel].child + 1;

        if (nbornw > 0) {
            knot02 = tree->cell[nbornw].swKnot;
            knot12 = tree->cell[nbornw].seKnot;
            knot22 = tree->cell[nborne].seKnot;
        }
    } else if (tree->cell[jcel].dtype == 2) {
        nbornw = tree->cell[jcel].child;
        nborne = tree->cell[jcel].child;

        if (nbornw > 0) {
            knot02 = tree->cell[nbornw].swKnot;
            knot22 = tree->cell[nborne].seKnot;
        }
    } else {
        nbornw = tree->cell[jcel].child;
        nborne = tree->cell[jcel].child + 1;

        if (nbornw > 0) {
            knot02 = tree->cell[nbornw].swKnot;
            knot12 = tree->cell[nbornw].seKnot;
            knot22 = tree->cell[nborne].seKnot;
        }
    }

    /*
     * number the new Knots (if any)
     */
    if        (dtype == 1) {
        if (knot00 < 0) knot00 = nknt++;
        if (knot10 < 0) knot10 = nknt++;
        if (knot20 < 0) knot20 = nknt++;
        if (knot02 < 0) knot02 = nknt++;
        if (knot12 < 0) knot12 = nknt++;
        if (knot22 < 0) knot22 = nknt++;
    } else if (dtype == 2) {
        if (knot00 < 0) knot00 = nknt++;
        if (knot20 < 0) knot20 = nknt++;
        if (knot01 < 0) knot01 = nknt++;
        if (knot21 < 0) knot21 = nknt++;
        if (knot02 < 0) knot02 = nknt++;
        if (knot22 < 0) knot22 = nknt++;
    } else {
        if (knot00 < 0) knot00 = nknt++;
        if (knot10 < 0) knot10 = nknt++;
        if (knot20 < 0) knot20 = nknt++;
        if (knot01 < 0) knot01 = nknt++;
        if (knot11 < 0) knot11 = nknt++;
        if (knot21 < 0) knot21 = nknt++;
        if (knot02 < 0) knot02 = nknt++;
        if (knot12 < 0) knot12 = nknt++;
        if (knot22 < 0) knot22 = nknt++;
    }

    /*
     * make sure that dtype is consistent with neighbors
     */
    if        (dtype == 1) {
        if (nborsw > 0 && nborsw == nborse) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
        if (nbornw > 0 && nbornw == nborne) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
    } else if (dtype == 2) {
        if (nborws > 0 && nborws == nborwn) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
        if (nbores > 0 && nbores == nboren) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
    } else {
        if (nborsw > 0 && nborsw == nborse) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
        if (nbornw > 0 && nbornw == nborne) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
        if (nborws > 0 && nborws == nborwn) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
        if (nbores > 0 && nbores == nboren) {
            status = PRM_BADDIVISION;
            goto cleanup;
        }
    }

    /*
     * tentative names of children Cells
     */
    if (dtype == 1) {
        isw = ncel;
        ise = ncel + 1;
        inw = ncel;
        ine = ncel + 1;
    } else if (dtype == 2) {
        isw = ncel;
        ise = ncel;
        inw = ncel + 1;
        ine = ncel + 1;
    } else {
        isw = ncel;
        ise = ncel + 1;
        inw = ncel + 2;
        ine = ncel + 3;
    }

//$$$    DPRINT2("          %5d     %5d",       nbornw, nborne);
//$$$    DPRINT3("     %5d     %5d     %5d",    knot02, knot12, knot22);
//$$$    DPRINT4("%5d     %5d     %5d     %5d", nborwn, inw, ine, nboren);
//$$$    DPRINT3("     %5d     %5d     %5d",    knot01, knot11, knot21);
//$$$    DPRINT4("%5d     %5d     %5d     %5d", nborws, isw, ise, nbores);
//$$$    DPRINT3("     %5d     %5d     %5d",    knot00, knot10, knot20);
//$$$    DPRINT2("          %5d     %5d",       nborsw, nborse);

    /*
     * extend the Cell array for the new Cells
     */
    if        (dtype == 1) {
        RALLOC(tree->cell, gridCell, tree->ncel+2);
    } else if (dtype == 2) {
        RALLOC(tree->cell, gridCell, tree->ncel+2);
    } else {
        RALLOC(tree->cell, gridCell, tree->ncel+4);
    }

    umin = tree->cell[icel].umin;
    umax = tree->cell[icel].umax;
    vmin = tree->cell[icel].vmin;
    vmax = tree->cell[icel].vmax;

    /*
     * make the children Cell
     */
    if        (dtype == 1) {
        tree->cell[ncel].swKnot = knot00;
        tree->cell[ncel].seKnot = knot10;
        tree->cell[ncel].nwKnot = knot02;
        tree->cell[ncel].neKnot = knot12;
        tree->cell[ncel].wCell  = nborws;
        tree->cell[ncel].eCell  = ise;
        tree->cell[ncel].sCell  = nborsw;
        tree->cell[ncel].nCell  = nbornw;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   =  umin;
        tree->cell[ncel].umax   = (umin + umax) / 2;
        tree->cell[ncel].vmin   =  vmin;
        tree->cell[ncel].vmax   =         vmax;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;

        tree->cell[ncel].swKnot = knot10;
        tree->cell[ncel].seKnot = knot20;
        tree->cell[ncel].nwKnot = knot12;
        tree->cell[ncel].neKnot = knot22;
        tree->cell[ncel].wCell  = isw;
        tree->cell[ncel].eCell  = nbores;
        tree->cell[ncel].sCell  = nborse;
        tree->cell[ncel].nCell  = nborne;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   = (umin + umax) / 2;
        tree->cell[ncel].umax   =         umax;
        tree->cell[ncel].vmin   =  vmin;
        tree->cell[ncel].vmax   =         vmax;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;
    } else if (dtype == 2) {
        tree->cell[ncel].swKnot = knot00;
        tree->cell[ncel].seKnot = knot20;
        tree->cell[ncel].nwKnot = knot01;
        tree->cell[ncel].neKnot = knot21;
        tree->cell[ncel].wCell  = nborws;
        tree->cell[ncel].eCell  = nbores;
        tree->cell[ncel].sCell  = nborsw;
        tree->cell[ncel].nCell  = inw;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   =  umin;
        tree->cell[ncel].umax   =         umax;
        tree->cell[ncel].vmin   =  vmin;
        tree->cell[ncel].vmax   = (vmin + vmax) / 2;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;

        tree->cell[ncel].swKnot = knot01;
        tree->cell[ncel].seKnot = knot21;
        tree->cell[ncel].nwKnot = knot02;
        tree->cell[ncel].neKnot = knot22;
        tree->cell[ncel].wCell  = nborwn;
        tree->cell[ncel].eCell  = nboren;
        tree->cell[ncel].sCell  = isw;
        tree->cell[ncel].nCell  = nbornw;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   =  umin;
        tree->cell[ncel].umax   =         umax;
        tree->cell[ncel].vmin   = (vmin + vmax) / 2;
        tree->cell[ncel].vmax   =         vmax;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;
    } else {
        tree->cell[ncel].swKnot = knot00;
        tree->cell[ncel].seKnot = knot10;
        tree->cell[ncel].nwKnot = knot01;
        tree->cell[ncel].neKnot = knot11;
        tree->cell[ncel].wCell  = nborws;
        tree->cell[ncel].eCell  = ise;
        tree->cell[ncel].sCell  = nborsw;
        tree->cell[ncel].nCell  = inw;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   =  umin;
        tree->cell[ncel].umax   = (umin + umax) / 2;
        tree->cell[ncel].vmin   =  vmin;
        tree->cell[ncel].vmax   = (vmin + vmax) / 2;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;

        tree->cell[ncel].swKnot = knot10;
        tree->cell[ncel].seKnot = knot20;
        tree->cell[ncel].nwKnot = knot11;
        tree->cell[ncel].neKnot = knot21;
        tree->cell[ncel].wCell  = isw;
        tree->cell[ncel].eCell  = nbores;
        tree->cell[ncel].sCell  = nborse;
        tree->cell[ncel].nCell  = ine;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   = (umin + umax) / 2;
        tree->cell[ncel].umax   =         umax;
        tree->cell[ncel].vmin   =  vmin;
        tree->cell[ncel].vmax   = (vmin + vmax) / 2;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;

        tree->cell[ncel].swKnot = knot01;
        tree->cell[ncel].seKnot = knot11;
        tree->cell[ncel].nwKnot = knot02;
        tree->cell[ncel].neKnot = knot12;
        tree->cell[ncel].wCell  = nborwn;
        tree->cell[ncel].eCell  = ine;
        tree->cell[ncel].sCell  = isw;
        tree->cell[ncel].nCell  = nbornw;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   =  umin;
        tree->cell[ncel].umax   = (umin + umax) / 2;
        tree->cell[ncel].vmin   = (vmin + vmax) / 2;
        tree->cell[ncel].vmax   =         vmax;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;

        tree->cell[ncel].swKnot = knot11;
        tree->cell[ncel].seKnot = knot21;
        tree->cell[ncel].nwKnot = knot12;
        tree->cell[ncel].neKnot = knot22;
        tree->cell[ncel].wCell  = inw;
        tree->cell[ncel].eCell  = nboren;
        tree->cell[ncel].sCell  = ise;
        tree->cell[ncel].nCell  = nborne;
        tree->cell[ncel].dtype  = 0;
        tree->cell[ncel].child  = 0;
        tree->cell[ncel].count  = 0;
        tree->cell[ncel].umin   = (umin + umax) / 2;
        tree->cell[ncel].umax   =         umax;
        tree->cell[ncel].vmin   = (vmin + vmax) / 2;
        tree->cell[ncel].vmax   =         vmax;
        tree->cell[ncel].rmserr = 0;
        tree->cell[ncel].maxerr = 0;
        ncel++;
    }

    /*
     * inform the neighbors of the new Cells
     */
    if (nborws > 0) tree->cell[nborws].eCell = isw;
    if (nborwn > 0) tree->cell[nborwn].eCell = inw;
    if (nbores > 0) tree->cell[nbores].wCell = ise;
    if (nboren > 0) tree->cell[nboren].wCell = ine;
    if (nborsw > 0) tree->cell[nborsw].nCell = isw;
    if (nborse > 0) tree->cell[nborse].nCell = ise;
    if (nbornw > 0) tree->cell[nbornw].sCell = inw;
    if (nborne > 0) tree->cell[nborne].sCell = ine;

    /*
     * inform icel of its new children
     */
    tree->cell[icel].dtype = dtype;
    tree->cell[icel].child = isw;

    /*
     * extend the Knot array
     */
    RALLOC(tree->knot, double, nvar4*nknt);

    for (iknt = tree->nknt; iknt < nknt; iknt++) {
        for (ivar = 0; ivar < nvar4; ivar++) {
            tree->knot[nvar4*iknt+ivar] = 0;
        }
    }

    /*
     * save the number of Cells and Knots in the Tree
     */
    tree->ncel = ncel;
    tree->nknt = nknt;

 cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * evalBicubic -- evaluate a bicubic                                            *
 *                                                                              *
 ********************************************************************************
 */
static void
evalBicubic(double   ss,                     /* (in)   U-fractional distance */
            double   tt,                     /* (in)   V-fractional distance */
            int      nvar,                   /* (in)   number of variables */
            double   sw[],                   /* (in)   Knot on southwest */
            double   se[],                   /* (in)   Knot on southeast */
            double   nw[],                   /* (in)   Knot on northwest */
            double   ne[],                   /* (in)   Knot on northeast */
            double   xyz[],                  /* (out)  dependent variables */
  /*@null@*/double   du[],                   /* (out)  optional d/dU */
  /*@null@*/double   dv[],                   /* (out)  optional d/dV */
  /*@null@*/double   duu[],                  /* (out)  optional d2/dU2 */
  /*@null@*/double   duv[],                  /* (out)  optional d2/dUdV */
  /*@null@*/double   dvv[])                  /* (out)  optional d2/dV2 */
{
    double      s0,  s1,  s2,  s3;
    double      t0,  t1,  t2,  t3;
    double      a11, a12, a13, a14;
    double      a21, a22, a23, a24;
    double      a31, a32, a33, a34;
    double      a41, a42, a43, a44;
    double      s10, s20, s30, s40, temp = 0;

    int         ivar;
    int         nvar2 = 2 * nvar;
    int         nvar3 = 3 * nvar;

//    ROUTINE(evalBicubic);
//    DPRINT3("%s(ss=%f, tt=%f) {",
//            routine, ss, tt);

    /* ----------------------------------------------------------------------- */

    for (ivar = 0; ivar < nvar; ivar++) {
        s0 = -3 * sw[     ivar] + 3 * nw[     ivar] -2 * sw[nvar2+ivar] - nw[nvar2+ivar];
        s1 = -3 * se[     ivar] + 3 * ne[     ivar] -2 * se[nvar2+ivar] - ne[nvar2+ivar];
        s2 = -3 * sw[nvar+ivar] + 3 * nw[nvar+ivar] -2 * sw[nvar3+ivar] - nw[nvar3+ivar];
        s3 = -3 * se[nvar+ivar] + 3 * ne[nvar+ivar] -2 * se[nvar3+ivar] - ne[nvar3+ivar];

        t0 =  2 * sw[     ivar] - 2 * nw[     ivar] +    sw[nvar2+ivar] + nw[nvar2+ivar];
        t1 =  2 * se[     ivar] - 2 * ne[     ivar] +    se[nvar2+ivar] + ne[nvar2+ivar];
        t2 =  2 * sw[nvar+ivar] - 2 * nw[nvar+ivar] +    sw[nvar3+ivar] + nw[nvar3+ivar];
        t3 =  2 * se[nvar+ivar] - 2 * ne[nvar+ivar] +    se[nvar3+ivar] + ne[nvar3+ivar];

        a11 =      sw[      ivar];
        a12 =      sw[nvar2+ivar];
        a13 =      s0;
        a14 =      t0;

        a21 =      sw[nvar +ivar];
        a22 =      sw[nvar3+ivar];
        a23 =      s2;
        a24 =      t2;

        a31 = -3 * sw[      ivar] + 3 * se[      ivar] - 2 * sw[nvar +ivar] - se[nvar +ivar];
        a32 = -3 * sw[nvar2+ivar] + 3 * se[nvar2+ivar] - 2 * sw[nvar3+ivar] - se[nvar3+ivar];
        a33 = -3 * s0             + 3 * s1             - 2 * s2             - s3;
        a34 = -3 * t0             + 3 * t1             - 2 * t2             - t3;

        a41 =  2 * sw[      ivar] - 2 * se[      ivar] +     sw[nvar +ivar] + se[nvar +ivar];
        a42 =  2 * sw[nvar2+ivar] - 2 * se[nvar2+ivar] +     sw[nvar3+ivar] + se[nvar3+ivar];
        a43 =  2 * s0             - 2 * s1             +     s2             + s3;
        a44 =  2 * t0              -2 * t1             +     t2             + t3;

        s10 = a11 + tt * (a12 + tt * (a13 + tt * a14));
        s20 = a21 + tt * (a22 + tt * (a23 + tt * a24));
        s30 = a31 + tt * (a32 + tt * (a33 + tt * a34));
        s40 = a41 + tt * (a42 + tt * (a43 + tt * a44));

        xyz[ivar] = s10 + ss * (s20 + ss * (s30 + ss * s40));

        if (du != NULL && dv != NULL) {
            du[ivar] = s20 + ss * (2 * s30 + 3 * ss * s40);
            temp     =             2 * s30 + 6 * ss * s40;

            s10 = a12 + tt * (2 * a13 + 3 * tt * a14);
            s20 = a22 + tt * (2 * a23 + 3 * tt * a24);
            s30 = a32 + tt * (2 * a33 + 3 * tt * a34);
            s40 = a42 + tt * (2 * a43 + 3 * tt * a44);

            dv[ivar] = s10 + ss * (s20 + ss * (s30 + ss * s40));

            if (duu != NULL && duv != NULL && dvv != NULL) {
                duu[ivar] = temp;
                duv[ivar] = s20 + ss * (2 * s30 + 3 * ss * s40);

                s10 = 2 * a13 + 6 * tt * a14;
                s20 = 2 * a23 + 6 * tt * a24;
                s30 = 2 * a33 + 6 * tt * a34;
                s40 = 2 * a43 + 6 * tt * a44;

                dvv[ivar] = s10 + ss * (s20 + ss * (s30 + ss * s40));
            }
        }
    }

// cleanup:
//    DPRINT1("%s --> }", routine);
}



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * evalToLevel -- evaluate a Grid to a given level                              *
 *                                                                              *
 ********************************************************************************
 */
extern double
evalToLevel(gridTree *tree,                  /* (in)   Grid Tree */
            int      nlvl,                   /* (in)   maximum evaluation levels */
            int      ivar,                   /* (in)   variable index */
            prmUV    uv)                     /* (in)   independent variables */
{
    double      answer;                      /* (out)  result of evaluation */

    double      xyz[99];

    int         icel, ilvl;
    double      ss, tt, umin, umax, vmin, vmax;

    int         nvar  =      tree->nvar;
    int         nvar4 =  4 * tree->nvar;

//    ROUTINE(evalToLevel);
//    DPRINT3("%s(uu=%f, vv=%f) {",
//            routine, uv.u, uv.v);

    /* ----------------------------------------------------------------------- */

    /*
     * initialize the dependent variables
     */
    answer = 0;

    /*
     * start evaluation at the root of the Grid Tree
     */
    icel = 0;

    /*
     * recursively add the appropriate evaluations
     */
    for (ilvl = 0; ilvl < nlvl; ilvl++) {

        /*
         * add in current Cell
         */
        umin = tree->cell[icel].umin;
        umax = tree->cell[icel].umax;
        vmin = tree->cell[icel].vmin;
        vmax = tree->cell[icel].vmax;

        ss =  (uv.u - umin) / (umax - umin);
        tt =  (uv.v - vmin) / (vmax - vmin);

        evalBicubic(ss, tt, nvar,
                    &(tree->knot[nvar4*(tree->cell[icel].swKnot)]),
                    &(tree->knot[nvar4*(tree->cell[icel].seKnot)]),
                    &(tree->knot[nvar4*(tree->cell[icel].nwKnot)]),
                    &(tree->knot[nvar4*(tree->cell[icel].neKnot)]),
                    xyz, NULL, NULL, NULL, NULL, NULL);

        answer += xyz[ivar];

        /*
         * if at bottom of hierarchy, stop
         */
        if (tree->cell[icel].dtype == 0) {
            break;

        /*
         * otherwise go to the appropriate child
         */
        } else if (tree->cell[icel].dtype == 1) {
            if (uv.u <= (umin+umax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 1;
            }

        } else if (tree->cell[icel].dtype == 2) {
            if (uv.v <= (vmin+vmax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 1;
            }

        } else if (uv.u <= (umin+umax)/2) {
            if (uv.v <= (vmin+vmax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 2;
            }
        } else {
            if (uv.v <= (vmin+vmax)/2) {
                icel  = tree->cell[icel].child + 1;
            } else {
                icel  = tree->cell[icel].child + 3;
            }
        }
    }

// cleanup:
//    DPRINT2("%s --> answer=%f}", routine, answer);
    return answer;
}
#endif



/*
 ********************************************************************************
 *                                                                              *
 * finestCell -- find the finest Cell containing (u,v)                          *
 *                                                                              *
 ********************************************************************************
 */
static int
finestCell(gridTree *tree,                   /* (in)   Grid Tree */
           double   u,                       /* (in)   first  parameter */
           double   v)                       /* (in)   second parameter */
{
    int         icel = 0;                    /* (out)  finest Cell containing u,v */

    double      umin, umax, vmin, vmax;

//    ROUTINE(finestCell);
//    DPRINT3("%s(u=%f, v=%f) {",
//            routine, u, v);

    /* ----------------------------------------------------------------------- */

    /*
     * traverse the Tree until we find an undivided Cell
     */
    while (1) {
        umin = tree->cell[icel].umin;
        umax = tree->cell[icel].umax;
        vmin = tree->cell[icel].vmin;
        vmax = tree->cell[icel].vmax;

        if        (tree->cell[icel].dtype == 0) {
            break;
        } else if (tree->cell[icel].dtype == 1) {
            if (u <= (umin+umax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 1;
            }

        } else if (tree->cell[icel].dtype == 2) {
            if (v <= (vmin+vmax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 1;
            }

        } else if (u <= (umin+umax)/2) {
            if (v <= (vmin+vmax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 2;
            }
        } else {
            if (v <= (vmin+vmax)/2) {
                icel  = tree->cell[icel].child + 1;
            } else {
                icel  = tree->cell[icel].child + 3;
            }
        }
    }

//cleanup:
//    DPRINT2("%s -> icel=%d}", routine, icel);
    return icel;
}



/*
 ********************************************************************************
 *                                                                              *
 * globalRefine2d -- refine a Tree globally                                     *
 *                                                                              *
 ********************************************************************************
 */
extern int
globalRefine2d(int      dtype,               /* (in)   division type */
               gridTree *tree,               /* (both) Grid Tree to refine */
               int      nvrt,                /* (in)   number of Vertices */
               int      periodic,            /* (in)   periodicity flag */
    /*@null@*/ int      ppnts[],             /* (in)   indices of periodic points */
               prmUV    uv[],                /* (in)   array  of Vertices */
               double   resid[],             /* (both) array  of residuals */
               double   *rmserr,             /* (out)  RMS     error at Vertices */
               double   *maxerr)             /* (out)  maximum error at Vertices */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      *xvrt   = NULL;
    double      *xspln  = NULL;
    double      *uspln  = NULL;
    double      *vspln  = NULL;
    double      *cspln  = NULL;
    double      *xyz    = NULL;

    int         ibeg, iend, ii, jj, ivrt, ivar, icel, iknt, nuv;
    int         ncel;
    double      umin, umax, vmin, vmax, uu, vv, ss, tt, err;

    int         nvar  =     tree->nvar;
    int         nvar2 = 2 * tree->nvar;
    int         nvar3 = 3 * tree->nvar;
    int         nvar4 = 4 * tree->nvar;

    ROUTINE(globalRefine2d);
    DPRINT4("%s(dtype=%d, nvrt=%d, periodic=%d) {",
            routine, dtype, nvrt, periodic);

    /* ----------------------------------------------------------------------- */

    *rmserr = 0;
    *maxerr = 0;

    MALLOC(xyz, double, nvar);

    /*
     * find the first and last Cell with dtype==0
     */
    ibeg = +99999;
    iend = -99999;

    for (icel = 0; icel < tree->ncel; icel++) {
        if (tree->cell[icel].dtype == 0) {
            ibeg = MIN(ibeg, icel);
            iend = MAX(iend, icel);
        }
    }

    /*
     * adjust number of Knots in Tree
     */
    if (dtype == 1) {
        tree->nu = 2 * tree->nu - 1;
    } else if (dtype == 2) {
        tree->nv = 2 * tree->nv - 1;
    } else {
        tree->nu = 2 * tree->nu - 1;
        tree->nv = 2 * tree->nv - 1;
    }

    /* ncel           first new Cell */
    /* tree->ncel-1   last  new Cell */
    ncel = tree->ncel;

    /*
     * create the new Cells and new Knots
     */
    for (icel = ibeg; icel <= iend; icel++) {
        status = divideCell(tree, icel, dtype);
        CHECK_STATUS;
    }

    /*
     * count number of Vertices in each new Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        icel = finestCell(tree, uv[ivrt].u, uv[ivrt].v);
        tree->cell[icel].count++;
    }

#ifdef DEBUG
    printTree(dbg_fp, tree);
#endif

    /*
     * allocate temp storage for the Vertices and the extra data
     */
    nuv =             (tree->nu    ) * (tree->nv    );
    DPRINT2("nu=%d, nv=%d", tree->nu, tree->nv);

    MALLOC(xvrt,  double, nvrt);
    MALLOC(xspln, double, nuv );
    MALLOC(uspln, double, nuv );
    MALLOC(vspln, double, nuv );
    MALLOC(cspln, double, nuv );

    for (ivar = 0; ivar < nvar; ivar++) {
        DPRINT1("working on ivar=%d", ivar);

        /*
         * copy the Vertex data into xvrt
         */
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            xvrt[ivrt] = resid[nvar*ivrt+ivar];
        }

        /*
         * spline fit for this group
         */
        status = splineWithGaps2d(nvrt, uv, xvrt, periodic, ppnts,
                                  tree->nu, tree->nv, xspln, uspln, vspln, cspln);
        CHECK_STATUS;

        /*
         * store the data in the Knots
         */
        for (icel = ncel; icel < tree->ncel; icel++) {
            ii = tree->nu * tree->cell[icel].umin;
            jj = tree->nv * tree->cell[icel].vmin;

            iknt = tree->cell[icel].swKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii  )+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii  )+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii  )+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii  )+tree->nu*(jj  )];

            iknt = tree->cell[icel].seKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii+1)+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii+1)+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii+1)+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii+1)+tree->nu*(jj  )];

            iknt = tree->cell[icel].nwKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii  )+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii  )+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii  )+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii  )+tree->nu*(jj+1)];

            iknt = tree->cell[icel].neKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii+1)+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii+1)+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii+1)+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii+1)+tree->nu*(jj+1)];
        }
    }

    /*
     * subtract the latest spline from "resid" and keep
     *    track of the max error associated with each Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uu = uv[ivrt].u;
        vv = uv[ivrt].v;

        for (icel = ncel; icel < tree->ncel; icel++) {
            umin = tree->cell[icel].umin;
            umax = tree->cell[icel].umax;
            vmin = tree->cell[icel].vmin;
            vmax = tree->cell[icel].vmax;

            if (uu >= umin && uu <= umax && vv >= vmin && vv <= vmax) {
                ss  = (uu - umin) / (umax - umin);
                tt  = (vv - vmin) / (vmax - vmin);
                evalBicubic(ss, tt, nvar,
                            &(tree->knot[nvar4*(tree->cell[icel].swKnot)]),
                            &(tree->knot[nvar4*(tree->cell[icel].seKnot)]),
                            &(tree->knot[nvar4*(tree->cell[icel].nwKnot)]),
                            &(tree->knot[nvar4*(tree->cell[icel].neKnot)]),
                            xyz, NULL, NULL, NULL, NULL, NULL);

                for (ivar = 0; ivar < nvar; ivar++) {
                    resid[nvar*ivrt+ivar] -= xyz[ivar];
                    err = fabs(resid[nvar*ivrt+ivar]);

                    tree->cell[icel].rmserr += SQR(err);
                    *rmserr                 += SQR(err);
                    *maxerr = MAX(*maxerr,         err);

                    if (err > tree->cell[icel].maxerr) tree->cell[icel].maxerr = err;
                }
                break;
            }
        }
    }

    *rmserr = sqrt(*rmserr / nvrt);

 cleanup:
    FREE(cspln);
    FREE(vspln);
    FREE(uspln);
    FREE(xspln);
    FREE(xvrt );
    FREE(xyz  );

    DPRINT4("%s --> rmserr=%f, maxerr=%f, status=%d}",
            routine, *rmserr, *maxerr, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * initialTree2d -- set up initial Tree and global Cell                         *
 *                                                                              *
 ********************************************************************************
 */
extern int
initialTree2d(gridTree *tree,                /* (both) Grid Tree to refine */
              int      nvar,                 /* (in)   number of variables */
              int      nvrt,                 /* (in)   number of Vertices */
              int      periodic,             /* (in)   periodicity flag */
   /*@null@*/ int      ppnts[],              /* (in)   indices of periodic points */
              prmUV    uv[],                 /* (in)   array  of Vertices */
              double   resid[],              /* (both) array  of residuals */
              double   *rmserr,              /* (out)  RMS     error at Vertices */
              double   *maxerr)              /* (out)  maximum error at Vertices */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      *xvrt   = NULL;
    double      *xspln  = NULL;
    double      *uspln  = NULL;
    double      *vspln  = NULL;
    double      *cspln  = NULL;
    double      *xyz    = NULL;

    int         ncel, ivrt, ivar, iknt, icel, ii, jj;
    double      err, uu, vv, umin, umax, vmin, vmax, ss, tt;

    int         nvar2 = 2 * nvar;
    int         nvar3 = 3 * nvar;
    int         nvar4 = 4 * nvar;

    ROUTINE(initialTree2d);
    DPRINT4("%s(nvar=%d, nvrt=%d, periodic=%d) {",
            routine, nvar, nvrt, periodic);

    /* ----------------------------------------------------------------------- */

    *rmserr = 0;
    *maxerr = 0;

    MALLOC(xyz, double, nvar);

    /*
     * initialize the Tree
     */
    tree->nvar = nvar;
    tree->nu   = 2;
    tree->nv   = 2;
    tree->ncel = 1;
    tree->cell = NULL;
    tree->nknt = 4;
    tree->knot = NULL;

    /*
     * allocate and initialize the Cell and Knot arrays
     */
    MALLOC(tree->cell, gridCell, tree->ncel      );
    MALLOC(tree->knot, double,   tree->nknt*nvar4);

    for (iknt = 0; iknt < tree->nknt; iknt++) {
        for (ivar = 0; ivar < nvar4; ivar++) {
            tree->knot[nvar4*iknt+ivar] = 0;
        }
    }

    /*
     * set up root (global) Cell
     */
    ncel = 0;
    tree->cell[ncel].swKnot = 0;
    tree->cell[ncel].seKnot = 1;
    tree->cell[ncel].nwKnot = 2;
    tree->cell[ncel].neKnot = 3;
    tree->cell[ncel].wCell = -1;
    tree->cell[ncel].eCell = -1;
    tree->cell[ncel].sCell = -1;
    tree->cell[ncel].nCell = -1;
    tree->cell[ncel].dtype  =  0;
    tree->cell[ncel].child  =  0;
    tree->cell[ncel].count  =  0;

    tree->cell[ncel].umin   =  0;
    tree->cell[ncel].umax   =  1;
    tree->cell[ncel].vmin   =  0;
    tree->cell[ncel].vmax   =  1;
    tree->cell[ncel].rmserr =  0;
    tree->cell[ncel].maxerr =  0;

    /*
     * if periodic, then refine the global Cell (twice).
     *
     *  note: ncel will be the first "fine" Cell
     */
    if        (periodic == 0) {
        tree->nu = 2;
        tree->nv = 2;
        ncel     = 0;
    } else if (periodic == 1) {
        tree->nu = 5;
        tree->nv = 2;
        ncel     = 3;
        for (icel = 0; icel < ncel; icel++) {
            status = divideCell(tree, icel, 1);
            CHECK_STATUS;
        }
    } else if (periodic == 2) {
        tree->nu = 2;
        tree->nv = 5;
        ncel     = 3;
        for (icel = 0; icel < ncel; icel++) {
            status = divideCell(tree, icel, 2);
            CHECK_STATUS;
        }
    } else if (periodic == 3) {
        tree->nu = 5;
        tree->nv = 5;
        ncel     = 5;
        for (icel = 0; icel < ncel; icel++) {
            status = divideCell(tree, icel, 3);
            CHECK_STATUS;
        }
    }

    /*
     * count number of Vertices in each Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        icel = finestCell(tree, uv[ivrt].u, uv[ivrt].v);
        tree->cell[icel].count++;
    }

    /*
     * allocate storage for (extended) Vertex table and spline data
     */
    MALLOC(xvrt,  double,       nvrt);
    MALLOC(xspln, double, tree->nknt);
    MALLOC(uspln, double, tree->nknt);
    MALLOC(vspln, double, tree->nknt);
    MALLOC(cspln, double, tree->nknt);

#ifdef DEBUG
    printTree(dbg_fp, tree);
#endif

    /*
     * create splines for each of the dependent variables
     */
    for (ivar = 0; ivar < nvar; ivar++) {

        /*
         * copy the Vertex data into xvrt
         */
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            xvrt[ivrt] = resid[nvar*ivrt+ivar];
        }

        /*
         * find the spline
         */
        status = splineWithGaps2d(nvrt, uv, xvrt, periodic, ppnts,
                                  tree->nu, tree->nv, xspln, uspln, vspln, cspln);
        CHECK_STATUS;

        /*
         * store the spline data in the knots
         */
        for (icel = ncel; icel < tree->ncel; icel++) {
            ii = tree->nu * tree->cell[icel].umin;
            jj = tree->nv * tree->cell[icel].vmin;

            iknt = tree->cell[icel].swKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii  )+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii  )+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii  )+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii  )+tree->nu*(jj  )];

            iknt = tree->cell[icel].seKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii+1)+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii+1)+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii+1)+tree->nu*(jj  )];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii+1)+tree->nu*(jj  )];

            iknt = tree->cell[icel].nwKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii  )+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii  )+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii  )+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii  )+tree->nu*(jj+1)];

            iknt = tree->cell[icel].neKnot;
            tree->knot[nvar4*iknt      +ivar] = xspln[(ii+1)+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar +ivar] = uspln[(ii+1)+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar2+ivar] = vspln[(ii+1)+tree->nu*(jj+1)];
            tree->knot[nvar4*iknt+nvar3+ivar] = cspln[(ii+1)+tree->nu*(jj+1)];
        }
    }

    /*
     * subtract the current spline from "resid" and keep
     *    track of the max error associated with the global Cell(s)
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uu = uv[ivrt].u;
        vv = uv[ivrt].v;

        icel = finestCell(tree, uu, vv);
        umin = tree->cell[icel].umin;
        umax = tree->cell[icel].umax;
        vmin = tree->cell[icel].vmin;
        vmax = tree->cell[icel].vmax;
        ss   = (uu - umin) / (umax - umin);
        tt   = (vv - vmin) / (vmax - vmin);

        evalBicubic(ss, tt, nvar,
                    &(tree->knot[nvar4*(tree->cell[icel].swKnot)]),
                    &(tree->knot[nvar4*(tree->cell[icel].seKnot)]),
                    &(tree->knot[nvar4*(tree->cell[icel].nwKnot)]),
                    &(tree->knot[nvar4*(tree->cell[icel].neKnot)]),
                    xyz, NULL, NULL, NULL, NULL, NULL);

        for (ivar = 0; ivar < nvar; ivar++) {
            resid[nvar*ivrt+ivar] -= xyz[ivar];
            err = fabs(resid[nvar*ivrt+ivar]);

            tree->cell[icel].rmserr += SQR(err);
            *rmserr                 += SQR(err);
            *maxerr = MAX(*maxerr,         err);
            DPRINT5("%5d %10.5f %10.5f  %12.5e %12.5e", ivrt, uu, vv, err, *maxerr);

            if (err > tree->cell[icel].maxerr) tree->cell[icel].maxerr = err;
        }
    }

    *rmserr = sqrt(*rmserr / nvrt);

 cleanup:
    FREE(cspln);
    FREE(vspln);
    FREE(uspln);
    FREE(xspln);
    FREE(xvrt );
    FREE(xyz  );

    DPRINT4("%s --> rmserr=%f, maxerr=%f, status=%d}",
            routine, *rmserr, *maxerr, status);
    return status;
}



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * plotGrid -- plot the Grid of Knots (using GRAFIC)                            *
 *                                                                              *
 ********************************************************************************
 */
static void
plotGrid(int      nu,                        /* (in)   number of U grid Knots */
         int      nv,                        /* (in)   number of V grid Knots */
         int      nvar,                      /* (in)   number of variables */
         int      nvrt,                      /* (in)   number of Vertices */
         prmUV    uv[],                      /* (in)   array  of Vertices */
         double   grid[],                    /* (in)   grid of Knots */
         double   ddu[],                     /* (in)   grid of Knot U-derivatives */
         double   ddv[])                     /* (in)   grid of Knot V-derivatives */
{
    int         indgr  = 1 + 4 + 16 + 64;
    int         ivar   = 0;
    int         icut   = 1;

    ROUTINE(plotGrid);
    DPRINT4("%s(nu=%d, nv=%d, nvar=%d) {",
            routine, nu, nv, nvar);

    /* ----------------------------------------------------------------------- */

    grctrl_(plotGridImage, &indgr, "~U~V~Contours of X",
            (void*)(&nu), (void*)(&nv), (void*)(&nvar), (void*)(grid),
            (void*)(&nvrt), (void*)(uv),
            (void*)(ddu), (void*)(ddv),
            (void*)(&ivar), (void*)(&icut),
            strlen("~U~V~Contours of X"));

// cleanup:
    DPRINT1("%s}", routine);
}
#endif



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * plotGridImage -- used by plotGrid                                            *
 *                                                                              *
 ********************************************************************************
 */
static void
plotGridImage(int   *ifunct,                 /* (in)   GRAFIC function indicator */
              void  *nuP,                    /* (in)   number of U grid Knots */
              void  *nvP,                    /* (in)   number of V grid Knots */
              void  *nvarP,                  /* (in)   number of variables */
              void  *gridP,                  /* (in)   grid of Knots */
              void  *nvrtP,                  /* (in)   number of Vertices */
              void  *uvP,                    /* (in)   array  of Vertices */
              void  *dduP,                   /* (in)   grid of Knot U-derivatives */
              void  *ddvP,                   /* (in)   grid of Knot V-derivatives */
              void  *ivarP,                  /* (both) variable number to plot */
              void  *icutP,                  /* (both) =0 plot contours */
                                             /*        =1 plot vs U */
                                             /*        =2 plot vs V */
              float *scale,                  /* (out)  array of scales */
              char  *text,                   /* (out)  help text */
              int   textlen)                 /* (in)   length of text */
{
    int      *nu   = (int    *) nuP;
    int      *nv   = (int    *) nvP;
    int      *nvar = (int    *) nvarP;
    double   *grid = (double *) gridP;
    int      *nvrt = (int    *) nvrtP;
    prmUV    *uv   = (prmUV  *) uvP;
    double   *ddu  = (double *) dduP;
    double   *ddv  = (double *) ddvP;
    int      *ivar = (int    *) ivarP;
    int      *icut = (int    *) icutP;

    int      ncol = 13;

    int      iu, iv, iuv, ivrt, icol, idash;
    double   xmin, xmax;
    float    uu[4], vv[4], ff[4], col[ncol], frac, dum;
    char     pltitl[255];

    int      ione    = 1;
    int      ifour   = 4;
    int      iblack  = GR_BLACK;
    int      icircle = GR_CIRCLE;

//    ROUTINE(plotGridImage);

    /* ----------------------------------------------------------------------- */

    /* ---------- return scales ----------*/
    if (*ifunct == 0) {
        if        (*icut == 0 || *icut == 1) {
            scale[0] = 0;
            scale[1] = 1;
            scale[2] = 0;
            scale[3] = 1;
        } else if (*icut == 2 || *icut == 3) {
            xmin = grid[(*nvar)*0+(*ivar)];
            xmax = grid[(*nvar)*0+(*ivar)];

            for (iv = 0; iv < *nv; iv++) {
                for (iu = 0; iu < *nu; iu++) {
                    iuv = (*nu) * iv + iu;
                    xmin = MIN(xmin, grid[(*nvar)*iuv+(*ivar)]);
                    xmax = MAX(xmax, grid[(*nvar)*iuv+(*ivar)]);
                }
            }

            scale[0] = 0;
            scale[1] = 1;
            scale[2] = xmin;
            scale[3] = xmax;
        } else if (*icut == 4 || *icut == 5) {
            xmin = ddu[(*nvar)*0+(*ivar)];
            xmax = ddu[(*nvar)*0+(*ivar)];

            for (iv = 0; iv < *nv; iv++) {
                for (iu = 0; iu < *nu; iu++) {
                    iuv = (*nu) * iv + iu;
                    xmin = MIN(xmin, ddu[(*nvar)*iuv+(*ivar)]);
                    xmax = MAX(xmax, ddu[(*nvar)*iuv+(*ivar)]);
                }
            }

            scale[0] = 0;
            scale[1] = 1;
            scale[2] = xmin;
            scale[3] = xmax;
        } else if (*icut == 6 || *icut == 7) {
            xmin = ddv[(*nvar)*0+(*ivar)];
            xmax = ddv[(*nvar)*0+(*ivar)];

            for (iv = 0; iv < *nv; iv++) {
                for (iu = 0; iu < *nu; iu++) {
                    iuv = (*nu) * iv + iu;
                    xmin = MIN(xmin, ddv[(*nvar)*iuv+(*ivar)]);
                    xmax = MAX(xmax, ddv[(*nvar)*iuv+(*ivar)]);
                }
            }

            scale[0] = 0;
            scale[1] = 1;
            scale[2] = xmin;
            scale[3] = xmax;
        }

        strncpy(text, "C G N", textlen-1);

    /* ---------- plot image ---------- */
    } else if (*ifunct == 1) {
        if (*icut == 0 || *icut == 1) {
            xmin = grid[(*nvar)*0+(*ivar)];
            xmax = grid[(*nvar)*0+(*ivar)];

            for (iv = 0; iv < *nv; iv++) {
                for (iu = 0; iu < *nu; iu++) {
                    iuv = (*nu) * iv + iu;
                    xmin = MIN(xmin, grid[(*nvar)*iuv+(*ivar)]);
                    xmax = MAX(xmax, grid[(*nvar)*iuv+(*ivar)]);
                }
            }

            for (icol = 0; icol < ncol; icol++) {
                frac      = (float)(icol) / (float)(ncol - 1);
                col[icol] = (1 - frac) * xmin + frac * xmax;
            }

            for (iv = 1; iv < *nv; iv++) {
                for (iu = 1; iu < *nu; iu++) {
                    iuv = (*nu) * iv + iu;

                    uu[0] = (float)(iu-1) / (float)(*nu-1);
                    vv[0] = (float)(iv-1) / (float)(*nv-1);
                    ff[0] = grid[(*nvar)*((*nu)*(iv-1)+(iu-1))+(*ivar)];

                    uu[1] = (float)(iu  ) / (float)(*nu-1);
                    vv[1] = (float)(iv-1) / (float)(*nv-1);
                    ff[1] = grid[(*nvar)*((*nu)*(iv-1)+(iu  ))+(*ivar)];

                    uu[2] = (float)(iu  ) / (float)(*nu-1);
                    vv[2] = (float)(iv  ) / (float)(*nv-1);
                    ff[2] = grid[(*nvar)*((*nu)*(iv  )+(iu  ))+(*ivar)];

                    uu[3] = (float)(iu-1) / (float)(*nu-1);
                    vv[3] = (float)(iv  ) / (float)(*nv-1);
                    ff[3] = grid[(*nvar)*((*nu)*(iv  )+(iu-1))+(*ivar)];

                    grcol2_(uu, vv, ff, &ifour, col, &ncol);
                }
            }

            if (*icut == 1) {
                grcolr_(&iblack);

                for (ivrt = 0; ivrt < *nvrt; ivrt++) {
                    uu[0] = uv[ivrt].u;
                    vv[0] = uv[ivrt].v;
                    grmov2_(uu, vv);
                    grsymb_(&icircle);
                }
            }

        } else if (*icut == 2) {
            for (iv = 0; iv < *nv; iv++) {
                idash = 1 + iv % 5;
                grdash_(&idash);

                iu    = 0;
                uu[0] = (float)(iu) / (float)(*nu-1);
                ff[0] = grid[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                grmov2_(uu, ff);
                grsymb_(&idash);

                for (iu = 1; iu < *nu; iu++) {
                    uu[0] = (float)(iu) / (float)(*nu-1);
                    ff[0] = grid[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                    grdrw2_(uu, ff);
                grsymb_(&idash);
                }
            }

            idash = 1;
            grdash_(&idash);
        } else if (*icut == 3) {
            for (iu = 0; iu < *nu; iu++) {
                idash = 1 + iu % 5;
                grdash_(&idash);

                iv    = 0;
                vv[0] = (float)(iv) / (float)(*nv-1);
                ff[0] = grid[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                grmov2_(vv, ff);
                grsymb_(&idash);

                for (iv = 1; iv < *nv; iv++) {
                    vv[0] = (float)(iv) / (float)(*nv-1);
                    ff[0] = grid[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                    grdrw2_(vv, ff);
                grsymb_(&idash);
                }
            }

            idash = 1;
            grdash_(&idash);
        } else if (*icut == 4) {
            for (iv = 0; iv < *nv; iv++) {
                idash = 1 + iv % 5;
                grdash_(&idash);

                iu    = 0;
                uu[0] = (float)(iu) / (float)(*nu-1);
                ff[0] = ddu[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                grmov2_(uu, ff);
                grsymb_(&idash);

                for (iu = 1; iu < *nu; iu++) {
                    uu[0] = (float)(iu) / (float)(*nu-1);
                    ff[0] = ddu[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                    grdrw2_(uu, ff);
                grsymb_(&idash);
                }
            }

            idash = 1;
            grdash_(&idash);
        } else if (*icut == 5) {
            for (iu = 0; iu < *nu; iu++) {
                idash = 1 + iu % 5;
                grdash_(&idash);

                iv    = 0;
                vv[0] = (float)(iv) / (float)(*nv-1);
                ff[0] = ddu[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                grmov2_(vv, ff);
                grsymb_(&idash);

                for (iv = 1; iv < *nv; iv++) {
                    vv[0] = (float)(iv) / (float)(*nv-1);
                    ff[0] = ddu[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                    grdrw2_(vv, ff);
                grsymb_(&idash);
                }
            }

            idash = 1;
            grdash_(&idash);
        } else if (*icut == 6) {
            for (iv = 0; iv < *nv; iv++) {
                idash = 1 + iv % 5;
                grdash_(&idash);

                iu    = 0;
                uu[0] = (float)(iu) / (float)(*nu-1);
                ff[0] = ddv[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                grmov2_(uu, ff);
                grsymb_(&idash);

                for (iu = 1; iu < *nu; iu++) {
                    uu[0] = (float)(iu) / (float)(*nu-1);
                    ff[0] = ddv[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                    grdrw2_(uu, ff);
                grsymb_(&idash);
                }
            }

            idash = 1;
            grdash_(&idash);
        } else if (*icut == 7) {
            for (iu = 0; iu < *nu; iu++) {
                idash = 1 + iu % 5;
                grdash_(&idash);

                iv    = 0;
                vv[0] = (float)(iv) / (float)(*nv-1);
                ff[0] = ddv[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                grmov2_(vv, ff);
                grsymb_(&idash);

                for (iv = 1; iv < *nv; iv++) {
                    vv[0] = (float)(iv) / (float)(*nv-1);
                    ff[0] = ddv[(*nvar)*((*nu)*(iv)+(iu))+(*ivar)];
                    grdrw2_(vv, ff);
                grsymb_(&idash);
                }
            }

            idash = 1;
            grdash_(&idash);
        }

   /* ---------- toggle icut ---------- */
    } else if (*ifunct == -3) {
        *icut = (*icut + 1) % 8;
        printf("setting icut=%d\n", *icut);

        if        (*icut == 0 || *icut == 1) {
            if (*ivar == 0) sprintf(pltitl, "~U~V~Contours of X");
            if (*ivar == 1) sprintf(pltitl, "~U~V~Contours of Y");
            if (*ivar == 2) sprintf(pltitl, "~U~V~Contours of Z");
        } else if (*icut == 2) {
            if (*ivar == 0) sprintf(pltitl, "~U~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~Z~ ");
        } else if (*icut == 3) {
            if (*ivar == 0) sprintf(pltitl, "~V~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~Z~ ");
        } else if (*icut == 4) {
            if (*ivar == 0) sprintf(pltitl, "~U~dXdU~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~dYdU~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~dZdU~ ");
        } else if (*icut == 5) {
            if (*ivar == 0) sprintf(pltitl, "~V~dXdU~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~dYdU~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~dZdU~ ");
        } else if (*icut == 6) {
            if (*ivar == 0) sprintf(pltitl, "~U~dXdV~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~dYdV~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~dZdV~ ");
        } else if (*icut == 7) {
            if (*ivar == 0) sprintf(pltitl, "~V~dXdV~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~dYdV~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~dZdV~ ");
        }
        grvalu_("LABLGR", &ione, &dum, pltitl, strlen("LABLGR"), strlen(pltitl));
        grscpt_(&ione, "O", strlen("O"));

   /* ---------- toggle ivar ---------- */
    } else if (*ifunct == -14) {
        *ivar = (*ivar + 1) % (*nvar);
        printf("setting ivar=%d\n", *ivar);

        if        (*icut == 0 || *icut == 1) {
            if (*ivar == 0) sprintf(pltitl, "~U~V~Contours of X");
            if (*ivar == 1) sprintf(pltitl, "~U~V~Contours of Y");
            if (*ivar == 2) sprintf(pltitl, "~U~V~Contours of Z");
        } else if (*icut == 2) {
            if (*ivar == 0) sprintf(pltitl, "~U~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~Z~ ");
        } else if (*icut == 3) {
            if (*ivar == 0) sprintf(pltitl, "~V~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~Z~ ");
        } else if (*icut == 4) {
            if (*ivar == 0) sprintf(pltitl, "~U~dXdU~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~dYdU~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~dZdU~ ");
        } else if (*icut == 5) {
            if (*ivar == 0) sprintf(pltitl, "~V~dXdU~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~dYdU~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~dZdU~ ");
        } else if (*icut == 6) {
            if (*ivar == 0) sprintf(pltitl, "~U~dXdV~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~dYdV~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~dZdV~ ");
        } else if (*icut == 7) {
            if (*ivar == 0) sprintf(pltitl, "~V~dXdV~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~dYdV~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~dZdV~ ");
        }
        grvalu_("LABLGR", &ione, &dum, pltitl, strlen("LABLGR"), strlen(pltitl));
        grscpt_(&ione, "O", strlen("O"));
    }
}
#endif



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * plotLevels -- plot the Grid at various levels (using GRAFIC)                 *
 *                                                                              *
 ********************************************************************************
 */
static void
plotLevels(gridTree *tree)                   /* (in)   pointer to Tree */
{
    int         indgr  = 1 + 4 + 16 + 64;
    int         ivar   = 0;
    int         icut   = 0;
    int         nlvl   = 1;
    double      xmin, xmax;

    ROUTINE(plotLevels);
    DPRINT1("%s() {",
            routine);

    /* ----------------------------------------------------------------------- */

    grctrl_(plotLevelsImage, &indgr, "~U~V~Contours of X",
            (void*)(tree), (void*)(&ivar), (void*)(&icut), (void*)(&nlvl),
            (void*)(&xmin), (void*)(&xmax), NULL, NULL, NULL, NULL,
            strlen("~U~V~Contours of X"));

// cleanup:
    DPRINT1("%s}", routine);
}
#endif



#ifdef GRAFIC
/*
 ********************************************************************************
 *                                                                              *
 * plotLevelsImage -- used by plotLevels                                        *
 *                                                                              *
 ********************************************************************************
 */
static void
plotLevelsImage(int   *ifunct,               /* (in)   GRAFIC function indicator */
                void  *treeP,                /* (in)   pointer to Tree */
                void  *ivarP,                /* (both) variable number to plot */
                void  *icutP,                /* (both) =0 plot contours */
                                             /*        =1 plot vs U */
                                             /*        =2 plot vs V */
                void  *nlvlP,                /* (both) maximum level to plot */
                void  *xminP,                   /* (in)   dummy GRAFIC argument */
                void  *xmaxP,                   /* (in)   dummy GRAFIC argument */
                void  *a6,                   /* (in)   dummy GRAFIC argument */
                void  *a7,                   /* (in)   dummy GRAFIC argument */
                void  *a8,                   /* (in)   dummy GRAFIC argument */
                void  *a9,                   /* (in)   dummy GRAFIC argument */
                float *scale,                /* (out)  array of scales */
                char  *text,                 /* (out)  help text */
                int   textlen)               /* (in)   length of text */
{
    gridTree *tree = (gridTree*) treeP;
    int      *ivar = (int    *) ivarP;
    int      *icut = (int    *) icutP;
    int      *nlvl = (int    *) nlvlP;
    double   *xmin = (double *) xminP;
    double   *xmax = (double *) xmaxP;

    int      ncol = 13;

    prmUV    uv;
    int      nvar, nu, nv, iu, iv, icol, idash, icolr;
    double   xx;
    float    uu[4], vv[4], ff[4], col[ncol], frac, dum;
    char     pltitl[255];

    int      ione    = 1;
    int      ifour   = 4;
    int      iblack  = GR_BLACK;

//    ROUTINE(plotLevelsImage);

    /* ----------------------------------------------------------------------- */

    nvar = tree->nvar;
    nu   = tree->nu;
    nv   = tree->nv;

    /* ---------- return scales ----------*/
    if (*ifunct == 0) {
        *xmin = +HUGEQ;
        *xmax = -HUGEQ;

        for (iv = 0; iv < nv; iv++) {
            for (iu = 0; iu < nu; iu++) {
                uv.u = (double)(iu) / (double)(nu-1);
                uv.v = (double)(iv) / (double)(nv-1);

                xx = evalToLevel(tree, 99, *ivar, uv);

                if (xx < *xmin) *xmin = xx;
                if (xx > *xmax) *xmax = xx;
            }
        }

        if (*icut == 0) {
            scale[0] = 0;
            scale[1] = 1;
            scale[2] = 0;
            scale[3] = 1;
        } else if (*icut == 1 || *icut == 2) {
            scale[0] = 0;
            scale[1] = 1;
            scale[2] = *xmin;
            scale[3] = *xmax;
        }

        strncpy(text, "C L N", textlen-1);

    /* ---------- plot image ---------- */
    } else if (*ifunct == 1) {
        if (*icut == 0) {
            for (icol = 0; icol < ncol; icol++) {
                frac      = (float)(icol) / (float)(ncol - 1);
                col[icol] = (1 - frac) * (*xmin) + frac * (*xmax);
            }

            for (iv = 1; iv < nv; iv++) {
                for (iu = 1; iu < nu; iu++) {
                    uv.u  = (double)(iu-1) / (double)(nu-1);
                    uv.v  = (double)(iv-1) / (double)(nv-1);
                    uu[0] = uv.u;
                    vv[0] = uv.v;
                    ff[0] = evalToLevel(tree, *nlvl, *ivar, uv);

                    uv.u  = (double)(iu  ) / (double)(nu-1);
                    uv.v  = (double)(iv-1) / (double)(nv-1);
                    uu[1] = uv.u;
                    vv[1] = uv.v;
                    ff[1] = evalToLevel(tree, *nlvl, *ivar, uv);

                    uv.u  = (double)(iu  ) / (double)(nu-1);
                    uv.v  = (double)(iv  ) / (double)(nv-1);
                    uu[2] = uv.u;
                    vv[2] = uv.v;
                    ff[2] = evalToLevel(tree, *nlvl, *ivar, uv);

                    uv.u  = (double)(iu-1) / (double)(nu-1);
                    uv.v  = (double)(iv  ) / (double)(nv-1);
                    uu[3] = uv.u;
                    vv[3] = uv.v;
                    ff[3] = evalToLevel(tree, *nlvl, *ivar, uv);

                    grcol2_(uu, vv, ff, &ifour, col, &ncol);
                }
            }
        } else if (*icut == 1) {
            for (iv = 0; iv < nv; iv++) {
                idash = 1 +      iv % 5;
                icolr = 1 + 2 * (iv % 5);
                grdash_(&idash);
                grcolr_(&icolr);

                uv.u  = 0;
                uv.v  = (double)(iv) / (double)(nv-1);
                uu[0] = uv.u;
                ff[0] = evalToLevel(tree, *nlvl, *ivar, uv);
                grmov2_(uu, ff);

                for (iu = 1; iu < nu; iu++) {
                    uv.u  = (double)(iu) / (double)(nu-1);
                    uu[0] = uv.u;
                    ff[0] = evalToLevel(tree, *nlvl, *ivar, uv);
                    grdrw2_(uu, ff);
                }
            }

            idash = 1;
            grdash_(&idash);
            grcolr_(&iblack);
        } else if (*icut == 2) {
            for (iu = 0; iu < nu; iu++) {
                idash = 1 +      iu % 5;
                icolr = 1 + 2 * (iu % 5);
                grdash_(&idash);
                grcolr_(&icolr);

                uv.u  = (double)(iu) / (double)(nu-1);
                uv.v  = 0;
                vv[0] = uv.v;
                ff[0] = evalToLevel(tree, *nlvl, *ivar, uv);
                grmov2_(vv, ff);

                for (iv = 1; iv < nv; iv++) {
                    uv.v  = (double)(iv) / (double)(nv-1);
                    vv[0] = uv.v;
                    ff[0] = evalToLevel(tree, *nlvl, *ivar, uv);
                    grdrw2_(vv, ff);
                }
            }

            idash = 1;
            grdash_(&idash);
            grcolr_(&iblack);
        }

   /* ---------- toggle icut ---------- */
    } else if (*ifunct == -3) {
        *icut = (*icut + 1) % 3;

        if        (*icut == 0) {
            if (*ivar == 0) sprintf(pltitl, "~U~V~Contours of X");
            if (*ivar == 1) sprintf(pltitl, "~U~V~Contours of Y");
            if (*ivar == 2) sprintf(pltitl, "~U~V~Contours of Z");
        } else if (*icut == 1) {
            if (*ivar == 0) sprintf(pltitl, "~U~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~Z~ ");
        } else if (*icut == 2) {
            if (*ivar == 0) sprintf(pltitl, "~V~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~Z~ ");
        }
        grvalu_("LABLGR", &ione, &dum, pltitl, strlen("LABLGR"), strlen(pltitl));
        grscpt_(&ione, "O", strlen("O"));

   /* ---------- toggle nlvl ---------- */
    } else if (*ifunct == -12) {
        *nlvl = (*nlvl + 1);

        if ((pow(2, *nlvl-1)+1 > nu) && (pow(2, *nlvl-1)+1 > nv)) {
            *nlvl = 1;
        }
        printf("setting *nlvl=%d\n", *nlvl);

        grscpt_(&ione, "O", strlen("O"));

   /* ---------- toggle ivar ---------- */
    } else if (*ifunct == -14) {
        *ivar = (*ivar + 1) % (nvar);

        if        (*icut == 0) {
            if (*ivar == 0) sprintf(pltitl, "~U~V~Contours of X");
            if (*ivar == 1) sprintf(pltitl, "~U~V~Contours of Y");
            if (*ivar == 2) sprintf(pltitl, "~U~V~Contours of Z");
        } else if (*icut == 1) {
            if (*ivar == 0) sprintf(pltitl, "~U~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~U~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~U~Z~ ");
        } else if (*icut == 2) {
            if (*ivar == 0) sprintf(pltitl, "~V~X~ ");
            if (*ivar == 1) sprintf(pltitl, "~V~Y~ ");
            if (*ivar == 2) sprintf(pltitl, "~V~Z~ ");
        }
        grvalu_("LABLGR", &ione, &dum, pltitl, strlen("LABLGR"), strlen(pltitl));
        grscpt_(&ione, "O", strlen("O"));
    }
}
#endif



#ifdef GRAFIC2
/*
 ********************************************************************************
 *                                                                              *
 * plotVrts -- plot the UV coordinates of Vertices (using GRAFIC)               *
 *                                                                              *
 ********************************************************************************
 */
static void
plotVrts(int      nvrt,                      /* (in)   number of Vertices */
         int      nvar,                      /* (in)   number of Variables */
         prmUV    uv[],                      /* (in)   array  of Parameters */
         double   xyz[])                     /* (in)   array  of Coordinates */
{
    int         io_kbd=5, io_scr=6, indgr=1+2+4+16+64;
    int         ivar, ivrt;
    double      umin, umax, vmin, vmax, xmin[99], xmax[99];

    ROUTINE(plotVrts);
    DPRINT2("%s(nvrt=%d) {",
            routine, nvrt);

    /* ----------------------------------------------------------------------- */

    umin = uv[0].u;
    umax = uv[0].u;
    vmin = uv[0].v;
    vmax = uv[0].v;
    for (ivar = 0; ivar < nvar; ivar++) {
        xmin[ivar] = xyz[ivar];
        xmax[ivar] = xyz[ivar];
    }

    for (ivrt = 1; ivrt < nvrt; ivrt++) {
        umin = MIN(umin, uv[ivrt].u);
        umax = MAX(umax, uv[ivrt].u);
        vmin = MIN(vmin, uv[ivrt].v);
        vmax = MAX(vmax, uv[ivrt].v);
        for (ivar = 0; ivar < nvar; ivar++) {
                xmin[ivar] = MIN(xmin[ivar], xyz[nvar*ivrt+ivar]);
                xmax[ivar] = MAX(xmax[ivar], xyz[nvar*ivrt+ivar]);
        }
    }

    printf("extrema of data being plotted:");
    printf("u   : %14.7f %14.7f\n", umin, umax);
    printf("v   : %14.7f %14.7f\n", vmin, vmax);
    for (ivar = 0; ivar < nvar; ivar++) {
        printf("x[%d]: %14.7f %14.7f\n", ivar, xmin[ivar], xmax[ivar]);
    }

    ivar = 0;

    grinit_(&io_kbd, &io_scr, "plotVrts", strlen("plotVrts"));
    grctrl_(plotVrtsImage, &indgr, "~U~V~Vertices colored by X",
            (void*)(&nvrt), (void*)(&nvar), (void*)(uv), (void*)(xyz),
            (void*)(&ivar), NULL, NULL, NULL, NULL, NULL,
            strlen("~U~V~Vertices colored by X"));

// cleanup:
    DPRINT1("%s}", routine);
}
#endif



#ifdef GRAFIC2
/*
 ********************************************************************************
 *                                                                              *
 * plotVrtsImage -- used by plotVrts                                            *
 *                                                                              *
 ********************************************************************************
 */
static void
plotVrtsImage(int   *ifunct,                 /* (in)   GRAFIC function indicator */
              void  *nvrtP,                  /* (in)   number of Vertices */
              void  *nvarP,                  /* (in)   number of variables */
              void  *uvP,                    /* (in)   array  of Parametric Coords */
              void  *xyzP,                   /* (in)   array  of physical Coords */
              void  *ivarP,                  /* (both) variable number */
              void  *a5,                     /* (in)   dummy GRAFIC argument */
              void  *a6,                     /* (in)   dummy GRAFIC argument */
              void  *a7,                     /* (in)   dummy GRAFIC argument */
              void  *a8,                     /* (in)   dummy GRAFIC argument */
              void  *a9,                     /* (in)   dummy GRAFIC argument */
              float *scale,                  /* (out)  array of scales */
              char  *text,                   /* (out)  help text */
              int   textlen)                 /* (in)   length of text */
{
    int    *nvar   = (int    *) nvarP;
    int    *nvrt   = (int    *) nvrtP;
    prmUV  *uv     = (prmUV  *) uvP;
    double *xyz    = (double *) xyzP;
    int    *ivar   = (int    *) ivarP;

    int    ivrt, icol;
    double umin, umax, vmin, vmax, xmin, xmax;
    float  uu, vv, dum;
    char   pltitl[255];

    double eps     = 1.0e-6;
    int    icircle = GR_CIRCLE;
    int    ione    = 1;

//    ROUTINE(plotVrtsImage);

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

        strncpy(text, "N", textlen-1);

    /* ---------- plot image ---------- */
    } else if (*ifunct == 1) {
        xmin = xyz[(*nvar)*0+(*ivar)];
        xmax = xyz[(*nvar)*0+(*ivar)];

        for (ivrt = 1; ivrt < *nvrt; ivrt++) {
            xmin = MIN(xmin, xyz[(*nvar)*ivrt+(*ivar)]);
            xmax = MAX(xmax, xyz[(*nvar)*ivrt+(*ivar)]);
        }

        for (ivrt = 0; ivrt < *nvrt; ivrt++) {
            icol = 1 + 13 * (xyz[(*nvar)*ivrt+(*ivar)] - xmin) / (xmax + eps - xmin);
            grcolr_(&icol);

            uu = uv[ivrt].u;
            vv = uv[ivrt].v;

            grmov2_(&uu, &vv);
            grsymb_(&icircle);
        }

   /* ---------- next variable ---------- */
    } else if (*ifunct == -14) {
        *ivar = (*ivar + 1) % (*nvar);
        printf("setting ivar=%d\n", *ivar);

        if (*ivar == 0) sprintf(pltitl, "~U~V~Vertices colored by X");
        if (*ivar == 1) sprintf(pltitl, "~U~V~Vertices colored by Y");
        if (*ivar == 2) sprintf(pltitl, "~U~V~Vertices colored by Z");

        grvalu_("LABLGR", &ione, &dum, pltitl, strlen("LABLGR"), strlen(pltitl));
        grscpt_(&ione, "R", strlen("R"));
    }
}
#endif



#ifdef DEBUG
/*
 ********************************************************************************
 *                                                                              *
 * printTree -- print the current Tree                                          *
 *                                                                              *
 ********************************************************************************
 */
static void
printTree(FILE     *fp,                      /* (in)   file pointer */
          gridTree *tree)                    /* (in)   Grid Tree */
{

    int         icel, iknt, ivar;

    int         nvar  =      tree->nvar;
    int         nvar2 =  2 * tree->nvar;
    int         nvar3 =  3 * tree->nvar;
    int         nvar4 =  4 * tree->nvar;

    ROUTINE(printTree);
    DPRINT1("%s() {",
            routine);

    /* ----------------------------------------------------------------------- */

    fprintf(fp, "        ------ Knot ------   ------ Cell ------\n");
    fprintf(fp, " icel    sw   se   nw   ne     w    e    s    n  dtype chld count     umin    umax    vmin    vmax     maxerr\n");
    for (icel = 0; icel < tree->ncel; icel++) {
        fprintf(fp, "%5d %5d%5d%5d%5d %5d%5d%5d%5d  %5d%5d%6d %8.5f%8.5f%8.5f%8.5f %12.7f\n", icel,
                tree->cell[icel].swKnot,
                tree->cell[icel].seKnot,
                tree->cell[icel].nwKnot,
                tree->cell[icel].neKnot,
                tree->cell[icel].wCell,
                tree->cell[icel].eCell,
                tree->cell[icel].sCell,
                tree->cell[icel].nCell,
                tree->cell[icel].dtype,
                tree->cell[icel].child,
                tree->cell[icel].count,
                tree->cell[icel].umin,
                tree->cell[icel].umax,
                tree->cell[icel].vmin,
                tree->cell[icel].vmax,
                tree->cell[icel].maxerr);
    }

    fprintf(fp, " iknt\n");
    for (iknt = 0; iknt < tree->nknt; iknt++) {
        fprintf(fp, "%5d ", iknt);
        for (ivar = 0; ivar < nvar; ivar++) {
            fprintf(fp, "%10.5f", tree->knot[nvar4*iknt      +ivar]);
        }
        fprintf(fp, "  ");
        for (ivar = 0; ivar < nvar; ivar++) {
            fprintf(fp, "%10.5f", tree->knot[nvar4*iknt+nvar +ivar]);
        }
        fprintf(fp, "  ");
        for (ivar = 0; ivar < nvar; ivar++) {
            fprintf(fp, "%10.5f", tree->knot[nvar4*iknt+nvar2+ivar]);
        }
        fprintf(fp, "  ");
        for (ivar = 0; ivar < nvar; ivar++) {
            fprintf(fp, "%10.5f", tree->knot[nvar4*iknt+nvar3+ivar]);
        }
        fprintf(fp, "\n");
    }

// cleanup:
    DPRINT1("%s}", routine);
}
#endif



/*
 ********************************************************************************
 *                                                                              *
 * splineWithGaps2d -- create Ferguson spline that least-square fits 2D data    *
 *                                                                              *
 ********************************************************************************
 */
static int
splineWithGaps2d(int      nvrt,              /* (in)   number of Vertices to fit */
                 prmUV    uv[],              /* (in)   UV for Vertices to fit */
                 double   xvrt[],            /* (in)   X (data) to fit */
                 int      periodic,          /* (in)   periodicity flag */
      /*@null@*/ int      ppnts[],           /* (in)   indices of periodic points */
                 int      nu,                /* (in)   number of U grid Knots */
                 int      nv,                /* (in)   number of V grid Knots */
                 double   x[],               /* (both) array of spline Knots */
                 double   du[],              /* (out)  array of u-derivs */
                 double   dv[],              /* (out)  array of v-derivs */
                 double   dc[])              /* (out)  array of cross derivs */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */

    double      *asmf = NULL;                /* sparse-matrix data */
    int         *ismf = NULL;                /* sparse-matrix indices */
    double      *rhs  = NULL;                /* right-hand sides */
    double      *xx   = NULL;                /* matrix solution: x, du, dv, dc */
    int         *nn   = NULL;                /* number of Vertices near each Knot */

    double      ss, As, Bs, Cs, Ds, uu;
    double      tt, At, Bt, Ct, Dt, vv;
    double      b00, b01, b02, b03;
    double      b10, b11, b12, b13;
    double      b20, b21, b22, b23;
    double      b30, b31, b32, b33;
    double      uleft, urite, vleft, vrite, xleft, xrite, uuu, vvv, xxx;

    int         i, j, ij, ii, jj, kk, ivrt, nuv, is, iw;
    int         isw, ise, inw, ine, jsw, jse, jnw, jne;
    int         nsmf, irow, icol, iuv, ileft, imidl, irite;

    double      amu    = 0.001;              /* smoothing coefficient */
    double      xtol   = 1.e-6;              /* tolerance on x for constant value */

    ROUTINE(splineWithGaps2d);
    DPRINT5("%s(nvrt=%d, periodic=%d, nu=%d, nv=%d) {",
            routine, nvrt, periodic, nu, nv);

    /* ----------------------------------------------------------------------- */

#ifdef GRAFIC2
    plotVrts(nvrt, 1, uv, xvrt);
#endif

    nuv = nu * nv;

    /*
     * find the count of Vertices adjacent to each Knot
     */
    MALLOC(nn, int, nuv);

    for (i = 0; i < nuv; i++) {
        nn[i] = 0;
    }

    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        i   = MIN(MAX(0, (nu-1)*uv[ivrt].u), nu-2);
        j   = MIN(MAX(0, (nv-1)*uv[ivrt].v), nv-2);
        iuv = (i) + nu * (j);

        nn[iuv     ]++;
        nn[iuv   +1]++;
        nn[iuv+nu  ]++;
        nn[iuv+nu+1]++;
    }

    /*
     * overall matrix equation is of form:
     *
     *     [      |      |      |      ] [     ]   [     ]
     *     [  M1  |  M2  |  M3  |  M4  ] [  x  ]   [ RHS ]  (minimize SSE)
     *     [      |      |      |      ] [     ]   [     ]
     *     [  - - - - - -| - - - - - - ] [ --- ]   [ --- ]
     *     [      |      |      |      ] [     ]   [     ]
     *     [  M5  |  M6  |  0   |  0   ] [ du  ]   [  0  ]  (spline x(u))
     *     [      |      |      |      ] [     ]   [     ]
     *     [  - - - - - -| - - - - - - ] [ --- ] = [ --- ]
     *     [      |      |      |      ] [     ]   [     ]
     *     [  M7  |  0   |  M8  |  0   ] [ dv  ]   [  0  ]  (spline x(v))
     *     [      |      |      |      ] [     ]   [     ]
     *     [  - - - - - -| - - - - - - ] [ --- ]   [ --- ]
     *     [      |      |      |      ] [     ]   [     ]
     *     [  0   |  0   |  M9  |  M10 ] [ duv ]   [  0  ]  (spline dv(u))
     *     [      |      |      |      ] [     ]   [     ]
     *
     *     set up the storage in sparse-matrix form: 51 = 9 + 9 + 9 + 9
     *                                                  + 2 + 3 + 0 + 0
     *                                                  + 2 + 0 + 3 + 0
     *                                                  + 0 + 0 + 2 + 3
     */

    MALLOC(asmf, double, 51*nuv);
    MALLOC(ismf, int,    51*nuv);
    MALLOC(rhs,  double,  4*nuv);
    MALLOC(xx,   double,  4*nuv);

    /*
     * start out with all zeros
     */
    for (i = 0; i < 51*nuv; i++) {
        asmf[i] = 0;
        ismf[i] = 0;
    }

    for (i = 0; i < 4*nuv; i++) {
        rhs[i] = 0;
        xx[ i] = 0;
    }

    /*
     * set up sparse-matrix index matrix as well as sparse-matrix
     *    data for the bottom blocks.  note that some of the
     *    off-diagonal matrix elements might turn out to be zero (which
     *    wastes a little space, but make the code below much simpler)
     */
      ismf[0] = 4*nuv + 1;
      nsmf    = 4*nuv;

    /*
     * ...blocks M1, M2, M3, and M4 implement least-squares fit -> x
     *    (set up indices only)
     */
    for (j = 0; j < nv; j++) {
        for (i = 0; i < nu; i++) {
            irow = (i) + nu * (j);

            for (kk = 0; kk < 4; kk++) {
                for (jj = MAX(0, j-1); jj <= MIN(nv-1, j+1); jj++) {
                    for (ii = MAX(0, i-1); ii <= MIN(nu-1, i+1); ii++) {
                        icol = (ii) + nu * (jj) + nuv * (kk);

                        if (irow != icol) {
                            nsmf++;
                            ismf[nsmf] = icol;
                            asmf[nsmf] = 0;
                        }
                    }
                }
            }

            ismf[irow+1] = nsmf + 1;
        }
    }

    /*
     *...M5 and M6 implement Ferguson spline for x(u) -> du
     *   (set up required matrix data and indices)
     */
    for (j = 0; j < nv; j++) {
        for (i = 0; i < nu; i++) {
            irow = (i) + nu * (j) + nuv;

            if (i == 0) {
                nsmf++;
                ismf[nsmf] = irow - nuv;
                asmf[nsmf] =  3;                             /* x( ij) */

                nsmf++;
                ismf[nsmf] = irow - nuv + 1;
                asmf[nsmf] = -3;                             /* x( ie) */

                nsmf++;
                ismf[nsmf] = irow + nu - 2;
                asmf[nsmf] =  0;                             /* periodic=1 */

                asmf[irow] =  2;                             /* du(ij) */

                nsmf++;
                ismf[nsmf] = irow       + 1;
                asmf[nsmf] =  1;                             /* du(ie) */
            } else if (i == nu-1) {
                nsmf++;
                ismf[nsmf] = irow - nuv - 1;
                asmf[nsmf] =  3;                             /* x( iw) */

                nsmf++;
                ismf[nsmf] = irow - nuv;
                asmf[nsmf] = -3;                             /* x( ij) */

                nsmf++;
                ismf[nsmf] = irow       - 1;
                asmf[nsmf] =  1;                             /* du(iw) */

                asmf[irow] =  2;                             /* du(ij) */

                nsmf++;
                ismf[nsmf] = irow - nu + 1;
                asmf[nsmf] =  0;                             /* periodic=1 */
            } else {
                nsmf++;
                ismf[nsmf] = irow - nuv - 1;
                asmf[nsmf] =  3;                             /* x( iw) */

                nsmf++;
                ismf[nsmf] = irow - nuv + 1;
                asmf[nsmf] = -3;                             /* x( ie) */

                nsmf++;
                ismf[nsmf] = irow       - 1;
                asmf[nsmf] =  1;                             /* du(iw) */

                asmf[irow] =  4;                             /* du(ij) */

                nsmf++;
                ismf[nsmf] = irow       + 1;
                asmf[nsmf] =  1;                             /* du(ie) */
            }

            ismf[irow+1] = nsmf + 1;
        }
    }

    /*
     *...M7 and M8 implement Ferguson spline for x(v) -> dv
     *   (set up required matrix data and indices)
     */
      for (j = 0; j < nv; j++) {
          for (i = 0; i < nu; i++) {
              irow = (i) + nu * (j) + 2 * nuv;

              if (j == 0) {
                  nsmf++;
                  ismf[nsmf] = irow - 2*nuv;
                  asmf[nsmf] =  3;                           /* x( ij) */

                  nsmf++;
                  ismf[nsmf] = irow - 2*nuv + nu;
                  asmf[nsmf] = -3;                           /* x( in) */

                  nsmf++;
                  ismf[nsmf] = irow +   nuv - 2*nu;
                  asmf[nsmf] =  0;                           /* periodic=2 */

                  asmf[irow] =  2;                           /* dv(ij) */

                  nsmf++;
                  ismf[nsmf] = irow         + nu;
                  asmf[nsmf] =  1;                           /* dv(in) */
              } else if (j == nv-1) {
                  nsmf++;
                  ismf[nsmf] = irow - 2*nuv - nu;
                  asmf[nsmf] =  3;                           /* x( is) */

                  nsmf++;
                  ismf[nsmf] = irow - 2*nuv;
                  asmf[nsmf] = -3;                           /* x( ij) */

                  nsmf++;
                  ismf[nsmf] = irow         - nu;
                  asmf[nsmf] =  1;                           /* dv(is) */

                  asmf[irow] =  2;                           /* dv(ij) */

                  nsmf++;
                  ismf[nsmf] = irow -   nuv + nu;
                  asmf[nsmf] =  0;                           /* periodic=2 */
              } else {
                  nsmf++;
                  ismf[nsmf] = irow - 2*nuv - nu;
                  asmf[nsmf] =  3;                           /* x( is) */

                  nsmf++;
                  ismf[nsmf] = irow - 2*nuv + nu;
                  asmf[nsmf] = -3;                           /* x( in) */

                  nsmf++;
                  ismf[nsmf] = irow         - nu;
                  asmf[nsmf] =  1;                           /* dv(is) */

                  asmf[irow] =  4;                           /* dv(ij) */

                  nsmf++;
                  ismf[nsmf] = irow         + nu;
                  asmf[nsmf] =  1;                           /* dv(in) */
              }

              ismf[irow+1] = nsmf + 1;
          }
      }

    /*
     *...M9 and M10 implement Ferguson spline for dv(u) -> duv
     *   (set up required matrix data and indices)
     */
    for (j = 0; j < nv; j++) {
        for (i = 0; i < nu; i++) {
            irow = (i) + nu * (j) + 3 * nuv;

            if (i == 0) {
                nsmf++;
                ismf[nsmf] = irow - nuv;
                asmf[nsmf] =  3;                            /* dv( ij) */

                nsmf++;
                ismf[nsmf] = irow - nuv + 1;
                asmf[nsmf] = -3;                            /* dv( ie) */

                nsmf++;
                ismf[nsmf] = irow + nu  - 2;
                asmf[nsmf] =  0;                            /* periodic=1 */

                asmf[irow] =  2;                            /* duv(ij) */

                nsmf++;
                ismf[nsmf] = irow       + 1;
                asmf[nsmf] =  1;                            /* duv(ie) */
            } else if (i == nu-1) {
                nsmf++;
                ismf[nsmf] = irow - nuv - 1;
                asmf[nsmf] =  3;                            /* dv( iw) */

                nsmf++;
                ismf[nsmf] = irow - nuv;
                asmf[nsmf] = -3;                            /* dv( ij) */

                nsmf++;
                ismf[nsmf] = irow       - 1;
                asmf[nsmf] =  1;                            /* duv(iw) */

                asmf[irow] =  2;                            /* duv(ij) */

                nsmf++;
                ismf[nsmf] = irow - nu  + 1;
                asmf[nsmf] =  0;                            /* periodic=1 */
            } else {
                nsmf++;
                ismf[nsmf] = irow - nuv - 1;
                asmf[nsmf] =  3;                            /* dv( iw) */

                nsmf++;
                ismf[nsmf] = irow - nuv + 1;
                asmf[nsmf] = -3;                            /* dv( ie) */

                nsmf++;
                ismf[nsmf] = irow       - 1;
                asmf[nsmf] =  1;                            /* duv(iw) */

                asmf[irow] =  4;                            /* duv(ij) */

                nsmf++;
                ismf[nsmf] = irow       + 1;
                asmf[nsmf] =  1;                            /* duv(ie) */
            }

            ismf[irow+1] = nsmf + 1;
        }
    }

    /*
     * blocks M1, M2, M3, and M4  (minimize the SSE for the Vertices)
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uu = uv[ivrt].u;
        vv = uv[ivrt].v;

        iw = (int)(uu * (nu - 1));
        iw = MIN(iw, nu-2);

        is = (int)(vv * (nv - 1));
        is = MIN(is, nv-2);

        isw = (iw  ) + nu * (is  );
        ise = (iw+1) + nu * (is  );
        inw = (iw  ) + nu * (is+1);
        ine = (iw+1) + nu * (is+1);

        ss = uu * (nu-1) - iw;
        As = 1 + ss * ( 0 + ss * (-3 + ss * ( 2)));
        Bs = 0 + ss * ( 0 + ss * ( 3 + ss * (-2)));
        Cs = 0 + ss * ( 1 + ss * (-2 + ss * ( 1)));
        Ds = 0 + ss * ( 0 + ss * (-1 + ss * ( 1)));

        tt = vv * (nv-1) - is;
        At = 1 + tt * ( 0 + tt * (-3 + tt * ( 2)));
        Bt = 0 + tt * ( 0 + tt * ( 3 + tt * (-2)));
        Ct = 0 + tt * ( 1 + tt * (-2 + tt * ( 1)));
        Dt = 0 + tt * ( 0 + tt * (-1 + tt * ( 1)));

        /*
         *       [ As*At  Bs*At  As*Bt  Bs*Bt ]
         *       [                            ]
         *       [ Cs*At  Ds*At  Cs*Bt  Ds*Bt ]
         *   b = [                            ]
         *       [ As*Ct  Bs*Ct  As*Dt  Bs*Dt ]
         *       [                            ]
         *       [ Cs*Ct  Ds*Ct  Cs*Dt  Ds*Dt ]
         */
        b00 = As * At;
        b01 = Bs * At;
        b02 = As * Bt;
        b03 = Bs * Bt;

        b10 = Cs * At;
        b11 = Ds * At;
        b12 = Cs * Bt;
        b13 = Ds * Bt;

        b20 = As * Ct;
        b21 = Bs * Ct;
        b22 = As * Dt;
        b23 = Bs * Dt;

        b30 = Cs * Ct;
        b31 = Ds * Ct;
        b32 = Cs * Dt;
        b33 = Ds * Dt;

        /*
         * diagonal elements
         */
        asmf[isw] += b00 * b00;                            /* x(  isw) */
        asmf[ise] += b01 * b01;                            /* x(  ise) */
        asmf[inw] += b02 * b02;                            /* x(  inw) */
        asmf[ine] += b03 * b03;                            /* x(  ine) */

        /*
         * off-diagonal elements (diagonals are commented out)
         */
        for (jsw = ismf[isw]; jsw < ismf[isw+1]; jsw++) {
//          if        (ismf[jsw] == isw      ) {asmf[jsw] += b00 * b00;  /* x(  isw) */
            if        (ismf[jsw] == ise      ) {asmf[jsw] += b00 * b01;  /* x(  ise) */
            } else if (ismf[jsw] == inw      ) {asmf[jsw] += b00 * b02;  /* x(  inw) */
            } else if (ismf[jsw] == ine      ) {asmf[jsw] += b00 * b03;  /* x(  ine) */
            } else if (ismf[jsw] == isw+  nuv) {asmf[jsw] += b00 * b10;  /* du( isw) */
            } else if (ismf[jsw] == ise+  nuv) {asmf[jsw] += b00 * b11;  /* du( ise) */
            } else if (ismf[jsw] == inw+  nuv) {asmf[jsw] += b00 * b12;  /* du( inw) */
            } else if (ismf[jsw] == ine+  nuv) {asmf[jsw] += b00 * b13;  /* du( ine) */
            } else if (ismf[jsw] == isw+2*nuv) {asmf[jsw] += b00 * b20;  /* dv( isw) */
            } else if (ismf[jsw] == ise+2*nuv) {asmf[jsw] += b00 * b21;  /* dv( ise) */
            } else if (ismf[jsw] == inw+2*nuv) {asmf[jsw] += b00 * b22;  /* dv( inw) */
            } else if (ismf[jsw] == ine+2*nuv) {asmf[jsw] += b00 * b23;  /* dv( ine) */
            } else if (ismf[jsw] == isw+3*nuv) {asmf[jsw] += b00 * b30;  /* duv(isw) */
            } else if (ismf[jsw] == ise+3*nuv) {asmf[jsw] += b00 * b31;  /* duv(ise) */
            } else if (ismf[jsw] == inw+3*nuv) {asmf[jsw] += b00 * b32;  /* duv(inw) */
            } else if (ismf[jsw] == ine+3*nuv) {asmf[jsw] += b00 * b33;  /* duv(ine) */
            }
         }

        for (jse = ismf[ise]; jse < ismf[ise+1]; jse++) {
            if        (ismf[jse] == isw      ) {asmf[jse] += b01 * b00;  /* x(  isw) */
//          } else if (ismf[jse] == ise      ) {asmf[jse] += b01 * b01;  /* x(  ise) */
            } else if (ismf[jse] == inw      ) {asmf[jse] += b01 * b02;  /* x(  inw) */
            } else if (ismf[jse] == ine      ) {asmf[jse] += b01 * b03;  /* x(  ine) */
            } else if (ismf[jse] == isw+  nuv) {asmf[jse] += b01 * b10;  /* du( isw) */
            } else if (ismf[jse] == ise+  nuv) {asmf[jse] += b01 * b11;  /* du( ise) */
            } else if (ismf[jse] == inw+  nuv) {asmf[jse] += b01 * b12;  /* du( inw) */
            } else if (ismf[jse] == ine+  nuv) {asmf[jse] += b01 * b13;  /* du( ine) */
            } else if (ismf[jse] == isw+2*nuv) {asmf[jse] += b01 * b20;  /* dv( isw) */
            } else if (ismf[jse] == ise+2*nuv) {asmf[jse] += b01 * b21;  /* dv( ise) */
            } else if (ismf[jse] == inw+2*nuv) {asmf[jse] += b01 * b22;  /* dv( inw) */
            } else if (ismf[jse] == ine+2*nuv) {asmf[jse] += b01 * b23;  /* dv( ine) */
            } else if (ismf[jse] == isw+3*nuv) {asmf[jse] += b01 * b30;  /* duv(isw) */
            } else if (ismf[jse] == ise+3*nuv) {asmf[jse] += b01 * b31;  /* duv(ise) */
            } else if (ismf[jse] == inw+3*nuv) {asmf[jse] += b01 * b32;  /* duv(inw) */
            } else if (ismf[jse] == ine+3*nuv) {asmf[jse] += b01 * b33;  /* duv(ine) */
            }
        }

        for (jnw = ismf[inw]; jnw <ismf[inw+1]; jnw++) {
            if        (ismf[jnw] == isw      ) {asmf[jnw] += b02 * b00;  /* x(  isw) */
            } else if (ismf[jnw] == ise      ) {asmf[jnw] += b02 * b01;  /* x(  ise) */
//          } else if (ismf[jnw] == inw      ) {asmf[jnw] += b02 * b02;  /* x(  inw) */
            } else if (ismf[jnw] == ine      ) {asmf[jnw] += b02 * b03;  /* x(  ine) */
            } else if (ismf[jnw] == isw+  nuv) {asmf[jnw] += b02 * b10;  /* du( isw) */
            } else if (ismf[jnw] == ise+  nuv) {asmf[jnw] += b02 * b11;  /* du( ise) */
            } else if (ismf[jnw] == inw+  nuv) {asmf[jnw] += b02 * b12;  /* du( inw) */
            } else if (ismf[jnw] == ine+  nuv) {asmf[jnw] += b02 * b13;  /* du( ine) */
            } else if (ismf[jnw] == isw+2*nuv) {asmf[jnw] += b02 * b20;  /* dv( isw) */
            } else if (ismf[jnw] == ise+2*nuv) {asmf[jnw] += b02 * b21;  /* dv( ise) */
            } else if (ismf[jnw] == inw+2*nuv) {asmf[jnw] += b02 * b22;  /* dv( inw) */
            } else if (ismf[jnw] == ine+2*nuv) {asmf[jnw] += b02 * b23;  /* dv( ine) */
            } else if (ismf[jnw] == isw+3*nuv) {asmf[jnw] += b02 * b30;  /* duv(isw) */
            } else if (ismf[jnw] == ise+3*nuv) {asmf[jnw] += b02 * b31;  /* duv(ise) */
            } else if (ismf[jnw] == inw+3*nuv) {asmf[jnw] += b02 * b32;  /* duv(inw) */
            } else if (ismf[jnw] == ine+3*nuv) {asmf[jnw] += b02 * b33;  /* duv(ine) */
            }
         }

        for (jne = ismf[ine]; jne < ismf[ine+1]; jne++) {
            if        (ismf[jne] == isw      ) {asmf[jne] += b03 * b00;  /* x(  isw) */
            } else if (ismf[jne] == ise      ) {asmf[jne] += b03 * b01;  /* x(  ise) */
            } else if (ismf[jne] == inw      ) {asmf[jne] += b03 * b02;  /* x(  inw) */
//          } else if (ismf[jne] == ine      ) {asmf[jne] += b03 * b03;  /* x(  ine) */
            } else if (ismf[jne] == isw+  nuv) {asmf[jne] += b03 * b10;  /* du( isw) */
            } else if (ismf[jne] == ise+  nuv) {asmf[jne] += b03 * b11;  /* du( ise) */
            } else if (ismf[jne] == inw+  nuv) {asmf[jne] += b03 * b12;  /* du( inw) */
            } else if (ismf[jne] == ine+  nuv) {asmf[jne] += b03 * b13;  /* du( ine) */
            } else if (ismf[jne] == isw+2*nuv) {asmf[jne] += b03 * b20;  /* dv( isw) */
            } else if (ismf[jne] == ise+2*nuv) {asmf[jne] += b03 * b21;  /* dv( ise) */
            } else if (ismf[jne] == inw+2*nuv) {asmf[jne] += b03 * b22;  /* dv( inw) */
            } else if (ismf[jne] == ine+2*nuv) {asmf[jne] += b03 * b23;  /* dv( ine) */
            } else if (ismf[jne] == isw+3*nuv) {asmf[jne] += b03 * b30;  /* duv(isw) */
            } else if (ismf[jne] == ise+3*nuv) {asmf[jne] += b03 * b31;  /* duv(ise) */
            } else if (ismf[jne] == inw+3*nuv) {asmf[jne] += b03 * b32;  /* duv(inw) */
            } else if (ismf[jne] == ine+3*nuv) {asmf[jne] += b03 * b33;  /* duv(ine) */
            }
         }

        /*
         * right-hand sides
         */
        rhs[isw] += b00 * xvrt[ivrt];
        rhs[ise] += b01 * xvrt[ivrt];
        rhs[inw] += b02 * xvrt[ivrt];
        rhs[ine] += b03 * xvrt[ivrt];

        /*
         * initial guesses
         */
        xx[isw] = xvrt[ivrt];
        xx[ise] = xvrt[ivrt];
        xx[inw] = xvrt[ivrt];
        xx[ine] = xvrt[ivrt];
    }

    /*
     * add smoothing
     */
    for (j = 0; j < nv; j++) {
        for (i = 1; i < nu-1; i++) {
            ij = (i) + nu * (j);

            asmf[ij] += 2 * amu;                          /* x(ij) */

            for (jj = ismf[ij]; jj < ismf[ij+1]; jj++) {
                if        (ismf[jj] == ij-1    ) {
                    asmf[jj] -= amu;                      /* x(iw) */
                } else if (ismf[jj] == ij+1    ) {
                    asmf[jj] -= amu;                      /* x(ie) */
                }
            }
        }
    }

    for (j = 1; j < nv-1; j++) {
        for (i = 0; i < nu; i++) {
            ij = (i) + nu * (j);

            asmf[ij] += 2 * amu;                          /* x(ij) */

            for (jj = ismf[ij]; jj < ismf[ij+1]; jj++) {
                if        (ismf[jj] == ij-nu   ) {
                    asmf[jj] -= amu;                      /* x(is) */
                } else if (ismf[jj] == ij+nu   ) {
                    asmf[jj] -= amu;                      /* x(in) */
                }
            }
        }
    }

    /*
     * modify Ferguson splines if periodic in U
     */
    if (periodic == 1 && ppnts != NULL) {

        /* blocks M1 to M4 */
        for (j = 0; j < nv; j++) {
            vvv   = (double)(j) / (double)(nv - 1);
            ileft = 1;
            irite = ppnts[0];
            while (irite > ileft+1) {
                imidl = (ileft + irite) / 2;
                if (vvv > uv[ppnts[imidl]].v) {
                    ileft = imidl;
                } else {
                    irite = imidl;
                }
            }

            vleft = uv[  ppnts[ileft]].v;
            vrite = uv[  ppnts[irite]].v;
            xleft = xvrt[ppnts[ileft]];
            xrite = xvrt[ppnts[irite]];
            xxx   = xleft + (vvv - vleft) * (xrite - xleft) / (vrite - vleft);

            for (i = 0; i < nu; i += nu-1) {
                irow = (i) + nu * (j);
                asmf[irow] = 1;
                rhs[ irow] = xxx;

                for (jj = ismf[irow]; jj < ismf[irow+1]; jj++) {
                    asmf[jj] = 0;
                }
            }
        }

        /* blocks M5 and M6 */
        for (j = 0; j < nv; j++) {
            irow = (0) + nu * (j) + nuv;
            asmf[irow] = 4;

            for (jj = ismf[irow]; jj < ismf[irow+1]; jj++) {
                if        (ismf[jj] == irow-nuv ) {
                    ismf[jj] += nu - 2;
                } else if (ismf[jj] == irow+nu-2) {
                    asmf[jj] = 1;
                }
            }

            for (jj = ismf[irow+nu-1]; jj < ismf[irow+nu]; jj++) {
                if (ismf[jj] == irow) {
                    asmf[jj] = -2;
                } else {
                    asmf[jj] =  0;
                }
            }
        }

        /* blocks M9 and M10 */
        for (j = 0; j < nv; j++) {
            irow = (0) + nu * (j) + 3*nuv;
            asmf[irow] = 4;

            for (jj = ismf[irow]; jj < ismf[irow+1]; jj++) {
                if        (ismf[jj] == irow-nuv ) {
                    ismf[jj] += nu - 2;
                } else if (ismf[jj] == irow+nu-2) {
                    asmf[jj] = 1;
                }
            }

            for (jj = ismf[irow+nu-1]; jj < ismf[irow+nu]; jj++) {
                if (ismf[jj] == irow) {
                    asmf[jj] = -2;
                } else {
                    asmf[jj] =  0;
                }
            }
        }

    /*
     * modify Ferguson splines if periodic in V
     */
    } else if (periodic == 2 && ppnts != NULL) {

        /* blocks M1 to M4 */
        for (i = 0; i < nu; i++) {
            uuu  = (double)(i) / (double)(nu - 1);
            ileft = 1;
            irite = ppnts[0];
            while (irite > ileft+1) {
                imidl = (ileft + irite) / 2;
                if (uuu > uv[ppnts[imidl]].u) {
                    ileft = imidl;
                } else {
                    irite = imidl;
                }
            }

            uleft = uv[  ppnts[ileft]].u;
            urite = uv[  ppnts[irite]].u;
            xleft = xvrt[ppnts[ileft]];
            xrite = xvrt[ppnts[irite]];
            xxx   = xleft + (uuu - uleft) * (xrite - xleft) / (urite - uleft);

            for (j = 0; j < nv; j += nv-1) {
                irow = (i) + nu * (j);
                asmf[irow] = 1;
                rhs[ irow] = xxx;

                for (jj = ismf[irow]; jj < ismf[irow+1]; jj++) {
                    asmf[jj] = 0;
                }

            }
        }

        /* blocks M7 and M8 */
        for (i = 0; i < nu; i++) {
            irow = i + nu * (0) + 2*nuv;
            asmf[irow] = 4;

            for (jj = ismf[irow]; jj < ismf[irow+1]; jj++) {
                if (ismf[jj] == irow-2*nuv) {
                    ismf[jj] += nuv -2*nu;
                } else if (ismf[jj] == irow+nuv-2*nu) {
                    asmf[jj] = 1;
                }
            }

            for (jj = ismf[irow+nuv-nu]; jj < ismf[irow+nuv-nu+1]; jj++) {
                if (ismf[jj] == irow) {
                    asmf[jj] = -2;
                } else {
                    asmf[jj] =  0;
                }
            }
        }
    }

//$$$    /*
//$$$     * make sure that we do not have a row of all zeros
//$$$     */
//$$$    for (j = 0; j < nv; j++) {
//$$$        for (i = 0; i < nu; i++) {
//$$$            ij = (i) + nu * (j);
//$$$
//$$$            sum = asmf[ij];
//$$$            for (jj = ismf[ij]; jj < ismf[ij+1]; jj++) {
//$$$                sum += asmf[jj];
//$$$            }
//$$$
//$$$            if (sum == 0) {
//$$$                asmf[ij] = 1;
//$$$                rhs[ ij] = 0;
//$$$            }
//$$$        }
//$$$    }

    /*
     * perturb initial guess so that we do not have a lot of zeroes
     */
    for (i = 0; i < 4*nuv; i++) {
        if (fabs(xx[i]) < xtol) {
            xx[i] = xtol;
        }
    }

    /*
     * pin the values (to zero) at Knots that have no adjacent Vertices
     */
    for (i = 0; i < nuv; i++) {
        if (nn[i] == 0) {
            DPRINT1("pinning %d", i);

            asmf[i] = 1;
            rhs[ i] = 0;

            for (jj = ismf[i]; jj < ismf[i+1]; jj++) {
                asmf[jj] = 0;
            }
        }
    }

    /*
     * solve for the new x-locations
     */
    status = sparseBCG(asmf, ismf, xx, rhs);
    DPRINT1("sparseBCG -> status=%d", status);
    CHECK_STATUS;

    /*
     * extract the new x values from the matrix solution
     */
    for (i = 0; i < nuv; i++) {
        x[ i] = xx[i      ];
        du[i] = xx[i+  nuv];
        dv[i] = xx[i+2*nuv];
        dc[i] = xx[i+3*nuv];
    }

#ifdef DEBUG
    DPRINT0("grid generated by splineWithGaps2D");
    DPRINT0x("i=   ");
    for (i = 0; i < nu; i++) {
        DPRINT1x(" %5d     ", i);
    }
    DPRINT0("");

    for (j = 0; j < nv; j++) {
        DPRINT1x("j=%3d", j);
        for (i = 0; i < nu; i++) {
            ij = (i) + nu * (j);
            DPRINT1x(" %10.5f", xx[ij]);
        }
        DPRINT0("");
    }
#endif

 cleanup:
    FREE(xx  );
    FREE(rhs );
    FREE(ismf);
    FREE(asmf);
    FREE(nn  );

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_LimitGridSize -- set a global Grid size limit                            *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_LimitGridSize(int      limit)            /* (in)   global size limit */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    ROUTINE(prm_LimitGridSize);
    DPRINT2("%s(limit=%d) {",
            routine, limit);

    /* ----------------------------------------------------------------------- */

    globalSizeLimit = limit;

// cleanup:
    DPRINT2("%s --> status=%d}",
            routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_FixedGrid -- set up a fixed Grid for a set of Vertices                   *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_FixedGrid(int      nvrt,                 /* (in)   number of Vertices */
              int      nvar,                 /* (in)   number of dependent vars */
              prmUV    uv[],                 /* (in)   array  of Parameters */
              double   var[],                /* (in)   array  of dependent vars */
              int      ntri,                 /* (in)   number of Triangles (optional) */
  /*@null@*/  prmTri   tri[],                /* (in)   array  of Triangles (optional) */
              int      periodic,             /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        = 2  periodic in V */
   /*@null@*/ int      ppnts[],              /* (in)   indices of periodic points */
              int      nu,                   /* (in)   number of Knots in U-dirn */
              int      nv,                   /* (in)   number of Knots in V-dirn */
              double   grid[],               /* (out)  grid of Knots */
              double   *rmserr,              /* (out)  RMS     error at Vertices */
              double   *maxerr,              /* (out)  maximum error at Vertices */
              double   *dotmin)              /* (out)  minimum dot product */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      *resid0 = NULL;
    double      *ddu    = NULL;
    double      *ddv    = NULL;

    gridTree    tree;
    int         iu, iv, iuv, ii, jj, dtype, ivar, nsize;
    double      xyz[100], dxyzdu[100], dxyzdv[100];
    prmUV       uuvv;

    ROUTINE(prm_FixedGrid);
    DPRINT6("%s(nvrt=%d, nvar=%d, periodic=%d, nu=%d, nv=%d) {",
            routine, nvrt, nvar, periodic, nu, nv);
    for (iv = 0; iv < nvrt; iv++) {
        DPRINT4("%5d  %10.5f %10.5f  %10.5f",
                iv, uv[iv].u, uv[iv].v, var[nvar*iv]);
    }

    /* ----------------------------------------------------------------------- */

    PPRINT4("   prm_FixedGrid(nvrt=%d, nvar=%d, ntri=%d, periodic=%d)",
            nvrt, nvar, ntri, periodic);

    *dotmin = 1;

    if (nvrt <= 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nvar <  0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nu <=  1) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nv <=  1) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

#ifdef GRAFIC2
    plotVrts(nvrt, nvar, uv, var);
#endif

    /*
     * make sure nu-1 and nv-1 are both powers of 2
     */
    ii = nu - 1;
    while (1) {
        if (ii == 1) {
            break;
        } else if (ii%2 == 0) {
            ii /= 2;
        } else {
            status = PRM_BADPARAM;
            goto cleanup;
        }
    }

    jj = nu - 1;
    while (1) {
        if (jj == 1) {
            break;
        } else if (jj%2 == 0) {
            jj /= 2;
        } else {
            status = PRM_BADPARAM;
            goto cleanup;
        }
    }

    /*
     * the initial residual is just the original data
     */
    nsize = nvrt * nvar * sizeof(double);
    MALLOC(resid0, double, nsize);
    MEMCPY(resid0, var,    nsize);

    /*
     * initialize the Tree and set up the first (global) Cell
     */
    status = initialTree2d(&tree, nvar, nvrt, periodic, ppnts,
                           uv, resid0, rmserr, maxerr);
    CHECK_STATUS;

    PPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree.nu, tree.nv, *rmserr, *maxerr);
    DPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree.nu, tree.nv, *rmserr, *maxerr);

    /*
     * keep refining until desired number of Knots exist
     */
    while (tree.nu < nu || tree.nv < nv) {
                          dtype =  0;
        if (tree.nu < nu) dtype += 1;
        if (tree.nv < nv) dtype += 2;

        status = globalRefine2d(dtype, &tree, nvrt, periodic, ppnts,
                                uv, resid0, rmserr, maxerr);
        CHECK_STATUS;

        PPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
               tree.nu, tree.nv, *rmserr, *maxerr);
        DPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
               tree.nu, tree.nv, *rmserr, *maxerr);
    }

    MALLOC(ddu, double, nvar*nu*nv);
    MALLOC(ddv, double, nvar*nu*nv);

    /*
     * evaluate Tree at every Knot so that ferguson array can be returned
     */
    for (iv = 0; iv < nv; iv++) {
        for (iu = 0; iu < nu; iu++) {
            uuvv.u  = (double)(iu) / (double)(nu - 1);
            uuvv.v  = (double)(iv) / (double)(nv - 1);
            iuv = (iu) + nu * (iv);

            prm_EvalGrid(tree, uuvv, xyz, dxyzdu, dxyzdv, NULL, NULL, NULL);

            for (ivar = 0; ivar < nvar; ivar++) {
                grid[nvar*iuv+ivar] = xyz[   ivar];
                ddu[ nvar*iuv+ivar] = dxyzdu[ivar];
                ddv[ nvar*iuv+ivar] = dxyzdv[ivar];
            }
        }
    }

#ifdef GRAFIC
    {
        int io_kbd=5, io_scr=6;
        grinit_(&io_kbd, &io_scr, "at end of prm_FixedGrid",
                           strlen("at end of prm_FixedGrid"));
        plotGrid(nu, nv, nvar, nvrt, uv, grid, ddu, ddv);
    }
#endif

    /*
     * check the quality of the Grid (if ntri > 0)
     */
    if (ntri > 0 && tri != NULL && nvar == 3) {
        *dotmin = computeDotmin(&tree, nvar, uv, var, ntri, tri);
    }

 cleanup:
    FREE(ddv   );
    FREE(ddu   );
    FREE(resid0);

    DPRINT4("%s --> rmserr=%f, maxerr=%f, status=%d}",
            routine, *rmserr, *maxerr, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_BestGrid -- set up the best (coarsest) Grid for a set of Vertices        *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_BestGrid(int      nvrt,                  /* (in)   number of Vertices */
             int      nvar,                  /* (in)   number of dependent vars */
             prmUV    uv[],                  /* (in)   array  of Parameters */
             double   var[],                 /* (in)   array  of dependent vars */
             int      ntri,                  /* (in)   number of Triangles (optional) */
  /*@null@*/ prmTri   tri[],                 /* (in)   array  of Triangles (optional) */
             double   tol,                   /* (in)   tolerance on maximum error */
             int      periodic,              /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        = 2  periodic in V */
  /*@null@*/ int      ppnts[],               /* (in)   indices of periodic points */
             int      *nu,                   /* (in)   limit on nu if > 0 */
                                             /* (out)  number of Knots in U-dirn */
             int      *nv,                   /* (in)   limit on nv if > 0 */
                                             /* (out)  number of Knots in V-dirn */
             double   *grid[],               /* (out)  pointer to grid of Knots */
             double   *rmserr,               /* (out)  RMS     error at Vertices */
             double   *maxerr,               /* (out)  maximum error at Vertices */
             double   *dotmin)               /* (out)  minimum dot product */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      *resid0 = NULL;
    double      *resid1 = NULL;
    double      *resid2 = NULL;
    double      *ddu    = NULL;
    double      *ddv    = NULL;

    gridTree    tree0;
    gridTree    tree1;
    gridTree    tree2;

    int         iu, iv, iuv, numax, nvmax, nuvmax, ivrt, ivar, nsize, i;
    int         ivrt_max, ivar_max;
    double      rmserr1, rmserr2, maxerr1, maxerr2, resid_max;
    double      xyz[100], dxyzdu[100], dxyzdv[100];
    prmUV       uuvv;
    char        tag[2];

    ROUTINE(prm_BestGrid);
    DPRINT7("%s(nvrt=%d, nvar=%d, tol=%f, periodic=%d, nu=%d, nv=%d) {",
            routine, nvrt, nvar, tol, periodic, *nu, *nv);
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        strcpy(tag, " ");
        if (ppnts != NULL) {
            for (i = 1; i <= ppnts[0]; i++) {
                if (ppnts[i] == ivrt) strcpy(tag, "*");
            }
        }
        DPRINT4x("%5d %s %10.5f %10.5f  ",
                 ivrt, tag, uv[ivrt].u, uv[ivrt].v);
        for (ivar = 0; ivar < nvar; ivar++) {
            DPRINT1x("%10.5f ", var[ivar+nvar*ivrt]);
        }
        DPRINT0("");
    }

    /* ----------------------------------------------------------------------- */

    PPRINT5("   prm_BestGrid(nvrt=%d, nvar=%d, ntri=%d, tol=%f, periodic=%d)",
            nvrt, nvar, ntri, tol, periodic);

    *dotmin = 1;

    if (nvrt <= 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nvar <  0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (tol <= 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

#ifdef GRAFIC2
    plotVrts(nvrt, nvar, uv, var);
#endif

#ifdef DEBUG
    /*
     * print periodic points (in order)
     */
    if (ppnts != NULL) {
        DPRINT1("periodic=%d", periodic);
        for (i = 1; i <= ppnts[0]; i++) {
            ivrt = ppnts[i];
            DPRINT7("ppnts[%3d]=%5d   %10.5f %10.5f   %10.5f %10.5f %10.5f",
                    i, ivrt, uv[ivrt].u, uv[ivrt].v,
                    var[nvar*ivrt], var[nvar*ivrt+1], var[nvar*ivrt+2]);
        }
    }
#endif

    /*
     * if nu and nv are given positives, use them as max alloable size
     */
    if (*nu > 0 && *nv == 0) {
        numax  = globalSizeLimit;
        nvmax  = globalSizeLimit;
        nuvmax = *nu;
    } else {
        if (*nu > 0) {
            numax = *nu;
        } else {
            numax = globalSizeLimit;
        }

        if (*nv > 0) {
            nvmax = *nv;
        } else {
            nvmax = globalSizeLimit;
        }
        nuvmax = numax * nvmax;
    }

    /*
     * initialize the pointers in the Trees
     */
    tree0.cell = NULL;
    tree0.knot = NULL;

    tree1.cell = NULL;
    tree1.knot = NULL;

    tree2.cell = NULL;
    tree2.knot = NULL;

    /*
     * the initial residual is just the original data
     */
    nsize = nvrt * nvar * sizeof(double);
    MALLOC(resid0, double, nsize);
    MALLOC(resid1, double, nsize);
    MALLOC(resid2, double, nsize);

    MEMCPY(resid0, var,    nsize);

    /*
     * initialize the Tree and set up the first (global) Cell(s)
     */
    status = initialTree2d(&tree0, nvar, nvrt, periodic, ppnts,
                           uv, resid0, rmserr, maxerr);
    CHECK_STATUS;

    PPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, tree0.nv, *rmserr, *maxerr);
    DPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, tree0.nv, *rmserr, *maxerr);

    /*
     * keep refining until desired tolerance is met (or max size is reached)
     */
    while (*maxerr > tol && (tree0.nu)*(tree0.nv) < nuvmax) {

        /*
         * try refinement type 1
         */
        if (tree0.nu < numax) {
            copyGrid(&tree1, &tree0);
            MEMCPY(resid1, resid0, nsize);

            status = globalRefine2d(1, &tree1, nvrt, periodic, ppnts,
                                    uv, resid1, &rmserr1, &maxerr1);
            CHECK_STATUS;

            PPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
                    tree1.nu, tree1.nv, rmserr1, maxerr1);
            DPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
                    tree1.nu, tree1.nv, rmserr1, maxerr1);
        } else {
            rmserr1 = HUGEQ;
            maxerr1 = HUGEQ;
        }

        /*
         * try refinement type 2
         */
        if (tree0.nv < nvmax) {
            copyGrid(&tree2, &tree0);
            MEMCPY(resid2, resid0, nsize);

            status = globalRefine2d(2, &tree2, nvrt, periodic, ppnts,
                                    uv, resid2, &rmserr2, &maxerr2);
            CHECK_STATUS;

            PPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
                    tree2.nu, tree2.nv, rmserr2, maxerr2);
            DPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
                    tree2.nu, tree2.nv, rmserr2, maxerr2);
        } else {
            rmserr2 = HUGEQ;
            maxerr2 = HUGEQ;
        }

        /*
         * if neither refinement is better than what we already have, stop
         *
         * note: the 10 here makes sure that we try at least a few levels
         *       of refinement since the errors might climb a little
         *       in the beginning
         */
        if (MIN(maxerr1, maxerr2) > *maxerr &&
            MIN(rmserr1, rmserr2) > *rmserr &&
            tree0.nu*tree0.nv     > 10        ) {

            status = PRM_TOLERANCEUNMET;
            break;

        /*
         * keep the better of refinement 1 or 2
         */
        } else if (rmserr1 < rmserr2) {
            copyGrid(&tree0, &tree1);
            MEMCPY(resid0, resid1, nsize);
            *rmserr = rmserr1;
            *maxerr = maxerr1;

        } else {
            copyGrid(&tree0, &tree2);
            MEMCPY(resid0, resid2, nsize);
            *rmserr = rmserr2;
            *maxerr = maxerr2;
        }

        /*
         * print the maximum residual
         */
        ivrt_max  = -1;
        ivar_max  = -1;
        resid_max =  0;

        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            for (ivar = 0; ivar < nvar; ivar++) {
                if (fabs(resid0[ivrt+ivar*nvrt]) > resid_max) {
                    ivrt_max  = ivrt;
                    ivar_max  = ivar;
                    resid_max = fabs(resid0[ivrt+ivar*nvrt]);
                }
            }
        }

        DPRINT4("(%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
               tree0.nu, tree0.nv, *rmserr, *maxerr);
        DPRINT2("Maximum resid @ ivrt=%d, ivar=%d: ", ivrt_max, ivar_max);
        DPRINT1("                u = %10.5f", uv[ivrt_max].u);
        DPRINT1("                v = %10.5f", uv[ivrt_max].v);
        for (ivar = 0; ivar < nvar; ivar++) {
            DPRINT3("                %d = %10.5f  resid = %10.3e",
                    ivar, var[ivrt_max+ivar*nvrt], resid0[ivrt_max+ivar*nvrt]);
        }
    }

    PPRINT0("   (---------)");
    PPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, tree0.nv, *rmserr, *maxerr);
    DPRINT4("   (%4d,%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, tree0.nv, *rmserr, *maxerr);

    /*
     * malloc the Grid
     */
    *nu = tree0.nu;
    *nv = tree0.nv;
    MALLOC(*grid, double, nvar*(*nu)*(*nv));
    MALLOC( ddu,  double, nvar*(*nu)*(*nv));
    MALLOC( ddv,  double, nvar*(*nu)*(*nv));

    /*
     * evaluate Tree at every Knot so that ferguson array can be returned
     */
    for (iv = 0; iv < *nv; iv++) {
        for (iu = 0; iu < *nu; iu++) {
            uuvv.u  = (double)(iu) / (double)(*nu - 1);
            uuvv.v  = (double)(iv) / (double)(*nv - 1);
            iuv = (iu) + (*nu) * (iv);

            prm_EvalGrid(tree0, uuvv, xyz, dxyzdu, dxyzdv, NULL, NULL, NULL);

            for (ivar = 0; ivar < nvar; ivar++) {
                (*grid)[nvar*iuv+ivar] = xyz[   ivar];
                ddu[    nvar*iuv+ivar] = dxyzdu[ivar];
                ddv[    nvar*iuv+ivar] = dxyzdv[ivar];
            }
        }
    }

#ifdef GRAFIC
    {
        int io_kbd=5, io_scr=6;
        grinit_(&io_kbd, &io_scr, "plotGrid at end of prm_BestGrid",
                           strlen("plotGrid at end of prm_BestGrid"));
        plotGrid(*nu, *nv, nvar, nvrt, uv, *grid, ddu, ddv);
    }
#endif
#ifdef GRAFIC
    {
        int io_kbd=5, io_scr=6;
        grinit_(&io_kbd, &io_scr, "plotLevels at end of prm_BestGrid",
                           strlen("plotLevels at end of prm_BestGrid"));
        plotLevels(&tree0);
    }
#endif

    /*
     * check the quality of the Grid (if ntri > 0)
     */
    if (ntri > 0 && tri != NULL && nvar == 3) {
        *dotmin = computeDotmin(&tree0, nvar, uv, var, ntri, tri);
    }

 cleanup:
    prm_FreeGrid(&tree2);
    prm_FreeGrid(&tree1);
    prm_FreeGrid(&tree0);

    FREE(ddv   );
    FREE(ddu   );
    FREE(resid2);
    FREE(resid1);
    FREE(resid0);

    DPRINT4("%s --> rmserr=%f, maxerr=%f, status=%d}",
            routine, *rmserr, *maxerr, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_EvalGrid -- evaluate the Grid at the given uu                            *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_EvalGrid(gridTree tree,                  /* (in)   Grid Tree */
             prmUV    uv,                    /* (in)   independent variables */
             double   var[],                 /* (out)  dependent variables */
  /*@null@*/ double   *du,                   /* (out)  optional first u-derivatives */
  /*@null@*/ double   *dv,                   /* (out)  optional first v-derivatives */
  /*@null@*/ double   *duu,                  /* (out)  optional second u-derivatives */
  /*@null@*/ double   *duv,                  /* (out)  optional cross derivatives */
  /*@null@*/ double   *dvv)                  /* (out)  optional second v-derivatives */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    double      *xyz = NULL;
    double      *Du  = NULL;
    double      *Dv  = NULL;
    double      *Duu = NULL;
    double      *Duv = NULL;
    double      *Dvv = NULL;

    int         ivar, icel;
    double      ss, tt, umin, umax, vmin, vmax;

    int         nvar  =      tree.nvar;
    int         nvar4 =  4 * tree.nvar;

    ROUTINE(prm_EvalGrid);
//    DPRINT3("%s(uu=%f, vv=%f) {",
//            routine, uv.u, uv.v);

    /* ----------------------------------------------------------------------- */

    MALLOC(xyz, double, nvar);
    MALLOC(Du,  double, nvar);
    MALLOC(Dv,  double, nvar);
    MALLOC(Duu, double, nvar);
    MALLOC(Duv, double, nvar);
    MALLOC(Dvv, double, nvar);

    /*
     * initialize the dependent variables
     */
    for (ivar = 0; ivar < nvar; ivar++) {
                         var[ivar] = 0;
        if (du  != NULL) du[ ivar] = 0;
        if (dv  != NULL) dv[ ivar] = 0;
        if (duu != NULL) duu[ivar] = 0;
        if (duv != NULL) duv[ivar] = 0;
        if (dvv != NULL) dvv[ivar] = 0;
    }

    /*
     * start evaluation at the root of the Grid Tree
     */
    icel = 0;

    /*
     * recursively add the appropriate evaluations
     */
    while (1) {
        /*
         * add in current Cell
         */
        umin = tree.cell[icel].umin;
        umax = tree.cell[icel].umax;
        vmin = tree.cell[icel].vmin;
        vmax = tree.cell[icel].vmax;

        ss =  (uv.u - umin) / (umax - umin);
        tt =  (uv.v - vmin) / (vmax - vmin);

        if (du == NULL || dv == NULL) {
            evalBicubic(ss, tt, nvar,
                        &(tree.knot[nvar4*(tree.cell[icel].swKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].seKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].nwKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].neKnot)]),
                        xyz, NULL, NULL, NULL, NULL, NULL);

            for (ivar = 0; ivar < nvar; ivar++) {
                var[ivar] += xyz[ivar];
            }
        } else if (duu == NULL || duv == NULL || dvv == NULL) {
            evalBicubic(ss, tt, nvar,
                        &(tree.knot[nvar4*(tree.cell[icel].swKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].seKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].nwKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].neKnot)]),
                        xyz, Du, Dv, NULL, NULL, NULL);

            for (ivar = 0; ivar < nvar; ivar++) {
                var[ivar] += xyz[ivar];
                du[ ivar] += Du[ ivar] / (umax - umin);
                dv[ ivar] += Dv[ ivar] / (vmax - vmin);
            }
        } else {
            evalBicubic(ss, tt, nvar,
                        &(tree.knot[nvar4*(tree.cell[icel].swKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].seKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].nwKnot)]),
                        &(tree.knot[nvar4*(tree.cell[icel].neKnot)]),
                        xyz, Du, Dv, Duu, Duv, Dvv);

            for (ivar = 0; ivar < nvar; ivar++) {
                var[ivar] += xyz[ivar];
                du[ ivar] += Du[ ivar] / (umax - umin);
                dv[ ivar] += Dv[ ivar] / (vmax - vmin);
                duu[ivar] += Duu[ivar] / (umax - umin) / (umax - umin);
                duv[ivar] += Duv[ivar] / (umax - umin) / (vmax - vmin);
                dvv[ivar] += Dvv[ivar] / (vmax - vmin) / (vmax - vmin);
            }
        }

        /*
         * if at bottom of hierarchy, stop
         */
        if (tree.cell[icel].dtype == 0) {
            break;

        /*
         * otherwise go to the appropriate child
         */
        } else if (tree.cell[icel].dtype == 1) {
            if (uv.u <= (umin+umax)/2) {
                icel  = tree.cell[icel].child;
            } else {
                icel  = tree.cell[icel].child + 1;
            }

        } else if (tree.cell[icel].dtype == 2) {
            if (uv.v <= (vmin+vmax)/2) {
                icel  = tree.cell[icel].child;
            } else {
                icel  = tree.cell[icel].child + 1;
            }

        } else if (uv.u <= (umin+umax)/2) {
            if (uv.v <= (vmin+vmax)/2) {
                icel  = tree.cell[icel].child;
            } else {
                icel  = tree.cell[icel].child + 2;
            }
        } else {
            if (uv.v <= (vmin+vmax)/2) {
                icel  = tree.cell[icel].child + 1;
            } else {
                icel  = tree.cell[icel].child + 3;
            }
        }
    }

 cleanup:
    FREE(Dvv);
    FREE(Duv);
    FREE(Duu);
    FREE(Dv );
    FREE(Du );
    FREE(xyz);
//    for (ivar = 0; ivar < nvar; ivar++) {
//        DPRINT1x("%12.6f ", var[ivar]);
//    }
//    DPRINT0("");

//    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_FreeGrid -- free Cells and Knots associated with Grid Tree               *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_FreeGrid(gridTree *tree)                 /* (in)   pointer to Grid Tree */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    ROUTINE(prm_FreeGrid);
    DPRINT1("%s() {",
            routine);

    /* ----------------------------------------------------------------------- */

    /*
     * delete the Cell and Knot arrays
     */
#ifndef __clang_analyzer__
    FREE(tree->cell);
    FREE(tree->knot);
#endif

    /*
     * reset Tree's information
     */
    DPRINT0("reseting Tree information");
    tree->nvar = 0;
    tree->nu   = 0;
    tree->nv   = 0;
    tree->ncel = 0;
    tree->nknt = 0;

// cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}
