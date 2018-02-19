/*
 ******************************************************************************
 * EGADS: Electronic Geometry Aircraft Design System                          *
 *                                                                            *
 * This module creates and manages curve fits (Cfits)                         *
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
static int    copyCfit                (cfitTree*, cfitTree);
static int    divideCell              (cfitTree*, int);
static void   evalCubic               (double, int, double[], double[], double[],
                                       /*@null@*/double[], /*@null@*/double[]);
static int    finestCell              (cfitTree*, double);
static int    globalRefine1d          (cfitTree*, int, int, double[], double[],
                                       double*, double*);
static int    initialTree1d           (cfitTree*, int, int, int, double[],
                                       double[], double*, double*);
#ifdef GRAFIC
static void   plotCfit                (int, int, double[], double[], int,
                                       double[], char[]);
#endif
#ifdef DEBUG
static void   printTree               (FILE*, cfitTree);
#endif
static int    splineWithGaps1d        (int, double[], double[], int, int,
                                       double[], double[]);

/*
 * external routines defined in prmUV
 */
extern int    sparseBCG               (double[], int[], double[], double[]);

/*
 * global constants
 */
#define  HUGEQ           1.0e+40
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
           {if (dbg_fp == NULL) dbg_fp = fopen("prmCfit.dbg", "w");}
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
 * copyCfit -- copy one Cfit Tree to another                                    *
 *                                                                              *
 ********************************************************************************
 */
static int
copyCfit(cfitTree *dest,                     /* (out)  pointer to destination Tree */
         cfitTree src)                       /* (in)   source Tree */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */

    ROUTINE(copyCfit);
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
    dest->nvar = src.nvar;
    dest->nu   = src.nu;
    dest->ncel = src.ncel;
    dest->nknt = src.nknt;

    /*
     * allocate new arrays
     */
    MALLOC(dest->cell, cfitCell,  dest->ncel                );
    MALLOC(dest->knot, double,   (dest->nknt)*2*(dest->nvar));

    memcpy(dest->cell, src.cell, (dest->ncel)               *sizeof(cfitCell));
    memcpy(dest->knot, src.knot, (dest->nknt)*2*(dest->nvar)*sizeof(double  ));

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
divideCell(cfitTree *tree,                   /* (in)   pointer to Cfit Tree */
           int      icel)                    /* (in)   Cell to divide */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_BADPARAM */
                                             /*        PRM_BADDIVISION */

    int         nvar2 = 2 * tree->nvar;

    int         iknt, nknt, ncel, ivar, jcel;
    int         iw, ie, knotw, knotc, knote, nborw, nbore;
    double      umin, umax;

    ROUTINE(divideCell);
    DPRINT2("%s(icel=%d) {",
            routine, icel);

    /* ----------------------------------------------------------------------- */

    /*
     * make sure that Cell is divisible
     */
    if (icel < 0 || icel >= tree->ncel) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (tree->cell[icel].child > 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

    ncel = tree->ncel;
    nknt = tree->nknt;

    knotw = -1;
    knotc = -1;
    knote = -1;

    jcel = tree->cell[icel].wCell;
    if (jcel <= 0) {
        nborw = jcel;
    } else if (tree->cell[jcel].child == 0) {
        nborw = 0;
    } else {
        nborw = tree->cell[jcel ].child + 1;
        knotw = tree->cell[nborw].eKnot;
    }

    jcel = tree->cell[icel].eCell;
    if (jcel <= 0) {
        nbore = jcel;
    } else if (tree->cell[jcel].child == 0) {
        nbore = 0;
    } else {
        nbore = tree->cell[jcel ].child;
        knote = tree->cell[nbore].wKnot;
    }

    /*
     * number the new Knots (if any)
     */
    if (knotw < 0) knotw = nknt++;
    if (knotc < 0) knotc = nknt++;
    if (knote < 0) knote = nknt++;

    /*
     * names of children Cells
     */
    iw = ncel;
    ie = ncel + 1;

    /*
     * extend the Cell array for the new Cells
     */
    RALLOC(tree->cell, cfitCell, tree->ncel+2);

    umin = tree->cell[icel].umin;
    umax = tree->cell[icel].umax;

    /*
     * make the children Cell
     */
    tree->cell[ncel].wKnot  = knotw;
    tree->cell[ncel].eKnot  = knotc;
    tree->cell[ncel].wCell  = nborw;
    tree->cell[ncel].eCell  = ie;
    tree->cell[ncel].child  = 0;
    tree->cell[ncel].count  = 0;
    tree->cell[ncel].umin   =  umin;
    tree->cell[ncel].umax   = (umin + umax) / 2;
    tree->cell[ncel].rmserr = 0;
    tree->cell[ncel].maxerr = 0;
    ncel++;

    tree->cell[ncel].wKnot  = knotc;
    tree->cell[ncel].eKnot  = knote;
    tree->cell[ncel].wCell  = iw;
    tree->cell[ncel].eCell  = nbore;
    tree->cell[ncel].child  = 0;
    tree->cell[ncel].count  = 0;
    tree->cell[ncel].umin   = (umin + umax) / 2;
    tree->cell[ncel].umax   =         umax;
    tree->cell[ncel].rmserr = 0;
    tree->cell[ncel].maxerr = 0;
    ncel++;

    /*
     * inform the neighbors of the new Cells
     */
    if (nborw > 0) tree->cell[nborw].eCell = iw;
    if (nbore > 0) tree->cell[nbore].wCell = ie;

    /*
     * inform icel of its new children
     */
    tree->cell[icel].child = iw;

    /*
     * extend the Knot array
     */
    RALLOC(tree->knot, double, nvar2*nknt);

    for (iknt = tree->nknt; iknt < nknt; iknt++) {
        for (ivar = 0; ivar < nvar2; ivar++) {
            tree->knot[nvar2*iknt+ivar] = 0;
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
 * evalCubic -- evaluate a cubic                                                *
 *                                                                              *
 ********************************************************************************
 */
static void
evalCubic(double   ss,                       /* (in)   fractional distance */
          int      nvar,                     /* (in)   number of variables */
          double   w[],                      /* (in)   Knot on west */
          double   e[],                      /* (in)   Knot on east */
          double   xyz[],                    /* (out)  dependent variables */
/*@null@*/double   dt[],                     /* (out)  optional d/dT */
/*@null@*/double   dtt[])                    /* (out)  optional d2/dT2 */
{
    double      a, b, c, d;
    int         ivar;

//    ROUTINE(evalCubic);
//    DPRINT2("%s(ss=%f) {",
//             routine, ss);

    /* ----------------------------------------------------------------------- */

    for (ivar = 0; ivar < nvar; ivar++) {
        a =       w[ivar];
        b =                                   w[nvar+ivar];
        c = - 3 * w[ivar] + 3 * e[ivar] - 2 * w[nvar+ivar] - e[nvar+ivar];
        d =   2 * w[ivar] - 2 * e[ivar] +     w[nvar+ivar] + e[nvar+ivar];

        xyz[ivar] = a + ss * (b + ss * (c + ss * d));

        if (dt != NULL) {
            dt[ivar] = b + ss * (2 * c + ss  * 3 * d);

            if (dtt != NULL) {
                dtt[ivar] = 2 * c + ss * 6 * d;
            }
        }
    }

// cleanup:
//    DPRINT1("%s --> }", routine);
}



/*
 ********************************************************************************
 *                                                                              *
 * finestCell -- find the finest Cell containing (u)                            *
 *                                                                              *
 ********************************************************************************
 */
static int
finestCell(cfitTree *tree,                   /* (in)   Cfit Tree */
           double   u)                       /* (in)   parameter */
{
    int         icel = 0;                    /* (out)  finest Cell containing u */

    double      umin, umax;

//    ROUTINE(finestCell);
//    DPRINT2("%s(u=%f) {",
//            routine, u);

    /* ----------------------------------------------------------------------- */

    /*
     * traverse the Tree until we find an undivided Cell
     */
    while (1) {
        umin = tree->cell[icel].umin;
        umax = tree->cell[icel].umax;

        if (tree->cell[icel].child == 0) {
            break;
        } else {
            if (u <= (umin+umax)/2) {
                icel  = tree->cell[icel].child;
            } else {
                icel  = tree->cell[icel].child + 1;
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
 * globalRefine1d -- refine a Tree globally                                     *
 *                                                                              *
 ********************************************************************************
 */
extern int
globalRefine1d(cfitTree *tree,               /* (both) Cfit Tree to refine */
               int      nvrt,                /* (in)   number of Vertices */
               int      periodic,            /* (in)   periodicity flag */
               double   u[],                 /* (in)   array  of Vertices */
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
    double      *xyz    = NULL;

    int         ibeg, iend, ii, ivrt, ivar, icel, iknt, nu;
    int         ncel;
    double      umin, umax, uu, ss, err;

    int         nvar  =     tree->nvar;
    int         nvar2 = 2 * tree->nvar;

    ROUTINE(globalRefine1d);
    DPRINT3("%s(nvrt=%d, periodic=%d) {",
            routine, nvrt, periodic);

    /* ----------------------------------------------------------------------- */

    *rmserr = 0;
    *maxerr = 0;

    MALLOC(xyz, double, nvar);

    /*
     * find the first and last Cell with child==0
     */
    ibeg = +99999;
    iend = -99999;

    for (icel = 0; icel < tree->ncel; icel++) {
        if (tree->cell[icel].child == 0) {
            ibeg = MIN(ibeg, icel);
            iend = MAX(iend, icel);
        }
    }

    /*
     * adjust number of Knots in Tree
     */
    tree->nu = 2 * tree->nu - 1;

    /* ncel           first new Cell */
    /* tree->ncel-1   last  new Cell */
    ncel = tree->ncel;

    /*
     * create the new Cells and new Knots
     */
    for (icel = ibeg; icel <= iend; icel++) {
        status = divideCell(tree, icel);
        CHECK_STATUS;
    }

    /*
     * count number of Vertices in each new Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        icel = finestCell(tree, u[ivrt]);
        tree->cell[icel].count++;
    }

#ifdef DEBUG
    printTree(dbg_fp, *tree);
#endif

    /*
     * allocate temp storage for the Vertices and the extra data
     */
    nu = tree->nu;
    DPRINT1("nu=%d", tree->nu);

    MALLOC(xvrt,  double, nvrt);
    MALLOC(xspln, double, nu  );
    MALLOC(uspln, double, nu  );

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
        status = splineWithGaps1d(nvrt, u, xvrt, periodic,
                                  tree->nu, xspln, uspln);
        CHECK_STATUS;

        /*
         * store the data in the Knots
         */
        for (icel = ncel; icel < tree->ncel; icel++) {
            ii = tree->nu * tree->cell[icel].umin;

            iknt = tree->cell[icel].wKnot;
            tree->knot[nvar2*iknt      +ivar] = xspln[(ii  )];
            tree->knot[nvar2*iknt+nvar +ivar] = uspln[(ii  )];

            iknt = tree->cell[icel].eKnot;
            tree->knot[nvar2*iknt      +ivar] = xspln[(ii+1)];
            tree->knot[nvar2*iknt+nvar +ivar] = uspln[(ii+1)];
        }
    }

    /*
     * subtract the latest spline from "resid" and keep
     *    track of the max error associated with each Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uu = u[ivrt];

        for (icel = ncel; icel < tree->ncel; icel++) {
            umin = tree->cell[icel].umin;
            umax = tree->cell[icel].umax;

            if (uu >= umin && uu <= umax) {
                ss  = (uu - umin) / (umax - umin);
                evalCubic(ss, nvar,
                          &(tree->knot[nvar2*(tree->cell[icel].wKnot)]),
                          &(tree->knot[nvar2*(tree->cell[icel].eKnot)]),
                          xyz, NULL, NULL);

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
 * initialTree1d -- set up initial Tree and global Cell                         *
 *                                                                              *
 ********************************************************************************
 */
extern int
initialTree1d(cfitTree *tree,                /* (both) Cfit Tree to refine */
              int      nvar,                 /* (in)   number of variables */
              int      nvrt,                 /* (in)   number of Vertices */
              int      periodic,             /* (in)   periodicity flag */
              double   u[],                  /* (in)   array  of Vertices */
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
    double      *xyz    = NULL;

    int         ncel, ivrt, ivar, iknt, icel, ii;
    double      err, uu, umin, umax, ss;

    int         nvar2 = 2 * nvar;

    ROUTINE(initialTree1d);
    DPRINT3("%s(nvrt=%d, periodic=%d) {",
            routine, nvrt, periodic);

    /* ----------------------------------------------------------------------- */

    *rmserr = 0;
    *maxerr = 0;

    MALLOC(xyz, double, nvar);

    /*
     * initialize the Tree
     */
    tree->nvar = nvar;
    tree->nu   = 2;
    tree->ncel = 1;
    tree->cell = NULL;
    tree->nknt = 2;
    tree->knot = NULL;

    /*
     * allocate the initialize the Cell and Knot arrays
     */
    MALLOC(tree->cell, cfitCell, tree->ncel      );
    MALLOC(tree->knot, double,   tree->nknt*nvar2);

    for (iknt = 0; iknt < tree->nknt; iknt++) {
        for (ivar = 0; ivar < nvar2; ivar++) {
            tree->knot[nvar2*iknt+ivar] = 0;
        }
    }

    /*
     * set up root (global) Cell
     */
    ncel = 0;
    tree->cell[ncel].wKnot  =  0;
    tree->cell[ncel].eKnot  =  1;
    tree->cell[ncel].wCell  = -1;
    tree->cell[ncel].eCell  = -1;
    tree->cell[ncel].child  =  0;
    tree->cell[ncel].count  =  0;

    tree->cell[ncel].umin   =  0;
    tree->cell[ncel].umax   =  1;
    tree->cell[ncel].rmserr =  0;
    tree->cell[ncel].maxerr =  0;
    ncel++;

    /*
     * if periodic, then refine the global Cell (twice).
     *
     *  note: ncel will be the first "fine" Cell
     */
    if (periodic != 1) {
        tree->nu = 2;
        ncel     = 0;
    } else {
        tree->nu = 5;
        ncel     = 3;
        for (icel = 0; icel < ncel; icel++) {
            status = divideCell(tree, icel);
            CHECK_STATUS;
        }
    }

    /*
     * count number of Vertices in each Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        icel = finestCell(tree, u[ivrt]);
        tree->cell[icel].count++;
    }

    /*
     * allocate storage for (extended) Vertex table and spline data
     */
    MALLOC(xvrt,  double,       nvrt);
    MALLOC(xspln, double, tree->nknt);
    MALLOC(uspln, double, tree->nknt);

#ifdef DEBUG
    printTree(dbg_fp, *tree);
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
        status = splineWithGaps1d(nvrt, u, xvrt, periodic,
                                  tree->nu, xspln, uspln);
        CHECK_STATUS;

        /*
         * store the spline data in the knots
         */
        for (icel = ncel; icel < tree->ncel; icel++) {
            ii = tree->nu * tree->cell[icel].umin;

            iknt = tree->cell[icel].wKnot;
            tree->knot[nvar2*iknt     +ivar] = xspln[ii  ];
            tree->knot[nvar2*iknt+nvar+ivar] = uspln[ii  ];

            iknt = tree->cell[icel].eKnot;
            tree->knot[nvar2*iknt     +ivar] = xspln[ii+1];
            tree->knot[nvar2*iknt+nvar+ivar] = uspln[ii+1];
        }
    }

    /*
     * subtract the current spline from "resid" and keep
     *    track of the max error associated with the global Cell
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uu = u[ivrt];

        icel = finestCell(tree, uu);
        umin = tree->cell[icel].umin;
        umax = tree->cell[icel].umax;
        ss   = (uu - umin) / (umax - umin);

        evalCubic(ss, nvar,
                  &(tree->knot[nvar2*(tree->cell[icel].wKnot)]),
                  &(tree->knot[nvar2*(tree->cell[icel].eKnot)]),
                  xyz, NULL, NULL);

        for (ivar = 0; ivar < nvar; ivar++) {
            resid[nvar*ivrt+ivar] -= xyz[ivar];
            err = fabs(resid[nvar*ivrt+ivar]);

            tree->cell[icel].rmserr += SQR(err);
            *rmserr                 += SQR(err);
            *maxerr = MAX(*maxerr,         err);

            if (err > tree->cell[icel].maxerr) tree->cell[icel].maxerr = err;
        }
    }

    *rmserr = sqrt(*rmserr / nvrt);

 cleanup:
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
 * prm_plotCfit -- plot the U coordinates of Vertices (using GRAFIC)            *
 *                                                                              *
 ********************************************************************************
 */
static void
plotCfit(int      nvrt,                      /* (in)   number of Vertices */
         int      nvar,                      /* (in)   number of Variables */
         double   u[],                       /* (in)   array  of Parameters */
         double   xyz[],                     /* (in)   array  of Coordinates */
         int      nu,                        /* (in)   number of Knots */
         double   cfit[],                    /* (in)   array  of Knots */
         char     pltitl[])                  /* (in)   plot title */
{
    int         io_kbd = 5, io_scr = 6, indgr = 1 + 4 + 16 + 64;
    int         nline, nplot, ilin[2*nvar], isym[2*nvar], nper[2*nvar];
    int         ivar, i;
    float       xplot[(nvrt+nu)*nvar], yplot[(nvrt+nu)*nvar];

    ROUTINE(plotCfit);
    DPRINT4("%s(nvrt=%d, nvar=%d, nu=%d) {",
            routine, nvrt, nvar, nu);

    /* ----------------------------------------------------------------------- */

    nplot = 0;
    nline = 0;

    /*
     * plot the (input) Vertices
     */
    for (ivar = 0; ivar < nvar; ivar++) {
        ilin[nline] = - (1 + ivar %  5);
        isym[nline] = + (1 + ivar % 15);
        nper[nline] = nvrt;
        nline++;

        for (i = 0; i < nvrt; i++) {
            xplot[nplot] = u[       i     ];
            yplot[nplot] = xyz[nvar*i+ivar];
            nplot++;
        }
    }

    /*
     * plot the (output) Knots
     */
    for (ivar = 0; ivar < nvar; ivar++) {
        ilin[nline] = + (1 + ivar %  5);
        isym[nline] = - (1 + ivar % 15);
        nper[nline] = nu;
        nline++;

        for (i = 0; i < nu; i++) {
            xplot[nplot] = (double)(i) / (double)(nu - 1);
            yplot[nplot] = cfit[nvar*i+ivar];
            nplot++;
        }
    }

    grinit_(&io_kbd, &io_scr, "plotCfit (X=sqr, Y=circ, Z=tri)",
                       strlen("plotCfit (X=sqr, Y=circ, Z=tri)"));
    grline_(ilin, isym, &nline, pltitl, &indgr, xplot, yplot, nper, strlen(pltitl));

// cleanup:
    DPRINT1("%s --> }", routine);
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
          cfitTree tree)                     /* (in)   Cfit Tree */
{

    int         nvar  =      tree.nvar;
    int         nvar2 =  2 * tree.nvar;

    int         icel, iknt, ivar;

    ROUTINE(printTree);
    DPRINT1("%s() {",
            routine);

    /* ----------------------------------------------------------------------- */

    fprintf(fp, " icel wKnot eKnot wCell eCell child count     umin    umax     maxerr\n");
    for (icel = 0; icel < tree.ncel; icel++) {
        fprintf(fp, "%5d %5d %5d %5d %5d %5d %5d %8.5f%8.5f %12.7f\n", icel,
                tree.cell[icel].wKnot,
                tree.cell[icel].eKnot,
                tree.cell[icel].wCell,
                tree.cell[icel].eCell,
                tree.cell[icel].child,
                tree.cell[icel].count,
                tree.cell[icel].umin,
                tree.cell[icel].umax,
                tree.cell[icel].maxerr);
    }

    fprintf(fp, " iknt\n");
    for (iknt = 0; iknt < tree.nknt; iknt++) {
        fprintf(fp, "%5d ", iknt);
        for (ivar = 0; ivar < nvar; ivar++) {
            fprintf(fp, "%10.5f", tree.knot[nvar2*iknt      +ivar]);
        }
        fprintf(fp, "  ");
        for (ivar = 0; ivar < nvar; ivar++) {
            fprintf(fp, "%10.5f", tree.knot[nvar2*iknt+nvar +ivar]);
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
 * splineWithGaps1d -- create Ferguson spline that least-square fits 1D data    *
 *                                                                              *
 ********************************************************************************
 */
static int
splineWithGaps1d(int      nvrt,              /* (in)   number of Vertices to fit */
                 double   uvrt[],            /* (in)   U for Vertices to fit */
                 double   xvrt[],            /* (in)   X (data) to fit */
                 int      periodic,          /* (in)   periodicity flag */
                 int      nu,                /* (in)   number of U grid Knots */
                 double   x[],               /* (out)  array of spline Knots */
                 double   d[])               /* (out)  array of Knot derivs */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */

    double      *asmf = NULL;                /* sparse-matrix data */
    int         *ismf = NULL;                /* sparse-matrix indices */
    double      *rhs  = NULL;                /* right-hand sides */
    double      *xx   = NULL;                /* matrix solution: x, du */
    double      *xmin = NULL;                /* minimum X near each Knot */
    double      *xmax = NULL;                /* maximum X near each Knot */
    int         *nn   = NULL;                /* number of Vertices near each Knot */

    double      uu, ss, As, Bs, Cs, Ds;
    double      b00, b01, b10, b11;
    double      umin, umax, xavg, xumin;

    int         i, ii, jj, kk, ivrt, iu, ipin, npin;
    int         iw, ie, jw, je;
    int         nsmf, irow, icol;

    int         maxpin = 10;                 /* maximum pinning passes */
    double      amu    = 0.001;              /* smoothing coefficient */
    double      xtol   = 1.e-6;              /* tolerance on x for constant value */

    ROUTINE(splineWithGaps1d);
    DPRINT4("%s(nvrt=%d, periodic=%d, nu=%d) {",
            routine, nvrt, periodic, nu);
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        DPRINT3("uvrt[%3d]=%10.6f, xvrt[]=%10.6f",
                ivrt, uvrt[ivrt], xvrt[ivrt]);
    }

    /* ----------------------------------------------------------------------- */

    /*
     * if periodic and nu <= 3, then the spline is simply a constant
     */
    if (periodic == 1 && nu <= 3) {
        xavg = 0;
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            xavg += xvrt[ivrt];
        }
        xavg /= nvrt;

        for (i = 0; i < nu; i++) {
            x[i] = xavg;
            d[i] = 0;
        }

        goto cleanup;
    }

    /*
     * find the extrema of the X-data in each cell
     */
    MALLOC(xmin, double, nu);
    MALLOC(xmax, double, nu);
    MALLOC(nn,   int,    nu);

    for (i = 0; i < nu; i++) {
        xmin[i] = +HUGEQ;
        xmax[i] = -HUGEQ;
        nn[  i] = 0;
    }

    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        iu = MIN(MAX(0, (nu-1)*uvrt[ivrt]), nu-2);

                                     nn[  iu  ]++;
        if (xvrt[ivrt] < xmin[iu  ]) xmin[iu  ] = xvrt[ivrt];
        if (xvrt[ivrt] > xmax[iu  ]) xmax[iu  ] = xvrt[ivrt];

                                     nn[  iu+1]++;
        if (xvrt[ivrt] < xmin[iu+1]) xmin[iu+1] = xvrt[ivrt];
        if (xvrt[ivrt] > xmax[iu+1]) xmax[iu+1] = xvrt[ivrt];
    }

    /*
     * overall matrix equation is of form:
     *
     *     [      |      ] [     ]   [     ]
     *     [  M1  |  M2  ] [  x  ]   [ RHS ]   (minimize SSE)
     *     [      |      ] [     ]   [     ]
     *     [ - - - - - --] [ --- ] = [ --- ]
     *     [      |      ] [     ]   [     ]
     *     [  M5  |  M6  ] [ du  ]   [  0  ]   (spline x(u))
     *     [      |      ] [     ]   [     ]
     *
     *     set up the storage in sparse-matrix form: 11 = 3 + 3
     *                                                  + 2 + 3
     */

    MALLOC(asmf, double, 12*nu); /* safety because of extra element */
    MALLOC(ismf, int,    12*nu);
    MALLOC(rhs,  double,  2*nu);
    MALLOC(xx,   double,  2*nu);

    /*
     * start out with all zeros
     */
    for (i = 0; i < 12*nu; i++) {
        asmf[i] = 0;
        ismf[i] = 0;
    }

    for (i = 0; i < 2*nu; i++) {
        rhs[i] = 0;
        xx[ i] = 0;
    }

    /*
     * set up sparse-matrix index matrix as well as sparse-matrix
     *    data for the bottom blocks.  note that some of the
     *    off-diagonal matrix elements might turn out to be zero (which
     *    wastes a little space, but make the code below much simpler)
     */
    ismf[0] = 2 * nu + 1;
    nsmf    = 2 * nu;

    /*
     * ...blocks M1 and M2 implement least-squares fit -> x
     *    (set up indices only - data filled in below)
     */
    for (i = 0; i < nu; i++) {
        irow = i;

        for (kk = 0; kk < 2; kk++) {
            for (ii = MAX(0, i-1); ii <= MIN(nu-1, i+1); ii++) {
                icol = (ii) + nu * (kk);

                if (irow != icol) {
                    nsmf++;
                    ismf[nsmf] = icol;
                    asmf[nsmf] = 0;
                }
            }
        }

        ismf[irow+1] = nsmf + 1;
    }

    /*
     * ...M5 and M6 implement Ferguson spline for x(u) -> du
     *    (set up required matrix data and indices)
     */
    for (i = 0; i < nu; i++) {
        irow = (i) + nu;

        if (i == 0) {
            nsmf++;
            ismf[nsmf] = irow - nu;
            asmf[nsmf] =  3;                                 /* x( i  ) */

            nsmf++;
            ismf[nsmf] = irow - nu + 1;
            asmf[nsmf] = -3;                                 /* x( i+1) */

            asmf[irow] =  2;                                 /* du(i  ) */

            nsmf++;
            ismf[nsmf] = irow      + 1;
            asmf[nsmf] =  1;                                 /* du(i+1) */

            nsmf++;
            ismf[nsmf] = 2 * nu - 2;
            asmf[nsmf] = 0;                                  /* du(n-2) */
        } else if (i == nu-1) {
            nsmf++;
            ismf[nsmf] = irow - nu - 1;
            asmf[nsmf] =  3;                                 /* x( i-1) */

            nsmf++;
            ismf[nsmf] = irow - nu;
            asmf[nsmf] = -3;                                 /* x( i  ) */

            nsmf++;
            ismf[nsmf] = irow      - 1;
            asmf[nsmf] =  1;                                 /* du(i-1) */

            asmf[irow] =  2;                                 /* du(i  ) */
        } else {
            nsmf++;
            ismf[nsmf] = irow - nu - 1;
            asmf[nsmf] =  3;                                 /* x( i-1) */

            nsmf++;
            ismf[nsmf] = irow - nu + 1;
            asmf[nsmf] = -3;                                 /* x( i+1) */

            nsmf++;
            ismf[nsmf] = irow      - 1;
            asmf[nsmf] =  1;                                 /* du(i-1) */

            asmf[irow] =  4;                                 /* du(i  ) */

            nsmf++;
            ismf[nsmf] = irow      + 1;
            asmf[nsmf] =  1;                                 /* du(i+1) */
        }

        ismf[irow+1] = nsmf + 1;
    }

    /*
     * blocks M1 and M2   (minimize the SSE for the Vertices)
     */
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        uu = uvrt[ivrt];
        iw = (int)(uu * (nu - 1));
        iw = MIN(iw, nu-2);
        ie = iw + 1;

        ss = uu * (nu-1) - iw;
        As = 1 + ss * (   + ss * (-3 + ss * ( 2)));
        Bs =   + ss * (   + ss * ( 3 + ss * (-2)));
        Cs =   + ss * ( 1 + ss * (-2 + ss       ));
        Ds =   + ss * (   + ss * (-1 + ss       ));

        /*
         *       [ As  Bs ]
         *   b = [        ]
         *       [ Cs  Ds ]
         */
        b00 = As;
        b01 = Bs;
        b10 = Cs;
        b11 = Ds;

        /*
         * diagonal elements
         */
        asmf[iw] += b00 * b00;                               /* x(  iw) */
        asmf[ie] += b01 * b01;                               /* x(  ie) */

        /*
         * off-diagonal elements
         */
        for (jw = ismf[iw]; jw < ismf[iw+1]; jw++) {
            if (       ismf[jw] == ie   ) {
                asmf[jw] += b00 * b01;                       /* x(  ie) */
            } else if (ismf[jw] == iw+nu) {
                asmf[jw] += b00 * b10;                       /* du( iw) */
            } else if (ismf[jw] == ie+nu) {
                asmf[jw] += b00 * b11;                       /* du( ie) */
            }
        }

        for (je = ismf[ie]; je < ismf[ie+1]; je++) {
            if (       ismf[je] == iw   ) {
                asmf[je] += b01 * b00;                       /* x(  iw) */
            } else if (ismf[je] == iw+nu) {
                asmf[je] += b01 * b10;                       /* du( iw) */
            } else if (ismf[je] == ie+nu) {
                asmf[je] += b01 * b11;                       /* du( ie) */
            }
        }

        /*
         * right-hand sides
         */
        rhs[iw] = rhs[iw] + b00 * xvrt[ivrt];
        rhs[ie] = rhs[ie] + b01 * xvrt[ivrt];

        /*
         * initial guesses
         */
        xx[iw] = xvrt[ivrt];
        xx[ie] = xvrt[ivrt];
    }

    /*
     * modify Ferguson splines if periodic
     */
    if (periodic == 1) {
        umin  = uvrt[0];
        xumin = xvrt[0];
        for (ivrt = 1; ivrt < nvrt; ivrt++) {
            if (uvrt[ivrt] < umin) {
                umin  = uvrt[ivrt];
                xumin = xvrt[ivrt];
            }
        }

        /*
         * row 0
         */
        asmf[0] = 1;                                     /* x(0   ) */
        rhs[ 0] = xumin;

        for (jj = ismf[0]; jj < ismf[1]; jj++) {
            asmf[jj] = 0;                                /* x(*   ) */
        }

        /*
         * row nu-1
         */
        asmf[nu-1] = 1;                                  /* x(nu-1) */
        rhs[ nu-1] = xumin;

        for (jj = ismf[nu-1]; jj < ismf[nu]; jj++) {
            asmf[jj] = 0;                                /* x(*   ) */
        }

        /*
         * row nu
         */
        asmf[nu] = 4;                                    /* du(0   ) */
        rhs[ nu] = 0;

        for (jj = ismf[nu]; jj < ismf[nu+1]; jj++) {
            if        (ismf[jj] == 0     ) {
                ismf[jj] = nu - 2;                       /* x(nu-2) */
            } else if (ismf[jj] == 2*nu-2) {
                asmf[jj] = 1;                            /* du(nu-2) */
            }
        }

        /*
         * row 2*nu-1
         */
        asmf[2*nu-1] = 1;                                /* du(nu-1) */
        rhs[ 2*nu-1] = 0;

        for (jj = ismf[2*nu-1]; jj < ismf[2*nu]; jj++) {
            if (ismf[jj] == 2*nu-2) {
                ismf[jj] = nu;
                asmf[jj] = -1;                           /* du(0   ) */
            } else {
                asmf[jj] =  0;                           /* x(*   ) */
            }
        }
    }

    /*
     * add smoothing
     */
    for (i = 1; i < nu-1; i++) {
        asmf[i] += 2 * amu;                                  /* x(i ) */

        for (jj = ismf[i]; jj < ismf[i+1]; jj++) {
            if        (ismf[jj] == i-1) {
                asmf[jj] -= amu;                             /* x(iw) */
            } else if (ismf[jj] == i+1) {
                asmf[jj] -= amu;                             /* x(ie) */
            }
        }
    }

    /*
     * modify Ferguson splines if value at left end is fixed
     */
    if (periodic == 11 || periodic == 12) {
        asmf[0] = 1;                                     /* x(0   ) */

        for (jj = ismf[0]; jj < ismf[1]; jj++) {
            asmf[jj] = 0;
        }

        umin = +HUGEQ;
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (uvrt[ivrt] >= 0 && uvrt[ivrt] < umin) {
                umin   = uvrt[ivrt];
                rhs[0] = xvrt[ivrt];
            }
        }
    }

    /*
     * modify Ferguson splines if value and slope at left end is zero
     */
    if (periodic == 21 || periodic == 22) {
        asmf[   0] = 1;                                  /* x(0   ) */
        asmf[nu+0] = 1;                                  /* d(0   ) */

        for (jj = ismf[   0]; jj < ismf[   1]; jj++) {
            asmf[jj] = 0;
        }
        for (jj = ismf[nu+0]; jj < ismf[nu+1]; jj++) {
            asmf[jj] = 0;
        }

        rhs[   0] = 0;
        rhs[nu+0] = 0;
    }

    /*
     * modify Ferguson splines if value at rite end is fixed
     */
    if (periodic == 11 || periodic == 21) {
        asmf[nu-1] = 1;                                  /* x(nu-1) */

        for (jj = ismf[nu-1]; jj < ismf[nu]; jj++) {
            asmf[jj] = 0;
        }

        umax = -HUGEQ;
        for (ivrt = 0; ivrt < nvrt; ivrt++) {
            if (uvrt[ivrt] <= 1 && uvrt[ivrt] > umax) {
                umax      = uvrt[ivrt];
                rhs[nu-1] = xvrt[ivrt];
            }
        }
    }

    /*
     * modify Ferguson splines if value and slope at rite end is zero
     */
    if (periodic == 12 || periodic == 22) {
        asmf[  nu-1] = 1;                                /* x(nu-1) */
        asmf[2*nu-1] = 1;                                /* d(nu-1) */

        for (jj = ismf[  nu-1]; jj < ismf[  nu]; jj++) {
            asmf[jj] = 0;
        }
        for (jj = ismf[2*nu-1]; jj < ismf[2*nu]; jj++) {
            asmf[jj] = 0;
        }

        rhs[  nu-1] = 0;
        rhs[2*nu-1] = 0;
    }

    /*
     * perturb initial guess so that we do not have a lot of zeroes
     */
    for (i = 0; i < 2*nu; i++) {
        if (fabs(xx[i]) < xtol) {
            xx[i] = xtol;
        }
    }

    /*
     * make up to maxpin loops to pin spline Knots to the
     *    extrema of the surrounding Vertices (if any)
     */
    ipin = 0;
    while (1) {
        ipin++;

        /*
         * solve for the new x-locations
         */
        status = sparseBCG(asmf, ismf, xx, rhs);
        DPRINT1("sparseBCG -> status=%d", status);
        CHECK_STATUS;

        if (ipin > maxpin) break;

        /*
         * if any of the Knot values exceed the values in the
         *    adjacent Cells, pin it by modifying the equations
         *    and re-solve
         */
        npin = 0;
        for (i = 0; i < nu; i++) {
            if (nn[i] == 0) {

            } else if (xx[i] < xmin[i]) {
                asmf[i] = 1;
                rhs[ i] = xmin[i];

                for (jj = ismf[i]; jj < ismf[i+1]; jj++) {
                    asmf[jj] = 0;
                }

                npin++;
            } else if (xx[i] > xmax[i]) {
                asmf[i] = 1;
                rhs[ i] = xmax[i];

                for (jj = ismf[i]; jj < ismf[i+1]; jj++) {
                    asmf[jj] = 0;
                }

                npin++;
            }
        }

        if (npin == 0) break;
        DPRINT2("ipin=%d  npin=%d", ipin, npin);
    }

    /*
     * extract the new x values from the matrix solution
     */
    for (i = 0; i < nu; i++) {
        x[i] = xx[i   ];
        d[i] = xx[i+nu];
    }

 cleanup:
    FREE(xx  );
    FREE(rhs );
    FREE(ismf);
    FREE(asmf);
    FREE(nn  );
    FREE(xmax);
    FREE(xmin);

    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_FixedCfit -- update the Cfit Knots to match the given Vertices           *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_FixedCfit(int      nvrt,                 /* (in)   number of Vertices */
              int      nvar,                 /* (in)   number of dependent vars */
              double   u[],                  /* (in)   array  of Parameters */
              double   var[],                /* (in)   array  of dependent vars */
              int      periodic,             /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        =11  left fixed, rite fixed */
                                             /*        =12  left fixed, rite zero  */
                                             /*        =21  left zero,  rite fixed */
                                             /*        =22  left zero,  rite zero  */
              int      nu,                   /* (in)   number of Knots */
              double   cfit[],               /* (out)  Cfit of Knots */
              double   *rmserr,              /* (out)  RMS     error at Vertices */
              double   *maxerr)              /* (out)  maximum error at Vertices */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      *resid  = NULL;              /* array of residuals */

    cfitTree    tree;
    int         iu, ii, ivar;
    double      uu, xyz[100];

    ROUTINE(prm_FixedCfit);
    DPRINT5("%s(nvrt=%d, nvar=%d, periodic=%d, nu=%d) {",
            routine, nvrt, nvar, periodic, nu);
    for (iu = 0; iu < nvrt; iu++) {
        DPRINT3("%5d  %10.5f   %10.5f", iu, u[iu], var[nvar*iu]);
    }

    /* ----------------------------------------------------------------------- */

    PPRINT3("   prm_FixedCfit(nvrt=%d, nvar=%d, periodic=%d)",
            nvrt, nvar, periodic);

    if        (nvrt <= 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nvar <  0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nu <= 1) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

    /*
     * make sure nu-1 is a power of 2
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

    /*
     * storage for residuals
     */
    MALLOC(resid, double, nvrt*nvar);

    /*
     * initialize the Tree and set up the first (global) Cell
     */
    memcpy(resid, var, nvrt*nvar*sizeof(double));
    status = initialTree1d(&tree, nvar, nvrt, periodic,
                           u, resid, rmserr, maxerr);
    CHECK_STATUS;

    PPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree.nu, *rmserr, *maxerr);
    DPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree.nu, *rmserr, *maxerr);

    /*
     * keep refining until desired number of Knots exist
     */
    while (tree.nu < nu) {
        status = globalRefine1d(&tree, nvrt, periodic,
                                u, resid, rmserr, maxerr);
        CHECK_STATUS;

        PPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f\n",
                tree.nu, *rmserr, *maxerr);
        DPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f\n",
                tree.nu, *rmserr, *maxerr);
    }

    /*
     * evaluate Tree at every Knot so that ferguson array can be returned
     */
    for (iu = 0; iu < nu; iu++) {
        uu = (double)(iu) / (double)(nu - 1);

        prm_EvalCfit(tree, uu, xyz, NULL, NULL);

        for (ivar = 0; ivar < nvar; ivar++) {
            cfit[nvar*iu+ivar] = xyz[ivar];
        }
    }

#ifdef GRAFIC
    plotCfit(nvrt, nvar, u, var, nu, cfit,
             "~U~var~Vertices and Knots in prm_FixedCfit");
#endif

cleanup:
    FREE(resid);

    DPRINT4("%s --> rmserr=%f, maxerr=%f, status=%d}",
            routine, *rmserr, *maxerr, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_BestCfit -- set up the best (ie, coarsest) Cfit for a set of Vertices    *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_BestCfit(int      nvrt,                  /* (in)   number of Vertices */
             int      nvar,                  /* (in)   number of dependent vars */
             double   u[],                   /* (in)   array  of Parameters */
             double   var[],                 /* (in)   array  of dependent vars */
             double   tol,                   /* (in)   tolerance on maximum error */
             int      periodic,              /* (in)   = 0  no periodicity */
                                             /*        = 1  periodic in U */
                                             /*        =11  left fixed, rite fixed */
                                             /*        =12  left fixed, rite zero  */
                                             /*        =21  left zero,  rite fixed */
                                             /*        =22  left zero,  rite zero  */
             int      *nu,                   /* (in)   limit on nu if > 0 */
                                             /* (out)  number of Knots */
             double   *cfit[],               /* (out)  pointer to Cfit of Knots
                                                       NOTE: user must EG_free after use */
             double   *rmserr,               /* (out)  RMS     error at Vertices */
             double   *maxerr)               /* (out)  maximum error at Vertices */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */
                                             /*        PRM_TOLERANCEUNMET */
                                             /*        PRM_NOTCONVERGED */
                                             /*        PRM_BADPARAM */

    double      *resid0 = NULL;              /* array of residuals assoc w/tree0 */
    double      *resid1 = NULL;              /* array of residuals assoc w/tree1 */

    cfitTree    tree0;
    cfitTree    tree1;

    int         iu, numax, ivrt, ivar, nsize;
    double      rmserr1, maxerr1, xyz[100], uu;

    ROUTINE(prm_BestCfit);
    DPRINT6("%s(nvrt=%d, nvar=%d, tol=%f, periodic=%d, nu=%d) {",
            routine, nvrt, nvar, tol, periodic, *nu);
    for (ivrt = 0; ivrt < nvrt; ivrt++) {
        DPRINT2x("%5d %10.5f  ",
                 ivrt, u[ivrt]);
        for (ivar = 0; ivar < nvar; ivar++) {
            DPRINT1x("%10.5f ", var[ivar+nvar*ivrt]);
        }
        DPRINT0("");
    }

    /* ----------------------------------------------------------------------- */

    PPRINT4("   prm_BestCfit(nvrt=%d, nvar=%d, tol=%f, periodic=%d)",
            nvrt, nvar, tol, periodic);

    if        (nvrt <= 0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (nvar <  0) {
        status = PRM_BADPARAM;
        goto cleanup;
    } else if (tol <=  0) {
        status = PRM_BADPARAM;
        goto cleanup;
    }

    /*
     * if nu is given positive, use it as maximum allowable size
     */
    if (*nu > 0) {
        numax = *nu;
    } else {
        numax = 513;
    }

    /*
     * initialize the pointers in the Trees
     */
    tree0.cell = NULL;
    tree0.knot = NULL;

    tree1.cell = NULL;
    tree1.knot = NULL;

    /*
     * the initial residual is just the original data
     */
    nsize = nvrt * nvar * sizeof(double);
    MALLOC(resid0, double, nsize);
    MALLOC(resid1, double, nsize);

    MEMCPY(resid0, var,    nsize);

    /*
     * initialize the Tree and set up the first (global) Cell(s)
     */
    status = initialTree1d(&tree0, nvar, nvrt, periodic,
                           u, resid0, rmserr, maxerr);
    CHECK_STATUS;

    PPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, *rmserr, *maxerr);
    DPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, *rmserr, *maxerr);

    /*
     * keep refining until desired tolerance is met (or max size is reached)
     */
    while (*maxerr > tol && tree0.nu < numax) {

        /*
         * try refinement
         */
        copyCfit(&tree1, tree0);
        MEMCPY(resid1, resid0, nsize);

        status = globalRefine1d(&tree1, nvrt, periodic,
                                u, resid1, &rmserr1, &maxerr1);
        CHECK_STATUS;

        PPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
                tree1.nu, rmserr1, maxerr1);
        DPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
                tree1.nu, rmserr1, maxerr1);

        /*
         * if refinement is not better than what we already have, stop
         *
         * note: the 5 here makes sure that we at least a few levels
         *       of refinement since the errors might climb a little
         *       in the beginning
         */
        if (maxerr1  > *maxerr &&
            rmserr1  > *rmserr &&
            tree0.nu > 5         ) {

            status = PRM_TOLERANCEUNMET;
            break;

        /*
         * keep the refinement
         */
        } else {
            copyCfit(&tree0, tree1);
            MEMCPY(resid0, resid1, nsize);
            *rmserr = rmserr1;
            *maxerr = maxerr1;
        }
    }

    PPRINT0("   (----)");
    PPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, *rmserr, *maxerr);
    DPRINT3("   (%4d)  rmserr=%12.6f  maxerr=%12.6f",
            tree0.nu, *rmserr, *maxerr);

    /*
     * malloc the Cfit
     */
    *nu = tree0.nu;
    MALLOC(*cfit, double, nvar*(*nu));

    /*
     * evaluate Tree at every Knot so that ferguson array can be returned
     */
    for (iu = 0; iu < *nu; iu++) {
        uu = (double)(iu) / (double)(*nu - 1);

        prm_EvalCfit(tree0, uu, xyz, NULL, NULL);

        for (ivar = 0; ivar < nvar; ivar++) {
            (*cfit)[nvar*iu+ivar] = xyz[ivar];
        }
    }

#ifdef GRAFIC
    plotCfit(nvrt, nvar, u, var, *nu, *cfit,
             "~U~var~Vertices and Knots in prm_BestCfit");
#endif

cleanup:
    prm_FreeCfit(&tree1);
    prm_FreeCfit(&tree0);

    FREE(resid1);
    FREE(resid0);

    DPRINT4("%s --> rmserr=%f, maxerr=%f, status=%d}",
            routine, *rmserr, *maxerr, status);
    return status;
}



/*
 ********************************************************************************
 *                                                                              *
 * prm_EvalCfit -- evaluate the Cfit at the given uu                            *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_EvalCfit(cfitTree tree,                  /* (in)   Cfit Tree */
             double   uu,                    /* (in)   independent variable */
             double   var[],                 /* (out)  dependent variables */
  /*@null@*/ double   *dt,                   /* (out)  optional first derivatives */
  /*@null@*/ double   *dtt)                  /* (out)  optional second derivatives */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    int         nvar  =     tree.nvar;
    int         nvar2 = 2 * tree.nvar;

    int         ivar, icel;
    double      ss, umin, umax;

    double      *xyz = NULL;
    double      *Dt  = NULL;
    double      *Dtt = NULL;

    ROUTINE(prm_EvalCfit);
//    DPRINT2("%s(uu=%f) {",
//            routine, uu);

    /* ----------------------------------------------------------------------- */

    MALLOC(xyz, double, nvar);
    MALLOC(Dt,  double, nvar);
    MALLOC(Dtt, double, nvar);

    /*
     * initialize the dependent variables
     */
    for (ivar = 0; ivar < nvar; ivar++) {
        var[ivar] = 0;
    }
    if (dt != NULL) {
        for (ivar = 0; ivar < nvar; ivar++) {
            dt[ivar] = 0;
        }
    }
    if (dtt != NULL) {
        for (ivar = 0; ivar < nvar; ivar++) {
            dtt[ivar] = 0;
        }
    }

    /*
     * start evaluation at the root ot the Cfit Tree
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

        ss = (uu - umin) / (umax - umin);

        if (dt == NULL) {
            evalCubic(ss, nvar,
                      &(tree.knot[nvar2*(tree.cell[icel].wKnot)]),
                      &(tree.knot[nvar2*(tree.cell[icel].eKnot)]),
                      xyz, NULL, NULL);

            for (ivar = 0; ivar < nvar; ivar++) {
                var[ivar] += xyz[ivar];
            }
        } else if (dtt == NULL) {
            evalCubic(ss, nvar,
                      &(tree.knot[nvar2*(tree.cell[icel].wKnot)]),
                      &(tree.knot[nvar2*(tree.cell[icel].eKnot)]),
                      xyz, Dt, NULL);

            for (ivar = 0; ivar < nvar; ivar++) {
                var[ivar] += xyz[ivar];
                dt[ ivar] += Dt[ ivar] / (umax - umin);
            }
        } else {
            evalCubic(ss, nvar,
                      &(tree.knot[nvar2*(tree.cell[icel].wKnot)]),
                      &(tree.knot[nvar2*(tree.cell[icel].eKnot)]),
                      xyz, Dt, Dtt);

            for (ivar = 0; ivar < nvar; ivar++) {
                var[ivar] += xyz[ivar];
                dt[ ivar] += Dt[ ivar] / (umax - umin);
                dtt[ivar] += Dtt[ivar] / (umax - umin) / (umax - umin);
            }
        }

        /*
         * if at bottom of hierarchy, stop
         */
        if (tree.cell[icel].child == 0) {
            break;

        /*
         * otherwise go to the appropriate child
         */
        } else if (uu <= (umin+umax)/2) {
            icel  = tree.cell[icel].child;
        } else {
            icel  = tree.cell[icel].child + 1;
        }
    }

 cleanup:
    FREE(Dtt);
    FREE(Dt );
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
 * prm_FreeCfit -- free Cells and Knots associated with Cfit Tree               *
 *                                                                              *
 ********************************************************************************
 */
extern int
prm_FreeCfit(cfitTree *tree)                 /* (in)   pointer to Cfit Tree */
{
    int         status = EGADS_SUCCESS;      /* (out)  return status */
                                             /*        EGADS_SUCCESS */

    ROUTINE(prm_FreeCfit);
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
    tree->ncel = 0;
    tree->nknt = 0;

// cleanup:
    DPRINT2("%s --> status=%d}", routine, status);
    return status;
}
