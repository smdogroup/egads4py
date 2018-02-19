/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Parameterization (prm)  Header
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#ifndef _PRM_H_
#define _PRM_H_

#define PRM_OK_ONEFACE        1
#define PRM_OK_SIMPLE         2
#define PRM_OK_JONES          3
#define PRM_OK_AXIAL          4
#define PRM_OK_PLANAR         5
#define PRM_OK_FLOATER        6
#define PRM_OK_UNROLLING      7

#define PRM_NOTCONVERGED   -401
#define PRM_BADNUMVERTICES -402
#define PRM_ZEROPIVOT      -403
#define PRM_NEGATIVEAREAS  -404
#define PRM_TOLERANCEUNMET -405
#define PRM_NOGLOBALUV     -406
#define PRM_BADPARAM       -407
#define PRM_BADDIVISION    -408
#define PRM_CANNOTFORMLOOP -409
#define PRM_WIGGLEDETECTED -410
#define PRM_INTERNAL       -499

/*
 * structures for holding Coordinate-type info
 */
typedef struct {
    double      x;			/* x-coordinate */
    double      y;                      /* y-coordinate */
    double      z;                      /* z-coordinate */
} prmXYZ;

typedef struct {
    double      u0;                     /* u-coordinate at point 0 of Tri */
    double      v0;                     /* v-coordinate at point 0 of Tri */
    double      u1;                     /* u-coordinate at point 1 of Tri */
    double      v1;                     /* v-coordinate at point 1 of Tri */
    double      u2;                     /* u-coordinate at point 2 of Tri */
    double      v2;                     /* v-coordinate at point 2 of Tri */
} prmUVF;

typedef struct {
    double      u;                      /* first  Parameter */
    double      v;                      /* second Parameter */
} prmUV;

typedef struct {
  int indices[3];			/* vert indices (1 bias) */
  int neigh[3];				/* neighboring tris (1 bias) */
  int own;				/* owning face index in source */
} prmTri;


/*
 * structure that holds Cells for a Cfit
 */
typedef struct {
    int              wKnot;             /* index of west Knot */
    int              eKnot;             /* index of east Knot */
    int              wCell;             /* west neighbor (or 0) */
    int              eCell;             /* east neighbor (or 0) */
    int              child;             /* index of west child (or 0) */
    int              count;             /* number of Vertices in Cell */
    double           umin;              /* minimum U */
    double           umax;              /* maximum U */
    double           rmserr;            /* RMS     error at Vertices */
    double           maxerr;            /* maximum error at Vertices */
} cfitCell;

/*
 * structure that holds a hierarchical Cfit
 */
typedef struct {
    int              nvar;              /* number of dependent vars */
    int              nu;                /* total Knots */
    int              ncel;              /* number of Cells */
    cfitCell         *cell;             /* array  of Cells */
    int              nknt;              /* number of Knots */
    double           *knot;             /* array  of Knots */
} cfitTree;

/*
 * structure that holds Cells for a Grid
 */
typedef struct {
    int              swKnot;            /* index of southwest Knot */
    int              seKnot;            /* index of southeast Knot */
    int              nwKnot;            /* index of northwest Knot */
    int              neKnot;            /* index of northeast Knot */
    int              wCell;             /* west  neighbor (or 0) */
    int              eCell;             /* east  neighbor (or 0) */
    int              sCell;             /* south neighbor (or 0) */
    int              nCell;             /* north neighbor (or 0) */
    int              dtype;             /* division type:
                                           0 = not divided
                                           1 = divided vertically (w and e)
                                           2 = divided horizontally (s and n)
                                           3 = quartered (sw, se, nw, and ne) */
    int              child;             /* index of southwest child (or 0) */
    int              count;             /* number of Vertices in Cell */
    double           umin;              /* minimum U */
    double           umax;              /* maximum U */
    double           vmin;              /* minimum V */
    double           vmax;              /* maximum V */
    double           rmserr;            /* RMS     error at Vertices */
    double           maxerr;            /* maximum error at Vertices */
} gridCell;

/*
 * structure that holds a hierarchical Grid
 */
typedef struct {
    int              nvar;              /* number of dependent vars */
    int              nu;                /* total Knots in U-direction */
    int              nv;                /* total Knots in V-direction */
    int              ncel;              /* number of Cells */
    gridCell         *cell;             /* array  of Cells */
    int              nknt;              /* number of Knots */
    double           *knot;             /* array  of Knots */
} gridTree;

/**********************************************************************/
/**********************************************************************/

/*
 * create a parameterization from a set of Vertices
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_CreateU(int      nvrt,              /* (in)   number of Vertices */
            double   u[],               /* (out)  array  of Parameters */
            prmXYZ   xyz[],             /* (in)   array  of Coordinates */
            double   tol,               /* (in)   tol for periodic check */
            int      *periodic);        /* (out)  = 0  no periodicity */
                                        /*        = 1  periodic in U */

/*
 * create a parameterization from a set of Triangles and Vertices
 *      returns:  0  NO parameterization found!
 *                1  Surface has one Face
 *                2  simple projection
 *                3  jones transformation
 *                4  axial projection
 *                5  planar projection
 *                6  floater's parameterization
 *                7  unrolling
 *               <0  error condition
 */
extern int
prm_CreateUV(int      UVtype,           /* (in)   type of transformation (see returns) */
             int      ntri,             /* (in)   number of Triangles */
             prmTri   tri[],            /* (in)   array  of Triangles */
  /*@null@*/ prmUVF   uvf[],            /* (in)   array  of UVF for Triangles */
             int      nvrt,             /* (in)   number of Vertices */
  /*@null@*/ int      ptype[],          /* (in)  point types */
  /*@null@*/ int      pindx[],          /* (in)  point indices */
             prmUV    uv[],             /* (out)  array  of Parameters */
             prmXYZ   xyz[],            /* (in)   array  of Coordinates */
             int      *periodic,        /* (out)  = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        = 2  periodic in V */
                                        /*        =-1  periodicity problem */
             int      *ppnts[]);        /* (out)  pointer to indices of periodic points
                                                  NOTE: user must EG_free after use */

/*
 * reparameterize based upon arc-lengths
 *      returns:  CAPS_SUCCESS
 *                PRM_NOTCONVERGED
 */
extern int
prm_SmoothU(int      periodic,          /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
            int      nvrt,              /* (in)   number of Vertices */
            int      nvar,              /* (in)   number of dependent vars */
            double   u[],               /* (both) array  of parameters */
            double   xyz[]);            /* (in)   array  of Coordinates */

/*
 * reparameterize based upon arc-lengths and angles
 *      returns:  CAPS_SUCCESS
 *                PRM_NOTCONVERGED
 *                PRM_NEGATIVEAREAS
 */
extern int
prm_SmoothUV(int      type,             /* (in)   smoothing type =1 for boundary only
                                                                 =2 for interior only
                                                                 =3 for both */
             int      periodic,         /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        = 2  periodic in V */
  /*@null@*/ int      ppnts[],          /* (in)   indices of periodic points */
             int      ntri,             /* (in)   number of Triangles */
             prmTri   tri[],            /* (in)   array  of Triangles */
             int      nvrt,             /* (in)   number of Vertices */
             int      nvar,             /* (in)   number of dependent vars */
             prmUV    uv[],             /* (both) array  of Parameters */
             double   xyz[]);           /* (in)   array  of Coordinates */

/*
 * normalize a parameterization
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_NormalizeU(double  halo,            /* (in)   halo size */
               int     periodic,        /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
               int     nvrt,            /* (in)   number of Vertices */
               double  u[]);            /* (both) array  of Parameters */

/*
 * normalize a parameterization
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_NormalizeUV(double  halo,           /* (in)   halo size */
                int     periodic,       /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        = 2  periodic in V */
                int     nvrt,           /* (in)   number of Vertices */
                prmUV   uv[]);          /* (both) array  of Parameters */

/*
 * set a global Grid size limit
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_LimitGridSize(int      limit);      /* (in)   global size limit */

/*
 * set up a fixed Cfit for a set of Vertices
 *      returns:  CAPS_SUCCESS
 *                PRM_NOTCONVERGED
 *                PRM_BADPARAM
 */
extern int
prm_FixedCfit(int      nvrt,            /* (in)   number of Vertices */
              int      nvar,            /* (in)   number of dependent vars */
              double   u[],             /* (in)   array  of Parameters */
              double   var[],           /* (in)   array  of dependent vars */
              int      periodic,        /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        =11  left fixed, rite fixed */
                                        /*        =12  left fixed, rite zero  */
                                        /*        =21  left zero,  rite fixed */
                                        /*        =22  left zero,  rite zero  */
              int      nu,              /* (in)   number of Knots */
              double   cfit[],          /* (out)  Cfit of Knots */
              double   *rmserr,         /* (out)  RMS     error at Vertices */
              double   *maxerr);        /* (out)  maximum error at Vertices */

/*
 * set up a fixed Grid for a set of Vertices
 *      returns:  CAPS_SUCCESS
 *                PRM_NOTCONVERGED
 *                PRM_BADPARAM
 */
extern int
prm_FixedGrid(int      nvrt,            /* (in)   number of Vertices */
              int      nvar,            /* (in)   number of dependent vars */
              prmUV    uv[],            /* (in)   array  of Parameters */
              double   var[],           /* (in)   array  of dependent vars */
              int      ntri,            /* (in)   number of Triangles (optional) */
   /*@null@*/ prmTri   tri[],           /* (in)   array  if Triangles (optional) */
              int      periodic,        /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        = 2  periodic in V */
   /*@null@*/ int      ppnts[],         /* (in)   indices of periodic points */
              int      nu,              /* (in)   number of Knots in U-dirn */
              int      nv,              /* (in)   number of Knots in V-dirn */
              double   grid[],          /* (out)  Grid of Knots */
              double   *rmserr,         /* (out)  RMS     error at Vertices */
              double   *maxerr,         /* (out)  maximum error at Vertices */
              double   *dotmin);        /* (out)  minimum dot product */

/*
 * set up the best (ie, coarsest) Cfit for a set of Vertices
 *      returns:  CAPS_SUCCESS
 *                PRM_TOLERANCEUNMET
 *                PRM_NOTCONVERGED
 *                PRM_BADPARAM
 */
extern int
prm_BestCfit(int      nvrt,             /* (in)   number of Vertices */
             int      nvar,             /* (in)   number of dependent vars */
             double   u[],              /* (in)   array  of Parameters */
             double   var[],            /* (in)   array  of dependent vars */
             double   tol,              /* (in)   tolerance on maximum error */
             int      periodic,         /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        =11  left fixed, rite fixed */
                                        /*        =12  left fixed, rite zero  */
                                        /*        =21  left zero,  rite fixed */
                                        /*        =22  left zero,  rite zero  */
             int      *nu,              /* (in)   limit on nu if > 0 */
                                        /* (out)  number of Knots */
             double   *cfit[],          /* (out)  pointer to Cfit of Knots
                                                  NOTE: user must EG_free after use */
             double   *rmserr,          /* (out)  RMS     error at Vertices */
             double   *maxerr);         /* (out)  maximum error at Vertices */

/*
 * set up the best (ie, coarsest) Grid for a set of Vertices
 *      returns:  CAPS_SUCCESS
 *                PRM_TOLERANCEUNMET
 *                PRM_NOTCONVERGED
 *                PRM_BADPARAM
 */
extern int
prm_BestGrid(int      nvrt,             /* (in)   number of Vertices */
             int      nvar,             /* (in)   number of dependent vars */
             prmUV    uv[],             /* (in)   array  of Parameters */
             double   var[],            /* (in)   array  of dependent vars */
             int      ntri,             /* (in)   number of Tris (optional) */
  /*@null@*/ prmTri   tri[],            /* (in)   array  if Tris (optional) */
             double   tol,              /* (in)   tolerance on maximum error */
             int      periodic,         /* (in)   = 0  no periodicity */
                                        /*        = 1  periodic in U */
                                        /*        = 2  periodic in V */
  /*@null@*/ int      ppnts[],          /* (in)   indices of periodic points */
             int      *nu,              /* (in)   limit on nu if > 0 */
                                        /* (out)  number of Knots in U-dirn */
             int      *nv,              /* (in)   limit on nv if > 0 */
                                        /* (out)  number of Knots in V-dirn */
             double   *grid[],          /* (out)  pointer to Grid of Knots
                                                  NOTE: user must EG_free after use */
             double   *rmserr,          /* (out)  RMS     error at Vertices */
             double   *maxerr,          /* (out)  maximum error at Vertices */
             double   *dotmin);         /* (out)  minimum dot product */

/*
 * evaluate the Cfit at the given uu
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_EvalCfit(cfitTree tree,             /* (in)   Cfit Tree */
             double   uu,               /* (in)   independent variable */
             double   *var,             /* (out)  dependent variables */
  /*@null@*/ double   *dt,              /* (out)  optional first derivatives */
  /*@null@*/ double   *dtt);            /* (out)  optional second derivatives */

/*
 * evaluate the Grid at the given uv
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_EvalGrid(gridTree tree,             /* (in)   Grid Tree */
             prmUV    uv,               /* (in)   independent variable */
             double   *var,             /* (out)  dependent variables */
  /*@null@*/ double   *du,              /* (out)  optional first u-derivatives */
  /*@null@*/ double   *dv,              /* (out)  optional first v-derivatives */
  /*@null@*/ double   *duu,             /* (out)  optional second u-derivatives */
  /*@null@*/ double   *duv,             /* (out)  optional cross derivatives */
  /*@null@*/ double   *dvv);            /* (out)  optional second v-derivatives */

/*
 * free Cells and Knots associated with Cfit Tree
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_FreeCfit(cfitTree *tree);           /* (in)   pointer to Cfit Tree */

/*
 * free Cells and Knots associated with Grid Tree
 *      returns:  CAPS_SUCCESS
 */
extern int
prm_FreeGrid(gridTree *tree);           /* (in)   pointer to Grid Tree */

#endif /*_PRM_H_*/
