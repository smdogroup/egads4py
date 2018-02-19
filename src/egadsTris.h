/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Triangulation Header
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#define CHUNK    256                   /* allocation chunk */

#define MAXELEN 2048                   /* max Edge length */
#define DEGENUV 1.e-14


  typedef struct {
    int keys[3];
  } KEY;

  typedef struct {
    int    close;
    double xyz[3];
  } DATA;

  typedef struct {
    KEY  key;
    DATA data;
  } ENTRY;

  typedef struct element {
    ENTRY          item;
    /*@null@*/
    struct element *next;
  } ELEMENT;
  

  typedef struct {
    int    type;			/* Topology type */
    int    edge;			/* Edge tessellation index */
    int    index;			/* index for Node or Edge */
    double xyz[3];
    double uv[2];
  } triVert;


/*                                 neighbors
 *               0            tri-side   vertices
 *              / \               0        1 2
 *             /   \              1        0 2
 *            /     \             2        0 1
 *           1-------2 
 */


  typedef struct {
    int    indices[3];		/* triVert indices for triangle */
    int    neighbors[3];	/* neighboring tri index (- seg) */
    double mid[3];              /* midpoint xyz */
    double area;                /* area of triangle */
    short  mark;                /* temp storage for marking tri */
    short  close;               /* mid marked too close to edge */
    short  hit;                 /* hit this before */
    short  count;
  } triTri;


  typedef struct {
    int indices[2];             /* indices for the bounding segment */
    int neighbor;		/* triangle neighbor index */
    int edge;                   /* owning Edge index (+/- sense) */
    int index;                  /* edge tessellation index */
  } triSeg;


  typedef struct {
    egObject *face;		/* face object */
    int      fIndex;            /* face index */
    int      orUV;              /* face sense */
    int      planar;            /* the face is a plane (== 1) */
    int      phase;             /* tessellation phase */
    double   VoverU;            /* UV ratio in physical coordinates */
    double   maxlen;		/* maximum length for side */
    double   chord;		/* sag for triangulation */
    double   dotnrm;            /* angle for the dihedral */
    double   accum;
    double   edist2;            /* largest edge segment */
    double   eps2;              /* smallest edge segment */
    double   devia2;            /* largest edge deviation */
    double   minlen;            /* smallest length for side */
    double   range[4];          /* face parameter range */
    double   qparm[3];          /* quadding parameters */
    double   *uvs;
    int      orCnt;             /* orientataion error count */
    int      maxPts;            /* max points for Face */
    int      mverts;		/* triangulation vert storage */
    int      nverts;
    triVert  *verts;
    int      mtris;		/* triangle storage */
    int      ntris;
    triTri   *tris;
    int      msegs;		/* bounding segment (edge) storage */
    int      nsegs;
    triSeg   *segs;
    int      nfrvrts;           /* number of verts in the frame */
    int      mframe;            /* starting triangle frame storage */
    int      nframe;
    int      *frame;
    int      mloop;             /* frame loop storage */
    int      nloop;
    int      *loop;
    int      lens[5];           /* quading sizes */
    int      numElem;		/* hash table -- number of elements */
    ELEMENT  **hashTab;
  } triStruct;


  typedef struct {
    int   sleft;		/* left  segment in front */
    int   i0;                   /* left  vertex index */
    int   i1;                   /* right vertex index */
    int   sright;		/* right segment in front */
    short snew;                 /* is this a new segment? */
    short mark;                 /* is this segment marked? */
  } Front;
  
  typedef struct {
    int    mfront;
    int    nfront;
    int    mpts;
    int    npts;
    int   *pts;
    int    nsegs;
    int   *segs;
    Front *front;
  } fillArea;

  typedef struct {
    int node1;                  /* 1nd node number for edge */
    int node2;                  /* 2nd node number for edge */
    int *tri;                   /* 1st triangle storage or NULL for match */
    int thread;                 /* thread to next face with 1st node number */
  } connect;


  /* structure to pass data to the thread for a block */
  typedef struct {
    void     *mutex;            /* the mutex or NULL for single thread */
    long     master;            /* master thread ID */
    int      end;               /* end of loop */
    int      index;             /* current loop index */
    /*@dependent@*/
    int      *mark;             /* do the index or NULL (for all) */
    egObject *tess;             /* Tessellation Object */
    egTessel *btess;            /* tessellation structure -- Edges */
    egObject *body;             /* Body Object to Tessellate */
    /*@dependent@*/
    egObject **faces;           /* Face Object list */
    /*@dependent@*/
    egObject **edges;           /* Edge Object list */
    double   *params;           /* Tessellation parameters */
    double   *tparam;
    double   qparam[3];         /* quading parameters */
    void     *ptr;              /* user pointer */
  } EMPtess;
