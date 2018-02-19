# Import string stuff
from libc.string cimport const_char

# Import numpy 
cimport numpy as np
import numpy as np

cdef extern from "egadsErrors.h":
    enum:
        _EGADS_TESSTATE"EGADS_TESSTATE"
        _EGADS_EXISTS"EGADS_EXISTS"
        _EGADS_ATTRERR"EGADS_ATTRERR"
        _EGADS_TOPOCNT"EGADS_TOPOCNT"
        _EGADS_OCSEGFLT"EGADS_OCSEGFLT"
        _EGADS_BADSCALE"EGADS_BADSCALE"
        _EGADS_NOTORTHO"EGADS_NOTORTHO"
        _EGADS_DEGEN"EGADS_DEGEN"
        _EGADS_CONSTERR"EGADS_CONSTERR"
        _EGADS_TOPOERR"EGADS_TOPOERR"
        _EGADS_GEOMERR"EGADS_GEOMERR"
        _EGADS_NOTBODY"EGADS_NOTBODY"
        _EGADS_WRITERR"EGADS_WRITERR"
        _EGADS_NOTMODEL"EGADS_NOTMODEL"
        _EGADS_NOLOAD"EGADS_NOLOAD"
        _EGADS_RANGERR"EGADS_RANGERR"
        _EGADS_NOTGEOM"EGADS_NOTGEOM"
        _EGADS_NOTTESS"EGADS_NOTTESS"
        _EGADS_EMPTY"EGADS_EMPTY"
        _EGADS_NOTTOPO"EGADS_NOTTOPO"
        _EGADS_REFERCE"EGADS_REFERCE"
        _EGADS_NOTXFORM"EGADS_NOTXFORM"
        _EGADS_NOTCNTX"EGADS_NOTCNTX"
        _EGADS_MIXCNTX"EGADS_MIXCNTX"
        _EGADS_NODATA"EGADS_NODATA"
        _EGADS_NONAME"EGADS_NONAME"
        _EGADS_INDEXERR"EGADS_INDEXERR"
        _EGADS_MALLOC"EGADS_MALLOC"
        _EGADS_NOTOBJ"EGADS_NOTOBJ"
        _EGADS_NULLOBJ"EGADS_NULLOBJ"
        _EGADS_NOTFOUND"EGADS_NOTFOUND"
        _EGADS_SUCCESS"EGADS_SUCCESS"
        _EGADS_OUTSIDE"EGADS_OUTSIDE"

cdef extern from "egadsTypes.h":
    cdef struct egObject:
        short   oclass
        short   mtype

    ctypedef egObject* ego

    enum:
        EGADS_NIL"NIL"
        EGADS_EMPTY"EMPTY"
        EGADS_REFERENCE"REFERENCE"
        EGADS_PCURVE"PCURVE"
        EGADS_CURVE"CURVE"
        EGADS_SURFACE"SURFACE"
        EGADS_NODE"NODE"
        EGADS_EDGE"EDGE"
        EGADS_LOOP"LOOP"
        EGADS_FACE"FACE"
        EGADS_SHELL"SHELL"
        EGADS_BODY"BODY"
        EGADS_MODEL"MODEL"

        # /* PCURVES & CURVES */
        EGADS_LINE"LINE"
        EGADS_CIRCLE"CIRCLE"
        EGADS_ELLIPSE"ELLIPSE"
        EGADS_PARABOLA"PARABOLA"
        EGADS_HYPERBOLA"HYPERBOLA"
        EGADS_TRIMMED"TRIMMED"
        EGADS_BEZIER"BEZIER"
        EGADS_BSPLINE"BSPLINE"
        EGADS_OFFSET"OFFSET"

        # /* SURFACES */
        EGADS_PLANE"PLANE"
        EGADS_SPHERICAL"SPHERICAL"
        EGADS_CYLINDRICAL"CYLINDRICAL"
        EGADS_REVOLUTION"REVOLUTION"
        EGADS_TOROIDAL"TOROIDAL"
        EGADS_CONICAL"CONICAL"
        EGADS_EXTRUSION"EXTRUSION"

        # /* TOPOLOGY */
        EGADS_SREVERSE"SREVERSE"
        EGADS_NOMTYPE"NOMTYPE"
        EGADS_SFORWARD"SFORWARD"
        EGADS_ONENODE"ONENODE"
        EGADS_TWONODE"TWONODE"
        EGADS_OPEN"OPEN"
        EGADS_CLOSED"CLOSED"
        EGADS_DEGENERATE"DEGENERATE"
        EGADS_WIREBODY"WIREBODY"
        EGADS_FACEBODY"FACEBODY"
        EGADS_SHEETBODY"SHEETBODY"
        EGADS_SOLIDBODY"SOLIDBODY"

        # # /* ATTRIBUTE TYPES */
        # EGADS_ATTRINT   
        # EGADS_ATTRREAL  
        # EGADS_ATTRSTRING
        # EGADS_ATTRCSYS  

        # # /* SOLID BOOLEAN OPERATIONS */
        # EGADS_SUBTRACTION
        # EGADS_INTERSECTION
        # EGADS_FUSION      

        # # /* SOLID BODY TYPES */
        # EGADS_BOX       
        # EGADS_SPHERE    
        # EGADS_CONE      
        # EGADS_CYLINDER  
        # EGADS_TORUS     

        # # /* ISOCLINE TYPES */
        # EGADS_UISO      
        # EGADS_VISO      

        # # /* FACE SOURCE TYPES */
        # EGADS_NODEOFF   
        # EGADS_EDGEOFF   
        # EGADS_FACEDUP   
        # EGADS_FACECUT   
        # EGADS_FACEOFF

cdef extern from "egads.h":
    void EG_revision( int *major, int *minor, const char **OCCrev )
    int  EG_open( ego *context )
    int  EG_loadModel( ego context, int bflg, const char *name, 
                       ego *model )
    int  EG_saveModel( const ego model, const char *name )
    int  EG_deleteObject( ego object )
    int  EG_makeTransform( ego context, const double *xform, 
                           ego *oform )
    int  EG_getTransformation( const ego oform, double *xform )
    int  EG_getContext( ego object, ego *context )
    int  EG_setOutLevel( ego context, int outLevel )
    int  EG_getInfo( const ego object, int *oclass, int *mtype, 
                     ego *topObj, ego *prev, ego *next )
    int  EG_copyObject( const ego object,void *oform,
                        ego *copy )
    int  EG_flipObject( const ego object, ego *flippedCopy )
    int  EG_close( ego context )

    int  EG_attributeAdd( ego obj, const char *name, int type, int len,
                          const int    *ints, 
                          const double *reals,
                          const char   *str )
    int  EG_attributeDel( ego object, const char *name )
    int  EG_attributeNum( const ego obj, int *num )
    int  EG_attributeGet( const ego obj, int index, const char **name,
                          int *atype, int *len, 
                          const int    **ints,
                          const double **reals, 
                          const char   **str )
    int  EG_attributeRet( const ego obj, const char *name, int *atype, 
                          int *len, const int    **ints,
                          const double **reals, 
                          const char   **str )
    int  EG_attributeDup( const ego src, ego dst )

    int  EG_getGeometry( const ego geom, int *oclass, int *mtype,
                         ego *refGeom, int **ivec, double **rvec )
    int  EG_makeGeometry( ego context, int oclass, int mtype,
                          ego refGeom,
                          const int *ivec, 
                          const double *rvec, ego *geom )
    int  EG_getRange( const ego geom, double *range, int *periodic )
    int  EG_evaluate( const ego geom, const double *param, 
                      double *results )
    int  EG_invEvaluate( const ego geom, double *xyz, double *param,
                         double *results )
    int  EG_invEvaluateGuess( const ego geom, double *xyz, 
                              double *param, double *results )
    int  EG_arcLength( const ego geom, double t1, double t2,
                       double *alen )
    int  EG_curvature( const ego geom, const double *param, 
                       double *results )
    int  EG_approximate( ego context, int maxdeg, double tol,
                         const int *sizes, const double *xyzs,
                         ego *bspline )
    int  EG_fitTriangles( ego context, int npts, double *xyzs,
                          int ntris, const int *tris,
                          const int *tric, double tol,
                          ego *bspline )
    int  EG_otherCurve( const ego surface, const ego curve,
                        double tol, ego *newcurve )
    int  EG_isSame( const ego geom1, const ego geom2 )
    int  EG_isoCline( const ego surface, int UV, double val,
                      ego *newcurve )
    int  EG_convertToBSpline( ego geom, ego *bspline )
    int  EG_convertToBSplineRange( ego geom, const double *range,
                                   ego *bspline )

    int  EG_getTolerance( const ego topo, double *tol )
    int  EG_setTolerance(       ego topo, double  tol )
    int  EG_getTopology( const ego topo, ego *geom, int *oclass, 
                         int *type, double *limits, 
                         int *nChildren, ego **children, int **sense )
    int  EG_makeTopology( ego context, ego geom, int oclass,
                          int mtype, double *limits,
                          int nChildren, ego *children,
                          int *senses, ego *topo )
    int  EG_makeLoop( int nedge, ego *edges, ego geom,
                      double toler, ego *result )
    int  EG_getArea( ego object, const double *limits,
                     double *area )
    int  EG_makeFace( ego object, int mtype,
                      const double *limits, ego *face )
    int  EG_getBodyTopos( const ego body, ego src,
                          int oclass, int *ntopo, ego **topos )
    int  EG_indexBodyTopo( const ego body, const ego src )
    int  EG_inTopology( const ego topo, const double *xyz )
    int  EG_inFace( const ego face, const double *uv )
    int  EG_getEdgeUV( const ego face, const ego edge, int sense,
                       double t, double *UV )
    int  EG_getBody( const ego object, ego *body )
    int  EG_makeSolidBody( ego context, int stype, const double *rvec,
                           ego *body )
    int  EG_getBoundingBox( const ego topo, double *bbox )
    int  EG_getMassProperties( const ego topo, double *result )
    int  EG_isEquivalent( const ego topo1, const ego topo2 )
    int  EG_sewFaces( int nobj, const ego *objs, double toler,
                     int flag, ego *result )
    int  EG_replaceFaces( const ego body, int nobj, ego *objs,
                          ego *result )
    int  EG_mapBody( const ego sBody,   const ego dBody,
                     const char *fAttr, ego *mapBody )
    int  EG_matchBodyFaces( const ego body1, const ego body2,
                            double toler, int *nmatch, int **match )

    int  EG_setTessParam( ego context, int iparam, double value,
                          double *oldvalue )
    int  EG_makeTessGeom( ego obj, double *params, int *sizes, 
                          ego *tess )
    int  EG_getTessGeom( const ego tess, int *sizes, double **xyz )

    int  EG_makeTessBody( ego object, double *params, ego *tess )
    int  EG_remakeTess( ego tess, int nobj, ego *objs, 
                        double *params )
    int  EG_mapTessBody( ego tess, ego body, ego *mapTess )
    int  EG_locateTessBody( const ego tess, int npt, const int *ifaces,
                            const double *uv, int *itri, 
                            double *results )

    int  EG_getTessEdge( const ego tess, int eIndex, int *len, 
                         const double **xyz, const double **t )
    int  EG_getTessFace( const ego tess, int fIndex, int *len, 
                         const double **xyz, const double **uv, 
                         const int **ptype, const int **pindex, 
                         int *ntri, const int **tris, 
                         const int **tric )
    int  EG_getTessLoops( const ego tess, int fIndex, int *nloop,
                          const int **lIndices )
    int  EG_getTessQuads( const ego tess, int *nquad,
                          int **fIndices )
    int  EG_makeQuads( ego tess, double *params, int fIndex )
    int  EG_getQuads( const ego tess, int fIndex, int *len, 
                      const double **xyz, const double **uv, 
                      const int **ptype, const int **pindex, 
                      int *npatch )
    int  EG_getPatch( const ego tess, int fIndex, int patch, 
                      int *nu, int *nv, const int **ipts, 
                      const int **bounds )
                               
    int  EG_insertEdgeVerts( ego tess, int eIndex, int vIndex, 
                             int npts, double *t )
    int  EG_deleteEdgeVert( ego tess, int eIndex, int vIndex, 
                            int dir )
    int  EG_moveEdgeVert( ego tess, int eIndex, int vIndex, 
                          double t )

    int  EG_solidBoolean( const ego src, const ego tool, int oper, 
                          ego *model )
    int  EG_intersection( const ego src, const ego tool, int *nedge, 
                          ego **facEdg, ego *model )
    int  EG_imprintBody( const ego src, int nedge, const ego *facEdg, 
                         ego *result )
    int  EG_filletBody( const ego src, int nedge, const ego *edges, 
                        double radius,
                        ego *result, int **facemap )
    int  EG_chamferBody( const ego src, int nedge, const ego *edges, 
                         const ego *faces, double dis1, double dis2, 
                         ego *result, int **facemap )
    int  EG_hollowBody( const ego src, int nface, const ego *faces, 
                        double offset, int join,
                        ego *result, int **facemap )
    int  EG_extrude( const ego src, double dist, const double *dir, 
                     ego *result )
    int  EG_rotate( const ego src, double angle, const double *axis, 
                    ego *result )
    int  EG_sweep( const ego src, const ego spine, int mode,
                   ego *result )
    int  EG_loft( int nsec, const ego *secs, int opt, ego *result )
    int  EG_blend( int nsec, const ego *secs, double *rc1,
                   double *rcN, ego *result )
    int  EG_ruled( int nsec, const ego *secs, ego *result )
