# Import numpy 
cimport numpy as np
import numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import the definitions
from egads4py cimport *

def revision():
    cdef int major, minor
    cdef const char *occrev
    EG_revision(&major, &minor, &occrev)
    return major, minor, occrev

cdef class Ego:
    cdef ego context
    cdef ego model

    def __init__(self, contxt=None):
        if contxt:
            self.context = <ego>contxt
        else:
            EG_open(&self.context)
        self.model = NULL

    def __dealloc__(self):
        if self.model:
            EG_deleteObject(self.model)

    def setOutLevel(self, int outlevel):
        EG_setOutLevel(self.context, outlevel)

    def loadModel(self, str filename, int bflag=1):
        if self.model:
            EG_deleteObject(self.model)
        stat = EG_loadModel(self.context, bflag, filename, &self.model)

    def saveModel(self, str filename):
        stat = EG_saveModel(self.model, filename)

    def makeTransform(self, np.ndarray[double, ndim=1, mode='c'] xform):
        new = Ego()
        EG_makeTransform(self.context, <double*>xform.data, &new.model)

    def getTransform(self):
        cdef np.ndarray xform = np.zeros(12, np.double)
        EG_getTransformation(self.model, <double*>xform.data)
        return xform

    def copy(self):
        new = Ego()
        EG_copyObject(self.model, NULL, &new.model)

    def flip(self):
        flip = Ego()
        EG_flipObject(self.model, &flip.model)

    # Set/retrieve attributes
    # int  EG_attributeAdd( ego obj, const char *name, int type, int len,
    #                       const int    *ints, 
    #                       const double *reals,
    #                       const char   *str )
    # int  EG_attributeDel( ego object, const char *name )
    # int  EG_attributeNum( const ego obj, int *num )
    # int  EG_attributeGet( const ego obj, int index, const char **name,
    #                                int *atype, int *len, 
    #                                const int    **ints,
    #                                const double **reals, 
    #                                const char   **str )
    # int  EG_attributeRet( const ego obj, const char *name, int *atype, 
    #                                int *len, const int    **ints,
    #                                          const double **reals, 
    #                                          const char   **str )
    # int  EG_attributeDup( const ego src, ego dst )

    def getRange(self):
        cdef double r[4]
        cdef int periodic
        cdef int oclass
        if (self.model.oclass == EGADS_PCURVE or 
            self.model.oclass == EGADS_CURVE or 
            self.model.oclass == EGADS_EDGE):
            EG_getRange(self.model, r, &periodic)
            return [r[0], r[1]], periodic
        elif (self.model.oclass == EGADS_SURFACE or 
              self.model.oclass == EGADS_EDGE):
            EG_getRange(self.model, r, &periodic)
            return [r[0], r[1], r[2], r[3]], periodic
        return None

    # def evaluate(self, np.ndarray[double, ndim=1, mode='c'] pt):
    #     cdef double result[3]
    #     stat = EG_evaluate(self.model, <double*>pt.data, result)

    # int  EG_invEvaluate( const ego geom, double *xyz, double *param,
    #                       double *results )
    # int  EG_invEvaluateGuess( const ego geom, double *xyz, 
    #                           double *param, double *results )
    # int  EG_arcLength( const ego geom, double t1, double t2,
    #                    double *alen )
    # int  EG_curvature( const ego geom, const double *param, 
    #                    double *results )
    # int  EG_approximate( ego context, int maxdeg, double tol,
    #                      const int *sizes, const double *xyzs,
    #                      ego *bspline )

    # int  EG_getGeometry( const ego geom, int *oclass, int *mtype,
    #                               ego *refGeom, int **ivec, double **rvec )
    # int  EG_makeGeometry( ego context, int oclass, int mtype,
    #                               ego refGeom,
    #                               const int *ivec, 
    #                               const double *rvec, ego *geom )
    # int  EG_getRange( const ego geom, double *range, int *periodic )
    # int  EG_evaluate( const ego geom, const double *param, 
    #                            double *results )
    # int  EG_invEvaluate( const ego geom, double *xyz, double *param,
    #                               double *results )
    # int  EG_invEvaluateGuess( const ego geom, double *xyz, 
    #                                    double *param, double *results )
    # int  EG_arcLength( const ego geom, double t1, double t2,
    #                             double *alen )
    # int  EG_curvature( const ego geom, const double *param, 
    #                             double *results )
    # int  EG_approximate( ego context, int maxdeg, double tol,
    #                               const int *sizes, const double *xyzs,
    #                               ego *bspline )
    # int  EG_fitTriangles( ego context, int npts, double *xyzs,
    #                                int ntris, const int *tris,
    #                                const int *tric, double tol,
    #                                ego *bspline )
    # int  EG_otherCurve( const ego surface, const ego curve,
    #                              double tol, ego *newcurve )
    # int  EG_isSame( const ego geom1, const ego geom2 )
    # int  EG_isoCline( const ego surface, int UV, double val,
    #                            ego *newcurve )
    # int  EG_convertToBSpline( ego geom, ego *bspline )
    # int  EG_convertToBSplineRange( ego geom, const double *range,
    #                                         ego *bspline )

    # int  EG_getTolerance( const ego topo, double *tol )
    # int  EG_setTolerance(       ego topo, double  tol )
    # int  EG_getTopology( const ego topo, ego *geom, int *oclass, 
    #                               int *type, double *limits, 
    #                               int *nChildren, ego **children, int **sense )
    # int  EG_makeTopology( ego context, ego geom, int oclass,
    #                                int mtype, double *limits,
    #                                int nChildren, ego *children,
    #                                int *senses, ego *topo )
    # int  EG_makeLoop( int nedge, ego *edges, ego geom,
    #                            double toler, ego *result )
    # int  EG_getArea( ego object, const double *limits,
    #                  double *area )
    # int  EG_makeFace( ego object, int mtype,
    #                   const double *limits, ego *face )
    # int  EG_getBodyTopos( const ego body, ego src,
    #                                int oclass, int *ntopo, ego **topos )
    # int  EG_indexBodyTopo( const ego body, const ego src )
    # int  EG_inTopology( const ego topo, const double *xyz )
    # int  EG_inFace( const ego face, const double *uv )
    # int  EG_getEdgeUV( const ego face, const ego edge, int sense,
    #                    double t, double *UV )
    # int  EG_getBody( const ego object, ego *body )
    # int  EG_makeSolidBody( ego context, int stype, const double *rvec,
    #                        ego *body )
    # int  EG_getBoundingBox( const ego topo, double *bbox )
    # int  EG_getMassProperties( const ego topo, double *result )
    # int  EG_isEquivalent( const ego topo1, const ego topo2 )
    # int  EG_sewFaces( int nobj, const ego *objs, double toler,
    #                  int flag, ego *result )
    # int  EG_replaceFaces( const ego body, int nobj, ego *objs,
    #                       ego *result )
    # int  EG_mapBody( const ego sBody,   const ego dBody,
    #                           const char *fAttr, ego *mapBody )
    # int  EG_matchBodyFaces( const ego body1, const ego body2,
    #                         double toler, int *nmatch, int **match )

    # int  EG_setTessParam( ego context, int iparam, double value,
    #                                double *oldvalue )
    # int  EG_makeTessGeom( ego obj, double *params, int *sizes, 
    #                                ego *tess )
    # int  EG_getTessGeom( const ego tess, int *sizes, double **xyz )

    # int  EG_makeTessBody( ego object, double *params, ego *tess )
    # int  EG_remakeTess( ego tess, int nobj, ego *objs, 
    #                              double *params )
    # int  EG_mapTessBody( ego tess, ego body, ego *mapTess )
    # int  EG_locateTessBody( const ego tess, int npt, const int *ifaces,
    #                                  const double *uv, int *itri, 
    #                                  double *results )

    # int  EG_getTessEdge( const ego tess, int eIndex, int *len, 
    #                               const double **xyz, const double **t )
    # int  EG_getTessFace( const ego tess, int fIndex, int *len, 
    #                               const double **xyz, const double **uv, 
    #                               const int **ptype, const int **pindex, 
    #                               int *ntri, const int **tris, 
    #                               const int **tric )
    # int  EG_getTessLoops( const ego tess, int fIndex, int *nloop,
    #                                const int **lIndices )
    # int  EG_getTessQuads( const ego tess, int *nquad,
    #                                int **fIndices )
    # int  EG_makeQuads( ego tess, double *params, int fIndex )
    # int  EG_getQuads( const ego tess, int fIndex, int *len, 
    #                               const double **xyz, const double **uv, 
    #                               const int **ptype, const int **pindex, 
    #                               int *npatch )
    # int  EG_getPatch( const ego tess, int fIndex, int patch, 
    #                            int *nu, int *nv, const int **ipts, 
    #                            const int **bounds )
                               
    # int  EG_insertEdgeVerts( ego tess, int eIndex, int vIndex, 
    #                                   int npts, double *t )
    # int  EG_deleteEdgeVert( ego tess, int eIndex, int vIndex, 
    #                                  int dir )
    # int  EG_moveEdgeVert( ego tess, int eIndex, int vIndex, 
    #                                double t )
 
    # int  EG_openTessBody( ego tess )
    # int  EG_initTessBody( ego object, ego *tess )
    # int  EG_statusTessBody( ego tess, ego *body, int *state, int *np )
    # int  EG_setTessEdge( ego tess, int eIndex, int len,
    #                               const double *xyz, const double *t )
    # int  EG_setTessFace( ego tess, int fIndex, int len,
    #                               const double *xyz, const double *uv,
    #                               int ntri, const int *tris )
    # int  EG_localToGlobal( const ego tess, int index, int local,
    #                                 int *global )
    # int  EG_getGlobal( const ego tess, int global, int *ptype,
    #                             int *pindex, double *xyz )

    # int  EG_solidBoolean( const ego src, const ego tool, int oper, 
    #                                ego *model )
    # int  EG_intersection( const ego src, const ego tool, int *nedge, 
    #                                ego **facEdg, ego *model )
    # int  EG_imprintBody( const ego src, int nedge, const ego *facEdg, 
    #                               ego *result )
    # int  EG_filletBody( const ego src, int nedge, const ego *edges, 
    #                              double radius,
    #                              ego *result, int **facemap )
    # int  EG_chamferBody( const ego src, int nedge, const ego *edges, 
    #                               const ego *faces, double dis1, double dis2, 
    #                               ego *result, int **facemap )
    # int  EG_hollowBody( const ego src, int nface, const ego *faces, 
    #                              double offset, int join,
    #                              ego *result, int **facemap )
    # int  EG_extrude( const ego src, double dist, const double *dir, 
    #                                 ego *result )
    # int  EG_rotate( const ego src, double angle, const double *axis, 
    #                                ego *result )
    # int  EG_sweep( const ego src, const ego spine, int mode,
    #                               ego *result )
    # int  EG_loft( int nsec, const ego *secs, int opt, ego *result )
    # int  EG_blend( int nsec, const ego *secs, double *rc1,
    #                         double *rcN, ego *result )
    # int  EG_ruled( int nsec, const ego *secs, ego *result )