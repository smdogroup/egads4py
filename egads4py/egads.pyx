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

# Class types
NIL = EGADS_NIL
EMPTY = EGADS_EMPTY
REFERENCE = EGADS_REFERENCE
PCURVE = EGADS_PCURVE
CURVE = EGADS_CURVE
SURFACE = EGADS_SURFACE
NODE = EGADS_NODE
EDGE = EGADS_EDGE
LOOP = EGADS_LOOP
FACE = EGADS_FACE
SHELL = EGADS_SHELL
BODY = EGADS_BODY
MODEL = EGADS_MODEL

# /* PCURVES & CURVES */
LINE = EGADS_LINE
CIRCLE = EGADS_CIRCLE
ELLIPSE = EGADS_ELLIPSE
PARABOLA = EGADS_PARABOLA
HYPERBOLA = EGADS_HYPERBOLA
TRIMMED = EGADS_TRIMMED
BEZIER = EGADS_BEZIER
BSPLINE = EGADS_BSPLINE
OFFSET = EGADS_OFFSET

# /* SURFACES */
PLANE = EGADS_PLANE
SPHERICAL = EGADS_SPHERICAL
CYLINDRICAL = EGADS_CYLINDRICAL
REVOLUTION = EGADS_REVOLUTION
TORODIAL = EGADS_TOROIDAL
CONICAL = EGADS_CONICAL
EXTRUSION = EGADS_EXTRUSION

# /* TOPOLOGY */
SREVERSED = EGADS_SREVERSE
NOMTYPE = EGADS_NOMTYPE
SFORWARD = EGADS_SFORWARD
ONENODE = EGADS_ONENODE
TWONODE = EGADS_TWONODE
OPEN = EGADS_OPEN
CLOSED = EGADS_CLOSED
DEGENERATE = EGADS_DEGENERATE
WIREBODY = EGADS_WIREBODY
FACEBODY = EGADS_FACEBODY
SHEETBODY = EGADS_SHEETBODY
SOLIDBODY = EGADS_SOLIDBODY

def _checkErr(int stat):
    # TODO: More descriptive error messages
    errmsg = 'An error occured with code %d'%(stat)
    raise RuntimeError(errmsg)

def revision():
    cdef int major, minor
    cdef const char *occrev
    EG_revision(&major, &minor, &occrev)
    return major, minor, occrev

cdef class pyego:
    cdef ego context
    cdef ego ptr

    def __init__(self, pyego contxt=None):
        cdef int stat
        if contxt:
            self.context = (<pyego>contxt).context
        else:
            stat = EG_open(&self.context)
            if stat:
                _checkErr(stat)
        self.ptr = NULL

    def __dealloc__(self):
        cdef int stat
        if self.ptr:
            stat = EG_deleteObject(self.ptr)
            if stat:
                self._checkErr(stat)

    def setOutLevel(self, int outlevel):
        cdef int stat
        stat = EG_setOutLevel(self.context, outlevel)
        if stat:
            _checkErr(stat)

    def loadModel(self, str filename, int bflag=1):
        cdef int stat
        if self.ptr:
            EG_deleteObject(self.ptr)
        stat = EG_loadModel(self.context, bflag, filename, &self.ptr)
        if stat:
            _checkErr(stat)

    def saveModel(self, str filename):
        cdef int stat
        stat = EG_saveModel(self.ptr, filename)
        if stat:
            _checkErr(stat)

    def makeTransform(self, np.ndarray[double, ndim=1, mode='c'] xform):
        cdef int stat
        new = pyego(self)
        stat = EG_makeTransform(self.context, <double*>xform.data, &new.ptr)
        if stat:
            _checkErr(stat)

    def getTransform(self):
        cdef int stat
        cdef np.ndarray xform = np.zeros(12, np.double)
        stat = EG_getTransformation(self.ptr, <double*>xform.data)
        if stat:
            _checkErr(stat)
        return xform

    def copy(self):
        cdef int stat
        new = pyego(self)
        stat = EG_copyObject(self.ptr, NULL, &new.ptr)
        if stat:
            _checkErr(stat)

    def flip(self):
        cdef int stat
        flip = pyego(self)
        stat = EG_flipObject(self.ptr, &flip.ptr)
        if stat:
            _checkErr(stat)

    def getGeometry(self):
        cdef int stat
        cdef int oclass
        cdef int mtype
        cdef int *ivec
        cdef double *rvec
        cdef ego refgeo
        stat = EG_getGeometry(self.ptr, &oclass, &mtype, &refgeo, 
                              &ivec, &rvec)
        if stat:
            _checkErr(stat)
        return oclass, mtype

    def makeGeometry(self, int oclass, int mtype, pyego refgeo, 
                     np.ndarray[int, ndim=1, mode='c'] ivec,
                     np.ndarray[double, ndim=1, mode='c'] rvec):
        cdef int stat
        stat = EG_makeGeometry(self.context, oclass, mtype, refgeo.ptr,
                               <int*>ivec.data, <double*>rvec.data, 
                               &self.ptr)
        if stat:
            _checkErr(stat)
        return

    def getRange(self):
        cdef int stat
        cdef double r[4]
        cdef int periodic
        cdef int oclass
        if (self.ptr.oclass == EGADS_PCURVE or 
            self.ptr.oclass == EGADS_CURVE or 
            self.ptr.oclass == EGADS_EDGE):
            stat = EG_getRange(self.ptr, r, &periodic)
            if stat:
                _checkErr(stat)
            return [r[0], r[1]], periodic
        elif (self.ptr.oclass == EGADS_SURFACE or 
              self.ptr.oclass == EGADS_EDGE):
            stat = EG_getRange(self.ptr, r, &periodic)
            if stat:
                _checkErr(stat)
            return [r[0], r[1], r[2], r[3]], periodic
        return None

    # def evaluate(self, np.ndarray[double, ndim=1, mode='c'] pt):
    #     cdef double result[3]
    #     stat = EG_evaluate(self.ptr, <double*>pt.data, result)

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

    def getTolerance(self):
        cdef int stat
        cdef double tol
        stat = EG_getTolerance(self.ptr, &tol)
        if stat:
            _checkErr(stat)
        return tol

    def setTolerance(self, double tol):
        cdef int stat
        stat = EG_setTolerance(self.ptr, tol)
        if stat:
            _checkErr(stat)

    # def getTopology(self):
    #     cdef int stat
    #     cdef ego geom
    #     cdef int oclass
    #     cdef int otype
    #     double limits[4]
    #     cdef int nchild
    #     cdef ego **nchild
    #     stat = EG_getTopology(self.ptr, &geom, &oclass, &otype, limits,
    #                           &nchild, &children, &sense)
    #     return geom, oclass, otype, limits, nchild, children, sense
    # int  EG_getTopology( const ego topo, ego *geom, int *oclass, 
    #                               int *type, double *limits, 
    #                               int *nChildren, ego **children, int **sense )

    # def makeTopology(self, pyego geom, int oclass, int otype,
    #                  limits, list children, list sense):

    # int  EG_makeTopology( ego context, ego geom, int oclass,
    #                                int mtype, double *limits,
    #                                int nChildren, ego *children,
    #                                int *senses, ego *topo )

    def makeLoop(self, list edges, pyego geom, double toler):
        cdef int stat 
        cdef int nedges
        cdef ego *edgs
        nedges = len(edges)
        edgs = <ego*>malloc(nedges*sizeof(ego))
        for i in range(nedges):
            edgs[i] = (<pyego>edges[i]).ptr
        new = pyego(self)
        stat = EG_makeLoop(nedges, edgs, geom.ptr, toler, &new.ptr)
        free(edgs)
        if stat:
            _checkErr(stat)
        return new

    # int  EG_getArea( ego object, const double *limits,
    #                  double *area )

    def makeFace(self, pyego obj, int mtype, 
                 np.ndarray[double, ndim=1, mode='c'] limits):
        cdef int stat
        new = pyego(self)
        stat = EG_makeFace(obj.ptr, mtype, <double*>limits.data, &new.ptr)
        if stat:
            _checkErr(stat)
        return new

    # int  EG_getBodyTopos( const ego body, ego src,
    #                                int oclass, int *ntopo, ego **topos )

    def solidBoolean(self, pyego tool, int oper):
        cdef int stat
        new = pyego(self)
        stat = EG_solidBoolean(self.ptr, tool.ptr, oper, &new.ptr)
        if stat:
            _checkErr(stat)

    def intersection(self, pyego tool):
        cdef int stat
        cdef int nedge
        new = pyego(self)
        stat = EG_intersection(self.ptr, tool.ptr, &nedge, NULL, &new.ptr)
        if stat:
            _checkErr(stat)
        return new

    def imprintBody(self, list edges):
        cdef int stat
        cdef nedges = 0
        cdef ego* edgs = NULL
        nedges = len(edges) 
        edgs = <ego*>malloc(nedges*sizeof(ego))
        for i in range(nedges):
            edgs[i] = (<pyego>edges[i]).ptr
        new = pyego(self)
        stat = EG_imprintBody(self.ptr, nedges, edgs, &new.ptr)
        free(edgs)
        if stat:
            _checkErr(stat)
        return new

    def filletBody(self, list edges, double radius):
        cdef int stat
        cdef nedges = 0
        cdef ego* edgs = NULL
        cdef int *facemap
        nedges = len(edges) 
        edgs = <ego*>malloc(nedges*sizeof(ego))
        for i in range(nedges):
            edgs[i] = (<pyego>edges[i]).ptr
        new = pyego(self)
        stat = EG_filletBody(self.ptr, nedges, edgs, radius, &new.ptr,
                             &facemap)
        free(edgs)
        if stat:
            _checkErr(stat)
        return new

    # int  EG_chamferBody( const ego src, int nedge, const ego *edges, 
    #                      const ego *faces, double dis1, double dis2, 
    #                      ego *result, int **facemap )
    # int  EG_hollowBody( const ego src, int nface, const ego *faces, 
    #                     double offset, int join,
    #                     ego *result, int **facemap )
    
    def extrude(self, double dist, 
                np.ndarray[double, ndim=1, mode='c'] _dir):
        cdef int stat
        new = pyego(self)
        stat = EG_extrude(self.ptr, dist, <double*>_dir.data, &new.ptr)
        if stat:
            _checkErr(stat)
        return new

    def rotate(self, double angle, 
               np.ndarray[double, ndim=1, mode='c'] axis):
        cdef int stat
        new = pyego(self)
        stat = EG_rotate(self.ptr, angle, <double*>axis.data, &new.ptr)
        if stat:
            _checkErr(stat)
        return new

    def sweep(self, pyego spline, int mode):
        cdef int stat
        new = pyego(self)
        stat = EG_sweep(self.ptr, spline.ptr, mode, &new.ptr)
        if stat:
            _checkErr(stat)
        return stat

    def blend(self, list sections, 
              np.ndarray[double, ndim=1, mode='c'] rc1,
              np.ndarray[double, ndim=1, mode='c'] rc2):
        cdef int stat
        cdef int nsec
        cdef ego *secs
        nsec = len(sections)
        secs = <ego*>malloc(nsec*sizeof(ego))
        for i in range(nsec):
            secs[i] = (<pyego>sections[i]).ptr
        new = pyego(self)
        stat = EG_blend(nsec, secs, <double*>rc1.data, <double*>rc2.data,
                        &new.ptr)
        if stat:
            _checkErr(stat)
        free(secs)
        return new

    def ruled(self, list sections):
        cdef int stat
        cdef int nsec
        cdef ego *secs
        if self.ptr:
            stat = EG_deleteObject(self.ptr)
            if stat:
                _checkErr(stat)
        nsec = len(sections)
        secs = <ego*>malloc(nsec*sizeof(ego))
        for i in range(nsec):
            secs[i] = (<pyego>sections[i]).ptr
        stat = EG_ruled(nsec, secs, &self.ptr)
        free(secs)
        if stat:
            _checkErr(stat)
        return




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

