# Import numpy 
cimport numpy as np
import numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import the definitions
from egads cimport *

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
TOROIDAL = EGADS_TOROIDAL
CONICAL = EGADS_CONICAL
EXTRUSION = EGADS_EXTRUSION

# /* TOPOLOGY */
SREVERSE = EGADS_SREVERSE
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

# /* ATTRIBUTE TYPES */
ATTRINT = EGADS_ATTRINT
ATTRREAL = EGADS_ATTRREAL
ATTRSTRING = EGADS_ATTRSTRING
ATTRCSYS = EGADS_ATTRCSYS

# /* SOLID BOOLEAN OPERATIONS */
SUBTRACTION = EGADS_SUBTRACTION
INTERSECTION = EGADS_INTERSECTION
FUSION = EGADS_FUSION

# /* SOLID BODY TYPES */
BOX = EGADS_BOX
SPHERE = EGADS_SPHERE
CONE = EGADS_CONE
CYLINDER = EGADS_CYLINDER
TORUS = EGADS_TORUS

# /* ISOCLINE TYPES */
UISO = EGADS_UISO
VISO = EGADS_VISO

# /* FACE SOURCE TYPES */
NODEOFF = EGADS_NODEOFF
EDGEOFF = EGADS_EDGEOFF
FACEDUP = EGADS_FACEDUP
FACECUT = EGADS_FACECUT
FACEOFF = EGADS_FACEOFF

oclass_str = {
    NIL : "EGADS_NIL",
    EMPTY : "EGADS_EMPTY",
    REFERENCE : "EGADS_REFERENCE",
    PCURVE : "EGADS_PCURVE",
    CURVE : "EGADS_CURVE",
    SURFACE : "EGADS_SURFACE",
    NODE : "EGADS_NODE",
    EDGE : "EGADS_EDGE",
    LOOP : "EGADS_LOOP",
    FACE : "EGADS_FACE",
    SHELL : "EGADS_SHELL",
    BODY : "EGADS_BODY",
    MODEL : "EGADS_MODEL" }

oclass_to_type_str = {
    PCURVE : {
        LINE : "EGADS_LINE",
        CIRCLE : "EGADS_CIRCLE",
        ELLIPSE : "EGADS_ELLIPSE",
        PARABOLA : "EGADS_PARABOLA",
        HYPERBOLA : "EGADS_HYPERBOLA",
        TRIMMED : "EGADS_TRIMMED",
        BEZIER : "EGADS_BEZIER",
        BSPLINE : "EGADS_BSPLINE",
        OFFSET : "EGADS_OFFSET" },
    SURFACE : {
        PLANE : "EGADS_PLANE",
        SPHERICAL : "EGADS_SPHERICAL",
        CYLINDRICAL : "EGADS_CYLINDRICAL",
        REVOLUTION : "EGADS_REVOLUTION",
        TOROIDAL : "EGADS_TOROIDAL",
        CONICAL : "EGADS_CONICAL",
        EXTRUSION : "EGADS_EXTRUSION" } }

# Error codes
EGADS_TESSTATE = _EGADS_TESSTATE
EGADS_EXISTS = _EGADS_EXISTS
EGADS_ATTRERR = _EGADS_ATTRERR
EGADS_TOPOCNT = _EGADS_TOPOCNT
EGADS_OCSEGFLT = _EGADS_OCSEGFLT
EGADS_BADSCALE = _EGADS_BADSCALE
EGADS_NOTORTHO = _EGADS_NOTORTHO
EGADS_DEGEN = _EGADS_DEGEN
EGADS_CONSTERR = _EGADS_CONSTERR
EGADS_TOPOERR = _EGADS_TOPOERR
EGADS_GEOMERR = _EGADS_GEOMERR
EGADS_NOTBODY = _EGADS_NOTBODY
EGADS_WRITERR = _EGADS_WRITERR
EGADS_NOTMODEL = _EGADS_NOTMODEL
EGADS_NOLOAD = _EGADS_NOLOAD
EGADS_RANGERR = _EGADS_RANGERR
EGADS_NOTGEOM = _EGADS_NOTGEOM
EGADS_NOTTESS = _EGADS_NOTTESS
EGADS_EMPTYERR = _EGADS_EMPTY
EGADS_NOTTOPO = _EGADS_NOTTOPO
EGADS_REFERCE = _EGADS_REFERCE
EGADS_NOTXFORM = _EGADS_NOTXFORM
EGADS_NOTCNTX = _EGADS_NOTCNTX
EGADS_MIXCNTX = _EGADS_MIXCNTX
EGADS_NODATA = _EGADS_NODATA
EGADS_NONAME = _EGADS_NONAME
EGADS_INDEXERR = _EGADS_INDEXERR
EGADS_MALLOC = _EGADS_MALLOC
EGADS_NOTOBJ = _EGADS_NOTOBJ
EGADS_NULLOBJ = _EGADS_NULLOBJ
EGADS_NOTFOUND = _EGADS_NOTFOUND
EGADS_SUCCESS = _EGADS_SUCCESS
EGADS_OUTSIDE = _EGADS_OUTSIDE

egads_error_codes = {
    EGADS_TESSTATE : "EGADS_TESSTATE",
    EGADS_EXISTS : "EGADS_EXISTS",
    EGADS_ATTRERR : "EGADS_ATTRERR",
    EGADS_TOPOCNT : "EGADS_TOPOCNT",
    EGADS_OCSEGFLT : "EGADS_OCSEGFLT",
    EGADS_BADSCALE : "EGADS_BADSCALE",
    EGADS_NOTORTHO : "EGADS_NOTORTHO",
    EGADS_DEGEN : "EGADS_DEGEN",
    EGADS_CONSTERR : "EGADS_CONSTERR",
    EGADS_TOPOERR : "EGADS_TOPOERR",
    EGADS_GEOMERR : "EGADS_GEOMERR",
    EGADS_NOTBODY : "EGADS_NOTBODY",
    EGADS_WRITERR : "EGADS_WRITERR",
    EGADS_NOTMODEL : "EGADS_NOTMODEL",
    EGADS_NOLOAD : "EGADS_NOLOAD",
    EGADS_RANGERR : "EGADS_RANGERR",
    EGADS_NOTGEOM : "EGADS_NOTGEOM",
    EGADS_NOTTESS : "EGADS_NOTTESS",
    EGADS_EMPTYERR : "EGADS_EMPTY",
    EGADS_NOTTOPO : "EGADS_NOTTOPO",
    EGADS_REFERCE : "EGADS_REFERCE",
    EGADS_NOTXFORM : "EGADS_NOTXFORM",
    EGADS_NOTCNTX : "EGADS_NOTCNTX",
    EGADS_MIXCNTX : "EGADS_MIXCNTX",
    EGADS_NODATA : "EGADS_NODATA",
    EGADS_NONAME : "EGADS_NONAME",
    EGADS_INDEXERR : "EGADS_INDEXERR",
    EGADS_MALLOC : "EGADS_MALLOC",
    EGADS_NOTOBJ : "EGADS_NOTOBJ",
    EGADS_NULLOBJ : "EGADS_NULLOBJ",
    EGADS_NOTFOUND : "EGADS_NOTFOUND",
    EGADS_SUCCESS : "EGADS_SUCCESS",
    EGADS_OUTSIDE : "EGADS_OUTSIDE" }

def _checkErr(int stat):
    # TODO: More descriptive error messages
    errmsg = 'Egads returned error code: %s'%(egads_error_codes[stat])
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
            self.context = contxt.context
        else:
            stat = EG_open(&self.context)
            if stat:
                _checkErr(stat)
        self.ptr = NULL

    def __dealloc__(self):
        # if self.ptr:
        #     EG_deleteObject(self.ptr)
        pass

    def setOutLevel(self, int outlevel):
        EG_setOutLevel(self.context, outlevel)

    def getInfo(self):
        cdef int stat
        cdef int oclass
        cdef int mtype
        cdef ego topObj
        cdef ego prev
        cdef ego next
        stat = EG_getInfo(self.ptr, &oclass, &mtype, &topObj, &prev, &next)
        if stat:
            _checkErr(stat)
        return oclass, mtype

    def getInfoStr(self):
        oclass, mtype = self.getInfo()
        return oclass_str[oclass]
        
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
        cdef int *ivec = NULL
        cdef double *rvec = NULL
        cdef ego refgeo
        stat = EG_getGeometry(self.ptr, &oclass, &mtype, &refgeo, 
                              &ivec, &rvec)
        if stat:
            _checkErr(stat)
        return oclass, mtype

    def makeGeometry(self, int oclass, int mtype,
                     idata=None, rdata=None,                   
                     pyego refgeo=None):
        cdef int stat
        cdef ego rgeo = NULL
        cdef int *iptr = NULL
        cdef double *rptr = NULL
        if self.ptr:
            stat = EG_deleteObject(self.ptr)
            if stat:
                _checkErr(stat)
        if not (oclass == CURVE or oclass == PCURVE or oclass == SURFACE):
            errmsg = 'makeGeometry only accepts CURVE, PCURVE or SURFACE'
            raise ValueError(errmsg)
        if refgeo is not None:
            rgeo = refgeo.ptr
        ivec, rvec = [], []
        try:
            if idata is not None:
                for item in idata:
                    if isinstance(item, int):
                        ivec.append(item)
                    else:
                        ivec.extend(item)
        except:
            errmsg = 'Failed to convert integer data'
            raise ValueError(errmsg)

        try:
            if rdata is not None:
                for item in rdata:
                    if isinstance(item, (int, float)):
                        rvec.append(item)
                    else:
                        rvec.extend(item)
        except:
            errmsg = 'Failed to convert real data'
            raise ValueError(errmsg)

        if len(ivec) > 0:
            iptr = <int*>malloc(len(ivec)*sizeof(int))
            for i in range(len(ivec)):
                iptr[i] = ivec[i]
        if len(rvec) > 0:
            rptr = <double*>malloc(len(rvec)*sizeof(double))
            for i in range(len(rvec)):
                rptr[i] = rvec[i]       
        stat = EG_makeGeometry(self.context, oclass, mtype, rgeo,
                               iptr, rptr, &self.ptr)
        if stat:
            _checkErr(stat)
        if len(ivec) > 0:
            free(iptr)
        if len(rvec) > 0:
            free(rptr)
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

    def evaluate(self, p):
        cdef int stat
        cdef double param[2]
        cdef double r[18]
        try:
            if isinstance(p, (int, float)):
                param[0] = p
            else:
                param[0] = p[0]
                param[1] = p[1]
        except:
            errmsg = 'Failed to convert parameter value'
            raise ValueError(errmsg)
        stat = EG_evaluate(self.ptr, param, r)
        if stat:
            _checkErr(stat)
        if self.ptr.oclass == EGADS_PCURVE:
            return [r[0], r[1]], [r[2], r[3]], [r[4], r[5]]
        elif (self.ptr.oclass == EGADS_EDGE or
              self.ptr.oclass == EGADS_CURVE):
            return [r[0], r[1], r[2]], [r[3], r[4], r[5]], [r[6], r[7], r[8]]
        elif (self.ptr.oclass == EGADS_FACE or
              self.ptr.oclass == EGADS_SURFACE):
            return [r[0], r[1], r[2]], \
                   [r[3], r[4], r[5], r[6], r[7], r[8]], \
                   [r[9], r[10], r[11],
                    r[12], r[13], r[14],
                    r[15], r[16], r[17]]
        return None

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

    def getBody(self):
        cdef int stat
        body = pyego(self)
        stat = EG_getBody(self.ptr, &body.ptr)
        if stat:
            _checkErr(stat)
        return body

    def getTopology(self):
        cdef int stat
        cdef ego geom = NULL
        cdef int oclass
        cdef int mtype
        cdef double limits[4]
        cdef int nchildren
        cdef ego *childarray
        cdef int *senses
        stat = EG_getTopology(self.ptr, &geom, &oclass, &mtype,
                              limits, &nchildren, &childarray, &senses)
        if stat:
            _checkErr(stat)
        children = []
        for i in range(nchildren):
            c = pyego(self)
            c.ptr = childarray[i]
            children.append(c)

        sens = None
        if oclass == LOOP or oclass == SHELL:
            sens = []
            for i in range(nchildren):
                sens.append(senses[i])
        lim = []
        if oclass == NODE:
            lim = [limits[0], limits[1], limits[2]]
        elif oclass == EDGE:
            lim = [limits[0], limits[1]]
        elif oclass == FACE:
            lim = [limits[0], limits[1], limits[2], limits[3]]
        geo = pyego(self)
        geo.ptr = geom
        return geo, oclass, mtype, lim, children, sens

    def makeTopology(self, int oclass, int mtype=0,
                     children=None, sens=None, rdata=None, refgeo=True):
        cdef int stat
        cdef ego refptr = NULL
        cdef double data[4]
        cdef int nchildren = 0
        cdef ego *childarray = NULL
        cdef int *senses = NULL
        errmsg = None
        if oclass == EDGE:
            if mtype != TWONODE and mtype != ONENODE and mtype != CLOSED:
                errmsg = 'EDGE must be TWONODE, ONENODE or CLOSED'
        elif oclass == LOOP:
            if mtype != OPEN and mtype != CLOSED:
                errmsg = 'LOOP must be OPEN or CLOSED'
        elif oclass == FACE:
            if mtype != SFORWARD and mtype != SREVERSE:
                errmsg = 'FACE must be SFORWARD or SREVERSE'
        elif oclass == SHELL:
            if mtype != OPEN and mtype != CLOSED:
                errmsg = 'SHELL must be OPEN or CLOSED'
        elif oclass == BODY:
            if (mtype != WIREBODY and mtype != FACEBODY and
                mtype != SHEETBODY and mtype != SOLIDBODY):
                errmsg = 'BODY must be WIREBODY, FACEBODY'
                errmsg += ', SHEETBODY or SOLIDBODY'
        if errmsg is not None:
            raise ValueError(errmsg)
        
        if ((oclass == SHELL or oclass == LOOP) and
            len(children) != len(sens)):
            errmsg = 'Children and senses list must be of equal length'
            raise ValueError(errmsg)

        if oclass == NODE:
            data[0] = rdata[0]
            data[1] = rdata[1]
            data[2] = rdata[2]
        elif oclass == EDGE:
            data[0] = rdata[0]
            data[1] = rdata[1]
        elif oclass == FACE:
            data[0] = rdata[0]
            data[1] = rdata[1]
            data[2] = rdata[2]
            data[3] = rdata[3]

        if children is not None:
            nchildren = len(children)
            childarray = <ego*>malloc(nchildren*sizeof(ego))
            for i in range(nchildren):
                childarray[i] = (<pyego>children[i]).ptr

        if oclass == FACE or oclass == LOOP:
            senses = <int*>malloc(nchildren*sizeof(int))
            for i in range(nchildren):
                senses[i] = sens[i]

        if refgeo:
            refptr = self.ptr
        else:
            refptr = NULL
            
        new = pyego(self)
        stat = EG_makeTopology(self.context, refptr, oclass, mtype,
                               data, nchildren, childarray,
                               senses, &new.ptr)
        if stat:
            _checkErr(stat)
        free(childarray)
        free(senses)
        return new

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

    def makeFace(self, pyego obj, int mtype, 
                 np.ndarray[double, ndim=1, mode='c'] limits):
        cdef int stat
        new = pyego(self)
        stat = EG_makeFace(obj.ptr, mtype, <double*>limits.data, &new.ptr)
        if stat:
            _checkErr(stat)
        return new

    def makeSolidBody(self, int stype, rdata=None):
        cdef int stat
        cdef double *rptr = NULL
        rvec = []
        try:
            if rdata is not None:
                for item in rdata:
                    if isinstance(item, (int, float)):
                        rvec.append(item)
                    else:
                        rvec.extend(item)
        except:
            errmsg = 'Failed to convert real data'
            raise ValueError(errmsg)
        
        if len(rvec) > 0:
            rptr = <double*>malloc(len(rvec)*sizeof(double))
            for i in range(len(rvec)):
                rptr[i] = rvec[i]

        body = pyego(self)
        stat = EG_makeSolidBody(self.context, stype, rptr, &body.ptr)
        if rptr:
            free(rptr)
        return body

    def sewFaces(self, list objlist, double toler=0.0, int flag=0):
        cdef ego context
        cdef int nobj
        cdef ego *objs

        if len(objlist) > 0:
            nobj = len(objlist)
            objs = <ego*>malloc(nobj*sizeof(ego))
            for i in range(nobj):
                objs[i] = (<pyego>objlist[i]).ptr
            new = pyego(self)
            stat = EG_sewFaces(nobj, objs, toler, flag, &new.ptr)
            if stat:
                _checkErr(stat)
            free(objs)
        return new

    def getBodyTopos(self, int oclass, pyego ref=None):
        cdef int stat
        cdef ego src = NULL
        cdef int ntopos = 0
        cdef ego *topos = NULL
        if ref is not None:
            src = ref.ptr
        if not (oclass == NODE or oclass == EDGE or oclass == LOOP or
                oclass == FACE or oclass == SHELL):
            errmsg = 'Object must be either NODE, EDGE, LOOP, FACE or SHELL'
            raise ValueError(errmsg)
        stat = EG_getBodyTopos(self.ptr, src, oclass, &ntopos, &topos)
        if stat:
            _checkErr(stat)
        tlist = []
        for i in range(ntopos):
            t = pyego(self)
            t.ptr = topos[i]
            tlist.append(t)
        free(topos)
        return tlist

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
