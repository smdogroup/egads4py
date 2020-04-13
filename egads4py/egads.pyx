# Import the definition required for const strings
from libc.string cimport const_char, strcpy, strlen
from libc.stdlib cimport malloc, free

# Import numpy
cimport numpy as np

# Import the definitions
from egads cimport *

# Import os for file I/O
import os

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

def _convert_rdata(rdata):
    rvec = []
    try:
        for item in rdata:
            if isinstance(item, (int, float)):
                rvec.append(item)
            elif isinstance(item, np.ndarray):
                rvec.extend(item.flatten())
            else:
                rvec.extend(item)
    except:
        errmsg = 'Failed to convert real data'
        raise ValueError(errmsg)
    return rvec

cdef class context:
    def __init__(self):
        cdef int stat
        stat = EG_open(&self.context)
        if stat:
            _checkErr(stat)
        self.refs = []
        return

    def __dealloc__(self):
        cdef int stat
        self.refs = None
        stat = EG_close(self.context)
        if stat:
            _checkErr(stat)
        return

    def setOutLevel(self, int outlevel):
        '''
        Set output level

        Parameters
        ----------
        outlevel: 0 <= outlevel <= 2
        '''
        EG_setOutLevel(self.context, outlevel)

    def makeTransform(self, xform):
        '''
        Creates a TRANSFORM object from the 12 values. The rotation
        portion [3][3] must be “scaled” orthonormal (orthogonal with a
        single scale factor).
        '''
        cdef int stat
        cdef double T[12]
        for i in range(12):
            T[i] = xform[i]
        new_obj = pyego(self)
        stat = EG_makeTransform(self.context, T, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj

    def makeGeometry(self, int oclass, int mtype,
                     pyego geom=None, idata=None, rdata=None):
        '''
        Creates a geometric object:

        Parameters
        ----------
        oclass: PCURVE, CURVE or SURFACE

        mtype:
        For PCURVE/CURVE:
        LINE, CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA, TRIMMED,
        BEZIER, BSPLINE, OFFSET

        For SURFACE:
        PLANE, SPHERICAL, CYLINDRICAL, REVOLUTION, TORIODAL,
        TRIMMED, BEZIER, BSPLINE, OFFSET, CONICAL, EXTRUSION

        geom:
        The reference geometry object (if none use None)

        idata:
        is a pointer to the block of integer information. Required for
        either BEZIER or BSPLINE.

        rdata:
        is the pointer to a block of double precision reals. The
        content and length depends on the oclass/mtype.

        - For mtype=LINE: rdata=[x0, dx], where x0 is the starting
          point of the line, and dx is the distance of the line in each
          direction
        - For mtype=PARABOLA: rdata=[O, XDir, YDir, F]
          The parametric equation takes the form:
          P(U) = O + U*U/(4.*F)*XDir + U*YDir
          where: P is the point of parameter U,
          O, XDir and YDir are respectively the origin, "X Direction"
          and "Y Direction" of its local coordinate system,
          F is the focal length of the parabola. Note that Xdir defines
          the axis of symmetry, not the standard "x-axis".

        returns:
        Resultant new geometry object
        '''
        cdef int stat
        cdef ego refptr = NULL
        cdef int *iptr = NULL
        cdef double *rptr = NULL

        if geom is not None:
            refptr = geom.ptr

        if not (oclass == CURVE or oclass == PCURVE or oclass == SURFACE):
            errmsg = 'makeGeometry only accepts CURVE, PCURVE or SURFACE'
            raise ValueError(errmsg)

        ivec, rvec = [], []
        try:
            if idata is not None:
                for item in idata:
                    if isinstance(item, int):
                        ivec.append(item)
                    elif isinstance(item, np.ndarray):
                        ivec.extend(item.flatten())
                    else:
                        ivec.extend(item)
        except:
            errmsg = 'Failed to convert integer data'
            raise ValueError(errmsg)

        if rdata is not None:
            rvec = _convert_rdata(rdata)

        if len(ivec) > 0:
            iptr = <int*>malloc(len(ivec)*sizeof(int))
            for i in range(len(ivec)):
                iptr[i] = ivec[i]
        if len(rvec) > 0:
            rptr = <double*>malloc(len(rvec)*sizeof(double))
            for i in range(len(rvec)):
                rptr[i] = rvec[i]
        new_obj = pyego(self)
        stat = EG_makeGeometry(self.context, oclass, mtype, refptr,
                               iptr, rptr, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        if len(ivec) > 0:
            free(iptr)
        if len(rvec) > 0:
            free(rptr)
        return new_obj

    def makeTopology(self, int oclass, int mtype=0, pyego geom=None,
                     children=None, sens=None, rdata=None):
        '''
        Creates and returns a topological object:

        Parameters
        ----------
        oclass:
        either NODE, EGDE, LOOP, FACE, SHELL, BODY or MODEL

        mtype:
        for EDGE is TWONODE, ONENODE or DEGENERATE
        for LOOP is OPEN or CLOSED
        for FACE is either SFORWARD or SREVERSE
        for SHELL is OPEN or CLOSED
        BODY is either WIREBODY, FACEBODY, SHEETBODY or SOLIDBODY

        geom:
        reference geometry object required for EDGEs and FACEs
        (optional for LOOP)

        rdata:
        may be None except for:
        NODE which contains the [x,y,z] location
        EDGE is the t-min and t-max (the parametric bounds)

        children:
        chldrn a list of children objects (nchild in length)
        if LOOP and has reference SURFACE, then 2*nchild in length
        (PCURVES follow)

        senses:
        a list of integer senses for the children (required for FACES
        & LOOPs only)

        returns:
        the resultant returned topological object
        '''
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

        if ((oclass == FACE or oclass == LOOP) and
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

        if geom is not None:
            refptr = geom.ptr

        new_obj = pyego(self)
        stat = EG_makeTopology(self.context, refptr, oclass, mtype,
                               data, nchildren, childarray,
                               senses, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        free(childarray)
        free(senses)
        return new_obj

    def makeLoop(self, list edges, pyego geom=None, double toler=0.0):
        '''
        Creates a LOOP from a list of EDGE Objects, where the EDGEs do
        not have to be topologically connected. The tolerance is used
        to build the NODEs for the LOOP. The orientation is set by the
        first non-NULL entry in the list, which is taken in the
        positive sense. This is designed to be executed until all list
        entries are exhausted.

        Parameters
        ----------
        edges:
        list of EDGEs, of which some may be NULL (nEdge in length)
        Note: list entries are None-ified when included in LOOPs

        geom:
        SURFACE Object for non-planar LOOPs to be used to bound
        FACEs (can be NULL)

        toler:
        tolerance used for the operation (0.0 - use EDGE tolerances)

        returns:
        obj: the resultant LOOP Object
        nedges: the number of non-None entries in edges when returned
        or error code
        '''
        cdef int nloop_edges = 0
        cdef int nedges
        cdef ego *edgs = NULL
        cdef ego geoptr = NULL
        if geom is not None:
            geoptr = geom.ptr
        nedges = len(edges)
        edgs = <ego*>malloc(nedges*sizeof(ego))
        for i in range(nedges):
            if edges[i] is None:
                edgs[i] = NULL
            else:
                edgs[i] = (<pyego>edges[i]).ptr
        new_obj = pyego(self)
        nloop_edges = EG_makeLoop(nedges, edgs, geoptr, toler, &new_obj.ptr)
        for i in range(nedges):
            if edgs[i] == NULL:
                edges[i] = None
        free(edgs)
        return new_obj, nloop_edges

    def makeFace(self, pyego obj, int mtype, rdata=None):
        '''
        Creates a simple FACE from a LOOP or a SURFACE. Also can be
        used to hollow a single LOOPed existing FACE. This function
        screates any required NODEs, EDGEs and LOOPs.

        Parameters
        ----------
        obj:
        Either a LOOP (for a planar cap), a SURFACE with [u,v] bounds,
        or a FACE to be hollowed out

        mtype:
        Is either SFORWARD or SREVERSE
        For LOOPs you may want to look at the orientation using
        getArea, ignored when the input object is a FACE

        rdata:
        May be None for LOOPs, but must be the limits for a SURFACE (4
        values), the hollow/offset distance and fillet radius (zero is
        for no fillets) for a FACE input object (2 values)

        returns:
        The resultant returned topological FACE object (a return of
        EGADS_OUTSIDE is the indication that offset distance was too
        large to produce any cutouts, and this result is the input
        object)
        '''
        cdef int stat
        cdef double data[4]
        cdef double *ptr = NULL
        if obj.ptr.oclass == SURFACE:
            # Limits of the surface
            data[0] = rdata[0]
            data[1] = rdata[1]
            data[2] = rdata[2]
            data[3] = rdata[3]
            ptr = data
        elif obj.ptr.oclass == FACE:
            data[0] = rdata[0]
            data[1] = rdata[1]
            ptr = data

        new_obj = pyego(self)
        stat = EG_makeFace(obj.ptr, mtype, ptr, &new_obj.ptr)
        if stat and stat != EGADS_OUTSIDE:
            _checkErr(stat)
        return new_obj

    def makeSolidBody(self, int stype, rdata=None):
        '''
        Creates a simple SOLIDBODY. Can be either a box, cylinder,
        sphere, cone, or torus.

        Parameters
        ----------

        stype: BOX, SPHERE, CONE, CYLINDER or TORUS

        rdata:
        Depends on stype
        For box [x,y,z] then [dx,dy,dz] for size of box
        For sphere [x,y,z] of center then radius
        For cone apex [x,y,z], base center [x,y,z], then radius
        For cylinder 2 axis points and the radius
        For torus [x,y,z] of center, direction of rotation, then major
        radius and minor radius

        returns:
        the resultant returned topological BODY object
        '''
        cdef int stat
        cdef double *rptr = NULL
        rvec = []
        if rdata is not None:
            rvec = _convert_rdata(rdata)

        if len(rvec) > 0:
            rptr = <double*>malloc(len(rvec)*sizeof(double))
            for i in range(len(rvec)):
                rptr[i] = rvec[i]

        body = pyego(self)
        stat = EG_makeSolidBody(self.context, stype, rptr, &body.ptr)
        if rptr:
            free(rptr)
        return body

    def sewFaces(self, list objlist, double toler=0.0, manifold=True):
        '''
        Creates a MODEL from a collection of Objects. The Objects can
        be either BODYs (not WIREBODY), SHELLs and/or FACEs. After the
        sewing operation, any unconnected Objects are returned as
        BODYs.

        Parameters
        ----------

        objlist:
        List of Objects to sew together

        toler:
        Tolerance used for the operation (0.0 - use Face tolerances)

        manifold:
        Indicates whether to produce manifold/non-manifold geometry

        returns:
        The resultant MODEL object
        '''
        cdef int nobj
        cdef ego *objs
        cdef int flag = 1
        if manifold:
            flag = 0

        if len(objlist) > 0:
            nobj = len(objlist)
            objs = <ego*>malloc(nobj*sizeof(ego))
            for i in range(nobj):
                objs[i] = (<pyego>objlist[i]).ptr
            new_obj = pyego(self)
            stat = EG_sewFaces(nobj, objs, toler, flag, &new_obj.ptr)
            if stat:
                _checkErr(stat)
            free(objs)
            return new_obj
        return None

    def loadModel(self, fname, split=False):
        '''
        Loads a MODEL object from a file

        Parameters
        ----------
        fname:     the file name to load
        split:     split closed/periodic entries

        returns:   the MODEL object
        '''
        cdef int stat
        cdef int bflag = 1
        cdef char *filename = NULL
        if split:
            bflag = 2
        new_obj = pyego(self)
        filename = egads_convert_to_chars(fname)
        stat = EG_loadModel(self.context, bflag, filename, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj

    def blend(self, list sections, list _rc1=None, list _rc2=None):
        '''
        Simply lofts the input Objects to create a BODY Object (that
        has the type SOLIDBODY or SHEETBODY). Cubic BSplines are
        used. All sections must have the same number of Edges (except
        for NODEs) and the Edge order in each (defined in a CCW
        manner) is used to specify the loft connectivity .

        Parameters
        ----------
        sections:
        list of WIREBODY or LOOP objects to Blend - nSection in len
        the first and last can be NODEs and/or FACEs (only one LOOP),
        if the first and last are NODEs and/or FACEs (and the
        intermediate sections are CLOSED) the result will be a
        SOLIDBODY otherwise a SHEETBODY will be constructed interior
        sections can be repeated once for C1 or twice for C0

        rc1: specifies treatment* at the first section (or None for no
        treatment)

        rc2: specifies treatment* at the last section (or None for no
        treatment)

        Returns
        -------
        the resultant BODY object

        * for NODEs -- elliptical treatment (8 in length): radius of
        curvature1, unit direction, rc2, orthogonal direction;
        nSection must be at least 3 (or 4 for treatments at both ends)
        for other sections -- setting tangency (4 in length):
        magnitude, unit direction for FACEs with 2 or 3 EDGEs -- make
        a Wing Tip-like cap: zero, growthFactor (len of 2)
        '''
        cdef int stat
        cdef int nsec
        cdef ego *secs
        cdef double rc1[8]
        cdef double rc2[8]
        cdef double *rc1ptr = NULL
        cdef double *rc2ptr = NULL
        nsec = len(sections)
        secs = <ego*>malloc(nsec*sizeof(ego))
        for i in range(nsec):
            secs[i] = (<pyego>sections[i]).ptr
        if _rc1 is not None:
            for i in range(min(len(_rc1), 8)):
                rc1[i] = _rc1[i]
            rc1ptr = rc1
        if _rc2 is not None:
            for i in range(min(len(_rc2), 8)):
                rc2[i] = _rc2[i]
            rc2ptr = rc2
        new_obj = pyego(self)
        stat = EG_blend(nsec, secs, rc1ptr, rc2ptr, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        free(secs)
        return new_obj

    def ruled(self, list sections):
        '''
        Produces a BODY Object (that has the type SOLIDBODY or SHEETBODY)
        that goes through the sections by ruled surfaces between each. All
        sections must have the same number of Edges (except for NODEs) and
        the Edge order in each is used to specify the connectivity.

        Parameters
        ----------
        sections:
        A list of NODE, WIREBODY, LOOP and/or FACE objects to operate
        upon. Any FACE objects must contain only a single LOOP, Only
        the first and last sections can be NODEs, If the first and
        last sections are NODEs and/or FACEs and all WIREBODY and LOOP
        objects are closed, the result will be a SOLIDBODY otherwise a
        SHEETBODY will be constructed result the resultant BODY object

        Note: for both blend and ruled all Loops must have their Edges
        ordered in a counterclockwise manner.
        '''
        cdef int stat
        cdef int nsec
        cdef ego *secs
        nsec = len(sections)
        secs = <ego*>malloc(nsec*sizeof(ego))
        for i in range(nsec):
            secs[i] = (<pyego>sections[i]).ptr
        new_obj = pyego(self)
        stat = EG_ruled(nsec, secs, &new_obj.ptr)
        free(secs)
        if stat:
            _checkErr(stat)
        return new_obj

    def approximate(self, xyz, int dim1=-1, int maxdeg=3,
                    double tol=0.0):
        '''
        Computes and returns the resultant geometry object created by
        approximating the data by a BSpline (OCC or EGADS method).

        Parameters
        ----------
        maxdeg:
        the maximum degree used by OCC [3-8], or cubic by EGADS [0-2]
        0 – fixes the bounds and uses natural end conditions
        1 – fixes the bounds and maintains the slope input at the bounds
        2 – fixes the bounds & quadratically maintains the slope at 2 nd order

        tol:
        is the tolerance to use for the BSpline approximation procedure,
        zero for a SURFACE fit (OCC).

        sizes:
        a vector of 2 integers that specifies the size and dimensionality of
        the data. If the second is zero, then a CURVE is fit and the first
        integer is the length of the number of [x,y,z] triads. If the second
        integer is nonzero then the input data reflects a 2D map.

        xyz:
        the data to fit (3 times the number of points in length)

        Returns
        -------
        the returned approximated (or fit) BSpline resultant object
        '''
        cdef int stat
        cdef int sizes[2]
        cdef double *xyz_array = NULL

        if dim1 <= 0:
            sizes[0] = len(xyz)
            sizes[1] = 0
            length = sizes[0]
        else:
            sizes[0] = dim1
            sizes[1] = len(xyz)//dim1
            length = sizes[0]*sizes[1]

        # Allocate the points
        xyz_array = <double*>malloc(3*length*sizeof(double))
        for i in range(length):
            xyz_array[3*i] = xyz[i][0]
            xyz_array[3*i+1] = xyz[i][1]
            xyz_array[3*i+2] = xyz[i][2]

        new_obj = pyego(self)
        stat = EG_approximate(self.context, maxdeg, tol, sizes,
                              xyz_array, &new_obj.ptr)
        free(xyz_array)
        if stat:
            _checkErr(stat)
        return new_obj

cdef class pyego:
    def __cinit__(self, context ctx):
        self.ctx = ctx
        self.ctx.refs.append(self)
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr:
            EG_deleteObject(self.ptr)
        return

    def isSame(self, pyego obj):
        '''
        Compares two objects for geometric equivalence.

        Parameters
        ----------
        obj: an object to compare with
        '''
        cdef int stat = 0
        stat = EG_isSame(self.ptr, obj.ptr)
        if stat == EGADS_SUCCESS:
            return True
        return False

    def getInfo(self):
        '''
        Return information about the object

        Returns
        -------
        oclass:  object class type
        mtype:   object sub-type
        '''
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
        '''Get a string representation of the class type'''
        oclass, mtype = self.getInfo()
        return oclass_str[oclass]

    def saveModel(self, fname, overwrite=False):
        '''
        Saves the model based on the filename extension
        '''
        cdef int stat
        cdef char* filename = egads_convert_to_chars(fname)
        if overwrite and os.path.exists(filename):
            os.remove(filename)
        stat = EG_saveModel(self.ptr, filename)
        if stat:
            _checkErr(stat)

    def getTransform(self):
        '''
        Returns the transformation information. This appears like is a
        column- major matrix that is 4 columns by 3 rows and could be
        thought of as [3][4] in C (though is flat) and in FORTRAN
        dimensioned as (4,3).
        '''
        cdef int stat
        cdef double T[12]
        stat = EG_getTransformation(self.ptr, T)
        xform = []
        for i in range(12):
            xform.append(T[i])
        if stat:
            _checkErr(stat)
        return xform

    def copy(self):
        '''
        Creates a new EGADS object by copying.
        '''
        cdef int stat
        new_obj = pyego(self.ctx)
        stat = EG_copyObject(self.ptr, NULL, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj

    def flip(self):
        '''
        Creates a new EGADS object by copying and reversing the input
        object.

        Parameters
        ----------
        returns:
        the input object: 3D geometry (flip the parameterization) or
        topology (reverse the sense). Not for NODE, BODY or
        MODEL. SURFACEs reverse only the u parameter.
        '''
        cdef int stat
        flip = pyego(self.ctx)
        stat = EG_flipObject(self.ptr, &flip.ptr)
        if stat:
            _checkErr(stat)
        return flip

    def getGeometry(self):
        '''
        Returns information about the geometric object:

        Returns
        -------
        oclass PCURVE, CURVE or SURFACE
        mtype PCURVE/CURVE
        . LINE, CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA, TRIMMED,
        BEZIER, BSPLINE, OFFSET
        SURFACE
        . PLANE, SPHERICAL, CYLINDRICAL, REVOLUTION, TORIODAL,
        TRIMMED, BEZIER, BSPLINE, OFFSET, CONICAL, EXTRUSION
        '''
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

    def getRange(self):
        '''
        Returns the valid range of the object: may be one of PCURVE,
        CURVE, SURFACE, EDGE or FACE

        Returns
        -------

        range:
        for PCURVE, CURVE or EDGE returns 2 values, t-start and t-end
        for SURFACE or FACE returns 4 values, u-min, u-max, v-min and v-max

        periodic:
        0 for non-periodic, 1 for periodic in t or u 2 for periodic in
        v (or-able)
        '''
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

    def getBoundingBox(self):
        '''
        Computes the Cartesian bounding box around the object:

        Returns
        -------
        the [x,y,z] min and [x,y,z] max
        '''
        cdef int stat
        cdef double data[6]
        stat = EG_getBoundingBox(self.ptr, data)
        if stat:
            _checkErr(stat)
        return [data[0], data[1], data[2]], [data[3], data[4], data[5]]

    def evaluate(self, p):
        '''
        Returns the result of evaluating the object at a parameter
        point. May be used for PCURVE, CURVE, SURFACE, EDGE or FACE.

        Parameters
        ----------
        p: The parametric location
        For PCURVE, CURVE, EDGE: The t-location
        For SURFACE, FACE: The (u, v) coordiantes

        Returns:
        For PCURVE, CURVE, EDGE: X, X,t and X,tt
        For SURFACE, FACE: X, [X,u, X,v], [X,uu, X,uv, X,vv]
        '''
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

    def invEvaluate(self, xyz):
        '''
        Returns the result of inverse evaluation on the object. For
        topology the result is limited to inside the EGDE/FACE valid
        bounds. The object may be one of PCURVE, CURVE, SURFACE, EDGE
        or FACE.

        Parameters
        ----------
        xyz: [u,v] for a PCURVE and [x,y,z] for all others

        Returns
        -------
        for PCURVE, CURVE or EDGE the one value is t
        for SURFACE or FACE the 2 values are u then v the closest
        position found is returned:
        [u,v] for a PCURVE (2) and [x,y,z] for all others (3)
        '''
        cdef int stat
        cdef double xp[3]
        cdef double params[2]
        cdef double result[3]

        if self.ptr.oclass == EGADS_PCURVE:
            xp[0] = xyz[0]
            xp[1] = xyz[1]
        else:
            xp[0] = xyz[0]
            xp[1] = xyz[1]
            xp[2] = xyz[2]

        stat = EG_invEvaluate(self.ptr, xp, params, result)
        if stat:
            _checkErr(stat)
        if self.ptr.oclass == EGADS_PCURVE:
            return params[0], [result[0], result[1]]
        elif (self.ptr.oclass == EGADS_EDGE or
              self.ptr.oclass == EGADS_CURVE):
            return params[0], [result[0], result[1], result[2]]
        elif (self.ptr.oclass == EGADS_FACE or
              self.ptr.oclass == EGADS_SURFACE):
            return [params[0], params[1]], [result[0], result[1], result[2]]
        return None, None

    def getTolerance(self):
        '''
        Get the tolerance associated with this object
        '''
        cdef int stat
        cdef double tol
        stat = EG_getTolerance(self.ptr, &tol)
        if stat:
            _checkErr(stat)
        return tol

    def setTolerance(self, double tol):
        '''
        Set the tolerance associated with this object
        '''
        cdef int stat
        stat = EG_setTolerance(self.ptr, tol)
        if stat:
            _checkErr(stat)

    def getBody(self):
        '''
        Get the body that this object is contained within
        '''
        cdef int stat
        body = pyego(self.ctx)
        stat = EG_getBody(self.ptr, &body.ptr)
        if stat:
            _checkErr(stat)
        return body

    def getArea(self, limits=None):
        '''
        Get the area of the geometric object
        '''
        cdef int stat
        cdef double lim[4]
        cdef double area
        if self.ptr.oclass == SURFACE:
            lim[0] = limits[0]
            lim[1] = limits[1]
            lim[2] = limits[2]
            lim[3] = limits[3]
        stat = EG_getArea(self.ptr, lim, &area)
        return area

    def attributeAdd(self, aname, int atype, data):
        '''
        Adds an attribute to the object. If an attribute exists with
        the name it is overwritten with the new information.

        Parameters
        ----------

        name:
        The name of the attribute. Must not contain a space or other
        special characters

        atype
        atype must be either:
        ATTRINT for integers
        ATTRREAL for double precision
        ATTRSTRING for a character string
        ATTRCSYS for a coordinate system (use reals for input)
        '''
        cdef int stat
        cdef int length = 0
        cdef int *ints = NULL
        cdef double *reals = NULL
        cdef const char *temp = NULL
        cdef char *chars = NULL
        cdef char *name = NULL
        temp = egads_convert_to_chars(aname)
        name = <char*>malloc(strlen(temp)+1)
        strcpy(name, temp)

        if atype == EGADS_ATTRINT:
            length = len(data)
            ints = <int*>malloc(length*sizeof(int))
            for i in range(length):
                ints[i] = <int>data[i]
        elif atype == EGADS_ATTRREAL or atype == EGADS_ATTRCSYS:
            length = len(data)
            reals = <double*>malloc(length*sizeof(double))
            for i in range(length):
                reals[i] = <double>data[i]
        elif atype == EGADS_ATTRSTRING:
            temp = egads_convert_to_chars(data)
            chars = <char*>malloc(strlen(temp)+1)
            strcpy(chars, temp)

        stat = EG_attributeAdd(self.ptr, name, atype, length,
                               ints, reals, chars)
        if stat:
            _checkErr(stat)

        free(name)
        if ints:
            free(ints)
        if reals:
            free(reals)
        if chars:
            free(chars)
        return

    def attributeDel(self, aname=None):
        '''
        Deletes an attribute from the object. If the name is NULL then
        all attributes are removed from this object
        '''
        cdef int stat
        cdef char *name = NULL
        if name is not None:
            name = egads_convert_to_chars(aname)
            stat = EG_attributeDel(self.ptr, name)
        else:
            stat = EG_attributeDel(self.ptr, NULL)
        if stat:
            _checkErr(stat)
        return

    def attributeNum(self):
        '''
        Returns the number of attributes found with this object
        '''
        cdef int stat
        cdef int num
        stat = EG_attributeNum(self.ptr, &num)
        if stat:
            _checkErr(stat)
        return num

    def attributeGet(self, int index):
        '''
        Retrieves a specific attribute from the object
        '''
        cdef int stat
        cdef const char *name
        cdef int atype
        cdef int length
        cdef const int *ints
        cdef const double *reals
        cdef const char *chars

        stat = EG_attributeGet(self.ptr, index, &name, &atype, &length,
                               &ints, &reals, &chars)
        if stat:
            _checkErr(stat)

        if atype == EGADS_ATTRINT:
            res = []
            for i in range(length):
                res.append(ints[i])
            return name, res
        elif atype == EGADS_ATTRREAL or atype == EGADS_ATTRCSYS:
            res = []
            for i in range(length):
                res.append(reals[i])
            return name, res
        elif atype == EGADS_ATTRSTRING:
            return egads_convert_chars_to_str(name),\
                egads_convert_chars_to_str(chars)
        return None

    def attributeRet(self, aname):
        '''
        Retrieves a specific attribute from the object
        '''
        cdef int stat
        cdef int atype
        cdef int length
        cdef const int *ints
        cdef const double *reals
        cdef const char *chars
        cdef char *name = egads_convert_to_chars(aname)

        stat = EG_attributeRet(self.ptr, name, &atype, &length,
                               &ints, &reals, &chars)
        if stat == EGADS_NOTFOUND:
            return None
        elif stat:
            _checkErr(stat)

        if atype == EGADS_ATTRINT:
            res = []
            for i in range(length):
                res.append(ints[i])
            return res
        elif atype == EGADS_ATTRREAL or atype == EGADS_ATTRCSYS:
            res = []
            for i in range(length):
                res.append(reals[i])
            return res
        elif atype == EGADS_ATTRSTRING:
            return egads_convert_chars_to_str(chars)
        return None

    def attributeDup(self, pyego dup):
        '''
        Removes all attributes from the destination object, then copies the
        attributes from the source
        '''
        cdef int stat
        stat = EG_attributeDup(self.ptr, dup.ptr)
        if stat:
            _checkErr(stat)
        return

    def getTopology(self):
        '''
        Returns information about the topological object

        Returns
        -------

        geo:
        The reference geometry object (if none this is returned as None)

        oclass:
        is NODE, EGDE, LOOP, FACE, SHELL, BODY or MODEL

        mtype:
        for EDGE is TWONODE, ONENODE or DEGENERATE
        for LOOP is OPEN or CLOSED
        for FACE is either SFORWARD or SREVERSE
        for SHELL is OPEN or CLOSED
        BODY is either WIREBODY, FACEBODY, SHEETBODY or SOLIDBODY

        lim:
        will retrieve at most 4 doubles:
        for NODE this contains the [x,y,z] location
        EDGE is the t-min and t-max (the parametric bounds)
        FACE returns the [u,v] box (the limits first for u then for v)
        number of children (lesser) topological objects

        children:
        is a returned pointer to the block of children objects.

        sens:
        is the returned pointer to a block of integer senses for
        the children.
        '''
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
            c = pyego(self.ctx)
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
        geo = pyego(self.ctx)
        geo.ptr = geom
        return geo, oclass, mtype, lim, children, sens

    def getChildren(self):
        '''
        Returns a list of the children from the object

        Returns
        -------
        children:
        is a returned pointer to the block of children objects.
        '''
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
            EG_referenceObject(childarray[i], self.ctx.context)
            c = pyego(self.ctx)
            c.ptr = childarray[i]
            children.append(c)
        return children

    def getBodyTopos(self, int oclass, pyego ref=None):
        '''
        Returns topologically connected objects:

        Parameters
        ----------
        oclass:
        is NODE, EGDE, LOOP, FACE or SHELL -- must not be the same
        class as ref

        ref:
        reference topological object or NULL. Sets the context for the
        returned objects (i.e. all objects of a class [oclass] in the
        tree looking towards that class from ref). None starts from
        the BODY (for example all NODEs in the BODY)

        Returns
        -------
        the returned number of requested topological objects is a
        returned pointer to the block of objects
        '''

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
            EG_referenceObject(topos[i], self.ctx.context)
            t = pyego(self.ctx)
            t.ptr = topos[i]
            tlist.append(t)
        EG_free(topos)
        return tlist

    def matchBodyFaces(self, pyego body, double toler):
        '''
        Examines the FACEs in one BODY against all of the FACEs in
        another. If the number of LOOPs, number of NODEs, the NODE
        locations, the number of EDGEs and the EDGE bounding boxes as
        well as the EDGE arc lengths match it is assumed that the
        FACEs match. A list of pairs of indices are returned.

        Parameters
        ----------
        body:
        body container object

        toler:
        the tolerance used (can be zero to use entity tolerances)

        Returns
        -------
        a list of the tuples of matching indices

        Note: This is useful for the situation where there are
        glancing FACEs and a UNION operation fails (or would
        fail). Simply find the matching FACEs and do not include them
        in a call to EG_sewFaces.
        '''
        cdef int stat
        cdef ego body1 = NULL
        cdef ego body2 = NULL
        cdef int nmatches = 0
        cdef int *match = NULL

        matches = None
        if body.ptr.oclass == BODY and self.ptr.oclass == BODY:
            stat = EG_matchBodyFaces(self.ptr, body.ptr, toler,
                                     &nmatches, &match)
            if stat:
                _checkErr(stat)
            if nmatches > 0:
                matches = []
                for i in range(nmatches):
                    matches.append((match[2*i], match[2*i+1]))

                EG_free(match)
        else:
            errmsg = 'Both objects must have object class BODY'
            raise ValueError(errmsg)

        return None

    def indexBodyTopo(self, pyego src):
        '''
        Return the (1-based) index of the topological object in the body
        '''
        return EG_indexBodyTopo(self.ptr, src.ptr)

    def solidBoolean(self, pyego tool, int oper):
        '''
        Performs the Solid Boolean Operations (SBOs) on the source
        BODY Object (that has the type SOLIDBODY). The tool object
        types depend on the operation. This supports Intersection,
        Subtraction and Union. The object must be a SOLIDBODY or
        MODEL.

        Note: This may be called with src being a MODEL. In this case
        tool may be a SOLIDBODY for Intersection or a FACE/FACEBODY
        for Fusion. The input MODEL may contain anything, but must not
        have duplicate topology.

        Parameters
        ----------
        tool: the tool object:
        either a SOLIDBODY for all operators or a FACE/FACEBODY for
        Subtraction.

        oper: the operation to perform

        Returns
        -------
        the resultant MODEL object (this is because there may be
        multiple bodies from either the subtraction or intersection
        operation).
        '''
        cdef int stat
        new_obj = pyego(self.ctx)
        stat = EG_solidBoolean(self.ptr, tool.ptr, oper, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj

    def intersection(self, pyego tool):
        '''
        Intersects the source BODY Object (that has the type
        SOLIDBODY, SHEETBODY or FACEBODY) with a surface or
        surfaces. The tool object contains the intersecting geometry
        in the form of a FACEBODY, SHEETBODY, SOLIDBODY or a single
        FACE.

        Parameters
        ----------

        tool:
        the FACE/FACEBODY/SHEETBODY/SOLIDBODY tool object

        Returns
        -------
        new:
        The resultant MODEL object which contains the set of WIREBODY
        BODY objects (this is because there may be multiple LOOPS as a
        result of the operation).  The EDGE objects contained within
        the LOOPS have the attributes of the FACE in src responsible
        for that EDGE.

        pairs:
        List of the FACE/EDGE object pairs - 2*nEdge in len.
        '''
        cdef int stat
        cdef int nobj
        cdef ego *faceEdgePairs = NULL
        new_obj = pyego(self.ctx)
        stat = EG_intersection(self.ptr, tool.ptr, &nobj, &faceEdgePairs,
                               &new_obj.ptr)
        if stat == EGADS_CONSTERR or nobj == 0:
            return None, []
        elif stat:
            _checkErr(stat)
        pairs = []
        for i in range(nobj):
            EG_referenceObject(faceEdgePairs[2*i], self.ctx.context)
            f = pyego(self.ctx)
            f.ptr = faceEdgePairs[2*i]
            pairs.append(f)

            EG_referenceObject(faceEdgePairs[2*i+1], self.ctx.context)
            e = pyego(self.ctx)
            e.ptr = faceEdgePairs[2*i+1]
            pairs.append(e)
        if faceEdgePairs:
            free(faceEdgePairs)
        return new_obj, pairs

    def imprintBody(self, list pairs):
        '''
        Imprints EDGE/LOOPs on the source BODY Object (that has the type
        SOLIDBODY, SHEETBODY or FACEBODY). The EDGE/LOOPs are
        paired with the FACEs in the source that will be scribed with the
        EDGE/LOOP.

        Parameter
        ---------
        pairs:
        list of FACE/EDGE and/or FACE/LOOP object pairs to scribe
        2*nObj in len -- can be the output from intersect
        result

        Returns
        -------
        the resultant BODY object (with the same type as the input source
        object, though the splitting of FACEBODY objects results in a
        SHEETBODY)
        '''
        cdef int stat
        cdef nobj = 0
        cdef ego* objs = NULL
        nobj = len(pairs)
        objs = <ego*>malloc(nobj*sizeof(ego))
        for i in range(nobj):
            objs[i] = (<pyego>pairs[i]).ptr
        new_obj = pyego(self.ctx)
        stat = EG_imprintBody(self.ptr, nobj/2, objs, &new_obj.ptr)
        free(objs)
        if stat:
            _checkErr(stat)
        return new_obj

    def filletBody(self, list edges, double radius):
        '''
        Fillets the EDGEs on the source BODY Object (that has the type
        SOLIDBODY or SHEETBODY).

        Parameters
        ----------
        edges:  list of EDGE objects to fillet
        radius: the radius of the fillets created

        Returns
        -------
        the resultant BODY object (with the same type as the input
        source object)
        '''
        cdef int stat
        cdef nedges = 0
        cdef ego* edgs = NULL
        cdef int *facemap
        nedges = len(edges)
        edgs = <ego*>malloc(nedges*sizeof(ego))
        for i in range(nedges):
            edgs[i] = (<pyego>edges[i]).ptr
        new_obj = pyego(self.ctx)
        stat = EG_filletBody(self.ptr, nedges, edgs, radius, &new_obj.ptr,
                             &facemap)
        free(edgs)
        free(facemap)
        if stat:
            _checkErr(stat)
        return new_obj

    def extrude(self, double dist, _dir):
        '''
        Extrudes the source Object through the distance specified.
        If the Object is either a LOOP or WIREBODY the result is a
        SHEETBODY. If the source is either a FACE or FACEBODY then
        the returned Object is a SOLIDBODY.

        Parameters
        ----------
        dist: the distance to extrude

        dir the vector that is the extrude direction (3 in length)

        Returns
        -------
        the resultant BODY object (type is one greater than the
        input source object)
        '''
        cdef int stat
        cdef double direction[3]
        for i in range(3):
            direction[i] = _dir[i]
        new_obj = pyego(self.ctx)
        stat = EG_extrude(self.ptr, dist, direction, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj

    def rotate(self, double angle, _axis):
        '''
        Rotates the source Object about the axis through the angle
        specified. If the Object is either a LOOP or WIREBODY the
        result is a SHEETBODY. If the source is either a FACE or
        FACEBODY then the returned Object is a SOLIDBODY.

        Parameters
        ----------
        angle: the angle to rotate the object through [0-360 Degrees]
        axis: a point (on the axis) and a direction (6 in length)

        Returns
        -------
        the resultant BODY object (type is one greater than the input
        source object)
        '''
        cdef int stat
        cdef double axis[6]
        for i in range(6):
            axis[i] = _axis[i]
        new_obj = pyego(self.ctx)
        stat = EG_rotate(self.ptr, angle, axis, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj

    def sweep(self, pyego spline, int mode):
        '''
        Sweeps the source Object through the “spine” specified. The
        spine can be either an EDGE, LOOP or WIREBODY. If the source
        Object is either a LOOP or WIREBODY the result is a
        SHEETBODY. If the source is either a FACE or FACEBODY then the
        returned Object is a SOLIDBODY.

        Parameters
        ----------
        spline: the Object used as guide curve segment(s) to sweep the
        source through

        mode: Integer indicating the mode
        '''
        cdef int stat
        new_obj = pyego(self.ctx)
        stat = EG_sweep(self.ptr, spline.ptr, mode, &new_obj.ptr)
        if stat:
            _checkErr(stat)
        return new_obj
