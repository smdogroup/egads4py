/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Geometry Functions
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <setjmp.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"

//#define MAPBSPLINE   /* also in egadsTopo.cpp */


#define PARAMACC 1.0e-4         // parameter accuracy
#define KNACC    1.0e-12	// knot accuracy


  extern "C" int  EG_destroyGeometry( egObject *geom );
  extern "C" int  EG_copyGeometry( const egObject *geom,
                                   /*@null@*/ double *xform, egObject **copy );
  extern "C" int  EG_flipGeometry( const egObject *geom, egObject **copy );
  extern "C" int  EG_spline1d( egObject *context, int endc, int imax,
                               const double *xyz, double tol, egObject **ecrv );
  extern "C" int  EG_spline2d( egObject *context, int endc,
                               /*@null@*/ const double **dr, int imax, int jmax,
                               const double *xyz, double tol, egObject **esrf );
  extern "C" int  EG_spline2dFit( egObject *context, const double *crosT,
                                  int imax, const double *uknot,
                                  const double *souT, const double *norT,
                                  int jmax, const double *vknot,
                                  const double *wesT, const double *easT,
                                  const double *xyz, double tol,
                                  egObject **esurf );
  extern     int  EG_copyAttrTopo( egadsBody *pbody, /*@null@*/ double *xform,
                                   gp_Trsf form, const egObject *src,
                                   egObject *dst, egObject *topObj );
  extern     void EG_checkStatus( const Handle_BRepCheck_Result tResult );

  extern "C" int  EG_getGeometry( const egObject *geom, int *oclass, int *type,
                                  egObject **rGeom, int **ivec, double **rvec );
  extern "C" int  EG_makeGeometry( egObject *context, int oclass, int mtype,
                                   /*@null@*/ egObject *refGeo, const int *ivec,
                                   const double *rvec, egObject **geom );
  extern "C" int  EG_getRange( const egObject *geom, double *range, int *pflg );
  extern "C" int  EG_curvature( const egObject *geom, const double *param,
                                double *result );
  extern "C" int  EG_evaluate( const egObject *geom, const double *param,
                               double *result );
  extern "C" int  EG_invEvaluate( const egObject *geom, double *xyz,
                                  double *param, double *result );
  extern "C" int  EG_invEvaluateGuess( const egObject *geom, double *xyz,
                                       double *param, double *result );
  extern "C" int  EG_arcLength( const egObject *geom, double t1, double t2,
                                double *alen );
  extern "C" int  EG_approximate( egObject *context, int maxdeg, double tol,
                                  const int *sizes, const double *xyzs,
                                  egObject **bspline );
  extern "C" int  EG_otherCurve( const egObject *surface, const egObject *curve,
                                 double tol, egObject **newcurve );
  extern "C" int  EG_isoCline( const egObject *surface, int UV, double value,
                               egObject **newcurve );
  extern "C" int  EG_convertToBSplineRange( egObject *geom, const double *range,
                                            egObject **bspline );
  extern "C" int  EG_convertToBSpline( egObject *geom, egObject **bspline );
  extern "C" int  EG_mapKnots( egObject *src, egObject *dst, egObject **rslt );
  extern "C" int  EG_mapSequen( egObject *src, egObject *dst, egObject **rslt );
  extern "C" void EG_mapTessTs( egTess1D src, egTess1D dst );
  extern "C" int  EG_relPosTs( egObject *geom, int n, const double *rel,
                               double *ts, double *xyzs );


// for trapping SegFaults
static jmp_buf jmpenv;

static void
segfault_handler(int x)
{
  longjmp(jmpenv, x);
}


int 
EG_destroyGeometry(egObject *geom)
{
  egObject *obj = NULL;
  
  if (geom->oclass == PCURVE) {
  
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    if (ppcurv != NULL) obj = ppcurv->basis;
    if (obj    != NULL)
      if (ppcurv->topFlg == 0) {
        EG_dereferenceObject(obj, geom);
      } else {
        EG_dereferenceTopObj(obj, geom);
      }
    if (ppcurv != NULL) delete ppcurv;  
  
  } else if (geom->oclass == CURVE) {

    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    if (pcurve != NULL) obj = pcurve->basis;
    if (obj    != NULL)
      if (pcurve->topFlg == 0) {
        EG_dereferenceObject(obj, geom);
      } else {
        EG_dereferenceTopObj(obj, geom);
      }
    if (pcurve != NULL) delete pcurve;

  } else {
  
    egadsSurface *psurf = (egadsSurface *) geom->blind;
    if (psurf != NULL) obj = psurf->basis;
    if (obj   != NULL)
      if (psurf->topFlg == 0) {
        EG_dereferenceObject(obj, geom);
      } else {
        EG_dereferenceTopObj(obj, geom);
      }
    if (psurf != NULL) delete psurf;

  }
  return EGADS_SUCCESS;
}


static int
EG_getPCurveType(Handle(Geom2d_Curve) &hCurve)
{

  Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) return LINE;

  Handle(Geom2d_Circle) hCirc = Handle(Geom2d_Circle)::DownCast(hCurve);
  if (!hCirc.IsNull()) return CIRCLE;

  Handle(Geom2d_Ellipse) hEllip = Handle(Geom2d_Ellipse)::DownCast(hCurve);
  if (!hEllip.IsNull()) return ELLIPSE;

  Handle(Geom2d_Parabola) hParab = Handle(Geom2d_Parabola)::DownCast(hCurve);
  if (!hParab.IsNull()) return PARABOLA;
 
  Handle(Geom2d_Hyperbola) hHypr = Handle(Geom2d_Hyperbola)::DownCast(hCurve);
  if (!hHypr.IsNull()) return HYPERBOLA;

  Handle(Geom2d_BezierCurve) hBezier = Handle(Geom2d_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) return BEZIER;

  Handle(Geom2d_BSplineCurve) hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) return BSPLINE;
  
  Handle(Geom2d_TrimmedCurve) hTrim = Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) return TRIMMED;

  Handle(Geom2d_OffsetCurve) hOffst = Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) return OFFSET;

  return 0;
}


void 
EG_completePCurve(egObject *geom, Handle(Geom2d_Curve) &hCurve)
{
  int                  stat;
  egObject             *obj;
  Handle(Geom2d_Curve) base;

  geom->oclass        = PCURVE;
  egadsPCurve *ppcurv = new egadsPCurve;
  ppcurv->handle      = hCurve;
  ppcurv->basis       = NULL;
  ppcurv->topFlg      = 0;
  geom->blind         = ppcurv;

  // stand alone geometry
  Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) {
    geom->mtype = LINE;
    return;
  }
  Handle(Geom2d_Circle) hCirc = Handle(Geom2d_Circle)::DownCast(hCurve);
  if (!hCirc.IsNull()) {
    geom->mtype = CIRCLE;
    return;
  }
  Handle(Geom2d_Ellipse) hEllip = Handle(Geom2d_Ellipse)::DownCast(hCurve);
  if (!hEllip.IsNull()) {
    geom->mtype = ELLIPSE;
    return;
  }
  Handle(Geom2d_Parabola) hParab = Handle(Geom2d_Parabola)::DownCast(hCurve);
  if (!hParab.IsNull()) {
    geom->mtype = PARABOLA;
    return;
  }
  Handle(Geom2d_Hyperbola) hHypr = Handle(Geom2d_Hyperbola)::DownCast(hCurve);
  if (!hHypr.IsNull()) {
    geom->mtype = HYPERBOLA;
    return;
  }
  Handle(Geom2d_BezierCurve) hBezier = Handle(Geom2d_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) {
    geom->mtype = BEZIER;
    return;
  }
  Handle(Geom2d_BSplineCurve) hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) {
    geom->mtype = BSPLINE;
    return;
  }
  
  // referencing geometry
  Handle(Geom2d_TrimmedCurve) hTrim = Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) {
    geom->mtype = TRIMMED;
    base = hTrim->BasisCurve();
  }
  Handle(Geom2d_OffsetCurve) hOffst = Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) {
    geom->mtype = OFFSET;
    base = hOffst->BasisCurve();
  }  
  if (geom->mtype == 0) {
    printf(" EGADS Error: Unknown PCurve Type!\n");
    return;
  }
  
  // make the reference curve
  stat = EG_makeObject(EG_context(geom), &obj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: make Object = %d (EG_completePCurve)!\n", stat);
    return;
  }
  ppcurv->basis = obj;
  if (geom->topObj == EG_context(geom)) {
    obj->topObj = geom;
  } else {
    obj->topObj = geom->topObj;
  }
  EG_completePCurve(obj,  base);
  EG_referenceObject(obj, geom);
}


void 
EG_completeCurve(egObject *geom, Handle(Geom_Curve) &hCurve)
{
  int                stat;
  egObject           *obj;
  Handle(Geom_Curve) base;

  geom->oclass       = CURVE;
  egadsCurve *pcurve = new egadsCurve;
  pcurve->handle     = hCurve;
  pcurve->basis      = NULL;
  pcurve->topFlg     = 0;
  geom->blind        = pcurve;

  // stand alone geometry
  Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) {
    geom->mtype = LINE;
    return;
  }
  Handle(Geom_Circle) hCirc = Handle(Geom_Circle)::DownCast(hCurve);
  if (!hCirc.IsNull()) {
    geom->mtype = CIRCLE;
    return;
  }
  Handle(Geom_Ellipse) hEllip = Handle(Geom_Ellipse)::DownCast(hCurve);
  if (!hEllip.IsNull()) {
    geom->mtype = ELLIPSE;
    return;
  }
  Handle(Geom_Parabola) hParab = Handle(Geom_Parabola)::DownCast(hCurve);
  if (!hParab.IsNull()) {
    geom->mtype = PARABOLA;
    return;
  }
  Handle(Geom_Hyperbola) hHypr = Handle(Geom_Hyperbola)::DownCast(hCurve);
  if (!hHypr.IsNull()) {
    geom->mtype = HYPERBOLA;
    return;
  }
  Handle(Geom_BezierCurve) hBezier = Handle(Geom_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) {
    geom->mtype = BEZIER;
    return;
  }
  Handle(Geom_BSplineCurve) hBSpline = Handle(Geom_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) {
    geom->mtype = BSPLINE;
    return;
  }
  
  // referencing geometry
  Handle(Geom_TrimmedCurve) hTrim = Handle(Geom_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) {
    geom->mtype = TRIMMED;
    base = hTrim->BasisCurve();
  }
  Handle(Geom_OffsetCurve) hOffst = Handle(Geom_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) {
    geom->mtype = OFFSET;
    base = hOffst->BasisCurve();
  }  
  if (geom->mtype == 0) {
    printf(" EGADS Error: Unknown Curve Type!\n");
    return;
  }
  
  // make the reference curve
  stat = EG_makeObject(EG_context(geom), &obj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: make Object = %d (EG_completeCurve)!\n", stat);
    return;
  }
  pcurve->basis = obj;
  if (geom->topObj == EG_context(geom)) {
    obj->topObj = geom;
  } else {
    obj->topObj = geom->topObj;
  }
  EG_completeCurve(obj, base);
  EG_referenceObject(obj, geom);
}


void EG_completeSurf(egObject *geom, Handle(Geom_Surface) &hSurf)
{
  int                  stat;
  egObject             *obj;
  Handle(Geom_Surface) base;
  Handle(Geom_Curve)   curve;

  geom->oclass        = SURFACE;
  egadsSurface *psurf = new egadsSurface;
  psurf->handle       = hSurf;
  psurf->basis        = NULL;
  psurf->topFlg       = 0;
  geom->blind         = psurf;
  
  // stand alone geometry
  Handle(Geom_Plane) hPlane = Handle(Geom_Plane)::DownCast(hSurf);
  if (!hPlane.IsNull()) {
    geom->mtype = PLANE;
    return;
  }
  Handle(Geom_SphericalSurface) hSphere = 
      Handle(Geom_SphericalSurface)::DownCast(hSurf);
  if (!hSphere.IsNull()) {
    geom->mtype = SPHERICAL;
    return;
  }
  Handle(Geom_ConicalSurface) hCone = 
    Handle(Geom_ConicalSurface)::DownCast(hSurf);
  if (!hCone.IsNull()) {
    geom->mtype = CONICAL;
    return;
  }
  Handle(Geom_CylindricalSurface) hCyl =
    Handle(Geom_CylindricalSurface)::DownCast(hSurf);
  if (!hCyl.IsNull()) {
    geom->mtype = CYLINDRICAL;
    return;
  }
  Handle(Geom_ToroidalSurface) hTorus =
    Handle(Geom_ToroidalSurface)::DownCast(hSurf);
  if (!hTorus.IsNull()) {
    geom->mtype = TOROIDAL;
    return;
  }
  Handle(Geom_BezierSurface) hBezier = 
    Handle(Geom_BezierSurface)::DownCast(hSurf);
  if (!hBezier.IsNull()) {
    geom->mtype = BEZIER;
    return;
  }
  Handle(Geom_BSplineSurface) hBSpline = 
    Handle(Geom_BSplineSurface)::DownCast(hSurf);
  if (!hBSpline.IsNull()) {
    geom->mtype = BSPLINE;
    return;
  }

  // referencing geometry -- surface
  Handle(Geom_OffsetSurface) hOffst = 
    Handle(Geom_OffsetSurface)::DownCast(hSurf);
  if (!hOffst.IsNull()) {
    geom->mtype = OFFSET;
    base = hOffst->BasisSurface();
    stat = EG_makeObject(EG_context(geom), &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Object = %d (EG_completeSurface)!\n", stat);
      return;
    }
    psurf->basis = obj;
    if (geom->topObj == EG_context(geom)) {
      obj->topObj = geom;
    } else {
      obj->topObj = geom->topObj;
    }
    EG_completeSurf(obj, base);
    EG_referenceObject(obj, geom);
    return;
  } 
  Handle(Geom_RectangularTrimmedSurface) hTrim = 
    Handle(Geom_RectangularTrimmedSurface)::DownCast(hSurf);
  if (!hTrim.IsNull()) {
    geom->mtype = TRIMMED;
    base = hTrim->BasisSurface();
    stat = EG_makeObject(EG_context(geom), &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Object = %d (EG_completeSurface)!\n", stat);
      return;
    }
    psurf->basis = obj;
    if (geom->topObj == EG_context(geom)) {
      obj->topObj = geom;
    } else {
      obj->topObj = geom->topObj;
    }
    EG_completeSurf(obj, base);
    EG_referenceObject(obj, geom);
    return;
  } 
  
  // referencing geometry -- curve
  Handle(Geom_SurfaceOfLinearExtrusion) hSLExtr = 
    Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(hSurf);
  if (!hSLExtr.IsNull()) {
    geom->mtype = EXTRUSION;
    curve       = hSLExtr->BasisCurve();
  }
  Handle(Geom_SurfaceOfRevolution) hSORev = 
    Handle(Geom_SurfaceOfRevolution)::DownCast(hSurf);
  if (!hSORev.IsNull()) {
    geom->mtype = REVOLUTION;
    curve       = hSORev->BasisCurve();
  }
  if (geom->mtype == 0) {
    printf(" EGADS Error: Unknown Surface Type!\n");
    return;
  }
  
  // make the reference curve
  stat = EG_makeObject(EG_context(geom), &obj);
  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Error: make Curve = %d (EG_completeSurface)!\n", stat);
    return;
  }
  psurf->basis = obj;
  if (geom->topObj == EG_context(geom)) {
    obj->topObj = geom;
  } else {
    obj->topObj = geom->topObj;
  }
  EG_completeCurve(obj, curve);
  EG_referenceObject(obj, geom);  
}


int
EG_copyGeometry(const egObject *geom, /*@null@*/ double *xform, egObject **copy)
{
  int      stat, outLevel;
  egObject *obj, *context;
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != CURVE) && (geom->oclass != SURFACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  context  = EG_context(geom);
          
  gp_Trsf form = gp_Trsf();
  if (xform != NULL)
    form.SetValues(xform[ 0], xform[ 1], xform[ 2], xform[ 3],
                   xform[ 4], xform[ 5], xform[ 6], xform[ 7],
                   xform[ 8], xform[ 9], xform[10], xform[11]
#if CASVER < 680
                   , Precision::Confusion(), Precision::Angular()
#endif
                   );
                   
  if (geom->oclass == CURVE) {

    egadsCurve           *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve)    hCurve = pcurve->handle;
    Handle(Geom_Geometry) nGeom  = hCurve->Transformed(form);
    Handle(Geom_Curve)    nCurve = Handle(Geom_Curve)::DownCast(nGeom);
    if (nCurve.IsNull()) {
      if (outLevel > 0) 
        printf(" EGADS Error: XForm Curve Failed (EG_copyGeometry)!\n");
      return EGADS_CONSTERR;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) 
        printf(" EGADS Error: makeObject = %d (EG_copyGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeCurve(obj, nCurve);
    
  } else {
  
    egadsSurface         *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface)  hSurf = psurf->handle;
    Handle(Geom_Geometry) nGeom = hSurf->Transformed(form);
    Handle(Geom_Surface)  nSurf = Handle(Geom_Surface)::DownCast(nGeom);
    if (nSurf.IsNull()) {
      if (outLevel > 0) 
        printf(" EGADS Error: XForm Surface Failed (EG_copyGeometry)!\n");
      return EGADS_CONSTERR;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) 
        printf(" EGADS Error: makeObject = %d (EG_copyGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeSurf(obj, nSurf);
    
  }

  EG_referenceObject(obj, context);
  *copy = obj;
  return EGADS_SUCCESS;
}


int
EG_flipGeometry(const egObject *geom, egObject **copy)
{
  int      stat, outLevel;
  egObject *obj, *context;
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE)     && (geom->oclass != CURVE) &&
      (geom->oclass != SURFACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  context  = EG_context(geom);
    
  if (geom->oclass == PCURVE) {
  
    egadsPCurve          *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve)  hPCurv = ppcurv->handle;
    Handle(Geom2d_Curve)  nPCurv = hPCurv->Reversed();

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) 
        printf(" EGADS Error: makeObject = %d (EG_flipGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completePCurve(obj, nPCurv);    
    
  } else if (geom->oclass == CURVE) {

    egadsCurve         *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve)  hCurve = pcurve->handle;
    Handle(Geom_Curve)  nCurve = hCurve->Reversed();

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) 
        printf(" EGADS Error: makeObject = %d (EG_flipGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeCurve(obj, nCurve);
    
  } else {
  
    egadsSurface         *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface)  hSurf = psurf->handle;
    Handle(Geom_Surface)  nSurf = hSurf->UReversed();

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0) 
        printf(" EGADS Error: makeObject = %d (EG_flipGeometry)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeSurf(obj, nSurf);
    
  }

  EG_referenceObject(obj, context);
  *copy = obj;
  return EGADS_SUCCESS;
}


int
EG_getGeometry(const egObject     *geom, int *oclass, int *type, 
                     egObject **refGeom, int **ivec,  double **rvec)
{
  int    *ints = NULL, i, j, len, outLevel;
  double *data = NULL;

  *ivec    = NULL;
  *rvec    = NULL;
  *refGeom = NULL;
  *oclass  = *type = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass < PCURVE) || (geom->oclass > SURFACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  *oclass = geom->oclass;
  *type   = geom->mtype;

  if (geom->oclass == PCURVE) {

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    *refGeom = ppcurv->basis;

    switch (geom->mtype) {

      case LINE:
        data = (double *) EG_alloc(4*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PLine (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
          gp_Dir2d direct = hLine->Direction();
          gp_Pnt2d locat  = hLine->Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = direct.X();
          data[3] = direct.Y();
          *rvec   = data;
        }
        break;
      
      case CIRCLE:
        data = (double *) EG_alloc(7*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PCircle (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Circle) hCirc = Handle(Geom2d_Circle)::DownCast(hCurve);
          gp_Circ2d circ  = hCirc->Circ2d();
          gp_Ax2d   xaxis = circ.XAxis();
          gp_Ax2d   yaxis = circ.YAxis();
          gp_Pnt2d  locat = circ.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = circ.Radius();
          *rvec   = data;
        }
        break;
      
      case ELLIPSE:
        data = (double *) EG_alloc(8*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PEllipse (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Ellipse) hEllip = Handle(Geom2d_Ellipse)::DownCast(hCurve);
          gp_Elips2d elips = hEllip->Elips2d();
          gp_Ax2d    xaxis = elips.XAxis();
          gp_Ax2d    yaxis = elips.YAxis();
          gp_Pnt2d   locat = elips.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = elips.MajorRadius();
          data[7] = elips.MinorRadius();
          *rvec   = data;
        }
        break;
      
      case PARABOLA:
        data = (double *) EG_alloc(7*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PParabola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Parabola) hParab=Handle(Geom2d_Parabola)::DownCast(hCurve);
          gp_Parab2d parab = hParab->Parab2d();
          gp_Ax22d   axes  = parab.Axis();
          gp_Ax2d    xaxis = axes.XAxis();
          gp_Ax2d    yaxis = axes.YAxis();
          gp_Pnt2d   locat = parab.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = parab.Focal();
          *rvec   = data;
        }
        break;
      
      case HYPERBOLA:
        data = (double *) EG_alloc(8*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PHyperbola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_Hyperbola) hHypr = 
            Handle(Geom2d_Hyperbola)::DownCast(hCurve);
          gp_Hypr2d hypr  = hHypr->Hypr2d();
          gp_Ax2d   xaxis = hypr.XAxis();
          gp_Ax2d   yaxis = hypr.YAxis();
          gp_Pnt2d  locat = hypr.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = xaxis.Direction().X();
          data[3] = xaxis.Direction().Y();
          data[4] = yaxis.Direction().X();
          data[5] = yaxis.Direction().Y();
          data[6] = hypr.MajorRadius();
          data[7] = hypr.MinorRadius();
          *rvec   = data;
        }
        break;
      
      case TRIMMED:
        data = (double *) EG_alloc(2*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PTrimmed Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_TrimmedCurve) hTrim = 
            Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
          data[0] = hTrim->FirstParameter();
          data[1] = hTrim->LastParameter();
          *rvec   = data;
        }
        break;
      
      case BEZIER:
        ints = (int *) EG_alloc(3*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PBezier Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_BezierCurve) hBezier = 
            Handle(Geom2d_BezierCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBezier->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBezier->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBezier->Degree();
          ints[2] = hBezier->NbPoles();
          len = ints[2]*2;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc on PBezier Data (EG_getGeometry)!\n");
              EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= ints[2]; i++, len+=2) {
            gp_Pnt2d P  = hBezier->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++) data[len] = hBezier->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;
      
      case BSPLINE:
        ints = (int *) EG_alloc(4*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on PBSpline Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_BSplineCurve) hBSpline = 
            Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBSpline->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBSpline->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBSpline->Degree();
          ints[2] = hBSpline->NbPoles();
          ints[3] = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++)
            ints[3] += hBSpline->Multiplicity(i);
          len = ints[3] + ints[2]*2;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc PBSpline Data (EG_getGeometry)!\n");
              EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++) {
            int km = hBSpline->Multiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->Knot(i);
          }
          for (i = 1; i <= ints[2]; i++, len+=2) {
            gp_Pnt2d P  = hBSpline->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++) 
              data[len] = hBSpline->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;
      
      case OFFSET:
        data = (double *) EG_alloc(1*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on POffset Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom2d_OffsetCurve) hOffst = 
            Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
          data[0] = hOffst->Offset();
          *rvec   = data;
        }
        break;
    }  
  
  } else if (geom->oclass == CURVE) {
    
    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    *refGeom = pcurve->basis;

    switch (geom->mtype) {

      case LINE:
        data = (double *) EG_alloc(6*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Line (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);
          gp_Lin line   = hLine->Lin();
          gp_Dir direct = line.Direction();
          gp_Pnt locat  = line.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = direct.X();
          data[4] = direct.Y();
          data[5] = direct.Z();
          *rvec   = data;
        }
        break;
      
      case CIRCLE:
        data = (double *) EG_alloc(10*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Circle (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Circle) hCirc = Handle(Geom_Circle)::DownCast(hCurve);
          gp_Circ circ  = hCirc->Circ();
          gp_Ax1  xaxis = circ.XAxis();
          gp_Ax1  yaxis = circ.YAxis();
          gp_Pnt  locat = circ.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          data[9] = circ.Radius();
          *rvec   = data;
        }
        break;
      
      case ELLIPSE:
        data = (double *) EG_alloc(11*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Ellipse (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Ellipse) hEllip = Handle(Geom_Ellipse)::DownCast(hCurve);
          gp_Elips elips = hEllip->Elips();
          gp_Ax1   xaxis = elips.XAxis();
          gp_Ax1   yaxis = elips.YAxis();
          gp_Pnt   locat = elips.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = xaxis.Direction().X();
          data[ 4] = xaxis.Direction().Y();
          data[ 5] = xaxis.Direction().Z();
          data[ 6] = yaxis.Direction().X();
          data[ 7] = yaxis.Direction().Y();
          data[ 8] = yaxis.Direction().Z();
          data[ 9] = elips.MajorRadius();
          data[10] = elips.MinorRadius();
          *rvec    = data;
        }
        break;
      
      case PARABOLA:
        data = (double *) EG_alloc(10*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Parabola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Parabola) hParab=Handle(Geom_Parabola)::DownCast(hCurve);
          gp_Parab parab = hParab->Parab();
          gp_Ax1   xaxis = parab.XAxis();
          gp_Ax1   yaxis = parab.YAxis();
          gp_Pnt   locat = parab.Location();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          data[9] = parab.Focal();
          *rvec   = data;
        }
        break;
      
      case HYPERBOLA:
        data = (double *) EG_alloc(11*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Hyperbola (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Hyperbola) hHypr = 
            Handle(Geom_Hyperbola)::DownCast(hCurve);
          gp_Hypr hypr  = hHypr->Hypr();
          gp_Ax1  xaxis = hypr.XAxis();
          gp_Ax1  yaxis = hypr.YAxis();
          gp_Pnt  locat = hypr.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = xaxis.Direction().X();
          data[ 4] = xaxis.Direction().Y();
          data[ 5] = xaxis.Direction().Z();
          data[ 6] = yaxis.Direction().X();
          data[ 7] = yaxis.Direction().Y();
          data[ 8] = yaxis.Direction().Z();
          data[ 9] = hypr.MajorRadius();
          data[10] = hypr.MinorRadius();
          *rvec    = data;
        }
        break;
      
      case TRIMMED:
        data = (double *) EG_alloc(2*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Trimmed Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_TrimmedCurve) hTrim = 
            Handle(Geom_TrimmedCurve)::DownCast(hCurve);
          data[0] = hTrim->FirstParameter();
          data[1] = hTrim->LastParameter();
          *rvec   = data;
        }
        break;
      
      case BEZIER:
        ints = (int *) EG_alloc(3*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Bezier Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BezierCurve) hBezier = 
            Handle(Geom_BezierCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBezier->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBezier->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBezier->Degree();
          ints[2] = hBezier->NbPoles();
          len = ints[2]*3;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc on Bezier Data (EG_getGeometry)!\n");
              EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= ints[2]; i++, len+=3) {
            gp_Pnt P    = hBezier->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
            data[len+2] = P.Z();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++) data[len] = hBezier->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;
      
      case BSPLINE:
        ints = (int *) EG_alloc(4*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on BSpline Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BSplineCurve) hBSpline = 
            Handle(Geom_BSplineCurve)::DownCast(hCurve);
          int rational = 0;
          if (hBSpline->IsRational()) rational = 1;
          ints[0] = rational*2;
          if (hBSpline->IsPeriodic()) ints[0] |= 4;
          ints[1] = hBSpline->Degree();
          ints[2] = hBSpline->NbPoles();
          ints[3] = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++)
            ints[3] += hBSpline->Multiplicity(i);
          len = ints[3] + ints[2]*3;
          if (rational == 1) len += ints[2];
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc BSpline Data (EG_getGeometry)!\n");
              EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= hBSpline->NbKnots(); i++) {
            int km = hBSpline->Multiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->Knot(i);
          }
          for (i = 1; i <= ints[2]; i++, len+=3) {
            gp_Pnt P    = hBSpline->Pole(i);
            data[len  ] = P.X();
            data[len+1] = P.Y();
            data[len+2] = P.Z();
          }
          if (rational == 1)
            for (i = 1; i <= ints[2]; i++,len++) 
              data[len] = hBSpline->Weight(i);
          *ivec = ints;
          *rvec = data;
        }
        break;
      
      case OFFSET:
        data = (double *) EG_alloc(4*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Offset Curve (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_OffsetCurve) hOffst = 
            Handle(Geom_OffsetCurve)::DownCast(hCurve);
          gp_Dir direct = hOffst->Direction();
          data[0] = direct.X();
          data[1] = direct.Y();
          data[2] = direct.Z();
          data[3] = hOffst->Offset();
          *rvec   = data;
        }
        break;
    }
    
  } else {
  
    egadsSurface *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    *refGeom = psurf->basis;

    switch (geom->mtype) {
    
      case PLANE:
        data = (double *) EG_alloc(9*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Plane (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_Plane) hPlane = Handle(Geom_Plane)::DownCast(hSurf);
          gp_Pln plane = hPlane->Pln();
          gp_Pnt locat = plane.Location();
          gp_Ax1 xaxis = plane.XAxis();
          gp_Ax1 yaxis = plane.YAxis();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          *rvec   = data;
        }
        break;
        
      case SPHERICAL:
        data = (double *) EG_alloc(10*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Sphere (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_SphericalSurface) hSphere = 
            Handle(Geom_SphericalSurface)::DownCast(hSurf);
          gp_Sphere sphere = hSphere->Sphere();
          gp_Pnt locat     = sphere.Location();
          gp_Ax1 xaxis     = sphere.XAxis();
          gp_Ax1 yaxis     = sphere.YAxis();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = xaxis.Direction().X();
          data[4] = xaxis.Direction().Y();
          data[5] = xaxis.Direction().Z();
          data[6] = yaxis.Direction().X();
          data[7] = yaxis.Direction().Y();
          data[8] = yaxis.Direction().Z();
          data[9] = sphere.Radius();
          *rvec   = data;
        }
        break;

      case CONICAL:
        data = (double *) EG_alloc(14*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Conical (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_ConicalSurface) hCone = 
            Handle(Geom_ConicalSurface)::DownCast(hSurf);
          gp_Cone cone  = hCone->Cone();
          gp_Ax3  axes  = cone.Position();
          gp_Pnt  locat = cone.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = axes.XDirection().X();
          data[ 4] = axes.XDirection().Y();
          data[ 5] = axes.XDirection().Z();
          data[ 6] = axes.YDirection().X();
          data[ 7] = axes.YDirection().Y();
          data[ 8] = axes.YDirection().Z();
          data[ 9] = axes.Direction().X();
          data[10] = axes.Direction().Y();
          data[11] = axes.Direction().Z();
          data[12] = cone.SemiAngle();
          data[13] = cone.RefRadius();
          *rvec    = data;
        }
        break;
        
      case CYLINDRICAL:
        data = (double *) EG_alloc(13*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Cylinder (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_CylindricalSurface) hCyl =
            Handle(Geom_CylindricalSurface)::DownCast(hSurf);
          gp_Cylinder cyl   = hCyl->Cylinder();
          gp_Ax3      axes  = cyl.Position();
          gp_Pnt      locat = cyl.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = axes.XDirection().X();
          data[ 4] = axes.XDirection().Y();
          data[ 5] = axes.XDirection().Z();
          data[ 6] = axes.YDirection().X();
          data[ 7] = axes.YDirection().Y();
          data[ 8] = axes.YDirection().Z();
          data[ 9] = axes.Direction().X();
          data[10] = axes.Direction().Y();
          data[11] = axes.Direction().Z();
          data[12] = cyl.Radius();
          *rvec    = data;
        }
        break;

      case TOROIDAL:
        data = (double *) EG_alloc(14*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Cylinder (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_ToroidalSurface) hTorus =
            Handle(Geom_ToroidalSurface)::DownCast(hSurf);
          gp_Torus torus = hTorus->Torus();
          gp_Ax3   axes  = torus.Position();
          gp_Pnt   locat = torus.Location();
          data[ 0] = locat.X();
          data[ 1] = locat.Y();
          data[ 2] = locat.Z();
          data[ 3] = axes.XDirection().X();
          data[ 4] = axes.XDirection().Y();
          data[ 5] = axes.XDirection().Z();
          data[ 6] = axes.YDirection().X();
          data[ 7] = axes.YDirection().Y();
          data[ 8] = axes.YDirection().Z();
          data[ 9] = axes.Direction().X();
          data[10] = axes.Direction().Y();
          data[11] = axes.Direction().Z();
          data[12] = torus.MajorRadius();
          data[13] = torus.MinorRadius();
          *rvec    = data;
        }
        break;
     
      case BEZIER:
        ints = (int *) EG_alloc(5*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc on Bezier header (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BezierSurface) hBezier = 
            Handle(Geom_BezierSurface)::DownCast(hSurf);
          int rational = 0;
          if (hBezier->IsURational()) rational = 1;
          if (hBezier->IsVRational()) rational = 1;
          ints[0] = rational*2;
          if (hBezier->IsUPeriodic()) ints[0] |=  4;
          if (hBezier->IsVPeriodic()) ints[0] |=  8;
          ints[1] = hBezier->UDegree();
          ints[3] = hBezier->VDegree();
          ints[2] = hBezier->NbUPoles();
          ints[4] = hBezier->NbVPoles();
          int nCP = ints[2];
          nCP    *= ints[4];
          len     = nCP*3;
          if (rational == 1) len += nCP;
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc on Bezier Surf (EG_getGeometry)!\n");
              EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (j = 1; j <= ints[4]; j++)
            for (i = 1; i <= ints[2]; i++, len+=3) {
              gp_Pnt P    = hBezier->Pole(i, j);
              data[len  ] = P.X();
              data[len+1] = P.Y();
              data[len+2] = P.Z();
            }
          if (rational == 1)
            for (j = 1; j <= ints[4]; j++)
              for (i = 1; i <= ints[2]; i++,len++) 
                data[len] = hBezier->Weight(i, j);
          *ivec = ints;
          *rvec = data;
        }
        break;

      case BSPLINE:
        ints = (int *) EG_alloc(7*sizeof(int));
        if (ints == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc BSpline header (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_BSplineSurface) hBSpline = 
            Handle(Geom_BSplineSurface)::DownCast(hSurf);
          int rational = 0;
          if (hBSpline->IsURational()) rational = 1;
          if (hBSpline->IsVRational()) rational = 1;
          ints[0] = rational*2;
          if (hBSpline->IsUPeriodic()) ints[0] |=  4;
          if (hBSpline->IsVPeriodic()) ints[0] |=  8;
          ints[1] = hBSpline->UDegree();
          ints[4] = hBSpline->VDegree();
          ints[2] = hBSpline->NbUPoles();
          ints[5] = hBSpline->NbVPoles();
          ints[3] = ints[6] = 0;
          for (i = 1; i <= hBSpline->NbUKnots(); i++)
            ints[3] += hBSpline->UMultiplicity(i);
          for (i = 1; i <= hBSpline->NbVKnots(); i++)
            ints[6] += hBSpline->VMultiplicity(i);
          int nCP = ints[2];
          nCP    *= ints[5];
          len     = ints[3] + ints[6] + nCP*3;
          if (rational == 1) len += nCP;
          data = (double *) EG_alloc(len*sizeof(double));
          if (data == NULL) {
            if (outLevel > 0)
              printf(" EGADS Error: Malloc BSpline Surf (EG_getGeometry)!\n");
              EG_free(ints);
            return EGADS_MALLOC;
          }
          len = 0;
          for (i = 1; i <= hBSpline->NbUKnots(); i++) {
            int km = hBSpline->UMultiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->UKnot(i);
          }
          for (i = 1; i <= hBSpline->NbVKnots(); i++) {
            int km = hBSpline->VMultiplicity(i);
            for (j = 1; j <= km; j++, len++) data[len] = hBSpline->VKnot(i);
          }
          for (j = 1; j <= ints[5]; j++)
            for (i = 1; i <= ints[2]; i++, len+=3) {
              gp_Pnt P    = hBSpline->Pole(i, j);
              data[len  ] = P.X();
              data[len+1] = P.Y();
              data[len+2] = P.Z();
            }
          if (rational == 1)
            for (j = 1; j <= ints[5]; j++)
              for (i = 1; i <= ints[2]; i++,len++) 
                data[len] = hBSpline->Weight(i, j);
          *ivec = ints;
          *rvec = data;
        }
        break;
        
      case OFFSET:
        data = (double *) EG_alloc(sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Offset Surface (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_OffsetSurface) hOffst = 
            Handle(Geom_OffsetSurface)::DownCast(hSurf);
          data[0] = hOffst->Offset();
          *rvec   = data;
        }
        break;
        
      case TRIMMED:
        data = (double *) EG_alloc(4*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Trimmed Surface (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_RectangularTrimmedSurface) hTrim = 
            Handle(Geom_RectangularTrimmedSurface)::DownCast(hSurf);
          hTrim->Bounds(data[0], data[1], data[2], data[3]);
          *rvec = data;
        }
        break;
        
      case EXTRUSION:
        data = (double *) EG_alloc(3*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Linear Extrusion (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_SurfaceOfLinearExtrusion) hSLExtr = 
            Handle(Geom_SurfaceOfLinearExtrusion)::DownCast(hSurf);
          gp_Dir direct = hSLExtr->Direction();
          data[0] = direct.X();
          data[1] = direct.Y();
          data[2] = direct.Z();
          *rvec   = data;
        }
        break;   

      case REVOLUTION:
        data = (double *) EG_alloc(6*sizeof(double));
        if (data == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Malloc Revolved Surface (EG_getGeometry)!\n");
          return EGADS_MALLOC;
        } else {
          Handle(Geom_SurfaceOfRevolution) hSORev = 
            Handle(Geom_SurfaceOfRevolution)::DownCast(hSurf);
          gp_Pnt locat     = hSORev->Location();
          gp_Ax1 axis      = hSORev->Axis();
          data[0] = locat.X();
          data[1] = locat.Y();
          data[2] = locat.Z();
          data[3] = axis.Direction().X();
          data[4] = axis.Direction().Y();
          data[5] = axis.Direction().Z();
          *rvec   = data;
        }
        break;

    }
  
  }

  return EGADS_SUCCESS;
}


int  
EG_makeGeometry(egObject *context, int oclass, int mtype, 
                /*@null@*/ egObject *refGeom, /*@null@*/ const int *ints,
                const double *data, egObject **geom)
{
  int      i, j, len, stat, outLevel, rational, nmult;
  egObject *obj = NULL, *basis = NULL;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  outLevel = EG_outLevel(context);

  if ((oclass < PCURVE) || (oclass > SURFACE)) {
    if (outLevel > 0)
      printf(" EGADS Error: oclass = %d (EG_makeGeometry)!\n", oclass);
    return EGADS_NOTGEOM;
  }
  
  if (oclass == PCURVE) {
  
    if ((mtype < LINE) || (mtype > OFFSET)) {
      if (outLevel > 0)
        printf(" EGADS Error: PCurve mtype = %d (EG_makeGeometry)!\n", 
               mtype);
      return EGADS_RANGERR;
    }
    if ((mtype == TRIMMED) || (mtype == OFFSET)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: PCrv mtype = %d Ref is NULL (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != PCURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: PCrv mtype = %d Ref is %d (EG_makeGeometry)!\n", 
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: PCrv mtype = %d Ref has no data (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NODATA;
      }
    }
    
    Handle(Geom2d_Curve) hCurve;
    try {
      switch (mtype) {

        case LINE:
          {
            gp_Pnt2d pntl(data[0], data[1]);
            gp_Dir2d dirl(data[2], data[3]);
            hCurve = new Geom2d_Line(pntl, dirl);
          }
          break;
    
        case CIRCLE:
          {
            gp_Pnt2d pntc(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pntc, dirx, diry);
            hCurve = new Geom2d_Circle(axi2, data[6]);
          }
          break;

        case ELLIPSE:
          {
            gp_Pnt2d pnte(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pnte, dirx, diry);
            hCurve = new Geom2d_Ellipse(axi2, data[6], data[7]);
          }
          break;
      
        case PARABOLA:
          {
            gp_Pnt2d pntp(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pntp, dirx, diry);
            hCurve = new Geom2d_Parabola(axi2, data[6]);
          }
          break;
     
        case HYPERBOLA:
          {
            gp_Pnt2d pnth(data[0], data[1]);
            gp_Dir2d dirx(data[2], data[3]);
            gp_Dir2d diry(data[4], data[5]);
            gp_Ax22d axi2(pnth, dirx, diry);
            hCurve = new Geom2d_Hyperbola(axi2, data[6], data[7]);
          } 
          break;
   
        case TRIMMED:
          {
            egadsPCurve *ppcurv = (egadsPCurve *) basis->blind;
            hCurve = new Geom2d_TrimmedCurve(ppcurv->handle, data[0], data[1]);
          }
          break;
     
        case BEZIER:
          {
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            len = ints[2];
            TColgp_Array1OfPnt2d aPoles(1, len);
            for (i = 1; i <= ints[2]; i++)
              aPoles(i) = gp_Pnt2d(data[2*i-2], data[2*i-1]);
            if (rational == 0) {
              hCurve = new Geom2d_BezierCurve(aPoles);
            } else {
              TColStd_Array1OfReal aWeights(1, len);
              len = 2*ints[2];
              for (i = 1; i <= ints[2]; i++, len++) 
                aWeights(i) = data[len];
              hCurve = new Geom2d_BezierCurve(aPoles, aWeights);
            }
          }
          break;
    
        case BSPLINE:
          {
            Standard_Boolean periodic = Standard_False;
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            if ((ints[0]&4) != 0) periodic = Standard_True;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) len++;
            TColStd_Array1OfReal    aKnots(1, len);
            TColStd_Array1OfInteger aMults(1, len);
            aKnots(1) = data[0];
            aMults(1) = nmult = 1;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) {
                len++;
                aKnots(len) = data[i];
                aMults(len) = nmult = 1;
              } else {
                nmult++;
                aMults(len) = nmult;
              }
            len = ints[3];
            TColgp_Array1OfPnt2d aPoles(1, ints[2]);
            for (i = 1; i <= ints[2]; i++, len+=2)
              aPoles(i) = gp_Pnt2d(data[len], data[len+1]);
            if (rational == 0) {
              hCurve = new Geom2d_BSplineCurve(aPoles, aKnots, aMults,
                                               ints[1], periodic);
            } else {
              TColStd_Array1OfReal aWeights(1, ints[2]);
              for (i = 1; i <= ints[2]; i++, len++) 
                aWeights(i) = data[len];
              hCurve = new Geom2d_BSplineCurve(aPoles, aWeights, aKnots, 
                                               aMults, ints[1], periodic);
            }
          }
          break;
    
        case OFFSET:
          {
            egadsPCurve *ppcurv = (egadsPCurve *) basis->blind;
            hCurve = new Geom2d_OffsetCurve(ppcurv->handle, data[0]);
          }
          break;
      }
    }
    catch (Standard_Failure) 
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      Handle_Standard_Failure e = Standard_Failure::Caught();
      printf("                %s\n", e->GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      return EGADS_GEOMERR;
    }
    
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make PCurve = %d (EG_makeGeometry)!\n", stat);
      return stat;
    }
    obj->oclass         = PCURVE;
    obj->mtype          = mtype;
    egadsPCurve *ppcurv = new egadsPCurve;
    ppcurv->handle      = hCurve;
    ppcurv->basis       = basis;
    ppcurv->topFlg      = 1;
    obj->blind          = ppcurv;
    EG_referenceObject(obj, context);
    if (basis != NULL) EG_referenceTopObj(basis, obj);
  
  } else if (oclass == CURVE) {

    if ((mtype < LINE) || (mtype > OFFSET)) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve mtype = %d (EG_makeGeometry)!\n", 
               mtype);
      return EGADS_RANGERR;
    }
    if ((mtype == TRIMMED) || (mtype == OFFSET)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Crv mtype = %d Ref is NULL (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != CURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: Crv mtype = %d Ref is %d (EG_makeGeometry)!\n", 
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Crv mtype = %d Ref has no data (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NODATA;
      }
    }
    
    Handle(Geom_Curve) hCurve;
    try {
      switch (mtype) {

        case LINE:
          {
            gp_Pnt pntl(data[0], data[1], data[2]);
            gp_Dir dirl(data[3], data[4], data[5]);
            hCurve = new Geom_Line(pntl, dirl);
          }
          break;
    
        case CIRCLE:
          {
            gp_Pnt pntc(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pntc, dirz, dirx);
            hCurve = new Geom_Circle(axi2, data[9]);
          }
          break;

        case ELLIPSE:
          {
            gp_Pnt pnte(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pnte, dirz, dirx);
            hCurve = new Geom_Ellipse(axi2, data[9], data[10]);
          }
          break;
      
        case PARABOLA:
          {
            gp_Pnt pntp(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pntp, dirz, dirx);
            hCurve = new Geom_Parabola(axi2, data[9]);
          }
          break;
     
        case HYPERBOLA:
          {
            gp_Pnt pnth(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2(pnth, dirz, dirx);
            hCurve = new Geom_Hyperbola(axi2, data[9], data[10]);
          } 
          break;
   
        case TRIMMED:
          {
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            hCurve = new Geom_TrimmedCurve(pcurve->handle, data[0], data[1]);
          }
          break;
     
        case BEZIER:
          {
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            len = ints[2];
            TColgp_Array1OfPnt aPoles(1, len);
            for (i = 1; i <= ints[2]; i++)
              aPoles(i) = gp_Pnt(data[3*i-3], data[3*i-2], data[3*i-1]);
            if (rational == 0) {
              hCurve = new Geom_BezierCurve(aPoles);
            } else {
              TColStd_Array1OfReal aWeights(1, len);
              len = 3*ints[2];
              for (i = 1; i <= ints[2]; i++, len++) 
                aWeights(i) = data[len];
              hCurve = new Geom_BezierCurve(aPoles, aWeights);
            }
          }
          break;
    
        case BSPLINE:
          {
            Standard_Boolean periodic = Standard_False;
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            if ((ints[0]&4) != 0) periodic = Standard_True;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) len++;
            TColStd_Array1OfReal    aKnots(1, len);
            TColStd_Array1OfInteger aMults(1, len);
            aKnots(1) = data[0];
            aMults(1) = nmult = 1;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) {
                len++;
                aKnots(len) = data[i];
                aMults(len) = nmult = 1;
              } else {
                nmult++;
                aMults(len) = nmult;
              }            
            len = ints[3];
            TColgp_Array1OfPnt aPoles(1, ints[2]);
            for (i = 1; i <= ints[2]; i++, len+=3)
              aPoles(i) = gp_Pnt(data[len], data[len+1], data[len+2]);
            if (rational == 0) {
              hCurve = new Geom_BSplineCurve(aPoles, aKnots, aMults,
                                             ints[1], periodic);
            } else {
              TColStd_Array1OfReal aWeights(1, ints[2]);
              for (i = 1; i <= ints[2]; i++, len++) 
                aWeights(i) = data[len];
              hCurve = new Geom_BSplineCurve(aPoles, aWeights, aKnots, 
                                             aMults, ints[1], periodic);
            }
          }
          break;
    
        case OFFSET:
          {
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            gp_Dir dir(data[0], data[1], data[2]);
            hCurve = new Geom_OffsetCurve(pcurve->handle, data[3], dir);
          }
          break;
      }
    }
    catch (Standard_Failure)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      Handle_Standard_Failure e = Standard_Failure::Caught();
      printf("                %s\n", e->GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      return EGADS_GEOMERR;
    }
    
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_makeGeometry)!\n", stat);
      return stat;
    }
    obj->oclass        = CURVE;
    obj->mtype         = mtype;
    egadsCurve *pcurve = new egadsCurve;
    pcurve->handle     = hCurve;
    pcurve->basis      = basis;
    pcurve->topFlg     = 1;
    obj->blind         = pcurve;
    EG_referenceObject(obj, context);
    if (basis != NULL) EG_referenceTopObj(basis, obj);

  } else {
  
    if ((mtype < PLANE) || (mtype > EXTRUSION)) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface mtype = %d (EG_makeGeometry)!\n", 
               mtype);
      return EGADS_RANGERR;
    }
    if ((mtype == EXTRUSION) || (mtype == REVOLUTION)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is NULL (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != CURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is %d (EG_makeGeometry)!\n", 
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref has no data (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NODATA;
      }
    }  
    if ((mtype == OFFSET) || (mtype == TRIMMED)) {
      if (refGeom == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is NULL (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NULLOBJ;
      }
      if (refGeom->oclass != SURFACE) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref is %d (EG_makeGeometry)!\n", 
                 mtype, refGeom->oclass);
        return EGADS_NOTGEOM;
      }
      basis = refGeom;
      if (basis->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Srf mtype = %d Ref has no data (EG_makeGeometry)!\n", 
                 mtype);
        return EGADS_NODATA;
      }
    }
    
    Handle(Geom_Surface) hSurf;
    try {
      switch (mtype) {

        case PLANE:
          {
            gp_Pnt pntp(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2;
            axi2.SetLocation(pntp);
            axi2.SetDirection(dirz);
            axi2.SetXDirection(dirx);
            axi2.SetYDirection(diry);
            gp_Ax3 axi3(axi2);
            hSurf = new Geom_Plane(axi3);
          }
          break;
          
        case SPHERICAL:
          {
            gp_Pnt pnts(data[0], data[1], data[2]);
            gp_Dir dirx(data[3], data[4], data[5]);
            gp_Dir diry(data[6], data[7], data[8]);
            gp_Dir dirz = dirx.Crossed(diry);
            gp_Ax2 axi2;
            axi2.SetLocation(pnts);
            axi2.SetDirection(dirz);
            axi2.SetXDirection(dirx);
            axi2.SetYDirection(diry);
            gp_Ax3 axi3(axi2);
            hSurf = new Geom_SphericalSurface(axi3, data[9]);
          }
          break;
          
        case CONICAL:
          {
            gp_Pnt pntc(data[0], data[1],  data[2]);
            gp_Dir dirx(data[3], data[4],  data[5]);
            gp_Dir diry(data[6], data[7],  data[8]);
            gp_Dir dirz(data[9], data[10], data[11]);
            gp_Ax3 axi3(pntc, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_ConicalSurface(axi3, data[12], data[13]);
          }
          break;
          
        case CYLINDRICAL:
          {
            gp_Pnt pntc(data[0], data[1],  data[2]);
            gp_Dir dirx(data[3], data[4],  data[5]);
            gp_Dir diry(data[6], data[7],  data[8]);
            gp_Dir dirz(data[9], data[10], data[11]);
            gp_Ax3 axi3(pntc, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_CylindricalSurface(axi3, data[12]);
          }
          break;

        case TOROIDAL:
          {
            gp_Pnt pntt(data[0], data[1],  data[2]);
            gp_Dir dirx(data[3], data[4],  data[5]);
            gp_Dir diry(data[6], data[7],  data[8]);
            gp_Dir dirz(data[9], data[10], data[11]);
            gp_Ax3 axi3(pntt, dirz);
            axi3.SetXDirection(dirx);
            axi3.SetYDirection(diry);
            hSurf = new Geom_ToroidalSurface(axi3, data[12], data[13]);
          }
          break;
        
        case BEZIER:
          {
            rational = 0;
            if ((ints[0]&2) != 0) rational = 1;
            TColgp_Array2OfPnt aPoles(1, ints[2], 1, ints[4]);
            len = 0;
            for (j = 1; j <= ints[4]; j++)
              for (i = 1; i <= ints[2]; i++, len+=3)
                aPoles(i, j) = gp_Pnt(data[len], data[len+1], data[len+2]);
            if (rational == 0) {
              hSurf = new Geom_BezierSurface(aPoles);
            } else {
              TColStd_Array2OfReal aWeights(1, ints[2], 1, ints[4]);
              for (j = 1; j <= ints[4]; j++)
                for (i = 1; i <= ints[2]; i++, len++)
                  aWeights(i, j) = data[len];
              hSurf = new Geom_BezierSurface(aPoles, aWeights);
            }
          }
          break;
          
        case BSPLINE:
          {
            Standard_Boolean uPeriodic = Standard_False;
            Standard_Boolean vPeriodic = Standard_False;
            rational = 0;
            if ((ints[0]&2) != 0) rational  = 1;
            if ((ints[0]&4) != 0) uPeriodic = Standard_True;
            if ((ints[0]&8) != 0) vPeriodic = Standard_True;
//          uKnots
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) len++;
            TColStd_Array1OfReal    uKnots(1, len);
            TColStd_Array1OfInteger uMults(1, len);
            uKnots(1) = data[0];
            uMults(1) = nmult = 1;
            for (len = i = 1; i < ints[3]; i++)
              if (fabs(data[i]-data[i-1]) > KNACC) {
                len++;
                uKnots(len) = data[i];
                uMults(len) = nmult = 1;
              } else {
                nmult++;
                uMults(len) = nmult;
              }
//          vKnots
            for (len = i = 1; i < ints[6]; i++)
              if (fabs(data[ints[3]+i]-data[ints[3]+i-1]) > KNACC) len++;
            TColStd_Array1OfReal    vKnots(1, len);
            TColStd_Array1OfInteger vMults(1, len);
            vKnots(1) = data[ints[3]];
            vMults(1) = nmult = 1;
            for (len = i = 1; i < ints[6]; i++)
              if (fabs(data[ints[3]+i]-data[ints[3]+i-1]) > KNACC) {
                len++;
                vKnots(len) = data[ints[3]+i];
                vMults(len) = nmult = 1;
              } else {
                nmult++;
                vMults(len) = nmult;
              }
            len = ints[3]+ints[6];
            TColgp_Array2OfPnt aPoles(1, ints[2], 1, ints[5]);
            for (j = 1; j <= ints[5]; j++)
              for (i = 1; i <= ints[2]; i++, len+=3)
                aPoles(i, j) = gp_Pnt(data[len], data[len+1], data[len+2]);
            if (rational == 0) {
              hSurf = new Geom_BSplineSurface(aPoles, uKnots, vKnots, 
                                              uMults, vMults, ints[1], ints[4], 
                                              uPeriodic, vPeriodic);
            } else {
              TColStd_Array2OfReal aWeights(1, ints[2], 1, ints[5]);
              for (j = 1; j <= ints[5]; j++)
                for (i = 1; i <= ints[2]; i++, len++)
                  aWeights(i, j) = data[len];
              hSurf = new Geom_BSplineSurface(aPoles, aWeights, uKnots, vKnots, 
                                              uMults, vMults, ints[1], ints[4], 
                                              uPeriodic, vPeriodic);
            }
          }
          break;

        case OFFSET:
          {
            egadsSurface *psurf = (egadsSurface *) basis->blind;
            hSurf = new Geom_OffsetSurface(psurf->handle, data[0]);
          }
          break;
        
        case TRIMMED:
          {
            egadsSurface *psurf = (egadsSurface *) basis->blind;
            hSurf = new Geom_RectangularTrimmedSurface(psurf->handle,
                                        data[0], data[1], data[2], data[3]);
          }
          break;
        
      case EXTRUSION:
          {
            gp_Dir dir(data[0], data[1],  data[2]);
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            hSurf = new Geom_SurfaceOfLinearExtrusion(pcurve->handle, dir);
          }
          break;   

        case REVOLUTION:
          {
            gp_Pnt pnt(data[0], data[1], data[2]);
            gp_Dir dir(data[3], data[4], data[5]);
            gp_Ax1 axi1(pnt, dir);
            egadsCurve *pcurve = (egadsCurve *) basis->blind;
            hSurf = new Geom_SurfaceOfRevolution(pcurve->handle, axi1);
          }
          break;

      }
    }
    catch (Standard_Failure)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      Handle_Standard_Failure e = Standard_Failure::Caught();
      printf("                %s\n", e->GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_makeGeometry)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Surface = %d (EG_makeGeometry)!\n", stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = mtype;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurf;
    psurf->basis        = basis;
    psurf->topFlg       = 1;
    obj->blind          = psurf;
    EG_referenceObject(obj, context);
    if (basis != NULL) EG_referenceTopObj(basis, obj);

  }

  *geom = obj;
  return EGADS_SUCCESS;
}


int
EG_getRange(const egObject *geom, double *range, int *periodic)
{
  int per;

  *periodic = 0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  if (geom->oclass == PCURVE) {
  
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    if (hCurve->IsPeriodic()) *periodic = 1;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
  
  } else if (geom->oclass == CURVE) {
    
    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    if (hCurve->IsPeriodic()) *periodic = 1;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();

  } else if (geom->oclass == SURFACE) {
  
    egadsSurface *psurf = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
    per = 0;
    if (hSurf->IsUPeriodic()) per  = 1;
    if (hSurf->IsVPeriodic()) per |= 2;
    *periodic = per;
    hSurf->Bounds(range[0],range[1], range[2],range[3]);
    
  } else if (geom->oclass == EDGE) {
  
    egadsEdge *pedge = (egadsEdge *) geom->blind;
    BRep_Tool::Range(pedge->edge, range[0],range[1]);
    BRepAdaptor_Curve aCurv(pedge->edge);
    if (aCurv.IsPeriodic()) *periodic = 1;
  
  } else {

    egadsFace *pface = (egadsFace *) geom->blind;
    BRepTools::UVBounds(pface->face, range[0],range[1], 
                                     range[2],range[3]);
    BRepAdaptor_Surface aSurf(pface->face, Standard_True);
    per = 0;
    if (aSurf.IsUPeriodic()) per  = 1;
    if (aSurf.IsVPeriodic()) per |= 2;
    *periodic = per;

  }
  
  return EGADS_SUCCESS;
}


int
EG_curvature(const egObject *geom, const double *param,
             double *result)
{
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  int outLevel = EG_outLevel(geom);
  
  if (geom->oclass == PCURVE) {
    
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    Geom2dLProp_CLProps2d aProp(hCurve, param[0], 2, Precision::Confusion());
    if (aProp.IsTangentDefined()) {
      gp_Dir2d tang;
      
      aProp.Tangent(tang);
      result[0] = aProp.Curvature();
      result[1] = tang.X();
      result[2] = tang.Y();
    } else {
      for (int i = 0; i < 3; i++) result[i] = 0.0;
      return EGADS_DEGEN;
    }
    
  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {
    
    // 1D -- curves & Edges
    Handle(Geom_Curve) hCurve;
    if (geom->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) geom->blind;
      hCurve = pcurve->handle;
    } else {
      egadsEdge *pedge = (egadsEdge *) geom->blind;
      egObject   *curv = pedge->curve;
      if (curv == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Object for Edge (EG_curvature)!\n");
        return EGADS_NULLOBJ;
      }
      egadsCurve *pcurve = (egadsCurve *) curv->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Data for Edge (EG_curvature)!\n");
        return EGADS_NODATA;
      }
      hCurve = pcurve->handle;
    }
    GeomLProp_CLProps aProp(hCurve, param[0], 2, Precision::Confusion());
    if (aProp.IsTangentDefined()) {
      gp_Dir tang;
  
      aProp.Tangent(tang);
      result[0] = aProp.Curvature();
      result[1] = tang.X();
      result[2] = tang.Y();
      result[3] = tang.Z();
    } else {
      for (int i = 0; i < 4; i++) result[i] = 0.0;
      return EGADS_DEGEN;
    }

  } else {
    
    // 2D -- surfaces & Faces
    Handle(Geom_Surface) hSurface;
    if (geom->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      hSurface = psurf->handle;
    } else {
      egadsFace *pface = (egadsFace *) geom->blind;
      egObject *surf = pface->surface;
      if (surf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Object for Face (EG_curvature)!\n");
        return EGADS_NULLOBJ;
      }
      egadsSurface *psurf = (egadsSurface *) surf->blind;
      if (psurf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Data for Face (EG_curvature)!\n");
        return EGADS_NODATA;
      }
      hSurface = psurf->handle;
    }
    GeomLProp_SLProps aProp(hSurface, param[0], param[1], 2,
                            Precision::Angular());
    if (aProp.IsCurvatureDefined()) {
      gp_Dir MaxD, MinD;
      
      aProp.CurvatureDirections(MaxD, MinD);
      result[0] = aProp.MaxCurvature();
      result[1] = MaxD.X();
      result[2] = MaxD.Y();
      result[3] = MaxD.Z();
      result[4] = aProp.MinCurvature();
      result[5] = MinD.X();
      result[6] = MinD.Y();
      result[7] = MinD.Z();
    } else {
      for (int i = 0; i < 8; i++) result[i] = 0.0;
      return EGADS_DEGEN;
    }
    
  }
  
  return EGADS_SUCCESS;
}


int
EG_evaluate(const egObject *geom, const double *param,
                                        double *result)
{
  int    outLevel;
  gp_Pnt P0;
  gp_Vec V1, V2, U1, U2, UV;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  
  if (geom->oclass == PCURVE) {
    gp_Pnt2d P2d;
    gp_Vec2d V12d, V22d;

    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    hCurve->D2(*param, P2d, V12d, V22d);
    result[0] = P2d.X();
    result[1] = P2d.Y();
    result[2] = V12d.X();
    result[3] = V12d.Y();
    result[4] = V22d.X();
    result[5] = V22d.Y();
  
  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {
  
    // 1D -- curves & Edges
    if (geom->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) geom->blind;
      Handle(Geom_Curve) hCurve = pcurve->handle;
      hCurve->D2(*param, P0, V1, V2);
    } else {
      egadsEdge *pedge = (egadsEdge *) geom->blind;
#ifdef ADAPTOR
      BRepAdaptor_Curve aCurv(pedge->edge);
      aCurv.D2(*param, P0, V1, V2);
#else
      egObject *curv = pedge->curve;
      if (curv == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Object for Edge (EG_evaluate)!\n");
        return EGADS_NULLOBJ;
      }
      egadsCurve *pcurve = (egadsCurve *) curv->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Data for Edge (EG_evaluate)!\n");
        return EGADS_NODATA;
      }
      Handle(Geom_Curve) hCurve = pcurve->handle;
      hCurve->D2(*param, P0, V1, V2);
#endif
    }
    result[0] = P0.X();
    result[1] = P0.Y();
    result[2] = P0.Z();
    result[3] = V1.X();
    result[4] = V1.Y();
    result[5] = V1.Z();
    result[6] = V2.X();
    result[7] = V2.Y();
    result[8] = V2.Z();
  
  } else {
  
    // 2D -- surfaces & Faces
    if (geom->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      Handle(Geom_Surface) hSurface = psurf->handle;
      hSurface->D2(param[0], param[1], P0, U1, V1, U2, V2, UV);
    } else {
      egadsFace *pface = (egadsFace *) geom->blind;
#ifdef ADAPTOR
      BRepAdaptor_Surface aSurf(pface->face, Standard_True);
      aSurf.D2(param[0], param[1], P0, U1, V1, U2, V2, UV);
#else
      egObject *surf = pface->surface;
      if (surf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Object for Face (EG_evaluate)!\n");
        return EGADS_NULLOBJ;
      }
      egadsSurface *psurf = (egadsSurface *) surf->blind;
      if (psurf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Data for Face (EG_evaluate)!\n");
        return EGADS_NODATA;
      }
      Handle(Geom_Surface) hSurface = psurf->handle;
      hSurface->D2(param[0], param[1], P0, U1, V1, U2, V2, UV);
#endif
    }
    result[ 0] = P0.X();
    result[ 1] = P0.Y();
    result[ 2] = P0.Z();
    result[ 3] = U1.X();
    result[ 4] = U1.Y();
    result[ 5] = U1.Z();
    result[ 6] = V1.X();
    result[ 7] = V1.Y();
    result[ 8] = V1.Z();
    result[ 9] = U2.X();
    result[10] = U2.Y();
    result[11] = U2.Z();
    result[12] = UV.X();
    result[13] = UV.Y();
    result[14] = UV.Z();
    result[15] = V2.X();
    result[16] = V2.Y();
    result[17] = V2.Z();
    
  }
  
  return EGADS_SUCCESS;
}


static void
EG_nearestPCurve(Handle_Geom2d_Curve hCurve, const double *coor,
                 double tmin, double tmax, int flag, double *t, double *uv)
{
  int      i;
  double   a, b, tx, pw[2];
  gp_Pnt2d pnt;
  gp_Vec2d t1, t2;
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};
  
  // sample and pick closest
  if (flag == 0) {
    b = 0.0;
    for (i = 0; i < 5; i++) {
      tx = (1.0-ratios[i])*tmin + ratios[i]*tmax;
      hCurve->D0(tx, pnt);
      a = (pnt.X()-coor[0])*(pnt.X()-coor[0]) +
          (pnt.Y()-coor[1])*(pnt.Y()-coor[1]);
      if (i == 0) {
        *t = tx;
        b  = a;
      } else {
        if (a < b) {
          *t = tx;
          b  = a;
        }
      }
    }
  }
  
  // netwon-raphson from picked position
  for (i = 0; i < 20; i++) {
    if ((*t < tmin) || (*t > tmax)) break;
    hCurve->D2(*t, pnt, t1, t2);
    pw[0] = pnt.X() - coor[0];
    pw[1] = pnt.Y() - coor[1];
    b     = -( pw[0]*t1.X() +  pw[1]*t1.Y());
    a     =  (t1.X()*t1.X() + t1.Y()*t1.Y()) +
             ( pw[0]*t2.X() +  pw[1]*t2.Y());
    if (a == 0.0) break;
    b  /= a;
//  if (fabs(b) < 1.e-10*(tmax-tmin)) break;
    *t += b;
  }
  if (*t < tmin) *t = tmin;
  if (*t > tmax) *t = tmax;
  
  hCurve->D0(*t, pnt);
  uv[0] = pnt.X();
  uv[1] = pnt.Y();
}


void
EG_nearestCurve(Handle_Geom_Curve hCurve, const double *coor,
                double tmin, double tmax, int flag, double *t, double *xyz)
{
  int    i;
  double a, b, tx, pw[3];
  gp_Pnt pnt;
  gp_Vec t1, t2;
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};

  // sample and pick closest
  if (flag == 0) {
    b = 0.0;
    for (i = 0; i < 5; i++) {
      tx = (1.0-ratios[i])*tmin + ratios[i]*tmax;
      hCurve->D0(tx, pnt);
      a = (pnt.X()-coor[0])*(pnt.X()-coor[0]) +
          (pnt.Y()-coor[1])*(pnt.Y()-coor[1]) +
          (pnt.Z()-coor[2])*(pnt.Z()-coor[2]);
      if (i == 0) {
        *t = tx;
        b  = a;
      } else {
        if (a < b) {
          *t = tx;
          b  = a;
        }
      }
    }
  }
  
  // netwon-raphson from picked position
  for (i = 0; i < 20; i++) {
    if ((*t < tmin) || (*t > tmax)) break;
    hCurve->D2(*t, pnt, t1, t2);
    pw[0] = pnt.X() - coor[0];
    pw[1] = pnt.Y() - coor[1];
    pw[2] = pnt.Z() - coor[2];
    b     = -( pw[0]*t1.X() +  pw[1]*t1.Y() +  pw[2]*t1.Z());
    a     =  (t1.X()*t1.X() + t1.Y()*t1.Y() + t1.Z()*t1.Z()) +
             ( pw[0]*t2.X() +  pw[1]*t2.Y() +  pw[2]*t2.Z());
    if (a == 0.0) break;
    b  /= a;
//  if (fabs(b) < 1.e-10*(tmax-tmin)) break;
    *t += b;
  }
  if (*t < tmin) *t = tmin;
  if (*t > tmax) *t = tmax;
  
  hCurve->D0(*t, pnt);
  xyz[0] = pnt.X();
  xyz[1] = pnt.Y();
  xyz[2] = pnt.Z();
}


static int
EG_nearestSurface(Handle_Geom_Surface hSurface, double *range,
                  const double *point, int flag, double *uv, double *coor)
{
  int    i, j, count;
  gp_Pnt P0;
  gp_Vec V1, V2, U1, U2, UV;
  double a00, a10, a11, b0, b1, det, dist, ldist, dx[3], uvs[2];
  static double ratios[5] = {0.02, 0.25, 0.5, 0.75, 0.98};
  
  if (flag == 0) {
    // find candidate starting point
    ldist = 1.e308;
    for (j = 0; j < 5; j++) {
      uvs[1] = (1.0-ratios[j])*range[2] + ratios[j]*range[3];
      for (i = 0; i < 5; i++) {
        uvs[0] = (1.0-ratios[i])*range[0] + ratios[i]*range[1];
        hSurface->D0(uv[0], uv[1], P0);
        dx[0] = P0.X() - point[0];
        dx[1] = P0.Y() - point[1];
        dx[2] = P0.Z() - point[2];
        dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        if (dist < ldist) {
          ldist = dist;
          uv[0] = uvs[0];
          uv[1] = uvs[1];
        }
      }
    }
  }

  // newton iteration
  for (count = 0; count < 10; count++) {
    hSurface->D2(uv[0], uv[1], P0, U1, V1, U2, V2, UV);
    dx[0] = P0.X() - point[0];
    dx[1] = P0.Y() - point[1];
    dx[2] = P0.Z() - point[2];
    dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    if (dist < Precision::Confusion()) break;
    if (count != 0) {
      if (fabs(dist-ldist) < Precision::Confusion()) break;
      if (dist > ldist) {
        uv[0] = uvs[0];
        uv[1] = uvs[1];
        hSurface->D0(uv[0], uv[1], P0);
        coor[0] = P0.X();
        coor[1] = P0.Y();
        coor[2] = P0.Z();
        return EGADS_EMPTY;
      }
    }

    b0  = -dx[0]*U1.X() -  dx[1]*U1.Y() -  dx[2]*U1.Z();
    b1  = -dx[0]*V1.X() -  dx[1]*V1.Y() -  dx[2]*V1.Z();
    a00 = U1.X()*U1.X() + U1.Y()*U1.Y() + U1.Z()*U1.Z() +
           dx[0]*U2.X() +  dx[1]*U2.Y() +  dx[2]*U2.Z();
    a10 = U1.X()*V1.X() + U1.Y()*V1.Y() + U1.Z()*V1.Z() +
           dx[0]*UV.X() +  dx[1]*UV.Y() +  dx[2]*UV.Z();
    a11 = V1.X()*V1.X() + V1.Y()*V1.Y() + V1.Z()*V1.Z() +
           dx[0]*V2.X() +  dx[1]*V2.Y() +  dx[2]*V2.Z();

    det    = a00*a11 - a10*a10;
    if (det == 0.0) return EGADS_DEGEN;
    det    = 1.0/det;
    uvs[0] = uv[0];
    uvs[1] = uv[1];
    uv[0] += det*(b0*a11 - b1*a10);
    uv[1] += det*(b1*a00 - b0*a10);
    ldist  = dist;
//  printf("   %d: %lf %lf   %le\n", count, uv[0], uv[1], ldist);
  }

  hSurface->D0(uv[0], uv[1], P0);
  coor[0] = P0.X();
  coor[1] = P0.Y();
  coor[2] = P0.Z();
  if (count == 10) return EGADS_EMPTY;

  return EGADS_SUCCESS;
}


int
EG_invEvaluate(const egObject *geom, double *xyz, 
               double *param, double *result)
{
  int           outLevel, stat;
  Standard_Real period, t, u, v, range[4], coor[3];

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  
  if (geom->oclass == PCURVE) {
  
    // 2D on PCurves
    gp_Pnt2d pnt(xyz[0], xyz[1]);
  
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
    Geom2dAPI_ProjectPointOnCurve projPnt(pnt, hCurve);
    if (projPnt.NbPoints() == 0) {
      EG_nearestPCurve(hCurve, xyz, range[0], range[1], 0, &t, result);
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
    } else {
      pnt = projPnt.NearestPoint();
      t   = projPnt.LowerDistanceParameter();
    }
    if (hCurve->IsPeriodic()) {
      period = hCurve->Period();
      if ((t+PARAMACC < range[0]) || (t-PARAMACC > range[1]))
        if (period != 0.0)
          if (t+PARAMACC < range[0]) {
            if (t+period-PARAMACC < range[1]) t += period;
          } else {
            if (t-period+PARAMACC > range[0]) t -= period;
          }
    }
    result[0] = pnt.X();
    result[1] = pnt.Y();
    *param    = t;
    return EGADS_SUCCESS;
  }
  
  // make the point
  gp_Pnt pnt(xyz[0], xyz[1], xyz[2]);
  
  if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {
  
    // 1D -- curves & Edges
    Handle(Geom_Curve) hCurve;
    if (geom->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) geom->blind;
      hCurve   = pcurve->handle;
      range[0] = hCurve->FirstParameter();
      range[1] = hCurve->LastParameter();
    } else {
      egadsEdge *pedge = (egadsEdge *) geom->blind;
      BRep_Tool::Range(pedge->edge, range[0], range[1]);
      egObject  *curv  = pedge->curve;
      if (curv == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Object for Edge (EG_invEvaluate)!\n");
        return EGADS_NULLOBJ;
      }
      egadsCurve *pcurve = (egadsCurve *) curv->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Data for Edge (EG_invEvaluate)!\n");
        return EGADS_NODATA;
      }
      hCurve = pcurve->handle;
    }

    GeomAPI_ProjectPointOnCurve projPnt(pnt, hCurve);
    if (projPnt.NbPoints() == 0) {
      EG_nearestCurve(hCurve, xyz, range[0], range[1], 0, &t, result);
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
      pnt.SetZ(result[2]);
    } else {
      pnt = projPnt.NearestPoint();
      t   = projPnt.LowerDistanceParameter();
    }

    if (hCurve->IsPeriodic()) {
      period = hCurve->Period();
      if ((t+PARAMACC < range[0]) || (t-PARAMACC > range[1])) {
        if (period != 0.0)
          if (t+PARAMACC < range[0]) {
            if (t+period-PARAMACC < range[1]) t += period;
          } else {
            if (t-period+PARAMACC > range[0]) t -= period;
          }
      }
    }
    
    /* clip it? */
    if (geom->oclass == EDGE)
      if ((t < range[0]) || (t > range[1])) {
/*      if (t < range[0]) t = range[0];
        if (t > range[1]) t = range[1];
        hCurve->D0(t, pnt);  */
        EG_nearestCurve(hCurve, xyz, range[0], range[1], 0, &t, result);
        pnt.SetX(result[0]);
        pnt.SetY(result[1]);
        pnt.SetZ(result[2]);
      }

    result[0] = pnt.X();
    result[1] = pnt.Y();
    result[2] = pnt.Z();
    *param    = t;

  } else {
  
    // 2D -- surfaces & Faces
    Handle(Geom_Surface) hSurface;
    if (geom->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      hSurface = psurf->handle;
      hSurface->Bounds(range[0],range[1], range[2],range[3]);
    } else {
      egadsFace *pface = (egadsFace *) geom->blind;
      BRepTools::UVBounds(pface->face, range[0],range[1], 
                                       range[2],range[3]);
      egObject  *surf  = pface->surface;
      if (surf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Object for Face (EG_invEvaluate)!\n");
        return EGADS_NULLOBJ;
      }
      egadsSurface *psurf = (egadsSurface *) surf->blind;
      if (psurf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Data for Face (EG_invEvaluate)!\n");
        return EGADS_NODATA;
      }
      hSurface = psurf->handle;
    }
    
    GeomAPI_ProjectPointOnSurf projPnt(pnt, hSurface);
    if (!projPnt.IsDone()) {
      stat = EG_nearestSurface(hSurface, range, xyz, 0, param, result);
      if (stat == EGADS_DEGEN) {
        if (outLevel > 0)
          printf(" EGADS Warning: Surf Proj Incomplete - DEGEN (EG_invEvaluate)!\n");
        return stat;
      } else if (stat == EGADS_EMPTY) {
        if (outLevel > 1)
          printf(" EGADS Warning: Surf Proj Incomplete (EG_invEvaluate)!\n");
      }
      u  = param[0];
      v  = param[1];
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
      pnt.SetZ(result[2]);
    } else {
      if (projPnt.NbPoints() == 0) {
        if (outLevel > 0)
          printf(" EGADS Warning: No projection on Surf (EG_invEvaluate)!\n");
        return EGADS_NOTFOUND;
      }
      pnt = projPnt.NearestPoint();
      projPnt.LowerDistanceParameters(u, v);
    }

    if (hSurface->IsUPeriodic()) {
      period = hSurface->UPeriod();
      if ((u+PARAMACC < range[0]) || (u-PARAMACC > range[1])) {
        if (period != 0.0)
          if (u+PARAMACC < range[0]) {
            if (u+period-PARAMACC < range[1]) u += period;
          } else {
            if (u-period+PARAMACC > range[0]) u -= period;
          }
      }
    }
    if (hSurface->IsVPeriodic()) {
      period = hSurface->VPeriod();
      if ((v+PARAMACC < range[2]) || (v-PARAMACC > range[3])) {
        if (period != 0.0)
          if (v+PARAMACC < range[2]) {
            if (v+period-PARAMACC < range[3]) v += period;
          } else {
            if (v-period+PARAMACC > range[2]) v -= period;
          }
      }
    }

    if (geom->oclass == FACE) {
      egadsFace *pface = (egadsFace *) geom->blind;
      Standard_Real tol = BRep_Tool::Tolerance(pface->face);
      gp_Pnt2d pnt2d(u, v);
      TopOpeBRep_PointClassifier pClass;
      pClass.Load(pface->face);
      if (pClass.Classify(pface->face, pnt2d, tol) == TopAbs_OUT) {
        /* find closest clipped point on the edges */
        double dist2 = 1.e308;
        gp_Pnt pnts(xyz[0], xyz[1], xyz[2]);
        gp_Pnt pntt(xyz[0], xyz[1], xyz[2]);
        TopExp_Explorer ExpW;
        for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
          TopoDS_Shape shapw = ExpW.Current();
          TopoDS_Wire  wire  = TopoDS::Wire(shapw);
          BRepTools_WireExplorer ExpWE;
          for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
            TopoDS_Shape shape = ExpWE.Current();
            TopoDS_Edge  wedge = TopoDS::Edge(shape);
            if (BRep_Tool::Degenerated(wedge)) continue;
            Standard_Real t1, t2;
            Handle(Geom_Curve) hCurve = BRep_Tool::Curve(wedge, t1, t2);
            GeomAPI_ProjectPointOnCurve projPnt(pnts, hCurve);
            if (projPnt.NbPoints() == 0) {
              EG_nearestCurve(hCurve, xyz, t1, t2, 0, &t, result);
              pnt.SetX(result[0]);
              pnt.SetY(result[1]);
              pnt.SetZ(result[2]);
            } else {
              pnt = projPnt.NearestPoint();
              t   = projPnt.LowerDistanceParameter();
            }
            if ((t < t1) || (t > t2)) {
/*            if (t < t1) t = t1;
              if (t > t2) t = t2;
              hCurve->D0(t, pnt);  */
              EG_nearestCurve(hCurve, xyz, t1, t2, 0, &t, result);
              pnt.SetX(result[0]);
              pnt.SetY(result[1]);
              pnt.SetZ(result[2]);
            }
            double d = (pnts.X()-pnt.X())*(pnts.X()-pnt.X()) +
                       (pnts.Y()-pnt.Y())*(pnts.Y()-pnt.Y()) +
                       (pnts.Z()-pnt.Z())*(pnts.Z()-pnt.Z());
            if (d < dist2) {
              pntt  = pnt;
              dist2 = d;
            }
          }
        }
        // project again but with clipped point
        GeomAPI_ProjectPointOnSurf projPnt(pntt, hSurface);
        if (!projPnt.IsDone()) {
          coor[0] = pntt.X();
          coor[1] = pntt.Y();
          coor[2] = pntt.Z();
          stat = EG_nearestSurface(hSurface, range, coor, 0, param, result);
          if (stat == EGADS_DEGEN) {
            if (outLevel > 0)
              printf(" EGADS Warning: Face Proj Incomplete - DEGEN (EG_invEvaluate)!\n");
            return stat;
          } else if (stat == EGADS_EMPTY) {
            if (outLevel > 1)
              printf(" EGADS Warning: Face Proj Incomplete (EG_invEvaluate)!\n");
          }
          u  = param[0];
          v  = param[1];
          pnt.SetX(result[0]);
          pnt.SetY(result[1]);
          pnt.SetZ(result[2]);
        } else {
          if (projPnt.NbPoints() == 0) {
            if (outLevel > 0)
              printf(" EGADS Warning: No projection on Face (EG_invEvaluate)!\n");
            return EGADS_NOTFOUND;
          }
          pnt = projPnt.NearestPoint();
          projPnt.LowerDistanceParameters(u, v);
        }

        if (hSurface->IsUPeriodic()) {
          period = hSurface->UPeriod();
          if ((u+PARAMACC < range[0]) || (u-PARAMACC > range[1])) {
            if (period != 0.0)
              if (u+PARAMACC < range[0]) {
                if (u+period-PARAMACC < range[1]) u += period;
              } else {
                if (u-period+PARAMACC > range[0]) u -= period;
              }
          }
        }
        if (hSurface->IsVPeriodic()) {
          period = hSurface->VPeriod();
          if ((v+PARAMACC < range[2]) || (v-PARAMACC > range[3])) {
            if (period != 0.0)
              if (v+PARAMACC < range[2]) {
                if (v+period-PARAMACC < range[3]) v += period;
              } else {
                if (v-period+PARAMACC > range[2]) v -= period;
              }
          }
        }
             
      }
    }

    result[0] = pnt.X();
    result[1] = pnt.Y();
    result[2] = pnt.Z();
    param[0]  = u;
    param[1]  = v;

  }
  
  return EGADS_SUCCESS;
}


int
EG_invEvaluateGuess(const egObject *geom, double *xyz,
                    double *param, double *result)
{
  int           outLevel;
  Standard_Real range[4];
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE) &&
      (geom->oclass != EDGE)   && (geom->oclass != FACE))
                                   return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(geom);
  
  if (geom->oclass == PCURVE) {
    
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
    EG_nearestPCurve(hCurve, xyz, range[0], range[1], 1, param, result);

  } else if ((geom->oclass == CURVE) || (geom->oclass == EDGE)) {
  
    Handle(Geom_Curve) hCurve;
    if (geom->oclass == CURVE) {
      egadsCurve *pcurve = (egadsCurve *) geom->blind;
      hCurve   = pcurve->handle;
      range[0] = hCurve->FirstParameter();
      range[1] = hCurve->LastParameter();
    } else {
      egadsEdge *pedge = (egadsEdge *) geom->blind;
      BRep_Tool::Range(pedge->edge, range[0], range[1]);
      egObject  *curv  = pedge->curve;
      if (curv == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Object for Edge (EG_invEvaluate)!\n");
        return EGADS_NULLOBJ;
      }
      egadsCurve *pcurve = (egadsCurve *) curv->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No curve Data for Edge (EG_invEvaluate)!\n");
        return EGADS_NODATA;
      }
      hCurve = pcurve->handle;
    }
    EG_nearestCurve(hCurve, xyz, range[0], range[1], 1, param, result);

  } else {
    // 2D -- surfaces & Faces
    Handle(Geom_Surface) hSurface;
    if (geom->oclass == SURFACE) {
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      hSurface = psurf->handle;
      hSurface->Bounds(range[0],range[1], range[2],range[3]);
    } else {
      egadsFace *pface = (egadsFace *) geom->blind;
      BRepTools::UVBounds(pface->face, range[0],range[1],
                          range[2],range[3]);
      egObject  *surf  = pface->surface;
      if (surf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Object for Face (EG_invEvaluate)!\n");
        return EGADS_NULLOBJ;
      }
      egadsSurface *psurf = (egadsSurface *) surf->blind;
      if (psurf == NULL) {
        if (outLevel > 0)
          printf(" EGADS Warning: No Surf Data for Face (EG_invEvaluate)!\n");
        return EGADS_NODATA;
      }
      hSurface = psurf->handle;
    }
    return EG_nearestSurface(hSurface, range, xyz, 1, param, result);
  }

  return EGADS_SUCCESS;
}



int
EG_arcLength(const egObject *geom, double t1, double t2, double *alen)
{

  *alen = 0.0;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != EDGE))      return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  if (geom->oclass == PCURVE) {
    
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    if (ppcurv == NULL) return EGADS_NULLOBJ;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    Geom2dAdaptor_Curve AC(hCurve);
    *alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);

  } else if (geom->oclass == CURVE) {
  
    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    if (pcurve == NULL) return EGADS_NULLOBJ;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAdaptor_Curve AC(hCurve);
    *alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);

  } else {

    if (geom->mtype == DEGENERATE) return EGADS_SUCCESS;
    egadsEdge *pedge  = (egadsEdge *) geom->blind;
    if (pedge  == NULL) return EGADS_NULLOBJ;
    egObject  *curvo  = pedge->curve;
    if (curvo  == NULL) return EGADS_NULLOBJ;
    egadsCurve *pcurve = (egadsCurve *) curvo->blind;
    if (pcurve == NULL) return EGADS_NULLOBJ;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAdaptor_Curve AC(hCurve);
    *alen = GCPnts_AbscissaPoint::Length(AC, t1, t2);

  }

  return EGADS_SUCCESS;
}


int
EG_approximate(egObject *context, int maxdeg, double tol, const int *sizes,
               const double *data, egObject **bspline)
{
  int      i, outLevel, stat, len = 0;
  egObject *obj;

  *bspline = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  outLevel = EG_outLevel(context);
  
  if ((maxdeg < 0) || (maxdeg > 8)) {
    if (outLevel > 0)
      printf(" EGADS Warning: maxDeg = %d (EG_approximate)!\n", maxdeg);
    return EGADS_RANGERR;
  }
  
  if (sizes[1] == -1)
    if ((maxdeg < 3) && (sizes[0] > 2))
      return EG_spline1d(context, maxdeg, -sizes[0], data, tol, bspline);

  if (sizes[1] == 0) {

    if ((maxdeg < 3) && (sizes[0] > 2))
      return EG_spline1d(context, maxdeg, sizes[0], data, tol, bspline);

    // curve
    Handle(Geom_BSplineCurve) hCurve;
    if (sizes[0] == 2) {
      TColStd_Array1OfReal    aKnots(1, sizes[0]);
      TColStd_Array1OfInteger aMults(1, sizes[0]);
      TColgp_Array1OfPnt      aPoles(1, sizes[0]);
      aKnots(1) = 0.0;
      aMults(1) = 2;
      aKnots(2) = 1.0;
      aMults(2) = 2;
      for (i = 0; i < sizes[0]; i++)
        aPoles(i+1) = gp_Pnt(data[3*i], data[3*i+1], data[3*i+2]);
      hCurve = new Geom_BSplineCurve(aPoles, aKnots, aMults, 1);
      if (hCurve.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Null Curve (EG_approximate)!\n");
        return EGADS_GEOMERR;
      }
    } else {
      try {
        TColgp_Array1OfPnt aPnts(1, sizes[0]);
        for (int i = 1; i <= sizes[0]; i++, len+=3)
          aPnts(i) = gp_Pnt(data[len], data[len+1], data[len+2]);
        hCurve = GeomAPI_PointsToBSpline(aPnts, 3, maxdeg, GeomAbs_C2,
                                         tol).Curve();
      }
      catch (Standard_Failure)
      {
        if (outLevel > 0) {
          printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
          Handle_Standard_Failure e = Standard_Failure::Caught();
          printf("                %s\n", e->GetMessageString());
        }
        return EGADS_GEOMERR;
      }
      catch (...)
      {
        if (outLevel > 0)
          printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
        return EGADS_GEOMERR;
      }
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: make Curve = %d (EG_approximate)!\n",
               stat);
      return stat;
    }
    obj->oclass        = CURVE;
    obj->mtype         = BSPLINE;
    egadsCurve *pcurve = new egadsCurve;
    pcurve->handle     = hCurve;
    pcurve->basis      = NULL;
    pcurve->topFlg     = 0;
    obj->blind         = pcurve;
    EG_referenceObject(obj, context);

  } else if ((sizes[0] <= 2) || (sizes[1] < 1)) {
    
    if (outLevel > 0)
      printf(" EGADS Error: Sizes = %d %d (EG_approximate)!\n",
             sizes[0], sizes[1]);
    return EGADS_RANGERR;

  } else {
  
    if (maxdeg < 3)
      return EG_spline2d(context, maxdeg, NULL, sizes[0], sizes[1], data, tol,
                         bspline);

    // surface
    Handle(Geom_BSplineSurface) hSurf;
    try {    
      TColgp_Array2OfPnt aPnts(1, sizes[0], 1, sizes[1]);
      for (int j = 1; j <= sizes[1]; j++)
        for (int i = 1; i <= sizes[0]; i++, len+=3)
          aPnts(i, j) = gp_Pnt(data[len], data[len+1], data[len+2]);
      if (tol != 0.0) {
        hSurf = GeomAPI_PointsToBSplineSurface(aPnts, 3, maxdeg, GeomAbs_C2, 
                                               tol).Surface();
      } else {
        GeomAPI_PointsToBSplineSurface P2BSpl;
        P2BSpl.Interpolate(aPnts);
        hSurf = P2BSpl.Surface();
      }
    }
    catch (Standard_Failure)
    {
      if (outLevel > 0) {
        printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
        Handle_Standard_Failure e = Standard_Failure::Caught();
        printf("                %s\n", e->GetMessageString());
      }
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      if (outLevel > 0)
        printf(" EGADS Warning: Internal Error (EG_approximate)!\n");
      return EGADS_GEOMERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: make Surface = %d (EG_approximate)!\n",
               stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = BSPLINE;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurf;
    psurf->basis        = NULL;
    psurf->topFlg       = 0;
    obj->blind          = psurf;
    EG_referenceObject(obj, context);

  }
  
  *bspline = obj;
  return EGADS_SUCCESS;
}


int
EG_otherCurve(const egObject *surface, const egObject *curve, 
              double tol, egObject **newcurve)
{
  int      outLevel, stat;
  egObject *context, *obj;

  *newcurve = NULL;
  if (surface == NULL)               return EGADS_NULLOBJ;
  if (surface->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (surface->oclass != SURFACE)    return EGADS_NOTGEOM;
  if (surface->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(surface);
  context  = EG_context(surface);

  if (curve == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Input Curve (EG_otherCurve)!\n"); 
    return EGADS_NULLOBJ;
  }
  if ((curve->oclass != PCURVE) && (curve->oclass != CURVE) &&
      (curve->oclass != EDGE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not a PCurve/Curve or Edge (EG_otherCurve)!\n"); 
    return EGADS_NOTGEOM;
  }
  if (curve->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: PCurve has no data (EG_otherCurve)!\n"); 
    return EGADS_NODATA;
  }
  if (EG_context(curve) != context) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_otherCurve)!\n");
    return EGADS_MIXCNTX;
  }

  egadsSurface *psurf           = (egadsSurface *) surface->blind;
  Handle(Geom_Surface) hSurface = psurf->handle;
  Standard_Real prec            = tol;
  if (prec < Precision::Confusion()) prec = Precision::Confusion();
  
  if (curve->oclass == PCURVE) {
  
    Standard_Real maxDev, aveDev;
  
    egadsPCurve *ppcurv         = (egadsPCurve *) curve->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    GeomAdaptor_Surface  aGAS   = hSurface;
    Handle(GeomAdaptor_HSurface) aHGAS = new GeomAdaptor_HSurface(aGAS);
    Handle(Geom2dAdaptor_HCurve) Crv   = new Geom2dAdaptor_HCurve(hCurve);
    Adaptor3d_CurveOnSurface ConS(Crv,aHGAS);
    
    Handle(Geom_Curve) newcrv;
    signal(SIGSEGV, segfault_handler);
    switch (stat = setjmp(jmpenv)) {
      case 0:
        GeomLib::BuildCurve3d(prec, ConS, hCurve->FirstParameter(),
                              hCurve->LastParameter(), newcrv, maxDev, aveDev);
        break;
      default:
        printf(" EGADS Fatal Error: OCC SegFault %d (EG_otherCurve)!\n", stat);
        signal(SIGSEGV, SIG_DFL);
        return EGADS_OCSEGFLT;
    }
    signal(SIGSEGV, SIG_DFL);

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_otherCurve)!\n", stat);
      return stat;
    }
    EG_completeCurve(obj, newcrv);

  } else {
  
    Handle(Geom2d_Curve) newcrv;

    if (curve->oclass == EDGE) {
    
      Standard_Real t1, t2;

      egadsEdge *pedge = (egadsEdge *) curve->blind;
      double toler     = BRep_Tool::Tolerance(pedge->edge);
      if (prec < toler) prec = toler;
      egObject *geom   = pedge->curve;
      if (geom->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: NULL Curve Data (EG_otherCurve)!\n");
        return EGADS_NODATA;
      }
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(pedge->edge, t1, t2);
      try {
        newcrv = GeomProjLib::Curve2d(hCurve, t1, t2, hSurface, prec);
      }
      catch (Standard_Failure)
      {
        printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
        Handle_Standard_Failure e = Standard_Failure::Caught();
        printf("                %s\n", e->GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...)
      {
        printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
        return EGADS_GEOMERR;
      }

    } else {
    
      egadsCurve *pcurve        = (egadsCurve *) curve->blind;
      Handle(Geom_Curve) hCurve = pcurve->handle;
      try {
        newcrv = GeomProjLib::Curve2d(hCurve, hCurve->FirstParameter(),
                                      hCurve->LastParameter(), hSurface, prec);
      }
      catch (Standard_Failure)
      {
        printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
        Handle_Standard_Failure e = Standard_Failure::Caught();
        printf("                %s\n", e->GetMessageString());
        return EGADS_GEOMERR;
      }
      catch (...)
      {
        printf(" EGADS Warning: Geometry Creation Error (EG_otherCurve)!\n");
        return EGADS_GEOMERR;
      }

    }

    if (EG_getPCurveType(newcrv) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot construct PCurve (EG_otherCurve)!\n");
      return EGADS_CONSTERR;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make PCurve = %d (EG_otherCurve)!\n", stat);
      return stat;
    }
    EG_completePCurve(obj, newcrv);
  }
  
  EG_referenceObject(obj, context);
  *newcurve = obj;
  return EGADS_SUCCESS;
}


int
EG_isoCline(const egObject *surface, int UV, double value, 
                  egObject **newcurve)
{
  int      stat, outLevel;
  egObject *context, *obj;

  *newcurve = NULL;
  if  (surface == NULL)               return EGADS_NULLOBJ;
  if  (surface->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((surface->oclass != SURFACE) &&
      (surface->oclass != LOOP))      return EGADS_NOTGEOM;
  if  (surface->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(surface);
  context  = EG_context(surface);

  // special loop fitting code
  if (surface->oclass == LOOP) {
    if (surface->mtype != CLOSED) {
      if (outLevel > 1)
        printf(" EGADS Info: Input Loop is NOT closed (EG_isocline)!\n");
      return EGADS_TOPOERR;
    }
    egadsLoop *ploop  = (egadsLoop *) surface->blind;
    if (ploop->surface != NULL) {
      if (outLevel > 1)
        printf(" EGADS Info: Input Loop has attached Surface (EG_isocline)!\n");
      return EGADS_TOPOERR;
    }
    if ((ploop->nedges < 3) || (ploop->nedges > 4)) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop has %d Edges (EG_isocline)!\n",
               ploop->nedges);
      return EGADS_GEOMERR;
    }
    
    Handle(Geom_BSplineCurve) bsc[4];
    ShapeConstruct_Curve      ShapeCC;
    Standard_Real             prec  = Precision::Confusion();
    if (value > prec)         prec  = value;
    GeomFill_FillingStyle     style = GeomFill_StretchStyle;
    if (UV < 0)               style = GeomFill_CurvedStyle;
    if (UV > 0)               style = GeomFill_CoonsStyle;
    for (int i = 0; i < ploop->nedges; i++) {
      obj = ploop->edges[i];
      egadsEdge *pedge = (egadsEdge *) obj->blind;
      if (pedge == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: %d/%d NULL pedge (EG_isocline)!\n",
                 i+1, ploop->nedges);
        return EGADS_NODATA;
      }
      double range[2];
      BRep_Tool::Range(pedge->edge, range[0], range[1]);
      obj = pedge->curve;
      egadsCurve *pcurve = (egadsCurve *) obj->blind;
      if (pcurve == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: %d/%d NULL pcurve (EG_isocline)!\n",
                 i+1, ploop->nedges);
        return EGADS_NODATA;
      }
      Handle(Geom_Curve) hCurve = pcurve->handle;
      bsc[i] = ShapeCC.ConvertToBSpline(hCurve, range[0], range[1], prec);
      if (bsc[i].IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to Convert %d/%d (EG_isocline)!\n",
                 i+1, ploop->nedges);
        return EGADS_GEOMERR;
      }
    }
    GeomFill_BSplineCurves fill;
    if (ploop->nedges == 3) {
      fill.Init(bsc[2], bsc[1], bsc[0], style);
    } else {
      fill.Init(bsc[3], bsc[2], bsc[1], bsc[0], style);
    }
    Handle(Geom_BSplineSurface) hSurface = fill.Surface();
    if (hSurface.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to Construct (EG_isocline)!\n");
      return EGADS_GEOMERR;
    }
    
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Surface = %d (EG_isocline)!\n", stat);
      return stat;
    }
    obj->oclass         = SURFACE;
    obj->mtype          = BSPLINE;
    egadsSurface *psurf = new egadsSurface;
    psurf->handle       = hSurface;
    psurf->basis        = NULL;
    psurf->topFlg       = 0;
    obj->blind          = psurf;
    EG_referenceObject(obj, context);
    *newcurve = obj;

    return EGADS_SUCCESS;
  }

  // normal isocline code
  egadsSurface *psurf           = (egadsSurface *) surface->blind;
  Handle(Geom_Surface) hSurface = psurf->handle;
  Handle(Geom_Curve)   newcrv;
  if (UV == UISO) {
    newcrv = hSurface->UIso(value);
  } else {
    newcrv = hSurface->VIso(value);
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 1)
      printf(" EGADS Error: make Curve = %d (EG_otherCurve)!\n", stat);
    return stat;
  }
  EG_completeCurve(obj, newcrv);
  EG_referenceObject(obj, context);
  *newcurve = obj;

  return EGADS_SUCCESS;
}


int
EG_convertToBSplineRange(egObject *object, const double *range,
                         egObject **bspline)
{
  int      outLevel, stat, header[4];
  double   data[8];
  gp_Pnt2d pnt;
  egObject *obj, *geom, *context;
  
  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != PCURVE) &&
      (object->oclass != CURVE)  && (object->oclass != SURFACE) &&
      (object->oclass != EDGE)   && (object->oclass != FACE))
    return EGADS_NOTGEOM;
  if (object->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  geom     = object;
  
  if (object->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) object->blind;
    geom = pedge->curve;
    if (geom->blind == NULL) return EGADS_NODATA;
  }
  if (object->oclass == FACE) {
    egadsFace *pface = (egadsFace *) object->blind;
    geom = pface->surface;
    if (geom->blind == NULL) return EGADS_NODATA;
  }
  if (geom->mtype == BSPLINE) {
    *bspline = geom;
    return EGADS_SUCCESS;
  }
  
  if (geom->oclass == PCURVE) {
    
    egadsPCurve *ppcurv         = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;

    if (geom->mtype == LINE) {
      header[0] = 0;
      header[1] = 1;
      header[2] = 2;
      header[3] = 4;
      data[0]   = data[1] = range[0];
      data[2]   = data[3] = range[1];
      hCurve->D0(range[0], pnt);
      data[4]   = pnt.X();
      data[5]   = pnt.Y();
      hCurve->D0(range[1], pnt);
      data[6]   = pnt.X();
      data[7]   = pnt.Y();
      return EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                             bspline);
    }

    try {
      Handle(Geom2d_BSplineCurve)
        hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                         range[1],
                                                         Precision::Confusion(),
                                                         GeomAbs_C2, 100, 20);
      if (hBSpline.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to convert (EG_convertToBSplineRange)!\n");
        return EGADS_GEOMERR;
      }
    
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: make PCurve = %d (EG_convertToBSplineRange)!\n",
               stat);
        return stat;
      }
      obj->oclass        = PCURVE;
      obj->mtype         = BSPLINE;
      egadsPCurve *ppcrv = new egadsPCurve;
      ppcrv->handle      = hBSpline;
      ppcrv->basis       = NULL;
      ppcrv->topFlg      = 0;
      obj->blind         = ppcrv;
    }
    catch (Standard_Failure)
    {
      if (outLevel > 0) {
        printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
        Handle_Standard_Failure e = Standard_Failure::Caught();
        printf("                %s\n", e->GetMessageString());
      }
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      if (outLevel > 0)
        printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
      return EGADS_GEOMERR;
    }
    
  } else if (geom->oclass == CURVE) {
    
    egadsCurve *pcurve        = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    Handle(Geom_BSplineCurve)
      hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                       range[1],
                                                       Precision::Confusion(),
                                                       GeomAbs_C2, 100, 20);
    if (hBSpline.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to convert (EG_convertToBSplineRange)!\n");
      return EGADS_GEOMERR;
    }
    
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_convertToBSplineRange)!\n",
             stat);
      return stat;
    }
    obj->oclass       = CURVE;
    obj->mtype        = BSPLINE;
    egadsCurve *pcurv = new egadsCurve;
    pcurv->handle     = hBSpline;
    pcurv->basis      = NULL;
    pcurv->topFlg     = 0;
    obj->blind        = pcurv;
    
  } else {
    
    egadsSurface *psurface        = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurface = psurface->handle;
    try {
      Handle(Geom_BSplineSurface)
        hBSpline = ShapeConstruct::ConvertSurfaceToBSpline(hSurface,
                                   range[0], range[1], range[2], range[3],
                                   Precision::Confusion(), GeomAbs_C2, 100, 20);
      if (hBSpline.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to Convert (EG_convertToBSplineRange)!\n");
        return EGADS_GEOMERR;
      }
      
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: make Surface = %d (EG_convertToBSplineRange)!\n",
               stat);
        return stat;
      }
      obj->oclass         = SURFACE;
      obj->mtype          = BSPLINE;
      egadsSurface *psurf = new egadsSurface;
      psurf->handle       = hBSpline;
      psurf->basis        = NULL;
      psurf->topFlg       = 0;
      obj->blind          = psurf;
    }
    catch (Standard_Failure)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
      Handle_Standard_Failure e = Standard_Failure::Caught();
      printf("                %s\n", e->GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSplineRange)!\n");
      return EGADS_GEOMERR;
    }
    
  }
  
  *bspline = obj;
  EG_referenceObject(obj, context);
  
  return EGADS_SUCCESS;
}


int
EG_convertToBSpline(egObject *object, egObject **bspline)
{
  int           outLevel, stat;
  egObject      *obj, *geom, *context;
  Standard_Real range[4];

  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != PCURVE) &&
      (object->oclass != CURVE)  && (object->oclass != SURFACE) &&
      (object->oclass != EDGE)   && (object->oclass != FACE))
                                     return EGADS_NOTGEOM;
  if (object->blind == NULL)         return EGADS_NODATA;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);
  geom     = object;
  
  if (object->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) object->blind;
    geom = pedge->curve;
    if (geom->blind == NULL) return EGADS_NODATA;
  }
  if (object->oclass == FACE) {
    egadsFace *pface = (egadsFace *) object->blind;
    geom = pface->surface;
    if (geom->blind == NULL) return EGADS_NODATA;
  }
  if (geom->mtype == BSPLINE) {
    *bspline = geom;
    return EGADS_SUCCESS;
  }
  
  if (geom->oclass == PCURVE) {
  
    egadsPCurve *ppcurv         = (egadsPCurve *) geom->blind;
    Handle(Geom2d_Curve) hCurve = ppcurv->handle;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
/*
    ShapeConstruct_Curve ShapeCC;
    Handle(Geom2d_BSplineCurve)
      hBSpline = ShapeCC.ConvertToBSpline(hCurve, range[0], range[1],
                                          Precision::Confusion());
*/
    Handle(Geom2d_BSplineCurve)
      hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                       range[1],
                                                       Precision::Confusion(),
                                                       GeomAbs_C2, 100, 20);
    if (hBSpline.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to convert (EG_convertToBSpline)!\n");
      return EGADS_GEOMERR;
    }
    
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make PCurve = %d (EG_convertToBSpline)!\n", 
             stat);
      return stat;
    }
    obj->oclass        = PCURVE;
    obj->mtype         = BSPLINE;
    egadsPCurve *ppcrv = new egadsPCurve;
    ppcrv->handle      = hBSpline;
    ppcrv->basis       = NULL;
    ppcrv->topFlg      = 0;
    obj->blind         = ppcrv;

  } else if (geom->oclass == CURVE) {
  
    egadsCurve *pcurve        = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    range[0] = hCurve->FirstParameter();
    range[1] = hCurve->LastParameter();
    if (object->oclass == EDGE) {
      egadsEdge *pedge = (egadsEdge *) object->blind;
      BRep_Tool::Range(pedge->edge, range[0], range[1]);
    }
/*
    ShapeConstruct_Curve ShapeCC;
    Handle(Geom_BSplineCurve)
      hBSpline = ShapeCC.ConvertToBSpline(hCurve, range[0], range[1],
                                          Precision::Confusion());
*/
    Handle(Geom_BSplineCurve)
      hBSpline = ShapeConstruct::ConvertCurveToBSpline(hCurve, range[0],
                                                       range[1],
                                                       Precision::Confusion(),
                                                       GeomAbs_C2, 100, 20);
    if (hBSpline.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Warning: Failure to convert (EG_convertToBSpline)!\n");
      return EGADS_GEOMERR;
    }
    
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: make Curve = %d (EG_convertToBSpline)!\n", 
             stat);
      return stat;
    }
    obj->oclass       = CURVE;
    obj->mtype        = BSPLINE;
    egadsCurve *pcurv = new egadsCurve;
    pcurv->handle     = hBSpline;
    pcurv->basis      = NULL;
    pcurv->topFlg     = 0;
    obj->blind        = pcurv;

  } else {
  
    egadsSurface *psurface        = (egadsSurface *) geom->blind;
    Handle(Geom_Surface) hSurface = psurface->handle;
    hSurface->Bounds(range[0],range[1], range[2],range[3]);
    if (object->oclass == FACE) {
      egadsFace *pface = (egadsFace *) object->blind;
      BRepTools::UVBounds(pface->face, range[0],range[1], range[2],range[3]);
    }
    try {
      Handle(Geom_BSplineSurface)
        hBSpline = ShapeConstruct::ConvertSurfaceToBSpline(hSurface,
                                   range[0], range[1], range[2], range[3],
                                   Precision::Confusion(), GeomAbs_C2, 100, 20);
      if (hBSpline.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Warning: Failure to Convert (EG_convertToBSpline)!\n");
        return EGADS_GEOMERR;
      }
    
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        printf(" EGADS Error: make Surface = %d (EG_convertToBSpline)!\n",
               stat);
        return stat;
      }
      obj->oclass         = SURFACE;
      obj->mtype          = BSPLINE;
      egadsSurface *psurf = new egadsSurface;
      psurf->handle       = hBSpline;
      psurf->basis        = NULL;
      psurf->topFlg       = 0;
      obj->blind          = psurf;
    }
    catch (Standard_Failure)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSpline)!\n");
      Handle_Standard_Failure e = Standard_Failure::Caught();
      printf("                %s\n", e->GetMessageString());
      return EGADS_GEOMERR;
    }
    catch (...)
    {
      printf(" EGADS Warning: Geometry Creation Error (EG_convertToBSpline)!\n");
      return EGADS_GEOMERR;
    }

  }
  
  *bspline = obj;
  EG_referenceObject(obj, context);
    
  return EGADS_SUCCESS;
}


#ifdef MAPBSPLINE
static int
EG_mapBSplines(egObject *src, egObject *dst, egObject **result)
{
  int      i, j, k, outLevel, stat, imax, jmax, oclass, mtype, per, *splInfo;
  double   data[18], crosT[12], param[2], *xyzs, *spl, *tang, *uk, *vk, dran[4];
  egObject *context, *ref;
  
  outLevel = EG_outLevel(dst);
  context  = EG_context(dst);
  
  dran[0] = dran[1] = dran[2] = dran[3] = 0.0;
  stat = EG_getRange(dst, dran, &per);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getRange = %d for dst (EG_mapKnots)!\n", stat);
    return stat;
  }
  stat = EG_getGeometry(src, &oclass, &mtype, &ref, &splInfo, &spl);
  if (stat != EGADS_SUCCESS) return stat;
  if ((splInfo[1] != 3) || (splInfo[4] != 3)) {
    if (outLevel > 0)
      printf(" EGADS Error: uDeg = %d vDeg = %d (EG_mapKnots)!\n", splInfo[1],
             splInfo[4]);
    EG_free(splInfo);
    EG_free(spl);
    return EGADS_CONSTERR;
  }
  imax = splInfo[2] - 2;
  jmax = splInfo[5] - 2;
  uk   = &spl[3];
  vk   = &spl[splInfo[3]+3];
  tang = (double *) EG_alloc(6*(imax+jmax)*sizeof(double));
  xyzs = (double *) EG_alloc(3*imax*jmax*sizeof(double));
  if ((xyzs == NULL) || (tang == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: CP/tan malloc %d %d (EG_mapKnots)!\n", imax, jmax);
    if (xyzs != NULL) EG_free(xyzs);
    if (tang != NULL) EG_free(tang);
    EG_free(splInfo);
    EG_free(spl);
    return EGADS_MALLOC;
  }
  
  /* fill the coordinate data based on the scaled source knots */
  for (k = j = 0; j < jmax; j++) {
    param[1] = dran[2] + (vk[j]-vk[0])*(dran[3]-dran[2])/(vk[jmax-1]-vk[0]);
    for (i = 0; i < imax; i++, k++) {
      param[0] = dran[0] + (uk[i]-uk[0])*(dran[1]-dran[0])/(uk[imax-1]-uk[0]);
      stat = EG_evaluate(dst, param, data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: point [%d,%d] evaluate = %d (EG_mapKnots)!\n",
                 i+1, j+1, stat);
        EG_free(xyzs);
        EG_free(tang);
        EG_free(splInfo);
        EG_free(spl);
        return stat;
      }
      xyzs[3*k  ] = data[0];
      xyzs[3*k+1] = data[1];
      xyzs[3*k+2] = data[2];
    }
  }
  
  /* fill in the tangent values */
  param[1] = dran[2];
  for (k = i = 0; i < imax; i++, k++) {
    param[0] = dran[0] + (uk[i]-uk[0])*(dran[1]-dran[0])/(uk[imax-1]-uk[0]);
    data[6]  = data[7] = data[8] = 0.0;
    stat = EG_evaluate(dst, param, data);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: south [%d] evaluate = %d (EG_mapKnots)!\n",
               i+1, stat);
      EG_free(xyzs);
      EG_free(tang);
      EG_free(splInfo);
      EG_free(spl);
      return stat;
    }
    tang[3*k  ] = data[6]*(dran[3]-dran[2]);
    tang[3*k+1] = data[7]*(dran[3]-dran[2]);
    tang[3*k+2] = data[8]*(dran[3]-dran[2]);
  }
  param[1] = dran[3];
  for (i = 0; i < imax; i++, k++) {
    param[0] = dran[0] + (uk[i]-uk[0])*(dran[1]-dran[0])/(uk[imax-1]-uk[0]);
    stat = EG_evaluate(dst, param, data);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: north [%d] evaluate = %d (EG_mapKnots)!\n",
               i+1, stat);
      EG_free(xyzs);
      EG_free(tang);
      EG_free(splInfo);
      EG_free(spl);
      return stat;
    }
    tang[3*k  ] = data[6]*(dran[3]-dran[2]);
    tang[3*k+1] = data[7]*(dran[3]-dran[2]);
    tang[3*k+2] = data[8]*(dran[3]-dran[2]);
  }
  param[0] = dran[0];
  for (j = 0; j < jmax; j++, k++) {
    param[1] = dran[2] + (vk[j]-vk[0])*(dran[3]-dran[2])/(vk[jmax-1]-vk[0]);
    stat = EG_evaluate(dst, param, data);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: west [%d] evaluate = %d (EG_mapKnots)!\n",
               j+1, stat);
      EG_free(xyzs);
      EG_free(tang);
      EG_free(splInfo);
      EG_free(spl);
      return stat;
    }
    tang[3*k  ] = data[3]*(dran[1]-dran[0]);
    tang[3*k+1] = data[4]*(dran[1]-dran[0]);
    tang[3*k+2] = data[5]*(dran[1]-dran[0]);
  }
  param[0] = dran[1];
  for (j = 0; j < jmax; j++, k++) {
    param[1] = dran[2] + (vk[j]-vk[0])*(dran[3]-dran[2])/(vk[jmax-1]-vk[0]);
    stat = EG_evaluate(dst, param, data);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: east [%d] evaluate = %d (EG_mapKnots)!\n",
               j+1, stat);
      EG_free(xyzs);
      EG_free(tang);
      EG_free(splInfo);
      EG_free(spl);
      return stat;
    }
    tang[3*k  ] = data[3]*(dran[1]-dran[0]);
    tang[3*k+1] = data[4]*(dran[1]-dran[0]);
    tang[3*k+2] = data[5]*(dran[1]-dran[0]);
  }

  /* the cross terms */
  param[0] = dran[0];
  param[1] = dran[2];
  data[12] = data[13] = data[14] = 0.0;
  stat = EG_evaluate(dst, param, data);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: point D evaluate = %d (EG_mapKnots)!\n", stat);
    EG_free(xyzs);
    EG_free(tang);
    EG_free(splInfo);
    EG_free(spl);
    return stat;
  }
  crosT[0] = data[12]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[1] = data[13]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[2] = data[14]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  
  param[0] = dran[1];
  stat = EG_evaluate(dst, param, data);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: point F evaluate = %d (EG_mapKnots)!\n", stat);
    EG_free(xyzs);
    EG_free(tang);
    EG_free(splInfo);
    EG_free(spl);
    return stat;
  }
  crosT[3] = data[12]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[4] = data[13]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[5] = data[14]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  
  param[1] = dran[3];
  stat = EG_evaluate(dst, param, data);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: point M evaluate = %d (EG_mapKnots)!\n", stat);
    EG_free(xyzs);
    EG_free(tang);
    EG_free(splInfo);
    EG_free(spl);
    return stat;
  }
  crosT[6] = data[12]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[7] = data[13]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[8] = data[14]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  
  param[0] = dran[0];
  stat = EG_evaluate(dst, param, data);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: point K evaluate = %d (EG_mapKnots)!\n", stat);
    EG_free(xyzs);
    EG_free(tang);
    EG_free(splInfo);
    EG_free(spl);
    return stat;
  }
  crosT[ 9] = data[12]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[10] = data[13]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  crosT[11] = data[14]*(dran[1]-dran[0])*(dran[3]-dran[2]);
  
  /* refit our surface using the source knots */
  stat = EG_spline2dFit(context, crosT, imax, uk, &tang[0], &tang[3*imax],
                        jmax, vk, &tang[6*imax], &tang[6*imax+3*jmax], xyzs,
                        1.e-8, result);
  if (outLevel > 0)
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Error: EG_spline2dFit = %d (EG_mapKnots)!\n", stat);
  EG_free(xyzs);
  EG_free(tang);
  EG_free(splInfo);
  EG_free(spl);

  return stat;
}


int
EG_mapKnots(egObject *src, egObject *dst, egObject **result)
{
  int      i, outLevel, nent, stat;
  egObject *context, *obj;

  *result = NULL;
  if  (src == NULL)                return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if ((src->oclass != SURFACE) &&
      (src->oclass != FACE))       return EGADS_NOTGEOM;
  if  (src->blind == NULL)         return EGADS_NODATA;
  if  (dst == NULL)                return EGADS_NULLOBJ;
  if  (dst->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (dst->oclass != src->oclass) return EGADS_NOTGEOM;
  if  (dst->mtype  != src->mtype)  return EGADS_GEOMERR;
  if  (dst->blind == NULL)         return EGADS_NODATA;

  if (src->oclass == SURFACE) {
    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    return EG_mapBSplines(src, dst, result);
  }
  
  outLevel = EG_outLevel(dst);
  context  = EG_context(dst);
  BRep_Builder    Builder;
  TopLoc_Location L;

  egadsFace *pface   = (egadsFace *) src->blind;
  egObject  *surfs   = (egObject *)  pface->surface;
  egadsFace *pfacd   = (egadsFace *) dst->blind;
  egObject  *surfd   = (egObject *)  pfacd->surface;
  if (surfs->mtype  != BSPLINE) return EGADS_GEOMERR;
  if (surfd->mtype  != BSPLINE) return EGADS_GEOMERR;
    
  gp_Trsf       form = gp_Trsf();
  TopoDS_Shape shape = pfacd->face;
  BRepBuilderAPI_Transform xForm(shape, form, Standard_True);
  if (!xForm.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Can't copy Face Topology (EG_mapKnots)!\n");
    return EGADS_CONSTERR;
  }
  TopoDS_Shape nTopo = xForm.ModifiedShape(shape);
  TopoDS_Face  Face  = TopoDS::Face(nTopo);

  egObject *surfr = NULL;
  stat = EG_mapBSplines(surfs, surfd, &surfr);
  if ((stat != EGADS_SUCCESS) || (surfr == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_mapBSplines = %d (EG_mapKnots)!\n", stat);
    return stat;
  }
  if (surfr->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_mapBSplines -> NULL blind (EG_mapKnots)!\n");
    return EGADS_NULLOBJ;
  }
  egadsSurface *psurface      = (egadsSurface *) surfr->blind;
  Handle(Geom_Surface) hSurfr = psurface->handle;
    
  // update the Face
  Builder.UpdateFace(Face, hSurfr, L, BRep_Tool::Tolerance(pfacd->face));
  BRepCheck_Analyzer fCheck(Face);
  if (!fCheck.IsValid()) {
    // try to fix the fault
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(Face);
    sfs->Perform();
    TopoDS_Shape fixedFace = sfs->Shape();
    if (fixedFace.IsNull()) {
      if (outLevel > 0) {
        printf(" EGADS Info: Invalid Face w/ NULL Fix (EG_mapKnots)!\n");
        EG_checkStatus(fCheck.Result(Face));
      }
      EG_deleteObject(surfr);
      return EGADS_CONSTERR;
    }
    BRepCheck_Analyzer fxCheck(fixedFace);
    if (!fxCheck.IsValid()) {
      if (outLevel > 0) {
        printf(" EGADS Info: Face is invalid (EG_mapKnots)!\n");
        EG_checkStatus(fxCheck.Result(fixedFace));
      }
      EG_deleteObject(surfr);
      return EGADS_CONSTERR;
    }
    Face = TopoDS::Face(fixedFace);
  }
    
  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Face makeObject = %d (EG_mapKnots)!\n", stat);
    EG_deleteObject(surfr);
    return EGADS_CONSTERR;
  }
  
  egadsBody       ebody;
  TopExp_Explorer Exp;
  TopExp::MapShapes(Face, TopAbs_VERTEX, ebody.nodes.map);
  nent = ebody.nodes.map.Extent();
  ebody.nodes.objs = new egObject*[nent];
  for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
  TopExp::MapShapes(Face, TopAbs_EDGE,   ebody.edges.map);
  nent = ebody.edges.map.Extent();
  ebody.edges.objs = new egObject*[nent];
  for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;
  TopExp::MapShapes(Face, TopAbs_WIRE,   ebody.loops.map);
  nent = ebody.loops.map.Extent();
  ebody.loops.objs = new egObject*[nent];
  for (i = 0; i < nent; i++) ebody.loops.objs[i] = NULL;
    
  egadsFace *pfacn = new egadsFace;
  pfacn->face      = Face;
  obj->blind       = pfacn;
  EG_copyAttrTopo(&ebody, NULL, form, dst, obj, obj);
    
  delete [] ebody.loops.objs;
  delete [] ebody.edges.objs;
  delete [] ebody.nodes.objs;
  EG_referenceObject(obj, context);
  EG_deleteObject(surfr);

  *result = obj;
  return EGADS_SUCCESS;
}
#endif


int
EG_mapSequen(egObject *src, egObject *dst, egObject **result)
{
  int      i, j, hit, outLevel, stat;
  egObject *context, *obj;
  
  if  (src == NULL)                return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if ((src->oclass != PCURVE)  && (src->oclass != CURVE) &&
      (src->oclass != SURFACE))    return EGADS_NOTGEOM;
  if  (src->blind == NULL)         return EGADS_NODATA;
  if  (dst == NULL)                return EGADS_NULLOBJ;
  if  (dst->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (dst->oclass != src->oclass) return EGADS_NOTGEOM;
  if  (dst->mtype  != src->mtype)  return EGADS_GEOMERR;
  if  (dst->blind == NULL)         return EGADS_NODATA;

  *result  = NULL;
  outLevel = EG_outLevel(dst);
  context  = EG_context(dst);

  if (src->oclass == PCURVE) {
    
    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    egadsPCurve *ppcurv            = (egadsPCurve *) src->blind;
    Handle(Geom2d_Curve)    hCurve = ppcurv->handle;
    Handle(Geom2d_BSplineCurve)
                          hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
    egadsPCurve *ppcurvd           = (egadsPCurve *) dst->blind;
    Handle(Geom2d_Curve)    hCurvt = ppcurvd->handle;
    gp_Trsf2d                 form = gp_Trsf2d();
    Handle(Geom2d_Geometry)  nGeom = hCurvt->Transformed(form);
    Handle(Geom2d_Curve)    hCurvd = Handle(Geom2d_Curve)::DownCast(nGeom);
    Handle(Geom2d_BSplineCurve)
                          hBSplind = Handle(Geom2d_BSplineCurve)::DownCast(hCurvd);
    if (hBSpline->NbKnots() != hBSplind->NbKnots()) return EGADS_INDEXERR;
    int len = hBSpline->NbKnots();
    TColStd_Array1OfInteger mults(1, len);
    TColStd_Array1OfInteger multd(1, len);
    hBSpline->Multiplicities(mults);
    hBSplind->Multiplicities(multd);
    for (i = 1; i <= len; i++)
      if (mults(i) != multd(i)) return EGADS_RANGERR;

    TColStd_Array1OfReal knots(1, len);
    TColStd_Array1OfReal knotd(1, len);
    hBSpline->Knots(knots);
    hBSplind->Knots(knotd);
    hit = 0;
    for (i = 2; i < len; i++) {
      double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
      scaledKnot *= knotd(len)-knotd(1);
      scaledKnot += knotd(1);
      for (j = 2; j < len; j++)
        if (fabs(scaledKnot-knotd(j)) <= 0.5*KNDIFF) break;
      if (j != len) continue;
      hBSplind->InsertKnot(scaledKnot);
      hit++;
      if (outLevel > 1)
        printf("   inserting knot = %lf (%lf)\n", scaledKnot, knots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->Knot(j);
        for (i = 2; i < len; i++) {
          double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
          scaledKnot *= knotd(len)-knotd(1);
          scaledKnot += knotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != len) continue;
        if (outLevel > 1)
          printf("   removing  knot = %lf\n", hBSplind->Knot(j));
        if (!hBSplind->RemoveKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove knot %lf (EG_mapSequen)!",
                   hBSplind->Knot(j));
      }
    }
#endif
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: PCurve makeObject = %d (EG_mapSequen)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completePCurve(obj, hCurvd);
    
  } else if (src->oclass == CURVE) {
    
    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    egadsCurve *pcurve           = (egadsCurve *) src->blind;
    Handle(Geom_Curve)    hCurve = pcurve->handle;
    Handle(Geom_BSplineCurve)
                        hBSpline = Handle(Geom_BSplineCurve)::DownCast(hCurve);
    egadsCurve *pcurvd           = (egadsCurve *) dst->blind;
    Handle(Geom_Curve)    hCurvt = pcurvd->handle;
    gp_Trsf                 form = gp_Trsf();
    Handle(Geom_Geometry)  nGeom = hCurvt->Transformed(form);
    Handle(Geom_Curve)    hCurvd = Handle(Geom_Curve)::DownCast(nGeom);
    Handle(Geom_BSplineCurve)
                        hBSplind = Handle(Geom_BSplineCurve)::DownCast(hCurvd);
    if (hBSpline->NbKnots() != hBSplind->NbKnots()) return EGADS_INDEXERR;
    int len = hBSpline->NbKnots();
    TColStd_Array1OfInteger mults(1, len);
    TColStd_Array1OfInteger multd(1, len);
    hBSpline->Multiplicities(mults);
    hBSplind->Multiplicities(multd);
    for (i = 1; i <= len; i++)
      if (mults(i) != multd(i)) return EGADS_RANGERR;
/*  GeomAbs_BSplKnotDistribution kDist = hBSplind->KnotDistribution();
    if (kDist == GeomAbs_NonUniform)      printf(" NonUniform!\n");
    if (kDist == GeomAbs_Uniform)         printf(" Uniform!\n");
    if (kDist == GeomAbs_QuasiUniform)    printf(" QuasiUniform!\n");
    if (kDist == GeomAbs_PiecewiseBezier) printf(" PiecewiseBezier!\n");  */
    
    TColStd_Array1OfReal knots(1, len);
    TColStd_Array1OfReal knotd(1, len);
    hBSpline->Knots(knots);
    hBSplind->Knots(knotd);
    hit = 0;
    for (i = 2; i < len; i++) {
      double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
      scaledKnot *= knotd(len)-knotd(1);
      scaledKnot += knotd(1);
      for (j = 2; j < len; j++)
        if (fabs(scaledKnot-knotd(j)) <= 0.5*KNDIFF) break;
      if (j != len) continue;
      hBSplind->InsertKnot(scaledKnot);
      hit++;
      if (outLevel > 1)
        printf("   inserting knot = %lf (%lf)\n", scaledKnot, knots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->Knot(j);
        for (i = 2; i < len; i++) {
          double scaledKnot = (knots(i)-knots(1))/(knots(len)-knots(1));
          scaledKnot *= knotd(len)-knotd(1);
          scaledKnot += knotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != len) continue;
        if (outLevel > 1)
          printf("   removing  knot = %lf\n", hBSplind->Knot(j));
        if (!hBSplind->RemoveKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove knot %lf (EG_mapSequen)!",
                   hBSplind->Knot(j));
      }
    }
#endif
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve makeObject = %d (EG_mapSequen)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeCurve(obj, hCurvd);
    
  } else {
    
    if (src->mtype != BSPLINE) return EGADS_GEOMERR;
    egadsSurface *psurface       = (egadsSurface *) src->blind;
    Handle(Geom_Surface)  hSurf  = psurface->handle;
    Handle(Geom_BSplineSurface)
                        hBSpline = Handle(Geom_BSplineSurface)::DownCast(hSurf);
    egadsSurface *psurfacd       = (egadsSurface *) dst->blind;
    Handle(Geom_Surface)  hSurft = psurfacd->handle;
    gp_Trsf                 form = gp_Trsf();
    Handle(Geom_Geometry)  nGeom = hSurft->Transformed(form);
    Handle(Geom_Surface)  hSurfd = Handle(Geom_Surface)::DownCast(nGeom);
    Handle(Geom_BSplineSurface)
                        hBSplind = Handle(Geom_BSplineSurface)::DownCast(hSurfd);
    if (hBSpline->NbUKnots() != hBSplind->NbUKnots()) return EGADS_INDEXERR;
    if (hBSpline->NbVKnots() != hBSplind->NbVKnots()) return EGADS_INDEXERR;
    int uLen = hBSpline->NbUKnots();
    int vLen = hBSpline->NbVKnots();
    TColStd_Array1OfInteger uMults(1, uLen);
    TColStd_Array1OfInteger uMultd(1, uLen);
    TColStd_Array1OfInteger vMults(1, vLen);
    TColStd_Array1OfInteger vMultd(1, vLen);
    hBSpline->UMultiplicities(uMults);
    hBSplind->UMultiplicities(uMultd);
    for (i = 1; i <= uLen; i++)
      if (uMults(i) != uMultd(i)) return EGADS_RANGERR;
    hBSpline->VMultiplicities(vMults);
    hBSplind->VMultiplicities(vMultd);
    for (i = 1; i <= vLen; i++)
      if (vMults(i) != vMultd(i)) return EGADS_RANGERR;
    
    TColStd_Array1OfReal uKnots(1, uLen);
    TColStd_Array1OfReal vKnots(1, vLen);
    TColStd_Array1OfReal uKnotd(1, uLen);
    TColStd_Array1OfReal vKnotd(1, vLen);
    hBSpline->UKnots(uKnots);
    hBSpline->VKnots(vKnots);
    hBSplind->UKnots(uKnotd);
    hBSplind->VKnots(vKnotd);

    // u knots
    hit = 0;
    for (i = 2; i < uLen; i++) {
      double scaledKnot = (uKnots(i)-uKnots(1))/(uKnots(uLen)-uKnots(1));
      scaledKnot *= uKnotd(uLen)-uKnotd(1);
      scaledKnot += uKnotd(1);
      for (j = 2; j < uLen; j++)
        if (fabs(scaledKnot-uKnotd(j)) <= 0.5*KNDIFF) break;
      if (j != uLen) continue;
      hBSplind->InsertUKnot(scaledKnot, 1, 0.5*KNDIFF);
      hit++;
      if (outLevel > 1)
        printf("   inserting u knot = %lf (%lf)\n", scaledKnot, uKnots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbUKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->UKnot(j);
        for (i = 2; i < uLen; i++) {
          double scaledKnot = (uKnots(i)-uKnots(1))/(uKnots(uLen)-uKnots(1));
          scaledKnot *= uKnotd(uLen)-uKnotd(1);
          scaledKnot += uKnotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != uLen) continue;
        if (outLevel > 1)
          printf("   removing  u knot = %lf\n", hBSplind->UKnot(j));
        if (!hBSplind->RemoveUKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove U knot %lf (EG_mapSequen)!",
                   hBSplind->UKnot(j));
      }
    }
#endif
    
    // v knots
    hit = 0;
    for (i = 2; i < vLen; i++) {
      double scaledKnot = (vKnots(i)-vKnots(1))/(vKnots(vLen)-vKnots(1));
      scaledKnot *= vKnotd(vLen)-vKnotd(1);
      scaledKnot += vKnotd(1);
      for (j = 2; j < vLen; j++)
        if (fabs(scaledKnot-vKnotd(j)) <= 0.5*KNDIFF) break;
      if (j != vLen) continue;
      hBSplind->InsertVKnot(scaledKnot, 1, 0.5*KNDIFF);
      hit++;
      if (outLevel > 1)
        printf("   inserting v knot = %lf (%lf)\n", scaledKnot, vKnots(i));
    }
#ifdef KNOTREMOVE
    if (hit != 0) {
      int nknot = hBSplind->NbVKnots();
      for (j = nknot-1; j >= 2; j--) {
        double knot = hBSplind->VKnot(j);
        for (i = 2; i < vLen; i++) {
          double scaledKnot = (vKnots(i)-vKnots(1))/(vKnots(vLen)-vKnots(1));
          scaledKnot *= vKnotd(vLen)-vKnotd(1);
          scaledKnot += vKnotd(1);
          if (fabs(scaledKnot-knot) <= 0.5*KNDIFF) break;
        }
        if (i != vLen) continue;
        if (outLevel > 1)
          printf("   removing  v knot = %lf\n", hBSplind->VKnot(j));
        if (!hBSplind->RemoveVKnot(j, 0, 1.0))
          if (outLevel > 0)
            printf(" EGADS Warning: Cannot remove V knot %lf (EG_mapSequen)!",
                   hBSplind->VKnot(j));
      }
    }
#endif
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface makeObject = %d (EG_mapSequen)!\n", stat);
      return EGADS_CONSTERR;
    }
    EG_completeSurf(obj, hSurfd);
    
  }

  *result = obj;
  return EGADS_SUCCESS;
}


void
EG_mapTessTs(egTess1D src, egTess1D dst)
{
  if (src.npts != dst.npts) {
    printf(" EGADS Warning: Len Mismatch src = %d, dst = %d (EG_mapTessTs)!\n",
           src.npts, dst.npts);
    return;
  }
  egadsEdge *pedgs = (egadsEdge *) src.obj->blind;
  if (pedgs == NULL) {
    printf(" EGADS Warning: NULL src Edge Object (EG_mapTessTs)!\n");
    return;
  }
  egObject  *curvs = pedgs->curve;
  if (curvs == NULL) {
    printf(" EGADS Warning: No curve Object for src Edge (EG_mapTessTs)!\n");
    return;
  }
  egadsCurve *pcurvs = (egadsCurve *) curvs->blind;
  if (pcurvs == NULL) {
    printf(" EGADS Warning: No curve Data for src Edge (EG_mapTessTs)!\n");
    return;
  }
  egadsEdge *pedgd = (egadsEdge *) dst.obj->blind;
  if (pedgd == NULL) {
    printf(" EGADS Warning: NULL dst Edge Object (EG_mapTessTs)!\n");
    return;
  }
  egObject  *curvd = pedgd->curve;
  if (curvd == NULL) {
    printf(" EGADS Warning: No curve Object for dst Edge (EG_mapTessTs)!\n");
    return;
  }
  egadsCurve *pcurvd = (egadsCurve *) curvd->blind;
  if (pcurvd == NULL) {
    printf(" EGADS Warning: No curve Data for dst Edge (EG_mapTessTs)!\n");
    return;
  }
  
  int n = src.npts;
  GeomAdaptor_Curve ACsrc(pcurvs->handle);
  GeomAdaptor_Curve ACdst(pcurvd->handle);
  double slen = GCPnts_AbscissaPoint::Length(ACsrc, src.t[0], src.t[n-1]);
  double dlen = GCPnts_AbscissaPoint::Length(ACdst, dst.t[0], dst.t[n-1]);
  
  //
  // have the relative arcLengths in the destination match the source
  for (int i = 1; i < n-1; i++) {
    double srcAlen = GCPnts_AbscissaPoint::Length(ACsrc, src.t[0], src.t[i]);
    double tgtAlen = dlen*srcAlen/slen;
    GCPnts_AbscissaPoint AP(ACdst, tgtAlen, dst.t[0], dst.t[i]);
    if (!AP.IsDone()) continue;
    double t       = AP.Parameter();
/*  printf("    %d:  %lf %lf\n", i+1, dst.t[i], t);  */
    gp_Pnt P0;
    pcurvd->handle->D0(t, P0);
    dst.t[i]       = t;
    dst.xyz[3*i  ] = P0.X();
    dst.xyz[3*i+1] = P0.Y();
    dst.xyz[3*i+2] = P0.Z();
  }
}


int
EG_relPosTs(egObject *geom, int n, const double *rel, double *ts, double *xyzs)
{
  egadsCurve *pcurve;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != CURVE)  &&
      (geom->oclass != EDGE))      return EGADS_NOTGEOM;
  if (geom->blind == NULL)         return EGADS_NODATA;
  
  if (geom->oclass == CURVE) {
    
    pcurve = (egadsCurve *) geom->blind;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    
  } else {
    
    if (geom->mtype == DEGENERATE) return EGADS_DEGEN;
    egadsEdge *pedge  = (egadsEdge *) geom->blind;
    if (pedge  == NULL) return EGADS_NULLOBJ;
    egObject  *curvo  = pedge->curve;
    if (curvo  == NULL) return EGADS_NULLOBJ;
    pcurve = (egadsCurve *) curvo->blind;
    
  }
  if (pcurve == NULL) return EGADS_NULLOBJ;
  GeomAdaptor_Curve AC(pcurve->handle);

  double alen = GCPnts_AbscissaPoint::Length(AC, ts[0], ts[n-1]);
  for (int i = 1; i < n-1; i++) {
    GCPnts_AbscissaPoint AP(AC, rel[i-1]*alen, ts[0]);
    if (!AP.IsDone()) continue;
    double t    = AP.Parameter();
/*  printf("    %d:  %lf %lf\n", i, ts[i], t);  */
    gp_Pnt P0;
    pcurve->handle->D0(t, P0);
    ts[i]       = t;
    xyzs[3*i  ] = P0.X();
    xyzs[3*i+1] = P0.Y();
    xyzs[3*i+2] = P0.Z();
  }
  
  return EGADS_SUCCESS;
}
