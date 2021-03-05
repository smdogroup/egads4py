/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Topology Functions
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <setjmp.h>


#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"

//#define MAPBSPLINE    /* also in egadsGeom.cpp */


/* structures */

  typedef struct {
    int    nFace;		/* the number of Faces in the Edge */
    int    seq;                 /* the sequence number */
    int    *fIndices;           /* the Face indices (nFace in len) */
    char   *ID;                 /* the ID including seq # */
    double CG[4];               /* Center of Gravity, then Length */
  } edgeID;


  extern "C" int  EG_destroyTopology( egObject *topo );
  extern "C" int  EG_copyObject( const egObject *object, /*@null@*/ void *ptr,
                                       egObject **copy );
  extern "C" int  EG_isSame( const egObject *obj1, const egObject *obj2 );
  extern "C" int  EG_attributeRet( const egObject *obj, const char *name,
                                   int *atype, int *len,
                                   /*@null@*/ const int **ints,
                                   /*@null@*/ const double **reals,
                                   /*@null@*/ const char **str );
  extern "C" int  EG_attributeAdd( egObject *obj, const char *name, int atype,
                                   int len, /*@null@*/ const int    *ints,
                                            /*@null@*/ const double *reals,
                                            /*@null@*/ const char   *str );
  extern "C" int  EG_mapKnots( egObject *src, egObject *dst, egObject **result );

  extern "C" int  EG_getGeometry( const egObject *geom, int *oclass, int *mtype,
                                  egObject **refGeom, int **ivec, double **rvec );
  extern "C" int  EG_otherCurve( const egObject *surface, const egObject *curve,
                                 double tol, egObject **newcurve );

  extern "C" int  EG_getTolerance( const egObject *topo, double *tol );
  extern "C" int  EG_setTolerance( const egObject *topo, double  tol );
  extern "C" int  EG_getTopology( const egObject *topo, egObject **geom,
                                  int *oclass, int *type,
                                  /*@null@*/ double *limits, int *nChildren,
                                  egObject ***children, int **senses );
  extern "C" int  EG_makeTopology( egObject *context, /*@null@*/ egObject *geom,
                                   int oclass, int mtype, /*@null@*/ double *limits,
                                   int nChildren, /*@null@*/ egObject **children,
                                   /*@null@*/ int *senses, egObject **topo );
  extern "C" int  EG_makeLoop( int nedge, egObject **edges,
                               /*@null@*/ egObject *geom, double toler,
                               egObject **result );
  extern "C" int  EG_makeFace( egObject *object, int mtype,
                               /*@null@*/ const double *limits, egObject **face );
  extern "C" int  EG_getPlane( const egObject *object, egObject **plane );
  extern "C" int  EG_getArea( egObject *object, /*@null@*/ const double *limits,
                              double *area );
  extern "C" int  EG_getBodyTopos( const egObject *body, /*@null@*/ egObject *src,
                                   int oclass, int *ntopo, egObject ***topos );
  extern "C" int  EG_indexBodyTopo( const egObject *body, const egObject *src );
  extern "C" int  EG_sameBodyTopo( const egObject *bod1, const egObject *bod2 );
  extern "C" int  EG_makeSolidBody( egObject *context, int stype,
                                    const double *rvec, egObject **body );
  extern "C" int  EG_getBoundingBox( const egObject *topo, double *box );
  extern "C" int  EG_getMassProperties( const egObject *topo, double *props );
  extern "C" int  EG_isEquivalent( const egObject *topo1, const egObject *topo2 );
  extern "C" int  EG_isPlanar( const egObject *topo );
  extern "C" int  EG_getEdgeUV( const egObject *face, const egObject *edge,
                                int sense, double t, double *result );
  extern "C" int  EG_getEdgeUVeval( const egObject *face, const egObject *edge,
                                    int sense, double t, double *result );
  extern "C" int  EG_getBody( const egObject *topo, egObject **result );
  extern "C" int  EG_inTopology( const egObject *topo, const double *xyz );
  extern "C" int  EG_inFace( const egObject *face, const double *uv );
  extern "C" int  EG_sewFaces( int nobj, const egObject **objs, double toler,
                               int opt, egObject **result );
  extern "C" int  EG_replaceFaces( const egObject *body,  int nobj,
                                         egObject **objs, egObject **result );
  extern "C" int  EG_mapBody( const egObject *sBody, const egObject *dBody,
                              const char *fAttr, egObject **mBody );

  extern     int  EG_attriBodyDup( const egObject *src, egObject *dst );
  extern     void EG_completePCurve( egObject *g, Handle(Geom2d_Curve) &hCurv );
  extern     void EG_completeCurve(  egObject *g, Handle(Geom_Curve)   &hCurv );
  extern     void EG_completeSurf(   egObject *g, Handle(Geom_Surface) &hSurf );


// for trapping SegFaults
static jmp_buf jmpenv;

static void
segfault_handler(int x)
{
  longjmp(jmpenv, x);
}


static void
EG_cleanMaps(egadsMap *map)
{
  if (map->objs == NULL) return;
  EG_free(map->objs);
  map->objs = NULL;
}


static void
EG_printStatus(BRepCheck_Status cStatus)
{
  if (cStatus == BRepCheck_NoError) return;

  if (cStatus == BRepCheck_InvalidPointOnCurve) {
    printf(" EGADS Fault: Invalid Point On Curve\n");
  } else if (cStatus == BRepCheck_InvalidPointOnCurveOnSurface) {
    printf(" EGADS Fault: Invalid Point On Curve On Surface\n");
  } else if (cStatus == BRepCheck_InvalidPointOnSurface) {
    printf(" EGADS Fault: Invalid Point On Surface\n");
  } else if (cStatus == BRepCheck_No3DCurve) {
    printf(" EGADS Fault: No 3D Curve\n");
  } else if (cStatus == BRepCheck_Multiple3DCurve) {
    printf(" EGADS Fault: Multiple 3D Curves\n");
  } else if (cStatus == BRepCheck_Invalid3DCurve) {
    printf(" EGADS Fault: Invalid 3D Curve\n");
  } else if (cStatus == BRepCheck_NoCurveOnSurface) {
    printf(" EGADS Fault: No Curve On Surface\n");
  } else if (cStatus == BRepCheck_InvalidCurveOnSurface) {
    printf(" EGADS Fault: Invalid Curve On Surface\n");
  } else if (cStatus == BRepCheck_InvalidCurveOnClosedSurface) {
    printf(" EGADS Fault: Invalid Curve On Closed Surface\n");
  } else if (cStatus == BRepCheck_InvalidSameRangeFlag) {
    printf(" EGADS Fault: Invalid SameRange Flag\n");
  } else if (cStatus == BRepCheck_InvalidSameParameterFlag) {
    printf(" EGADS Fault: Invalid Same Parameter Flag\n");
  } else if (cStatus == BRepCheck_InvalidDegeneratedFlag) {
    printf(" EGADS Fault: Invalid Degenerated Flag\n");
  } else if (cStatus == BRepCheck_FreeEdge) {
    printf(" EGADS Fault: Free Edge\n");
  } else if (cStatus == BRepCheck_InvalidMultiConnexity) {
    printf(" EGADS Fault: Invalid Multi Connexity\n");
  } else if (cStatus == BRepCheck_InvalidRange) {
    printf(" EGADS Fault: Invalid Range\n");
  } else if (cStatus == BRepCheck_EmptyWire) {
    printf(" EGADS Fault: Empty Wire\n");
  } else if (cStatus == BRepCheck_RedundantEdge) {
    printf(" EGADS Fault: Redundant Edge\n");
  } else if (cStatus == BRepCheck_SelfIntersectingWire) {
    printf(" EGADS Fault: Self Intersecting Wire\n");
  } else if (cStatus == BRepCheck_NoSurface) {
    printf(" EGADS Fault: No Surface\n");
  } else if (cStatus == BRepCheck_InvalidWire) {
    printf(" EGADS Fault: Invalid Wire\n");
  } else if (cStatus == BRepCheck_RedundantWire) {
    printf(" EGADS Fault: Redundant Wire\n");
  } else if (cStatus == BRepCheck_IntersectingWires) {
    printf(" EGADS Fault: Intersecting Wires\n");
  } else if (cStatus == BRepCheck_InvalidImbricationOfWires) {
    printf(" EGADS Fault: Invalid Imbrication Of Wires\n");
  } else if (cStatus == BRepCheck_EmptyShell) {
    printf(" EGADS Fault: Empty Shell\n");
  } else if (cStatus == BRepCheck_RedundantFace) {
    printf(" EGADS Fault: Redundant Face\n");
  } else if (cStatus == BRepCheck_UnorientableShape) {
    printf(" EGADS Fault: Unorientable Shape\n");
  } else if (cStatus == BRepCheck_NotClosed) {
    printf(" EGADS Fault: Not Closed\n");
  } else if (cStatus == BRepCheck_NotConnected) {
    printf(" EGADS Fault: Not Connected\n");
  } else if (cStatus == BRepCheck_SubshapeNotInShape) {
    printf(" EGADS Fault: Subshape Not In Shape\n");
  } else if (cStatus == BRepCheck_BadOrientation) {
    printf(" EGADS Fault: Bad Orientation\n");
  } else if (cStatus == BRepCheck_BadOrientationOfSubshape) {
    printf(" EGADS Fault: Bad Orientation Of Subshape\n");
#if CASVER >= 680
  } else if (cStatus == BRepCheck_InvalidPolygonOnTriangulation) {
    printf(" EGADS Fault: Invalid Polygon On Triangulation\n");
#endif
  } else if (cStatus == BRepCheck_InvalidToleranceValue) {
    printf(" EGADS Fault: Invalid Tolerance Value\n");
  } else if (cStatus == BRepCheck_CheckFail) {
    printf(" EGADS Fault: Check Fail\n");
  } else {
    printf(" EGADS Unknown Fault = %d\n", cStatus);
  }
}


void
EG_checkStatus(const Handle_BRepCheck_Result tResult)
{

  const BRepCheck_ListOfStatus& tList = tResult->Status();
  if (tList.Extent() <= 0) return;

  BRepCheck_Status cStatus = tList.First();
  EG_printStatus(cStatus);
}


int
EG_destroyTopology(egObject *topo)
{
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;

  if (topo->blind == NULL) return EGADS_SUCCESS;

  if (topo->oclass == MODEL) {
    egadsModel *mshape = (egadsModel *) topo->blind;
    if (mshape->bodies != NULL) {
      for (int i = 0; i < mshape->nbody; i++)
        EG_dereferenceObject(mshape->bodies[i], topo);
      delete [] mshape->bodies;
    }
    mshape->shape.Nullify();
    delete mshape;

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL) {
      if (topo->mtype == WIREBODY) {
        int nwire = pbody->loops.map.Extent();
        for (int i = 0; i < nwire; i++)
          EG_dereferenceObject(pbody->loops.objs[i], topo);
      } else if (topo->mtype == FACEBODY) {
        int nface = pbody->faces.map.Extent();
        for (int i = 0; i < nface; i++)
          EG_dereferenceObject(pbody->faces.objs[i], topo);
      } else {
        int nshell = pbody->shells.map.Extent();
        for (int i = 0; i < nshell; i++)
          EG_dereferenceObject(pbody->shells.objs[i], topo);
        if (topo->mtype == SOLIDBODY) delete [] pbody->senses;
      }
      EG_cleanMaps(&pbody->shells);
      EG_cleanMaps(&pbody->faces);
      EG_cleanMaps(&pbody->loops);
      EG_cleanMaps(&pbody->edges);
      EG_cleanMaps(&pbody->nodes);
      delete pbody;
    }

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL) {
      if (pshell->topFlg == 0) {
        for (int i = 0; i < pshell->nfaces; i++)
          EG_dereferenceObject(pshell->faces[i], topo);
      } else {
        for (int i = 0; i < pshell->nfaces; i++)
          EG_dereferenceTopObj(pshell->faces[i], topo);
      }
      delete [] pshell->faces;
      delete pshell;
    }

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL) {
      if (pface->topFlg == 0) {
        for (int i = 0; i < pface->nloops; i++)
          EG_dereferenceObject(pface->loops[i], topo);
        EG_dereferenceObject(pface->surface, topo);
      } else {
        for (int i = 0; i < pface->nloops; i++)
          EG_dereferenceTopObj(pface->loops[i], topo);
        EG_dereferenceTopObj(pface->surface, topo);
      }
      delete [] pface->senses;
      delete [] pface->loops;
      delete pface;
    }

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL) {
      if (ploop->topFlg == 0) {
        for (int i = 0; i < ploop->nedges; i++) {
          EG_dereferenceObject(ploop->edges[i], topo);
          if (ploop->surface != NULL)
            EG_dereferenceObject(ploop->edges[i+ploop->nedges], topo);
        }
        if (ploop->surface != NULL)
          EG_dereferenceObject(ploop->surface, topo);
      } else {
        for (int i = 0; i < ploop->nedges; i++) {
          EG_dereferenceTopObj(ploop->edges[i], topo);
          if (ploop->surface != NULL)
            EG_dereferenceTopObj(ploop->edges[i+ploop->nedges], topo);
        }
        if (ploop->surface != NULL)
          EG_dereferenceTopObj(ploop->surface, topo);
      }
      delete [] ploop->senses;
      delete [] ploop->edges;
      delete ploop;
    }

  } else if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL) {
      int degen = 0;

      if ((pedge->curve == NULL) &&
          (topo->mtype  == DEGENERATE)) degen = 1;
      if (pedge->topFlg == 0) {
        if (degen == 0) EG_dereferenceObject(pedge->curve, topo);
        EG_dereferenceObject(pedge->nodes[0], topo);
        EG_dereferenceObject(pedge->nodes[1], topo);
      } else {
        if (degen == 0) EG_dereferenceTopObj(pedge->curve, topo);
        EG_dereferenceTopObj(pedge->nodes[0], topo);
        EG_dereferenceTopObj(pedge->nodes[1], topo);
      }
      delete pedge;
    }

  } else {

    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode != NULL) delete pnode;

  }

  return EGADS_SUCCESS;
}


void
EG_splitPeriodics(egadsBody *body)
{
  int          hit    = 0;
  TopoDS_Shape bshape = body->shape;
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(bshape, TopAbs_EDGE, MapE);
  for (int i = 1; i <= MapE.Extent(); i++) {
    TopoDS_Shape shape = MapE(i);
    TopoDS_Edge  Edge  = TopoDS::Edge(shape);
    if (Edge.Closed()) hit++;
  }
  if (hit == 0) {
    TopTools_IndexedMapOfShape MapF;
    TopExp::MapShapes(bshape, TopAbs_FACE, MapF);
    for (int i = 1; i <= MapF.Extent(); i++) {
      TopoDS_Shape shape = MapF(i);
      TopoDS_Face  Face  = TopoDS::Face(shape);
      BRepAdaptor_Surface aSurf(Face, Standard_True);
      if (aSurf.IsUClosed()) hit++;
      if (aSurf.IsVClosed()) hit++;
    }
  }
  if (hit == 0) return;

  // use the OpenCASCADE method ->

  TopoDS_Shape solid = bshape;
  Handle(ShapeBuild_ReShape) reShape = new ShapeBuild_ReShape();
  ShapeUpgrade_ShapeDivideClosed aShape(bshape);
  aShape.SetNbSplitPoints(1);
  aShape.SetContext(reShape);
  if (aShape.Perform(Standard_False)) {
    solid = reShape->Apply(bshape);
    if (solid.IsNull()) {
      printf(" EGADS Warning: Can't Split Periodics!\n");
      solid = bshape;
    } else {
      BRepCheck_Analyzer fCheck(solid);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(solid);
        sfs->Perform();
        TopoDS_Shape fixedSolid = sfs->Shape();
        if (fixedSolid.IsNull()) {
          printf(" EGADS Warning: Periodic Split is Invalid!\n");
          solid = bshape;
        } else {
          BRepCheck_Analyzer sfCheck(fixedSolid);
          if (!sfCheck.IsValid()) {
            printf(" EGADS Warning: Periodic Split is InValid!\n");
            solid = bshape;
          } else {
            solid = fixedSolid;
          }
        }
      }
    }
  }

  body->shape = solid;
}


void
EG_splitMultiplicity(egadsBody *body, int outLevel)
{
  TopoDS_Shape bshape = body->shape;

  int hite = 0;
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(bshape, TopAbs_EDGE, MapE);
  for (int i = 1; i <= MapE.Extent(); i++) {
    Standard_Real t1, t2;
    TopoDS_Shape shape = MapE(i);
    TopoDS_Edge  Edge  = TopoDS::Edge(shape);
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
    Handle(Geom_BSplineCurve) hBSpline =
                                    Handle(Geom_BSplineCurve)::DownCast(hCurve);
    if  (hBSpline.IsNull()) continue;
    if  (hBSpline->Continuity() == GeomAbs_C0)  hite++;
/*  if ((hBSpline->Degree() > 1) &&
        (hBSpline->Continuity() == GeomAbs_C1)) hite++;  */
/*  int c = -1;
    if (hBSpline->Continuity() == GeomAbs_C0) c = 0;
    if (hBSpline->Continuity() == GeomAbs_C1) c = 1;
    if (hBSpline->Continuity() == GeomAbs_C2) c = 2;
    if (hBSpline->Continuity() == GeomAbs_C3) c = 3;
    if (hBSpline->Continuity() == GeomAbs_CN) c = 99;
    printf(" Edge %d: degree = %d, c = %d!\n", i, hBSpline->Degree(), c);  */
  }
  int hitf = 0;
  TopTools_IndexedMapOfShape MapF;
  TopExp::MapShapes(bshape, TopAbs_FACE, MapF);
  for (int i = 1; i <= MapF.Extent(); i++) {
    TopoDS_Shape shape = MapF(i);
    TopoDS_Face  Face  = TopoDS::Face(shape);
    Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
    Handle(Geom_BSplineSurface) hBSpline =
                                Handle(Geom_BSplineSurface)::DownCast(hSurface);
    if   (hBSpline.IsNull()) continue;
    if   (hBSpline->Continuity() == GeomAbs_C0)  hitf++;
/*  if (((hBSpline->UDegree() > 1) || (hBSpline->VDegree() > 1)) &&
         (hBSpline->Continuity() == GeomAbs_C1)) hitf++;  */
  }
  if (hite+hitf == 0) return;

  // use the OpenCASCADE method ->

  TopoDS_Shape shape = bshape;
  ShapeUpgrade_ShapeDivideContinuity aShape(bshape);
/*  aShape.SetBoundaryCriterion(GeomAbs_C2);
    aShape.SetPCurveCriterion(GeomAbs_C2);
    aShape.SetSurfaceCriterion(GeomAbs_C2);  */
  if (aShape.Perform()) {
    if (aShape.Status(ShapeExtend_OK)) {
      printf(" EGADS Warning: No Splitting %d %d (EG_splitMultiplicity)!\n",
             hite, hitf);
      return;
    }
    if (outLevel > 1) {
      if (aShape.Status(ShapeExtend_DONE1))
        printf(" EGADS Info: Some Edges Split!\n");
      if (aShape.Status(ShapeExtend_DONE2))
        printf(" EGADS Info: Surfaces were Split!\n");
      if (aShape.Status(ShapeExtend_DONE3))
        printf(" EGADS Info: Surfaces were modified without splitting!\n");
      if (aShape.Status(ShapeExtend_FAIL1))
        printf(" EGADS Info: Some errors occured splitting Wires!\n");
      if (aShape.Status(ShapeExtend_FAIL2))
        printf(" EGADS Info: Faces could not be split!\n");
    }
    if (!aShape.Status(ShapeExtend_DONE)) {
      printf(" EGADS Warning: Not Done (EG_splitMultiplicity)!\n");
      return;
    }
    shape = aShape.Result();
    if (shape.IsNull()) {
      printf(" EGADS Warning: Can't Do Continuity Split!\n");
      shape = bshape;
    } else {
      BRepCheck_Analyzer fCheck(shape);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          printf(" EGADS Warning: Continuity Split is Invalid!\n");
          shape = bshape;
        } else {
          BRepCheck_Analyzer sfCheck(fixedShape);
          if (!sfCheck.IsValid()) {
            printf(" EGADS Warning: Continuity Split is InValid!\n");
            shape = bshape;
          } else {
            shape = fixedShape;
          }
        }
      }
    }
  } else {
    printf(" EGADS Warning: Perform is False (EG_splitMultiplicity)!\n");
    return;
  }

  body->shape = shape;
}


void
EG_fillPCurves(TopoDS_Face face, egObject *surfo, egObject *loopo,
                                 egObject *topObj)
{
  int           i = 0;
  egObject      *geom;
  Standard_Real f, l;

  egadsLoop *ploop = (egadsLoop *) loopo->blind;
  if (ploop->surface == NULL) return;
  if (ploop->surface != surfo) {
    printf(" EGADS Internal: Loop/Face mismatch on Surface!\n");
    return;
  }

  TopoDS_Wire wire = ploop->loop;
  BRepTools_WireExplorer ExpWE;
  for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
    if (ploop->edges[ploop->nedges+i] != NULL) {
      printf(" EGADS Internal: PCurve already Filled!\n");
      return;
    }
    TopoDS_Shape shape = ExpWE.Current();
    TopoDS_Edge  edge  = TopoDS::Edge(shape);
    if (EG_makeObject(EG_context(surfo),
          &ploop->edges[ploop->nedges+i]) == EGADS_SUCCESS) {
      geom = ploop->edges[ploop->nedges+i];
      Handle(Geom2d_Curve) hCurve = BRep_Tool::
                                    CurveOnSurface(edge, face, f, l);
      if (hCurve.IsNull()) {
        Handle(Geom_Curve)   hCurv    = BRep_Tool::Curve(edge, f, l);
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(face);
        double toler = BRep_Tool::Tolerance(edge);
        try {
          hCurve = GeomProjLib::Curve2d(hCurv, f, l, hSurface, toler);
        }
        catch (Standard_Failure)
        {
          printf(" EGADS Info: Geometry Creation Error (EG_fillPCurves)!\n");
        }
        catch (...)
        {
          printf(" EGADS Info: Geometry Creation Error (EG_fillPCurves)!\n");
        }
      }
      geom->topObj = topObj;
      EG_completePCurve(geom,  hCurve);
      EG_referenceObject(geom, loopo);
    }
    i++;
  }

}


int
EG_shellClosure(egadsShell *pshell, int mtype)
{
  int ret, i, *hits;

  TopoDS_Shell Shell = pshell->shell;
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(Shell, TopAbs_EDGE, MapE);
  if (MapE.Extent() == 0) return CLOSED;

  hits = new int[MapE.Extent()];
  for (i = 0; i < MapE.Extent(); i++) hits[i] = 0;

  TopExp_Explorer ExpW;
  for (ExpW.Init(Shell, TopAbs_EDGE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shape = ExpW.Current();
    TopoDS_Edge  edge  = TopoDS::Edge(shape);
    if (BRep_Tool::Degenerated(edge)) continue;
    i = MapE.FindIndex(edge);
    if (i == 0) {
      printf(" EGADS Internal: Edge not found (EG_shellClosure)!\n");
      continue;
    }
    hits[i-1]++;
  }

  ret = CLOSED;
  for (i = 0; i < MapE.Extent(); i++)
    if ((hits[i] != 2) && (hits[i] != 0)) ret = OPEN;
  if ((mtype == DEGENERATE) && (ret == OPEN))
    for (i = 0; i < MapE.Extent(); i++)
      printf(" EGADS Info: Edge %d: hits = %d\n", i+1, hits[i]);

  delete [] hits;
  return ret;
}


static void
EG_fillTopoObjs(egObject *object, egObject *topObj)
{
  int           outLevel, stat, degen = 0;
  egObject      *context;
  Standard_Real t1, t2;
  TopoDS_Vertex V1, V2;

  outLevel = EG_outLevel(object);
  context  = EG_context(object);

  if (object->oclass == EDGE) {

    egObject *geom   = NULL;
    egObject *pn1    = NULL;
    egObject *pn2    = NULL;
    egadsEdge *pedge = (egadsEdge *) object->blind;
    TopoDS_Edge Edge = pedge->edge;
    if (BRep_Tool::Degenerated(Edge)) {
      degen = 1;
    } else {
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
      stat = EG_makeObject(context, &geom);
      if (stat == EGADS_SUCCESS) {
        geom->topObj = topObj;
        EG_completeCurve(geom, hCurve);
      }
    }

    TopExp::Vertices(Edge, V2, V1, Standard_True);
    EG_makeObject(context, &pn1);
    if (pn1 != NULL) {
      egadsNode *pnode = new egadsNode;
      gp_Pnt pv        = BRep_Tool::Pnt(V1);
      pnode->node      = V1;
      pnode->xyz[0]    = pv.X();
      pnode->xyz[1]    = pv.Y();
      pnode->xyz[2]    = pv.Z();
      pn1->oclass      = NODE;
      pn1->blind       = pnode;
      pn1->topObj      = topObj;
      BRepCheck_Analyzer v1Check(V1);
      if(!v1Check.IsValid())
        if (outLevel > 0)
        printf(" EGADS Info: Node1 may be invalid (EG_fillTopoObjs)!\n");
    }
    if (V1.IsSame(V2)) {
      object->mtype = ONENODE;
      pn2           = pn1;
    } else {
      object->mtype = TWONODE;
      EG_makeObject(context, &pn2);
      if (pn2 != NULL) {
        egadsNode *pnode = new egadsNode;
        gp_Pnt pv        = BRep_Tool::Pnt(V2);
        pnode->node      = V2;
        pnode->xyz[0]    = pv.X();
        pnode->xyz[1]    = pv.Y();
        pnode->xyz[2]    = pv.Z();
        pn2->oclass      = NODE;
        pn2->blind       = pnode;
        pn2->topObj      = topObj;
        BRepCheck_Analyzer v2Check(V2);
        if(!v2Check.IsValid())
          if (outLevel > 0)
          printf(" EGADS Info: Node2 may be invalid (EG_fillTopoObjs)!\n");
      }
    }
    if (Edge.Orientation() != TopAbs_REVERSED) {
      pedge->nodes[0] = pn2;
      pedge->nodes[1] = pn1;
    } else {
      pedge->nodes[0] = pn1;
      pedge->nodes[1] = pn2;
    }

    pedge->curve   = geom;
    pedge->topFlg  = 0;
    object->topObj = topObj;
    if (degen == 1) {
      object->mtype = DEGENERATE;
    } else {
      EG_referenceObject(geom, object);
    }
    EG_referenceObject(pn1,  object);
    EG_referenceObject(pn2,  object);
    BRepCheck_Analyzer eCheck(Edge);
    if (!eCheck.IsValid())
      if (outLevel > 0)
        printf(" EGADS Info: Edge may be invalid (EG_fillTopoObjs)!\n");

  } else if (object->oclass == LOOP) {

    int      *senses = NULL;
    egObject **edgeo = NULL;
    egadsLoop *ploop = (egadsLoop *) object->blind;
    TopoDS_Wire Wire = ploop->loop;
    int            n = 1;
    int           ne = 0;
    int       closed = 0;
    if (Wire.Closed()) closed = 1;
    // more reliable for checking closure of Wires
    TopExp::Vertices(Wire, V1, V2);
    if (!V1.IsNull() && !V2.IsNull())
      if (V1.IsSame(V2)) {
        closed = 1;
      } else {
        closed = 0;
      }
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) ne++;
    if (ploop->surface != NULL) n = 2;

    if (ne > 0) {
      edgeo  = new egObject*[n*ne];
      senses = new int[ne];
    }
    int k = 0;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) {
      if (edgeo == NULL) continue;
      TopoDS_Shape shapW = ExpWE.Current();
      TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
      edgeo[k]           = NULL;
      senses[k]          = 1;
      if (n == 2) edgeo[k+ne] = NULL;
      if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
      for (int j = 0; j < k; j++) {
        egadsEdge *pedg = (egadsEdge *) edgeo[j]->blind;
        if (Edge.IsSame(pedg->edge)) {
          edgeo[k] = edgeo[j];
          break;
        }
      }
      if (edgeo[k] == NULL) {
        stat = EG_makeObject(context, &edgeo[k]);
        if (stat != EGADS_SUCCESS) continue;
        edgeo[k]->oclass = EDGE;
        egadsEdge *pedge = new egadsEdge;
        pedge->edge      = Edge;
        pedge->curve     = NULL;
        pedge->nodes[0]  = NULL;
        pedge->nodes[1]  = NULL;
        edgeo[k]->blind  = pedge;
        EG_fillTopoObjs(edgeo[k], topObj);
      }
      EG_referenceObject(edgeo[k], object);
      k++;
    }

    ploop->nedges  = ne;
    ploop->edges   = edgeo;
    ploop->senses  = senses;
    ploop->topFlg  = 0;
    object->topObj = topObj;
    object->mtype  = OPEN;
    if (closed == 1) object->mtype = CLOSED;
    BRepCheck_Analyzer wCheck(Wire);
    if (!wCheck.IsValid())
      if (outLevel > 0)
        printf(" EGADS Info: Loop may be invalid (EG_fillTopoObjs)!\n");

  } else {

    printf(" EGADS Internal: Not Implemented (EG_fillTopoObjs)!\n");

  }

}


int
EG_traverseBody(egObject *context, int i, egObject *bobj,
                egObject *topObj, egadsBody *body, int *nerr)
{
  int          ii, j, k, outLevel, stat, hit = 0, solid = 0;
  TopoDS_Shape shape;
  egObject     *obj, *geom;

  *nerr    = 0;
  outLevel = EG_outLevel(context);
  if (body->shape.ShapeType() == TopAbs_SOLID) solid = 1;

  TopExp_Explorer Exp;
  TopExp::MapShapes(body->shape, TopAbs_VERTEX, body->nodes.map);
  TopExp::MapShapes(body->shape, TopAbs_EDGE,   body->edges.map);
  TopExp::MapShapes(body->shape, TopAbs_WIRE,   body->loops.map);
  TopExp::MapShapes(body->shape, TopAbs_FACE,   body->faces.map);
  TopExp::MapShapes(body->shape, TopAbs_SHELL,  body->shells.map);
  int nNode  = body->nodes.map.Extent();
  int nEdge  = body->edges.map.Extent();
  int nLoop  = body->loops.map.Extent();
  int nFace  = body->faces.map.Extent();
  int nShell = body->shells.map.Extent();
  bobj->oclass = BODY;
  bobj->mtype  = WIREBODY;
  if (nFace > 0) {
    bobj->mtype = FACEBODY;
    if (nShell > 0) {
      bobj->mtype = SHEETBODY;
      if (solid == 1) bobj->mtype = SOLIDBODY;
    }
  }

  if (outLevel > 1)
    printf(" EGADS Info: Shape %d has %d Nodes, %d Edges, %d Loops, %d Faces and %d Shells\n",
           i+1, nNode, nEdge, nLoop, nFace, nShell);

  // allocate ego storage

  if (nNode > 0) {
    body->nodes.objs = (egObject **) EG_alloc(nNode*sizeof(egObject *));
    if (body->nodes.objs == NULL) return EGADS_MALLOC;
    for (j = 0; j < nNode; j++) {
      stat = EG_makeObject(context, &body->nodes.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nEdge > 0) {
    body->edges.objs = (egObject **) EG_alloc(2*nEdge*sizeof(egObject *));
    if (body->edges.objs == NULL) {
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 2*nEdge; j++) {
      stat = EG_makeObject(context, &body->edges.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nLoop > 0) {
    body->loops.objs = (egObject **) EG_alloc(nLoop*sizeof(egObject *));
    if (body->loops.objs == NULL) {
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nLoop; j++) {
      stat = EG_makeObject(context, &body->loops.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nFace > 0) {
    body->faces.objs = (egObject **) EG_alloc(2*nFace*sizeof(egObject *));
    if (body->faces.objs == NULL) {
      EG_cleanMaps(&body->loops);
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < 2*nFace; j++) {
      stat = EG_makeObject(context, &body->faces.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->faces);
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }
  if (nShell > 0) {
    body->shells.objs = (egObject **) EG_alloc(nShell*sizeof(egObject *));
    if (body->shells.objs == NULL) {
      EG_cleanMaps(&body->faces);
      EG_cleanMaps(&body->loops);
      EG_cleanMaps(&body->edges);
      EG_cleanMaps(&body->nodes);
      return EGADS_MALLOC;
    }
    for (j = 0; j < nShell; j++) {
      stat = EG_makeObject(context, &body->shells.objs[j]);
      if (stat != EGADS_SUCCESS) {
        EG_cleanMaps(&body->shells);
        EG_cleanMaps(&body->faces);
        EG_cleanMaps(&body->loops);
        EG_cleanMaps(&body->edges);
        EG_cleanMaps(&body->nodes);
        return stat;
      }
    }
  }

  // fill our stuff

  for (j = 0; j < nNode; j++) {
    egadsNode *pnode   = new egadsNode;
    obj                = body->nodes.objs[j];
    shape              = body->nodes.map(j+1);
    TopoDS_Vertex Vert = TopoDS::Vertex(shape);
    gp_Pnt pv          = BRep_Tool::Pnt(Vert);
    pnode->node        = Vert;
    pnode->xyz[0]      = pv.X();
    pnode->xyz[1]      = pv.Y();
    pnode->xyz[2]      = pv.Z();
    obj->oclass        = NODE;
    obj->blind         = pnode;
    obj->topObj        = topObj;
  }

  for (j = 0; j < nEdge; j++) {
    int           n1, n2, degen = 0;
    TopoDS_Vertex V1, V2;
    Standard_Real t1, t2;

    egadsEdge *pedge = new egadsEdge;
    obj              = body->edges.objs[j];
    geom             = body->edges.objs[j+nEdge];
    shape            = body->edges.map(j+1);
    geom->topObj     = topObj;
    TopoDS_Edge Edge = TopoDS::Edge(shape);
    if (BRep_Tool::Degenerated(Edge)) {
      degen        = 1;
      geom->oclass = CURVE;
      geom->mtype  = DEGENERATE;
      geom->blind  = NULL;
    } else {
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
      EG_completeCurve(geom, hCurve);
    }

    TopExp::Vertices(Edge, V2, V1, Standard_True);
    if (Edge.Orientation() != TopAbs_REVERSED) {
      n1 = body->nodes.map.FindIndex(V2);
      n2 = body->nodes.map.FindIndex(V1);
    } else {
      n1 = body->nodes.map.FindIndex(V1);
      n2 = body->nodes.map.FindIndex(V2);
    }
    if (outLevel > 2)
      printf(" Edge %d:  nodes = %d %d  degen = %d (%lf, %lf)\n",
             j+1, n1, n2, degen, t1, t2);

    if (((n1 == 0) || (n2 == 0)) && (degen == 0)) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Node(s) not found for Edge (%d %d)!\n", n1, n2);
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
      if (n1 == 0) {
        gp_Pnt pnt1;
        hCurve->D0(t1, pnt1);
        if (outLevel > 0)
          printf("                beg: %lf %lf %lf",
                 pnt1.X(), pnt1.Y(), pnt1.Z());
        double dist = 1.e308;
        for (ii = 0; ii < nNode; ii++) {
          TopoDS_Shape  shapv = body->nodes.map(ii+1);
          TopoDS_Vertex vertv = TopoDS::Vertex(shapv);
          gp_Pnt pv           = BRep_Tool::Pnt(vertv);
          double d            = sqrt((pnt1.X()-pv.X())*(pnt1.X()-pv.X()) +
                                     (pnt1.Y()-pv.Y())*(pnt1.Y()-pv.Y()) +
                                     (pnt1.Z()-pv.Z())*(pnt1.Z()-pv.Z()));
          if (d >= dist) continue;
          dist = d;
          n1   = ii+1;
        }
        if (outLevel > 0)
          printf(" vert = %d dist = %le\n", n1, dist);
      }
      if (n2 == 0) {
        gp_Pnt pnt2;
        hCurve->D0(t2, pnt2);
        if (outLevel > 0)
          printf("                end: %lf %lf %lf",
                 pnt2.X(), pnt2.Y(), pnt2.Z());
        double dist = 1.e308;
        for (ii = 0; ii < nNode; ii++) {
          TopoDS_Shape  shapv = body->nodes.map(ii+1);
          TopoDS_Vertex vertv = TopoDS::Vertex(shapv);
          gp_Pnt pv           = BRep_Tool::Pnt(vertv);
          double d            = sqrt((pnt2.X()-pv.X())*(pnt2.X()-pv.X()) +
                                     (pnt2.Y()-pv.Y())*(pnt2.Y()-pv.Y()) +
                                     (pnt2.Z()-pv.Z())*(pnt2.Z()-pv.Z()));
          if (d >= dist) continue;
          dist = d;
          n2   = ii+1;
        }
        if (outLevel > 0)
          printf(" vert = %d dist = %le\n", n2, dist);
      }
    } else if ((n1 == 0) && (n2 == 0) && (degen == 1)) {
      hit++;
      if (outLevel > 0)
        printf(" EGADS Warning: Node not found for Degen Edge (%d)!\n", n1);
      Bnd_Box Box;
      double bbox[6];
      BRepBndLib::Add(Edge, Box);
      Box.Get(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
      gp_Pnt pnt1(bbox[0], bbox[1], bbox[2]);
      if (outLevel > 0)
        printf("                pnt: %lf %lf %lf",
               pnt1.X(), pnt1.Y(), pnt1.Z());
      double dist = 1.e308;
      for (ii = 0; ii < nNode; ii++) {
        TopoDS_Shape  shapv = body->nodes.map(ii+1);
        TopoDS_Vertex vertv = TopoDS::Vertex(shapv);
        gp_Pnt pv           = BRep_Tool::Pnt(vertv);
        double d            = sqrt((pnt1.X()-pv.X())*(pnt1.X()-pv.X()) +
                                   (pnt1.Y()-pv.Y())*(pnt1.Y()-pv.Y()) +
                                   (pnt1.Z()-pv.Z())*(pnt1.Z()-pv.Z()));
        if (d >= dist) continue;
        dist = d;
        n1   = n2 = ii+1;
      }
      if (outLevel > 0)
        printf(" vert = %d dist = %le\n", n1, dist);
    } else if (((n1 == 0) || (n2 == 0)) && (degen == 1)) {
      if (n1 == 0) n1 = n2;
      if (n2 == 0) n2 = n1;
    }
    egObject *pn1   = body->nodes.objs[n1-1];
    egObject *pn2   = body->nodes.objs[n2-1];

    pedge->edge     = Edge;
    pedge->curve    = geom;
    pedge->nodes[0] = pn1;
    pedge->nodes[1] = pn2;
    pedge->topFlg   = 0;
    obj->oclass     = EDGE;
    obj->blind      = pedge;
    obj->topObj     = topObj;
    obj->mtype      = TWONODE;
    if (n1 == n2) obj->mtype = ONENODE;
    if (degen == 1) {
      obj->mtype = DEGENERATE;
    } else {
      EG_referenceObject(geom, obj);
    }
    EG_referenceObject(pn1, obj);
    EG_referenceObject(pn2, obj);
  }

  for (j = 0; j < nLoop; j++) {
    int      *senses = NULL, closed = 0, ne = 0;
    egObject **edgeo = NULL;

    egadsLoop *ploop = new egadsLoop;
    obj              = body->loops.objs[j];
    shape            = body->loops.map(j+1);
    obj->oclass      = LOOP;
    if (shape.Closed()) closed = 1;
    TopoDS_Wire Wire = TopoDS::Wire(shape);
    // more reliable for checking closure of Wires
    TopoDS_Vertex V1, V2;
    TopExp::Vertices(Wire, V1, V2);
    if (!V1.IsNull() && !V2.IsNull())
      if (V1.IsSame(V2)) {
        closed = 1;
      } else {
        closed = 0;
      }
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) ne++;
    if (outLevel > 2)
      printf(" Loop %d: # edges = %d, closed = %d\n", j+1, ne, closed);

    // find the Face
    TopoDS_Face Face;
    int         hit;
    for (hit = k = 0; k < nFace; k++) {
      TopoDS_Shape shapf = body->faces.map(k+1);
      Face = TopoDS::Face(shapf);
      TopExp_Explorer ExpW;
      for (ExpW.Init(shapf, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
        TopoDS_Shape shapw = ExpW.Current();
        TopoDS_Wire  fwire = TopoDS::Wire(shapw);
        if (fwire.IsSame(Wire)) {
          hit++;
          break;
        }
      }
      if (hit != 0) break;
    }
    if ((hit == 0) && (outLevel > 0) && (nFace != 0))
      printf(" EGADS Internal: Loop without a Face!\n");
    if (hit != 0) {
      geom = body->faces.objs[k+nFace];
      if (geom->oclass != SURFACE) {
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
        geom->topObj = topObj;
        EG_completeSurf(geom, hSurface);
      }
      hit = 2;
      if (geom->mtype == PLANE) hit = 1;
    } else {
      hit = 1;
    }
    if (hit == 1) {
      geom = NULL;
    } else {
      EG_referenceObject(geom, obj);
    }

    if (ne > 0) {
      edgeo  = new egObject*[hit*ne];
      senses = new int[ne];
    }
    k = 0;
#ifdef FACEWIRE
    if (hit != 0) {
      for (ExpWE.Init(Wire, Face); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Shape shapW = ExpWE.Current();
        TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
        int          ed    = body->edges.map.FindIndex(Edge);
        edgeo[k]           = NULL;
        senses[k]          = 1;
        if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
        if (ed != 0) {
          egObject *eobj = body->edges.objs[ed-1];
          edgeo[k]       = eobj;
          if (hit == 2) edgeo[k+ne] = NULL;
          EG_referenceObject(eobj, obj);
        } else {
          printf(" EGADS Warning: Edge not found for Loop!\n");
        }
        if (outLevel > 2)
          printf("        %d  edge = %d   sense = %d\n", k, ed, senses[k]);
        k++;
      }
    } else {
#endif
      for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Shape shapW = ExpWE.Current();
        TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
        int          ed    = body->edges.map.FindIndex(Edge);
        edgeo[k]           = NULL;
        senses[k]          = 1;
        if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
        if (ed != 0) {
          egObject *eobj = body->edges.objs[ed-1];
          edgeo[k]       = eobj;
          if (hit == 2) edgeo[k+ne] = NULL;
          EG_referenceObject(eobj, obj);
        } else {
          printf(" EGADS Warning: Edge not found for Loop!\n");
        }
        if (outLevel > 2)
          printf("        %d  edge = %d   sense = %d\n", k, ed, senses[k]);
        k++;
      }
#ifdef FACEWIRE
    }
#endif
    ploop->loop    = Wire;
    ploop->surface = geom;
    ploop->nedges  = ne;
    ploop->edges   = edgeo;
    ploop->senses  = senses;
    ploop->topFlg  = 0;
    obj->blind     = ploop;
    obj->topObj    = topObj;
    obj->mtype     = OPEN;
    if (closed == 1) obj->mtype = CLOSED;
    if (bobj->mtype == WIREBODY) EG_referenceObject(obj, bobj);
  }

  int *lcnt = (int *) EG_alloc(nLoop*sizeof(int));
  if (lcnt != NULL) for (j = 0; j < nLoop; j++) lcnt[j] = 0;
  for (j = 0; j < nFace; j++) {
    int      *senses = NULL;
    egObject **loopo = NULL;

    egadsFace *pface = new egadsFace;
    obj              = body->faces.objs[j];
    geom             = body->faces.objs[j+nFace];
    shape            = body->faces.map(j+1);
    obj->oclass      = FACE;
    TopoDS_Face Face = TopoDS::Face(shape);
    if (geom->oclass != SURFACE) {
      Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
      geom->topObj = topObj;
      EG_completeSurf(geom, hSurface);
    }
    EG_referenceObject(geom, obj);

    int nl = 0;
    TopExp_Explorer ExpW;
    for (ExpW.Init(shape, TopAbs_WIRE); ExpW.More(); ExpW.Next()) nl++;
    if (outLevel > 2)
      printf(" Face %d: # loops = %d\n", j+1, nl);
    TopoDS_Wire oWire = BRepTools::OuterWire(Face);

    if (nl > 0) {
      loopo  = new egObject*[nl];
      senses = new int[nl];
    }
    k = 0;
    for (ExpW.Init(shape, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  Wire  = TopoDS::Wire(shapw);
      loopo[k]           = NULL;
      senses[k]          = -1;
      if (Wire.IsSame(oWire)) senses[k] = 1;
      int lp = body->loops.map.FindIndex(Wire);
      if (lp != 0) {
        if (lcnt != NULL) lcnt[lp-1]++;
        loopo[k] = body->loops.objs[lp-1];
        if (lcnt != NULL) {
          if (lcnt[lp-1] == 1) {
            EG_fillPCurves(Face, geom, loopo[k], topObj);
          } else {
            senses[k] *= 2;
            if (outLevel > 0)
              printf(" EGADS Info: Loop %d found in 2 Faces (sense = %d in Face %d)!\n",
                   lp, senses[k], j+1);
          }
        } else {
          EG_fillPCurves(Face, geom, loopo[k], topObj);
        }
        EG_referenceObject(loopo[k], obj);
      } else {
        printf(" EGADS Warning: Loop not found for Face!\n");
      }
      if (outLevel > 2)
        printf("        %d  loop = %d     outer = %d\n", k, lp, senses[k]);
      k++;
    }
    pface->face    = Face;
    pface->surface = geom;
    pface->nloops  = nl;
    pface->loops   = loopo;
    pface->senses  = senses;
    pface->topFlg  = 0;
    obj->blind     = pface;
    obj->topObj    = topObj;
    obj->mtype     = SFORWARD;
    if (Face.Orientation() == TopAbs_REVERSED) obj->mtype = SREVERSE;
    if (bobj->mtype == FACEBODY) EG_referenceObject(obj, bobj);
  }
  if (lcnt != NULL) EG_free(lcnt);

  if (nShell > 0) {
    TopoDS_Shell oShell;
    if (solid == 1) {
      TopoDS_Solid Solid = TopoDS::Solid(body->shape);
#if CASVER >= 660
      oShell = BRepClass3d::OuterShell(Solid);
#else
      oShell = BRepTools::OuterShell(Solid);
#endif
      body->senses = new int[nShell];
    }

    for (j = 0; j < nShell; j++) {
      egObject   **faceo = NULL;
      egadsShell *pshell = new egadsShell;
      obj                = body->shells.objs[j];
      shape              = body->shells.map(j+1);
      obj->oclass        = SHELL;
      TopoDS_Shell Shell = TopoDS::Shell(shape);
      if (solid == 1) {
        body->senses[j] = -1;
        if (Shell.IsSame(oShell)) body->senses[j] = 1;
      }

      int nf = 0;
      TopExp_Explorer ExpF;
      for (ExpF.Init(shape, TopAbs_FACE); ExpF.More(); ExpF.Next()) nf++;

      if (nf > 0) faceo = new egObject*[nf];

      k = 0;
      for (ExpF.Init(shape, TopAbs_FACE); ExpF.More(); ExpF.Next()) {
        if (faceo == NULL) continue;
        TopoDS_Shape shapf = ExpF.Current();
        TopoDS_Face  Face  = TopoDS::Face(shapf);
        faceo[k]           = NULL;
        int fa = body->faces.map.FindIndex(Face);
        if (fa != 0) {
          faceo[k] = body->faces.objs[fa-1];
          EG_referenceObject(faceo[k], obj);
        } else {
          printf(" EGADS Warning: Face not found for Shell!\n");
        }
        if (outLevel > 2)
          printf(" Shell %d/%d: Face = %d\n", k, j+1, fa);
        k++;
      }
      pshell->shell  = Shell;
      pshell->nfaces = nf;
      pshell->faces  = faceo;
      pshell->topFlg = 0;
      obj->blind     = pshell;
      obj->topObj    = topObj;
      obj->mtype     = EG_shellClosure(pshell, 0);
      if (bobj->mtype >= SHEETBODY) EG_referenceObject(obj, bobj);
    }
  }

  *nerr = hit;
  return EGADS_SUCCESS;
}


int
EG_getTolerance(const egObject *topo, double *tol)
{
  *tol = 0.0;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass <  NODE)       return EGADS_NOTTOPO;
  if (topo->oclass >= MODEL)      return EGADS_NOTTOPO;

  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode != NULL) *tol = BRep_Tool::Tolerance(pnode->node);
  } else if (topo->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL) *tol = BRep_Tool::Tolerance(pedge->edge);
  } else if (topo->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL)
      for (int i = 0; i < ploop->nedges; i++) {
        egadsEdge *pedge = (egadsEdge *) ploop->edges[i]->blind;
        if (pedge == NULL) continue;
        double toler = BRep_Tool::Tolerance(pedge->edge);
        if (toler > *tol) *tol = toler;
      }
  } else if (topo->oclass == FACE) {
    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL) *tol = BRep_Tool::Tolerance(pface->face);
  } else if (topo->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL)
      for (int i = 0; i < pshell->nfaces; i++) {
        egadsFace *pface = (egadsFace *) pshell->faces[i]->blind;
        if (pface == NULL) continue;
        double toler = BRep_Tool::Tolerance(pface->face);
        if (toler > *tol) *tol = toler;
      }
  } else {
    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL)
      if (topo->mtype == WIREBODY) {
        int nedge = pbody->edges.map.Extent();
        for (int i = 0; i < nedge; i++) {
          TopoDS_Edge Edge = TopoDS::Edge(pbody->edges.map(i+1));
          double toler     = BRep_Tool::Tolerance(Edge);
          if (toler > *tol) *tol = toler;
        }
      } else {
        int nface = pbody->faces.map.Extent();
        for (int i = 0; i < nface; i++) {
          TopoDS_Face Face = TopoDS::Face(pbody->faces.map(i+1));
          double toler     = BRep_Tool::Tolerance(Face);
          if (toler > *tol) *tol = toler;
        }
      }
  }

  return EGADS_SUCCESS;
}


int
EG_setTolerance(const egObject *topo, double tol)
{
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass <  NODE)       return EGADS_NOTTOPO;
  if (topo->oclass >= MODEL)      return EGADS_NOTTOPO;

  ShapeFix_ShapeTolerance sTol;

  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode != NULL)  sTol.SetTolerance(pnode->node, tol);
  } else if (topo->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL)  sTol.SetTolerance(pedge->edge, tol);
  } else if (topo->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL)  sTol.SetTolerance(ploop->loop, tol, TopAbs_EDGE);
  } else if (topo->oclass == FACE) {
    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL)  sTol.SetTolerance(pface->face, tol);
  } else if (topo->oclass == SHELL) {
    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL) sTol.SetTolerance(pshell->shell, tol, TopAbs_FACE);
  } else {
    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL)  sTol.SetTolerance(pbody->shape, tol, TopAbs_SHELL);
  }

  return EGADS_SUCCESS;
}


int
EG_getTopology(const egObject *topo, egObject **geom, int *oclass,
               int *type, /*@null@*/ double *limits, int *nChildren,
               egObject ***children, int **senses)
{
  *geom      = NULL;
  *oclass    = *type = 0;
  *nChildren = 0;
  *children  = NULL;
  *senses    = NULL;
  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass < NODE)        return EGADS_NOTTOPO;
  *oclass = topo->oclass;
  *type   = topo->mtype;

  if (topo->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) topo->blind;
    if ((limits != NULL) && (pnode != NULL)) {
      limits[0] = pnode->xyz[0];
      limits[1] = pnode->xyz[1];
      limits[2] = pnode->xyz[2];
    }

  } else if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    if (pedge != NULL) {
      *geom      = pedge->curve;
      *nChildren = 1;
      if (topo->mtype == TWONODE) *nChildren = 2;
      *children = pedge->nodes;
      if (limits != NULL)
        BRep_Tool::Range(pedge->edge, limits[0], limits[1]);
    }

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    if (ploop != NULL) {
      *geom      = ploop->surface;
      *nChildren = ploop->nedges;
      *children  = ploop->edges;
      *senses    = ploop->senses;
    }

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    if (pface != NULL) {
      *geom      = pface->surface;
      *nChildren = pface->nloops;
      *children  = pface->loops;
      *senses    = pface->senses;
      if (limits != NULL) {
        Standard_Real umin, umax, vmin, vmax;

        BRepTools::UVBounds(pface->face, umin, umax, vmin, vmax);
        limits[0] = umin;
        limits[1] = umax;
        limits[2] = vmin;
        limits[3] = vmax;
      }
    }

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    if (pshell != NULL) {
      *nChildren = pshell->nfaces;
      *children  = pshell->faces;
    }

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    if (pbody != NULL)
      if (topo->mtype == WIREBODY) {
        *nChildren = pbody->loops.map.Extent();
        *children  = pbody->loops.objs;
      } else if (topo->mtype == FACEBODY) {
        *nChildren = pbody->faces.map.Extent();
        *children  = pbody->faces.objs;
      } else {
        *nChildren = pbody->shells.map.Extent();
        *children  = pbody->shells.objs;
        if (topo->mtype == SOLIDBODY) *senses = pbody->senses;
      }

  } else {

    egadsModel *pmodel = (egadsModel *) topo->blind;
    if (pmodel != NULL) {
      *nChildren = pmodel->nbody;
      *children  = pmodel->bodies;
    }

  }

  return EGADS_SUCCESS;
}


static void
EG_makePCurves(TopoDS_Face face, egObject *surfo, egObject *loopo,
               Standard_Real prec, int flag)
{
  int      i = 0;
  egObject *geom;

  egadsLoop *ploop = (egadsLoop *) loopo->blind;
  if (ploop->surface == NULL) return;
  if (ploop->surface != surfo) {
    printf(" EGADS Internal: Loop/Face mismatch on Surface (EG_makePCurves)!\n");
    return;
  }

  int outLevel     = EG_outLevel(surfo);
  TopoDS_Wire wire = ploop->loop;
  BRep_Builder           Builder;
  BRepTools_WireExplorer ExpWE;
  for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next(), i++) {
    geom                = ploop->edges[ploop->nedges+i];
    egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
    TopoDS_Shape shape  = ExpWE.Current();
    TopoDS_Edge  edge   = TopoDS::Edge(shape);
    if (ppcurv == NULL) continue;
    Handle(Geom2d_Curve) hCurv2d = ppcurv->handle;

    if ((flag == 0) || (outLevel > 2)) {
      TopoDS_Vertex V1, V2;
      Standard_Real t, t1, t2, delta, mdelta;
      gp_Pnt        pnt, pnte, pv1, pv2;
      gp_Pnt2d      uv;

      egadsSurface *psurf = (egadsSurface *) surfo->blind;
      Handle(Geom_Surface) hSurface = psurf->handle;
      if (edge.Orientation() == TopAbs_REVERSED) {
        TopExp::Vertices(edge, V2, V1, Standard_True);
      } else {
        TopExp::Vertices(edge, V1, V2, Standard_True);
      }
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
      pv1 = BRep_Tool::Pnt(V1);
      pv2 = BRep_Tool::Pnt(V2);
      if (outLevel > 2)
        printf(" PCurve #%d: Limits = %lf %lf    prec = %le\n",
               i, t1, t2, prec);
      hCurv2d->D0(t1, uv);
      hSurface->D0(uv.X(), uv.Y(), pnt);
      mdelta = sqrt((pnt.X()-pv1.X())*(pnt.X()-pv1.X()) +
                    (pnt.Y()-pv1.Y())*(pnt.Y()-pv1.Y()) +
                    (pnt.Z()-pv1.Z())*(pnt.Z()-pv1.Z()));
      if (outLevel > 2)
        printf("            delta for 1st Node     = %le  %lf %lf %lf\n",
               mdelta, pv1.X(), pv1.Y(), pv1.Z());
      delta = 0.0;
      for (int j = 1; j < 36; j++) {
        t = t1 + j*(t2-t1)/36.0;
        hCurv2d->D0(t, uv);
        hSurface->D0(uv.X(), uv.Y(), pnt);
        if (BRep_Tool::Degenerated(edge)) {
          pnte = pv1;
        } else {
          hCurve->D0(t, pnte);
        }
        delta += sqrt((pnt.X()-pnte.X())*(pnt.X()-pnte.X()) +
                      (pnt.Y()-pnte.Y())*(pnt.Y()-pnte.Y()) +
                      (pnt.Z()-pnte.Z())*(pnt.Z()-pnte.Z()));
      }
      delta /= 35;
      if (outLevel > 2)
        printf("            ave delta against Edge = %le\n", delta);
      if (delta > mdelta) mdelta = delta;
      hCurv2d->D0(t2, uv);
      hSurface->D0(uv.X(), uv.Y(), pnt);
      delta = sqrt((pnt.X()-pv2.X())*(pnt.X()-pv2.X()) +
                   (pnt.Y()-pv2.Y())*(pnt.Y()-pv2.Y()) +
                   (pnt.Z()-pv2.Z())*(pnt.Z()-pv2.Z()));
      if (outLevel > 2)
        printf("            delta for 2nd Node     = %le  %lf %lf %lf\n",
               delta, pv2.X(), pv2.Y(), pv2.Z());
      if (delta > mdelta) mdelta = delta;
      if ((flag == 0) && (mdelta*1.001 > prec))
        Builder.SameParameter(edge, Standard_False);
    }

    Builder.UpdateEdge(edge, hCurv2d, face, prec);
  }

}


int
EG_makeTopology(egObject *context, /*@null@*/ egObject *geom,
                int oclass, int mtype, /*@null@*/ double *limits,
                int nChildren, /*@null@*/ egObject **children,
                /*@null@*/ int *senses, egObject **topo)
{
  int      i, n, stat, outLevel, nerr;
  egCntxt  *cntx;
  egObject *obj;

  *topo = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                  return EGADS_NODATA;
  outLevel = cntx->outLevel;

  if ((oclass < NODE) || (oclass > MODEL)) {
    if (outLevel > 0)
      printf(" EGADS Error: oclass = %d (EG_makeTopology)!\n", oclass);
    return EGADS_NOTTOPO;
  }

  if (oclass == NODE) {

    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Node with no Data (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }

    gp_Pnt pnt(limits[0], limits[1], limits[2]);
    TopoDS_Vertex vert = BRepBuilderAPI_MakeVertex(pnt);
    BRepCheck_Analyzer vCheck(vert);
    if (!vCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Node is invalid (EG_makeTopology)!\n");
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Node object (EG_makeTopology)!\n");
      return stat;
    }
    egadsNode *pnode   = new egadsNode;
    pnode->node        = vert;
    pnode->xyz[0]      = limits[0];
    pnode->xyz[1]      = limits[1];
    pnode->xyz[2]      = limits[2];
    obj->oclass        = NODE;
    obj->blind         = pnode;
    obj->topObj        = context;
    EG_referenceObject(obj, context);

  } else if (oclass == EDGE) {
    TopoDS_Vertex V1, V2;
    Standard_Real P1, P2;

    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Limits is NULL (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    if (limits[0] >= limits[1]) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Tmin (%lf) >= Tmax (%lf) (EG_makeTopology)!\n",
               limits[0], limits[1]);
      return EGADS_RANGERR;
    }

    if (mtype == DEGENERATE) {
      if (nChildren != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge with %d Verts (EG_makeTopology)!\n",
                 nChildren);
        return EGADS_TOPOERR;
      }
      if (children[0] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ Vert NULL (EG_makeTopology)!\n");
        return EGADS_NULLOBJ;
      }
      if (children[0]->oclass != NODE) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ nonNode Child (EG_makeTopology)!\n");
        return EGADS_NOTTOPO;
      }
      if (children[0]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Degen Edge w/ NULL Node Child (EG_makeTopology)!\n");
        return EGADS_NODATA;
      }

      // make a degenerate Edge
      egadsNode *pnode = (egadsNode *) children[0]->blind;
      V1 = pnode->node;
      P1 = limits[0];
      P2 = limits[1];
#if CASVER > 690
      gp_Pnt pv = BRep_Tool::Pnt(V1);
      // construct BSpline of (approx) circle with radius 1.e-6
      TColStd_Array1OfReal    aKnots(1, 2);
      TColStd_Array1OfInteger aMults(1, 2);
      TColgp_Array1OfPnt      aPoles(1, 7);
      aMults(1) = 7;
      aKnots(1) = P1;
      aMults(2) = 7;
      aKnots(2) = P2;
      aPoles(1) = gp_Pnt(pv.X(), pv.Y(), pv.Z());
      aPoles(2) = gp_Pnt(pv.X(), pv.Y()+1.04719755119660e-06, pv.Z());
      aPoles(3) = gp_Pnt(pv.X()-1.20996765359363e-06,
                         pv.Y()+2.40661999743566e-06, pv.Z());
      aPoles(4) = gp_Pnt(pv.X()-4.57635593428032e-06, pv.Y(), pv.Z());
      aPoles(5) = gp_Pnt(pv.X()-1.20996765359363e-06,
                         pv.Y()-2.40661999743566e-06, pv.Z());
      aPoles(6) = gp_Pnt(pv.X(), pv.Y()-1.04719755119660e-06, pv.Z());
      aPoles(7) = gp_Pnt(pv.X(), pv.Y(), pv.Z());
      Handle(Geom_Curve) hCurve = new Geom_BSplineCurve(aPoles, aKnots, aMults,
                                                        6, Standard_False);
      BRepBuilderAPI_MakeEdge MEdge(hCurve, V1, V1, P1, P2);
      TopoDS_Edge Edge = MEdge.Edge();
#else
      // construct circle with radius zero
      gp_Pnt pv = BRep_Tool::Pnt(V1);
      gp_Ax2 axi2(pv, gp_Dir(1.0, 0.0, 0.0), gp_Dir(0.0, 1.0, 0.0));
      Handle(Geom_Curve) hCurve = new Geom_Circle(axi2, 0.0);
      BRepBuilderAPI_MakeEdge MEdge(hCurve, V1, V1, P1, P2);
      TopoDS_Edge Edge = MEdge.Edge();
#if CASVER < 680
      BRep_Builder Builder;
      Builder.Degenerated(Edge, Standard_True);
      if (!BRep_Tool::Degenerated(Edge))
        printf(" EGADS Info: Degenerate Edge NOT Degenerate!\n");
#endif
#endif
      BRepCheck_Analyzer eCheck(Edge);
      if (!eCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Info: Degen Edge is invalid (EG_makeTopology)!\n");
        EG_checkStatus(eCheck.Result(Edge));
        return EGADS_CONSTERR;
      }
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Degen Edge obj (EG_makeTopology)!\n");
        return stat;
      }
      egadsEdge *pedge = new egadsEdge;
      pedge->edge      = Edge;
      pedge->curve     = NULL;
      pedge->nodes[0]  = children[0];
      pedge->nodes[1]  = children[0];
      pedge->topFlg    = 1;
      obj->oclass      = EDGE;
      obj->blind       = pedge;
      obj->topObj      = context;
      obj->mtype       = DEGENERATE;
      EG_referenceTopObj(pedge->nodes[0],  obj);
      EG_referenceTopObj(pedge->nodes[1],  obj);
      EG_referenceObject(obj,          context);

      *topo = obj;
      return EGADS_SUCCESS;
    }

    if (geom == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with NULL Geom (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (geom->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with No Geom (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    if (geom->oclass != CURVE) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge Geom not CURVE (EG_makeTopology)!\n");
      return EGADS_NOTGEOM;
    }
    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if ((nChildren != 1) && (nChildren != 2)) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with %d Verts (EG_makeTopology)!\n",
               nChildren);
      return EGADS_TOPOERR;
    }
    if (children[0] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with Vert[0] NULL (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (children[0]->oclass != NODE) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge with nonNode Child[0] (EG_makeTopology)!\n");
      return EGADS_NOTTOPO;
    }
    if (children[0]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge w/ NULL Node Child[0] (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    if (nChildren == 2) {
      if (children[1] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge with Vert[1] NULL (EG_makeTopology)!\n");
        return EGADS_NULLOBJ;
      }
      if (children[1]->oclass != NODE) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge w/ nonNode Child[1] (EG_makeTopology)!\n");
        return EGADS_NOTTOPO;
      }
      if (children[1]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge w/ NULL Node Child[1] (EG_makeTopology)!\n");
        return EGADS_NODATA;
      }
    }
    egadsCurve *pcurve = (egadsCurve *) geom->blind;
    egadsNode  *pnode1 = (egadsNode *)  children[0]->blind;
    egadsNode  *pnode2 = (egadsNode *)  children[0]->blind;
    if (nChildren == 2) pnode2 = (egadsNode *) children[1]->blind;

    P1 = limits[0];
    P2 = limits[1];
    V1 = pnode1->node;
    V2 = pnode2->node;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    gp_Pnt pnt1, pnt2, pv1, pv2;
    hCurve->D0(P1, pnt1);
    hCurve->D0(P2, pnt2);
    pv1 = BRep_Tool::Pnt(V1);
    pv2 = BRep_Tool::Pnt(V2);
    if (outLevel > 2) {
      printf(" P1 = %lf %lf %lf  %lf %lf %lf\n", pnt1.X(), pnt1.Y(), pnt1.Z(),
                pnode1->xyz[0], pnode1->xyz[1], pnode1->xyz[2]);
      printf("      vert = %lf %lf %lf\n", pv1.X(), pv1.Y(), pv1.Z());
      printf(" P2 = %lf %lf %lf  %lf %lf %lf\n", pnt2.X(), pnt2.Y(), pnt2.Z(),
              pnode2->xyz[0], pnode2->xyz[1], pnode2->xyz[2]);
      printf("      vert = %lf %lf %lf\n", pv2.X(), pv2.Y(), pv2.Z());
    }
    double delta1 = sqrt((pnt1.X()-pv1.X())*(pnt1.X()-pv1.X()) +
                         (pnt1.Y()-pv1.Y())*(pnt1.Y()-pv1.Y()) +
                         (pnt1.Z()-pv1.Z())*(pnt1.Z()-pv1.Z()));
    double delta2 = sqrt((pnt2.X()-pv2.X())*(pnt2.X()-pv2.X()) +
                         (pnt2.Y()-pv2.Y())*(pnt2.Y()-pv2.Y()) +
                         (pnt2.Z()-pv2.Z())*(pnt2.Z()-pv2.Z()));

    Standard_Real old  = BRepBuilderAPI::Precision();
    Standard_Real prec = old;
    if (outLevel > 2)
      printf("   Limits = %f %lf, Tol = %le %le   %le\n",
             P1, P2, delta1, delta2, old);
    if (delta1*1.001 > prec) prec = 1.001*delta1;
    if (delta2*1.001 > prec) prec = 1.001*delta2;
    BRepBuilderAPI::Precision(prec);
    BRepBuilderAPI_MakeEdge MEdge;
    MEdge.Init(hCurve, V1, V2, P1, P2);
    BRepBuilderAPI::Precision(old);
    if (!MEdge.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with the Edge (EG_makeTopology)!\n");
       return EGADS_NODATA;
    }
    TopoDS_Edge Edge = MEdge.Edge();
    BRepCheck_Analyzer eCheck(Edge);
    if (!eCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Edge is invalid (EG_makeTopology)!\n");
      EG_checkStatus(eCheck.Result(Edge));
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Edge object (EG_makeTopology)!\n");
      return stat;
    }
    egadsEdge *pedge = new egadsEdge;
    pedge->edge      = Edge;
    pedge->curve     = geom;
    pedge->nodes[0]  = children[0];
    pedge->nodes[1]  = children[0];
    pedge->topFlg    = 1;
    obj->oclass      = EDGE;
    obj->blind       = pedge;
    obj->topObj      = context;
    obj->mtype       = ONENODE;
    if (nChildren == 2) {
      obj->mtype      = TWONODE;
      pedge->nodes[1] = children[1];
    }
    EG_referenceTopObj(geom,             obj);
    EG_referenceTopObj(pedge->nodes[0],  obj);
    EG_referenceTopObj(pedge->nodes[1],  obj);
    EG_referenceObject(obj,          context);

  } else if (oclass == LOOP) {

    if ((children == NULL) || (senses == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop with NULL Input (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop with %d Edges (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    n = 1;
    if (geom != NULL) {
      if (geom->oclass != SURFACE) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop Geom not SURFACE (EG_makeTopology)!\n");
        return EGADS_NOTGEOM;
      }
      if (geom->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop with No Geom Data (EG_makeTopology)!\n");
        return EGADS_NODATA;
      }
      if (geom->mtype == PLANE) {
        if (outLevel > 0)
          printf(" EGADS Info: Loop with Planar Surface (EG_makeTopology)!\n");
      } else {
        n = 2;
      }
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop with Edge[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != EDGE) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ nonEdge Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ NULL Edge Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if (n == 1) continue;
      if (children[i+nChildren] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop with PCurve[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i+nChildren]->oclass != PCURVE) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ nonPCurve Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i+nChildren]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop w/ NULL PCurve Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
    }
    BRepBuilderAPI_MakeWire MW;
    for (i = 0; i < nChildren; i++) {
      egadsEdge *pedge = (egadsEdge *) children[i]->blind;
      if ((children[i]->mtype == DEGENERATE) && (geom != NULL))
        if (!BRep_Tool::Degenerated(pedge->edge)) {
          // fix up degenerate edge if not already "Degenerated"
          egadsSurface *psurf = (egadsSurface *) geom->blind;
          egadsPCurve  *pcurv = (egadsPCurve *)  children[i+nChildren]->blind;
          TopLoc_Location L;
          BRep_Builder    Builder;
          Builder.Degenerated(pedge->edge, Standard_True);
          // remove 3D curve
          Builder.UpdateEdge(pedge->edge, NULL, Precision::Confusion());
          // set PCurve & Surface
          Builder.UpdateEdge(pedge->edge, pcurv->handle, psurf->handle, L,
                             Precision::Confusion());
          if (!BRep_Tool::Degenerated(pedge->edge))
            printf(" EGADS Info: Degenerate Edge %d in Loop NOT Degenerate!\n",
                   i+1);
          BRepCheck_Analyzer eCheck(pedge->edge);
          if (!eCheck.IsValid()) {
            if (outLevel > 0)
              printf(" EGADS Info: Degen Edge %d in Loop is invalid!\n", i+1);
            EG_checkStatus(eCheck.Result(pedge->edge));
            return EGADS_CONSTERR;
          }
        }
      TopoDS_Edge edge = pedge->edge;

      // may only be required for the first Edge, must be in order!
      if (edge.Orientation() == TopAbs_REVERSED) {
        if (senses[i] ==  1) edge.Orientation(TopAbs_FORWARD);
      } else {
        if (senses[i] == -1) edge.Orientation(TopAbs_REVERSED);
      }
/*
      if (edge.Orientation() == TopAbs_REVERSED) {
        if (senses[i] ==  1) {
          TopoDS_Shape shape = edge.Oriented(TopAbs_FORWARD);
          edge = TopoDS::Edge(shape);
        }
      } else {
        if (senses[i] == -1) {
          TopoDS_Shape shape = edge.Oriented(TopAbs_REVERSED);
          edge = TopoDS::Edge(shape);
        }
      }
*/
      MW.Add(edge);
      if (MW.Error()) {
        if (outLevel > 0)
          printf(" EGADS Error: Problem with Edge %d (EG_makeTopology)!\n",
                 i+1);
        return EGADS_NODATA;
      }
    }
    if (!MW.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with Loop (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    TopoDS_Wire wire = MW.Wire();

    // validate against the senses
    if (outLevel > 2) {
      i = 0;
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Shape shape = ExpWE.Current();
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        int sense = 1;
        if (shape.Orientation() == TopAbs_REVERSED) sense = -1;
        egadsEdge *pedge = (egadsEdge *) children[i]->blind;
        if (Edge.IsSame(pedge->edge)) {
          printf("  %d: Edges same senses = %d %d\n",
                 i, senses[i], sense);
        } else {
          printf("  %d: Edges NOT the same senses = %d %d\n",
                 i, senses[i], sense);
        }
        i++;
      }
    }

    BRepCheck_Analyzer wCheck(wire);
    if (!wCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Wire is invalid (EG_makeTopology)!\n");
      EG_checkStatus(wCheck.Result(wire));
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Loop object (EG_makeTopology)!\n");
      return stat;
    }
    obj->oclass = LOOP;

    egadsLoop *ploop  = new egadsLoop;
    egObject  **edgeo = new egObject*[n*nChildren];
    int       *esense = new int[nChildren];
    int       closed  = 0;
    if (wire.Closed()) closed = 1;
    // more reliable for checking closure of Wires
    TopoDS_Vertex V1, V2;
    TopExp::Vertices(wire, V1, V2);
    if (!V1.IsNull() && !V2.IsNull())
      if (V1.IsSame(V2)) {
        closed = 1;
      } else {
        closed = 0;
      }
    for (i = 0; i < nChildren; i++) {
      edgeo[i]  = children[i];
      esense[i] = senses[i];
      EG_referenceTopObj(children[i], obj);
      if (n == 1) continue;
      edgeo[i+nChildren] = children[i+nChildren];
      EG_referenceTopObj(children[i+nChildren], obj);
    }
    ploop->loop    = wire;
    ploop->surface = geom;
    ploop->nedges  = nChildren;
    ploop->edges   = edgeo;
    ploop->senses  = esense;
    ploop->topFlg  = 1;
    obj->blind     = ploop;
    obj->topObj    = context;
    obj->mtype     = OPEN;
    if (closed == 1) obj->mtype = CLOSED;
    EG_referenceObject(obj, context);
    if (geom != NULL)
      if (geom->mtype == PLANE) {
        ploop->surface = NULL;
      } else {
        EG_referenceTopObj(geom, obj);
      }
    if (mtype == CLOSED)
      if ((outLevel > 0) && (closed == 0))
        printf(" EGADS Info: Wire is Open (EG_makeTopology)!\n");
    if (mtype == OPEN)
      if ((outLevel > 0) && (closed == 1))
        printf(" EGADS Info: Wire is Closed (EG_makeTopology)!\n");

  } else if (oclass == FACE) {

    if ((mtype != SFORWARD) && (mtype != SREVERSE)) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with MType = %d (EG_makeTopology)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    if (geom == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with NULL Geom (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (geom->oclass != SURFACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Geom not SURFACE (EG_makeTopology)!\n");
      return EGADS_NOTGEOM;
    }
    if (geom->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with No Geom (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    egadsSurface *psurf = (egadsSurface *) geom->blind;
    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Face with %d Loops (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Face with Loop[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != LOOP) {
        if (outLevel > 0)
          printf(" EGADS Error: Face w/ nonLoop Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->mtype != CLOSED) {
        if (outLevel > 0)
          printf(" EGADS Error: Face with OPEN Loop[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Face w/ NULL Loop Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      egadsLoop *ploop = (egadsLoop *) children[i]->blind;
      if ((ploop->surface != geom) && (geom->mtype != PLANE)) {
        if (outLevel > 0)
          printf(" EGADS Error: Face/Loop[%d] Geom Mismatch (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTGEOM;
      }
    }

    TopoDS_Face face;
    BRepBuilderAPI_MakeFace MFace;
    Standard_Real old = BRepBuilderAPI::Precision();
    Standard_Real prc = old;
    int nTrys = 5;              // # of tol attempts before setting SameParam = F
    for (int itry = 0; itry < nTrys; itry++) {
      if (itry != 0) face.Nullify();
      if (itry == nTrys-1) {    // last attempt -- set tol back
        prc = old;
        BRepBuilderAPI::Precision(old);
      }
//    this can change the underlying surface!
//             a copy prevents this when the surface is shared between Faces
/*    gp_Trsf               form  = gp_Trsf();
      Handle(Geom_Surface)  hSurf = psurf->handle;
      Handle(Geom_Geometry) nGeom = hSurf->Transformed(form);
      Handle(Geom_Surface)  nSurf = Handle(Geom_Surface)::DownCast(nGeom);
 */
      Handle(Geom_Surface)  nSurf = psurf->handle;
#if CASVER >= 652
      MFace.Init(nSurf, Standard_False, prc);
#else
      MFace.Init(nSurf, Standard_False);
#endif
      for (i = 0; i < nChildren; i++) {
        egadsLoop *ploop = (egadsLoop *) children[i]->blind;
        TopoDS_Wire wire = ploop->loop;
        if (mtype == SREVERSE) wire.Reverse();
        MFace.Add(wire);
        if (MFace.Error()) {
          if (outLevel > 0)
            printf(" EGADS Error: Problem with Loop %d (EG_makeTopology)!\n",
                   i+1);
          return EGADS_NODATA;
        }
      }
      if (MFace.IsDone()) {
        face = MFace.Face();
        if (mtype == SREVERSE) {
          face.Orientation(TopAbs_REVERSED);
        } else {
          face.Orientation(TopAbs_FORWARD);
        }
        for (i = 0; i < nChildren; i++)
          EG_makePCurves(face, geom, children[i], prc, nTrys-itry-1);
        BRepLib::SameParameter(face);
        BRepCheck_Analyzer oCheck(face);
        if (oCheck.IsValid()) break;
      }
      prc *= 10.0;
      BRepBuilderAPI::Precision(prc);
      if (outLevel > 1)
        printf(" EGADS Info: Adjusting Precision for Face - itry = %d  prec = %lf\n",
               itry, prc);
    }
    BRepBuilderAPI::Precision(old);
    if (!MFace.IsDone()) {
      if (outLevel > 0)
        printf(" EGADS Error: Problem with the Face (EG_makeTopology)!\n");
      return EGADS_NODATA;
    }
    for (i = 0; i < nChildren; i++) {
      egadsLoop *ploop = (egadsLoop *) children[i]->blind;
      BRepCheck_Wire wireCheck(ploop->loop);
      wireCheck.InContext(face);
      BRepCheck_ListOfStatus stati;
      stati = wireCheck.Status();
      if (stati.Extent() > 0) {
        BRepCheck_Status cStatus = stati.First();
        if (outLevel > 0)
          EG_printStatus(cStatus);
        if (cStatus != BRepCheck_NoError) return EGADS_GEOMERR;
      }
    }
    BRepCheck_Analyzer fCheck(face);
    if (!fCheck.IsValid()) {
      // try to fix the fault
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(face);
      sfs->Perform();
      TopoDS_Shape fixedFace = sfs->Shape();
      if (fixedFace.IsNull()) {
        if (outLevel > 0) {
          printf(" EGADS Info: Invalid Face w/ NULL Fix (EG_makeTopology)!\n");
          EG_checkStatus(fCheck.Result(face));
        }
        return  EGADS_CONSTERR;
      }
      BRepCheck_Analyzer fxCheck(fixedFace);
      if (!fxCheck.IsValid()) {
        if (outLevel > 0) {
          printf(" EGADS Info: Face is invalid (EG_makeTopology)!\n");
          EG_checkStatus(fxCheck.Result(fixedFace));
        }
        return  EGADS_CONSTERR;
      }
      face = TopoDS::Face(fixedFace);
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeTopology)!\n");
      return stat;
    }
    obj->oclass = FACE;

    egadsFace *pface = new egadsFace;
    egObject **loopo = new egObject*[nChildren];
    int      *lsense = new int[nChildren];
    for (i = 0; i < nChildren; i++) {
      loopo[i]  = children[i];
      lsense[i] = senses[i];
      EG_referenceTopObj(children[i], obj);
    }
    pface->face    = face;
    pface->surface = geom;
    pface->nloops  = nChildren;
    pface->loops   = loopo;
    pface->senses  = lsense;
    pface->topFlg  = 1;
    obj->blind     = pface;
    obj->topObj    = context;
    obj->mtype     = mtype;
    EG_referenceTopObj(geom, obj);
    EG_referenceObject(obj,  context);

  } else if (oclass == SHELL) {

    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Shell with NULL Input (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Shell with %d Faces (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell with Face[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != FACE) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell w/ nonFace Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Shell w/ NULL Face Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
    }
    BRep_Builder builder3D;
    TopoDS_Shell shell;
    builder3D.MakeShell(shell);
    for (i = 0; i < nChildren; i++) {
      egadsFace *pface = (egadsFace *) children[i]->blind;
      builder3D.Add(shell, pface->face);
    }
    BRepLib::SameParameter(shell);
    BRepCheck_Analyzer shCheck(shell);
    if (!shCheck.IsValid()) {
      if (outLevel > 0) {
        printf(" EGADS Info: Shell is invalid (EG_makeTopology)!\n");
        EG_checkStatus(shCheck.Result(shell));
      }
      if (mtype != DEGENERATE) return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Shell object (EG_makeTopology)!\n");
      return stat;
    }
   obj->oclass = SHELL;

    egadsShell *pshell = new egadsShell;
    egObject   **faceo = new egObject*[nChildren];
    for (i = 0; i < nChildren; i++) {
      faceo[i] = children[i];
      EG_referenceTopObj(children[i], obj);
    }
    pshell->shell  = shell;
    pshell->nfaces = nChildren;
    pshell->faces  = faceo;
    pshell->topFlg = 1;
    obj->blind     = pshell;
    obj->topObj    = context;
    obj->mtype     = EG_shellClosure(pshell, mtype);
    EG_referenceObject(obj, context);
    if (mtype == CLOSED)
      if ((outLevel > 0) && (obj->mtype == OPEN))
        printf(" EGADS Info: Shell is Open (EG_makeTopology)!\n");
    if (mtype == OPEN)
      if ((outLevel > 0) && (obj->mtype == CLOSED))
        printf(" EGADS Info: Shell is Closed (EG_makeTopology)!\n");
    if (mtype == DEGENERATE)
      if (obj->mtype == OPEN) {
        printf(" EGADS Info: Shell is Open (EG_makeTopology)!\n");
      } else {
        printf(" EGADS Info: Shell is Closed (EG_makeTopology)!\n");
      }

  } else if (oclass == BODY) {

    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Body with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Body with %d Children (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    if ((mtype < WIREBODY) || (mtype > SOLIDBODY)) {
      if (outLevel > 0)
        printf(" EGADS Error: Body with mtype = %d (EG_makeTopology)!\n",
               mtype);
      return EGADS_RANGERR;
    }
    int                    cclass = SHELL;
    if (mtype == FACEBODY) cclass = FACE;
    if (mtype == WIREBODY) cclass = LOOP;
    if ((mtype != SOLIDBODY) && (nChildren != 1)) {
      if (outLevel > 0)
        printf(" EGADS Error: non SolidBody w/ %d children (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Body with child[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != cclass) {
        if (outLevel > 0)
          printf(" EGADS Error: Body w/ invalid Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Body with NULL Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
      if ((children[i]->mtype != CLOSED) && (mtype == SOLIDBODY)) {
        if (outLevel > 0)
          printf(" EGADS Error: Solid w/ nonClosed Shell[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_RANGERR;
      }
    }
    int          fixedSld = 0;
    TopoDS_Shape shape;
    if (mtype == WIREBODY) {
      egadsLoop *ploop = (egadsLoop *) children[0]->blind;
      shape = ploop->loop;
    } else if (mtype == FACEBODY) {
      egadsFace *pface = (egadsFace *) children[0]->blind;
      shape = pface->face;
    } else if (mtype == SHEETBODY) {
      egadsShell *pshell = (egadsShell *) children[0]->blind;
      shape = pshell->shell;
    } else {
      BRep_Builder builder3D;
      TopoDS_Solid solid;
      builder3D.MakeSolid(solid);
      for (i = 0; i < nChildren; i++) {
        egadsShell *pshell = (egadsShell *) children[i]->blind;
        builder3D.Add(solid, pshell->shell);
      }
      try {
        BRepLib::OrientClosedSolid(solid);
      }
      catch (Standard_Failure)
      {
        printf(" EGADS Warning: Cannot Orient Solid (EG_makeTopology)!\n");
        Handle_Standard_Failure e = Standard_Failure::Caught();
        printf("                %s\n", e->GetMessageString());
        return EGADS_TOPOERR;
      }
      catch (...) {
        printf(" EGADS Warning: Cannot Orient Solid (EG_makeTopology)!\n");
        return EGADS_TOPOERR;
      }
      BRepCheck_Analyzer sCheck(solid);
      if (!sCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(solid);
        sfs->Perform();
        TopoDS_Shape fixedSolid = sfs->Shape();
        if (fixedSolid.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Solid is invalid (EG_makeTopology)!\n");
          solid.Nullify();
          return EGADS_CONSTERR;
        } else {
          if (outLevel > 0)
            printf(" EGADS Warning: Fixing Solid (EG_makeTopology)!\n");
          solid = TopoDS::Solid(fixedSolid);
          fixedSld = 1;
        }
      }
      shape = solid;
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Body object (EG_makeTopology)!\n");
      return stat;
    }
    obj->oclass        = oclass;
    obj->mtype         = mtype;
    egadsBody *pbody   = new egadsBody;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pbody->shape       = shape;
    obj->blind         = pbody;
    stat = EG_traverseBody(context, 0, obj, obj, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      delete pbody;
      return stat;
    }
    for (i = 0; i < nChildren; i++)
      EG_attriBodyDup(children[i], obj);
    /* patch up Face Attributes for fixed solid */
    if (fixedSld == 1) {
      for (i = 0; i < nChildren; i++) {
        egadsShell *pshell = (egadsShell *) children[i]->blind;
        for (int j = 0; j < pshell->nfaces; j++) {
          egadsFace *pface = (egadsFace *) pshell->faces[j]->blind;
          int index = pbody->faces.map.FindIndex(pface->face);
          if (index != 0) continue;
          for (int k = 0; k < pbody->faces.map.Extent(); k++) {
            stat = EG_isSame(pshell->faces[j], pbody->faces.objs[k]);
            if (stat != EGADS_SUCCESS) continue;
            EG_attributeDup(pshell->faces[j], pbody->faces.objs[k]);
            break;
          }
        }
      }
    }
    EG_referenceObject(obj, context);

  } else {

    if (children == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Model with NULL Children (EG_makeTopology)!\n");
      return EGADS_NULLOBJ;
    }
    if (nChildren <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Model with %d Bodies (EG_makeTopology)!\n",
               nChildren);
      return EGADS_RANGERR;
    }
    for (i = 0; i < nChildren; i++) {
      if (children[i] == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Model with Body[%d] NULL (EG_makeTopology)!\n",
                 i);
        return EGADS_NULLOBJ;
      }
      if (children[i]->oclass != BODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ nonBody Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NOTTOPO;
      }
      if (children[i]->topObj != context) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ body[%d] reference (EG_makeTopology)!\n",
                 i);
        return EGADS_REFERCE;
      }
      if (children[i]->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Model w/ NULL Body Child[%d] (EG_makeTopology)!\n",
                 i);
        return EGADS_NODATA;
      }
    }
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Model object (EG_makeTopology)!\n");
      return stat;
    }
    egadsModel *pmodel = new egadsModel;
    pmodel->bodies     = new egObject*[nChildren];
    pmodel->nbody      = nChildren;
    obj->oclass = MODEL;
    obj->blind  = pmodel;
    TopoDS_Compound compound;
    BRep_Builder    builder3D;
    builder3D.MakeCompound(compound);
    for (i = 0; i < nChildren; i++) {
      pmodel->bodies[i] = children[i];
      egadsBody *pbody  = (egadsBody *) children[i]->blind;
      builder3D.Add(compound, pbody->shape);
      EG_referenceObject(children[i], obj);
      EG_removeCntxtRef(children[i]);
    }
    pmodel->shape = compound;
    EG_referenceObject(obj, context);

  }

  *topo = obj;
  return EGADS_SUCCESS;
}


int
EG_makeLoop(int nedge, egObject **edges, /*@null@*/ egObject *surf,
            double toler, egObject **result)
{
  int      i, first, stat, outLevel, oclass, mtype, mtypc, close, nn, nc, nnew;
  int      hit, *lsenses, *senses;
  double   tol, nodetol, tlim[2], xyz0[3], xyz1[3], xyz2[3], xyz3[3];
  egObject *context, *ref, *curv, **children, **nodes, *fnode, *nnode, *geom;
  egObject **newedges, *nds[2];

  *result = NULL;
  geom    = surf;
  if (nedge < 1)                          return EGADS_EMPTY;
  if (edges == NULL)                      return EGADS_NULLOBJ;
  for (first = 0; first < nedge; first++)
    if (edges[first] != NULL) break;
  if (first == nedge)                     return EGADS_NULLOBJ;
  if (edges[first]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (surf != NULL) {
    if (surf->oclass != SURFACE)          return EGADS_NOTGEOM;
    if (surf->mtype  == PLANE) geom = NULL;
  }
  outLevel = EG_outLevel(edges[first]);
  context  = EG_context(edges[first]);

  tol = nodetol = 0.0;
  for (i = 0; i < nedge; i++) {
    if (edges[i] == NULL) continue;
    if (edges[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not an EGO (EG_makeLoop)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (edges[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d has no data (EG_makeLoop)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(edges[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d Context Mismatch (EG_makeLoop)!\n",
               i+1);
      return EGADS_MIXCNTX;
    }
    if (edges[i]->oclass != EDGE) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is Not an EDGE (EG_makeLoop)!\n", i+1);
      return EGADS_NOTTOPO;
    }
    egadsEdge *pedge = (egadsEdge *) edges[i]->blind;
    double toltmp = BRep_Tool::Tolerance(pedge->edge);
    if (toltmp > tol) tol = toltmp;
    BRep_Tool::Range(pedge->edge, tlim[0], tlim[1]);
    egObject  *curvo  = pedge->curve;
    if (curvo  == NULL) continue;
    egadsCurve *pcurve = (egadsCurve *) curvo->blind;
    if (pcurve == NULL) continue;
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAdaptor_Curve AC(hCurve);
    toltmp = GCPnts_AbscissaPoint::Length(AC, tlim[0], tlim[1]);
    if (nodetol == 0.0) {
      nodetol = toltmp;
    } else {
      if (nodetol > toltmp) nodetol = toltmp;
    }
  }
  nodetol *= 0.25;
  if (nodetol > tol) nodetol = tol;
  if (toler   > tol) tol     = toler;
  lsenses = (int *) EG_alloc(nedge*sizeof(int));
  if (lsenses == NULL) return EGADS_MALLOC;
  newedges = (egObject **) EG_alloc(2*nedge*sizeof(egObject *));
  if (newedges == NULL) {
    EG_free(lsenses);
    return EGADS_MALLOC;
  }
  if (outLevel > 1)
    printf(" EG_makeLoop: Nedge = %d  Tolerance = %le (%le)  nodetol = %le\n",
           nedge, tol, toler, nodetol);

  /* set the first edge */
  xyz0[0] = xyz0[1] = xyz0[2] = 0.0;
  stat = EG_getTopology(edges[first], &ref, &oclass, &mtype, tlim, &nn, &nodes,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    EG_free(newedges);
    EG_free(lsenses);
    return stat;
  }
  stat = EG_getTopology(nodes[0], &ref, &oclass, &mtype, xyz0, &nc, &children,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    EG_free(newedges);
    EG_free(lsenses);
    return stat;
  }
  xyz1[0] = xyz0[0];
  xyz1[1] = xyz0[1];
  xyz1[2] = xyz0[2];
  fnode   = nodes[0];
  if (nn != 1) {
    nnode = nodes[1];
    stat  = EG_getTopology(nodes[1], &ref, &oclass, &mtype, xyz1, &nc, &children,
                           &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(newedges);
      EG_free(lsenses);
      return stat;
    }
  }
  newedges[0]  = edges[first];
  lsenses[0]   = 1;
  edges[first] = NULL;
  nnew         = 1;
  if (nn == 1) {
    if (geom != NULL) {
      stat = EG_otherCurve(geom, ref, tol, &newedges[1]);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
    }
    /* make the loop */
    stat = EG_makeTopology(context, geom, LOOP, CLOSED, NULL, nnew, newedges,
                           lsenses, result);
    EG_free(newedges);
    EG_free(lsenses);
    if (stat == EGADS_SUCCESS)
      for (i = 0; i < nedge; i++)
        if (edges[i] != NULL) stat++;
    return stat;
  }
  close = OPEN;

  /* serach for next edges (forward) */
  do {
    for (hit = i = 0; i < nedge; i++) {
      if (edges[i] == NULL) continue;
      stat = EG_getTopology(edges[i], &curv, &oclass, &mtypc, tlim, &nn, &nodes,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (nn == 1) continue;
      stat = EG_getTopology(nodes[0], &ref, &oclass, &mtype, xyz2, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      stat = EG_getTopology(nodes[1], &ref, &oclass, &mtype, xyz3, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) +
               (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) +
               (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2])) <= nodetol) {
        lsenses[nnew] =  1;
        xyz1[0] = xyz3[0];
        xyz1[1] = xyz3[1];
        xyz1[2] = xyz3[2];
        nds[0]  = nnode;
        nds[1]  = nnode = nodes[1];
      } else if (sqrt((xyz1[0]-xyz3[0])*(xyz1[0]-xyz3[0]) +
                      (xyz1[1]-xyz3[1])*(xyz1[1]-xyz3[1]) +
                      (xyz1[2]-xyz3[2])*(xyz1[2]-xyz3[2])) <= nodetol) {
        lsenses[nnew] = -1;
        xyz1[0] = xyz2[0];
        xyz1[1] = xyz2[1];
        xyz1[2] = xyz2[2];
        nds[1]  = nnode;
        nds[0]  = nnode = nodes[0];
      } else {
        continue;
      }
      hit++;
      /* are we closed -- xyz1 == xyz0? */
      if (sqrt((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]) +
               (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1]) +
               (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2])) <= nodetol) {
        close  = CLOSED;
        if (lsenses[nnew] == 1) {
          nds[1] = fnode;
        } else {
          nds[0] = fnode;
        }
      }
      if ((nodes[0] == nds[0]) && (nodes[1] == nds[1])) {
        newedges[nnew] = edges[i];
      } else {
        if (outLevel > 1)
          printf(" EG_makeLoop: New Edge for %d\n", i);
        stat = EG_makeTopology(context, curv, EDGE, mtypc, tlim, 2, nds,
                               senses, &newedges[nnew]);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
      }
      edges[i] = NULL;
      nnew++;
    }
  } while ((close == OPEN) && (hit != 0));

  if (close == CLOSED) {
    if (geom != NULL) {
      for (i = 0; i < nnew; i++) {
        stat = EG_getTopology(newedges[i], &curv, &oclass, &mtypc, tlim, &nn,
                              &nodes, &senses);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
        stat = EG_otherCurve(geom, curv, tol, &newedges[nnew+i]);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
      }
    }
    /* make the loop */
    stat = EG_makeTopology(context, geom, LOOP, CLOSED, NULL, nnew, newedges,
                           lsenses, result);
    EG_free(newedges);
    EG_free(lsenses);
    if (stat == EGADS_SUCCESS)
      for (i = 0; i < nedge; i++)
        if (edges[i] != NULL) stat++;
    return stat;
  }

  /* start from the first and look back -- will be open */
  for (i = 0; i < nnew/2; i++) {
    nn           =  nnew-i-1;
    ref          =  newedges[i];
    newedges[i]  =  newedges[nn];
    newedges[nn] =  ref;
    nc           =  lsenses[i];
    lsenses[i]   = -lsenses[nn];
    lsenses[nn]  = -nc;
  }
  xyz1[0] = xyz0[0];
  xyz1[1] = xyz0[1];
  xyz1[2] = xyz0[2];
  nnode   = fnode;

  do {
    for (hit = i = 0; i < nedge; i++) {
      if (edges[i] == NULL) continue;
      stat = EG_getTopology(edges[i], &curv, &oclass, &mtypc, tlim, &nn, &nodes,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (nn == 1) continue;
      stat = EG_getTopology(nodes[0], &ref, &oclass, &mtype, xyz2, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      stat = EG_getTopology(nodes[1], &ref, &oclass, &mtype, xyz3, &nc,
                            &children, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      if (sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) +
               (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) +
               (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2])) <= nodetol) {
        lsenses[nnew] =  1;
        xyz1[0] = xyz3[0];
        xyz1[1] = xyz3[1];
        xyz1[2] = xyz3[2];
        nds[0]  = nnode;
        nds[1]  = nnode = nodes[1];
      } else if (sqrt((xyz1[0]-xyz3[0])*(xyz1[0]-xyz3[0]) +
                      (xyz1[1]-xyz3[1])*(xyz1[1]-xyz3[1]) +
                      (xyz1[2]-xyz3[2])*(xyz1[2]-xyz3[2])) <= nodetol) {
        lsenses[nnew] = -1;
        xyz1[0] = xyz2[0];
        xyz1[1] = xyz2[1];
        xyz1[2] = xyz2[2];
        nds[1]  = nnode;
        nds[0]  = nnode = nodes[0];
      } else {
        continue;
      }
      hit++;
      if ((nodes[0] == nds[0]) && (nodes[1] == nds[1])) {
        newedges[nnew] = edges[i];
      } else {
        if (outLevel > 1)
          printf(" EG_makeLoop: New Edge for %d\n", i);
        stat = EG_makeTopology(context, curv, EDGE, mtypc, tlim, 2, nds,
                               senses, &newedges[nnew]);
        if (stat != EGADS_SUCCESS) {
          EG_free(newedges);
          EG_free(lsenses);
          return stat;
        }
      }
      edges[i] = NULL;
      nnew++;
    }
  } while (hit != 0);

  if (geom != NULL) {
    for (i = 0; i < nnew; i++) {
      stat = EG_getTopology(newedges[i], &curv, &oclass, &mtypc, tlim, &nn,
                            &nodes, &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
      stat = EG_otherCurve(geom, curv, tol, &newedges[nnew+i]);
      if (stat != EGADS_SUCCESS) {
        EG_free(newedges);
        EG_free(lsenses);
        return stat;
      }
    }
  }
  /* make the loop */
  stat = EG_makeTopology(context, geom, LOOP, OPEN, NULL, nnew, newedges,
                         lsenses, result);
  EG_free(newedges);
  EG_free(lsenses);
  if (stat == EGADS_SUCCESS)
    for (i = 0; i < nedge; i++)
      if (edges[i] != NULL) stat++;
  return stat;
}


int
EG_getPlane(const egObject *object, egObject **plane)
{
  *plane = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass != LOOP)       return EGADS_GEOMERR;
  if (object->blind == NULL)        return EGADS_NODATA;
  egadsLoop *ploop  = (egadsLoop *) object->blind;
  if (ploop->surface != NULL)       return EGADS_GEOMERR;
  int outLevel = EG_outLevel(object);
  egObject *context = EG_context(object);

  // try to fit a plane
  Standard_Real tol = Precision::Confusion();
  Handle(Geom_Surface) hSurface;
  int nTrys = 4;              // # of tol attempts
  for (int itry = 0; itry < nTrys; itry++) {
    BRepBuilderAPI_FindPlane FPlane(ploop->loop, tol);
    if (FPlane.Found()) {
      hSurface = FPlane.Plane();
      break;
    }
    tol *= 10.0;
    if (outLevel > 1)
      printf(" EGADS Info: Adjusting Prec for getPlane - itry = %d  prec = %le\n",
             itry, tol);
  }
  if (hSurface.IsNull()) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Planar Surface (EG_getPlane)!\n");
    return EGADS_GEOMERR;
  }

  egObject *obj;
  int stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: make Surface = %d (EG_getPlane)!\n", stat);
    return stat;
  }
  obj->oclass         = SURFACE;
  obj->mtype          = PLANE;
  egadsSurface *psurf = new egadsSurface;
  psurf->handle       = hSurface;
  psurf->basis        = NULL;
  psurf->topFlg       = 1;
  obj->blind          = psurf;
  EG_referenceObject(obj, context);

  *plane = obj;
  return EGADS_SUCCESS;
}


int
EG_getArea(egObject *object, /*@null@*/ const double *limits,
           double *area)
{
  int         stat;
  double      sense = 1.0;
  TopoDS_Face Face;

  *area = 0.0;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != SURFACE) && (object->oclass != LOOP) &&
      (object->oclass != FACE))     return EGADS_GEOMERR;
  if (object->blind == NULL)        return EGADS_NODATA;
  int outLevel = EG_outLevel(object);

  if (object->oclass == FACE) {

    egadsFace *pface = (egadsFace *) object->blind;
    Face = pface->face;

  } else if (object->oclass == SURFACE) {

    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface with NULL Limits (EG_getArea)!\n");
      return EGADS_NODATA;
    }
    egadsSurface *psurf = (egadsSurface *) object->blind;
#if CASVER >= 652
    Standard_Real tol = BRepLib::Precision();
    BRepLib_MakeFace MFace(psurf->handle, limits[0], limits[1],
                                          limits[2], limits[3], tol);
#else
    BRepLib_MakeFace MFace(psurf->handle, limits[0], limits[1],
                                          limits[2], limits[3]);
#endif
    Face = MFace.Face();

  } else {

    egadsLoop *ploop  = (egadsLoop *) object->blind;
    if (ploop->surface == NULL) {

      // try to fit a plane
      Standard_Real tol = Precision::Confusion();
      Handle(Geom_Surface) hSurface;
      int nTrys = 4;              // # of tol attempts
      for (int itry = 0; itry < nTrys; itry++) {
        BRepBuilderAPI_FindPlane FPlane(ploop->loop, tol);
        if (FPlane.Found()) {
          hSurface = FPlane.Plane();
          break;
        }
        tol *= 10.0;
        if (outLevel > 1)
          printf(" EGADS Info: Adjusting Prec for makeFace - itry = %d  prec = %le\n",
                 itry, tol);
      }
      if (hSurface.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Planar Surface (EG_getArea)!\n");
        return EGADS_GEOMERR;
      }

      signal(SIGSEGV, segfault_handler);
      switch (stat = setjmp(jmpenv)) {
      case 0:
        {
          try {
            BRepLib_MakeFace MFace(hSurface, ploop->loop);
            Face = MFace.Face();
          }
          catch (Standard_Failure) {
            printf(" EGADS Warning: Cannot makeFace (EG_getArea)!\n");
            Handle_Standard_Failure e = Standard_Failure::Caught();
            printf("                %s\n", e->GetMessageString());
            return EGADS_CONSTERR;
          }
          catch (...) {
            printf(" EGADS Warning: Cannot makeFace (EG_getArea)!\n");
            return EGADS_CONSTERR;
          }
        }
        break;
      default:
        printf(" EGADS Fatal Error: OCC SegFault %d (EG_getArea)!\n", stat);
        signal(SIGSEGV, SIG_DFL);
        return EGADS_OCSEGFLT;
      }
      signal(SIGSEGV, SIG_DFL);

      // did making the Face flip the Loop?
      TopExp_Explorer ExpW;
      for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
        TopoDS_Shape shapw = ExpW.Current();
        TopoDS_Wire  wire  = TopoDS::Wire(shapw);
        if (!wire.IsSame(ploop->loop)) continue;
        if (!wire.IsEqual(ploop->loop)) sense = -1.0;
      }

    } else {

      // make a standard Face
      egObject *geom = ploop->surface;
      if (geom->blind == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop had NULL Ref Surface (EG_getArea)!\n");
        return EGADS_NOTGEOM;
      }
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      BRepBuilderAPI_MakeFace MFace;
      Standard_Real old = BRepBuilderAPI::Precision();
      Standard_Real prc = old;
      int nTrys = 5;              // # of tol attempts before setting SameParam = F
      for (int itry = 0; itry < nTrys; itry++) {
        if (itry != 0) Face.Nullify();
        if (itry == nTrys-1) {    // last attempt -- set tol back
          prc = old;
          BRepBuilderAPI::Precision(old);
        }
#if CASVER >= 652
        MFace.Init(psurf->handle, Standard_False, prc);
#else
        MFace.Init(psurf->handle, Standard_False);
#endif
        MFace.Add(ploop->loop);
        if (MFace.Error()) {
          if (outLevel > 0)
            printf(" EGADS Error: Problem with Loop (EG_getArea)!\n");
          return EGADS_NODATA;
        }
        if (MFace.IsDone()) {
          Face = MFace.Face();
          EG_makePCurves(Face, ploop->surface, object, prc, nTrys-itry-1);
          BRepLib::SameParameter(Face);
          BRepCheck_Analyzer oCheck(Face);
          if (oCheck.IsValid()) break;
        }
        prc *= 10.0;
        BRepBuilderAPI::Precision(prc);
        if (outLevel > 1)
          printf(" EGADS Info: Adjusting Precision for Face - itry = %d  prec = %lf\n",
                 itry, prc);
      }
      BRepBuilderAPI::Precision(old);
      if (!MFace.IsDone()) {
        if (outLevel > 0)
          printf(" EGADS Error: Problem making the Face (EG_getArea)!\n");
        return EGADS_NODATA;
      }

    }

  }

  BRepGProp    BProps;
  GProp_GProps SProps;
  BProps.SurfaceProperties(Face, SProps);
  *area = sense*SProps.Mass();

  return EGADS_SUCCESS;
}


int
EG_makeFace(egObject *object, int mtype,
            /*@null@*/ const double *limits, egObject **face)
{
  int      stat, outLevel, nl = 1;
  egObject *obj, *context, *loop = NULL, *geom = NULL;

  *face = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass != SURFACE) && (object->oclass != LOOP) &&
      (object->oclass != FACE))     return EGADS_GEOMERR;
  if (object->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(object);
  context  = EG_context(object);

  if (object->oclass != FACE)
    if ((mtype != SFORWARD) && (mtype != SREVERSE)) {
      if (outLevel > 0)
        printf(" EGADS Error: Mtype = %d (EG_makeFace)!\n", mtype);
      return EGADS_TOPOERR;
    }
  if (object->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) object->blind;
    if (ploop->surface != NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Loop had Ref Surface (EG_makeFace)!\n");
      return EGADS_NOTGEOM;
    }
  } else {
    if (limits == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Surface/Face with NULL Limits (EG_makeFace)!\n");
      return EGADS_NODATA;
    }
  }

  TopoDS_Face Face;
  if (object->oclass == SURFACE) {

    egadsSurface *psurf = (egadsSurface *) object->blind;
    Handle(Geom_Surface) hSurf = psurf->handle;
#if CASVER >= 652
    Standard_Real tol = BRepLib::Precision();
    BRepLib_MakeFace MFace(hSurf, limits[0], limits[1],
                                  limits[2], limits[3], tol);
#else
    BRepLib_MakeFace MFace(hSurf, limits[0], limits[1],
                                  limits[2], limits[3]);
#endif
    Face = MFace.Face();
    if (mtype == SREVERSE) {
      Face.Orientation(TopAbs_REVERSED);
    } else {
      Face.Orientation(TopAbs_FORWARD);
    }
    BRepLib::SameParameter(Face);
    BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
      EG_checkStatus(fCheck.Result(Face));
      return EGADS_CONSTERR;
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeFace)!\n");
      return stat;
    }
    obj->oclass = FACE;
    geom        = object;
    TopExp_Explorer ExpW;
    for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  Wire  = TopoDS::Wire(shapw);
      stat = EG_makeObject(context, &loop);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Surface object (EG_makeFace)!\n");
        obj->oclass = NIL;
        EG_deleteObject(obj);
        return stat;
      }
      loop->oclass      = LOOP;
      egadsLoop *ploop  = new egadsLoop;
      loop->blind       = ploop;
      ploop->loop       = Wire;
      ploop->nedges     = 0;
      ploop->edges      = NULL;
      ploop->senses     = NULL;
      ploop->topFlg     = 0;
      ploop->surface    = NULL;
      if (object->mtype != PLANE) {
        ploop->surface = geom;
        EG_referenceObject(geom, loop);
      }
      EG_fillTopoObjs(loop, obj);
      EG_fillPCurves(Face, geom, loop, obj);
      break;
    }

  } else if (object->oclass == FACE) {

    egadsFace *pfasrc = (egadsFace *) object->blind;
    geom = pfasrc->surface;
    if (geom->mtype != PLANE) {
      if (outLevel > 0)
        printf(" EGADS Error: Face is NOT Planar (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }
    if (pfasrc->nloops != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face has multiple Loops (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }
    BRepOffsetAPI_MakeOffset mOffset(pfasrc->face);
    mOffset.Perform(-limits[0]);
    mOffset.Build();
    TopoDS_Shape shapf = mOffset.Shape();

    int nWire = 0;
    TopExp_Explorer Exp;
    if (shapf.ShapeType() == TopAbs_WIRE) {
      TopoDS_Wire wire = TopoDS::Wire(shapf);
      BRepCheck_Analyzer sCheck(wire);
      if (!sCheck.IsValid()) {
        if (outLevel > 0)
          printf(" EGADS Info: Offset Wire is NOT valid (EG_makeFace)!\n");
        return EGADS_CONSTERR;
      }
      nWire++;
//    printf(" Is a Wire!\n");
    } else if (shapf.ShapeType() == TopAbs_COMPOUND) {
      for (Exp.Init(shapf, TopAbs_WIRE); Exp.More(); Exp.Next()) {
        TopoDS_Wire wire = TopoDS::Wire(Exp.Current());
        BRepCheck_Analyzer wCheck(wire);
        if (!wCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Info: Offset Wire is NOT valid (EG_makeFace)!\n");
          return EGADS_CONSTERR;
        }
        nWire++;
      }
//    printf(" Is a Compound -- with %d Wires!\n", nWire);
    } else {
      printf(" EGADS Info: What is this shape (EG_makeFace)?!\n");
    }
    if (nWire == 0) {                             // not really an error
      *face = object;                             //  just no hole!
      return EGADS_OUTSIDE;
    }
    TopoDS_Wire *wires = new TopoDS_Wire[nWire];
    if (shapf.ShapeType() == TopAbs_WIRE) {
      wires[0] = TopoDS::Wire(shapf);
    } else {
      int i = 0;
      for (Exp.Init(shapf, TopAbs_WIRE); Exp.More(); Exp.Next()) {
        wires[i] = TopoDS::Wire(Exp.Current());
        i++;
      }
    }

    // fillet the new Wires (holes)?
    if (limits[1] > 0.0) {
      BRepFilletAPI_MakeFillet2d f2d;
      egadsSurface *psurf = (egadsSurface *) geom->blind;
      for (int i = 0; i < nWire; i++) {
        BRepLib_MakeFace MFace(psurf->handle, wires[i]);
        Face = MFace.Face();
        f2d.Init(Face);
        TopTools_IndexedMapOfShape vmap;
        TopExp::MapShapes(wires[i], TopAbs_VERTEX, vmap);
        if (vmap.Extent() == 0) continue;
        for (int j = 0; j < vmap.Extent(); j++) {
          TopoDS_Shape  shapv = vmap(j+1);
          TopoDS_Vertex vert  = TopoDS::Vertex(shapv);
          f2d.AddFillet(vert, limits[1]);
/*
          ChFi2d_ConstructionError err = f2d.Status();
          if (err == ChFi2d_NotPlanar)            printf(" status = NotPlanar!\n");
          if (err == ChFi2d_NoFace)               printf(" status = NoFace!\n");
          if (err == ChFi2d_InitialisationError)  printf(" status = NoInit!\n");
          if (err == ChFi2d_ParametersError)      printf(" status = ParamErr!\n");
          if (err == ChFi2d_Ready)                printf(" status = Ready!\n");
          if (err == ChFi2d_IsDone)               printf(" status = IsDone!\n");
          if (err == ChFi2d_ComputationError)     printf(" status = ComputErr!\n");
          if (err == ChFi2d_ConnexionError)       printf(" status = Connexion!\n");
          if (err == ChFi2d_TangencyError)        printf(" status = TangenErr!\n");
          if (err == ChFi2d_FirstEdgeDegenerated) printf(" status = 1Degen!\n");
          if (err == ChFi2d_LastEdgeDegenerated)  printf(" status = LastDegen!\n");
          if (err == ChFi2d_BothEdgesDegenerated) printf(" status = BothDegen!\n");
          if (err == ChFi2d_NotAuthorized)        printf(" status = NotAuthor!\n");
*/
        }
        try {
          TopoDS_Shape shape = f2d.Shape();
//        if (shape.ShapeType() == TopAbs_FACE) printf(" filletted shape: Face!\n");
          Face = TopoDS::Face(shape);
        }
        catch (Standard_Failure) {
          printf(" EGADS Error: Cannot get Filletted Face (EG_makeFace)!\n");
          Handle_Standard_Failure e = Standard_Failure::Caught();
          printf("                %s\n", e->GetMessageString());
          delete [] wires;
          return EGADS_CONSTERR;
        }
        catch (...) {
          printf(" EGADS Error: Cannot get Filletted Face (EG_makeFace)!\n");
          delete [] wires;
          return EGADS_CONSTERR;
        }

        // update the wires with the filletted version
        TopTools_IndexedMapOfShape wmap;
        TopExp::MapShapes(Face, TopAbs_WIRE, wmap);
        if (wmap.Extent() != 1)
          if (outLevel > 0)
            printf(" EGADS Warning: Fillet %d #Loops = %d (EG_makeFace)!\n",
                   i+1, wmap.Extent());
        wires[i] = TopoDS::Wire(wmap(1));
      }
    }

    // make the new Face
    BRepBuilderAPI_MakeFace mFace(pfasrc->face);
    for (int i = 0; i < nWire; i++) {
      if (object->mtype == SFORWARD) wires[i].Reverse();   // reverse holes
      mFace.Add(wires[i]);
    }
    if (!mFace.IsDone()) {
      delete [] wires;
      printf(" EGADS Warning: BRepBuilderAPI_MakeFace Error (EG_imakeFace)!\n");
      return EGADS_GEOMERR;
    }
    Face = mFace.Face();
    BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      delete [] wires;
      if (outLevel > 0)
        printf(" EGADS Info: Offset Face is NOT valid (EG_makeFace)!\n");
      return EGADS_CONSTERR;
    }

    // make the EGADS objects
    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      delete [] wires;
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeFace)!\n");
      return stat;
    }
    obj->oclass      = FACE;
    egObject **loopo = new egObject*[nWire+1];
    int *senses      = new int[nWire+1];
    loopo[0]         = pfasrc->loops[0];
    senses[0]        = 1;
    for (int i = 0; i < nWire; i++) {
      stat = EG_makeObject(context, &loop);
      if (stat != EGADS_SUCCESS) {
        delete [] loopo;
        delete [] senses;
        delete [] wires;
        if (outLevel > 0)
          printf(" EGADS Error: Cannot make Loop object (EG_makeFace)!\n");
        obj->oclass = NIL;
        EG_deleteObject(obj);
        return stat;
      }
      loop->oclass     = LOOP;
      egadsLoop *ploop = new egadsLoop;
      loop->blind      = ploop;
      ploop->loop      = wires[i];
      ploop->nedges    = 0;
      ploop->edges     = NULL;
      ploop->senses    = NULL;
      ploop->topFlg    = 0;
      ploop->surface   = NULL;
      EG_fillTopoObjs(loop, obj);
      EG_fillPCurves(Face, geom, loop, obj);
      loopo[i+1]       = loop;
      senses[i+1]      = -1;
    }
    delete [] wires;

    egadsFace *pface = new egadsFace;
    pface->face      = Face;
    pface->surface   = geom;
    pface->nloops    = nWire+1;
    pface->loops     = loopo;
    pface->senses    = senses;
    pface->topFlg    = 0;
    obj->blind       = pface;
    obj->mtype       = object->mtype;

    EG_referenceObject(geom, obj);
    for (int i = 0; i <= nWire; i++) EG_referenceObject(loopo[i], obj);
    EG_referenceObject(obj,  context);
    EG_attributeDup(object,  obj);

    *face = obj;
    return EGADS_SUCCESS;

  } else {

    egadsLoop *ploop  = (egadsLoop *) object->blind;
    Standard_Real tol = Precision::Confusion();
    Handle(Geom_Surface) hSurface;
    int nTrys = 4;              // # of tol attempts
    for (int itry = 0; itry < nTrys; itry++) {
      BRepBuilderAPI_FindPlane FPlane(ploop->loop, tol);
      if (FPlane.Found()) {
        hSurface = FPlane.Plane();
        break;
      }
      tol *= 10.0;
      if (outLevel > 1)
        printf(" EGADS Info: Adjusting Prec for makeFace - itry = %d  prec = %le\n",
               itry, tol);
    }
    if (hSurface.IsNull()) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Planar Surface (EG_makeFace)!\n");
      return EGADS_GEOMERR;
    }

    signal(SIGSEGV, segfault_handler);
    switch (stat = setjmp(jmpenv)) {
    case 0:
      {
        try {
          BRepLib_MakeFace MFace(hSurface, ploop->loop);
          Face = MFace.Face();
        }
        catch (Standard_Failure) {
          printf(" EGADS Warning: Cannot makeFace (EG_makeFace)!\n");
          Handle_Standard_Failure e = Standard_Failure::Caught();
          printf("                %s\n", e->GetMessageString());
          return EGADS_CONSTERR;
        }
        catch (...) {
          printf(" EGADS Warning: Cannot makeFace (EG_makeFace)!\n");
          return EGADS_CONSTERR;
        }
      }
      break;
    default:
      printf(" EGADS Fatal Error: OCC SegFault %d (EG_makeFace)!\n", stat);
      signal(SIGSEGV, SIG_DFL);
      return EGADS_OCSEGFLT;
    }
    signal(SIGSEGV, SIG_DFL);

    if (mtype == SREVERSE) {
      Face.Orientation(TopAbs_REVERSED);
    } else {
      Face.Orientation(TopAbs_FORWARD);
    }
    BRepLib::SameParameter(Face);
/*  BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      if (outLevel > 0)
        printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
      return EGADS_CONSTERR;
    }  */
    BRepCheck_Wire wireCheck(ploop->loop);
    wireCheck.InContext(Face);
    BRepCheck_ListOfStatus stati;
    stati = wireCheck.Status();
    if (stati.Extent() > 0) {
      BRepCheck_Status cStatus = stati.First();
      if (outLevel > 0)
        EG_printStatus(cStatus);
      if (cStatus != BRepCheck_NoError) return EGADS_GEOMERR;
    }

    BRepCheck_Analyzer fCheck(Face);
    if (!fCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(Face);
      sfs->Perform();
      TopoDS_Shape fixedFace = sfs->Shape();
      if (fixedFace.IsNull()) {
        if (outLevel > 0)
          printf(" EGADS Info: Fixing Face Failed (EG_makeFace)!\n");
        return EGADS_CONSTERR;
      } else {
        BRepCheck_Analyzer sfCheck(fixedFace);
        if (!sfCheck.IsValid()) {
          if (outLevel > 0)
            printf(" EGADS Info: Face may be invalid (EG_makeFace)!\n");
          return EGADS_CONSTERR;
        } else {
          Face = TopoDS::Face(fixedFace);
        }
      }
    }

    stat = EG_makeObject(context, &obj);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Face object (EG_makeFace)!\n");
      return stat;
    }
    obj->oclass = FACE;

    loop = object;
    stat = EG_makeObject(context, &geom);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Surface object (EG_makeFace)!\n");
      obj->oclass = NIL;
      EG_deleteObject(obj);
      return stat;
    }
    geom->topObj = obj;
    EG_completeSurf(geom, hSurface);

  }

  egadsFace *pface = new egadsFace;
  egObject **loopo = new egObject*[nl];
  int *senses      = new int[nl];
  loopo[0]         = loop;
  senses[0]        = 1;
  pface->face      = Face;
  pface->surface   = geom;
  pface->nloops    = nl;
  pface->loops     = loopo;
  pface->senses    = senses;
  pface->topFlg    = 0;
  obj->blind       = pface;
  obj->mtype       = mtype;

  EG_referenceObject(geom, obj);
  EG_referenceObject(loop, obj);
  EG_referenceObject(obj,  context);

  *face = obj;
  return EGADS_SUCCESS;
}


int
EG_getBodyTopos(const egObject *body, /*@null@*/ egObject *src,
                int oclass, int *ntopo, egObject ***topos)
{
  int      outLevel, i, n, index;
  egadsMap *map;
  egObject **objs;

  *ntopo = 0;
  *topos = NULL;
  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(body);

  if ((oclass < NODE) || (oclass > SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: oclass = %d (EG_getBodyTopos)!\n",
             oclass);
    return EGADS_NOTTOPO;
  }

  egadsBody *pbody = (egadsBody *) body->blind;
  if (oclass == NODE) {
    map = &pbody->nodes;
  } else if (oclass == EDGE) {
    map = &pbody->edges;
  } else if (oclass == LOOP) {
    map = &pbody->loops;
  } else if (oclass == FACE) {
    map = &pbody->faces;
  } else {
    map = &pbody->shells;
  }

  if (src == NULL) {

    n = map->map.Extent();
    if (n == 0) return EGADS_SUCCESS;
    objs = (egObject **) EG_alloc(n*sizeof(egObject *));
    if (objs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc oclass = %d, n = %d (EG_getBodyTopos)!\n",
               oclass, n);
      return EGADS_MALLOC;
    }
    for (i = 0; i < n; i++) objs[i] = map->objs[i];

  } else {

    if (src->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: src not an EGO (EG_getBodyTopos)!\n");
      return EGADS_NOTOBJ;
    }
    if ((src->oclass < NODE) || (src->oclass > SHELL)) {
      if (outLevel > 0)
        printf(" EGADS Error: src not a Topo (EG_getBodyTopos)!\n");
      return EGADS_NOTTOPO;
    }
    if (src->oclass == oclass) {
      if (outLevel > 0)
        printf(" EGADS Error: src Topo is oclass (EG_getBodyTopos)!\n");
      return EGADS_TOPOERR;
    }
    if (EG_context(body) != EG_context(src)) {
      if (outLevel > 0)
        printf(" EGADS Error: Context mismatch (EG_getBodyTopos)!\n");
      return EGADS_MIXCNTX;
    }
    if (src->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL src pointer (EG_getBodyTopos)!\n");
      return EGADS_NODATA;
    }

    TopoDS_Shape     shape;
    TopAbs_ShapeEnum senum;
    if (src->oclass == NODE) {
      egadsNode *pnode = (egadsNode *) src->blind;
      shape = pnode->node;
      senum = TopAbs_VERTEX;
    } else if (src->oclass == EDGE) {
      egadsEdge *pedge = (egadsEdge *) src->blind;
      shape = pedge->edge;
      senum = TopAbs_EDGE;
    } else if (src->oclass == LOOP) {
      egadsLoop *ploop = (egadsLoop *) src->blind;
      shape = ploop->loop;
      senum = TopAbs_WIRE;
    } else if (src->oclass == FACE) {
      egadsFace *pface = (egadsFace *) src->blind;
      shape = pface->face;
      senum = TopAbs_FACE;
    } else {
      egadsShell *pshell = (egadsShell *) src->blind;
      shape = pshell->shell;
      senum = TopAbs_SHELL;
    }

    if (src->oclass > oclass) {

      // look down the tree (get sub-shapes)
      if (oclass == NODE) {
        senum = TopAbs_VERTEX;
      } else if (oclass == EDGE) {
        senum = TopAbs_EDGE;
      } else if (oclass == LOOP) {
        senum = TopAbs_WIRE;
      } else if (oclass == FACE) {
        senum = TopAbs_FACE;
      } else {
        senum = TopAbs_SHELL;
      }
      TopTools_IndexedMapOfShape smap;
      TopExp::MapShapes(shape, senum, smap);
      n = smap.Extent();
      if (n == 0) return EGADS_SUCCESS;
      objs = (egObject **) EG_alloc(n*sizeof(egObject *));
      if (objs == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc oclass = %d, n = %d (EG_getBodyTopos)!\n",
                 oclass, n);
        return EGADS_MALLOC;
      }
      for (i = 0; i < n; i++) {
        objs[i] = NULL;
        TopoDS_Shape shapo = smap(i+1);
        index = map->map.FindIndex(shapo);
        if (index == 0) {
          if (outLevel > 0)
            printf(" EGADS Warning: %d/%d NotFound oclass = %d (EG_getBodyTopos)!\n",
                   i+1, n, oclass);
        } else {
          objs[i] = map->objs[index-1];
        }
      }

    } else {

      // look up (get super-shapes)
      for (n = i = 0; i < map->map.Extent(); i++) {
        TopoDS_Shape shapo = map->map(i+1);
        TopTools_IndexedMapOfShape smap;
        TopExp::MapShapes(shapo, senum, smap);
        if (smap.FindIndex(shape) != 0) n++;
      }
      if (n == 0) return EGADS_SUCCESS;
      objs = (egObject **) EG_alloc(n*sizeof(egObject *));
      if (objs == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Malloc oclass = %d, N = %d (EG_getBodyTopos)!\n",
                 oclass, n);
        return EGADS_MALLOC;
      }
      for (n = i = 0; i < map->map.Extent(); i++) {
        TopoDS_Shape shapo = map->map(i+1);
        TopTools_IndexedMapOfShape smap;
        TopExp::MapShapes(shapo, senum, smap);
        if (smap.FindIndex(shape) != 0) {
          objs[n] = map->objs[i];
          n++;
        }
      }
    }
  }

  *ntopo = n;
  *topos = objs;

  return EGADS_SUCCESS;
}


int
EG_indexBodyTopo(const egObject *body, const egObject *src)
{
  int outLevel, index;

  if (body == NULL)               return EGADS_NULLOBJ;
  if (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body->oclass != BODY)       return EGADS_NOTBODY;
  if (body->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(body);

  if (src->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: src not an EGO (EG_indexBodyTopo)!\n");
    return EGADS_NOTOBJ;
  }
  if ((src->oclass < NODE) || (src->oclass > SHELL)) {
    if (outLevel > 0)
      printf(" EGADS Error: src not a Topo (EG_indexBodyTopo)!\n");
    return EGADS_NOTTOPO;
  }
  if (EG_context(body) != EG_context(src)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_indexBodyTopo)!\n");
    return EGADS_MIXCNTX;
  }
  if (src->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL src pointer (EG_indexBodyTopo)!\n");
    return EGADS_NODATA;
  }

  egadsBody *pbody = (egadsBody *) body->blind;
  if (src->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) src->blind;
    index = pbody->nodes.map.FindIndex(pnode->node);
  } else if (src->oclass == EDGE) {
    egadsEdge *pedge = (egadsEdge *) src->blind;
    index = pbody->edges.map.FindIndex(pedge->edge);
  } else if (src->oclass == LOOP) {
    egadsLoop *ploop = (egadsLoop *) src->blind;
    index = pbody->loops.map.FindIndex(ploop->loop);
  } else if (src->oclass == FACE) {
    egadsFace *pface = (egadsFace *) src->blind;
    index = pbody->faces.map.FindIndex(pface->face);
  } else {
    egadsShell *pshell = (egadsShell *) src->blind;
    index = pbody->shells.map.FindIndex(pshell->shell);
  }

  if (index == 0) index = EGADS_NOTFOUND;
  return index;
}


int
EG_sameBodyTopo(const egObject *bod1, const egObject *bod2)
{
  int      i, j, outLevel, ind1, ind2, nc, err = 0;
  egObject *obj1, *obj2;

  if (bod1 == NULL)               return EGADS_NULLOBJ;
  if (bod1->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (bod1->oclass != BODY)       return EGADS_NOTBODY;
  if (bod1->blind == NULL)        return EGADS_NODATA;
  if (bod2 == NULL)               return EGADS_NULLOBJ;
  if (bod2->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (bod2->oclass != BODY)       return EGADS_NOTBODY;
  if (bod2->blind == NULL)        return EGADS_NODATA;
  if (bod1->mtype != bod2->mtype) return EGADS_TOPOERR;
  egadsBody *pbod1 = (egadsBody *) bod1->blind;
  egadsBody *pbod2 = (egadsBody *) bod2->blind;
  outLevel = EG_outLevel(bod1);
//  outLevel = 2;

  /* sizes must match */
  if ((bod1->mtype == SHEETBODY) || (bod1->mtype == SOLIDBODY))
    if (pbod1->shells.map.Extent() != pbod2->shells.map.Extent()) {
      if (outLevel > 1)
        printf(" EGADS Warning: #Shells = %d %d (EG_sameBodyTopo)!\n",
               pbod1->shells.map.Extent(), pbod2->shells.map.Extent());
      return EGADS_TOPOCNT;
    }

  if (bod1->mtype != WIREBODY)
    if (pbod1->faces.map.Extent() != pbod2->faces.map.Extent()) {
      if (outLevel > 1)
        printf(" EGADS Warning: #Faces = %d %d (EG_sameBodyTopo)!\n",
               pbod1->faces.map.Extent(), pbod2->faces.map.Extent());
      return EGADS_TOPOCNT;
    }

  if (pbod1->loops.map.Extent() != pbod2->loops.map.Extent()) {
    if (outLevel > 1)
      printf(" EGADS Warning: #Loops = %d %d (EG_sameBodyTopo)!\n",
             pbod1->loops.map.Extent(), pbod2->loops.map.Extent());
    return EGADS_TOPOCNT;
  }
  if (pbod1->edges.map.Extent() != pbod2->edges.map.Extent()) {
    if (outLevel > 1)
      printf(" EGADS Warning: #Edges = %d %d (EG_sameBodyTopo)!\n",
             pbod1->edges.map.Extent(), pbod2->edges.map.Extent());
    return EGADS_TOPOCNT;
  }
  if (pbod1->nodes.map.Extent() != pbod2->nodes.map.Extent()) {
    if (outLevel > 1)
      printf(" EGADS Warning: #Nodes = %d %d (EG_sameBodyTopo)!\n",
             pbod1->nodes.map.Extent(), pbod2->nodes.map.Extent());
    return EGADS_TOPOCNT;
  }

  /* look at children */

  for (i = 0; i < pbod1->edges.map.Extent(); i++) {
    obj1 = pbod1->edges.objs[i];
    obj2 = pbod2->edges.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Edge %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      err++;
      continue;
    }
    egadsEdge *pedg1 = (egadsEdge *) obj1->blind;
    egadsEdge *pedg2 = (egadsEdge *) obj2->blind;
                                nc = 1;
    if (obj1->mtype == TWONODE) nc = 2;
    for (j = 0; j < nc; j++) {
      egObject  *nod1  = pedg1->nodes[j];
      egObject  *nod2  = pedg2->nodes[j];
      egadsNode *pnod1 = (egadsNode *) nod1->blind;
      egadsNode *pnod2 = (egadsNode *) nod2->blind;
      ind1 = pbod1->nodes.map.FindIndex(pnod1->node);
      ind2 = pbod2->nodes.map.FindIndex(pnod2->node);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Edge %d - nIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;

  for (i = 0; i < pbod1->loops.map.Extent(); i++) {
    obj1 = pbod1->loops.objs[i];
    obj2 = pbod2->loops.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Loop %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      return EGADS_TOPOERR;
    }
    egadsLoop *ploo1 = (egadsLoop *) obj1->blind;
    egadsLoop *ploo2 = (egadsLoop *) obj2->blind;
    nc = ploo1->nedges;
    if (nc != ploo2->nedges) {
      if (outLevel > 1)
        printf(" EGADS Warning: Loop %d - nEdges = %d %d (EG_sameBodyTopo)!\n",
               i+1, nc, ploo2->nedges);
      return EGADS_TOPOERR;
    }
    for (j = 0; j < nc; j++) {
      egObject  *edg1  = ploo1->edges[j];
      egObject  *edg2  = ploo2->edges[j];
      egadsEdge *pedg1 = (egadsEdge *) edg1->blind;
      egadsEdge *pedg2 = (egadsEdge *) edg2->blind;
      ind1 = pbod1->edges.map.FindIndex(pedg1->edge);
      ind2 = pbod2->edges.map.FindIndex(pedg2->edge);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Loop %d - eIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;
  if (bod1->mtype == WIREBODY) return EGADS_SUCCESS;

  for (i = 0; i < pbod1->faces.map.Extent(); i++) {
    obj1 = pbod1->faces.objs[i];
    obj2 = pbod2->faces.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Face %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      return EGADS_TOPOERR;
    }
    egadsFace *pfac1 = (egadsFace *) obj1->blind;
    egadsFace *pfac2 = (egadsFace *) obj2->blind;
    nc = pfac1->nloops;
    if (nc != pfac2->nloops) {
      if (outLevel > 1)
        printf(" EGADS Warning: Face %d - nLoops = %d %d (EG_sameBodyTopo)!\n",
               i+1, nc, pfac2->nloops);
      return EGADS_TOPOERR;
    }
    for (j = 0; j < nc; j++) {
      egObject  *loo1  = pfac1->loops[j];
      egObject  *loo2  = pfac2->loops[j];
      egadsLoop *ploo1 = (egadsLoop *) loo1->blind;
      egadsLoop *ploo2 = (egadsLoop *) loo2->blind;
      ind1 = pbod1->loops.map.FindIndex(ploo1->loop);
      ind2 = pbod2->loops.map.FindIndex(ploo2->loop);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Face %d - lIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;
  if (bod1->mtype == FACEBODY) return EGADS_SUCCESS;

  for (i = 0; i < pbod1->shells.map.Extent(); i++) {
    obj1 = pbod1->shells.objs[i];
    obj2 = pbod2->shells.objs[i];
    if (obj1->mtype != obj2->mtype) {
      if (outLevel > 1)
        printf(" EGADS Warning: Shell %d - Types = %d %d (EG_sameBodyTopo)!\n",
               i+1, obj1->mtype, obj2->mtype);
      return EGADS_TOPOERR;
    }
    egadsShell *pshl1 = (egadsShell *) obj1->blind;
    egadsShell *pshl2 = (egadsShell *) obj2->blind;
    nc = pshl1->nfaces;
    if (nc != pshl2->nfaces) {
      if (outLevel > 1)
        printf(" EGADS Warning: Shell %d - nFaces = %d %d (EG_sameBodyTopo)!\n",
               i+1, nc, pshl2->nfaces);
      return EGADS_TOPOERR;
    }
    for (j = 0; j < nc; j++) {
      egObject  *fac1  = pshl1->faces[j];
      egObject  *fac2  = pshl2->faces[j];
      egadsFace *pfac1 = (egadsFace *) fac1->blind;
      egadsFace *pfac2 = (egadsFace *) fac2->blind;
      ind1 = pbod1->faces.map.FindIndex(pfac1->face);
      ind2 = pbod2->faces.map.FindIndex(pfac2->face);
      if (ind1 != ind2) {
        if (outLevel > 1)
          printf(" EGADS Warning: Shell %d - fIndices = %d %d (EG_sameBodyTopo)!\n",
                 i+1, ind1, ind2);
        err++;
      }
    }
  }
  if (err != 0) return EGADS_TOPOERR;

  return EGADS_SUCCESS;
}


int
EG_makeSolidBody(egObject *context, int stypx, const double *data,
                 egObject **body)
{
  int           outLevel, stype, stat, nerr;
  egObject      *obj;
  TopoDS_Shape  solid;
  Standard_Real height;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  outLevel = EG_outLevel(context);
  stype    = abs(stypx);

  if ((stype < BOX) || (stype > TORUS)) {
    if (outLevel > 0)
      printf(" EGADS Error: stype = %d (EG_makeSolidBody)!\n", stype);
    return EGADS_RANGERR;
  }

  switch (stype) {

  /* box=1 */
  case BOX:
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeBox(gp_Pnt(data[0], data[1], data[2]),
                                data[3], data[4], data[5]);
    break;

  /* sphere=2 */
  case SPHERE:
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeSphere(gp_Pnt(data[0], data[1], data[2]),
                                   data[3]);
    break;

  /* cone=3 */
  case CONE:
    height = sqrt( (data[3]-data[0])*(data[3]-data[0]) +
                   (data[4]-data[1])*(data[4]-data[1]) +
                   (data[5]-data[2])*(data[5]-data[2]) );
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(data[0], data[1], data[2]),
                                 gp_Dir(data[3]-data[0], data[4]-data[1],
                                        data[5]-data[2])),
                                 0.0, data[6], height);
    break;

  /* cylinder=4 */
  case CYLINDER:
    height = sqrt( (data[3]-data[0])*(data[3]-data[0]) +
                   (data[4]-data[1])*(data[4]-data[1]) +
                   (data[5]-data[2])*(data[5]-data[2]) );
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(data[0], data[1], data[2]),
                                     gp_Dir(data[3]-data[0], data[4]-data[1],
                                            data[5]-data[2])),
                                     data[6], height);
    break;

  /* torus=5 */
  case TORUS:
    solid = (TopoDS_Solid)
            BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(data[0], data[1], data[2]),
                                  gp_Dir(data[3], data[4], data[5])),
                                  data[6], data[7]);
    break;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_makeSolidBody)!\n");
    return stat;
  }
  egadsBody *pbody   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = SOLIDBODY;
  pbody->nodes.objs  = NULL;
  pbody->edges.objs  = NULL;
  pbody->loops.objs  = NULL;
  pbody->faces.objs  = NULL;
  pbody->shells.objs = NULL;
  pbody->senses      = NULL;
  pbody->shape       = solid;
  obj->blind         = pbody;
  if (stypx > 0) EG_splitPeriodics(pbody);
  stat = EG_traverseBody(context, 0, obj, obj, pbody, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete pbody;
    return stat;
  }

  EG_referenceObject(obj, context);
  *body = obj;
  return EGADS_SUCCESS;
}


int
EG_getBoundingBox(const egObject *topo, double *bbox)
{
  int      i, n;
  egObject *obj;
  Bnd_Box  Box;

  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->oclass < NODE)        return EGADS_NOTTOPO;
  if (topo->blind == NULL)        return EGADS_NODATA;

  try {

    if (topo->oclass == NODE) {

      egadsNode *pnode = (egadsNode *) topo->blind;
      BRepBndLib::Add(pnode->node, Box);

    } else if (topo->oclass == EDGE) {

      egadsEdge *pedge = (egadsEdge *) topo->blind;
      BRepBndLib::Add(pedge->edge, Box);

    } else if (topo->oclass == LOOP) {

      egadsLoop *ploop = (egadsLoop *) topo->blind;
      BRepBndLib::Add(ploop->loop, Box);

    } else if (topo->oclass == FACE) {

      egadsFace *pface = (egadsFace *) topo->blind;
      BRepBndLib::Add(pface->face, Box);

    } else if (topo->oclass == SHELL) {

      egadsShell *pshell = (egadsShell *) topo->blind;
      BRepBndLib::Add(pshell->shell, Box);

    } else if (topo->oclass == BODY) {

      egadsBody *pbody = (egadsBody *) topo->blind;
      BRepBndLib::Add(pbody->shape, Box);

    } else {

      egadsModel *pmodel = (egadsModel *) topo->blind;
      if (pmodel != NULL) {
        n = pmodel->nbody;
        for (i = 0; i < n; i++) {
          obj = pmodel->bodies[i];
          if (obj == NULL) continue;
          egadsBody *pbody = (egadsBody *) obj->blind;
          if (pbody == NULL) continue;
          BRepBndLib::Add(pbody->shape, Box);
        }
      }

    }

    Box.Get(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);

  }
  catch (Standard_Failure)
  {
    printf(" EGADS Warning: BoundingBox failure (EG_getBoundingBox)!\n");
    Handle_Standard_Failure e = Standard_Failure::Caught();
    printf("                %s\n", e->GetMessageString());
    return EGADS_TOPOERR;
  }
  catch (...) {
    printf(" EGADS Warning: BoundingBox failure (EG_getBoundingBox)!\n");
    return EGADS_TOPOERR;
  }

  return EGADS_SUCCESS;
}


int
EG_getMassProperties(const egObject *topo, double *data)
{
  gp_Pnt       CofG;
  gp_Mat       Inert;
  BRepGProp    BProps;
  GProp_GProps SProps, VProps;

  if  (topo == NULL)               return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((topo->oclass < EDGE) ||
      (topo->oclass > BODY))       return EGADS_NOTTOPO;
  if  (topo->blind == NULL)        return EGADS_NODATA;

  data[0] = 0.0;
  if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    BProps.LinearProperties(pedge->edge, SProps);
    CofG  = SProps.CentreOfMass();
    Inert = SProps.MatrixOfInertia();

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    BProps.LinearProperties(ploop->loop, SProps);
    CofG  = SProps.CentreOfMass();
    Inert = SProps.MatrixOfInertia();

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    BProps.SurfaceProperties(pface->face, SProps);
    CofG  = SProps.CentreOfMass();
    Inert = SProps.MatrixOfInertia();

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    BProps.SurfaceProperties(pshell->shell, SProps);
    CofG  = SProps.CentreOfMass();
    Inert = SProps.MatrixOfInertia();

  } else {

    egadsBody *pbody = (egadsBody *) topo->blind;
    if (topo->mtype == SOLIDBODY) {
      BProps.SurfaceProperties(pbody->shape, SProps);
      BProps.VolumeProperties( pbody->shape, VProps);
      CofG    = VProps.CentreOfMass();
      Inert   = VProps.MatrixOfInertia();
      data[0] = VProps.Mass();
    } else if (topo->mtype == WIREBODY) {
      BProps.LinearProperties(pbody->shape, SProps);
      CofG    = SProps.CentreOfMass();
      Inert   = SProps.MatrixOfInertia();
    } else {
      BProps.SurfaceProperties(pbody->shape, SProps);
      CofG    = SProps.CentreOfMass();
      Inert   = SProps.MatrixOfInertia();
    }

  }

  data[ 1] = SProps.Mass();
  data[ 2] = CofG.X();
  data[ 3] = CofG.Y();
  data[ 4] = CofG.Z();
  data[ 5] = Inert.Value(1,1);
  data[ 6] = Inert.Value(1,2);
  data[ 7] = Inert.Value(1,3);
  data[ 8] = Inert.Value(2,1);
  data[ 9] = Inert.Value(2,2);
  data[10] = Inert.Value(2,3);
  data[11] = Inert.Value(3,1);
  data[12] = Inert.Value(3,2);
  data[13] = Inert.Value(3,3);

  return EGADS_SUCCESS;
}


int
EG_isEquivalent(const egObject *topo1, const egObject *topo2)
{
  int          i, stat;
  TopoDS_Shape shape1, shape2;

  if (topo1 == topo2)                 return EGADS_SUCCESS;
  if (topo1 == NULL)                  return EGADS_NULLOBJ;
  if (topo1->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (topo1->oclass < NODE)           return EGADS_NOTTOPO;
  if (topo1->blind == NULL)           return EGADS_NODATA;
  if (topo2 == NULL)                  return EGADS_NULLOBJ;
  if (topo2->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (topo2->oclass != topo1->oclass) return EGADS_NOTTOPO;
  if (topo2->blind == NULL)           return EGADS_NODATA;

  if (topo1->oclass == NODE) {

    egadsNode *pnode1 = (egadsNode *) topo1->blind;
    egadsNode *pnode2 = (egadsNode *) topo2->blind;
    shape1            = pnode1->node;
    shape2            = pnode2->node;

  } else if (topo1->oclass == EDGE) {

    egadsEdge *pedge1 = (egadsEdge *) topo1->blind;
    egadsEdge *pedge2 = (egadsEdge *) topo2->blind;
    shape1            = pedge1->edge;
    shape2            = pedge2->edge;

  } else if (topo1->oclass == LOOP) {

    egadsLoop *ploop1 = (egadsLoop *) topo1->blind;
    egadsLoop *ploop2 = (egadsLoop *) topo2->blind;
    shape1            = ploop1->loop;
    shape2            = ploop2->loop;

  } else if (topo1->oclass == FACE) {

    egadsFace *pface1 = (egadsFace *) topo1->blind;
    egadsFace *pface2 = (egadsFace *) topo2->blind;
    shape1            = pface1->face;
    shape2            = pface2->face;

  } else if (topo1->oclass == SHELL) {

    egadsShell *pshell1 = (egadsShell *) topo1->blind;
    egadsShell *pshell2 = (egadsShell *) topo2->blind;
    shape1              = pshell1->shell;
    shape2              = pshell2->shell;

  } else if (topo1->oclass == BODY) {

    egadsBody *pbody1 = (egadsBody *) topo1->blind;
    egadsBody *pbody2 = (egadsBody *) topo2->blind;
    shape1            = pbody1->shape;
    shape2            = pbody2->shape;

  } else {

    egadsModel *pmodel1 = (egadsModel *) topo1->blind;
    egadsModel *pmodel2 = (egadsModel *) topo2->blind;
    shape1              = pmodel1->shape;
    shape2              = pmodel2->shape;
  }
  if (shape1.IsSame(shape2)) return EGADS_SUCCESS;

  /* try to match via geometry */

  if (topo1->oclass == NODE) {

    return EG_isSame(topo1, topo2);

  } else if (topo1->oclass == EDGE) {

    stat = EG_isSame(topo1, topo2);
    if (stat != EGADS_SUCCESS) return stat;
    egadsEdge *pedge1 = (egadsEdge *) topo1->blind;
    egadsEdge *pedge2 = (egadsEdge *) topo2->blind;
    stat = EG_isSame(pedge1->nodes[0], pedge2->nodes[0]);
    if (stat != EGADS_SUCCESS) return stat;
    return EG_isSame(pedge1->nodes[1], pedge2->nodes[1]);

  } else if (topo1->oclass == LOOP) {

    egadsLoop *ploop1 = (egadsLoop *) topo1->blind;
    egadsLoop *ploop2 = (egadsLoop *) topo2->blind;

    if (ploop1->nedges != ploop2->nedges) return EGADS_OUTSIDE;
    for (i = 0; i < ploop1->nedges; i++) {
      if (ploop1->senses[i] != ploop2->senses[i]) return EGADS_OUTSIDE;
      stat = EG_isEquivalent(ploop1->edges[i], ploop2->edges[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;

  } else if (topo1->oclass == FACE) {

    egadsFace *pface1 = (egadsFace *) topo1->blind;
    egadsFace *pface2 = (egadsFace *) topo2->blind;

    stat = EG_isSame(topo1, topo2);
    if (stat != EGADS_SUCCESS) return stat;
    if (pface1->nloops != pface2->nloops) return EGADS_OUTSIDE;
    for (i = 0; i < pface1->nloops; i++) {
      stat = EG_isEquivalent(pface1->loops[i], pface2->loops[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;

  } else if (topo1->oclass == SHELL) {

    egadsShell *pshell1 = (egadsShell *) topo1->blind;
    egadsShell *pshell2 = (egadsShell *) topo2->blind;

    if (pshell1->nfaces != pshell2->nfaces) return EGADS_OUTSIDE;
    for (i = 0; i < pshell1->nfaces; i++) {
      stat = EG_isEquivalent(pshell1->faces[i], pshell2->faces[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;

  } else if (topo1->oclass == BODY) {

    egadsBody *pbody1 = (egadsBody *) topo1->blind;
    egadsBody *pbody2 = (egadsBody *) topo2->blind;

    if (pbody1->shells.map.Extent() != pbody2->shells.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->faces.map.Extent()  != pbody2->faces.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->loops.map.Extent()  != pbody2->loops.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->edges.map.Extent()  != pbody2->edges.map.Extent())
      return EGADS_OUTSIDE;
    if (pbody1->nodes.map.Extent()  != pbody2->nodes.map.Extent())
      return EGADS_OUTSIDE;

    if (pbody1->shells.map.Extent() != 0) {
      for (i = 0; i < pbody1->shells.map.Extent(); i++) {
        stat = EG_isEquivalent(pbody1->shells.objs[i], pbody2->shells.objs[i]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    } else if (pbody1->faces.map.Extent() != 0) {
      for (i = 0; i < pbody1->faces.map.Extent(); i++) {
        stat = EG_isEquivalent(pbody1->faces.objs[i],  pbody2->faces.objs[i]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    } else {
      for (i = 0; i < pbody1->edges.map.Extent(); i++) {
        stat = EG_isEquivalent(pbody1->edges.objs[i],  pbody2->edges.objs[i]);
        if (stat != EGADS_SUCCESS) return stat;
      }
    }
    return EGADS_SUCCESS;

  }

  return EGADS_OUTSIDE;
}


int
EG_isPlanar(const egObject *topo)
{
  TopoDS_Shape shape;

  if (topo == NULL)                  return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (topo->oclass <= NODE)          return EGADS_NOTTOPO;
  if (topo->blind == NULL)           return EGADS_NODATA;

  if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    shape            = pedge->edge;

  } else if (topo->oclass == LOOP) {

    egadsLoop *ploop = (egadsLoop *) topo->blind;
    shape            = ploop->loop;

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    shape            = pface->face;

  } else if (topo->oclass == SHELL) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    shape              = pshell->shell;

  } else if (topo->oclass == BODY) {

    egadsBody *pbody = (egadsBody *) topo->blind;
    shape            = pbody->shape;

  } else {

    egadsModel *pmodel = (egadsModel *) topo->blind;
    shape              = pmodel->shape;

  }

  // note: single linear edges return not planar!
  BRepBuilderAPI_FindPlane planar(shape);
  if (planar.Found()) return EGADS_SUCCESS;
  return EGADS_OUTSIDE;
}


int
EG_getEdgeUV(const egObject *face, const egObject *topo, int sense, double t,
             double *uv)
{
  int             outLevel, found;
  gp_Pnt2d        P2d;
  TopoDS_Edge     edge;
  TopExp_Explorer ExpW;

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (topo == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_getEdgeUV)!\n");
    return EGADS_NULLOBJ;
  }
  if (topo->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_getEdgeUV)!\n");
    return EGADS_NOTOBJ;
  }
  if ((topo->oclass != EDGE) && (topo->oclass != NODE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_getEdgeUV)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -1) || (sense > 1)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_getEdgeUV)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(topo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_getEdgeUV)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_getEdgeUV)!\n");
    return EGADS_NODATA;
  }

  // undocumented -- topo is a node!
  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Node pointer (EG_getEdgeUV)!\n");
      return EGADS_NODATA;
    }
    TopoDS_Vertex vert = TopoDS::Vertex(pnode->node);

    double tt = 0.0;
    found     = 0;
    for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Wire wire = TopoDS::Wire(ExpW.Current());
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        edge = TopoDS::Edge(ExpWE.Current());
        if (BRep_Tool::Degenerated(edge)) continue;
        TopoDS_Vertex V1, V2;
        if (edge.Orientation() == TopAbs_REVERSED) {
          TopExp::Vertices(edge, V2, V1, Standard_True);
        } else {
          TopExp::Vertices(edge, V1, V2, Standard_True);
        }
        if ((vert.IsSame(V1)) || (vert.IsSame(V2))) {
          double t1, t2;
          Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
          if (vert.IsSame(V1)) {
            tt = t1;
          } else {
            tt = t2;
          }
          found = 1;
          break;
        }
      }
      if (found != 0) break;
    }
    if (found == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Node not in Face (EG_getEdgeUV)!\n");
      return EGADS_NOTFOUND;
    }

    BRepAdaptor_Curve2d Curve2d(edge, pface->face);
    Curve2d.D0(tt, P2d);
    uv[0] = P2d.X();
    uv[1] = P2d.Y();

    return EGADS_SUCCESS;
  }

  // topo is an edge
  egadsEdge *pedge = (egadsEdge *) topo->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_getEdgeUV)!\n");
    return EGADS_NODATA;
  }
  edge = pedge->edge;

  // is edge in the face?

  found = 0;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
      if (!wedge.IsSame(edge)) continue;
      if (sense == 0) {
        found++;
      } else if (sense == -1) {
        if (shape.Orientation() == TopAbs_REVERSED) found++;
      } else {
        if (shape.Orientation() != TopAbs_REVERSED) found++;
      }
      if (found == 0) continue;
      edge = wedge;
      break;
    }
    if (found != 0) break;
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_getEdgeUV)!\n");
    return EGADS_NOTFOUND;
  }
  if (!BRep_Tool::SameRange(edge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge & PCurve not SameRange (EG_getEdgeUV)!\n");
    return EGADS_GEOMERR;
  }
/*
  if (!BRep_Tool::SameParameter(edge)) {
    printf(" EGADS Info: Edge has SameParameter NOT set!\n");
    printf(" Edge tol = %le\n", BRep_Tool::Tolerance(edge));
  }
 */

  // get and evaluate the pcurve

  BRepAdaptor_Curve2d Curve2d(edge, pface->face);
  Curve2d.D0(t, P2d);
  uv[0] = P2d.X();
  uv[1] = P2d.Y();

  return EGADS_SUCCESS;
}


int
EG_getEdgeUVeval(const egObject *face, const egObject *topo, int sense,
                 double t, double *result)
{
  int             outLevel, found;
  gp_Pnt2d        P2d;
  gp_Vec2d        V12d, V22d;
  TopoDS_Edge     edge;
  TopExp_Explorer ExpW;

  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(face);

  if (topo == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Reference (EG_getEdgeUV)!\n");
    return EGADS_NULLOBJ;
  }
  if (topo->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: topo not an EGO (EG_getEdgeUV)!\n");
    return EGADS_NOTOBJ;
  }
  if ((topo->oclass != EDGE) && (topo->oclass != NODE)) {
    if (outLevel > 0)
      printf(" EGADS Error: Not an Edge (EG_getEdgeUV)!\n");
    return EGADS_NOTTOPO;
  }
  if ((sense < -1) || (sense > 1)) {
    if (outLevel > 0)
      printf(" EGADS Error: Sense = %d (EG_getEdgeUV)!\n", sense);
    return EGADS_RANGERR;
  }
  if (EG_context(face) != EG_context(topo)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_getEdgeUV)!\n");
    return EGADS_MIXCNTX;
  }
  egadsFace *pface = (egadsFace *) face->blind;
  if (pface == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face pointer (EG_getEdgeUV)!\n");
    return EGADS_NODATA;
  }

  // undocumented -- topo is a node!
  if (topo->oclass == NODE) {
    egadsNode *pnode = (egadsNode *) topo->blind;
    if (pnode == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Node pointer (EG_getEdgeUV)!\n");
      return EGADS_NODATA;
    }
    TopoDS_Vertex vert = TopoDS::Vertex(pnode->node);

    double tt = 0.0;
    found     = 0;
    for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Wire wire = TopoDS::Wire(ExpW.Current());
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
        edge = TopoDS::Edge(ExpWE.Current());
        if (BRep_Tool::Degenerated(edge)) continue;
        TopoDS_Vertex V1, V2;
        if (edge.Orientation() == TopAbs_REVERSED) {
          TopExp::Vertices(edge, V2, V1, Standard_True);
        } else {
          TopExp::Vertices(edge, V1, V2, Standard_True);
        }
        if ((vert.IsSame(V1)) || (vert.IsSame(V2))) {
          double t1, t2;
          Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
          if (vert.IsSame(V1)) {
            tt = t1;
          } else {
            tt = t2;
          }
          found = 1;
          break;
        }
      }
      if (found != 0) break;
    }
    if (found == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Node not in Face (EG_getEdgeUV)!\n");
      return EGADS_NOTFOUND;
    }

    BRepAdaptor_Curve2d Curve2d(edge, pface->face);
    Curve2d.D2(tt, P2d, V12d, V22d);
    result[0] = P2d.X();
    result[1] = P2d.Y();
    result[2] = V12d.X();
    result[3] = V12d.Y();
    result[4] = V22d.X();
    result[5] = V22d.Y();

    return EGADS_SUCCESS;
  }

  // topo is an edge
  egadsEdge *pedge = (egadsEdge *) topo->blind;
  if (pedge == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Edge pointer (EG_getEdgeUV)!\n");
    return EGADS_NODATA;
  }
  edge = pedge->edge;

  // is edge in the face?

  found = 0;
  for (ExpW.Init(pface->face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
    TopoDS_Shape shapw = ExpW.Current();
    TopoDS_Wire  wire  = TopoDS::Wire(shapw);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shape = ExpWE.Current();
      TopoDS_Edge  wedge = TopoDS::Edge(shape);
      if (!wedge.IsSame(edge)) continue;
      if (sense == 0) {
        found++;
      } else if (sense == -1) {
        if (shape.Orientation() == TopAbs_REVERSED) found++;
      } else {
        if (shape.Orientation() != TopAbs_REVERSED) found++;
      }
      if (found == 0) continue;
      edge = wedge;
      break;
    }
    if (found != 0) break;
  }
  if (found == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge/Sense not in Face (EG_getEdgeUV)!\n");
    return EGADS_NOTFOUND;
  }
  if (!BRep_Tool::SameRange(edge)) {
    if (outLevel > 0)
      printf(" EGADS Error: Edge & PCurve not SameRange (EG_getEdgeUV)!\n");
    return EGADS_GEOMERR;
  }
/*
   if (!BRep_Tool::SameParameter(edge)) {
     printf(" EGADS Info: Edge has SameParameter NOT set!\n");
     printf(" Edge tol = %le\n", BRep_Tool::Tolerance(edge));
   }
*/

  // get and evaluate the pcurve

  BRepAdaptor_Curve2d Curve2d(edge, pface->face);
  Curve2d.D2(t, P2d, V12d, V22d);
  result[0] = P2d.X();
  result[1] = P2d.Y();
  result[2] = V12d.X();
  result[3] = V12d.Y();
  result[4] = V22d.X();
  result[5] = V22d.Y();

  return EGADS_SUCCESS;
}


int
EG_getBody(const egObject *obj, egObject **body)
{
  int i, index;

  *body = NULL;
  if (obj == NULL)                  return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC)    return EGADS_NOTOBJ;
  if (obj->blind == NULL)           return EGADS_NODATA;
  if ((obj->oclass < NODE) ||
      (obj->oclass > SHELL))        return EGADS_NOTTOPO;
  egObject *topObj = obj->topObj;
  if (topObj == NULL)               return EGADS_NULLOBJ;
  if (topObj->magicnumber != MAGIC) return EGADS_NOTOBJ;

  if (topObj->oclass == BODY) {
    *body = topObj;
  } else if (topObj->oclass == MODEL) {
    egadsModel *pmodel = (egadsModel *) topObj->blind;
    if (pmodel != NULL) {
      for (i = 0; i < pmodel->nbody; i++) {
        egObject *bod = pmodel->bodies[i];
        if (bod == NULL) continue;
        egadsBody *pbody = (egadsBody *) bod->blind;
        if (obj->oclass == NODE) {
          egadsNode *pnode = (egadsNode *) obj->blind;
          index = pbody->nodes.map.FindIndex(pnode->node);
        } else if (obj->oclass == EDGE) {
          egadsEdge *pedge = (egadsEdge *) obj->blind;
          index = pbody->edges.map.FindIndex(pedge->edge);
        } else if (obj->oclass == LOOP) {
          egadsLoop *ploop = (egadsLoop *) obj->blind;
          index = pbody->loops.map.FindIndex(ploop->loop);
        } else if (obj->oclass == FACE) {
          egadsFace *pface = (egadsFace *) obj->blind;
          index = pbody->faces.map.FindIndex(pface->face);
        } else {
          egadsShell *pshell = (egadsShell *) obj->blind;
          index = pbody->shells.map.FindIndex(pshell->shell);
        }
        if (index != 0) {
          *body = bod;
          break;
        }
      }
    }
  }

  return EGADS_SUCCESS;
}

int
EG_inTopology(const egObject *topo, const double *xyz)
{
  int           outLevel;
  Standard_Real tol, t, u, v, range[2];

  if (topo == NULL)               return EGADS_NULLOBJ;
  if (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (topo->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(topo);

  gp_Pnt pnt(xyz[0], xyz[1], xyz[2]);

  if (topo->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) topo->blind;
    egObject  *curv  = pedge->curve;
    if (curv == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Object for Edge (EG_inTopology)!\n");
      return EGADS_NULLOBJ;
    }
    egadsCurve *pcurve = (egadsCurve *) curv->blind;
    if (pcurve == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No curve Data for Edge (EG_inTopology)!\n");
      return EGADS_NODATA;
    }
    Handle(Geom_Curve) hCurve = pcurve->handle;
    GeomAPI_ProjectPointOnCurve projPnt(pnt, hCurve);
    if (projPnt.NbPoints() == 0) {
      if (outLevel > 0)
        printf(" EGADS Warning: No projection on Curve (EG_inTopology)!\n");
      return EGADS_NOTFOUND;
    }
    tol = BRep_Tool::Tolerance(pedge->edge);
    if (projPnt.LowerDistance() > tol) return EGADS_OUTSIDE;
    t = projPnt.LowerDistanceParameter();
    BRep_Tool::Range(pedge->edge, range[0],range[1]);
    if ((t < range[0]) || (t > range[1])) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  } else if (topo->oclass == FACE) {

    egadsFace *pface = (egadsFace *) topo->blind;
    egObject  *surf  = pface->surface;
    if (surf == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No Surf Object for Face (EG_inTopology)!\n");
      return EGADS_NULLOBJ;
    }
    egadsSurface *psurf = (egadsSurface *) surf->blind;
    if (psurf == NULL) {
      if (outLevel > 0)
        printf(" EGADS Warning: No Surf Data for Face (EG_inTopology)!\n");
      return EGADS_NODATA;
    }
    GeomAPI_ProjectPointOnSurf projPnt(pnt, psurf->handle);
    if (!projPnt.IsDone()) {
      printf(" EGADS Warning: GeomAPI_ProjectPointOnSurf (EG_inTopology)!\n");
      return EGADS_GEOMERR;
    }
    tol = BRep_Tool::Tolerance(pface->face);
    if (projPnt.LowerDistance() > tol) return EGADS_OUTSIDE;
    projPnt.LowerDistanceParameters(u, v);
    gp_Pnt2d pnt2d(u, v);
    TopOpeBRep_PointClassifier pClass;
    pClass.Load(pface->face);
    if (pClass.Classify(pface->face, pnt2d, tol) == TopAbs_OUT)
      return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  } else if ((topo->oclass == SHELL) && (topo->mtype == CLOSED)) {

    egadsShell *pshell = (egadsShell *) topo->blind;
    TopOpeBRepTool_SolidClassifier sClass;
    sClass.LoadShell(pshell->shell);
    if (sClass.Classify(pshell->shell, pnt,
          Precision::Confusion()) == TopAbs_OUT) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  } else if ((topo->oclass == BODY) && (topo->mtype == SOLIDBODY)) {

    egadsBody   *pbody = (egadsBody *) topo->blind;
    TopoDS_Solid solid = TopoDS::Solid(pbody->shape);
    TopOpeBRepTool_SolidClassifier sClass;
    sClass.LoadSolid(solid);
    if (sClass.Classify(solid, pnt, Precision::Confusion()) == TopAbs_OUT)
      return EGADS_OUTSIDE;
    return EGADS_SUCCESS;

  }

  return EGADS_NOTTOPO;
}


int
EG_inFace(const egObject *face, const double *uv)
{
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;

  egadsFace  *pface = (egadsFace *) face->blind;
  Standard_Real tol = BRep_Tool::Tolerance(pface->face);
  gp_Pnt2d pnt2d(uv[0], uv[1]);
  TopOpeBRep_PointClassifier pClass;
  pClass.Load(pface->face);
  if (pClass.Classify(pface->face, pnt2d, tol) == TopAbs_OUT)
    return EGADS_OUTSIDE;

  return EGADS_SUCCESS;
}


int
EG_sewFaces(int nobj, const egObject **objs, double toler, int opt,
            egObject **result)
{
  int      i, j, k, n, outLevel, stat, nerr, *amap;
  double   tol;
  egObject *context, *omodel;

  *result = NULL;
  if (nobj <= 1)                     return EGADS_EMPTY;
  if (objs == NULL)                  return EGADS_NULLOBJ;
  if (objs[0] == NULL)               return EGADS_NULLOBJ;
  if (objs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  outLevel = EG_outLevel(objs[0]);
  context  = EG_context(objs[0]);

  tol = 0.0;
  for (i = 0; i < nobj; i++) {
    if (objs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Object %d (EG_sewFaces)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (objs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is not an EGO (EG_sewFaces)!\n", i+1);
      return EGADS_NOTOBJ;
    }
    if (objs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d has no data (EG_sewFaces)!\n", i+1);
      return EGADS_NODATA;
    }
    if (EG_context(objs[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d Context Mismatch (EG_sewFaces)!\n",
               i+1);
      return EGADS_MIXCNTX;
    }
    if (objs[i]->oclass == BODY) {
      if (objs[i]->mtype == WIREBODY) {
        if (outLevel > 0)
          printf(" EGADS Error: Object %d is a WireBody (EG_sewFaces)!\n",
                 i+1);
        return EGADS_NOTTOPO;
      }
      egadsBody *pbody = (egadsBody *) objs[i]->blind;
      int nface = pbody->faces.map.Extent();
      for (i = 0; i < nface; i++) {
        TopoDS_Face facet = TopoDS::Face(pbody->faces.map(i+1));
        double toltmp = BRep_Tool::Tolerance(facet);
        if (toltmp > tol) tol = toltmp;
      }
    } else if (objs[i]->oclass == SHELL) {
      egadsShell *pshell = (egadsShell *) objs[i]->blind;
      for (i = 0; i < pshell->nfaces; i++) {
        egadsFace *pface = (egadsFace *) pshell->faces[i]->blind;
        if (pface == NULL) continue;
        double toltmp = BRep_Tool::Tolerance(pface->face);
        if (toltmp > tol) tol = toltmp;
      }
    } else if (objs[i]->oclass == FACE) {
      egadsFace *pface = (egadsFace *) objs[i]->blind;
      double toltmp = BRep_Tool::Tolerance(pface->face);
      if (toltmp > tol) tol = toltmp;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d Does Not have Faces (EG_sewFaces)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
  }
  if (toler > tol) tol = toler;

  Standard_Boolean flag = Standard_False;
  if (opt == 1)    flag = Standard_True;
  BRepBuilderAPI_Sewing sew(tol, Standard_True, Standard_True, Standard_True,
                            flag);
  for (i = 0; i < nobj; i++) {
    TopoDS_Shape shape;
    if (objs[i]->oclass == BODY) {
      egadsBody *pbody = (egadsBody *) objs[i]->blind;
      shape = pbody->shape;
    } else if (objs[i]->oclass == SHELL) {
      egadsShell *pshell = (egadsShell *) objs[i]->blind;
      shape = pshell->shell;
    } else {
      egadsFace *pface = (egadsFace *) objs[i]->blind;
      shape = pface->face;
    }
    sew.Add(shape);
  }
  sew.Perform();
  TopoDS_Shape sewShape = sew.SewedShape();

  // map the faces
  TopTools_IndexedMapOfShape smap;
  TopExp::MapShapes(sewShape, TopAbs_FACE, smap);
  amap = NULL;
  n    = smap.Extent();
  if (n > 0) {
    amap = (int *) EG_alloc(2*n*sizeof(int));
    if (amap != NULL) {
      for (j = 0; j < 2*n; j++) amap[j] = 0;
      for (i = 0; i < nobj; i++) {
        if (objs[i]->oclass == BODY) {
          egadsBody *pbodf = (egadsBody *) objs[i]->blind;
          int nface = pbodf->faces.map.Extent();
          for (j = 0; j < nface; j++) {
            TopoDS_Face faci =
                     TopoDS::Face(sew.ModifiedSubShape(pbodf->faces.map(j+1)));
            k = smap.FindIndex(faci) - 1;
            if (k < 0) continue;
            amap[2*k]   = i+1;
            amap[2*k+1] = j+1;
          }
        } else if (objs[i]->oclass == SHELL) {
          egadsShell *pshelf = (egadsShell *) objs[i]->blind;
          for (j = 0; j < pshelf->nfaces; j++) {
            egadsFace *pfacf = (egadsFace *) pshelf->faces[j]->blind;
            if (pfacf == NULL) continue;
            TopoDS_Face faci = TopoDS::Face(sew.ModifiedSubShape(pfacf->face));
            k = smap.FindIndex(faci) - 1;
            if (k < 0) continue;
            amap[2*k]   = i+1;
            amap[2*k+1] = j+1;
          }
        } else {
          egadsFace *pface = (egadsFace *) objs[i]->blind;
          TopoDS_Face faci = TopoDS::Face(sew.Modified(pface->face));
          k = smap.FindIndex(faci) - 1;
          if (k < 0) continue;
          amap[2*k]   = i+1;
          amap[2*k+1] = 1;
        }
      }
    }
  }

  // check for promoting sheets to solids
  TopExp_Explorer Exp;
  for (Exp.Init(sewShape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) {
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Shell shell = TopoDS::Shell(shape);
    BRep_Builder builder3D;
    TopoDS_Solid solid;
    builder3D.MakeSolid(solid);
    builder3D.Add(solid, shell);
    try {
      BRepGProp    BProps;
      GProp_GProps VProps;

      BRepLib::OrientClosedSolid(solid);
      BProps.VolumeProperties(solid, VProps);
      if (VProps.Mass() < 0.0) solid.Reverse();
    }
    catch (Standard_Failure)
    {
/*    printf(" EGADS Warning: Cannot Orient Solid (EG_sewFaces)!\n");
      Handle_Standard_Failure e = Standard_Failure::Caught();
      printf("                %s\n", e->GetMessageString()); */
      solid.Nullify();
    }
    catch (...) {
//    printf(" EGADS Warning: Cannot Orient Solid (EG_sewFaces)!\n");
      solid.Nullify();
    }
    if (!solid.IsNull()) {
      BRepCheck_Analyzer sCheck(solid);
      if (!sCheck.IsValid()) {
//      printf(" EGADS Warning: Invalid Solid (EG_sewFaces)!\n");
        solid.Nullify();
      }
    }
    if (!solid.IsNull()) {
      if (sewShape.ShapeType() == TopAbs_SHELL) {
        sewShape = solid;
        break;
      } else {
        builder3D.Add(sewShape, solid);
        builder3D.Remove(sewShape, shape);
      }
    }
  }

  // count our bodies
  int nFace  = 0;
  int nSheet = 0;
  int nSolid = 0;
  for (Exp.Init(sewShape, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) nFace++;
  for (Exp.Init(sewShape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) nSheet++;
  for (Exp.Init(sewShape, TopAbs_SOLID); Exp.More(); Exp.Next()) nSolid++;

  if (outLevel > 1)
    printf(" EGADS Info: Sewn Result has %d Solids, %d Sheets and %d Faces\n",
           nSolid, nSheet, nFace);

  int nBody = nFace+nSheet+nSolid;
  if (nBody == 0) {
    if (amap != NULL) EG_free(amap);
    sewShape.Nullify();
    if (outLevel > 0)
      printf(" EGADS Warning: Nothing found in Result (EG_sewFaces)!\n");
    return EGADS_NODATA;
  }

  egadsModel *mshape = new egadsModel;
  mshape->shape      = sewShape;
  mshape->nbody      = nBody;
  mshape->bodies     = new egObject*[nBody];
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (int j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      if (amap != NULL) EG_free(amap);
      return stat;
    }
    egObject  *pobj    = mshape->bodies[i];
    egadsBody *pbody   = new egadsBody;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pobj->blind        = pbody;
  }

  i = 0;
  for (Exp.Init(mshape->shape, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }
  for (Exp.Init(mshape->shape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
  }
  for (Exp.Init(mshape->shape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
    egObject  *obj     = mshape->bodies[i++];
    egadsBody *pbody   = (egadsBody *) obj->blind;
    TopoDS_Shape shape = Exp.Current();
    TopoDS_Solid solid = TopoDS::Solid(shape);
    BRepCheck_Analyzer sCheck(solid);
    if (!sCheck.IsValid()) {
      Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(solid);
      sfs->Perform();
      TopoDS_Shape fixedSolid = sfs->Shape();
      if (fixedSolid.IsNull()) {
        printf(" EGADS Warning: Cannot Fix Solid (EG_sewFaces)!\n");
      } else {
        if (outLevel > 0)
          printf(" EGADS Warning: Fixing Solid (EG_sewFaces)!\n");
        solid = TopoDS::Solid(fixedSolid);
      }
    }
    pbody->shape = solid;
  }

  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    sewShape.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    if (amap != NULL) EG_free(amap);
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  EG_referenceObject(omodel, context);

  for (i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      mshape->nbody = i;
      EG_destroyTopology(omodel);
      delete [] mshape->bodies;
      delete mshape;
      if (amap != NULL) EG_free(amap);
      return stat;
    }
  }

  // set the attributes
  if (amap != NULL) {
    TopTools_IndexedMapOfShape mmap;
    TopExp::MapShapes(sewShape, TopAbs_FACE, mmap);
    for (i = 0; i < nBody; i++) {
      egObject  *pobj  = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      for (j = 0; j < pbody->faces.map.Extent(); j++) {
        TopoDS_Face facb = TopoDS::Face(pbody->faces.map(j+1));
        int index = mmap.FindIndex(facb) - 1;
        if (index < 0) continue;
//      printf(" body %d: face %d -> %d %d\n", i+1, j+1, index, amap[2*index]);
        k     = amap[2*index+1] - 1;
        index = amap[2*index  ] - 1;
        if (index < 0) continue;
        if (objs[index]->oclass == BODY) {
          egadsBody *pbodf = (egadsBody *) objs[index]->blind;
          if ((k < 0) || (k >= pbodf->faces.map.Extent())) continue;
          EG_attributeDup(pbodf->faces.objs[k], pbody->faces.objs[j]);
        } else if (objs[index]->oclass == SHELL) {
          egadsShell *pshelf = (egadsShell *) objs[index]->blind;
          if ((k < 0) || (k >= pshelf->nfaces)) continue;
          EG_attributeDup(pshelf->faces[k], pbody->faces.objs[j]);
        } else {
          EG_attributeDup(objs[index], pbody->faces.objs[j]);
        }
      }
    }
    EG_free(amap);
  }

  *result = omodel;
  return EGADS_SUCCESS;
}


int
EG_replaceFaces(const egObject *body, int nobj, egObject **objs,
                      egObject **result)
{
  int      i, j, outLevel, stat, mtype, nerr;
  egObject *context, *obj, *src;

  *result = NULL;
  if  (body == NULL)               return EGADS_NULLOBJ;
  if  (body->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if  (body->oclass != BODY)       return EGADS_NOTBODY;
  if  (body->blind == NULL)        return EGADS_NODATA;
  if ((body->mtype != SHEETBODY) &&
      (body->mtype != SOLIDBODY))  return EGADS_TOPOERR;
  if  (nobj < 1)                   return EGADS_EMPTY;
  if  (objs == NULL)               return EGADS_NULLOBJ;
  outLevel = EG_outLevel(body);
  context  = EG_context(body);

  // check the input objects
  egadsBody *pbody = (egadsBody *) body->blind;
  if (body->mtype == SOLIDBODY)
    if (pbody->shells.map.Extent() != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: SolidBody with > 1 Shell (EG_replaceFaces)!\n");
      return EGADS_TOPOERR;
    }
  for (i = 0; i < nobj; i++) {
    if (objs[2*i  ] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Src Object %d (EG_replaceFaces)!\n", i+1);
      return EGADS_NULLOBJ;
    }
    if (objs[2*i  ]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Src Object %d is not an EGO (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (objs[2*i  ]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Src Object %d has no data (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NODATA;
    }
    if (objs[2*i  ]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Src Object %d is Not a Face (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
    egadsFace *pface = (egadsFace *) objs[2*i  ]->blind;
    if (pbody->faces.map.FindIndex(pface->face) == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Object %d is NOT in Body (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTBODY;
    }

    if (objs[2*i+1] == NULL) continue;
    if (objs[2*i+1]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Object %d is not an EGO (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTOBJ;
    }
    if (objs[2*i+1]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Object %d has no data (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NODATA;
    }
    if (EG_context(objs[2*i+1]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Obj %d Context Mismatch (EG_replaceFaces)!\n",
               i+1);
      return EGADS_MIXCNTX;
    }
    if (objs[2*i+1]->oclass != FACE) {
      if (outLevel > 0)
        printf(" EGADS Error: Repl Object %d is Not a Face (EG_replaceFaces)!\n",
               i+1);
      return EGADS_NOTTOPO;
    }
  }

  // make the changes
  BRepTools_ReShape reshape;
  for (i = 0; i < nobj; i++) {
    egadsFace *pface = (egadsFace *) objs[2*i  ]->blind;
    if (objs[2*i+1] == NULL) {
      reshape.Remove(pface->face);
    } else {
      egadsFace *pfacn = (egadsFace *) objs[2*i+1]->blind;
      reshape.Replace(pface->face, pfacn->face);
    }
  }

  // get mappings for attributes
  int nNull = 0;
  int *amap = new int[pbody->faces.map.Extent()];
  for (i = 1; i <= pbody->faces.map.Extent(); i++) {
    TopoDS_Face face = TopoDS::Face(pbody->faces.map(i));
    TopoDS_Face facn = TopoDS::Face(reshape.Value(face));
    if (facn.IsNull()) {
      amap[i-1] = 0;
      nNull++;
    } else if (facn.IsSame(face)) {
      amap[i-1] = i;
    } else {
      for (j = 0; j < nobj; j++) {
        if (objs[2*j+1] == NULL) continue;
        egadsFace *pfacn = (egadsFace *) objs[2*j+1]->blind;
        if (facn.IsSame(pfacn->face)) {
          amap[i-1] = -(j+1);
          break;
        }
      }
      if (j == nobj) {
        if (outLevel > 0)
          printf(" EGADS Error: New Face not Found (EG_replaceFaces)!\n");
        delete [] amap;
        return EGADS_NOTFOUND;
      }
    }
  }

  // apply the changes
  TopoDS_Shape newShape;
  if (body->mtype == SOLIDBODY) {
    egadsShell *pshell = (egadsShell *) pbody->shells.objs[0]->blind;
    newShape = reshape.Apply(pshell->shell, TopAbs_FACE);
  } else {
    newShape = reshape.Apply(pbody->shape, TopAbs_FACE);
  }

  // check our result -- sheet body
  if (newShape.ShapeType() != TopAbs_SHELL) {
    if (outLevel > 0)
      printf(" EGADS Error: Result Not a Sheet (EG_replaceFaces)!\n");
    delete [] amap;
    return EGADS_CONSTERR;
  }
  mtype = SHEETBODY;

  // promote to a solid?
  TopoDS_Shell shell = TopoDS::Shell(newShape);
  BRep_Builder builder3D;
  TopoDS_Solid solid;
  builder3D.MakeSolid(solid);
  builder3D.Add(solid, shell);
  try {
    BRepLib::OrientClosedSolid(solid);
  }
  catch (Standard_Failure) {
/*  printf(" EGADS Warning: Cannot Orient Solid (EG_replaceFaces)!\n");
    Handle_Standard_Failure e = Standard_Failure::Caught();
    printf("                %s\n", e->GetMessageString()); */
    solid.Nullify();
  }
  catch (...) {
//  printf(" EGADS Warning: Cannot Orient Solid (EG_replaceFaces)!\n");
    solid.Nullify();
  }
  if (!solid.IsNull()) {
    BRepCheck_Analyzer sCheck(solid);
    if (!sCheck.IsValid()) {
      printf(" EGADS Warning: Invalid Solid (EG_replaceFaces)!\n");
      solid.Nullify();
    }
  }
  if (!solid.IsNull()) {
    newShape = solid;
    mtype    = SOLIDBODY;
  }

  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_replaceFaces)!\n");
    delete [] amap;
    return stat;
  }
  egadsBody *pbodn   = new egadsBody;
  obj->oclass        = BODY;
  obj->mtype         = mtype;
  pbodn->nodes.objs  = NULL;
  pbodn->edges.objs  = NULL;
  pbodn->loops.objs  = NULL;
  pbodn->faces.objs  = NULL;
  pbodn->shells.objs = NULL;
  pbodn->senses      = NULL;
  pbodn->shape       = newShape;
  obj->blind         = pbodn;
  stat = EG_traverseBody(context, 0, obj, obj, pbodn, &nerr);
  if (stat != EGADS_SUCCESS) {
    delete [] amap;
    delete pbodn;
    return stat;
  }

  // assign attributes
  if (pbody->faces.map.Extent() != pbodn->faces.map.Extent()+nNull) {
    if (outLevel > 0)
      printf(" EGADS Warning: Attribute length Mismatch (EG_replaceFaces)!\n");
  } else {
    for (j = i = 0; i < pbody->faces.map.Extent(); i++) {
      if (amap[i] == 0) continue;
      if (amap[i] >  0) {
        src   = pbody->faces.objs[amap[i]-1];
      } else {
        int k = -amap[i]-1;
        src   = objs[2*k+1];
      }
      EG_attributeDup(src, pbodn->faces.objs[j]);
      j++;
    }
  }
  delete [] amap;

  EG_referenceObject(obj, context);
  *result = obj;

  return EGADS_SUCCESS;
}


static int
EG_edgeCmp(double *cg0, double *cg1)
{
  if (fabs(cg0[0]-cg1[0]) > 1.e-6) {
    if (cg0[0] < cg1[0]) return 1;
    return 0;
  }
  if (fabs(cg0[1]-cg1[1]) > 1.e-6) {
    if (cg0[1] < cg1[1]) return 1;
    return 0;
  }
  if (fabs(cg0[2]-cg1[2]) > 1.e-6) {
    if (cg0[2] < cg1[2]) return 1;
    return 0;
  }

  if (cg0[3] < cg1[3]) return 1;
  return 0;
}


static int
EG_getEdgeIDs(egadsBody *pbody, const char *fAttr, edgeID **IDs)
{
  int          i, j, k, m, n, len, index, cnt, hit, stat, aType, aLen;
  char         **facAttr, line[1025];
  const int    *ints;
  const char   *str;
  const double *reals;
  edgeID       *edgeIDs;
  gp_Pnt       CofG;
  BRepGProp    BProps;
  GProp_GProps SProps;

  *IDs    = NULL;
  len     = pbody->edges.map.Extent();
  edgeIDs = (edgeID *) EG_alloc(len*sizeof(edgeID));
  if (edgeIDs == NULL) return EGADS_MALLOC;
  for (i = 0; i < len; i++) {
    edgeIDs[i].nFace    = 0;
    edgeIDs[i].seq      = 0;
    edgeIDs[i].fIndices = NULL;
    edgeIDs[i].ID       = NULL;
  }

  // find # of Faces per Edge
  for (i = 1; i <= pbody->faces.map.Extent(); i++) {
    TopTools_IndexedMapOfShape smap;
    TopExp::MapShapes(pbody->faces.map(i), TopAbs_EDGE, smap);
    n = smap.Extent();
    for (j = 1; j <= n; j++) {
      index = pbody->edges.map.FindIndex(smap(j));
      if (index == 0) {
        printf(" EGADS Internal: Face %d -- cannot find edge %d!\n", i, j);
        continue;
      }
      edgeIDs[index-1].nFace++;
    }
  }

  // allocate our space per Edge
  for (i = 0; i < len; i++) {
    if (edgeIDs[i].nFace == 0) continue;
    edgeIDs[i].fIndices = (int *) EG_alloc((edgeIDs[i].nFace+1)*sizeof(int));
    if (edgeIDs[i].fIndices == NULL) {
      for (j = 0; j < i; j++)
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
      EG_free(edgeIDs);
      return EGADS_MALLOC;
    }
    for (j = 0; j < edgeIDs[i].nFace+1; j++) edgeIDs[i].fIndices[j] = 0;
    edgeIDs[i].nFace = 0;
  }

  // store Face indices
  for (i = 1; i <= pbody->faces.map.Extent(); i++) {
    TopTools_IndexedMapOfShape smap;
    TopExp::MapShapes(pbody->faces.map(i), TopAbs_EDGE, smap);
    n = smap.Extent();
    for (j = 1; j <= n; j++) {
      index = pbody->edges.map.FindIndex(smap(j));
      edgeIDs[index-1].fIndices[edgeIDs[index-1].nFace] = i;
      edgeIDs[index-1].nFace++;
    }
  }

  // look for duplicates
  for (i = 0; i < len-1; i++) {
    cnt = 0;
    if (edgeIDs[i].fIndices[edgeIDs[i].nFace] != 0) continue;
    for (j = i+1; j < len; j++) {
      if (edgeIDs[i].nFace != edgeIDs[j].nFace) continue;
      for (k = 0; k < edgeIDs[i].nFace; k++)
        if (edgeIDs[i].fIndices[k] != edgeIDs[j].fIndices[k]) break;
      if (k == edgeIDs[i].nFace) {
        if (edgeIDs[i].fIndices[k] == 0) {
          BProps.LinearProperties(pbody->edges.map(i+1), SProps);
          CofG = SProps.CentreOfMass();
          edgeIDs[i].seq         = 1;
          edgeIDs[i].fIndices[k] = -i-1;
          edgeIDs[i].CG[0]       = CofG.X();
          edgeIDs[i].CG[1]       = CofG.Y();
          edgeIDs[i].CG[2]       = CofG.Z();
          edgeIDs[i].CG[3]       = SProps.Mass();
/*        printf("  %d: %lf %lf %lf %lf\n", i+1, CofG.X(), CofG.Y(), CofG.Z(),
                 SProps.Mass());  */
        }
        BProps.LinearProperties(pbody->edges.map(j+1), SProps);
        CofG = SProps.CentreOfMass();
        edgeIDs[j].seq         = cnt+2;
        edgeIDs[j].fIndices[k] = -i-1;
        edgeIDs[j].CG[0]       = CofG.X();
        edgeIDs[j].CG[1]       = CofG.Y();
        edgeIDs[j].CG[2]       = CofG.Z();
        edgeIDs[j].CG[3]       = SProps.Mass();
/*      printf("  %d: %lf %lf %lf %lf\n", j+1, CofG.X(), CofG.Y(), CofG.Z(),
               SProps.Mass());  */
        cnt++;
      }
    }
    if (cnt != 0) {
      do {
        hit = 0;
        for (m = 1; m <= cnt; m++) {
          for (j = 0; j < len; j++)
            if ((edgeIDs[j].fIndices[edgeIDs[j].nFace] == -i-1) &&
                (edgeIDs[j].seq == m)) break;
          if (j == len) {
            printf("  Not found ERROR  %d!\n", m);
            continue;
          }
          for (k = 0; k < len; k++)
            if ((edgeIDs[k].fIndices[edgeIDs[k].nFace] == -i-1) &&
                (edgeIDs[k].seq == m+1)) break;
          if (k == len) {
            printf("  Not found ERROR  %d!\n", m+1);
            continue;
          }
          if (EG_edgeCmp(edgeIDs[j].CG, edgeIDs[k].CG) == 1) continue;
          edgeIDs[k].seq = m;
          edgeIDs[j].seq = m+1;
          hit++;
        }
      } while (hit != 0);
    }
  }

  // convert FaceIDs to strings

  facAttr = (char **) EG_alloc(pbody->faces.map.Extent()*sizeof(char *));
  if (facAttr == NULL) {
    for (j = 0; j < len; j++)
      if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
    EG_free(edgeIDs);
    return EGADS_MALLOC;
  }
  for (i = 0; i < pbody->faces.map.Extent(); i++) facAttr[i] = NULL;
  for (i = 0; i < pbody->faces.map.Extent(); i++) {
    stat = EG_attributeRet(pbody->faces.objs[i], fAttr, &aType, &aLen,
                           &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) {
      printf(" EGADS Error: Face Attr = %d for Face %d (EG_mapBody)!\n",
             stat, i+1);
      for (j = 0; j < i; j++) EG_free(facAttr[i]);
      EG_free(facAttr);
      for (j = 0; j < len; j++)
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
      EG_free(edgeIDs);
      return stat;
    }
    if (aType == ATTRSTRING) {
      facAttr[i] = EG_strdup(str);
    } else if (aType == ATTRREAL) {
      sprintf(line, "%20.13le", reals[0]);
      for (j = 1; j < aLen; j++) {
        hit = strlen(line);
        if (hit+21 > 1024) {
          printf(" EGADS Warning: FaceID %d truncated to %d ints (EG_mapBody)!\n",
                 i+1, j);
          break;
        }
        sprintf(&line[hit], ":%20.13le", reals[j]);
      }
      facAttr[i] = EG_strdup(line);
    } else {
      sprintf(line, "%d", ints[0]);
      for (j = 1; j < aLen; j++) {
        hit = strlen(line);
        if (hit+10 > 1024) {
          printf(" EGADS Warning: FaceID %d truncated to %d ints (EG_mapBody)!\n",
                 i+1, j);
          break;
        }
        sprintf(&line[hit], ":%d", ints[j]);
      }
      facAttr[i] = EG_strdup(line);
    }
  }

  // do we have all of the Face strings?
  for (i = 0; i < pbody->faces.map.Extent(); i++)
    if (facAttr[i] == NULL) {
      printf(" EGADS Error: FaceID %d converted to NULL (EG_mapBody)!\n",
             i+1);
      for (i = 0; i < pbody->faces.map.Extent(); i++) EG_free(facAttr[i]);
      EG_free(facAttr);
      for (j = 0; j < len; j++)
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
      EG_free(edgeIDs);
      return EGADS_MALLOC;
    }

  // make the Edge ID strings from the Faces

  for (i = 0; i < len; i++) {
    // order the Faces based on strcmps
    do {
      hit = 0;
      for (m = 0; m < edgeIDs[i].nFace-1; m++) {
        if (strcmp(facAttr[edgeIDs[i].fIndices[m  ]-1],
                   facAttr[edgeIDs[i].fIndices[m+1]-1]) <= 0) continue;
        k                        = edgeIDs[i].fIndices[m  ];
        edgeIDs[i].fIndices[m  ] = edgeIDs[i].fIndices[m+1];
        edgeIDs[i].fIndices[m+1] = k;
        hit++;
      }
    } while (hit != 0);
    for (cnt = j = 0; j < edgeIDs[i].nFace; j++)
      cnt += strlen(facAttr[edgeIDs[i].fIndices[j]-1]) + 1;
    edgeIDs[i].ID = (char *) EG_alloc((cnt+11)*sizeof(char));
    if (edgeIDs[i].ID == NULL) continue;
    // concatinate the ordered Face strings
    for (cnt = j = 0; j < edgeIDs[i].nFace; j++) {
      sprintf(&edgeIDs[i].ID[cnt], "%s|", facAttr[edgeIDs[i].fIndices[j]-1]);
      cnt += strlen(facAttr[edgeIDs[i].fIndices[j]-1]) + 1;
    }
    // append the sequence number
    cnt--;
    sprintf(&edgeIDs[i].ID[cnt], "\\%d", edgeIDs[i].seq);
  }
  for (i = 0; i < pbody->faces.map.Extent(); i++) EG_free(facAttr[i]);
  EG_free(facAttr);

  // do we have them all?
  for (i = 0; i < pbody->edges.map.Extent(); i++) {
/*  printf("  Edge %d: %d Faces ->", i+1, edgeIDs[i].nFace);
    for (j = 0; j < edgeIDs[i].nFace; j++)
      printf("  %d", edgeIDs[i].fIndices[j]);
    printf(" : %d\n", edgeIDs[i].seq);  */
    if (edgeIDs[i].ID == NULL) {
      printf(" EGADS Error: edgeID %d is NULL (EG_mapBody)!\n", i+1);
      for (j = 0; j < len; j++) {
        if (edgeIDs[j].fIndices != NULL) EG_free(edgeIDs[j].fIndices);
        if (edgeIDs[j].ID       != NULL) EG_free(edgeIDs[j].ID);
      }
      EG_free(edgeIDs);
      return EGADS_MALLOC;
    }
    if (edgeIDs[i].fIndices != NULL) EG_free(edgeIDs[i].fIndices);
/*  printf(" Edge %d: %s\n", i+1, edgeIDs[i].ID);  */
  }

  *IDs = edgeIDs;
  return EGADS_SUCCESS;
}


static int
EG_sameAttrs(int aTypes, int lens, const int *ints, const double *reals,
             const char *strs, int aTyped, int lend, const int *intd,
             const double *reald, const char *strd)
{
  int i;

  if (aTypes != aTyped) return EGADS_OUTSIDE;

  if (aTypes == ATTRINT) {
    if (lens != lend) return EGADS_OUTSIDE;
    for (i = 0; i < lens; i++)
      if (ints[i] != intd[i]) return EGADS_OUTSIDE;
  } else if ((aTypes == ATTRREAL) || (aTypes == ATTRCSYS)) {
    if (lens != lend) return EGADS_OUTSIDE;
    for (i = 0; i < lens; i++)
      if (reals[i] != reald[i]) return EGADS_OUTSIDE;
  } else {
    if (strcmp(strs,strd) != 0) return EGADS_OUTSIDE;
  }

  return EGADS_SUCCESS;
}


static int
EG_bodyMapping(const egObject  *sBody, const egObject *dBody, const char *fAttr,
               int **pnMap, int **peMap, int **pfMap)
{
  int          i, j, in0, in1, jn0, jn1, outLevel, stat;
  int          aTypes, aTyped, lens, lend;
  int          *nMap = NULL, *eMap = NULL, *fMap = NULL;
  edgeID       *edgeIDs, *edgeIDd;
  const int    *ints,  *intd;
  const char   *strs,  *strd;
  const double *reals, *reald;

  egadsBody *pbods = (egadsBody *) sBody->blind;
  egadsBody *pbodd = (egadsBody *) dBody->blind;
  outLevel = EG_outLevel(sBody);

  if (fAttr == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Face Attribute (EG_mapBody)!\n");
    return EGADS_NONAME;
  }

  // Faces first
  fMap = (int *) EG_alloc(pbods->faces.map.Extent()*sizeof(int));
  if (fMap == NULL) return EGADS_MALLOC;
  for (i = 0; i < pbods->faces.map.Extent(); i++) fMap[i] = 0;
  for (i = 0; i < pbods->faces.map.Extent(); i++) {
    stat = EG_attributeRet(pbods->faces.objs[i], fAttr, &aTypes, &lens,
                           &ints, &reals, &strs);
    if (stat == EGADS_NOTFOUND) {
      if (outLevel > 0)
        printf(" EGADS Error: Face Attr Not in Face %d (EG_mapBody)!\n",
               i+1);
      EG_free(fMap);
      return stat;
    }
    for (j = 0; j < pbodd->faces.map.Extent(); j++) {
      stat = EG_attributeRet(pbodd->faces.objs[j], fAttr, &aTyped, &lend,
                             &intd, &reald, &strd);
      if (stat != EGADS_SUCCESS) continue;
      if (EG_sameAttrs(aTypes, lens, ints, reals, strs,
                       aTyped, lend, intd, reald, strd) == EGADS_SUCCESS) {
        if (fMap[i] == 0) {
          fMap[i] = j+1;
        } else {
          if (outLevel > 0)
            printf(" EGADS Error: Face Attr Multip in %d (EG_mapBody)!\n",
                   i+1);
          EG_free(fMap);
          return EGADS_ATTRERR;
        }
      }
    }
  }
/*
  for (i = 0; i < pbods->faces.map.Extent(); i++)
    printf("  fMap %d: %d\n", i+1, fMap[i]);
 */
  for (i = 0; i < pbods->faces.map.Extent(); i++)
    if (fMap[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: No Mapping for src Face %d (EG_mapBody)!\n",
               i+1);
      EG_free(fMap);
      return EGADS_ATTRERR;
    }

  // now Edges & Nodes
  eMap = (int *) EG_alloc(pbods->edges.map.Extent()*sizeof(int));
  if (eMap == NULL) {
    EG_free(fMap);
    return EGADS_MALLOC;
  }
  for (i = 0; i < pbods->edges.map.Extent(); i++) eMap[i] = 0;
  nMap = (int *) EG_alloc(pbods->nodes.map.Extent()*sizeof(int));
  if (nMap == NULL) {
    EG_free(eMap);
    EG_free(fMap);
    return EGADS_MALLOC;
  }
  for (i = 0; i < pbods->nodes.map.Extent(); i++) nMap[i] = 0;

  // make the Edge IDs
  stat = EG_getEdgeIDs(pbods, fAttr, &edgeIDs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getEdgeIDs src ret = %d (EG_mapBody)!\n", stat);
    EG_free(nMap);
    EG_free(eMap);
    EG_free(fMap);
    return stat;
  }
  stat = EG_getEdgeIDs(pbodd, fAttr, &edgeIDd);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getEdgeIDs dst ret = %d (EG_mapBody)!\n", stat);
    for (i = 0; i < pbods->edges.map.Extent(); i++)
      if (edgeIDs[i].ID != NULL) EG_free(edgeIDs[i].ID);
    EG_free(edgeIDs);
    EG_free(nMap);
    EG_free(eMap);
    EG_free(fMap);
    return stat;
  }

  for (i = 0; i < pbods->edges.map.Extent(); i++) {
    egadsEdge     *pedge = (egadsEdge *) pbods->edges.objs[i]->blind;
    egadsNode     *pnode = (egadsNode *) pedge->nodes[0]->blind;
    TopoDS_Vertex vert   = TopoDS::Vertex(pnode->node);
    in0 = in1            = pbods->nodes.map.FindIndex(vert);
    if (in0 <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot find N0 in Edge %d (EG_mapBody)!\n",
               i+1);
      for (j = 0; j < pbods->edges.map.Extent(); j++)
        if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
      EG_free(edgeIDs);
      for (j = 0; j < pbodd->edges.map.Extent(); j++)
        if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
      EG_free(edgeIDd);
      EG_free(nMap);
      EG_free(eMap);
      EG_free(fMap);
      return EGADS_TOPOERR;
    }
    if (pbods->edges.objs[i]->mtype == TWONODE) {
      pnode = (egadsNode *) pedge->nodes[1]->blind;
      vert  = TopoDS::Vertex(pnode->node);
      in1   = pbods->nodes.map.FindIndex(vert);
      if (in1 <= 0) {
        if (outLevel > 0)
          printf(" EGADS Error: Cannot find N1 in Edge %d (EG_mapBody)!\n",
                 i+1);
        for (j = 0; j < pbods->edges.map.Extent(); j++)
          if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
        EG_free(edgeIDs);
        for (j = 0; j < pbodd->edges.map.Extent(); j++)
          if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
        EG_free(edgeIDd);
        EG_free(nMap);
        EG_free(eMap);
        EG_free(fMap);
        return EGADS_TOPOERR;
      }
    }
    for (j = 0; j < pbodd->edges.map.Extent(); j++) {
      if (strcmp(edgeIDs[i].ID, edgeIDd[j].ID) == 0) {
        if (eMap[i] == 0) {
          if (pbods->edges.objs[i]->mtype != pbodd->edges.objs[j]->mtype) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge Type Mismatch %d %d (EG_mapBody)!\n",
                     i+1, j+1);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_TOPOERR;
          }
          eMap[i]   = j+1;
          pedge     = (egadsEdge *) pbodd->edges.objs[j]->blind;
          pnode     = (egadsNode *) pedge->nodes[0]->blind;
          vert      = TopoDS::Vertex(pnode->node);
          jn0 = jn1 = pbodd->nodes.map.FindIndex(vert);
          if (jn0 <= 0) {
            if (outLevel > 0)
              printf(" EGADS Error: Cannot find N0 dstE %d (EG_mapBody)!\n",
                     j+1);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_TOPOERR;
          }
          if (pbodd->edges.objs[j]->mtype == TWONODE) {
            pnode = (egadsNode *) pedge->nodes[1]->blind;
            vert  = TopoDS::Vertex(pnode->node);
            jn1   = pbodd->nodes.map.FindIndex(vert);
            if (jn1 <= 0) {
              if (outLevel > 0)
                printf(" EGADS Error: Cannot find N1 dstE %d (EG_mapBody)!\n",
                       j+1);
              for (j = 0; j < pbods->edges.map.Extent(); j++)
                if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
              EG_free(edgeIDs);
              for (j = 0; j < pbodd->edges.map.Extent(); j++)
                if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
              EG_free(edgeIDd);
              EG_free(nMap);
              EG_free(eMap);
              EG_free(fMap);
              return EGADS_TOPOERR;
            }
          }
          if (nMap[in0-1] == 0) {
            nMap[in0-1] = jn0;
          } else if (nMap[in0-1] != jn0) {
            if (outLevel > 0)
              printf(" EGADS Error: Node attr Multip %d: %d %d (EG_mapBody)!\n",
                     in0, nMap[in0-1], jn0);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_ATTRERR;
          }
          if (nMap[in1-1] == 0) {
            nMap[in1-1] = jn1;
          } else if (nMap[in1-1] != jn1) {
            if (outLevel > 0)
              printf(" EGADS Error: Node Attr Multip %d: %d %d (EG_mapBody)!\n",
                     in1, nMap[in1-1], jn1);
            for (j = 0; j < pbods->edges.map.Extent(); j++)
              if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
            EG_free(edgeIDs);
            for (j = 0; j < pbodd->edges.map.Extent(); j++)
              if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
            EG_free(edgeIDd);
            EG_free(nMap);
            EG_free(eMap);
            EG_free(fMap);
            return EGADS_ATTRERR;
          }
        } else {
          if (outLevel > 0)
            printf(" EGADS Error: Edge Attr Multip %d: %d %d (EG_mapBody)!\n",
                   i+1, eMap[i], j+1);
          for (j = 0; j < pbods->edges.map.Extent(); j++)
            if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
          EG_free(edgeIDs);
          for (j = 0; j < pbodd->edges.map.Extent(); j++)
            if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
          EG_free(edgeIDd);
          EG_free(nMap);
          EG_free(eMap);
          EG_free(fMap);
          return EGADS_ATTRERR;
        }
      }
    }
  }
/*
  for (i = 0; i < pbods->edges.map.Extent(); i++) {
    j = eMap[i] - 1;
    printf("  eMap %d: %d  %d %d\n", i+1, eMap[i],
    pbods->edges.objs[i]->mtype, pbodd->edges.objs[j]->mtype);
  }
 */
  for (i = 0; i < pbods->edges.map.Extent(); i++)
    if (eMap[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: No Mapping for src Edge %d (EG_mapBody)!\n",
               i+1);
      for (j = 0; j < pbods->edges.map.Extent(); j++)
        if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
      EG_free(edgeIDs);
      for (j = 0; j < pbodd->edges.map.Extent(); j++)
        if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
      EG_free(edgeIDd);
      EG_free(nMap);
      EG_free(eMap);
      EG_free(fMap);
      return EGADS_ATTRERR;
    }
/*
  for (i = 0; i < pbods->nodes.map.Extent(); i++)
    printf("  nMap %d: %d\n", i+1, nMap[i]);
 */
  for (i = 0; i < pbods->nodes.map.Extent(); i++)
    if (nMap[i] == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: No Mapping for src Node %d (EG_mapBody)!\n",
               i+1);
      for (j = 0; j < pbods->edges.map.Extent(); j++)
        if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
      EG_free(edgeIDs);
      for (j = 0; j < pbodd->edges.map.Extent(); j++)
        if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
      EG_free(edgeIDd);
      EG_free(nMap);
      EG_free(eMap);
      EG_free(fMap);
      return EGADS_ATTRERR;
    }

  for (j = 0; j < pbods->edges.map.Extent(); j++)
    if (edgeIDs[j].ID != NULL) EG_free(edgeIDs[j].ID);
  EG_free(edgeIDs);
  for (j = 0; j < pbodd->edges.map.Extent(); j++)
    if (edgeIDd[j].ID != NULL) EG_free(edgeIDd[j].ID);
  EG_free(edgeIDd);

  *pnMap = nMap;
  *peMap = eMap;
  *pfMap = fMap;
  return EGADS_SUCCESS;
}


int
EG_mapBody(const egObject  *sBody, const egObject *dBody, const char *fAttr,
                 egObject **mBody)
{
  int outLevel, stat, *nMap = NULL, *eMap = NULL, *fMap = NULL, exacTopo = 0;
#ifdef MAPBSPLINE
  int i, j, k, cnt;
#endif

  *mBody = NULL;
  if  (sBody == NULL)                return EGADS_NULLOBJ;
  if  (sBody->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (sBody->oclass != BODY)        return EGADS_NOTBODY;
  if  (sBody->blind == NULL)         return EGADS_NODATA;
  if ((sBody->mtype != FACEBODY)  &&
      (sBody->mtype != SHEETBODY) &&
      (sBody->mtype != SOLIDBODY))   return EGADS_TOPOERR;
  if  (dBody == NULL)                return EGADS_NULLOBJ;
  if  (dBody->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if  (dBody->oclass != BODY)        return EGADS_NOTBODY;
  if  (dBody->blind == NULL)         return EGADS_NODATA;
  if  (dBody->mtype != sBody->mtype) return EGADS_TOPOERR;
  egadsBody *pbodd = (egadsBody *) dBody->blind;
  outLevel = EG_outLevel(sBody);

  stat = EG_sameBodyTopo(sBody, dBody);
  if (stat == EGADS_TOPOCNT) return stat;
  if (stat == EGADS_SUCCESS) exacTopo = 1;
  if (sBody->mtype == FACEBODY) {
    if (exacTopo != 1) return EGADS_TOPOERR;
    /* fix spline surface? */
    goto fixBSplSur;
  }

  // differing topology -- make mappings
  if (exacTopo == 0) {
    stat = EG_bodyMapping(sBody, dBody, fAttr, &nMap, &eMap, &fMap);
    if (stat != EGADS_SUCCESS) return stat;
  }

  // look at BSpline surfaces for potential Face updates

fixBSplSur:
#ifdef MAPBSPLINE
  egadsBody *pbods = (egadsBody *) sBody->blind;
  int *bMap = (int *) EG_alloc(pbods->faces.map.Extent()*sizeof(int));
  if (bMap == NULL) {
    if (nMap != NULL) EG_free(nMap);
    if (eMap != NULL) EG_free(eMap);
    if (fMap != NULL) EG_free(fMap);
    return EGADS_MALLOC;
  }
  for (cnt = j = 0; j < pbods->faces.map.Extent(); j++) {
    int      oclass, mtype, OK, *senses;
    double   tr, range[4];
    egObject *geom, **child;

    bMap[j] = 0;
    k       = j;
    if (fMap != NULL) k = fMap[j] - 1;

    /* check BSpline surfaces */
    OK   = 0;
    stat = EG_getTopology(pbodd->faces.objs[k], &geom, &oclass, &mtype, range,
                          &i, &child, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTopology dFace %d = %d (EG_mapBody)!\n",
               k+1, stat);
      EG_free(bMap);
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return stat;
    }
    if (geom->mtype == BSPLINE) {
      int      *ivec = NULL, *ivecb = NULL;
      double   *rvec = NULL, *rvecb = NULL, sknot, sknotb;
      egObject *rGeom, *geomm;

      OK   = -1;
      stat = EG_getGeometry(geom, &oclass, &mtype, &rGeom, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getGeometry dFace %d = %d (EG_mapBody)!\n",
                 k+1, stat);
        EG_free(bMap);
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        return stat;
      }
      tr   = 0.0;
      stat = EG_getTopology(pbods->faces.objs[j], &geomm, &oclass, &mtype,
                            range, &i, &child, &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTopology sFace %d = %d (EG_mapBody)!\n",
                 j+1, stat);
        EG_free(ivec);
        EG_free(rvec);
        EG_free(bMap);
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        return stat;
      }
      if (geomm->mtype != BSPLINE) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d BSpline Mismatch (EG_mapBody)!\n",
                 j+1);
        EG_free(ivec);
        EG_free(rvec);
        EG_free(bMap);
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        return EGADS_GEOMERR;
      }
      stat = EG_getGeometry(geomm, &oclass, &mtype, &rGeom, &ivecb, &rvecb);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getGeometry sFace %d = %d (EG_mapBody)!\n",
                 j+1, stat);
        EG_free(ivec);
        EG_free(rvec);
        EG_free(bMap);
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        return stat;
      }
      if ((ivec[3] == ivecb[3]) && (ivec[6] == ivecb[6])) {
        for (i = 0; i < ivec[3]; i++) {
          sknot  = (rvec[i] -rvec[0]) /(rvec[ivec[3]-1] -rvec[0]);
          sknotb = (rvecb[i]-rvecb[0])/(rvecb[ivec[3]-1]-rvecb[0]);
          if (fabs(sknot-sknotb) > KNDIFF) {
            if (fabs(sknot-sknotb) > tr) tr = fabs(sknot-sknotb);
            OK = 1;
          }
        }
        for (i = ivec[3]; i < ivec[3]+ivec[6]; i++) {
          sknot  = (rvec[i] -rvec[ivec[3]]) /
                   (rvec[ivec[3]+ivec[6]-1] -rvec[ivec[3]]);
          sknotb = (rvecb[i]-rvecb[ivec[3]])/
                   (rvecb[ivec[3]+ivec[6]-1]-rvecb[ivec[3]]);
          if (fabs(sknot-sknotb) > KNDIFF) {
            if (fabs(sknot-sknotb) > tr) tr = fabs(sknot-sknotb);
            if (OK == -1) {
              OK = 2;
            } else {
              if (OK == 1) OK = 3;
            }
          }
        }
        if (OK == -1) OK = 0;
      } else {
        if (outLevel > 0)
          printf(" EGADS Warning: Face %d %d BSpline #knots (EG_mapBody)!\n",
                 j+1, k+1);
      }
      EG_free(ivecb);
      EG_free(rvecb);
      EG_free(ivec);
      EG_free(rvec);
      if (OK == -1) {
        EG_free(bMap);
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        return EGADS_GEOMERR;
      } else if (OK > 0) {
        bMap[j] = 1;
        cnt++;
        if (outLevel > 0)
          printf(" EGADS Info: Face %d BSpline knots UV=%d  %le (EG_mapBody)!\n",
                 j+1, OK, tr);
      }
    }
  }

  // need to set knots and make new body

  if (cnt != 0) {
    egObject *model;

    i = pbodd->faces.map.Extent();
    egObject **faces = (egObject **) EG_alloc(i*sizeof(egObject *));
    if (faces == NULL) {
      EG_free(bMap);
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return EGADS_MALLOC;
    }
    for (j = 0; j < pbods->faces.map.Extent(); j++) {
      k = j;
      if (fMap != NULL) k = fMap[j] - 1;

      faces[k] = pbodd->faces.objs[k];
      if (bMap[j] != 0) {
        stat = EG_mapKnots(pbods->faces.objs[j], pbodd->faces.objs[k],
                           &faces[k]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_mapKnots sFace %d = %d (EG_mapBody)!\n",
                   j+1, stat);
          EG_free(faces);
          EG_free(bMap);
          if (nMap != NULL) EG_free(nMap);
          if (eMap != NULL) EG_free(eMap);
          if (fMap != NULL) EG_free(fMap);
          return stat;
        }
      }
    }

    stat = EG_sewFaces(i, (const egObject **) faces, 0.0, 0, &model);
    for (j = 0; j < i; j++) {
      k = j;
      if (fMap != NULL) k = fMap[j] - 1;
      if (bMap[j] != 0) EG_deleteObject(faces[k]);
    }
    EG_free(faces);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_sewFaces = %d (EG_mapBody)!\n", stat);
      EG_free(bMap);
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return stat;
    }
    egadsModel *pmdl = (egadsModel *) model->blind;
    if (pmdl->nbody != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_sewFaces returns %d Bodies (EG_mapBody)!\n",
               pmdl->nbody);
      EG_deleteObject(model);
      EG_free(bMap);
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return EGADS_CONSTERR;
    }
    stat = EG_copyObject(pmdl->bodies[0], NULL, mBody);
    EG_deleteObject(model);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_copyObject returns %d (EG_mapBody)!\n", stat);
      EG_free(bMap);
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return stat;
    }

    // is the mapping still OK?
    if (exacTopo == 1) {
      stat = EG_sameBodyTopo(sBody, *mBody);
      /* remap? */
      if (stat == EGADS_TOPOCNT) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_sameBodyTopo returns %d (EG_mapBody)!\n",
                 stat);
        EG_free(bMap);
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        return stat;
      } else if (stat != EGADS_SUCCESS) {
        exacTopo = 0;
        if (nMap != NULL) EG_free(nMap);
        if (eMap != NULL) EG_free(eMap);
        if (fMap != NULL) EG_free(fMap);
        stat = EG_bodyMapping(sBody, *mBody, fAttr, &nMap, &eMap, &fMap);
        if (stat != EGADS_SUCCESS) {
          EG_free(bMap);
          return stat;
        }
      }
    } else {
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      stat = EG_bodyMapping(sBody, *mBody, fAttr, &nMap, &eMap, &fMap);
      if (stat != EGADS_SUCCESS) {
        EG_free(bMap);
        return stat;
      }
    }
  }
  EG_free(bMap);
#endif

  // can we just use the destination Body?
  if ((*mBody == NULL) && (exacTopo == 1)) return EGADS_SUCCESS;

  // make a new body if not already done
  if (*mBody == NULL) {
    stat = EG_copyObject(dBody, NULL, mBody);
    if (stat != EGADS_SUCCESS) {
      if (nMap != NULL) EG_free(nMap);
      if (eMap != NULL) EG_free(eMap);
      if (fMap != NULL) EG_free(fMap);
      return stat;
    }
  }

  // put the mapping attributes on the body
  if (nMap != NULL) {
    stat = EG_attributeAdd(*mBody, ".nMap", ATTRINT, pbodd->nodes.map.Extent(),
                           nMap, NULL, NULL);
    EG_free(nMap);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(*mBody);
      *mBody = NULL;
      if (outLevel > 0)
        printf(" EGADS Warning: Attribute Node Mapping is %d (EG_mapBody)!\n",
               stat);
      EG_free(eMap);
      EG_free(fMap);
      return stat;
    }
  }
  if (eMap != NULL) {
    stat = EG_attributeAdd(*mBody, ".eMap", ATTRINT, pbodd->edges.map.Extent(),
                           eMap, NULL, NULL);
    EG_free(eMap);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(*mBody);
      *mBody = NULL;
      if (outLevel > 0)
        printf(" EGADS Warning: Attribute Edge Mapping is %d (EG_mapBody)!\n",
               stat);
      EG_free(fMap);
      return stat;
    }
  }
  if (fMap != NULL) {
    stat = EG_attributeAdd(*mBody, ".fMap", ATTRINT, pbodd->faces.map.Extent(),
                           fMap, NULL, NULL);
    EG_free(fMap);
    if (stat != EGADS_SUCCESS) {
      EG_deleteObject(*mBody);
      *mBody = NULL;
      if (outLevel > 0)
        printf(" EGADS Warning: Attribute Face Mapping is %d (EG_mapBody)!\n",
               stat);
      return stat;
    }
  }

  return EGADS_SUCCESS;
}
