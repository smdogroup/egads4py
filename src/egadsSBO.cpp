/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Solid Boolean Operator Functions
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

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"

#include "egadsSBO.h"


#define PARAMACC 1.0e-4         // parameter accuracy


extern void EG_nearestCurve(Handle_Geom_Curve hCurve, const double *coor,
                            double tmin, double tmax, int flag,
                            double *t, double *xyz);
extern void EG_checkStatus(const Handle_BRepCheck_Result tResult);


class egadsInterLoop
{
public:
  TopoDS_Wire oldLoop;
  TopoDS_Face oldFace;
  TopoDS_Wire newLoop;
};



#ifdef DEBUG
static void showPCurve(Handle(Geom2d_Curve) hCurve)
{
  Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) {
    printf("  PCurve is a Line!\n");
    gp_Dir2d direct = hLine->Direction();
    gp_Pnt2d locat  = hLine->Location();
    printf("      %lf %lf   %lf %lf\n",
           locat.X(), locat.Y(), direct.X(), direct.Y());
    return;
  }
  
  Handle(Geom2d_Circle) hCircle = Handle(Geom2d_Circle)::DownCast(hCurve);
  if (!hCircle.IsNull()) {
    printf("   PCurve is a Circle!\n");
    return;
  }
  
  Handle(Geom2d_Ellipse) hEllipse = Handle(Geom2d_Ellipse)::DownCast(hCurve);
  if (!hEllipse.IsNull()) {
    printf("   PCurve is a Circle!\n");
    return;
  }
  
  Handle(Geom2d_Parabola) hParabola = Handle(Geom2d_Parabola)::DownCast(hCurve);
  if (!hParabola.IsNull()) {
    printf("   PCurve is a Parabola!\n");
    return;
  }
  
  Handle(Geom2d_Hyperbola) hHyperbola = Handle(Geom2d_Hyperbola)::DownCast(hCurve);
  if (!hHyperbola.IsNull()) {
    printf("   PCurve is a Hyperbola!\n");
    return;
  }

  Handle(Geom2d_BezierCurve) hBezier = Handle(Geom2d_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) {
    printf("   PCurve is a Bezier!\n");
    return;
  }

  Handle(Geom2d_BSplineCurve) hBSpline = Handle(Geom2d_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) {
    printf("   PCurve is a BSpline!\n");
    return;
  }
  
  Handle(Geom2d_TrimmedCurve) hTrim = Handle(Geom2d_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) {
    printf("   PCurve is a Trim!\n");
    return;
  }
  
  Handle(Geom2d_OffsetCurve) hOffst = Handle(Geom2d_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) {
    printf("   PCurve is an Offset!\n");
    return;
  }
  
  printf("   PCurve is of unknown type!\n");
}


static void showCurve(Handle(Geom_Curve) hCurve)
{
  Handle(Geom_Line) hLine = Handle(Geom_Line)::DownCast(hCurve);
  if (!hLine.IsNull()) {
    printf("  Curve is a Line!\n");
    return;
  }
  
  Handle(Geom_Circle) hCircle = Handle(Geom_Circle)::DownCast(hCurve);
  if (!hCircle.IsNull()) {
    printf("   Curve is a Circle!\n");
    gp_Circ circle = hCircle->Circ();
    gp_Pnt  locate = circle.Location();
    gp_Ax1  xaxis  = circle.XAxis();
    gp_Ax1  yaxis  = circle.YAxis();
    printf("         %lf %lf %lf  Radius = %lf\n",
           locate.X(), locate.Y(), locate.Z(), circle.Radius());
    printf("         %lf %lf %lf   %lf %lf %lf\n",
           xaxis.Direction().X(), xaxis.Direction().Y(), xaxis.Direction().Z(),
           yaxis.Direction().X(), yaxis.Direction().Y(), yaxis.Direction().Z());
    return;
  }
  
  Handle(Geom_Ellipse) hEllipse = Handle(Geom_Ellipse)::DownCast(hCurve);
  if (!hEllipse.IsNull()) {
    printf("   Curve is a Circle!\n");
    return;
  }
  
  Handle(Geom_Parabola) hParabola = Handle(Geom_Parabola)::DownCast(hCurve);
  if (!hParabola.IsNull()) {
    printf("   Curve is a Parabola!\n");
    return;
  }
  
  Handle(Geom_Hyperbola) hHyperbola = Handle(Geom_Hyperbola)::DownCast(hCurve);
  if (!hHyperbola.IsNull()) {
    printf("   Curve is a Hyperbola!\n");
    return;
  }
  
  Handle(Geom_BezierCurve) hBezier = Handle(Geom_BezierCurve)::DownCast(hCurve);
  if (!hBezier.IsNull()) {
    printf("   Curve is a Bezier!\n");
    return;
  }
  
  Handle(Geom_BSplineCurve) hBSpline = Handle(Geom_BSplineCurve)::DownCast(hCurve);
  if (!hBSpline.IsNull()) {
    printf("   Curve is a BSpline!\n");
    return;
  }
  
  Handle(Geom_TrimmedCurve) hTrim = Handle(Geom_TrimmedCurve)::DownCast(hCurve);
  if (!hTrim.IsNull()) {
    printf("   Curve is a Trim!\n");
    return;
  }
  
  Handle(Geom_OffsetCurve) hOffst = Handle(Geom_OffsetCurve)::DownCast(hCurve);
  if (!hOffst.IsNull()) {
    printf("   Curve is an Offset!\n");
    return;
  }
  
  printf("   Curve is of unknown type!\n");
}
#endif


static bool sameCurve(Handle(Geom_Curve) hCurv1, Handle(Geom_Curve) hCurv2)
{
  
  // line
  Handle(Geom_Line) hLin1 = Handle(Geom_Line)::DownCast(hCurv1);
  Handle(Geom_Line) hLin2 = Handle(Geom_Line)::DownCast(hCurv2);
  if (!hLin1.IsNull() && !hLin2.IsNull()) {
    gp_Lin line1   = hLin1->Lin();
    gp_Lin line2   = hLin2->Lin();
    gp_Dir direct1 = line1.Direction();
    gp_Dir direct2 = line2.Direction();
    if (!direct1.IsEqual(direct2, 0.0)) return false;
    gp_Pnt locat1  = line1.Location();
    gp_Pnt locat2  = line2.Location();
    if (!locat1.IsEqual(locat2, 0.0))   return false;
    return true;
  }
  
  // circle
  Handle(Geom_Circle) hCirc1 = Handle(Geom_Circle)::DownCast(hCurv1);
  Handle(Geom_Circle) hCirc2 = Handle(Geom_Circle)::DownCast(hCurv2);
  if (!hCirc1.IsNull() && !hCirc2.IsNull()) {
    gp_Circ circ1  = hCirc1->Circ();
    gp_Circ circ2  = hCirc2->Circ();
    if (circ1.Radius() != circ2.Radius()) return false;
    gp_Ax1 xaxis1  = circ1.XAxis();
    gp_Ax1 xaxis2  = circ2.XAxis();
    if (!xaxis1.IsParallel(xaxis2, 0.0))  return false;
    gp_Ax1 yaxis1  = circ1.YAxis();
    gp_Ax1 yaxis2  = circ2.YAxis();
    if (!yaxis1.IsParallel(yaxis2, 0.0))  return false;
    gp_Pnt locat1  = circ1.Location();
    gp_Pnt locat2  = circ2.Location();
    if (!locat1.IsEqual(locat2, 0.0))     return false;
    return true;
  }
  
  // ellipse
  Handle(Geom_Ellipse) hEllip1 = Handle(Geom_Ellipse)::DownCast(hCurv1);
  Handle(Geom_Ellipse) hEllip2 = Handle(Geom_Ellipse)::DownCast(hCurv2);
  if (!hEllip1.IsNull() && !hEllip2.IsNull()) {
    gp_Elips elips1 = hEllip1->Elips();
    gp_Elips elips2 = hEllip2->Elips();
    if (elips1.MajorRadius() != elips2.MajorRadius()) return false;
    if (elips1.MinorRadius() != elips2.MinorRadius()) return false;
    gp_Ax1 xaxis1   = elips1.XAxis();
    gp_Ax1 xaxis2   = elips2.XAxis();
    if (!xaxis1.IsParallel(xaxis2, 0.0))              return false;
    gp_Ax1 yaxis1   = elips1.YAxis();
    gp_Ax1 yaxis2   = elips2.YAxis();
    if (!yaxis1.IsParallel(yaxis2, 0.0))              return false;
    gp_Pnt locat1   = elips1.Location();
    gp_Pnt locat2   = elips2.Location();
    if (!locat1.IsEqual(locat2, 0.0))                 return false;
    return true;
  }

  // parabola
  Handle(Geom_Parabola) hParab1 = Handle(Geom_Parabola)::DownCast(hCurv1);
  Handle(Geom_Parabola) hParab2 = Handle(Geom_Parabola)::DownCast(hCurv2);
  if (!hParab1.IsNull() && !hParab2.IsNull()) {
    gp_Parab parab1 = hParab1->Parab();
    gp_Parab parab2 = hParab2->Parab();
    if (parab1.Focal() != parab2.Focal()) return false;
    gp_Ax1 xaxis1   = parab1.XAxis();
    gp_Ax1 xaxis2   = parab2.XAxis();
    if (!xaxis1.IsParallel(xaxis2, 0.0))  return false;
    gp_Ax1 yaxis1   = parab1.YAxis();
    gp_Ax1 yaxis2   = parab2.YAxis();
    if (!yaxis1.IsParallel(yaxis2, 0.0))  return false;
    gp_Pnt locat1   = parab1.Location();
    gp_Pnt locat2   = parab2.Location();
    if (!locat1.IsEqual(locat2, 0.0))     return false;
    return true;
  }
  
  // hyperbola
  Handle(Geom_Hyperbola) hHypr1 = Handle(Geom_Hyperbola)::DownCast(hCurv1);
  Handle(Geom_Hyperbola) hHypr2 = Handle(Geom_Hyperbola)::DownCast(hCurv2);
  if (!hHypr1.IsNull() && !hHypr2.IsNull()) {
    gp_Hypr hypr1  = hHypr1->Hypr();
    gp_Hypr hypr2  = hHypr2->Hypr();
    if (hypr1.MajorRadius() != hypr2.MajorRadius()) return false;
    if (hypr1.MinorRadius() != hypr2.MinorRadius()) return false;
    gp_Ax1 xaxis1  = hypr1.XAxis();
    gp_Ax1 xaxis2  = hypr2.XAxis();
    if (!xaxis1.IsParallel(xaxis2, 0.0))            return false;
    gp_Ax1 yaxis1  = hypr1.YAxis();
    gp_Ax1 yaxis2  = hypr2.YAxis();
    if (!yaxis1.IsParallel(yaxis2, 0.0))            return false;
    gp_Pnt locat1  = hypr1.Location();
    gp_Pnt locat2  = hypr2.Location();
    if (!locat1.IsEqual(locat2, 0.0))               return false;
    return true;
  }

  // bezier
  Handle(Geom_BezierCurve) hBezier1 = Handle(Geom_BezierCurve)::DownCast(hCurv1);
  Handle(Geom_BezierCurve) hBezier2 = Handle(Geom_BezierCurve)::DownCast(hCurv2);
  if (!hBezier1.IsNull() && !hBezier2.IsNull()) {
    if (hBezier1->IsRational() != hBezier2->IsRational()) return false;
    if (hBezier1->Degree()     != hBezier2->Degree())     return false;
    if (hBezier1->NbPoles()    != hBezier2->NbPoles())    return false;
    for (int i = 1; i <= hBezier1->NbPoles(); i++) {
      gp_Pnt P1 = hBezier1->Pole(i);
      gp_Pnt P2 = hBezier2->Pole(i);
      if (!P1.IsEqual(P2, 0.0))                           return false;
    }
    if (!hBezier1->IsRational())                          return true;
    for (int i = 1; i <= hBezier1->NbPoles(); i++)
      if (hBezier1->Weight(i) != hBezier2->Weight(i))     return false;
    return true;
  }

  // bspline
  Handle(Geom_BSplineCurve) hBSplin1 = Handle(Geom_BSplineCurve)::DownCast(hCurv1);
  Handle(Geom_BSplineCurve) hBSplin2 = Handle(Geom_BSplineCurve)::DownCast(hCurv2);
  if (!hBSplin1.IsNull() && !hBSplin2.IsNull()) {
    if (hBSplin1->IsRational() != hBSplin2->IsRational())         return false;
    if (hBSplin1->Degree()     != hBSplin2->Degree())             return false;
    if (hBSplin1->NbPoles()    != hBSplin2->NbPoles())            return false;
    if (hBSplin1->NbKnots()    != hBSplin2->NbKnots())            return false;
    for (int i = 1; i <= hBSplin1->NbKnots(); i++) {
      if (hBSplin1->Multiplicity(i) != hBSplin2->Multiplicity(i)) return false;
      if (hBSplin1->Knot(i)         != hBSplin2->Knot(i))         return false;
    }
    for (int i = 1; i <= hBSplin1->NbPoles(); i++) {
      gp_Pnt P1 = hBSplin1->Pole(i);
      gp_Pnt P2 = hBSplin2->Pole(i);
      if (!P1.IsEqual(P2, 0.0))                                   return false;
    }
    if (!hBSplin1->IsRational())                                  return true;
    for (int i = 1; i <= hBSplin1->NbPoles(); i++)
      if (hBSplin1->Weight(i) != hBSplin2->Weight(i))             return false;
    return true;
  }
  
  // referencing geometry -- trimmed
  Handle(Geom_TrimmedCurve) hTrim1 = Handle(Geom_TrimmedCurve)::DownCast(hCurv1);
  Handle(Geom_TrimmedCurve) hTrim2 = Handle(Geom_TrimmedCurve)::DownCast(hCurv2);
  if (!hTrim1.IsNull() && !hTrim2.IsNull()) {
    if (hTrim1->FirstParameter() != hTrim2->FirstParameter()) return false;
    if (hTrim1->LastParameter()  != hTrim2->LastParameter())  return false;
    return sameCurve(hTrim1->BasisCurve(), hTrim2->BasisCurve());
  }
  
  // referencing geometry -- offset
  Handle(Geom_OffsetCurve) hOffst1 = Handle(Geom_OffsetCurve)::DownCast(hCurv1);
  Handle(Geom_OffsetCurve) hOffst2 = Handle(Geom_OffsetCurve)::DownCast(hCurv2);
  if (!hOffst1.IsNull() && !hOffst2.IsNull()) {
    gp_Dir direct1 = hOffst1->Direction();
    gp_Dir direct2 = hOffst2->Direction();
    if (!direct1.IsEqual(direct2, 0.0)) return false;
    return sameCurve(hOffst1->BasisCurve(), hOffst2->BasisCurve());
  }

  return false;
}


static void adjPeriodics(double Uper, double umin, double umax, double &u)
{
  if (Uper == 0.0) return;
  
  if ((u+PARAMACC < umin) || (u-PARAMACC > umax))
    if (u+PARAMACC < umin) {
      if (u+Uper-PARAMACC < umax) u += Uper;
    } else {
      if (u-Uper+PARAMACC > umin) u -= Uper;
    }
}


#ifdef DEBUG
static void compareShapes(TopoDS_Shape shap1, TopoDS_Shape shap2)
{
  TopTools_IndexedMapOfShape vMap1, vMap2, eMap1, eMap2, lMap1, lMap2, fMap1, fMap2;
  TopExp::MapShapes(shap1, TopAbs_VERTEX, vMap1);
  TopExp::MapShapes(shap2, TopAbs_VERTEX, vMap2);
  TopExp::MapShapes(shap1, TopAbs_EDGE,   eMap1);
  TopExp::MapShapes(shap2, TopAbs_EDGE,   eMap2);
  TopExp::MapShapes(shap1, TopAbs_WIRE,   lMap1);
  TopExp::MapShapes(shap2, TopAbs_WIRE,   lMap2);
  TopExp::MapShapes(shap1, TopAbs_FACE,   fMap1);
  TopExp::MapShapes(shap2, TopAbs_FACE,   fMap2);

  printf("  Vertex Maps = %d %d\n", vMap1.Extent(), vMap2.Extent());
  for (int i = 1; i <= vMap2.Extent(); i++) {
    TopoDS_Vertex node2 = TopoDS::Vertex(vMap2(i));
    printf("    %d %d\n", vMap1.FindIndex(node2), i);
  }
  
  printf("  Edge   Maps = %d %d\n", eMap1.Extent(), eMap2.Extent());
  int degen2 = 0;
  for (int i = 1; i <= eMap2.Extent(); i++) {
    TopoDS_Edge edge2 = TopoDS::Edge(eMap2(i));
    if (BRep_Tool::Degenerated(edge2)) {
      degen2++;
    } else {
      printf("    %d %d\n", eMap1.FindIndex(edge2), i);
    }
  }
  int degen1 = 0;
  for (int i = 1; i <= eMap1.Extent(); i++) {
    TopoDS_Edge edge1 = TopoDS::Edge(eMap1(i));
    if (BRep_Tool::Degenerated(edge1)) degen1++;
  }
  printf("  Edge Degens = %d %d\n", degen1, degen2);
  
  printf("  Loop   Maps = %d %d\n", lMap1.Extent(), lMap2.Extent());
  for (int i = 1; i <= lMap2.Extent(); i++) {
    TopoDS_Wire loop2 = TopoDS::Wire(lMap2(i));
    printf("    %d %d\n", lMap1.FindIndex(loop2), i);
  }
  
  printf("  Face   Maps = %d %d\n", fMap1.Extent(), fMap2.Extent());
  for (int i = 1; i <= fMap2.Extent(); i++) {
    TopoDS_Face face2 = TopoDS::Face(fMap2(i));
    printf("    %d %d\n", fMap1.FindIndex(face2), i);
  }
}
#endif


static void parseLoops(TopoDS_Shape shape)
{
  TopTools_IndexedMapOfShape lmap, emap, nmap;
  TopExp::MapShapes(shape, TopAbs_WIRE,   lmap);
  TopExp::MapShapes(shape, TopAbs_EDGE,   emap);
  TopExp::MapShapes(shape, TopAbs_VERTEX, nmap);
  for (int i = 1; i <= nmap.Extent(); i++) {
    TopoDS_Vertex vert = TopoDS::Vertex(nmap(i));
    gp_Pnt pv = BRep_Tool::Pnt(vert);
    printf("   Node %d: %lf %lf %lf   %le\n", i, pv.X(), pv.Y(), pv.Z(),
           BRep_Tool::Tolerance(vert));
  }
  for (int i = 1; i <= lmap.Extent(); i++) {
    TopoDS_Wire loop = TopoDS::Wire(lmap(i));
    printf("   Loop #%d:\n", i);
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(loop); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Edge edge = TopoDS::Edge(ExpWE.Current());
      TopoDS_Vertex V1, V2;
      TopExp::Vertices(edge, V1, V2, Standard_True);
      if (BRep_Tool::Degenerated(edge)) {
        printf("      Edge %d: Degenerate @ %d\n", emap.FindIndex(edge),
               nmap.FindIndex(V1));
      } else {
        double frst, last;
        BRep_Tool::Range(edge, frst, last);
        int sense = 1;
        if (edge.Orientation() == TopAbs_REVERSED) sense = -1;
        printf("      Edge %d: %lf - %lf   %d %d   %d   %le\n",
               emap.FindIndex(edge), frst, last, nmap.FindIndex(V1),
               nmap.FindIndex(V2), sense, BRep_Tool::Tolerance(edge));
      }
    }
  }
}


static TopoDS_Shape
lookupShapeByIndex(TopoDS_Shape shap1, TopoDS_Shape shap2,
                   const TopAbs_ShapeEnum sType, TopoDS_Shape find)
{
  TopTools_IndexedMapOfShape map1, map2;
  TopExp::MapShapes(shap1, sType, map1);
  TopExp::MapShapes(shap2, sType, map2);
  
  int k = map1.FindIndex(find);
  if (k > 0) return map2(k);
  
  printf("  ERROR: lookupShapeByIndex -- Shape not found!\n");
  TopoDS_Shape shape;
  return shape;
}


static TopoDS_Face getFace(TopoDS_Shape shape, TopoDS_Wire loop)
{
  TopTools_IndexedMapOfShape fMap;
  TopExp::MapShapes(shape, TopAbs_FACE, fMap);
  for (int i = 1; i <= fMap.Extent(); i++) {
    TopTools_IndexedMapOfShape lMap;
    TopExp::MapShapes(fMap(i), TopAbs_WIRE, lMap);
    for (int j = 1; j <= lMap.Extent(); j++) {
      TopoDS_Wire wire = TopoDS::Wire(lMap(j));
      if (loop.IsSame(wire)) return TopoDS::Face(fMap(i));
    }
  }
  
  printf("  ERROR: getFace -- Face not found!\n");
  TopoDS_Face face;
  return face;
}


static void setPCurve(TopoDS_Edge &edge, TopoDS_Face Face)
{
  double pfrst, plast, efrst, elast, toler;

  BRep_Builder Builder;
  toler = BRep_Tool::Tolerance(edge);
  Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
  Handle(Geom_Plane)   hPlane   = Handle(Geom_Plane)::DownCast(hSurface);
  if (!hPlane.IsNull()) return;
  Handle(Geom_Curve)   hCurve   = BRep_Tool::Curve(edge, efrst, elast);
  Handle(Geom2d_Curve) hCurv2d  = BRep_Tool::CurveOnSurface(edge,  Face,
                                                            pfrst, plast);
  if (hCurv2d.IsNull()) {
/*  printf("     Ranges: %lf %lf   %lf %lf", efrst, elast, pfrst, plast);
    if (hCurve  == NULL) printf("     C = NULL!");
    if (hCurv2d == NULL) printf("    pC = NULL!");
    printf("\n");  */
    try {
      hCurv2d = GeomProjLib::Curve2d(hCurve, efrst, elast, hSurface, toler);
      if (hCurv2d.IsNull()) {
        printf(" EGADS Error: Cannot get PCurve (setPCurve)!\n");
        return;
      }
    }
    catch (Standard_Failure)
    {
      printf(" EGADS Info: Geometry Creation Error (setPCurve)!\n");
      return;
    }
    catch (...)
    {
      printf(" EGADS Info: Geometry Creation ERROR (setPCurve)!\n");
      return;
    }
  }
  Builder.UpdateEdge(edge, hCurv2d, Face, toler);
  Builder.SameParameter(edge, Standard_True);
  Builder.SameRange(edge, Standard_True);
  BRepLib::SameParameter(Face);
  
  // check for out-of-range
  if (hCurve->IsPeriodic()) {
    gp_Pnt2d p2f, p2l;
    double   u1, u2, v1, v2, t1, t2, ts, umin, umax, vmin, vmax, per;
  
    BRepTools::UVBounds(Face, umin, umax, vmin, vmax);
    BRep_Tool::Range(edge, t1, t2);
    if (edge.Orientation() == TopAbs_REVERSED) {
      ts = t1;
      t1 = t2;
      t2 = ts;
    }
    per = hCurve->Period();
    BRepAdaptor_Curve2d Curve2d(edge, Face);
    Curve2d.D0(t1, p2f);
    Curve2d.D0(t2, p2l);
    u1 = p2f.X();
    u2 = p2l.X();
    v1 = p2f.Y();
    v2 = p2l.Y();
    if ((u1 > umin-PARAMACC) && (u1 < umax+PARAMACC) && (u2 > umin-PARAMACC) &&
        (u2 < umax+PARAMACC) && (v1 > vmin-PARAMACC) && (v1 < vmax+PARAMACC) &&
        (v2 > vmin-PARAMACC) && (v2 < vmax+PARAMACC)) return;
    printf("  PCurve limits out of range -- u = [%lf %lf] v = [%lf %lf]!\n",
           umin, umax, vmin, vmax);
    printf("  EndPts = %lf %lf    %lf %lf\n", u1, v1, u2, v2);
    Handle(Geom2d_Line) hLine = Handle(Geom2d_Line)::DownCast(hCurv2d);
    if (!hLine.IsNull()) {
      gp_Dir2d direct = hLine->Direction();
      gp_Pnt2d locat  = hLine->Location();
      printf("  Fixing Linear PCurve  %lf %lf  %lf %lf  per = %lf!\n",
             locat.X(), locat.Y(), direct.X(), direct.Y(), per);
      if (u2 < u1) {
        ts = u1;
        u1 = u2;
        u2 = ts;
      }
      if (u1+PARAMACC < umin) locat.SetX(locat.X()+per);
      if (u2-PARAMACC > umax) locat.SetX(locat.X()-per);
      if (v2 < v1) {
        ts = v1;
        v1 = v2;
        v2 = ts;
      }
      if (v1+PARAMACC < vmin) locat.SetY(locat.Y()+per);
      if (v2-PARAMACC > vmax) locat.SetY(locat.Y()-per);
      Handle(Geom2d_Curve) hCurv2dFix = new Geom2d_Line(locat, direct);
      Builder.UpdateEdge(edge, hCurv2dFix, Face, toler);
      Builder.SameParameter(edge, Standard_True);
      Builder.SameRange(edge, Standard_True);
      BRepLib::SameParameter(Face);
    }
  }

}


static void fixPCurves(TopoDS_Face &Face)
{
  TopTools_IndexedMapOfShape MapL;
  TopExp::MapShapes(Face, TopAbs_WIRE, MapL);
  for (int iloop = 1; iloop <= MapL.Extent(); iloop++) {
    TopoDS_Wire loop = TopoDS::Wire(MapL(iloop));
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(loop); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Edge edge = TopoDS::Edge(ExpWE.Current());
      if (BRep_Tool::Degenerated(edge)) continue;
      setPCurve(edge, Face);
    }
  }
  BRepLib::SameParameter(Face);

  BRepCheck_Analyzer fCheck(Face);
  if (!fCheck.IsValid()) {
    printf("   Face is still not OK!\n");
    EG_checkStatus(fCheck.Result(Face));
    parseLoops(Face);

    // try to fix the fault
    Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(Face);
    sfs->Perform();
    TopoDS_Shape fixedFace = sfs->Shape();
    if (fixedFace.IsNull()) {
      printf(" EGADS Info: Invalid Face w/ NULL Fix (fixPCurves)!\n");
      return;
    }
    BRepCheck_Analyzer fxCheck(fixedFace);
    if (!fxCheck.IsValid()) {
      printf(" EGADS Info: Face is invalid (fixPCurves)!\n");
      return;
    }
    Face = TopoDS::Face(fixedFace);
    printf("   fixed Face -- toler = %le\n", BRep_Tool::Tolerance(Face));
  }
}


static void coincidentFace(int ie, TopoDS_Shape body, egadsInterEdge &edgx,
                           int index)
{
  int                  i, j;
  double               tol, d, t1, t2, te1, te2;
  Handle(Geom2d_Curve) pcurve;
  Handle(Geom_Curve)   hCurvE;
  
  hCurvE = BRep_Tool::Curve(edgx.newEdge, t1, t2);
  tol    = BRep_Tool::Tolerance(edgx.newEdge);
  
  // look for matching Edge
  TopTools_IndexedMapOfShape emap;
  TopExp::MapShapes(body, TopAbs_EDGE, emap);
  for (i = 1; i <= emap.Extent(); i++) {
    TopoDS_Edge edge = TopoDS::Edge(emap(i));
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, te1, te2);
    if (sameCurve(hCurvE, hCurve)) {
      printf(" %d/%d: sameCurve ts cut = %lf %lf, orig = %lf %lf --",
             ie, index, t1, t2, te1, te2);
      if ((t1 >= te1) && (t2 <= te2)) {
        printf(" OK!\n");
        edgx.Edge[index] = edge;
        return;
      }
      printf(" not used!\n");
    }
  }

  // store away Face with Edge on that Face
  TopTools_IndexedMapOfShape fmap;
  TopExp::MapShapes(body, TopAbs_FACE, fmap);
  for (i = 1; i <= fmap.Extent(); i++) {
    TopoDS_Face Face = TopoDS::Face(fmap(i));
    Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
    pcurve = BRep_Tool::CurveOnSurface(edgx.newEdge, Face, te1, te2);
    if (pcurve.IsNull()) {
      try {
        pcurve = GeomProjLib::Curve2d(hCurvE, t1, t2, hSurface, tol);
        if (pcurve.IsNull()) {
          printf(" EGADS Error: Cannot get PCurve (coincidentFace)!\n");
          continue;
        }
      }
      catch (Standard_Failure)
      {
        printf(" EGADS Info: Geometry Creation Error (coincidentFace)!\n");
        Handle_Standard_Failure e = Standard_Failure::Caught();
        printf("                %s\n", e->GetMessageString());
        continue;
      }
      catch (...)
      {
        printf(" EGADS Info: Geometry Creation Error (coincidentFace)!\n");
        continue;
      }
    }
    for (j = 0; j < 11; j++) {
      double t = t1 + j*(t2-t1)/10.0;
      gp_Pnt pe, ps;
      hCurvE->D0(t, pe);
      gp_Pnt2d p2d;
      pcurve->D0(t, p2d);
      hSurface->D0(p2d.X(), p2d.Y(), ps);
      d = sqrt((pe.X()-ps.X())*(pe.X()-ps.X()) +
               (pe.Y()-ps.Y())*(pe.Y()-ps.Y()) +
               (pe.Z()-ps.Z())*(pe.Z()-ps.Z()));
      if (d > tol) break;
    }
    if (j == 11) {
      printf(" %d/%d: Found Face!\n", ie, index);
      edgx.oFace[index] = Face;
      return;
    }
  }
}


static void fillInterNode(TopoDS_Shape shapes[2], TopoDS_Face Face[2],
                          egadsInterNode &node)
{
  gp_Pnt pv  = BRep_Tool::Pnt(node.newNode);
  double tol = BRep_Tool::Tolerance(node.newNode);
  node.t[0]  = node.t[1] = 0.0;
  
  // look at existing Nodes
  for (int index = 0; index < 2; index++) {
    TopTools_IndexedMapOfShape vmap;
    if (Face[index].IsNull()) {
      TopExp::MapShapes(shapes[index], TopAbs_VERTEX, vmap);
    } else {
      TopExp::MapShapes(Face[index],   TopAbs_VERTEX, vmap);
    }
    for (int i = 0; i < vmap.Extent(); i++) {
      TopoDS_Vertex vert = TopoDS::Vertex(vmap(i+1));
      gp_Pnt pvx   = BRep_Tool::Pnt(vert);
      double toler = BRep_Tool::Tolerance(vert);
      if (toler < tol) toler = tol;
      double d = sqrt((pv.X()-pvx.X())*(pv.X()-pvx.X()) +
                      (pv.Y()-pvx.Y())*(pv.Y()-pvx.Y()) +
                      (pv.Z()-pvx.Z())*(pv.Z()-pvx.Z()));
      if (d <= toler) {
        node.Node[index] = vert;
        printf("   %d", index);
        if (Face[index].IsNull()) {
          printf(" ");
        } else {
          printf("F");
        }
        printf(": at Node %d -- %le (%le)\n", i+1, d, toler);
        break;
      }
    }
  }
  
  // look at cutting existing Edges
  for (int index = 0; index < 2; index++) {
    if (!node.Node[index].IsNull()) continue;
    TopTools_IndexedMapOfShape emap;
    if (Face[index].IsNull()) {
      TopExp::MapShapes(shapes[index], TopAbs_EDGE, emap);
    } else {
      TopExp::MapShapes(Face[index],   TopAbs_EDGE, emap);
    }
    for (int i = 0; i < emap.Extent(); i++) {
      TopoDS_Edge edge = TopoDS::Edge(emap(i+1));
      if (BRep_Tool::Degenerated(edge)) continue;
      Standard_Real te1, te2;
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, te1, te2);
      double toler = BRep_Tool::Tolerance(edge);
      if (toler < tol) toler = tol;
      GeomAPI_ProjectPointOnCurve projPnt(pv, hCurve);
      gp_Pnt pnt;
      double t;
      if (projPnt.NbPoints() == 0) {
        double t, result[3], xyz[3];
        xyz[0] = pv.X();
        xyz[1] = pv.Y();
        xyz[2] = pv.Z();
        EG_nearestCurve(hCurve, xyz, te1, te2, 0, &t, result);
        pnt.SetX(result[0]);
        pnt.SetY(result[1]);
        pnt.SetZ(result[2]);
      } else {
        pnt = projPnt.NearestPoint();
        t   = projPnt.LowerDistanceParameter();
      }
      if (hCurve->IsPeriodic()) {
        double period = hCurve->Period();
        adjPeriodics(period, te1, te2, t);
      }
      if ((t <= te1) || (t >= te2)) {
        printf("      t out of range %lf [%lf-%lf] %d (fillInterNode)!\n",
               t, te1, te2, projPnt.NbPoints());
        continue;
      }
      double d = sqrt((pv.X()-pnt.X())*(pv.X()-pnt.X()) +
                      (pv.Y()-pnt.Y())*(pv.Y()-pnt.Y()) +
                      (pv.Z()-pnt.Z())*(pv.Z()-pnt.Z()));
      if (d <= toler) {
        node.Edge[index] = edge;
        node.t[index]    = t;
        printf("   %d", index);
        if (Face[index].IsNull()) {
          printf(" ");
        } else {
          printf("F");
        }
        printf(": at Edge %d -- %le (%le) t = %lf [%lf-%lf]\n",
               i+1, d, toler, node.t[index], te1, te2);
        break;
      }
    }
  }
}


static void fillInterEdge(TopoDS_Shape body, egadsInterEdge &edgx, int index,
                          egadsInterNode &node)
{
  Standard_Real t1, t2;
  
  TopoDS_Edge Edge          = edgx.newEdge;
  gp_Pnt pv                 = BRep_Tool::Pnt(node.newNode);
  double tol                = BRep_Tool::Tolerance(node.newNode);
  Handle(Geom_Curve) hCurvE = BRep_Tool::Curve(Edge, t1, t2);
  node.t[index]             = 0.0;
  
  // look at existing Nodes
  TopTools_IndexedMapOfShape vmap;
  TopExp::MapShapes(body, TopAbs_VERTEX, vmap);
  for (int i = 0; i < vmap.Extent(); i++) {
    TopoDS_Vertex vert = TopoDS::Vertex(vmap(i+1));
    gp_Pnt pvx   = BRep_Tool::Pnt(vert);
    double toler = BRep_Tool::Tolerance(vert);
    if (toler < tol) toler = tol;
    double d = sqrt((pv.X()-pvx.X())*(pv.X()-pvx.X()) +
                    (pv.Y()-pvx.Y())*(pv.Y()-pvx.Y()) +
                    (pv.Z()-pvx.Z())*(pv.Z()-pvx.Z()));
    if (d <= toler) {
      node.Node[index] = vert;
      printf("   %dX: at Node %d -- %le (%le)\n", index, i+1, d, toler);
      return;
    }
  }
  
  // look at cutting existing Edges
  TopTools_IndexedMapOfShape emap;
  TopExp::MapShapes(body, TopAbs_EDGE, emap);
  for (int i = 0; i < emap.Extent(); i++) {
    TopoDS_Edge edge = TopoDS::Edge(emap(i+1));
    if (BRep_Tool::Degenerated(edge)) continue;
    Standard_Real te1, te2;
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, te1, te2);
    double toler = BRep_Tool::Tolerance(edge);
    if (toler < tol) toler = tol;
    GeomAPI_ProjectPointOnCurve projPnt(pv, hCurve);
    gp_Pnt pnt;
    double t;
    if (projPnt.NbPoints() == 0) {
      double result[3], xyz[3];
      xyz[0] = pv.X();
      xyz[1] = pv.Y();
      xyz[2] = pv.Z();
      EG_nearestCurve(hCurve, xyz, te1, te2, 0, &t, result);
      pnt.SetX(result[0]);
      pnt.SetY(result[1]);
      pnt.SetZ(result[2]);
    } else {
      pnt = projPnt.NearestPoint();
      t   = projPnt.LowerDistanceParameter();
    }
    if (hCurve->IsPeriodic()) {
      double period = hCurve->Period();
      adjPeriodics(period, te1, te2, t);
    }
    if ((t <= te1) || (t >= te2)) {
      printf("      t out of range %lf [%lf-%lf] %d (fillInterEdge)!\n",
             t, te1, te2, projPnt.NbPoints());
      continue;
    }
    double d = sqrt((pv.X()-pnt.X())*(pv.X()-pnt.X()) +
                    (pv.Y()-pnt.Y())*(pv.Y()-pnt.Y()) +
                    (pv.Z()-pnt.Z())*(pv.Z()-pnt.Z()));
    if (d <= toler) {
      node.Edge[index] = edge;
      node.t[index]    = t;
      printf("   %dX: at Edge %d -- %le (%le) t = %lf [%lf-%lf]\n",
             index, i+1, d, toler, node.t[index], te1, te2);
      if ( edgx.Edge[index].IsNull()) break;
      if (!edgx.Edge[index].IsSame(edge))
        printf(" Warning: Not same Edge for Coincidence endPnt!\n");
      break;
    }
  }
}


static int insertNodes(TopoDS_Shape &shape, int index, int nInter,
                       egadsInterEdge *inters)
{
  int i, j, k, inode, iedge, iloop, iface;
  
  // find and map all of the newNodes to intersect
  TopTools_IndexedMapOfShape MapNN;
  for (i = 0; i < nInter; i++) {
    if (!inters[i].start.Edge[index].IsNull())
      MapNN.Add(inters[i].start.newNode);
    if (!inters[i].end.Edge[index].IsNull())
      MapNN.Add(inters[i].end.newNode);
  }
  int n = MapNN.Extent();
  printf("   # new Edge-splitting Nodes = %d\n", n);
  if (n == 0) return EGADS_SUCCESS;
  int    *newEdgeNodes = new int[n];
  double *edgeHit      = new double[n];
  for (i = 0; i < n; i++) newEdgeNodes[i] = 0;

  // major loop matching existing Edges with the intersections
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(shape, TopAbs_EDGE, MapE);
  for (iedge = 1; iedge <= MapE.Extent(); iedge++) {
    TopoDS_Edge Edge = TopoDS::Edge(MapE(iedge));
    if (BRep_Tool::Degenerated(Edge)) continue;
    for (i = 0; i < n; i++) edgeHit[i] = -1.e101;
    for (i = 0; i < nInter; i++) {
      if (inters[i].start.Edge[index].IsNull())      continue;
      if (!Edge.IsSame(inters[i].start.Edge[index])) continue;
      printf("   insert start Node on Edge %d at %lf",
             MapE.FindIndex(inters[i].start.Edge[index]),
             inters[i].start.t[index]);
      j = MapNN.FindIndex(inters[i].start.newNode);
      if (newEdgeNodes[j-1] == 0) printf(" -- first for newNode %d", j);
      newEdgeNodes[j-1]++;
      edgeHit[j-1] = inters[i].start.t[index];
      printf("\n");
    }
    for (i = 0; i < nInter; i++) {
      if (inters[i].end.Edge[index].IsNull())      continue;
      if (!Edge.IsSame(inters[i].end.Edge[index])) continue;
      printf("   insert  end  Node on Edge %d at %lf",
             MapE.FindIndex(inters[i].end.Edge[index]), inters[i].end.t[index]);
      j = MapNN.FindIndex(inters[i].end.newNode);
      if (newEdgeNodes[j-1] == 0) printf(" -- first for newNode %d", j);
      newEdgeNodes[j-1]++;
      edgeHit[j-1] = inters[i].end.t[index];
      printf("\n");
    }
    
    // get all of the newNodes and where they intersect the Edge (in t)
    for (j = i = 0; i < n; i++) if (edgeHit[i] > -1.e100) j++;
    if (j == 0) continue;
    printf("   Edge intersections: ");
    for (i = 0; i < n; i++)
      if (edgeHit[i] > -1.e100) {
        printf(" %lf", edgeHit[i]);
      } else {
        printf(" *");
      }
    printf("\n");

    // make the split Edges
    TopTools_IndexedMapOfShape MapL;
    TopExp::MapShapes(shape, TopAbs_WIRE, MapL);
    double        t1, t2, ts;
    TopoDS_Vertex V1, V2, Vs;
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
    TopExp::Vertices(Edge, V2, V1, Standard_True);
    if (Edge.Orientation() != TopAbs_REVERSED) {
      Vs = V1;
      V1 = V2;
      V2 = Vs;
    }
    TopoDS_Edge *splitEdges = new TopoDS_Edge[j+1];
    Standard_Real old = BRepBuilderAPI::Precision();
    BRepBuilderAPI::Precision(BRep_Tool::Tolerance(Edge));
    for (i = 0; i <= j; i++) {
      if (i == j) {
        Vs = V2;
        ts = t2;
      } else {
        ts    = t2;
        inode = n;
        for (k = 0; k < n; k++) {
          if (edgeHit[k] < -1.e100) continue;
          if ((edgeHit[k] > t1) && (edgeHit[k] < ts)) {
            inode = k;
            ts    = edgeHit[k];
          }
        }
        if (inode == n) {
          printf(" Error: edgeHit Not Found (insertNodes)!\n");
//        return EGADS_NOTFOUND;
        }
        Vs = TopoDS::Vertex(MapNN(inode+1));
      }
      printf("     making Edge from %lf to %lf\n", t1, ts);
      BRepBuilderAPI_MakeEdge MEdge;
      MEdge.Init(hCurve, V1, Vs, t1, ts);
      if (!MEdge.IsDone()) {
        printf(" EGADS Error: Problem with the Edge!\n");
//      return EGADS_NODATA;
      }
      splitEdges[i] = MEdge.Edge();
      BRepCheck_Analyzer eCheck(splitEdges[i]);
      if (!eCheck.IsValid()) {
        printf(" EGADS Info: Edge %d is invalid!\n", i+1);
        EG_checkStatus(eCheck.Result(splitEdges[i]));
//      return EGADS_CONSTERR;
      }
      V1 = Vs;
      t1 = ts;
    }
    BRepBuilderAPI::Precision(old);
    
    // find loops to insert the split Edges
    BRepTools_ReShape reshape;
    egadsInterLoop *loopData = new egadsInterLoop[MapL.Extent()];
    for (iloop = 1; iloop <= MapL.Extent(); iloop++) {
      TopoDS_Wire Loop = TopoDS::Wire(MapL(iloop));
      TopTools_IndexedMapOfShape MapLE;
      TopExp::MapShapes(Loop, TopAbs_EDGE, MapLE);
      k = MapLE.FindIndex(Edge);
      if (k <= 0) continue;
      loopData[iloop-1].oldLoop = Loop;
      loopData[iloop-1].oldFace = getFace(shape, Loop);
#ifdef DEBUG
      int sense[3] = {1, 1, 1};
      if (Loop.Orientation() == TopAbs_REVERSED) sense[2] = -1;
      if (Edge.Orientation() == TopAbs_REVERSED) sense[0] = -1;
      TopoDS_Shape lEdg = MapLE(j);
      if (lEdg.Orientation() == TopAbs_REVERSED) sense[1] = -1;
      printf("    Found in Loop %d: nseg = %d (nEdge = %d) -- %d %d  %d!\n",
             iloop, j+1, MapLE.Extent(), sense[0], sense[1], sense[2]);
#endif
      // reinsert updated loops
      BRepBuilderAPI_MakeWire MW;
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(Loop); ExpWE.More(); ExpWE.Next()) {
        TopoDS_Edge edge = TopoDS::Edge(ExpWE.Current());
        if (Edge.IsSame(edge)) {
          if (edge.Orientation() == TopAbs_REVERSED) {
            for (i = j; i >= 0; i--) {
              TopoDS_Edge edg = splitEdges[i];
              edg.Orientation(TopAbs_REVERSED);
              MW.Add(edg);
            }
          } else {
            for (i = 0; i <= j; i++) {
              TopoDS_Edge edg = splitEdges[i];
              edg.Orientation(TopAbs_FORWARD);
              MW.Add(edg);
            }
          }
        } else {
          MW.Add(edge);
        }
      }
      if (!MW.IsDone()) {
        printf(" EGADS Error: Problem with Loop %d!\n", iloop);
//      return EGADS_NODATA;
      }
      TopoDS_Wire wire = MW.Wire();
      BRepCheck_Analyzer wCheck(wire);
      if (!wCheck.IsValid()) {
        printf(" EGADS Info: Wire is invalid %d!\n", iloop);
        EG_checkStatus(wCheck.Result(wire));
//      return EGADS_CONSTERR;
      }
      reshape.Replace(Loop, wire);
      loopData[iloop-1].newLoop = wire;
    }
    TopoDS_Shape newShape = reshape.Apply(shape, TopAbs_WIRE);
    
    // fix up the PCurves on the affected Faces
    BRepCheck_Analyzer sCheck(newShape);
    TopTools_IndexedMapOfShape MapNF;
    TopExp::MapShapes(newShape, TopAbs_FACE, MapNF);
    for (iface = 1; iface <= MapNF.Extent(); iface++) {
      TopoDS_Face Face = TopoDS::Face(MapNF(iface));
      if (!sCheck.IsValid(Face)) {
#ifdef DEBUG
        printf("   Face %d is Not OK!\n", iface);
#endif
        BRepTools_ReShape reshapeFace;
        TopoDS_Face newFace = Face;
        fixPCurves(newFace);
        reshapeFace.Replace(Face, newFace);
        newShape = reshapeFace.Apply(newShape, TopAbs_FACE);
        // update the loopData info
        for (iloop = 1; iloop <= MapL.Extent(); iloop++)
          if (!loopData[iloop-1].oldLoop.IsNull()) {
            TopoDS_Shape oflu = lookupShapeByIndex(shape, newShape, TopAbs_FACE,
                                                   loopData[iloop-1].oldFace);
            TopoDS_Face faclu = TopoDS::Face(oflu);
            if (faclu.IsSame(newFace)) {
              TopExp_Explorer ExpW;
              ExpW.Init(newFace, TopAbs_WIRE);
              loopData[iloop-1].newLoop = TopoDS::Wire(ExpW.Current());
#ifdef DEBUG
              printf("   Patching up Face %d (%d)!\n", iface, iloop);
#endif
            }
          }
      }
    }
    sCheck.Init(newShape);
    if (!sCheck.IsValid()) {
      printf(" EGADS Info: newShape %d is InValid!\n", index);
//    EG_checkStatus(sCheck.Result(newShape));
//    BRepTools::Write(shape, "Output.BRep");
    } else {
      // update the modified Face objects
      for (iloop = 1; iloop <= MapL.Extent(); iloop++)
        if (!loopData[iloop-1].oldLoop.IsNull()) {
          TopoDS_Face newFace = getFace(newShape, loopData[iloop-1].newLoop);
          for (i = 0; i < nInter; i++)
            if (inters[i].mFace[index].IsSame(loopData[iloop-1].oldFace))
              inters[i].mFace[index] = newFace;
        }
      // update the working shape
      shape = newShape;
    }
    delete [] loopData;
    delete [] splitEdges;
  }
  
  delete [] newEdgeNodes;
  return EGADS_SUCCESS;
}


static int insertEdges(TopoDS_Shape oshape, TopoDS_Shape &shape, int index,
                       int nInter, egadsInterEdge *inters,
                       TopTools_ListOfShape *faceList)
{
  int      i, j, k, prev, next, nCut, nEdge, nFace, nwe, iface, iLoop, oIndex;
  double   t1, t2, ts, umin, umax, vmin, vmax, Uper, Vper, u1, u2, v1, v2, toler;
  gp_Pnt2d p2f, p2l;

  // get our Mappings
  TopTools_IndexedMapOfShape MapF, oMapF, MapS;
  TopExp::MapShapes( shape, TopAbs_FACE,  MapF);
  TopExp::MapShapes(oshape, TopAbs_FACE, oMapF);
  
  // major Loop matching existing Edges with the intersections
  for (iface = 1; iface <= MapF.Extent(); iface++) {
    TopoDS_Face Face = TopoDS::Face(MapF(iface));
    for (nCut = i = 0; i < nInter; i++) {
      if (inters[i].mFace[index].IsNull())      continue;
      if (!Face.IsSame(inters[i].mFace[index])) continue;
      oIndex = oMapF.FindIndex(inters[i].oFace[index]);
      nCut++;
    }
    if (nCut == 0) continue;
    TopExp::MapShapes(shape, TopAbs_SHELL, MapS);
    if (MapS.Extent() == 0) {
      printf(" Error: No Shell!\n");
      return EGADS_NOTBODY;
    }
    TopoDS_Wire          outer    = BRepTools::OuterWire(Face);
    Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
    Handle(Geom_Plane)   hPlane   = Handle(Geom_Plane)::DownCast(hSurface);
    BRepTools::UVBounds(Face, umin, umax, vmin, vmax);
    Uper = Vper = 0.0;
    if (hPlane.IsNull()) {
      Uper = hSurface->UPeriod();
      Vper = hSurface->VPeriod();
    }

    // get all of the Edges
    //     * only include those from inner loops that are cut
    //     * or non-cut loops firure out what subFace they belong in
    //     * leave outer loop alone if only inners are cut
    //
    TopTools_IndexedMapOfShape MapL;
    TopExp::MapShapes(Face, TopAbs_WIRE, MapL);
    nEdge = 2*nCut;
    for (iLoop = 1; iLoop <= MapL.Extent(); iLoop++) {
      TopoDS_Wire Loop = TopoDS::Wire(MapL(iLoop));
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(Loop); ExpWE.More(); ExpWE.Next()) nEdge++;
    }
    TopoDS_Edge *edges = new TopoDS_Edge[nEdge];
    j = 0;
    for (iLoop = 1; iLoop <= MapL.Extent(); iLoop++) {
      TopoDS_Wire Loop = TopoDS::Wire(MapL(iLoop));
      BRepTools_WireExplorer ExpWE;
      for (ExpWE.Init(Loop); ExpWE.More(); ExpWE.Next()) {
        edges[j] = TopoDS::Edge(ExpWE.Current());
        j++;
      }
    }
    for (i = 0; i < nInter; i++) {
      if (inters[i].mFace[index].IsNull())      continue;
      if (!Face.IsSame(inters[i].mFace[index])) continue;
      int           hit = 0, degen = 0;
      double        t1, t2;
      TopoDS_Vertex V1, V2, Vs;
      TopoDS_Edge   Edge = inters[i].newEdge;
      Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
      TopExp::Vertices(Edge, V2, V1, Standard_True);
      if (!inters[i].start.Node[index].IsNull()) {
        V2 = inters[i].start.Node[index];
        t1 = inters[i].start.t[index];
        hit++;
      }
      if (!inters[i].end.Node[index].IsNull()) {
        V1 = inters[i].end.Node[index];
        t2 = inters[i].end.t[index];
        hit++;
      }
      if (hit != 0) {
        printf(" ts = %lf %lf", t1, t2);
        if (inters[i].newEdge.Orientation() != TopAbs_REVERSED) {
          Vs = V1;
          V1 = V2;
          V2 = Vs;
        }
        if (V1.IsSame(V2)) {
          printf("  Same Vertex!\n");
          degen = 1;
        }
        printf("\n");
        BRepBuilderAPI_MakeEdge MEdge;
        MEdge.Init(hCurve, V1, V2, t1, t2);
        if (!MEdge.IsDone()) {
          printf(" EGADS Error: Problem with the Edge @ Intersection %d!\n", i+1);
//        return EGADS_NODATA;
        }
        Edge = MEdge.Edge();
      }
      edges[j  ] = Edge;
      edges[j+1] = Edge;
      edges[j  ].Orientation(TopAbs_FORWARD);
      setPCurve(edges[j  ], Face);
      edges[j+1].Orientation(TopAbs_REVERSED);
      setPCurve(edges[j+1], Face);
/*    if (degen == 1) {
        BRep_Tool::Degenerated(edges[j  ]);
        BRep_Tool::Degenerated(edges[j+1]);
      }  */
      j += 2;
    }
//#ifdef DEBUG
    printf(" %2d: nEdge = %d, nCut = %d\n", iface, nEdge, nCut);
    printf("   Urange = %lf %lf, %lf    Vrange = %lf %lf, %lf\n",
           umin, umax, Uper, vmin, vmax, Vper);
    TopTools_IndexedMapOfShape MapN;
    for (i = 0; i < nEdge; i++) {
      TopoDS_Vertex last  = TopExp::LastVertex(edges[i],  Standard_True);
      TopoDS_Vertex first = TopExp::FirstVertex(edges[i], Standard_True);
      gp_Pnt pvf = BRep_Tool::Pnt(first);
      gp_Pnt pvl = BRep_Tool::Pnt(last);
      MapN.Add(first);
      MapN.Add(last);
      BRep_Tool::Range(edges[i], t1, t2);
      int sense = 1;
      if (edges[i].Orientation() == TopAbs_REVERSED) {
        ts    = t1;
        t1    = t2;
        t2    = ts;
        sense = -1;
      }
      BRepAdaptor_Curve2d Curve2d(edges[i], Face);
      Curve2d.D0(t1, p2f);
      Curve2d.D0(t2, p2l);
      u1 = p2f.X();
      u2 = p2l.X();
      v1 = p2f.Y();
      v2 = p2l.Y();
      printf(" %2d: %d  %lf %lf %lf   %lf %lf\n", i,
             MapN.FindIndex(first), pvf.X(), pvf.Y(), pvf.Z(), u1, v1);
      printf( "     %d  %lf %lf %lf   %lf %lf  %d",
             MapN.FindIndex(last),  pvl.X(), pvl.Y(), pvl.Z(), u2, v2, sense);
      if (BRep_Tool::Degenerated(edges[i])) printf("   Degenerate");
      printf("\n");
/*
      if (BRep_Tool::Degenerated(edges[i])) {
        printf("     t1 = %lf    t2 = %lf  \n", t1, t2);
      } else {
        double te1, te2;
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edges[i], te1, te2);
        hCurve->D0(t1, pvf);
        hCurve->D0(t2, pvl);
        printf("     t1 = %lf    %lf %lf %lf\n", t1, pvf.X(), pvf.Y(), pvf.Z());
        printf("     t2 = %lf    %lf %lf %lf\n", t2, pvl.X(), pvl.Y(), pvl.Z());
      }
*/
    }
//#endif
    
    // build Loops for each subFace
    next = nFace = 0;
    TopTools_ListOfShape subFaces;
    do {
      int degen = 0;
      if (BRep_Tool::Degenerated(edges[next])) degen = 1;
      printf("  Starting with %d  degen = %d!\n", next, degen);
      BRepBuilderAPI_MakeWire MW;
      MW.Add(edges[next]);
      TopoDS_Vertex first = TopExp::FirstVertex(edges[next], Standard_True);
      TopoDS_Vertex workn = TopExp::LastVertex(edges[next],  Standard_True);
      BRep_Tool::Range(edges[next], t1, t2);
      if (edges[next].Orientation() == TopAbs_REVERSED) {
        ts = t1;
        t1 = t2;
        t2 = ts;
      }
      BRepAdaptor_Curve2d Curve2d(edges[next], Face);
      Curve2d.D0(t2, p2l);
      edges[next].Nullify();
      nwe = 1;
      
      // follow the loop around and make new Wire
      while (!workn.IsSame(first) || (degen == 1)) {
        degen = 0;
#ifdef DEBUG
        printf(" %d: looking for Vertex %d\n", next, MapN.FindIndex(workn));
#endif
        prev = next;
        next = -1;
        for (i = nEdge-2*nCut; i < nEdge; i++) {
          if (edges[i].IsNull()) continue;
          if ((nEdge-i-1)/2 == (nEdge-prev-1)/2) continue;
          TopoDS_Vertex v = TopExp::FirstVertex(edges[i], Standard_True);
          if (v.IsSame(workn)) {
            toler = 2.0*BRep_Tool::Tolerance(edges[i]);
            BRep_Tool::Range(edges[i], t1, t2);
            if (edges[i].Orientation() == TopAbs_REVERSED) {
              ts = t1;
              t1 = t2;
              t2 = ts;
            }
            BRepAdaptor_Curve2d Curve2x(edges[i], Face);
            Curve2x.D0(t1, p2f);
            printf("  %d: dist = %le (%le)\n", i, p2l.Distance(p2f), toler);
            if (p2l.Distance(p2f) > toler) continue;
            Curve2x.D0(t2, p2l);
            next = i;
            break;
          }
        }
        if (next == -1)
          for (i = 0; i < nEdge-2*nCut; i++) {
            if (edges[i].IsNull()) continue;
            TopoDS_Vertex v = TopExp::FirstVertex(edges[i], Standard_True);
            if (v.IsSame(workn)) {
              toler = 2.0*BRep_Tool::Tolerance(edges[i]);
              BRep_Tool::Range(edges[i], t1, t2);
              if (edges[i].Orientation() == TopAbs_REVERSED) {
                ts = t1;
                t1 = t2;
                t2 = ts;
              }
              BRepAdaptor_Curve2d Curve2x(edges[i], Face);
              Curve2x.D0(t1, p2f);
              printf("  %d: dist = %le (%le)\n", i, p2l.Distance(p2f), toler);
              if (p2l.Distance(p2f) > toler) continue;
              Curve2x.D0(t2, p2l);
              next = i;
              break;
            }
          }
        if (next == -1) {
          printf(" ERROR: can not find next Edge!\n");
          break;
        }
        MW.Add(edges[next]);
        workn = TopExp::LastVertex(edges[next], Standard_True);
        edges[next].Nullify();
        nwe++;
      }
      nFace++;
      if (!MW.IsDone()) {
        printf(" EGADS Error: Problem with Loop for subFace %d!\n", nFace);
//      return EGADS_NODATA;
      }
      TopoDS_Wire wire = MW.Wire();
      BRepCheck_Analyzer wCheck(wire);
      if (!wCheck.IsValid()) {
        printf(" EGADS Info: Wire is invalid for subFace %d!\n", nFace);
        EG_checkStatus(wCheck.Result(wire));
//      return EGADS_CONSTERR;
      }
      
      toler = BRep_Tool::Tolerance(Face);
      BRepBuilderAPI_MakeFace MFace;
      Standard_Real old = BRepBuilderAPI::Precision();
      BRepBuilderAPI::Precision(toler);
#if CASVER >= 652
      MFace.Init(hSurface, Standard_False, toler);
#else
      MFace.Init(hSurface, Standard_False);
#endif
      MFace.Add(wire);
      TopoDS_Face face = MFace.Face();
      if (Face.Orientation() == TopAbs_REVERSED) {
        face.Orientation(TopAbs_REVERSED);
      } else {
        face.Orientation(TopAbs_FORWARD);
      }
      BRepBuilderAPI::Precision(old);
      if (!MFace.IsDone()) {
        printf(" EGADS Error: Problem with the Face for subFace %d!\n", nFace);
//      return EGADS_NODATA;
      }
//    if (hPlane.IsNull()) fixPCurves(face);
/*    {
        TopTools_IndexedMapOfShape MapNewL;
        TopExp::MapShapes(face, TopAbs_WIRE, MapNewL);
        for (iLoop = 1; iLoop <= MapNewL.Extent(); iLoop++) {
          printf(" Loop %d:\n", iLoop);
          TopoDS_Wire Loop = TopoDS::Wire(MapNewL(iLoop));
          BRepTools_WireExplorer ExpWE;
          for (ExpWE.Init(Loop); ExpWE.More(); ExpWE.Next()) {
            TopoDS_Edge   edge  = TopoDS::Edge(ExpWE.Current());
            TopoDS_Vertex last  = TopExp::LastVertex(edge,  Standard_True);
            TopoDS_Vertex first = TopExp::FirstVertex(edge, Standard_True);
            gp_Pnt pvf = BRep_Tool::Pnt(first);
            gp_Pnt pvl = BRep_Tool::Pnt(last);
            BRep_Tool::Range(edge, t1, t2);
            if (edge.Orientation() == TopAbs_REVERSED) {
              ts = t1;
              t1 = t2;
              t2 = ts;
            }
            BRepAdaptor_Curve2d Curve2d(edge, face);
            Curve2d.D0(t1, p2f);
            Curve2d.D0(t2, p2l);
            u1 = p2f.X();
            u2 = p2l.X();
            v1 = p2f.Y();
            v2 = p2l.Y();
            printf("  %lf %lf %lf   %lf %lf\n",
                   pvf.X(), pvf.Y(), pvf.Z(), u1, v1);
            printf("  %lf %lf %lf   %lf %lf",
                   pvl.X(), pvl.Y(), pvl.Z(), u2, v2);
            if (BRep_Tool::Degenerated(edge)) printf("   Degenerate");
            printf("\n");
          }
        }
      }  */

      int OK = 0;
      BRepCheck_Analyzer fCheck(face);
      if (!fCheck.IsValid()) {
        // try to fix the fault
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(face);
        sfs->Perform();
        TopoDS_Shape fixedFace = sfs->Shape();
        if (fixedFace.IsNull()) {
          printf(" EGADS Info: Invalid Face w/ NULL Fix for subFace %d!\n",
                 nFace);
          EG_checkStatus(fCheck.Result(face));
//        return  EGADS_CONSTERR;
        }
        BRepCheck_Analyzer fxCheck(fixedFace);
        if (!fxCheck.IsValid()) {
          printf(" EGADS Info: Face is invalid for subFace %d!\n", nFace);
          EG_checkStatus(fxCheck.Result(fixedFace));
//        BRepTools::Write(face, "Output.BRep");
          parseLoops(face);
//        return  EGADS_CONSTERR;
          OK = 1;
        }
        face = TopoDS::Face(fixedFace);
      }
      if (OK == 0) {
        subFaces.Append(face);
        printf("       subFace %d: #Edges = %d\n", nFace, nwe);
      }
      
      for (next = 0; next < nEdge; next++)
        if (!edges[next].IsNull()) break;
    } while (next != nEdge);
    
    delete [] edges;
    
    // get Shell and replace the Face with subFaces
    for (i = 1; i <= MapS.Extent(); i++) {
      TopoDS_Shell Shell = TopoDS::Shell(MapS(i));
      TopTools_IndexedMapOfShape MapSF;
      TopExp::MapShapes(Shell, TopAbs_FACE, MapSF);
      k = MapSF.FindIndex(Face);
      if (k <= 0) continue;

      // rebuild the Shell replacing the Face with the subFaces
      BRep_Builder builder3D;
      TopoDS_Shell shell;
      builder3D.MakeShell(shell);
      for (j = 1; j <= MapSF.Extent(); j++) {
        TopoDS_Face face = TopoDS::Face(MapSF(j));
        if (Face.IsSame(face)) {
          TopTools_ListIteratorOfListOfShape it(subFaces);
          for (; it.More(); it.Next())
            builder3D.Add(shell, TopoDS::Face(it.Value()));
        } else {
          builder3D.Add(shell, face);
        }
      }
      BRepLib::SameParameter(shell);
      BRepCheck_Analyzer shCheck(shell);
      if (!shCheck.IsValid()) {
        printf(" EGADS Info: Shell is invalid for subFace %d!\n", nFace);
//       EG_checkStatus(shCheck.Result(shell));
      } else {
        BRepTools_ReShape reshape;
        reshape.Replace(Shell, shell);
        // orientation?
        TopoDS_Shape newShape = reshape.Apply(shape, TopAbs_SHELL);
        BRepCheck_Analyzer sCheck(newShape);
        if (!sCheck.IsValid()) {
          printf(" EGADS Info: newShape %d is invalid for subFace %d!\n",
                 index, nFace);
        } else {
          shape = newShape;
        }
      }
    }
    
    // append subFaces to the faceList
    faceList[oIndex-1] = subFaces;
    printf("   Cut Face %d  %d,  nCut = %d,  nEdge = %d,   nFaces = %d\n",
           iface, oIndex, nCut, nEdge, nFace);
  }
 
  return EGADS_SUCCESS;
}


int
EG_intersect(TopoDS_Shape src, TopoDS_Shape tool, int outLevel, int *nInter,
             egadsInterEdge **inters)
{
  TopoDS_Shape shapes[2];
  shapes[0] = src;
  shapes[1] = tool;
  *nInter   = 0;
  *inters   = NULL;
  
  // perform the intersection
  BRepAlgoAPI_Section Sec(src, tool, Standard_False);
  Sec.Approximation(Standard_True);
  Sec.Build();
  if (!Sec.IsDone()) {
    if (outLevel > 0)
      printf(" EGADS Error: Can't Intersect (EG_intersect)!\n");
    return EGADS_GEOMERR;
  }
  TopoDS_Shape result = Sec.Shape();
 
  // get our intersection Edges & Nodes
  TopTools_IndexedMapOfShape MapN;
  TopExp::MapShapes(result, TopAbs_VERTEX, MapN);
  int nnode = MapN.Extent();
  
  TopTools_IndexedMapOfShape MapE;
  TopExp::MapShapes(result, TopAbs_EDGE, MapE);
  int nedge = MapE.Extent();
  if ((nedge == 0) && (nnode == 0)) {
    if (outLevel > 0)
      printf(" EGADS Error: No Intersection (EG_intersect)!\n");
    return EGADS_CONSTERR;
  }
  
  // Nodes only!
  if (nedge == 0) {
    egadsInterEdge *interx = new egadsInterEdge[nnode];
    for (int i = 0; i < nnode; i++) {
      interx[i].prev          = 0;
      interx[i].next          = 0;
      interx[i].start.newNode = TopoDS::Vertex(MapN(i+1));
      fillInterNode(shapes, interx[i].oFace, interx[i].start);
    }
    *nInter = nnode;
    *inters = interx;
    return EGADS_SUCCESS;
  }
  
  // Edges
  egadsInterEdge *interx = new egadsInterEdge[nedge];
  for (int i = 0; i < nedge; i++) {
    TopoDS_Shape shape = MapE(i+1);
    TopoDS_Edge  Edge  = TopoDS::Edge(shape);
    interx[i].newEdge  = Edge;
    Sec.HasAncestorFaceOn1(Edge, interx[i].oFace[0]);
    Sec.HasAncestorFaceOn2(Edge, interx[i].oFace[1]);
    if (interx[i].oFace[0].IsNull()) coincidentFace(i, src,  interx[i], 0);
    if (interx[i].oFace[1].IsNull()) coincidentFace(i, tool, interx[i], 1);
    interx[i].mFace[0] = interx[i].oFace[0];
    interx[i].mFace[1] = interx[i].oFace[1];
    interx[i].prev     = 0;
    interx[i].next     = 0;
    // note: Edge can be TopAbs_REVERSED!
    TopExp::Vertices(Edge, interx[i].end.newNode, interx[i].start.newNode,
                     Standard_True);
    if (interx[i].start.newNode.IsSame(interx[i].end.newNode))
      interx[i].prev = interx[i].next = i+1;
  }
  
  // fill in the Node intersection information
  for (int i = 0; i < nedge; i++) {
    printf(" Intersection %d:\n", i);
    fillInterNode(shapes, interx[i].oFace, interx[i].start);
    fillInterNode(shapes, interx[i].oFace, interx[i].end);
    if (interx[i].oFace[0].IsNull()) {
      fillInterEdge(src,  interx[i], 0, interx[i].start);
      fillInterEdge(src,  interx[i], 0, interx[i].end);
    }
    if (interx[i].oFace[1].IsNull()) {
      fillInterEdge(tool, interx[i], 1, interx[i].start);
      fillInterEdge(tool, interx[i], 1, interx[i].end);
    }
    // check -- at least 1 hit is required
    if (interx[i].start.Node[0].IsNull() && interx[i].start.Edge[0].IsNull() &&
        interx[i].start.Node[1].IsNull() && interx[i].start.Edge[1].IsNull())
      printf(" Error: No start Intersection Found!\n");
    if (interx[i].end.Node[0].IsNull()   && interx[i].end.Edge[0].IsNull() &&
        interx[i].end.Node[1].IsNull()   && interx[i].end.Edge[1].IsNull())
      printf(" Error: No end Intersection Found!\n");
  }
  printf("\n");
  
  // connect the Edges
  for (int i = 0; i < nedge; i++)
    for (int j = 0; j < nedge; j++) {
      if (i == j) continue;
      if (interx[i].start.newNode.IsSame(interx[j].start.newNode)) {
        interx[i].prev = j+1;
        if ((interx[j].prev != 0) && (interx[j].prev != i+1))
          printf(" Edge %d: previous hit again with %d %d\n",
                 j+1, interx[j].prev, i+1);
        interx[j].prev = i+1;
      }
      if (interx[i].start.newNode.IsSame(interx[j].end.newNode)) {
        interx[i].prev = j+1;
        if ((interx[j].next != 0) && (interx[j].next != i+1))
          printf(" Edge %d: next hit again with %d %d\n",
                 j+1, interx[j].next, i+1);
        interx[j].next = i+1;
      }
      if (interx[i].end.newNode.IsSame(interx[j].end.newNode)) {
        interx[i].next = j+1;
        if ((interx[j].next != 0) && (interx[j].next != i+1))
          printf(" Edge %d: next hit again with %d %d\n",
                 j+1, interx[j].next, i+1);
        interx[j].next = i+1;
      }
      if (interx[i].end.newNode.IsSame(interx[j].start.newNode)) {
        interx[i].next = j+1;
        if ((interx[j].prev != 0) && (interx[j].prev != i+1))
          printf(" Edge %d: previous hit again with %d %d\n",
                 j+1, interx[j].prev, i+1);
        interx[j].prev = i+1;
      }
    }
  
  *nInter = nnode;
  *inters = interx;
  return EGADS_SUCCESS;
}


extern "C" int
EG_testIntersect(const egObject *srco, const egObject *toolo, const char *name)
{
  int            i, j, stat, start, sense, nInter;
  egadsInterEdge *inters;
  TopoDS_Shape   src, tool, mods, modt;
  
  if ((srco == NULL) || (toolo == NULL)) return EGADS_NULLOBJ;
  if  (srco->magicnumber != MAGIC)       return EGADS_NOTOBJ;
  if ((srco->oclass != MODEL) &&
      (srco->oclass != BODY))            return EGADS_NOTBODY;
  if  (srco->blind  == NULL)             return EGADS_NODATA;
  if ((toolo->oclass != BODY) &&
      (toolo->oclass != FACE))           return EGADS_NOTBODY;
  if  (toolo->blind == NULL)             return EGADS_NODATA;
  int outLevel = EG_outLevel(srco);

  if (srco->oclass == BODY) {
    if (srco->mtype == WIREBODY)         return EGADS_TOPOERR;
    egadsBody *pbods = (egadsBody *) srco->blind;
    src  = pbods->shape;
  } else {
    egadsModel *pmodel = (egadsModel *) srco->blind;
    src  = pmodel->shape;
  }
  if (toolo->oclass == BODY) {
    if (srco->mtype == WIREBODY)         return EGADS_TOPOERR;
    egadsBody *pbodt = (egadsBody *) toolo->blind;
    tool = pbodt->shape;
  } else {
    egadsFace *pface = (egadsFace *) toolo->blind;
    tool = pface->face;
  }
  
  stat = EG_intersect(src, tool, outLevel, &nInter, &inters);
  if (stat != EGADS_SUCCESS) return stat;
/*
  for (i = 0; i < nInter; i++)
    printf(" %d:  %d %d\n", i+1, inters[i].prev, inters[i].next);
 */
  
  // split Edges with New Nodes, then split the Faces with the Edges
  mods = src;
  TopTools_IndexedMapOfShape MapFs;
  TopExp::MapShapes(src,  TopAbs_FACE, MapFs);
  TopTools_ListOfShape *sMapping = new TopTools_ListOfShape[MapFs.Extent()];
  printf(" Insert src  Nodes = %d\n", insertNodes(mods, 0, nInter, inters));
  printf(" Insert src  Edges = %d\n",
         insertEdges(src,  mods, 0, nInter, inters, sMapping));

  modt = tool;
  TopTools_IndexedMapOfShape MapFt;
  TopExp::MapShapes(tool, TopAbs_FACE, MapFt);
  TopTools_ListOfShape *tMapping = new TopTools_ListOfShape[MapFt.Extent()];
  printf(" Insert tool Nodes = %d\n", insertNodes(modt, 1, nInter, inters));
  printf(" Insert tool Edges = %d\n",
         insertEdges(tool, modt, 1, nInter, inters, tMapping));
  
  // remove Face mappings
  for (i = 0; i < MapFs.Extent(); i++)
    if (sMapping[i].Extent() != 0)
      printf(" sMapping %3d: %d\n", i+1, sMapping[i].Extent());
  for (i = 0; i < MapFt.Extent(); i++)
    if (tMapping[i].Extent() != 0)
      printf(" tMapping %3d: %d\n", i+1, tMapping[i].Extent());
  delete [] sMapping;
  delete [] tMapping;
  
  TopoDS_Compound compound;
  BRep_Builder    builder3D;
  builder3D.MakeCompound(compound);
  builder3D.Add(compound, mods);
  builder3D.Add(compound, modt);
  
  int nWire = 0;
  // find and make the loops
  for (;;) {
    BRepBuilderAPI_MakeWire MW;
    // find open beginnings
    for (i = 0; i < nInter; i++) {
      if (inters[i].newEdge.IsNull()) continue;
      if (inters[i].prev == 0) break;
    }
    if (i != nInter) {
      MW.Add(inters[i].newEdge);
      inters[i].prev = -(nInter+1);
      start = i+1;
      sense = 1;
      while (((sense ==  1) && (inters[start-1].next > 0)) ||
             ((sense == -1) && (inters[start-1].prev > 0))) {
        if (sense == 1) {
          j = inters[start-1].next;
          inters[start-1].next = -j;
        } else {
          j = inters[start-1].prev;
          inters[start-1].prev = -j;
        }
        if (inters[j-1].prev == start) {
          MW.Add(inters[j-1].newEdge);
          sense =  1;
        } else if (inters[j-1].next == start) {
          MW.Add(TopoDS::Edge(inters[j-1].newEdge.Reversed()));
          sense = -1;
        } else {
          printf(" ERROR: Cannot find orientation!\n");
          break;
        }
        start = j;
      }
      builder3D.Add(compound, MW.Wire());
      nWire++;
      continue;
    }
    // find open endings
    for (i = 0; i < nInter; i++) {
      if (inters[i].newEdge.IsNull()) continue;
      if (inters[i].next == 0) break;
    }
    if (i != nInter) {
      MW.Add(TopoDS::Edge(inters[i].newEdge.Reversed()));
      inters[i].next = -(nInter+1);
      start = i+1;
      sense = -1;
      while (((sense ==  1) && (inters[start-1].next > 0)) ||
             ((sense == -1) && (inters[start-1].prev > 0))) {
        if (sense == 1) {
          j = inters[start-1].next;
          inters[start-1].next = -j;
        } else {
          j = inters[start-1].prev;
          inters[start-1].prev = -j;
        }
        if (inters[j-1].prev == start) {
          MW.Add(inters[j-1].newEdge);
          sense =  1;
        } else if (inters[j-1].next == start) {
          MW.Add(TopoDS::Edge(inters[j-1].newEdge.Reversed()));
          sense = -1;
        } else {
          printf(" ERROR: Cannot find orientation!\n");
          break;
        }
        start = j;
      }
      builder3D.Add(compound, MW.Wire());
      nWire++;
      continue;
    }
    
    // find any other start
    for (i = 0; i < nInter; i++) {
      if (inters[i].newEdge.IsNull()) continue;
      if (inters[i].prev > 0) break;
    }
    if (i == nInter) break;
    MW.Add(inters[i].newEdge);
    if (inters[i].prev == i+1) {
      inters[i].prev = -inters[i].prev;
      inters[i].next = -inters[i].next;
    } else {
      start = i+1;
      sense = 1;
      while (((sense ==  1) && (inters[start-1].next > 0)) ||
             ((sense == -1) && (inters[start-1].prev > 0))) {
        if (sense == 1) {
          j = inters[start-1].next;
          inters[start-1].next = -j;
          inters[start-1].prev = -inters[start-1].prev;
        } else {
          j = inters[start-1].prev;
          inters[start-1].prev = -j;
          inters[start-1].next = -inters[start-1].next;
        }
        if (inters[j-1].prev == start) {
          MW.Add(inters[j-1].newEdge);
          sense =  1;
        } else if (inters[j-1].next == start) {
          MW.Add(TopoDS::Edge(inters[j-1].newEdge.Reversed()));
          sense = -1;
        } else {
          break;
        }
        start = j;
      }
    }
/*
    printf(" done with wire %d -- start = %d, sense = %d\n",
           nWire+1, start, sense);
    for (int k = 0; k < nInter; k++)
      printf(" %d:  %d %d\n", k+1, inters[k].prev, inters[k].next);  */
    builder3D.Add(compound, MW.Wire());
    nWire++;
  }
  
  printf(" nInter = %d, nWire = %d\n", nInter, nWire);

  if (!BRepTools::Write(compound, name)) {
    printf(" EGADS Warning: OpenCASCADE Write!\n");
    delete [] inters;
    return EGADS_WRITERR;
  }
  
  delete [] inters;
  
  return EGADS_SUCCESS;
}
