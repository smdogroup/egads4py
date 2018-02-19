/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Copy-based Topology Functions
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

  extern "C" int  EG_destroyTopology( egObject *topo );

  extern "C" int  EG_copyTopology( const egObject *topo,
                                   /*@null@*/ double *xform, egObject **copy );
  extern "C" int  EG_flipTopology( const egObject *topo, egObject **copy );
  extern "C" int  EG_matchBodyFaces( const egObject *bod1, const egObject *bod2,
                                     double toler, int *nmatch, int **match );

  extern     int  EG_traverseBody( egObject *context, int i, egObject *bobj, 
                                   egObject *topObj, egadsBody *body,
                                   int *nerr );
  extern     int  EG_shellClosure( egadsShell *pshell, int flag );
  extern     int  EG_attriBodyCopy( const egObject *src,
                                    /*@null@*/ double *xform, egObject *dst );
  extern     void EG_fillPCurves( TopoDS_Face face, egObject *surfo, 
                                  egObject *loopo, egObject *topObj );
  extern     void EG_completePCurve( egObject *g, Handle(Geom2d_Curve) &hCurv );
  extern     void EG_completeCurve(  egObject *g, Handle(Geom_Curve)   &hCurv );
  extern     void EG_completeSurf(   egObject *g, Handle(Geom_Surface) &hSurf );  


int
EG_copyAttrTopo(egadsBody *pbody, /*@null@*/ double *xform, gp_Trsf form,
                const egObject *src, egObject *dst, egObject *topObj)
{
  int      index, stat;
  egObject *context;
  
  context     = EG_context(dst);
  dst->topObj = topObj;
  if (dst == topObj) dst->topObj = context;

  if (src->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) dst->blind;
    gp_Pnt pv        = BRep_Tool::Pnt(pnode->node);
    pnode->xyz[0]    = pv.X();
    pnode->xyz[1]    = pv.Y();
    pnode->xyz[2]    = pv.Z();
    dst->oclass      = NODE;

  } else if (src->oclass == EDGE) {
  
    TopoDS_Vertex V1, V2;
    Standard_Real t1, t2;

    egObject *geom   = NULL;
    egObject *pn1    = NULL;
    egObject *pn2    = NULL;
    int      degen   = 0;
    egadsEdge *sedge = (egadsEdge *) src->blind;
    egadsEdge *pedge = (egadsEdge *) dst->blind;
    TopoDS_Edge Edge = pedge->edge;
    dst->oclass      = EDGE;
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
    index = pbody->nodes.map.FindIndex(V1);
    if (index > 0) pn1 = pbody->nodes.objs[index-1];
    if (pn1 == NULL) {
      EG_makeObject(context, &pn1);
      if (pn1 != NULL) {
        egadsNode *pnode = new egadsNode;
        pnode->node      = V1;
        pn1->blind       = pnode;
        egObject *snode  = sedge->nodes[0];
        if (Edge.Orientation() != TopAbs_REVERSED) 
          snode          = sedge->nodes[1];
        EG_copyAttrTopo(pbody, xform, form, snode, pn1, topObj);
        if (index > 0) pbody->nodes.objs[index-1] = pn1;
      }
    }
    if (V1.IsSame(V2)) {
      dst->mtype = ONENODE;
      pn2        = pn1;
    } else {
      dst->mtype = TWONODE;
      index = pbody->nodes.map.FindIndex(V2);
      if (index > 0) pn2 = pbody->nodes.objs[index-1];
      if (pn2 == NULL) {
        EG_makeObject(context, &pn2);
        if (pn2 != NULL) {
          egadsNode *pnode = new egadsNode;
          pnode->node      = V2;
          pn2->blind       = pnode;
          egObject *snode  = sedge->nodes[1];
          if (Edge.Orientation() != TopAbs_REVERSED) 
            snode          = sedge->nodes[0];
          EG_copyAttrTopo(pbody, xform, form, snode, pn2, topObj);
          if (index > 0) pbody->nodes.objs[index-1] = pn2;
        }
      }
    }
    if (Edge.Orientation() != TopAbs_REVERSED) {
      pedge->nodes[0] = pn2;
      pedge->nodes[1] = pn1; 
    } else {
      pedge->nodes[0] = pn1;
      pedge->nodes[1] = pn2;
    }

    pedge->curve  = geom;
    pedge->topFlg = 0;
    if (degen == 1) {
      dst->mtype = DEGENERATE;
    } else {
      EG_referenceObject(geom, dst);
    }
    EG_referenceObject(pn1, dst);
    EG_referenceObject(pn2, dst);

  } else if (src->oclass == LOOP) {
  
    int      *senses = NULL, closed = 0, ne = 0;
    egObject **edgeo = NULL;
    egadsLoop *sloop = (egadsLoop *) src->blind;
    egadsLoop *ploop = (egadsLoop *) dst->blind;
    TopoDS_Wire Wire = ploop->loop;
    dst->oclass      = LOOP;
    if (Wire.Closed()) closed = 1;
    if ((ploop->surface == NULL) && (dst == topObj))
      if (sloop->surface != NULL) {
        // top of the hierarchy -- use surface from source and transform
        egObject *geom = sloop->surface;
        if (geom->blind != NULL) {
          egadsSurface         *psurf = (egadsSurface *) geom->blind;
          Handle(Geom_Surface)  hSurf = psurf->handle;
          Handle(Geom_Geometry) nGeom = hSurf->Transformed(form);
          Handle(Geom_Surface)  nSurf = Handle(Geom_Surface)::DownCast(nGeom);
          if (!nSurf.IsNull()) {
            egObject *ngeom;
            stat = EG_makeObject(context, &ngeom);
            if (stat == EGADS_SUCCESS) {
              ngeom->topObj = topObj;
              EG_completeSurf(ngeom, nSurf);
              ploop->surface = ngeom;
            }
          }
        }
      }
    int hit = 1;
    if (ploop->surface != NULL) {
      EG_referenceObject(ploop->surface, dst);
      hit = 2;
    }
    
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) ne++;
    if (ne > 0) {
      edgeo  = new egObject*[hit*ne];
      senses = new int[ne];
    }
    int k = 0;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shapW = ExpWE.Current();
      TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
      edgeo[k]           = NULL;
      senses[k]          = 1;
      if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
      if (hit == 2) edgeo[k+ne] = NULL;
      index = pbody->edges.map.FindIndex(Edge);
      if (index > 0) edgeo[k] = pbody->edges.objs[index-1];
      if (edgeo[k] == NULL) {
        EG_makeObject(context, &edgeo[k]);
        if (edgeo[k] != NULL) {
          egadsEdge *pedge = new egadsEdge;
          pedge->edge      = Edge;
          edgeo[k]->blind  = pedge;
          EG_copyAttrTopo(pbody, xform, form, sloop->edges[k], edgeo[k], topObj);
          if (index > 0) pbody->edges.objs[index-1] = edgeo[k];
        }
      }
      if (edgeo[k] != NULL) EG_referenceObject(edgeo[k], dst);
      k++;
    }
    // deal with Loop at top w/ surface -- copy PCurves from src
    if ((ploop->surface != NULL) && (dst == topObj))
      for (k = 0; k < ne; k++) {
        EG_makeObject(context, &edgeo[k+ne]);
        if (edgeo[k+ne] != NULL) {
          edgeo[k+ne]->topObj = topObj;
          egObject *geom      = sloop->edges[k+ne];
          egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
          Handle(Geom2d_Curve) hCurve = ppcurv->handle;
          EG_completePCurve(edgeo[k+ne],  hCurve);
          EG_referenceObject(edgeo[k+ne], dst);
        }
      }
    ploop->nedges  = ne;
    ploop->edges   = edgeo;
    ploop->senses  = senses;
    ploop->topFlg  = 0;
    dst->mtype     = OPEN;
    if (closed == 1) dst->mtype = CLOSED;
  
  } else if (src->oclass == FACE) {
  
    int      *senses = NULL;
    egObject **loopo = NULL;
    egObject *geom   = NULL;
    egadsFace *sface = (egadsFace *) src->blind;
    egadsFace *pface = (egadsFace *) dst->blind;
    TopoDS_Face Face = pface->face;
    dst->oclass      = FACE;
    stat = EG_makeObject(context, &geom);
    if (stat == EGADS_SUCCESS) {
      geom->topObj = topObj;
      Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
      EG_completeSurf(geom, hSurface);
      EG_referenceObject(geom, dst);
    }
    
    int nl = 0;
    TopExp_Explorer ExpW;
    for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) nl++;
    TopoDS_Wire oWire = BRepTools::OuterWire(Face);
    if (nl > 0) {
      loopo  = new egObject*[nl];
      senses = new int[nl];
    }
    int k = 0;
    for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  Wire  = TopoDS::Wire(shapw);
      loopo[k]           = NULL;
      senses[k]          = -1;
      if (Wire.IsSame(oWire)) senses[k] = 1;
      index = pbody->loops.map.FindIndex(Wire);
      if (index > 0) loopo[k] = pbody->loops.objs[index-1];
      if (loopo[k] == NULL) {
        EG_makeObject(context, &loopo[k]);
        if (loopo[k] != NULL) {
          egadsLoop *ploop = new egadsLoop;
          ploop->loop      = Wire;
          ploop->surface   = geom;
          if (geom != NULL)
            if (geom->mtype == PLANE) ploop->surface = NULL;
          loopo[k]->blind  = ploop;
          EG_copyAttrTopo(pbody, xform, form, sface->loops[k], loopo[k], topObj);
          if (index > 0) pbody->loops.objs[index-1] = loopo[k];
          EG_fillPCurves(Face, geom, loopo[k], topObj);
        }
      }
      if (loopo[k] != NULL) EG_referenceObject(loopo[k], dst);
      k++;
    }
    pface->surface = geom;
    pface->nloops  = nl;
    pface->loops   = loopo;
    pface->senses  = senses;
    pface->topFlg  = 0;
    dst->mtype     = SFORWARD;
    if (Face.Orientation() == TopAbs_REVERSED) dst->mtype = SREVERSE;

  } else {
  
    egObject   **faceo = NULL;
    egadsShell *sshell = (egadsShell *) src->blind;
    egadsShell *pshell = (egadsShell *) dst->blind;
    dst->oclass        = SHELL;
    TopoDS_Shell Shell = pshell->shell;
    int             nf = 0;
    TopExp_Explorer ExpF;
    for (ExpF.Init(Shell, TopAbs_FACE); ExpF.More(); ExpF.Next()) nf++;
      
    if (nf > 0) faceo = new egObject*[nf];

    int k = 0;
    for (ExpF.Init(Shell, TopAbs_FACE); ExpF.More(); ExpF.Next()) {
      if (faceo == NULL) continue;
      TopoDS_Shape shapf = ExpF.Current();
      TopoDS_Face  Face  = TopoDS::Face(shapf);
      faceo[k]           = NULL;
      index = pbody->faces.map.FindIndex(Face);
      if (index > 0) faceo[k] = pbody->faces.objs[index-1];
      if (faceo[k] == NULL) {
        EG_makeObject(context, &faceo[k]);
        if (faceo[k] != NULL) {
          egadsFace *pface = new egadsFace;
          pface->face      = Face;
          faceo[k]->blind  = pface;
          EG_copyAttrTopo(pbody, xform, form, sshell->faces[k], faceo[k], topObj);
          if (index > 0) pbody->faces.objs[index-1] = faceo[k];
        }
      }
      if (faceo[k] != NULL) EG_referenceObject(faceo[k], dst);
      k++;
    }
    pshell->nfaces = nf;
    pshell->faces  = faceo;
    pshell->topFlg = 0;
    dst->mtype     = EG_shellClosure(pshell, 0);
  
  }

  EG_attributeXDup(src, xform, dst);

  return EGADS_SUCCESS;
}


int
EG_copyTopology(const egObject *topo, /*@null@*/ double *xform, egObject **copy)
{
  int             i, stat, nent, outLevel, nerr;
  egObject        *obj, *context;
  egadsBody       ebody;
  TopoDS_Shape    shape, nTopo;
  TopExp_Explorer Exp;

  if  (topo == NULL)               return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((topo->oclass < NODE) || (topo->oclass > MODEL))
                                   return EGADS_NOTTOPO;
  if  (topo->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(topo);
  context  = EG_context(topo);
          
  gp_Trsf form = gp_Trsf();
  if (xform != NULL)
    form.SetValues(xform[ 0], xform[ 1], xform[ 2], xform[ 3],
                   xform[ 4], xform[ 5], xform[ 6], xform[ 7],
                   xform[ 8], xform[ 9], xform[10], xform[11]
#if CASVER < 680
                   , Precision::Confusion(), Precision::Angular()
#endif
                   );
                   
  if (topo->oclass == NODE) {
  
    egadsNode *pnode = (egadsNode *) topo->blind;
    shape = pnode->node;

  } else if (topo->oclass == EDGE) {
  
    egadsEdge *pedge = (egadsEdge *) topo->blind;
    shape = pedge->edge;

  } else if (topo->oclass == LOOP) {
  
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    shape = ploop->loop;
  
  } else if (topo->oclass == FACE) {
  
    egadsFace *pface = (egadsFace *) topo->blind;
    shape = pface->face;
    
  } else if (topo->oclass == SHELL) {
  
    egadsShell *pshell = (egadsShell *) topo->blind;
    shape = pshell->shell;
    
  } else if (topo->oclass == BODY) {
  
    egadsBody *pbody = (egadsBody *) topo->blind;
    shape = pbody->shape;
    
  } else {
  
    egadsModel *pmodel = (egadsModel *) topo->blind;
    shape = pmodel->shape;

  }

  // got the OCC topology -- now transform
  BRepBuilderAPI_Transform xForm(shape, form, Standard_True);
  if (!xForm.IsDone()) {
    printf(" EGADS Error: Can't copy Topology (EG_copyTopology)!\n");
    return EGADS_CONSTERR;
  }
  nTopo = xForm.ModifiedShape(shape);

  // got the new shape -- parse and fill
  stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_copyTopology)!\n");
    return stat;
  }
  if (topo->oclass == NODE) {
 
    egadsNode *pnode   = new egadsNode;
    TopoDS_Vertex Vert = TopoDS::Vertex(nTopo);
    pnode->node        = Vert;
    obj->blind         = pnode;
    EG_copyAttrTopo(NULL, xform, form, topo, obj, obj);

  } else if (topo->oclass == EDGE) {
  
    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;

    egadsEdge *pedge = new egadsEdge;
    TopoDS_Edge Edge = TopoDS::Edge(nTopo);
    pedge->edge      = Edge;
    obj->blind       = pedge;
    EG_copyAttrTopo(&ebody, xform, form, topo, obj, obj);
    
    delete [] ebody.nodes.objs;

  } else if (topo->oclass == LOOP) {
  
    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_EDGE,   ebody.edges.map);
    nent = ebody.edges.map.Extent();
    ebody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;

    egadsLoop *ploop = new egadsLoop;
    TopoDS_Wire Loop = TopoDS::Wire(nTopo);
    ploop->loop      = Loop;
    ploop->surface   = NULL;
    obj->blind       = ploop;
    EG_copyAttrTopo(&ebody, xform, form, topo, obj, obj);
    
    delete [] ebody.edges.objs;
    delete [] ebody.nodes.objs;

  } else if (topo->oclass == FACE) {
  
    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_EDGE,   ebody.edges.map);
    nent = ebody.edges.map.Extent();
    ebody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_WIRE,   ebody.loops.map);
    nent = ebody.loops.map.Extent();
    ebody.loops.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.loops.objs[i] = NULL;

    egadsFace *pface = new egadsFace;
    TopoDS_Face Face = TopoDS::Face(nTopo);
    pface->face      = Face;
    obj->blind       = pface;
    EG_copyAttrTopo(&ebody, xform, form, topo, obj, obj);
    
    delete [] ebody.loops.objs;
    delete [] ebody.edges.objs;
    delete [] ebody.nodes.objs;
    
  } else if (topo->oclass == SHELL) {
  
    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_EDGE,   ebody.edges.map);
    nent = ebody.edges.map.Extent();
    ebody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_WIRE,   ebody.loops.map);
    nent = ebody.loops.map.Extent();
    ebody.loops.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.loops.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_FACE,   ebody.faces.map);
    nent = ebody.faces.map.Extent();
    ebody.faces.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.faces.objs[i] = NULL;

    egadsShell *pshell = new egadsShell;
    TopoDS_Shell Shell = TopoDS::Shell(nTopo);
    pshell->shell      = Shell;
    obj->blind         = pshell;
    EG_copyAttrTopo(&ebody, xform, form, topo, obj, obj);
    
    delete [] ebody.faces.objs;
    delete [] ebody.loops.objs;
    delete [] ebody.edges.objs;
    delete [] ebody.nodes.objs;
    
  } else if (topo->oclass == BODY) {
  
    egadsBody *pbody   = new egadsBody;
    obj->oclass        = BODY;
    obj->mtype         = topo->mtype;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pbody->shape       = nTopo;
    obj->blind         = pbody;
    stat = EG_traverseBody(context, 0, obj, obj, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      delete pbody;
      return stat;
    }
    EG_attriBodyCopy(topo, xform, obj);

  } else {

    egadsModel *pmodel = (egadsModel *) topo->blind;
    egadsModel *mshape = new egadsModel;
    mshape->shape  = nTopo;
    mshape->nbody  = pmodel->nbody;
    mshape->bodies = new egObject*[pmodel->nbody];
    for (i = 0; i < pmodel->nbody; i++) {
      stat = EG_makeObject(context, &mshape->bodies[i]);
      if (stat != EGADS_SUCCESS) {
        for (int j = 0; j < i; j++) {
          egObject  *bobj  = mshape->bodies[j];
          egadsBody *pbody = (egadsBody *) bobj->blind;
          delete pbody;
          EG_deleteObject(mshape->bodies[j]);
        }
        delete [] mshape->bodies;
        delete mshape;
        EG_deleteObject(obj);
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
    TopExp_Explorer Exp;
    for (Exp.Init(mshape->shape, TopAbs_WIRE,  TopAbs_FACE); 
         Exp.More(); Exp.Next()) {
      egObject  *pobj  = mshape->bodies[i++];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      pbody->shape     = Exp.Current();
    }
    for (Exp.Init(mshape->shape, TopAbs_FACE,  TopAbs_SHELL);
         Exp.More(); Exp.Next()) {
      egObject  *pobj  = mshape->bodies[i++];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      pbody->shape     = Exp.Current();
    }
    for (Exp.Init(mshape->shape, TopAbs_SHELL, TopAbs_SOLID);
         Exp.More(); Exp.Next()) {
      egObject  *pobj  = mshape->bodies[i++];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      pbody->shape     = Exp.Current();
    }
    for (Exp.Init(mshape->shape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
      egObject  *pobj  = mshape->bodies[i++];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      pbody->shape     = Exp.Current();
    }

    obj->oclass = MODEL;
    obj->blind  = mshape;
    for (i = 0; i < pmodel->nbody; i++) {
      egObject  *pobj  = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      pobj->topObj     = obj;
      stat = EG_traverseBody(context, i, pobj, obj, pbody, &nerr);
      if (stat != EGADS_SUCCESS) {
        mshape->nbody = i;
        EG_destroyTopology(obj);
        delete [] mshape->bodies;
        delete mshape;
        return stat;
      }
      egObject *sobj = pmodel->bodies[i];
      EG_attriBodyCopy(sobj, xform, pobj);
    }
    EG_attributeXDup(topo, xform, obj);

  }  

  EG_referenceObject(obj, context);
  *copy = obj;
  return EGADS_SUCCESS;
}


static void
EG_fillObjTopo(egadsBody *pbody, const egObject *obj)
{
  int             i, index;
  egObject        *src;
  TopExp_Explorer Exp;
  
  src = (egObject *) obj;
  
  if (src->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) src->blind;
    index = pbody->nodes.map.FindIndex(pnode->node);
    if (index > 0)
      if (pbody->nodes.objs[index-1] == NULL) 
        pbody->nodes.objs[index-1] = src;

  } else if (src->oclass == EDGE) {

    egadsEdge *pedge = (egadsEdge *) src->blind;
    TopoDS_Edge Edge = pedge->edge;
    index = pbody->edges.map.FindIndex(Edge);
    if (index <= 0) {
      EG_fillObjTopo(pbody, pedge->nodes[0]);
      if (src->mtype == TWONODE) EG_fillObjTopo(pbody, pedge->nodes[1]);
    } else {
      if (pbody->edges.objs[index-1] == NULL) {
        EG_fillObjTopo(pbody, pedge->nodes[0]);
        if (src->mtype == TWONODE) EG_fillObjTopo(pbody, pedge->nodes[1]);
        pbody->edges.objs[index-1] = src;
      }
    }

  } else if (src->oclass == LOOP) {
  
    egadsLoop *ploop = (egadsLoop *) src->blind;
    TopoDS_Wire Wire = ploop->loop;
    index = pbody->loops.map.FindIndex(Wire);
    if (index <= 0) {
      for (i = 0; i < ploop->nedges; i++)
        EG_fillObjTopo(pbody, ploop->edges[i]);
    } else {
      if (pbody->loops.objs[index-1] == NULL) {
        for (i = 0; i < ploop->nedges; i++)
          EG_fillObjTopo(pbody, ploop->edges[i]);
        pbody->loops.objs[index-1] = src;
      }
    }
      
  } else if (src->oclass == FACE) {
  
    egadsFace *pface = (egadsFace *) src->blind;
    TopoDS_Face Face = pface->face;
    index = pbody->faces.map.FindIndex(Face);
    if (index <= 0) {
      for (i = 0; i < pface->nloops; i++)
        EG_fillObjTopo(pbody, pface->loops[i]);
    } else {
      if (pbody->faces.objs[index-1] == NULL) {
        for (i = 0; i < pface->nloops; i++)
          EG_fillObjTopo(pbody, pface->loops[i]);
        pbody->faces.objs[index-1] = src;
      }
    }

  } else {
  
    egadsShell *pshell = (egadsShell *) src->blind;
    for (i = 0; i < pshell->nfaces; i++)
      EG_fillObjTopo(pbody, pshell->faces[i]);
  
  }

}


static int
EG_flipAttrTopo(egadsBody *pbody, egadsBody *tbody, const egObject *src, 
                egObject *dst, egObject *topObj)
{
  int      index, stat;
  egObject *context, *obj;
  
  context     = EG_context(dst);
  dst->topObj = topObj;
  if (dst == topObj) dst->topObj = context;

  if (dst->oclass == NODE) {

    egadsNode *pnode = (egadsNode *) dst->blind;
    gp_Pnt pv        = BRep_Tool::Pnt(pnode->node);
    pnode->xyz[0]    = pv.X();
    pnode->xyz[1]    = pv.Y();
    pnode->xyz[2]    = pv.Z();
    index = tbody->nodes.map.FindIndex(pnode->node);
    if (index > 0) EG_attributeDup(tbody->nodes.objs[index-1], dst);

  } else if (dst->oclass == EDGE) {
  
    TopoDS_Vertex V1, V2;
    Standard_Real t1, t2;

    egObject *geom   = NULL;
    egObject *pn1    = NULL;
    egObject *pn2    = NULL;
    int      degen   = 0;
    egadsEdge *pedge = (egadsEdge *) dst->blind;
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
    index = pbody->nodes.map.FindIndex(V1);
    if (index > 0) pn1 = pbody->nodes.objs[index-1];
    if (pn1 == NULL) {
      EG_makeObject(context, &pn1);
      if (pn1 != NULL) {
        egadsNode *pnode = new egadsNode;
        pnode->node      = V1;
        pn1->blind       = pnode;
        pn1->oclass      = NODE;
        EG_flipAttrTopo(pbody, tbody, src, pn1, topObj);
        if (index > 0) pbody->nodes.objs[index-1] = pn1;
      }
    }
    if (V1.IsSame(V2)) {
      dst->mtype = ONENODE;
      pn2        = pn1;
    } else {
      dst->mtype = TWONODE;
      index = pbody->nodes.map.FindIndex(V2);
      if (index > 0) pn2 = pbody->nodes.objs[index-1];
      if (pn2 == NULL) {
        EG_makeObject(context, &pn2);
        if (pn2 != NULL) {
          egadsNode *pnode = new egadsNode;
          pnode->node      = V2;
          pn2->blind       = pnode;
          pn2->oclass      = NODE;
          EG_flipAttrTopo(pbody, tbody, src, pn2, topObj);
          if (index > 0) pbody->nodes.objs[index-1] = pn2;
        }
      }
    }
    if (Edge.Orientation() != TopAbs_REVERSED) {
      pedge->nodes[0] = pn2;
      pedge->nodes[1] = pn1; 
    } else {
      pedge->nodes[0] = pn1;
      pedge->nodes[1] = pn2;
    }

    pedge->curve  = geom;
    pedge->topFlg = 0;
    if (degen == 1) {
      dst->mtype = DEGENERATE;
    } else {
      EG_referenceObject(geom, dst);
    }
    EG_referenceObject(pn1, dst);
    EG_referenceObject(pn2, dst);
    index = tbody->edges.map.FindIndex(Edge);
    if (index > 0) EG_attributeDup(tbody->edges.objs[index-1], dst);

  } else if (dst->oclass == LOOP) {
  
    int      *senses = NULL, closed = 0, ne = 0;
    egObject **edgeo = NULL;
    egadsLoop *ploop = (egadsLoop *) dst->blind;
    TopoDS_Wire Wire = ploop->loop;
    if (Wire.Closed()) closed = 1;
    if ((ploop->surface == NULL) && (dst == topObj)) {
      egadsLoop *sloop = (egadsLoop *) src->blind;
      if (sloop->surface != NULL) {
        // top of the hierarchy -- use surface from source
        egObject *geom = sloop->surface;
        if (geom->blind != NULL) {
          egadsSurface *psurf = (egadsSurface *) geom->blind;
          egObject     *ngeom = NULL;
          EG_makeObject(context, &ngeom);
          if (ngeom != NULL) {
            ngeom->topObj = topObj;
            EG_completeSurf(ngeom, psurf->handle);
            ploop->surface = ngeom;
          }
        } 
      }
    }
    int hit = 1;
    if (ploop->surface != NULL) {
      EG_referenceObject(ploop->surface, dst);
      hit = 2;
    }
    
    BRepTools_WireExplorer ExpWE;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) ne++;
    if (ne > 0) {
      edgeo  = new egObject*[hit*ne];
      senses = new int[ne];
    }
    int k = 0;
    for (ExpWE.Init(Wire); ExpWE.More(); ExpWE.Next()) {
      TopoDS_Shape shapW = ExpWE.Current();
      TopoDS_Edge  Edge  = TopoDS::Edge(shapW);
      edgeo[k]           = NULL;
      senses[k]          = 1;
      if (shapW.Orientation() == TopAbs_REVERSED) senses[k] = -1;
      if (hit == 2) edgeo[k+ne] = NULL;
      index = pbody->edges.map.FindIndex(Edge);
      if (index > 0) edgeo[k] = pbody->edges.objs[index-1];
      if (edgeo[k] == NULL) {
        EG_makeObject(context, &edgeo[k]);
        if (edgeo[k] != NULL) {
          egadsEdge *pedge = new egadsEdge;
          pedge->edge      = Edge;
          edgeo[k]->blind  = pedge;
          edgeo[k]->oclass = EDGE;
          EG_flipAttrTopo(pbody, tbody, src, edgeo[k], topObj);
          if (index > 0) pbody->edges.objs[index-1] = edgeo[k];
        }
      }
      if (edgeo[k] != NULL) EG_referenceObject(edgeo[k], dst);
      k++;
    }
    // deal with Loop at top w/ surface -- copy PCurves from src
    if ((ploop->surface != NULL) && (dst == topObj)) {
      egadsLoop *sloop = (egadsLoop *) src->blind;
      for (k = 0; k < ne; k++) {
        if (edgeo[k] == NULL) continue;
        egadsEdge *pedge = (egadsEdge *) edgeo[k]->blind;
        TopoDS_Edge Edge = TopoDS::Edge(pedge->edge);
        index = -1;
        for (int j = 0; j < ne; j++) {
          obj   = sloop->edges[j];
          pedge = (egadsEdge *) obj->blind;
          // reversed -- pick up correct periodic Edge
          if (( Edge.IsSame( pedge->edge)) && 
              (!Edge.IsEqual(pedge->edge))) {
            index = j;
            break;
          }
        }
        if (index != -1) EG_makeObject(context, &edgeo[k+ne]);
        if (edgeo[k+ne] != NULL) {
          edgeo[k+ne]->topObj = topObj;
          egObject *geom      = sloop->edges[index+ne];
          egadsPCurve *ppcurv = (egadsPCurve *) geom->blind;
          Handle(Geom2d_Curve) hCurve = ppcurv->handle;
          EG_completePCurve(edgeo[k+ne],  hCurve);
          EG_referenceObject(edgeo[k+ne], dst);
        }
      }
    }
    
    ploop->nedges = ne;
    ploop->edges  = edgeo;
    ploop->senses = senses;
    ploop->topFlg = 0;
    dst->mtype    = OPEN;
    if (closed == 1) dst->mtype = CLOSED;
    if (dst != topObj) {
      index = tbody->loops.map.FindIndex(Wire);
      if (index > 0) EG_attributeDup(tbody->loops.objs[index-1], dst);
    }
  
  } else if (dst->oclass == FACE) {
  
    int      *senses = NULL;
    egObject **loopo = NULL;
    egObject *geom   = NULL;
    egadsFace *pface = (egadsFace *) dst->blind;
    TopoDS_Face Face = pface->face;
    stat = EG_makeObject(context, &geom);
    if (stat == EGADS_SUCCESS) {
      geom->topObj = topObj;
      Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
      EG_completeSurf(geom, hSurface);
      EG_referenceObject(geom, dst);
    }
    
    int nl = 0;
    TopExp_Explorer ExpW;
    for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) nl++;
    TopoDS_Wire oWire = BRepTools::OuterWire(Face);
    if (nl > 0) {
      loopo  = new egObject*[nl];
      senses = new int[nl];
    }
    int k = 0;
    for (ExpW.Init(Face, TopAbs_WIRE); ExpW.More(); ExpW.Next()) {
      TopoDS_Shape shapw = ExpW.Current();
      TopoDS_Wire  Wire  = TopoDS::Wire(shapw);
      loopo[k]           = NULL;
      senses[k]          = -1;
      if (Wire.IsSame(oWire)) senses[k] = 1;
      index = pbody->loops.map.FindIndex(Wire);
      if (index > 0) loopo[k] = pbody->loops.objs[index-1];
      if (loopo[k] == NULL) {
        EG_makeObject(context, &loopo[k]);
        if (loopo[k] != NULL) {
          egadsLoop *ploop = new egadsLoop;
          ploop->loop      = Wire;
          ploop->surface   = geom;
          if (geom != NULL)
            if (geom->mtype == PLANE) ploop->surface = NULL;
          loopo[k]->blind  = ploop;
          loopo[k]->oclass = LOOP;
          EG_flipAttrTopo(pbody, tbody, src, loopo[k], topObj);
          if (index > 0) pbody->loops.objs[index-1] = loopo[k];
          EG_fillPCurves(Face, geom, loopo[k], topObj);
        }
      }
      if (loopo[k] != NULL) EG_referenceObject(loopo[k], dst);
      k++;
    }
    pface->surface = geom;
    pface->nloops  = nl;
    pface->loops   = loopo;
    pface->senses  = senses;
    pface->topFlg  = 0;
    dst->mtype     = SFORWARD;
    if (Face.Orientation() == TopAbs_REVERSED) dst->mtype = SREVERSE;
    if (dst != topObj) {
      index = tbody->faces.map.FindIndex(Face);
      if (index > 0) EG_attributeDup(tbody->faces.objs[index-1], dst);
    }

  } else {
  
    egObject   **faceo = NULL;
    egadsShell *pshell = (egadsShell *) dst->blind;
    dst->oclass        = SHELL;
    TopoDS_Shell Shell = pshell->shell;
    int             nf = 0;
    TopExp_Explorer ExpF;
    for (ExpF.Init(Shell, TopAbs_FACE); ExpF.More(); ExpF.Next()) nf++;
      
    if (nf > 0) faceo = new egObject*[nf];

    int k = 0;
    for (ExpF.Init(Shell, TopAbs_FACE); ExpF.More(); ExpF.Next()) {
      if (faceo == NULL) continue;
      TopoDS_Shape shapf = ExpF.Current();
      TopoDS_Face  Face  = TopoDS::Face(shapf);
      faceo[k]           = NULL;
      index = pbody->faces.map.FindIndex(Face);
      if (index > 0) faceo[k] = pbody->faces.objs[index-1];
      if (faceo[k] == NULL) {
        EG_makeObject(context, &faceo[k]);
        if (faceo[k] != NULL) {
          egadsFace *pface = new egadsFace;
          pface->face      = Face;
          faceo[k]->blind  = pface;
          faceo[k]->oclass = FACE;
          EG_flipAttrTopo(pbody, tbody, src, faceo[k], topObj);
          if (index > 0) pbody->faces.objs[index-1] = faceo[k];
        }
      }
      if (faceo[k] != NULL) EG_referenceObject(faceo[k], dst);
      k++;
    }
    pshell->nfaces = nf;
    pshell->faces  = faceo;
    pshell->topFlg = 0;
    dst->mtype     = EG_shellClosure(pshell, 0);
  
  }

  return EGADS_SUCCESS;
}


int
EG_flipTopology(const egObject *topo, egObject **copy)
{
  int             i, nent, outLevel;
  egObject        *context, *obj, *geom;
  egadsBody       ebody, tbody;
  TopoDS_Shape    shape, nTopo;
  egadsEdge       *pedges;
  TopExp_Explorer Exp;
  
  if  (topo == NULL)               return EGADS_NULLOBJ;
  if  (topo->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((topo->oclass < EDGE) || (topo->oclass > SHELL))
                                   return EGADS_NOTTOPO;
  if  (topo->blind == NULL)        return EGADS_NODATA;
  outLevel = EG_outLevel(topo);
  context  = EG_context(topo);

  if (topo->oclass == EDGE) {

    if (topo->mtype == DEGENERATE) return EGADS_DEGEN;
    pedges = (egadsEdge *) topo->blind;
    shape  = pedges->edge;
    
  } else if (topo->oclass == LOOP) {
  
    egadsLoop *ploop = (egadsLoop *) topo->blind;
    shape = ploop->loop;
  
  } else if (topo->oclass == FACE) {
  
    egadsFace *pface = (egadsFace *) topo->blind;
    shape = pface->face;
    
  } else {
  
    egadsShell *pshell = (egadsShell *) topo->blind;
    shape = pshell->shell;

  }

  // got the OCC topology -- now copy and flip
  nTopo = shape.Reversed();
  if (topo->oclass == EDGE) {
    double t1, t2;
    int stat = EG_makeObject(context, &geom);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Cannot make Geom object (EG_flipTopology)!\n");
      return stat;
    }
    TopoDS_Edge Edge = TopoDS::Edge(nTopo);
    Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
    EG_completeCurve(geom, hCurve);
  }
  // got the new shape -- parse and fill
  int stat = EG_makeObject(context, &obj);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Cannot make Body object (EG_flipTopology)!\n");
    return stat;
  }
  
  // allocate our Maps for attribute retrieval and minimal copy
  
  if (topo->oclass == EDGE) {

    egadsEdge *pedge = new egadsEdge;
    TopoDS_Edge Edge = TopoDS::Edge(nTopo);
    pedge->edge      = Edge;
    pedge->curve     = geom;
    pedge->nodes[0]  = pedges->nodes[1];
    pedge->nodes[1]  = pedges->nodes[0];
    pedge->topFlg    = 0;
    obj->oclass      = EDGE;
    obj->blind       = pedge;
    obj->mtype       = topo->mtype;
    EG_referenceObject(geom,            obj);
    EG_referenceObject(pedge->nodes[0], obj);
    EG_referenceObject(pedge->nodes[1], obj);
    EG_attributeDup(topo, obj);

  } else if (topo->oclass == LOOP) {
  
    TopExp::MapShapes(shape, TopAbs_VERTEX, tbody.nodes.map);
    nent = tbody.nodes.map.Extent();
    tbody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.nodes.objs[i] = NULL;
    TopExp::MapShapes(shape, TopAbs_EDGE,   tbody.edges.map);
    nent = tbody.edges.map.Extent();
    tbody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.edges.objs[i] = NULL;
    EG_fillObjTopo(&tbody, topo);

    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_EDGE,   ebody.edges.map);
    nent = ebody.edges.map.Extent();
    ebody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;
    
    egadsLoop *ploop = new egadsLoop;
    TopoDS_Wire Loop = TopoDS::Wire(nTopo);
    ploop->loop      = Loop;
    ploop->surface   = NULL;
    obj->blind       = ploop;
    obj->oclass      = LOOP;
    EG_flipAttrTopo(&ebody, &tbody, topo, obj, obj);
    EG_attributeDup(topo, obj);
  
    delete [] ebody.edges.objs;
    delete [] ebody.nodes.objs;

    delete [] tbody.edges.objs;
    delete [] tbody.nodes.objs;
  
  } else if (topo->oclass == FACE) {
  
    TopExp::MapShapes(shape, TopAbs_VERTEX, tbody.nodes.map);
    nent = tbody.nodes.map.Extent();
    tbody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.nodes.objs[i] = NULL;
    TopExp::MapShapes(shape, TopAbs_EDGE,   tbody.edges.map);
    nent = tbody.edges.map.Extent();
    tbody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.edges.objs[i] = NULL;
    TopExp::MapShapes(shape, TopAbs_WIRE,   tbody.loops.map);
    nent = tbody.loops.map.Extent();
    tbody.loops.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.loops.objs[i] = NULL;
    EG_fillObjTopo(&tbody, topo);

    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_EDGE,   ebody.edges.map);
    nent = ebody.edges.map.Extent();
    ebody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_WIRE,   ebody.loops.map);
    nent = ebody.loops.map.Extent();
    ebody.loops.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.loops.objs[i] = NULL;
    
    egadsFace *pface = new egadsFace;
    TopoDS_Face Face = TopoDS::Face(nTopo);
    pface->face      = Face;
    obj->blind       = pface;
    obj->oclass      = FACE;
    EG_flipAttrTopo(&ebody, &tbody, topo, obj, obj);
    EG_attributeDup(topo, obj);
  
    delete [] ebody.loops.objs;
    delete [] ebody.edges.objs;
    delete [] ebody.nodes.objs;

    delete [] tbody.loops.objs;
    delete [] tbody.edges.objs;
    delete [] tbody.nodes.objs;

  } else {
  
    TopExp::MapShapes(shape, TopAbs_VERTEX, tbody.nodes.map);
    nent = tbody.nodes.map.Extent();
    tbody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.nodes.objs[i] = NULL;
    TopExp::MapShapes(shape, TopAbs_EDGE,   tbody.edges.map);
    nent = tbody.edges.map.Extent();
    tbody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.edges.objs[i] = NULL;
    TopExp::MapShapes(shape, TopAbs_WIRE,   tbody.loops.map);
    nent = tbody.loops.map.Extent();
    tbody.loops.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.loops.objs[i] = NULL;
    TopExp::MapShapes(shape, TopAbs_FACE,   tbody.faces.map);
    nent = tbody.faces.map.Extent();
    tbody.faces.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) tbody.faces.objs[i] = NULL;
    EG_fillObjTopo(&tbody, topo);

    TopExp::MapShapes(nTopo, TopAbs_VERTEX, ebody.nodes.map);
    nent = ebody.nodes.map.Extent();
    ebody.nodes.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.nodes.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_EDGE,   ebody.edges.map);
    nent = ebody.edges.map.Extent();
    ebody.edges.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.edges.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_WIRE,   ebody.loops.map);
    nent = ebody.loops.map.Extent();
    ebody.loops.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.loops.objs[i] = NULL;
    TopExp::MapShapes(nTopo, TopAbs_FACE,   ebody.faces.map);
    nent = ebody.faces.map.Extent();
    ebody.faces.objs = new egObject*[nent];
    for (i = 0; i < nent; i++) ebody.faces.objs[i] = NULL;
    
    egadsShell *pshell = new egadsShell;
    TopoDS_Shell Shell = TopoDS::Shell(nTopo);
    pshell->shell      = Shell;
    obj->blind         = pshell;
    obj->oclass        = SHELL;
    EG_flipAttrTopo(&ebody, &tbody, topo, obj, obj);
    EG_attributeDup(topo, obj);
  
    delete [] ebody.faces.objs;
    delete [] ebody.loops.objs;
    delete [] ebody.edges.objs;
    delete [] ebody.nodes.objs;

    delete [] tbody.faces.objs;
    delete [] tbody.loops.objs;
    delete [] tbody.edges.objs;
    delete [] tbody.nodes.objs;

  }

  EG_referenceObject(obj, context);
  *copy = obj;
  return EGADS_SUCCESS;
}


#ifndef BBOX
static void
EG_edgeBBox(TopoDS_Edge edge, double *ebx)
{
  TopoDS_Vertex V1, V2;
  Standard_Real t1, t2, t;
  gp_Pnt        xyz;

  if (edge.Orientation() == TopAbs_REVERSED) {
    TopExp::Vertices(edge, V2, V1, Standard_True);
  } else {
    TopExp::Vertices(edge, V1, V2, Standard_True);
  }
  xyz = BRep_Tool::Pnt(V1);
  ebx[0] = ebx[3] = xyz.X();
  ebx[1] = ebx[4] = xyz.Y();
  ebx[2] = ebx[5] = xyz.Z();
  xyz = BRep_Tool::Pnt(V2);
  if (xyz.X() < ebx[0]) ebx[0] = xyz.X();
  if (xyz.X() > ebx[3]) ebx[3] = xyz.X();
  if (xyz.Y() < ebx[1]) ebx[1] = xyz.Y();
  if (xyz.Y() > ebx[4]) ebx[4] = xyz.Y();
  if (xyz.Z() < ebx[2]) ebx[2] = xyz.Z();
  if (xyz.Z() > ebx[5]) ebx[5] = xyz.Z();
  Handle(Geom_Curve) hCurve = BRep_Tool::Curve(edge, t1, t2);
  
  for (int i = 1; i <= 12; i++) {
    t = t1 + i*(t2-t1)/13;
    hCurve->D0(t, xyz);
    if (xyz.X() < ebx[0]) ebx[0] = xyz.X();
    if (xyz.X() > ebx[3]) ebx[3] = xyz.X();
    if (xyz.Y() < ebx[1]) ebx[1] = xyz.Y();
    if (xyz.Y() > ebx[4]) ebx[4] = xyz.Y();
    if (xyz.Z() < ebx[2]) ebx[2] = xyz.Z();
    if (xyz.Z() > ebx[5]) ebx[5] = xyz.Z();
  }
}
#endif


int
EG_matchBodyFaces(const egObject *body1, const egObject *body2, double toler,
                  int *nmatch, int **match)
{
  int i, j, k, l, n, hit;

  *nmatch = 0;
  *match  = NULL;
  if (body1 == NULL)               return EGADS_NULLOBJ;
  if (body1->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body1->oclass != BODY)       return EGADS_NOTBODY;
  if (body1->blind == NULL)        return EGADS_NODATA;
  if (body2 == NULL)               return EGADS_NULLOBJ;
  if (body2->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (body2->oclass != BODY)       return EGADS_NOTBODY;
  if (body2->blind == NULL)        return EGADS_NODATA;
  
  int outLevel = EG_outLevel(body1);
  if (EG_context(body1) != EG_context(body2)) {
    if (outLevel > 0)
      printf(" EGADS Error: Context mismatch (EG_matchBodyFaces)!\n");
    return EGADS_MIXCNTX;
  }
  egadsBody *pbod1 = (egadsBody *) body1->blind;
  egadsBody *pbod2 = (egadsBody *) body2->blind;
  int       nface1 = pbod1->faces.map.Extent();
  int       nface2 = pbod2->faces.map.Extent();

  int       *map1  = (int *) EG_alloc(nface1*sizeof(int));
  if (map1 == NULL) return EGADS_MALLOC;
  for (i = 0; i < nface1; i++) map1[i] = -1;
  
  /* check each face in body1 against all of body2 */
  for (i = 0; i < nface1; i++) {
    TopTools_IndexedMapOfShape nmap1, emap1, lmap1;
    TopoDS_Shape shape1 = pbod1->faces.map(i+1);
    TopoDS_Face  face1  = TopoDS::Face(shape1);
    TopExp::MapShapes(shape1, TopAbs_VERTEX, nmap1);
    TopExp::MapShapes(shape1, TopAbs_EDGE,   emap1);
    TopExp::MapShapes(shape1, TopAbs_WIRE,   lmap1);
#ifdef BBOX
    double ftol1        = BRep_Tool::Tolerance(face1);
    Bnd_Box fbox1;
    BRepBndLib::Add(shape1, fbox1);
    double fbx1[6];
    fbox1.Get(fbx1[0], fbx1[1], fbx1[2], fbx1[3], fbx1[4], fbx1[5]);
#endif

    for (j = 0; j < nface2; j++) {
      TopTools_IndexedMapOfShape nmap2, emap2, lmap2;
      TopoDS_Shape shape2 = pbod2->faces.map(j+1);
      TopoDS_Face  face2  = TopoDS::Face(shape2);
#ifdef BBOX
      double ftol         = BRep_Tool::Tolerance(face2);
      Bnd_Box fbox2;
      BRepBndLib::Add(shape2, fbox2);
      double fbx2[6];
      fbox2.Get(fbx2[0], fbx2[1], fbx2[2], fbx2[3], fbx2[4], fbx2[5]);
  
      if (ftol < ftol1) ftol = ftol1;
      if (toler != 0.0) ftol = toler;
      double ll = sqrt((fbx2[0]-fbx1[0])*(fbx2[0]-fbx1[0]) +
                       (fbx2[1]-fbx1[1])*(fbx2[1]-fbx1[1]) +
                       (fbx2[2]-fbx1[2])*(fbx2[2]-fbx1[2]));
      double ur = sqrt((fbx2[3]-fbx1[3])*(fbx2[3]-fbx1[3]) +
                       (fbx2[4]-fbx1[4])*(fbx2[4]-fbx1[4]) +
                       (fbx2[5]-fbx1[5])*(fbx2[5]-fbx1[5]));
      if ((ll > ftol) || (ur > ftol)) continue;
      if (outLevel > 1)
        printf(" EGADS Info: Faces %d and %d pass  Face  check!\n",
               i+1, j+1);
#endif

      /* loops */
      TopExp::MapShapes(shape2, TopAbs_WIRE,   lmap2);
      if (lmap1.Extent() != lmap2.Extent()) continue;
      if (outLevel > 1)
        printf(" EGADS Info: Faces %d and %d pass #Loops check!\n",
               i+1, j+1);

      /* nodes */
      TopExp::MapShapes(shape2, TopAbs_VERTEX, nmap2);
      if (nmap1.Extent() != nmap2.Extent()) continue;
      if (outLevel > 1)
        printf(" EGADS Info: Faces %d and %d pass #Nodes check!\n",
               i+1, j+1);
      for (hit = k = 0; k < nmap1.Extent(); k++) {
        TopoDS_Shape  shape = nmap1(k+1);
        TopoDS_Vertex vert  = TopoDS::Vertex(shape);
        gp_Pnt pv1          = BRep_Tool::Pnt(vert);
        double tol1         = BRep_Tool::Tolerance(vert);
        for (l = 0; l < nmap2.Extent(); l++) {
          TopoDS_Shape  shapv = nmap2(l+1);
          TopoDS_Vertex vrt   = TopoDS::Vertex(shapv);
          gp_Pnt pv2          = BRep_Tool::Pnt(vrt);
          double tol          = BRep_Tool::Tolerance(vrt);
          if (tol   < tol1) tol = tol1;
          if (toler != 0.0) tol = toler;
          double dist = sqrt((pv2.X()-pv1.X())*(pv2.X()-pv1.X()) +
                             (pv2.Y()-pv1.Y())*(pv2.Y()-pv1.Y()) +
                             (pv2.Z()-pv1.Z())*(pv2.Z()-pv1.Z()));
          if (dist <= tol) hit++;
        }
      }
      if (hit != nmap1.Extent()) continue;
      if (outLevel > 1)
        printf(" EGADS Info: Faces %d and %d pass  Node  check!\n", i+1, j+1);
      
      /* edges */
      TopExp::MapShapes(shape2, TopAbs_EDGE,   emap2);
      if (emap1.Extent() != emap2.Extent()) continue;
      if (outLevel > 1)
        printf(" EGADS Info: Faces %d and %d pass #Edges check!\n",
               i+1, j+1);
      for (hit = k = 0; k < emap1.Extent(); k++) {
        TopoDS_Shape shape = emap1(k+1);
        TopoDS_Edge  edge1 = TopoDS::Edge(shape);
        if (BRep_Tool::Degenerated(edge1)) {
          hit++;
          continue;
        }
        double etol1       = BRep_Tool::Tolerance(edge1);
        BRepGProp    BProps;
        GProp_GProps SProps1;
        BProps.LinearProperties(edge1, SProps1);
        double ebx1[6];
#ifdef BBOX
        Bnd_Box ebox1;
        BRepBndLib::Add(shape, ebox1);
        ebox1.Get(ebx1[0], ebx1[1], ebx1[2], ebx1[3], ebx1[4], ebx1[5]);
#else
        EG_edgeBBox(edge1, ebx1);
#endif
        for (l = 0; l < emap2.Extent(); l++) {
          TopoDS_Shape shapx = emap2(l+1);
          TopoDS_Edge  edge2 = TopoDS::Edge(shapx);
          if (BRep_Tool::Degenerated(edge2)) continue;
          double etol        = BRep_Tool::Tolerance(edge2);
          if (etol < etol1) etol = etol1;
          if (toler != 0.0) etol = toler;
          GProp_GProps SProps2;
          BProps.LinearProperties(edge2, SProps2);
          double ebx2[6];
#ifdef BBOX
          Bnd_Box ebox2;
          BRepBndLib::Add(shapx, ebox2);
          ebox2.Get(ebx2[0], ebx2[1], ebx2[2], ebx2[3], ebx2[4], ebx2[5]);
#else
          EG_edgeBBox(edge2, ebx2);
#endif
          double ll = sqrt((ebx2[0]-ebx1[0])*(ebx2[0]-ebx1[0]) +
                           (ebx2[1]-ebx1[1])*(ebx2[1]-ebx1[1]) +
                           (ebx2[2]-ebx1[2])*(ebx2[2]-ebx1[2]));
          double ur = sqrt((ebx2[3]-ebx1[3])*(ebx2[3]-ebx1[3]) +
                           (ebx2[4]-ebx1[4])*(ebx2[4]-ebx1[4]) +
                           (ebx2[5]-ebx1[5])*(ebx2[5]-ebx1[5]));
          if ((ll > etol) || (ur > etol)) continue;
          if (fabs(SProps1.Mass()-SProps2.Mass()) > etol) continue;
          hit++;
          break;
        }
      }
      if (hit != emap1.Extent()) continue;
      if (outLevel > 1)
        printf(" EGADS Info: Faces %d and %d pass  Edge  checks!\n", i+1, j+1);
      
      map1[i] = j;
      break;
    }
  }
  
  /* collect the results and return */
  for (n = i = 0; i < nface1; i++)
    if (map1[i] != -1) n++;
  
  if (n != 0) {
    int *fill = (int *) EG_alloc(2*n*sizeof(int));
    if (fill == NULL) {
      EG_free(map1);
      return EGADS_MALLOC;
    }
    for (n = i = 0; i < nface1; i++)
      if (map1[i] != -1) {
        fill[2*n  ] = i+1;
        fill[2*n+1] = map1[i]+1;
        n++;
      }
    *nmatch = n;
    *match  = fill;
  }
  EG_free(map1);

  return EGADS_SUCCESS;
}
