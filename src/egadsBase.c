/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Base Object Functions
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef WIN32
#define LONG long
#include <execinfo.h>
#else
#ifdef _OCC64
#define LONG long long
#else
#define LONG long
#endif
#endif

#include "egadsTypes.h"
#include "egadsInternals.h"


#define ZERO            1.e-5           /* allow for float-like precision */
#define STRING(a)       #a
#define STR(a)          STRING(a)


  static char *EGADSprop[2] = {STR(EGADSPROP),
                               "\nEGADSprop: Copyright 2011-2016 MIT. All Rights Reserved."};


  extern void EG_initOCC( );
  extern void EG_exactInit( );
  extern int  EG_destroyGeometry( egObject *geom );
  extern int  EG_destroyTopology( egObject *topo );
  extern int  EG_copyGeometry( const egObject *geom, /*@null@*/ double *xform,
                               egObject **copy );
  extern int  EG_copyTopology( const egObject *topo, /*@null@*/ double *xform,
                               egObject **copy );
  extern int  EG_flipGeometry( const egObject *geom, egObject **copy );
  extern int  EG_flipTopology( const egObject *topo, egObject **copy );
  extern int  EG_getTopology( const egObject *topo, egObject **geom, 
                              int *ocls, int *type, /*@null@*/ double *limits, 
                              int *nobjs, egObject ***objs, int **senses );
  extern int  EG_getGeometry( const egObject *geom, int *oclass, int *mtype,
                              egObject **refGeom, int **ivec, double **rvec );
  extern int  EG_getTolerance( const egObject *topo, double *tol );

  extern int  EG_computeTessMap( egTessel *btess, int outLevel );
  extern int  EG_initTessBody( egObject *object, egObject **tess );
  extern int  EG_getTessEdge( const egObject *tess, int eIndex, int *len,
                              const double **xyz, const double **t );
  extern int  EG_setTessEdge( egObject *tess, int eIndex, int len,
                              const double *xyz, const double *t );
  extern int  EG_getTessFace( const egObject *tess, int fIndex, int *len,
                              const double **xyz, const double **uv,
                              const int **ptype, const int **pindex, int *ntri,
                              const int **tris, const int **tric );
  extern int  EG_setTessFace( egObject *tess, int fIndex, int len,
                              const double *xyz, const double *uv, int ntri,
                              const int *tris );



static void
EG_traceback()
{
#ifndef WIN32
  int    i;
  void   *array[100];
  size_t size;
  char   **strings;

  size = backtrace(array, 100);
  strings = backtrace_symbols (array, size);
  i = size;
  printf ("\nObtained %d stack frames:\n", i);
  for (i = 0; i < size; i++)
    printf ("%s\n", strings[i]);
  free (strings);
  printf("\n");
#endif
}


/*@kept@*/ /*@null@*/ egObject *
EG_context(const egObject *obj)
{
  int      cnt;
  egObject *object, *topObj;

  if (obj == NULL) {
    printf(" EGADS Internal: EG_context called with NULL!\n");
    return NULL;
  }
  if (obj->magicnumber != MAGIC) {
    printf(" EGADS Internal: EG_context Object NOT an ego!\n");
    return NULL;
  }
  if (obj->oclass == CONTXT) return (egObject *) obj;
  
  object = obj->topObj;
  if (object == NULL) {
    printf(" EGADS Internal: EG_context topObj is NULL!\n");
    return NULL;
  }
  if (object->magicnumber != MAGIC) {
    printf(" EGADS Internal: EG_context topObj NOT an ego!\n");
    return NULL;
  }
  if (object->oclass == CONTXT) return object;
  
  cnt = 0;
  do {
    topObj = object->topObj;
    if (topObj == NULL) {
      printf(" EGADS Internal: %d EG_context contents of topObj is NULL!\n",
             cnt);
      return NULL;
    }
    if (topObj->magicnumber != MAGIC) {
      printf(" EGADS Internal: %d EG_context contents of topObj NOT an ego!\n",
             cnt);
      return NULL;
    }
    if (topObj->oclass == CONTXT) return topObj;
    object = topObj;
    cnt++;
  } while (object != NULL);

  printf(" EGADS Internal: Cannot find context -- depth = %d!\n", cnt);
  EG_traceback();
  return NULL;
}


int
EG_outLevel(const egObject *obj)
{
  egObject *context;
  egCntxt  *cntxt;

  if (obj == NULL)               return 0;
  if (obj->magicnumber != MAGIC) return 0;
  context = EG_context(obj);
  if (context == NULL)           return 0;
  
  cntxt = (egCntxt *) context->blind;
  return cntxt->outLevel;
}


int
EG_setOutLevel(egObject *context, int outLevel)
{
  int     old;
  egCntxt *cntx;

  if  (context == NULL)                 return EGADS_NULLOBJ;
  if  (context->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (context->oclass != CONTXT)       return EGADS_NOTCNTX;
  if ((outLevel < 0) || (outLevel > 3)) return EGADS_RANGERR;
  cntx = (egCntxt *) context->blind;
  if  (cntx == NULL)                    return EGADS_NODATA;
  old            = cntx->outLevel;
  cntx->outLevel = outLevel;
  
  return old;
}


int
EG_setTessParam(egObject *context, int iParam, double value, double *oldValue)
{
  egCntxt *cntx;
  
  *oldValue = 0.0;
  if  (context == NULL)                      return EGADS_NULLOBJ;
  if  (context->magicnumber != MAGIC)        return EGADS_NOTOBJ;
  if  (context->oclass != CONTXT)            return EGADS_NOTCNTX;
  if ((iParam < 1) || (iParam > MTESSPARAM)) return EGADS_RANGERR;
  cntx = (egCntxt *) context->blind;
  if  (cntx == NULL)                          return EGADS_NODATA;
  *oldValue            = cntx->tess[iParam-1];
  cntx->tess[iParam-1] = value;
  
  return EGADS_SUCCESS;
}


int
EG_makeObject(/*@null@*/ egObject *context, egObject **obj)
{
  int      outLevel;
  egObject *object, *prev;
  egCntxt  *cntx;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                  return EGADS_NODATA;
  outLevel = cntx->outLevel;

  /* any objects in the pool? */
  object = cntx->pool;
  if (object == NULL) {
    object = (egObject *) EG_alloc(sizeof(egObject));
    if (object == NULL) {
      if (outLevel > 0) 
        printf(" EGADS Error: Malloc on Object (EG_makeObject)!\n");
      return EGADS_MALLOC;
    }
  } else {
    cntx->pool   = object->next;
    object->prev = NULL;
  }
  
  prev                = cntx->last;
  object->magicnumber = MAGIC;
  object->oclass      = NIL;
  object->mtype       = 0;
  object->tref        = NULL;
  object->attrs       = NULL;
  object->blind       = NULL;
  object->topObj      = context;
  object->prev        = prev;
  object->next        = NULL;
  prev->next          = object;

  *obj = object;
  cntx->last = *obj;
  return EGADS_SUCCESS;
}


int
EG_open(egObject **context)
{
  int      i;
  egObject *object;
  egCntxt  *cntx;

  cntx   = (egCntxt *) EG_alloc(sizeof(egCntxt));
  if (cntx == NULL) return EGADS_MALLOC;
  object = (egObject *) EG_alloc(sizeof(egObject));
  if (object == NULL) {
    EG_free(cntx);
    return EGADS_MALLOC;
  }
  for (i = 0; i < MTESSPARAM; i++) cntx->tess[i] = 0.0;
  cntx->outLevel  = 1;
  cntx->signature = EGADSprop;
  cntx->pool      = NULL;
  cntx->last      = object;
  
  object->magicnumber = MAGIC;
  object->oclass      = CONTXT;
  object->mtype       = 0;
  object->tref        = NULL;
  object->attrs       = NULL;
  object->blind       = cntx;
  object->topObj      = NULL;
  object->prev        = NULL;
  object->next        = NULL;

  EG_initOCC();
  EG_exactInit();
  *context = object;
  return EGADS_SUCCESS;
}


int
EG_referenceObject(egObject *object, /*@null@*/ const egObject *ref)
{
  int      cnt, stat, outLevel;
  egObject *next, *last, *obj, *ocontext, *rcontext;
  
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == NIL)        return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  outLevel = EG_outLevel(object);

  if (ref == NULL) {
    if (outLevel > 0) 
      printf(" EGADS Error: NULL Reference (EG_referenceObject)!\n");
    return EGADS_NULLOBJ;
  }
  if (ref->magicnumber != MAGIC) {
    if (outLevel > 0) 
      printf(" EGADS Error: Reference not an EGO (EG_referenceObject)!\n");
    return EGADS_NOTOBJ;
  }
  if ((ref->oclass == EMPTY) || (ref->oclass == NIL)) {
    if (outLevel > 0) 
      printf(" EGADS Error: Reference is Empty (EG_refreneceObject)!\n");
    return EGADS_EMPTY;
  }
  ocontext = EG_context(object);
  rcontext = EG_context(ref);
  if (rcontext != ocontext) {
    if (outLevel > 0) 
      printf(" EGADS Error: Context mismatch (EG_referenceObject)!\n");
    return EGADS_MIXCNTX;
  }

  cnt = 1;
  obj = NULL;
  if (object->tref == NULL) {
    stat = EG_makeObject(ocontext, &obj);
    if (outLevel > 2) 
      printf(" 0 makeRef oclass %d for rclass %d = %d\n", 
             object->oclass, ref->oclass, stat);
    if (stat != EGADS_SUCCESS) return stat;
    if (obj != NULL) {
      obj->oclass  = REFERENCE;
      obj->attrs   = (void *) ref;
      object->tref = obj;
    }
  } else {
    next = object->tref;
    while (next != NULL) {
    if (outLevel > 2) {
      if (next->magicnumber != MAGIC)
        printf(" %d: Thread not an EGO!\n", cnt);
      if (next->oclass != REFERENCE)
        printf(" %d: Not a Reference - class = %d!\n", cnt, next->oclass);
      }
      obj = (egObject *) next->attrs;
      if (outLevel > 2) {
        if (obj == NULL) {
          printf(" %d: Reference is NULL!\n", cnt);
        } else {
          if (obj->magicnumber != MAGIC)
            printf(" %d: Reference not an EGO!\n", cnt);
          if ((obj->oclass == EMPTY) || (obj->oclass == NIL))
            printf(" %d: Reference is Empty!\n", cnt);
        }
      }
/*    single node edges do double reference!
      if (obj == ref) return EGADS_SUCCESS;  */
      last = next;
      next = (egObject *) last->blind;          /* next reference */
      cnt++;
    }
    stat = EG_makeObject(ocontext, &obj);
    if (outLevel > 2) 
      printf(" %d makeRef oclass %d for rclass %d = %d\n", 
             cnt, object->oclass, ref->oclass, stat);
    if (stat != EGADS_SUCCESS) return stat;
    if (obj != NULL) {
      obj->oclass = REFERENCE;
      obj->attrs  = (void *) ref;
      last->blind = obj;
    }
  }

  return cnt;
}


int
EG_referenceObjects(egObject *object, int *nobj, egObject ***objs)
{
  int      i;
  egObject *next, **objects;
  
  *nobj = 0;
  *objs = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;

  i    = 0;
  next = object->tref;
  while (next != NULL) {
    next = (egObject *) next->blind;          /* next reference */
    i++;
  }
  if (i == 0) return EGADS_SUCCESS;
  
  objects = (egObject **) EG_alloc(i*sizeof(egObject *));
  if (objects == NULL) return EGADS_MALLOC;

  i    = 0;
  next = object->tref;
  while (next != NULL) {
    objects[i] = (egObject *) next->attrs;
    next       = (egObject *) next->blind;
    i++;
  }
  
  *nobj = i;
  *objs = objects;

  return EGADS_SUCCESS;
}


int
EG_referenceTopObj(egObject *object, /*@null@*/ const egObject *ref)
{
  egObject *context, *obj;
  
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == NIL)        return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  context = EG_context(object);
  obj     = object;
  if (object->topObj != context) obj = object->topObj;
  
  return EG_referenceObject(obj, ref);
}


static int
EG_derefObj(egObject *object, /*@null@*/ const egObject *refx, int flg)
{
  int      i, j, stat, outLevel;
  LONG     ptr1, ptr2;
  egObject *pobj, *nobj, *obj, *context;
  egCntxt  *cntx;
  egTessel *tess;
  const egObject *ref;

  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  context = EG_context(object);
  if (context == NULL)              return EGADS_NOTCNTX;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                 return EGADS_NODATA;
  outLevel = cntx->outLevel;
  ref      = refx;

  /* context is an attempt to delete */
  
  if ((ref == context) && (object->tref != NULL)) {
    nobj = object->tref;
    pobj = NULL;
    i = 0;
    while (nobj != NULL) {
      obj = (egObject *) nobj->attrs;
      if (obj != ref) i++;
      pobj = nobj;
      nobj = (egObject *) pobj->blind;          /* next reference */
    }
    if (object->topObj == context)
      if (i > 0) {
        if (outLevel > 0)
          printf(" EGADS Info: dereference with %d active objects!\n", i); 
        return i;
      }
  }
  if (ref == NULL) ref = context;
  
  /* we should never see a NULL reference! */
  if (object->tref != NULL) {
    nobj = object->tref;
    pobj = NULL;
    while (nobj != NULL) {
      obj = (egObject *) nobj->attrs;
      if (obj == ref) break;
      pobj = nobj;
      nobj = (egObject *) pobj->blind;          /* next reference */
    }
    if (nobj == NULL) {
      if (refx != NULL) {
        ptr1 = (LONG) object;
        ptr2 = (LONG) ref;
        printf(" EGADS Internal: Ref Not Found (EG_dereferenceObject)!\n");
#if defined(WIN32) && defined(_OCC64)
        printf("                 Object %llx = %d/%d,  ref %llx = %d/%d\n",
               ptr1, object->oclass, object->mtype,
               ptr2, ref->oclass, ref->mtype);
#else
        printf("                 Object %lx = %d/%d,  ref %lx = %d/%d\n",
               ptr1, object->oclass, object->mtype,
               ptr2, ref->oclass, ref->mtype);
#endif
      }
      return EGADS_NOTFOUND;
    }
    if (pobj == NULL) {
      object->tref = (egObject *) nobj->blind;
    } else {
      pobj->blind  = nobj->blind;
    }
    obj  = nobj;
    pobj = obj->prev;
    nobj = obj->next;
    if (nobj != NULL) {
      nobj->prev = pobj;
    } else {
      cntx->last = pobj;
    }
    pobj->next   = nobj;
    obj->mtype   = REFERENCE;
    obj->oclass  = EMPTY;
    obj->blind   = NULL;
    obj->topObj  = context;
    obj->prev    = NULL;
    obj->next    = cntx->pool;
    cntx->pool   = obj;  
  }
  if (object->tref != NULL) return EGADS_SUCCESS;

  stat = EG_attributeDel(object, NULL);
  if (stat != EGADS_SUCCESS)
    if (outLevel > 0)
      printf(" EGADS Warning: Del Attributes = %d (EG_destroyObject)!\n",
             stat);

  stat = EGADS_SUCCESS;
  if (object->oclass == TRANSFORM) {
  
    EG_free(object->blind);

  } else if (object->oclass == TESSELLATION) {

    tess = (egTessel *) object->blind;
    if (tess != NULL) {
      EG_dereferenceTopObj(tess->src, object);
      if (tess->xyzs != NULL) EG_free(tess->xyzs);
      if (tess->tess1d != NULL) {
        for (i = 0; i < tess->nEdge; i++) {
          if (tess->tess1d[i].faces[0].faces != NULL)
            EG_free(tess->tess1d[i].faces[0].faces);
          if (tess->tess1d[i].faces[1].faces != NULL)
            EG_free(tess->tess1d[i].faces[1].faces);
          if (tess->tess1d[i].faces[0].tric  != NULL)
            EG_free(tess->tess1d[i].faces[0].tric);
          if (tess->tess1d[i].faces[1].tric  != NULL)
            EG_free(tess->tess1d[i].faces[1].tric);
          if (tess->tess1d[i].xyz    != NULL)
            EG_free(tess->tess1d[i].xyz);
          if (tess->tess1d[i].t      != NULL)
            EG_free(tess->tess1d[i].t);
          if (tess->tess1d[i].global != NULL)
            EG_free(tess->tess1d[i].global);
        }
        EG_free(tess->tess1d);
      }
      if (tess->tess2d != NULL) {
        for (i = 0; i < 2*tess->nFace; i++) {
          if (tess->tess2d[i].mKnots != NULL) 
            EG_deleteObject(tess->tess2d[i].mKnots);
          if (tess->tess2d[i].xyz    != NULL) 
            EG_free(tess->tess2d[i].xyz);
          if (tess->tess2d[i].uv     != NULL) 
            EG_free(tess->tess2d[i].uv);
          if (tess->tess2d[i].global != NULL)
            EG_free(tess->tess2d[i].global);
          if (tess->tess2d[i].ptype  != NULL) 
            EG_free(tess->tess2d[i].ptype);
          if (tess->tess2d[i].pindex != NULL) 
            EG_free(tess->tess2d[i].pindex);
          if (tess->tess2d[i].bary   != NULL)
            EG_free(tess->tess2d[i].bary);
          if (tess->tess2d[i].frame != NULL)
            EG_free(tess->tess2d[i].frame);
          if (tess->tess2d[i].frlps != NULL)
            EG_free(tess->tess2d[i].frlps);
          if (tess->tess2d[i].tris   != NULL) 
            EG_free(tess->tess2d[i].tris);
          if (tess->tess2d[i].tric   != NULL) 
            EG_free(tess->tess2d[i].tric);
          if (tess->tess2d[i].patch  != NULL) {
            for (j = 0; j < tess->tess2d[i].npatch; j++) {
              if (tess->tess2d[i].patch[j].ipts != NULL) 
                EG_free(tess->tess2d[i].patch[j].ipts);
              if (tess->tess2d[i].patch[j].bounds != NULL) 
                EG_free(tess->tess2d[i].patch[j].bounds);
            }
            EG_free(tess->tess2d[i].patch);
          }
        }
        EG_free(tess->tess2d);
      }
      if (tess->globals != NULL) EG_free(tess->globals);
      EG_free(tess);
    }

  } else if (object->oclass <= SURFACE) {
  
    if ((object->oclass != NIL) && (flg == 0))
      stat = EG_destroyGeometry(object);
  
  } else {
  
    if (flg == 0) stat = EG_destroyTopology(object);

  }
  object->mtype  = object->oclass;
  object->oclass = EMPTY;
  object->blind  = NULL;
  
  /* patch up the lists & put the object in the pool */

  pobj = object->prev;          /* always have a previous -- context! */
  nobj = object->next;
  if (nobj == NULL) {
    if (object != cntx->last)
      printf(" EGADS Info: Context Last NOT Object Next w/ NULL!\n");
    cntx->last = pobj;
  } else {
    nobj->prev = pobj;
  }
  if (pobj == NULL) {
    printf(" EGADS Info: PrevObj is NULL (EG_destroyObject)!\n");
  } else {
    pobj->next = nobj;
  }
  object->prev = NULL;
  object->next = cntx->pool;
  cntx->pool   = object;

  return stat;
}


int
EG_dereferenceTopObj(egObject *object, /*@null@*/ const egObject *ref)
{
  egObject *context, *obj;
  
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == NIL)        return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  context = EG_context(object);
  obj     = object;
  if (object->topObj != context) obj = object->topObj;
  
  return EG_derefObj(obj, ref, 1);
}


int
EG_dereferenceObject(egObject *object, /*@null@*/ const egObject *ref)
{
  return EG_derefObj(object, ref, 0);
}


int
EG_deleteObject(egObject *object)
{
  int      outLevel, total, cnt, nref, stat, oclass, mtype, nbody;
  int      i, *senses;
  egObject *context, *obj, *next, *pobj, **bodies;
  egCntxt  *cntx;

  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  if (object->oclass != CONTXT) {
    /* normal dereference */
    context = EG_context(object);
    if (context == NULL)            return EGADS_NOTCNTX;
    outLevel = EG_outLevel(object);
    
    /* special check for body references in models */
    if (object->oclass == MODEL) {
      stat = EG_getTopology(object, &next, &oclass, &mtype, NULL, 
                            &nbody, &bodies, &senses);
      if (stat != EGADS_SUCCESS) return stat;
      for (cnt = i = 0; i < nbody; i++) {
        if (bodies[i]->tref == NULL) continue;
        next = bodies[i]->tref;
        pobj = NULL;
        while (next != NULL) {
          obj = (egObject *) next->attrs;
          if (obj != object) cnt++;
          pobj = next;
          next = (egObject *) pobj->blind;      /* next reference */
        }
      }
      if (cnt > 0) {
        if (outLevel > 0)
          printf(" EGADS Info: Model delete w/ %d active Body Refs!\n", 
                 cnt); 
        return cnt;
      }
    }
  
    return EG_dereferenceObject(object, context);
  }
  
  /* delete all non-body attached topology and geometry */
 
  context  = object; 
  cntx     = (egCntxt *) context->blind;
  if (cntx == NULL) return EGADS_NODATA;
  outLevel = cntx->outLevel;
  
  nref = 0;
  obj  = context->next;
  while (obj != NULL) {
    next = obj->next;
    if (obj->oclass == REFERENCE) nref++;
    obj  = next;
  }

  cntx->outLevel = total = 0;
  do {
    cnt = 0;
    obj = context->next;
    while (obj != NULL) {
      next = obj->next;
      if ((obj->oclass >= PCURVE) && (obj->oclass <= SHELL) &&
          (obj->topObj == context)) {
        stat = EG_dereferenceObject(obj, context);
        if (stat == EGADS_SUCCESS) {
          cnt++;
          break;
        }
      }
      obj = next;
    }
    total += cnt;
  } while (cnt != 0);
  cntx->outLevel = outLevel;
  
  cnt = 0;
  obj = context->next;
  while (obj != NULL) {
    next = obj->next;
    if (obj->oclass == REFERENCE) cnt++;
    obj  = next;
  }
  
  if ((outLevel > 0) && (total != 0))
    printf(" EGADS Info: %d unattached Objects (%d References) removed!\n",
           total, nref-cnt);

  return EGADS_SUCCESS;
}


int
EG_removeCntxtRef(egObject *object)
{
  egObject *context, *nobj, *pobj, *obj;
  egCntxt  *cntx;
  
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == NIL)        return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  if (object->tref   == NULL)       return EGADS_SUCCESS;
  context = EG_context(object);
  if (context == NULL)              return EGADS_NULLOBJ;
  cntx    = (egCntxt *) context->blind;
  if (cntx == NULL)                 return EGADS_NODATA;

  nobj = object->tref;
  pobj = NULL;
  while (nobj != NULL) {
    obj = (egObject *) nobj->attrs;
    if (obj == context) break;
    pobj = nobj;
    nobj = (egObject *) pobj->blind;          /* next reference */
  }
  if (nobj == NULL) return EGADS_NOTFOUND;
  if (pobj == NULL) {
    object->tref = (egObject *) nobj->blind;
  } else {
    pobj->blind  = nobj->blind;
  }
  obj  = nobj;
  pobj = obj->prev;
  nobj = obj->next;
  if (nobj != NULL) {
    nobj->prev = pobj;
  } else {
    cntx->last = pobj;
  }
  pobj->next   = nobj;
  obj->mtype   = REFERENCE;
  obj->oclass  = EMPTY;
  obj->blind   = NULL;
  obj->topObj  = context;
  obj->prev    = NULL;
  obj->next    = cntx->pool;
  /*@ignore@*/ 
  cntx->pool   = obj;
  /*@end@*/

  return EGADS_SUCCESS;
}


int
EG_makeTransform(egObject *context, const double *xform, egObject **oform)
{
  int      i, stat, outLevel;
  double   *reals, dotXY, dotXZ, dotYZ, dotXX, dotYY, dotZZ;
  egObject *object;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (xform == NULL)                 return EGADS_NODATA;
  outLevel = EG_outLevel(context);

  /* check for "scaled" orthonormality */

  dotXX = xform[0]*xform[0] + xform[1]*xform[1] + xform[ 2]*xform[ 2];
  dotXY = xform[0]*xform[4] + xform[1]*xform[5] + xform[ 2]*xform[ 6];
  dotXZ = xform[0]*xform[8] + xform[1]*xform[9] + xform[ 2]*xform[10];

  dotYY = xform[4]*xform[4] + xform[5]*xform[5] + xform[ 6]*xform[ 6];
  dotYZ = xform[4]*xform[8] + xform[5]*xform[9] + xform[ 6]*xform[10];

  dotZZ = xform[8]*xform[8] + xform[9]*xform[9] + xform[10]*xform[10];

  if (sqrt(dotXX) < ZERO) {
    if (outLevel > 0) 
      printf(" EGADS Error: No Length on Transform  (EG_makeTransform)!\n");
    return EGADS_DEGEN;
  }
  if ((fabs((dotXX-dotYY)/dotXX) > ZERO) ||
      (fabs((dotXX-dotZZ)/dotXX) > ZERO)) {
    if (outLevel > 0) 
      printf(" EGADS Error: Skew Scaling in Transform  (EG_makeTransform)!\n");
    return EGADS_BADSCALE;
  }
  if ((fabs(dotXY/dotXX) > ZERO) || (fabs(dotXZ/dotXX) > ZERO) ||
      (fabs(dotYZ/dotXX) > ZERO)) {
    if (outLevel > 0) 
      printf(" EGADS Error: Transform not Orthogonal (EG_makeTransform)!\n");
    return EGADS_NOTORTHO;
  }

  reals = (double *) EG_alloc(12*sizeof(double));
  if (reals == NULL) {
      if (outLevel > 0) 
        printf(" EGADS Error: Malloc on transform (EG_makeTransform)!\n");
    return EGADS_MALLOC;
  }
  stat = EG_makeObject(context, &object);
  if (stat != EGADS_SUCCESS) {
    EG_free(reals);
    return stat;
  }
  
  object->oclass = TRANSFORM;
  object->blind  = reals;
  for (i = 0; i < 12; i++) reals[i] = xform[i];
  
  *oform = object;

  return EGADS_SUCCESS;
}


int
EG_getTransformation(const egObject *oform, double *xform)
{
  int    i;
  double *reals;

  for (i = 1; i < 12; i++) xform[i] = 0.0;
  xform[0] = xform[5] = xform[10]   = 1.0;
  if (oform == NULL)               return EGADS_NULLOBJ;
  if (oform->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (oform->oclass != TRANSFORM)  return EGADS_NOTXFORM;
  reals = (double *) oform->blind;
  if (reals == NULL)               return EGADS_NOTFOUND;
  
  for (i = 0; i < 12; i++) xform[i] = reals[i];
  
  return EGADS_SUCCESS;
}


int
EG_getContext(egObject *object, egObject **context)
{
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  
  *context = EG_context(object);
  return EGADS_SUCCESS;
}


int
EG_getInfo(const egObject *object, int *oclass, int *mtype, egObject **top,
           egObject **prev, egObject **next)
{
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;

  *oclass = object->oclass;
  *mtype  = object->mtype;
  *top    = object->topObj;
  *prev   = object->prev;
  *next   = object->next;

  return EGADS_SUCCESS;
}


int
EG_copyObject(const egObject *object, /*@null@*/ void *ptr, egObject **copy)
{
  int          i, j, k, stat, outLevel, npts, ntri, nn, *inode;
  double       *xform = NULL;
  const int    *ptype, *pindex, *tris, *tric;
  const double *xyz, *prms;
  egObject     *ocontext, *xcontext, *oform, *sobj, *obj = NULL;
  egTessel     *btess, *ctess;

  *copy = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == NIL)        return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == TRANSFORM)  return EGADS_NOTXFORM;
  outLevel = EG_outLevel(object);

  /* special tessellation object treatment */
  if (object->oclass == TESSELLATION) {
    
    btess = (egTessel *) object->blind;
    if (btess == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Blind Object (EG_copyObject)!\n");
      return EGADS_NOTFOUND;
    }
    sobj = btess->src;
    if (sobj == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Source Object (EG_copyObject)!\n");
      return EGADS_NULLOBJ;
    }
    if (sobj->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Source Not an Object (EG_copyObject)!\n");
      return EGADS_NOTOBJ;
    }
    if (sobj->oclass != BODY) {
      if (outLevel > 0)
        printf(" EGADS Error: Source Not Body (EG_copyObject)!\n");
      return EGADS_NOTBODY;
    }
    if (btess->done == 0) return EGADS_TESSTATE;

    xform = (double *) ptr;
    if (xform != NULL) {
      if (btess->globals == NULL) {
        stat = EG_computeTessMap(btess, outLevel);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_computeTessMap = %d (EG_copyObject)!\n",
                   stat);
          return stat;
        }
      }
    }
    
    stat = EG_initTessBody(sobj, copy);
    if ((stat != EGADS_SUCCESS) || (*copy == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_initTessBody = %d (EG_copyObject)!\n", stat);
      return stat;
    }
    for (i = 1; i <= btess->nEdge; i++) {
      stat = EG_getTessEdge(object, i, &npts, &xyz, &prms);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTessEdge %d = %d (EG_copyObject)!\n",
                 i, stat);
        EG_deleteObject(*copy);
        *copy = NULL;
        return stat;
      }
      stat = EG_setTessEdge(*copy, i, npts, xyz, prms);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_setTessEdge %d = %d (EG_copyObject)!\n",
                 i, stat);
        EG_deleteObject(*copy);
        *copy = NULL;
        return stat;
      }
    }
    for (i = 1; i <= btess->nFace; i++) {
      stat = EG_getTessFace(object, i, &npts, &xyz, &prms, &ptype, &pindex,
                            &ntri, &tris, &tric);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTessFace %d = %d (EG_copyObject)!\n",
                 i, stat);
        EG_deleteObject(*copy);
        *copy = NULL;
        return stat;
      }
      stat = EG_setTessFace(*copy, i, npts, xyz, prms, ntri, tris);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_setTessFace %d = %d (EG_copyObject)!\n",
                 i, stat);
        EG_deleteObject(*copy);
        *copy = NULL;
        return stat;
      }
    }
    ctess = (egTessel *) (*copy)->blind;
    ctess->done = 1;
    if ((btess->done != 2) && (xform == NULL)) return EGADS_SUCCESS;

    /* deal with displaced tessellation */
    if (ctess->globals == NULL) {
      stat = EG_computeTessMap(ctess, outLevel);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_computeTessMap Copy = %d (EG_copyObject)!\n",
                 stat);
        EG_deleteObject(*copy);
        *copy = NULL;
        return stat;
      }
    }
    if (btess->done == 2) ctess->done = 2;
    if (xform == NULL) return EGADS_SUCCESS;
    
    for (nn = i = 0; i < ctess->nEdge; i++) {
      if (ctess->tess1d[i].global == NULL) continue;
      if (nn < ctess->tess1d[i].nodes[0]) nn = ctess->tess1d[i].nodes[0];
      if (nn < ctess->tess1d[i].nodes[1]) nn = ctess->tess1d[i].nodes[1];
      for (j = 0; j < ctess->tess1d[i].npts; j++) {
        k = ctess->tess1d[i].global[j] - 1;
        if (k < 0) continue;
        ctess->tess1d[i].xyz[3*j  ] += xform[3*k  ];
        ctess->tess1d[i].xyz[3*j+1] += xform[3*k+1];
        ctess->tess1d[i].xyz[3*j+2] += xform[3*k+2];
      }
    }
    for (i = 0; i < ctess->nFace; i++) {
      if (ctess->tess2d[i].global == NULL) continue;
      for (j = 0; j < ctess->tess2d[i].npts; j++) {
        k = ctess->tess2d[i].global[j] - 1;
        if (k < 0) continue;
        ctess->tess2d[i].xyz[3*j  ] += xform[3*k  ];
        ctess->tess2d[i].xyz[3*j+1] += xform[3*k+1];
        ctess->tess2d[i].xyz[3*j+2] += xform[3*k+2];
      }
    }

    /* correct Node coordinates */
    inode = (int *) EG_alloc(nn*sizeof(int));
    if (inode == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation of %d Node markers (EG_copyObject)!\n",
               nn);
      EG_deleteObject(*copy);
      *copy = NULL;
      return EGADS_MALLOC;
    }
    for (j = 0; j < nn; j++) inode[j] = 0;
    for (i = 0; i < ctess->nEdge; i++) {
      if (ctess->tess1d[i].global == NULL) continue;
      j = ctess->tess1d[i].nodes[0] - 1;
      k = ctess->tess1d[i].global[0] - 1;
      if ((k >= 0) && (inode[j] == 0)) {
        ctess->xyzs[3*j  ] += xform[3*k  ];
        ctess->xyzs[3*j+1] += xform[3*k+1];
        ctess->xyzs[3*j+2] += xform[3*k+2];
        inode[j]++;
      }
      j = ctess->tess1d[i].nodes[1] - 1;
      k = ctess->tess1d[i].global[ctess->tess1d[i].npts-1] - 1;
      if ((k >= 0) && (inode[j] == 0)) {
        ctess->xyzs[3*j  ] += xform[3*k  ];
        ctess->xyzs[3*j+1] += xform[3*k+1];
        ctess->xyzs[3*j+2] += xform[3*k+2];
        inode[j]++;
      }
    }
    ctess->done = 2;
    EG_free(inode);
    return EGADS_SUCCESS;
  }
  
  /* get transform object (if any) */
  oform = (egObject *) ptr;
  if (oform != NULL) {
    if (oform->magicnumber != MAGIC) {
      if (outLevel > 0) 
        printf(" EGADS Error: 2nd argument not an EGO (EG_copyObject)!\n");
      return EGADS_NOTOBJ;
    }
    if (oform->oclass != TRANSFORM) {
      if (outLevel > 0) 
        printf(" EGADS Error: 2nd argument not an XForm (EG_copyObject)!\n");
      return EGADS_NOTXFORM;
    }
    ocontext = EG_context(object);
    xcontext = EG_context(oform);
    if (xcontext != ocontext) {
      if (outLevel > 0) 
        printf(" EGADS Error: Context mismatch (EG_copyObject)!\n");
      return EGADS_MIXCNTX;
    }
    xform = (double *) oform->blind;
  }
  
  if (object->oclass == PCURVE) {
  
    /* this does not make sense -- 2D! */
    if (outLevel > 0) 
      printf(" EGADS Error: PCurve is 2D (EG_copyObject)!\n");
    stat = EGADS_CONSTERR;

  } else if (object->oclass <= SURFACE) {
  
    stat = EG_copyGeometry(object, xform, &obj);
  
  } else {
  
    stat = EG_copyTopology(object, xform, &obj);

  }

  if (obj != NULL) {
    stat  = EG_attributeXDup(object, xform, obj);
    *copy = obj;
  }
  return stat;
}


int
EG_flipObject(const egObject *object, egObject **copy)
{
  int      stat;
  egObject *obj = NULL;

  *copy = NULL;
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == NIL)        return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;
  if (object->oclass == CONTXT)     return EGADS_NOTCNTX;
  if (object->oclass == TRANSFORM)  return EGADS_NOTXFORM;
  
  if (object->oclass == TESSELLATION) {
  
    /* what do we do here? */
    stat = EGADS_NOTTESS;
    
  } else if (object->oclass <= SURFACE) {
  
    stat = EG_flipGeometry(object, &obj);
  
  } else {
  
    stat = EG_flipTopology(object, &obj);

  }

  if (obj != NULL) {
    stat  = EG_attributeDup(object, obj);
    *copy = obj;
  }
  return stat;
}


int
EG_close(egObject *context)
{
  int      outLevel, cnt, ref, total, stat;
  egObject *obj, *next, *last;
  egCntxt  *cntx;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                  return EGADS_NODATA;
  outLevel = cntx->outLevel;

  /* count all active objects */
  
  cnt = ref = 0;
  obj = context->next;
  while (obj != NULL) {
    if (obj->magicnumber != MAGIC) {
      printf(" EGADS Info: Found BAD Object in cleanup (EG_close)!\n");
      printf("             Class = %d\n", obj->oclass);
      return EGADS_NOTFOUND;
    }
    if (obj->oclass == REFERENCE) {
      ref++;
    } else {
      if (outLevel > 2)
        printf(" EGADS Info: Object oclass = %d, mtype = %d Found!\n",
               obj->oclass, obj->mtype);
      cnt++;
    }
    obj = obj->next;
  }
  total = ref+cnt;
  obj   = cntx->pool;
  while (obj != NULL) {
    next = obj->next;
    obj = next;
    total++;
  }
  if (outLevel > 0)
    printf(" EGADS Info: %d Objects, %d Reference in Use (of %d) at Close!\n",
           cnt, ref, total);

  /* delete unattached geometry and topology objects */
  
  EG_deleteObject(context);
  
  /* delete tessellation objects */

  obj  = context->next;
  last = NULL;
  while (obj != NULL) {
    next = obj->next;
    if (obj->oclass == TESSELLATION)
      if (EG_deleteObject(obj) == EGADS_SUCCESS) {
        obj = last;
        if (obj == NULL) {
          next = context->next;
        } else {
          next = obj->next;
        }
      }
    last = obj;
    obj  = next;
  }

  /* delete all models */

  obj  = context->next;
  last = NULL;
  while (obj != NULL) {
    next = obj->next;
    if (obj->oclass == MODEL)
      if (EG_deleteObject(obj) == EGADS_SUCCESS) {
        obj = last;
        if (obj == NULL) {
          next = context->next;
        } else {
          next = obj->next;
        }
      }
    last = obj;
    obj  = next;
  }
  
  /* delete all bodies that are left */
  
  obj  = context->next;
  last = NULL;
  while (obj != NULL) {
    next = obj->next;
    if (obj->oclass == BODY)
      if (EG_deleteObject(obj) == EGADS_SUCCESS) {
        obj = last;
        if (obj == NULL) {
          next = context->next;
        } else {
          next = obj->next;
        }
      }
    last = obj;
    obj  = next;
  }
  
  /* dereference until nothing is left (which should be the case) */

  do {
    cnt = 0;
    obj = context->next;
    while (obj != NULL) {
      next = obj->next;
      if (obj->oclass != REFERENCE) {
        stat = EG_dereferenceTopObj(obj, NULL);
        if (stat == EGADS_SUCCESS) {
          cnt++;
          break;
        }
      }
      obj = next;
    }
    if (cnt != 0) continue;
    obj = context->next;
    while (obj != NULL) {
      next = obj->next;
      if (obj->oclass != REFERENCE) {
        stat = EG_dereferenceObject(obj, NULL);
        if (stat == EGADS_SUCCESS) {
          cnt++;
          break;
        }
      }
      obj = next;
    }    
  } while (cnt != 0);

  ref = cnt = 0;
  obj = context->next;
  while (obj != NULL) {
    if (cnt == 0)
      if (outLevel > 1)
        printf(" EGADS Info: Undeleted Object(s) in cleanup (EG_close):\n");
    if (obj->oclass == REFERENCE) {
      ref++;
    } else {
      if (outLevel > 1)
        printf("             %d: Class = %d, Type = %d\n", 
               cnt, obj->oclass, obj->mtype);
    }
    obj = obj->next;
    cnt++;
  }
  if (outLevel > 1)
    if ((cnt != 0) && (ref != 0))
      printf("             In Addition to %d Refereces\n", ref);  

  /* clean up the pool */
  
  obj = cntx->pool;
  while (obj != NULL) {
    if (obj->magicnumber != MAGIC) {
      printf(" EGADS Info: Found BAD Object in Cleanup (EG_close)!\n");
      printf("             Class = %d\n", obj->oclass);
      break;
    }
    next = obj->next;
    EG_free(obj);
    obj = next;
  }
  EG_attributeDel(context, NULL);
  EG_free(context);
  EG_free(cntx);
    
  return EGADS_SUCCESS;
}


int
EG_isSame(const egObject *obj1, const egObject *obj2)
{
  int      i, stat, oclass1, mtype1, oclass2, mtype2, nc = 0;
  int      *info1, *info2, *senses;
  double   limits[4], xyz1[3], xyz2[3], tol, toler, *rv1, *rv2;
  egObject *geom1, *geom2, *ref1, *ref2, **children;
  
  if (obj1 == obj2)                return EGADS_SUCCESS;
  if (obj1 == NULL)                return EGADS_NULLOBJ;
  if (obj1->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if ((obj1->oclass != NODE)  &&
      (obj1->oclass != EDGE)  &&
      (obj1->oclass != CURVE) &&
      (obj1->oclass != FACE)  &&
      (obj1->oclass != SURFACE))   return EGADS_NOTGEOM;
  if (obj1->blind == NULL)         return EGADS_NODATA;
  if (obj2 == NULL)                return EGADS_NULLOBJ;
  if (obj2->magicnumber != MAGIC)  return EGADS_NOTOBJ;
  if (obj1->oclass == NODE)
    if (obj2->oclass != NODE)      return EGADS_GEOMERR;
  if ((obj1->oclass == EDGE) ||
      (obj1->oclass == CURVE))
    if ((obj2->oclass != EDGE) &&
        (obj2->oclass != CURVE))   return EGADS_GEOMERR;
  if ((obj1->oclass == FACE) ||
      (obj1->oclass == SURFACE))
    if ((obj2->oclass != FACE) &&
        (obj2->oclass != SURFACE)) return EGADS_GEOMERR;
  if (obj2->blind == NULL)         return EGADS_NODATA;
  
  /* special Node checking */
  if (obj1->oclass == NODE) {
    stat = EG_getTolerance(obj1, &tol);
    if (stat != EGADS_SUCCESS) return stat;
    stat = EG_getTolerance(obj2, &toler);
    if (stat != EGADS_SUCCESS) return stat;
    if (toler > tol) toler = tol;
    stat = EG_getTopology(obj1, &geom1, &oclass1, &mtype1, xyz1,
                          &nc, &children, &senses);
    if (stat != EGADS_SUCCESS) return stat;
    stat = EG_getTopology(obj2, &geom2, &oclass2, &mtype2, xyz2,
                          &nc, &children, &senses);
    if (stat != EGADS_SUCCESS) return stat;
    if (sqrt((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) +
             (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) +
             (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2])) > toler) return EGADS_OUTSIDE;
    return EGADS_SUCCESS;
  }

  /* get the geometric objects */
  geom1  = (egObject *) obj1;
  mtype1 = 0;
  if ((obj1->oclass == EDGE) || (obj1->oclass == FACE)) {
    stat = EG_getTopology(obj1, &geom1, &oclass1, &mtype1, limits,
                          &nc, &children, &senses);
    if (stat != EGADS_SUCCESS)       return stat;
    if (geom1 == NULL)               return EGADS_NULLOBJ;
    if (geom1->magicnumber != MAGIC) return EGADS_NOTOBJ;
  }
  geom2  = (egObject *) obj1;
  mtype2 = 0;
  if ((obj2->oclass == EDGE) || (obj2->oclass == FACE)) {
    stat = EG_getTopology(obj2, &geom2, &oclass2, &mtype2, limits,
                          &nc, &children, &senses);
    if (stat != EGADS_SUCCESS)       return stat;
    if (geom2 == NULL)               return EGADS_NULLOBJ;
    if (geom2->magicnumber != MAGIC) return EGADS_NOTOBJ;
  }
  if (obj1->oclass == EDGE) {
    if ((mtype1 == DEGENERATE) && (mtype2 == DEGENERATE)) return EGADS_SUCCESS;
    if ((mtype1 == DEGENERATE) || (mtype2 == DEGENERATE)) return EGADS_OUTSIDE;
  }
  stat = EG_getGeometry(geom1, &oclass1, &mtype1, &ref1, &info1, &rv1);
  if (stat != EGADS_SUCCESS) return stat;
  stat = EG_getGeometry(geom2, &oclass2, &mtype2, &ref2, &info2, &rv2);
  if (stat != EGADS_SUCCESS) {
    if (info1 != NULL) EG_free(info1);
    EG_free(rv1);
    return stat;
  }
  
  /* are we the same type? */
  if (mtype1 != mtype2) {
    if (info1 != NULL) EG_free(info1);
    EG_free(rv1);
    if (info2 != NULL) EG_free(info2);
    EG_free(rv2);
    return EGADS_OUTSIDE;
  }
  
  if (mtype1 == TRIMMED) {
    /* ignore the trimming -- look at reference */
    EG_free(rv1);
    EG_free(rv2);
    return EG_isSame(ref1, ref2);
  }

  if (oclass1 == CURVE) {
    if (mtype1 == LINE) {
      nc = 6;
    } else if (mtype1 == CIRCLE) {
      nc = 10;
    } else if (mtype1 == ELLIPSE) {
      nc = 11;
    } else if (mtype1 == PARABOLA) {
      nc = 10;
    } else if (mtype1 == HYPERBOLA) {
      nc = 11;
    } else if (mtype1 == OFFSET) {
      for (stat = i = 0; i < 4; i++)
        if (rv1[i] != rv2[i]) {
          stat++;
          break;
        }
      EG_free(rv1);
      EG_free(rv2);
      if (stat != 0) return EGADS_OUTSIDE;
      return EG_isSame(ref1, ref2);
    } else if (mtype1 == BEZIER) {
      if ((info1[0] != info2[0]) || (info1[1] != info2[1]) ||
          (info1[2] != info2[2])) {
        EG_free(info1);
        EG_free(rv1);
        EG_free(info2);
        EG_free(rv2);
        return EGADS_OUTSIDE;
      }
      nc = info1[2]*3;
      if ((info1[0]&2) != 0) nc += info1[2];
    } else if (mtype1 == BSPLINE) {
      if ((info1[0] != info2[0]) || (info1[1] != info2[1]) ||
          (info1[2] != info2[2]) || (info1[3] != info2[3])) {
        EG_free(info1);
        EG_free(rv1);
        EG_free(info2);
        EG_free(rv2);
        return EGADS_OUTSIDE;
      }
      nc = info1[3] + info1[2]*3;
      if ((info1[0]&2) != 0) nc += info1[2];
    }
  } else {
    if (mtype1 == PLANE) {
      nc = 9;
    } else if (mtype1 == SPHERICAL) {
      nc = 10;
    } else if (mtype1 == CONICAL) {
      nc = 14;
    } else if (mtype1 == CYLINDRICAL) {
      nc = 13;
    } else if (mtype1 == TOROIDAL) {
      nc = 14;
    } else if (mtype1 == OFFSET) {
      stat = 0;
      if (rv1[0] != rv2[0]) stat++;
      EG_free(rv1);
      EG_free(rv2);
      if (stat != 0) return EGADS_OUTSIDE;
      return EG_isSame(ref1, ref2);
    } else if (mtype1 == EXTRUSION) {
      /* ref geometry is curve! */
      for (stat = i = 0; i < 3; i++)
        if (rv1[i] != rv2[i]) {
          stat++;
          break;
        }
      EG_free(rv1);
      EG_free(rv2);
      if (stat != 0) return EGADS_OUTSIDE;
      return EG_isSame(ref1, ref2);
    } else if (mtype1 == REVOLUTION) {
      /* ref geometry is curve! */
      for (stat = i = 0; i < 6; i++)
        if (rv1[i] != rv2[i]) {
          stat++;
          break;
        }
      EG_free(rv1);
      EG_free(rv2);
      if (stat != 0) return EGADS_OUTSIDE;
      return EG_isSame(ref1, ref2);
    } else if (mtype1 == BEZIER) {
      if ((info1[0] != info2[0]) || (info1[1] != info2[1]) ||
          (info1[2] != info2[2]) || (info1[3] != info2[3]) ||
          (info1[4] != info2[4])) {
        EG_free(info1);
        EG_free(rv1);
        EG_free(info2);
        EG_free(rv2);
        return EGADS_OUTSIDE;
      }
      nc = info1[2]*info1[4]*3;
      if ((info1[0]&2) != 0) nc += info1[2]*info1[4];
    } else if (mtype1 == BSPLINE) {
      if ((info1[0] != info2[0]) || (info1[1] != info2[1]) ||
          (info1[2] != info2[2]) || (info1[3] != info2[3]) ||
          (info1[4] != info2[4]) || (info1[5] != info2[5]) ||
          (info1[6] != info2[6])) {
        EG_free(info1);
        EG_free(rv1);
        EG_free(info2);
        EG_free(rv2);
        return EGADS_OUTSIDE;
      }
      nc = info1[3] + info1[6] + info1[2]*info1[5]*3;
      if ((info1[0]&2) != 0) nc += info1[2]*info1[5];

    }
  }
  if (info1 != NULL) EG_free(info1);
  if (info2 != NULL) EG_free(info2);
  
  for (stat = i = 0; i < nc; i++)
    if (rv1[i] != rv2[i]) {
      stat++;
      break;
    }
  
  EG_free(rv1);
  EG_free(rv2);
  if (stat != 0) return EGADS_OUTSIDE;
  
  return EGADS_SUCCESS;
}
