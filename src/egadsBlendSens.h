#ifndef EGADSBLENDSENS_H
#define EGADSBLENDSENS_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Blend Derivative Functions Header
 *
 *      Copyright 2011-2016, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"


#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif

__ProtoExt__ int  EG_blend_init( int nSection, const ego *sections,
                                 /*@null@*/ const double *rc1,
                                 /*@null@*/ const double *rc1_dot,
                                 /*@null@*/ const double *rcN,
                                 /*@null@*/ const double *rcN_dot,
                                 int *nStrips, void **cache );
  
__ProtoExt__ int  EG_blend_pos( void *cache, int sIndex, ego **objs,
                                int *imax, double **ts );

__ProtoExt__ int  EG_blend_sens( void *cache, int sIndex,
                                 const double *xyzs, const double *xyzs_dot,
                                 /*@null@*/ const double *tbeg,
                                 /*@null@*/ const double *tbeg_dot,
                                 /*@null@*/ const double *tend,
                                 /*@null@*/ const double *tend_dot);
  
__ProtoExt__ int  EG_blend_seval( void *cache, int sIndex, const double *uv,
                                  double *xyz, double *vel );
  
__ProtoExt__ void EG_sens_free( /*@null@*/ /*@only@*/ void *cache );
  
#ifdef __cplusplus
}
#endif

#endif