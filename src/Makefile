# ============================================
#
# Make file for EGADS_DIR/src/
#
# ============================================

include ../Makefile.in
include ../EGADS_Common.mk

CXX_OBJS = egadsAttrs.o \
	egadsBlendSens.o \
	egadsFit.o \
	egadsHLevel.o \
	egadsMemory.o \
	egadsRobust.o \
	egadsSpline.o \
	egadsTess.o \
	egadsTopo.o \
	emp.o \
	prmGrid.o \
	egadsBase.o \
	egadsCopy.o \
	egadsGeom.o \
	egadsIO.o \
	egadsQuads.o \
	egadsSBO.o \
	egadsSplineFit.o \
	egadsTessInp.o \
	egadsTris.o \
	prmCfit.o \
	prmUV.o \

DIR=${EGADS_DIR}/src

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${EGADS_LIB} ${CXX_OBJS}

debug: EGADS_CC_FLAGS=${EGADS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
