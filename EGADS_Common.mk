EGADS_LIB = ${EGADS_DIR}/lib/libegads.a

EGADS_INCLUDE = -I${EGADS_DIR}/include ${OPENCASCADE_INCLUDE} ${NETGEN_INCLUDE}

# Set the compiler flags for EGADS
EGADS_CC_FLAGS = ${EGADS_FLAGS} ${EGADS_INCLUDE}
EGADS_DEBUG_CC_FLAGS = ${EGADS_DEBUG_FLAGS} ${EGADS_INCLUDE}

# Set the compiler flags
EGADS_EXTERN_LIBS = ${OPENCASCADE_LIB_PATH} ${OPENCASCADE_LIBS}
EGADS_LD_FLAGS = ${EGADS_LD_CMD} ${EGADS_EXTERN_LIBS}

# This is the one rule that is used to compile all the source
%.o: %.cpp
	${CXX} ${EGADS_CC_FLAGS} -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.cpp successfully ---"
	@echo

# This is the one rule that is used to compile all the source
%.o: %.c
	${CC} ${EGADS_CC_FLAGS} -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.c successfully ---"
	@echo
