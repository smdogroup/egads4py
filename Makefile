include Makefile.in
include EGADS_Common.mk

EGADS_SUBDIRS = src

EGADS_OBJS := $(addsuffix /*.o, ${EGADS_SUBDIRS})

default:
	echo "Building Real EGADS"
	@for subdir in $(EGADS_SUBDIRS) ; do \
	    echo "making $@ in $$subdir"; \
	    echo; (cd $$subdir && $(MAKE) EGADS_DIR=${EGADS_DIR}) || exit 1; \
	done
	${CXX} ${SO_LINK_FLAGS} ${EGADS_OBJS} ${EGADS_EXTERN_LIBS} -o ${EGADS_DIR}/lib/libegads.${SO_EXT} 

debug:
	echo "Building Real EGADS"
	@for subdir in $(EGADS_SUBDIRS) ; do \
	    echo "making $@ in $$subdir"; \
	    echo; (cd $$subdir && $(MAKE) debug EGADS_DIR=${EGADS_DIR}) || exit 1; \
	done
	${CXX} ${SO_LINK_FLAGS} ${EGADS_OBJS} ${EGADS_EXTERN_LIBS} -o ${EGADS_DIR}/lib/libegads.${SO_EXT}

interface:
	pip install -e .\[all\]

clean:
	${RM} lib/*.a lib/*.so
	@for subdir in $(EGADS_SUBDIRS) ; do \
	    echo "making $@ in $$subdir"; \
	    echo; (cd $$subdir && $(MAKE) clean EGADS_DIR=${EGADS_DIR}) || exit 1; \
	done
