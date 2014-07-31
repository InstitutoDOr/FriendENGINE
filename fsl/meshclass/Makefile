# $Id: Makefile,v 1.8 2007/11/02 11:06:19 flitney Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = meshclass

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB}

LIBS=-lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lz

MESHOBJS=point.o mpoint.o triangle.o mesh.o pt_special.o profile.o 
DRAWOBJS=drawmesh.o mpoint.o triangle.o mesh.o point.o pt_special.o
INTEROBJS=selfintersection.o mpoint.o triangle.o mesh.o point.o pt_special.o

XFILES=drawmesh selfintersection

DBGFLAGS=-g

all: meshclass.a ${XFILES}

meshclass.a: ${MESHOBJS}
	${AR} -r libmeshclass.a ${MESHOBJS}

drawmesh:${DRAWOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DRAWOBJS} ${LIBS} 

selfintersection:${INTEROBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${INTEROBJS} ${LIBS} 



