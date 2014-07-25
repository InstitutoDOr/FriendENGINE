# $Id: Makefile,v 1.12 2007/06/14 14:50:20 duncan Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = bet2

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS=-lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

BETOBJS=bet2.o 
BETSURFOBJS=betsurf.o 

XFILES=bet2 betsurf

SCRIPTS = old_betall bet

RUNTCLS = Bet

all: ${XFILES}

bet2: ${BETOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BETOBJS} ${LIBS} 

betsurf: ${BETSURFOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${BETSURFOBJS} ${LIBS} 


