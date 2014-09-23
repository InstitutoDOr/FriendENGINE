include ${FSLCONFDIR}/default.mk

PROJNAME=fnirt

USRINCFLAGS=-I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS=-L${LIB_NEWMAT} 

FNIRT_OBJS=fnirt_costfunctions.o intensity_mappers.o matching_points.o fnirtfns.o fnirt.o
LIBS=-lwarpfns -lmeshclass -lbasisfield -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lutils -lnewmat -lm -lz

XFILES=fnirt

all: ${XFILES} config

config:
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	@if [ ! -d ${DESTDIR}/etc/flirtsch ] ; then ${MKDIR} ${DESTDIR}/etc/flirtsch ; ${CHMOD} g+w ${DESTDIR}/etc/flirtsch ; fi
	${CP} -rf fnirtcnf/* ${DESTDIR}/etc/flirtsch/.

fnirt: ${FNIRT_OBJS}
	${CXX} ${LDFLAGS} -o fnirt ${FNIRT_OBJS} ${LIBS}
