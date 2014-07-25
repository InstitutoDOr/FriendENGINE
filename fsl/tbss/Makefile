include ${FSLCONFDIR}/default.mk

PROJNAME = tbss

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

XFILES = distancemap tbss_skeleton

SCRIPTS = tbss_1_preproc tbss_2_reg tbss_3_postreg tbss_4_prestats fsl_reg tbss_non_FA tbss_fill tbss_3_postreg_new tbss_non_FA_new tbss_deproject tbss_sym

all: ${XFILES} 

distancemap: distancemap.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ distancemap.cc ${LIBS} 

tbss_skeleton: tbss_skeleton.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ tbss_skeleton.cc ${LIBS} 

