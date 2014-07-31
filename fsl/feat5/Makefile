include $(FSLCONFDIR)/default.mk

PROJNAME = feat5

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} 
USRLDFLAGS =  -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG}

IOLIBS = -lfslio -lniftiio -lznz
LIBS = -lutils -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lprob -lnewmat $(IOLIBS) -lm -lz  -lgdc -lgd -lpng

XFILES = cutoffcalc feat_model tsplot prewhiten

RUNTCLS = Feat Renderhighres Featquery Glm

SCRIPTS = easythresh feat mccutup renderhighres featregapply featquery fslerrorreport invfeatreg updatefeatreg estnoise perfusion_subtract match_smoothing fsl_motion_outliers feat_gm_prepare

all: cutoffcalc feat_model tsplot prewhiten defaultfiles

defaultfiles:
	sed s/VVV/`grep "set fmri(version)" featlib.tcl | grep "\." | awk '{print $$3}'`/g doc/ORIGtop.html > doc/top.html
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	${CP} -rf default_flobs.flobs ${DESTDIR}/etc
	${CP} -f featconf doc/feat.tcl
	@if [ ! -d ${DESTDIR}/etc/fslconf ] ; then ${MKDIR} ${DESTDIR}/etc/fslconf ; ${CHMOD} g+w ${DESTDIR}/etc/fslconf ; fi
	${CP} -f featconf ${DESTDIR}/etc/fslconf/feat.tcl

featlib.o: featlib.cc featlib.h
	$(CXX) $(CXXFLAGS) -c -o featlib.o featlib.cc

tsplot: tsplot.cc featlib.o featlib.h
	$(CXX) $(CXXFLAGS) -o tsplot featlib.o tsplot.cc $(LDFLAGS) $(LIBS)

prewhiten: prewhiten.cc featlib.o featlib.h
	$(CXX) $(CXXFLAGS) -o prewhiten featlib.o prewhiten.cc $(LDFLAGS) $(LIBS)

cutoffcalc: cutoffcalc.cc 
	$(CXX) $(CXXFLAGS) -o cutoffcalc cutoffcalc.cc $(LDFLAGS) $(LIBS)

feat_model: feat_model.cc featlib.o
	$(CXX) $(CXXFLAGS) -o feat_model featlib.o feat_model.cc $(LDFLAGS) $(LIBS)

retroicor_evs: retroicor_evs.cc 
	$(CXX) $(CXXFLAGS) -o retroicor_evs retroicor_evs.cc $(LDFLAGS) $(LIBS)

