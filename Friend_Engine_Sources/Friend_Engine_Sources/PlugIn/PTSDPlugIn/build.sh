sudo g++ -fPIC -g -c \
PTSDPlugIn.cpp ../vardb.cpp ../masks.cpp ../../FRIEND_Engine/utils.cpp ../../Friend_Engine/filefuncs.cpp ../../Friend_Engine/parser.cpp ../../Friend_Engine/intervals.cpp ../../Friend_Engine/session.cpp ../../Friend_Engine/fslfuncs.cpp ../../Friend_Engine/socket.cxx ../../Friend_Engine/socket2.cpp ../../Friend_Engine/defs.cpp \
../../alglib/statistics.cpp ../../alglib/alglibinternal.cpp ../../alglib/ap.cpp ../../alglib/specialfunctions.cpp ../../alglib/linalg.cpp ../../alglib/alglibmisc.cpp \
-DUNIX -DDARWIN -DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I. \
-I.. \
-I../../alglib \
-I../../Friend_Engine \
-I$FSLDIR/src \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/libgd \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/src/libpng \
-I$FSLDIR/extras/src/libprob \
-I../../simpleini


g++ -dynamiclib -o ../../Application/libPTSD.dylib PTSDPlugIn.o vardb.o masks.o intervals.o utils.o fslfuncs.o filefuncs.o defs.o parser.o socket.o socket2.o statistics.o alglibinternal.o ap.o specialfunctions.o linalg.o alglibmisc.o session.o \
-L$FSLDIR/extras/lib \
-L$FSLDIR/lib \
-lfslio \
-lnewimage \
-lmiscmaths \
-lcprob \
-lnewmat \
-lniftiio \
-lz \
-lznz \
-lmiscplot \
-lgdc \
-lgd \
-lpng
