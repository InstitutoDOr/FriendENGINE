g++ -w -fPIC -g -c \
ROIPlugIn.cpp ../vardb.cpp ../masks.cpp  ../../FRIEND_Engine/logObject.cpp ../../FRIEND_Engine/utils.cpp ../../FRIEND_Engine/filefuncs.cpp ../../FRIEND_Engine/parser.cpp ../../FRIEND_Engine/intervals.cpp ../../FRIEND_Engine/fslfuncs.cpp ../../FRIEND_Engine/socket.cxx ../../FRIEND_Engine/socket2.cpp ../../FRIEND_Engine/defs.cpp \
../../alglib/statistics.cpp ../../alglib/alglibinternal.cpp ../../alglib/ap.cpp ../../alglib/specialfunctions.cpp ../../alglib/linalg.cpp ../../alglib/alglibmisc.cpp \
-DUNIX -DLINUX -DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I. \
-I.. \
-I../../alglib \
-I../../FRIEND_Engine \
-I$FSLDIR/src \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/libgd \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/libpng \
-I$FSLDIR/extras/src/libprob \
-I../../simpleini

g++ -shared -o ../../Application/libROI.so ROIPlugIn.o vardb.o masks.o logObject.o intervals.o utils.o fslfuncs.o filefuncs.o defs.o parser.o socket.o socket2.o statistics.o alglibinternal.o ap.o specialfunctions.o linalg.o alglibmisc.o  \
-L../../libFiles \
-L$FSLDIR/extras/lib \
-L$FSLDIR/lib \
-lm \
-lfslio \
-lnewimage \
-lmiscmaths \
-lcprob \
-lprob \
-lnewmat \
-lniftiio \
-lz \
-lznz \
-lmiscplot \
-lgdc \
-lgd \
-lpng

