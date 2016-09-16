sudo g++ -fPIC -g -c \
connectivityPlugIn.cpp ../vardb.cpp ../masks.cpp ../../FRIEND_Engine/parser.cpp ../../FRIEND_Engine/utils.cpp ../../FRIEND_Engine/filefuncs.cpp ../../FRIEND_Engine/intervals.cpp ../../FRIEND_Engine/fslfuncs.cpp ../../FRIEND_Engine/defs.cpp ../../FRIEND_Engine/socket.cxx ../../FRIEND_Engine/socket2.cpp \
../../alglib/statistics.cpp ../../alglib/alglibinternal.cpp ../../alglib/ap.cpp ../../alglib/specialfunctions.cpp ../../alglib/linalg.cpp ../../alglib/alglibmisc.cpp \
-DUNIX -DDARWIN -DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I. \
-I.. \
-I../../alglib \
-I../../FRIEND_Engine \
-I$FSLDIR/src \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/libgd \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/src/libpng \
-I$FSLDIR/extras/src/libprob \
-I../../simpleini

g++ -dynamiclib -o ../../Application/libconnectivity.dylib connectivityPlugIn.o vardb.o masks.o intervals.o socket.o socket2.o fslfuncs.o utils.o filefuncs.o defs.o parser.o statistics.o alglibinternal.o ap.o specialfunctions.o linalg.o alglibmisc.o \
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

