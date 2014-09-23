sudo gcc -fPIC -g -c -w \
svmPlugIn.cpp ../vardb.cpp ../masks.cpp ../../Friend_Engine/filefuncs.cpp ../../Friend_Engine/parser.cpp ../../Friend_Engine/intervals.cpp ../../Friend_Engine/fslfuncs.cpp svmobj.cpp svmfuncs.cpp ../../Friend_Engine/socket.cxx ../../Friend_Engine/socket2.cpp ../../Friend_Engine/defs.cpp ../../libsvm/svm.cpp \
../../alglib/statistics.cpp ../../alglib/alglibinternal.cpp ../../alglib/ap.cpp ../../alglib/specialfunctions.cpp ../../alglib/linalg.cpp ../../alglib/alglibmisc.cpp \
-DUNIX -DLINUX -DDARWIN -DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I. \
-I.. \
-I../../alglib \
-I../../Friend_Engine \
-I../../libsvm \
-I$FSLDIR/src \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/libgd \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/src/libpng \
-I$FSLDIR/extras/src/libprob \
-I../../simpleini
