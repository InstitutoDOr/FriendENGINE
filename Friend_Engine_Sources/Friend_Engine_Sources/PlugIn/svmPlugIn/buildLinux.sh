sudo g++ -fPIC -g -c \
svmPlugIn.cpp ../vardb.cpp ../masks.cpp ../../FRIEND_Engine/confusionmatrix.cpp ../../FRIEND_Engine/filefuncs.cpp ../../FRIEND_Engine/parser.cpp ../../FRIEND_Engine/intervals.cpp ../../FRIEND_Engine/fslfuncs.cpp svmobj.cpp svmfuncs.cpp ../../FRIEND_Engine/socket.cxx ../../FRIEND_Engine/socket2.cpp ../../FRIEND_Engine/defs.cpp ../../libsvm/svm.cpp \
../../alglib/statistics.cpp ../../alglib/alglibinternal.cpp ../../alglib/ap.cpp ../../alglib/specialfunctions.cpp ../../alglib/linalg.cpp ../../alglib/alglibmisc.cpp \
-DUNIX -DLINUX -DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I. \
-I.. \
-I../../alglib \
-I../../FRIEND_Engine \
-I../../libsvm \
-I$FSLDIR/src \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/libgd \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/libpng \
-I$FSLDIR/extras/src/libprob \
-I../../simpleini

g++ -shared -o ../../Application/libBrainDecoding.so svmPlugIn.o confusionmatrix.o vardb.o masks.o intervals.o fslfuncs.o svmobj.o svmfuncs.o defs.o svm.o parser.o statistics.o socket.o socket2.o alglibinternal.o filefuncs.o ap.o specialfunctions.o linalg.o alglibmisc.o \
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
