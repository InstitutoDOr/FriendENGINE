echo recreating newmat library
sudo g++ -fPIC -g -c \
$FSLDIR/extras/src/newmat/bandmat.cpp $FSLDIR/extras/src/newmat/cholesky.cpp $FSLDIR/extras/src/newmat/evalue.cpp $FSLDIR/extras/src/newmat/fft.cpp $FSLDIR/extras/src/newmat/hholder.cpp $FSLDIR/extras/src/newmat/jacobi.cpp $FSLDIR/extras/src/newmat/myexcept.cpp $FSLDIR/extras/src/newmat/newmat1.cpp $FSLDIR/extras/src/newmat/newmat2.cpp $FSLDIR/extras/src/newmat/newmat3.cpp $FSLDIR/extras/src/newmat/newmat4.cpp $FSLDIR/extras/src/newmat/newmat5.cpp $FSLDIR/extras/src/newmat/newmat6.cpp $FSLDIR/extras/src/newmat/newmat7.cpp $FSLDIR/extras/src/newmat/newmat8.cpp $FSLDIR/extras/src/newmat/newmat9.cpp $FSLDIR/extras/src/newmat/newmatex.cpp $FSLDIR/extras/src/newmat/newmatnl.cpp $FSLDIR/extras/src/newmat/newmatrm.cpp $FSLDIR/extras/src/newmat/solution.cpp $FSLDIR/extras/src/newmat/sort.cpp $FSLDIR/extras/src/newmat/submat.cpp $FSLDIR/extras/src/newmat/svd.cpp $FSLDIR/extras/src/newmat/newfft.cpp \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/newmat 


sudo ar rcs libnewmat.a bandmat.o cholesky.o evalue.o fft.o hholder.o jacobi.o myexcept.o newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o newmat6.o newmat7.o newmat8.o newmat9.o newmatex.o newmatnl.o newmatrm.o solution.o sort.o submat.o svd.o newfft.o


echo recreating niftiio library
sudo gcc -fPIC -g -c $FSLDIR/src/niftiio/nifti1_io.c -I$FSLDIR/src/niftiio \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/include

sudo ar rcs libniftiio.a nifti1_io.o



echo recreating newimage library
sudo g++ -fPIC -g -c \
$FSLDIR/src/newimage/lazy.cc $FSLDIR/src/newimage/newimage.cc $FSLDIR/src/newimage/generalio.cc $FSLDIR/src/newimage/newimagefns.cc $FSLDIR/src/newimage/complexvolume.cc $FSLDIR/src/newimage/imfft.cc $FSLDIR/src/newimage/costfns.cc \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/include/newmat

sudo ar rcs libnewimage.a lazy.o newimage.o generalio.o newimagefns.o complexvolume.o imfft.o costfns.o



echo recreating miscmaths library
sudo g++ -fPIC -g -c $FSLDIR/src/miscmaths/miscmaths.cc $FSLDIR/src/miscmaths/optimise.cc $FSLDIR/src/miscmaths/miscprob.cc $FSLDIR/src/miscmaths/kernel.cc $FSLDIR/src/miscmaths/histogram.cc $FSLDIR/src/miscmaths/base2z.cc $FSLDIR/src/miscmaths/t2z.cc $FSLDIR/src/miscmaths/f2z.cc $FSLDIR/src/miscmaths/minimize.cc $FSLDIR/src/miscmaths/cspline.cc $FSLDIR/src/miscmaths/sparse_matrix.cc $FSLDIR/src/miscmaths/sparsefn.cc $FSLDIR/src/miscmaths/rungekutta.cc \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/include/newmat 

sudo ar rcs libmiscmaths.a miscmaths.o optimise.o miscprob.o kernel.o histogram.o base2z.o t2z.o f2z.o minimize.o cspline.o sparse_matrix.o sparsefn.o rungekutta.o



echo recreating fslio library
sudo gcc -fPIC -g -c $FSLDIR/src/fslio/fslio.c \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/include/newmat 

sudo ar rcs libfslio.a fslio.o

echo recreating zlib library
sudo gcc -fPIC -g -c $FSLDIR/extras/src/zlib/adler32.c $FSLDIR/extras/src/zlib/compress.c \
 $FSLDIR/extras/src/zlib/crc32.c $FSLDIR/extras/src/zlib/deflate.c $FSLDIR/extras/src/zlib/gzclose.c \
 $FSLDIR/extras/src/zlib/gzlib.c $FSLDIR/extras/src/zlib/gzread.c \
 $FSLDIR/extras/src/zlib/gzwrite.c $FSLDIR/extras/src/zlib/infback.c $FSLDIR/extras/src/zlib/inffast.c \
 $FSLDIR/extras/src/zlib/inflate.c $FSLDIR/extras/src/zlib/inftrees.c $FSLDIR/extras/src/zlib/trees.c \
 $FSLDIR/extras/src/zlib/uncompr.c $FSLDIR/extras/src/zlib/zutil.c \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/include/newmat 

sudo ar rcs libz.a adler32.o compress.o crc32.o deflate.o gzclose.o gzlib.o gzread.o \
 gzwrite.o infback.o inffast.o inflate.o inftrees.o trees.o uncompr.o zutil.o


echo recreating znzlib library
sudo gcc -fPIC -g -c $FSLDIR/src/znzlib/znzlib.c \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/include/newmat \
-L. \
-lz -lm

sudo ar rcs libznz.a znzlib.o adler32.o compress.o crc32.o deflate.o gzclose.o gzlib.o gzread.o gzwrite.o infback.o inffast.o inflate.o inftrees.o trees.o uncompr.o zutil.o

echo recreating libgdc library
sudo gcc -w -fPIC -g -c $FSLDIR/extras/src/libgdc/gifencode.c \
 $FSLDIR/extras/src/libgdc/price_conv.c $FSLDIR/extras/src/libgdc/gdc.c \
 $FSLDIR/extras/src/libgdc/gdc_pie.c $FSLDIR/extras/src/libgdc/gdchart.c \
 $FSLDIR/extras/src/libgdc/array_alloc.c \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/include/newmat \
-lgd -lpng -lz -lm

sudo ar rcs libgdc.a gifencode.o price_conv.o gdc.o gdc_pie.o gdchart.o array_alloc.o


echo recreating miscplot library
sudo gcc -fPIC -g -c $FSLDIR/src/libvis/miscplot.cc \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/libprob \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/include/newmat \
-lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm -lgdc -lgd -lpng -lz

sudo ar rcs libmiscplot.a miscplot.o




echo recreating libprob library
sudo gcc -fPIC -g -c $FSLDIR/extras/src/cprob/bdtr.c $FSLDIR/extras/src/cprob/btdtr.c \
 $FSLDIR/extras/src/cprob/chdtr.c $FSLDIR/extras/src/cprob/drand.c \
 $FSLDIR/extras/src/cprob/expx2.c $FSLDIR/extras/src/cprob/fdtr.c \
 $FSLDIR/extras/src/cprob/gamma.c $FSLDIR/extras/src/cprob/gdtr.c \
 $FSLDIR/extras/src/cprob/igam.c $FSLDIR/extras/src/cprob/igami.c \
 $FSLDIR/extras/src/cprob/incbet.c $FSLDIR/extras/src/cprob/incbi.c \
 $FSLDIR/extras/src/cprob/mtherr.c $FSLDIR/extras/src/cprob/nbdtr.c \
 $FSLDIR/extras/src/cprob/ndtr.c $FSLDIR/extras/src/cprob/ndtri.c \
 $FSLDIR/extras/src/cprob/pdtr.c $FSLDIR/extras/src/cprob/stdtr.c \
 $FSLDIR/extras/src/cprob/unity.c $FSLDIR/extras/src/cprob/polevl.c \
 $FSLDIR/extras/src/cprob/const.c $FSLDIR/extras/src/cprob/xmath.c \
-I$FSLDIR/extras/src/cprob \

sudo ar rcs libprob.a bdtr.o btdtr.o chdtr.o drand.o expx2.o fdtr.o gamma.o gdtr.o \
igam.o igami.o incbet.o incbi.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o \
stdtr.o unity.o polevl.o const.o xmath.o

ln -s libprob.a libcprob.a


