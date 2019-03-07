echo recreating newmat library
export BUILDSTRING="\"5.0.10\"" 
sudo g++ -fPIC -g -c \
$FSLDIR/extras/src/newmat/bandmat.cpp $FSLDIR/extras/src/newmat/cholesky.cpp $FSLDIR/extras/src/newmat/evalue.cpp $FSLDIR/extras/src/newmat/fft.cpp $FSLDIR/extras/src/newmat/hholder.cpp $FSLDIR/extras/src/newmat/jacobi.cpp $FSLDIR/extras/src/newmat/myexcept.cpp $FSLDIR/extras/src/newmat/newmat1.cpp $FSLDIR/extras/src/newmat/newmat2.cpp $FSLDIR/extras/src/newmat/newmat3.cpp $FSLDIR/extras/src/newmat/newmat4.cpp $FSLDIR/extras/src/newmat/newmat5.cpp $FSLDIR/extras/src/newmat/newmat6.cpp $FSLDIR/extras/src/newmat/newmat7.cpp $FSLDIR/extras/src/newmat/newmat8.cpp $FSLDIR/extras/src/newmat/newmat9.cpp $FSLDIR/extras/src/newmat/newmatex.cpp $FSLDIR/extras/src/newmat/newmatnl.cpp $FSLDIR/extras/src/newmat/newmatrm.cpp $FSLDIR/extras/src/newmat/solution.cpp $FSLDIR/extras/src/newmat/sort.cpp $FSLDIR/extras/src/newmat/submat.cpp $FSLDIR/extras/src/newmat/svd.cpp $FSLDIR/extras/src/newmat/newfft.cpp \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/extras/src/newmat 

sudo ar rcs libnewmat.a bandmat.o cholesky.o evalue.o fft.o hholder.o jacobi.o myexcept.o newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o newmat6.o newmat7.o newmat8.o newmat9.o newmatex.o newmatnl.o newmatrm.o solution.o sort.o submat.o svd.o newfft.o


echo recreating niftiio
sudo gcc -fPIC -g -c $FSLDIR/src/niftiio/nifti1_io.c -I$FSLDIR/src/niftiio \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include

sudo ar rcs libniftiio.a nifti1_io.o



echo recreating newimage library
sudo g++ -fPIC -g -c \
$FSLDIR/src/newimage/lazy.cc $FSLDIR/src/newimage/newimage.cc $FSLDIR/src/newimage/generalio.cc $FSLDIR/src/newimage/newimagefns.cc $FSLDIR/src/newimage/complexvolume.cc $FSLDIR/src/newimage/imfft.cc $FSLDIR/src/newimage/costfns.cc \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/include/newmat

sudo ar rcs libnewimage.a lazy.o newimage.o generalio.o newimagefns.o complexvolume.o imfft.o costfns.o



echo recreating miscmaths library
sudo g++ -fPIC -g -c $FSLDIR/src/miscmaths/miscmaths.cc $FSLDIR/src/miscmaths/optimise.cc $FSLDIR/src/miscmaths/miscprob.cc $FSLDIR/src/miscmaths/kernel.cc $FSLDIR/src/miscmaths/histogram.cc $FSLDIR/src/miscmaths/base2z.cc $FSLDIR/src/miscmaths/t2z.cc $FSLDIR/src/miscmaths/f2z.cc $FSLDIR/src/miscmaths/minimize.cc $FSLDIR/src/miscmaths/cspline.cc $FSLDIR/src/miscmaths/sparse_matrix.cc $FSLDIR/src/miscmaths/sparsefn.cc $FSLDIR/src/miscmaths/rungekutta.cc $FSLDIR/src/miscmaths/nonlin.cpp $FSLDIR/src/miscmaths/Simplex.cpp $FSLDIR/src/miscmaths/bfmatrix.cpp \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/src/newmat \
-I$FSLDIR/extras/src/zlib \
-I$FSLDIR/extras/include/boost \

sudo ar rcs libmiscmaths.a miscmaths.o optimise.o miscprob.o kernel.o histogram.o base2z.o t2z.o f2z.o minimize.o cspline.o sparse_matrix.o sparsefn.o rungekutta.o



echo recreating fslio library
sudo gcc -fPIC -g -c $FSLDIR/src/fslio/fslio.c \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB -DBUILDSTRING=$BUILDSTRING \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/include/newmat 

sudo ar rcs libfslio.a fslio.o

echo recreating znzlib library
sudo gcc -fPIC -g -c $FSLDIR/src/znzlib/znzlib.c \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/include/newmat 

sudo ar rcs libznz.a znzlib.o

echo recreating warpfns library
sudo gcc -fPIC -g -c $FSLDIR/src/warpfns/warpfns.cc $FSLDIR/src/warpfns/fnirt_file_reader.cpp $FSLDIR/src/warpfns/fnirt_file_writer.cpp $FSLDIR/src/warpfns/point_list.cpp \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/include/newmat \
-I$FSLDIR/extras/include/boost/

sudo ar rcs libwarpfns.a warpfns.o fnirt_file_reader.o fnirt_file_writer.o point_list.o

echo recreating meshclass library
sudo gcc -fPIC -g -c $FSLDIR/src/meshclass/point.cpp $FSLDIR/src/meshclass/mpoint.cpp $FSLDIR/src/meshclass/triangle.cpp $FSLDIR/src/meshclass/mesh.cpp $FSLDIR/src/meshclass/pt_special.cpp $FSLDIR/src/meshclass/profile.cpp \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/include/newmat \
-I$FSLDIR/extras/include/boost/

sudo ar rcs libmeshclass.a point.o mpoint.o triangle.o mesh.o pt_special.o profile.o 

echo recreating basisfield library
sudo gcc -fPIC -g -c $FSLDIR/src/basisfield/dctfield.cpp $FSLDIR/src/basisfield/splinefield.cpp $FSLDIR/src/basisfield/basisfield.cpp $FSLDIR/src/basisfield/splines.c \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/include/newmat \
-I$FSLDIR/extras/include/boost/

sudo ar rcs libbasisfield.a dctfield.o splinefield.o basisfield.o splines.o 

echo recreating utils library
sudo gcc -fPIC -g -c $FSLDIR/src/utils/matches.cc $FSLDIR/src/utils/functions.cc $FSLDIR/src/utils/usage.cc $FSLDIR/src/utils/check.cc $FSLDIR/src/utils/parse.cc $FSLDIR/src/utils/log.cc $FSLDIR/src/utils/time_tracer.cc \
-DEXPOSE_TREACHEROUS -DHAVE_LIBPNG -DHAVE_ZLIB -DBUILDSTRING=$BUILDSTRING \
-I$FSLDIR/src \
-I$FSLDIR/extras/include \
-I$FSLDIR/extras/src/libgdc \
-I$FSLDIR/extras/include/libprob \
-I$FSLDIR/extras/include/newmat \
-I$FSLDIR/extras/include/boost/

sudo ar rcs libutils.a matches.o functions.o usage.o check.o parse.o log.o time_tracer.o 
