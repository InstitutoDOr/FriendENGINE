cd Application
rm -f *.so
rm -f engine

cd ..
cd libFiles
rm -f *.o
rm -f *.a

cd ..
cd FRIEND_Engine
rm -f *.o
rm -f engine

cd ..
cd PlugIn
cd ROIbasedPlugIn

rm -f *.o
rm -f *.dylib
rm -f *.so

cd ..
cd connectivityPlugIn

rm -f *.o
rm -f *.dylib
rm -f *.so

cd ..
cd svmPlugIn

rm -f *.o
rm -f *.dylib
rm -f *.so 

cd ..
cd MotorPlugIn

rm -f *.o
rm -f *.dylib
rm -f *.so

cd ..
cd ..

cd TestData/SUBJ002
rm -f -R input
rm -f -R glm
rm -f -R log
rm -f -R preproc_RUN01
rm -f -R preproc_RUN02
rm -f -R preproc_RUN03
rm -f -R preproc_RUN04
rm -f -R copied_RUN01
rm -f -R copied_RUN02
rm -f -R copied_RUN03
rm -f -R copied_RUN04
rm -f study_params_RUN01.txt
rm -f study_params_RUN02.txt
rm -f study_params_RUN03.txt
rm -f study_params_RUN04.txt
rm -f MNI2RFI_RFI2.mat
rm -f RFI2MNI_RFI2.*

cd ..
rm -f hmat_spm_final_std.nii.gz
cd ..

rm -f Application/dcm2niiserver.ini
