mkdir -p Application
./clean.sh
cd libFiles
rm -f *.a
./buildLinux.sh
cd ..
cd FRIEND_Engine
./buildLinux.sh
cd ..
cd PlugIn
cd connectivityPlugIn
./buildLinux.sh
cd ..
cd MotorPlugIn
./buildLinux.sh
cd ..
cd ROIbasedPlugIn
./buildLinux.sh
cd ..
cd svmPlugIn
./buildLinux.sh
cd ..
cd ..
