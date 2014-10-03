mkdir -p Application
./clean.sh
cd libFiles
rm -f *.a
chmod +x buildUbuntu.sh
./buildUbuntu.sh
cd ..
cd FRIEND_Engine
chmod +x buildUbuntu.sh
./buildUbuntu.sh
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
