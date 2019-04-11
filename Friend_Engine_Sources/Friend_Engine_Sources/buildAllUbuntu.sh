mkdir -p Application
./clean.sh
cd libFiles
rm -f *.a
chmod +x buildUbuntu.sh
./buildUbuntu.sh

echo compiling engine Application
cd ..
cd FRIEND_Engine
chmod +x buildUbuntu.sh
./buildUbuntu.sh

echo compiling connectivity PlugIn
cd ..
cd PlugIn
cd connectivityPlugIn
chmod +x buildLinux.sh
./buildLinux.sh

echo compiling Motor PlugIn
cd ..
cd MotorPlugIn
chmod +x buildLinux.sh
./buildLinux.sh

echo compiling Roi PlugIn
cd ..
cd ROIbasedPlugIn
chmod +x buildLinux.sh
./buildLinux.sh

echo compiling SVM PlugIn
cd ..
cd svmPlugIn
chmod +x buildLinux.sh
./buildLinux.sh
cd ..
cd ..
