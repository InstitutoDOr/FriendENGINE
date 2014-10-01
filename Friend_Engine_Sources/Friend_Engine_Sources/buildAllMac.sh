mkdir -p Application
./clean.sh
cd FRIEND_Engine
chmod +x build.sh
./build.sh
cd ..
cd PlugIn
cd connectivityPlugIn
chmod +x build.sh
./build.sh
cd ..
cd MotorPlugIn
chmod +x build.sh
./build.sh
cd ..
cd ROIbasedPlugIn
chmod +x build.sh
./build.sh
cd ..
cd svmPlugIn
chmod +x build.sh
./build.sh
cd ..
cd ..
