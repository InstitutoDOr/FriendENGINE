mkdir -p Application
./clean.sh
cd libFiles
rm -f *.a
echo chmod +x buildMavericks.sh
echo ./buildMavericks.sh
cd ..
cd FRIEND_Engine
chmod +x buildMavericks.sh
./buildMavericks.sh
cd ..
cd PlugIn
cd connectivityPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh
cd ..
cd MotorPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh
cd ..
cd ROIbasedPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh
cd ..
cd svmPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh
cd ..
cd ..
