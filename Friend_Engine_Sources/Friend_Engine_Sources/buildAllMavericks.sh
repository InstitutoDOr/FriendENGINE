mkdir -p Application
./clean.sh
cd libFiles
rm -f *.a
chmod +x buildMavericks.sh
./buildMavericks.sh

echo compiling engine Application
cd ..
cd FRIEND_Engine
chmod +x buildMavericks.sh
./buildMavericks.sh

echo compiling connectivity PlugIn
cd ..
cd PlugIn
cd connectivityPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh

echo compiling Motor PlugIn
cd ..
cd MotorPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh

echo compiling Roi PlugIn
cd ..
cd ROIbasedPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh
echo compiling SVM PlugIn
cd ..
cd svmPlugIn
chmod +x buildMavericks.sh
./buildMavericks.sh
cd ..
cd ..
