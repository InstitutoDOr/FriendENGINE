## Usage - FRIEND Engine

1. Download and installation

Users can download the Friend Engine software through:

* a zip file at [https://github.com/InstitutoDOr/FriendENGINE](https://github.com/InstitutoDOr/FriendENGINE).
* the command:

```git
git clone [https://github.com/InstitutoDOr/FriendENGINE.git](https://github.com/InstitutoDOr/FriendENGINE.gitE)
```

> Preferred option as the users can keep track of the changes made to the software through git commands.

After unpacking or cloning the FRIEND Engine distribution, the following (reduced for presentation issues) structure will be shown:

* [**Friend\_Engine\_Sources**](https://github.com/InstitutoDOr/FriendENGINE) **(aka:FES)**
* **Friend\_Engine\_Sources (aka:FES)**
* Alglib
* Application
* study\_params.txt
* FRIEND\_Engine
* build.sh
* buildLinux.sh
    ...
* libFiles
* liblinear
* libsvm
* PlugIn
* connectivityPlugIn
* MotorPlugIn
* ROIbasedPlugIn
* svmPlugIn
* ManualPlugIn.rtf
...
* vardb.cpp
* Simpleini
* TestData
* buildAllLinux.sh
* buildUbuntu.sh
* buildAllMac.sh
* clean.sh
* Windows
* Application
* FriendEngine
* libBrainDecoding
* libConnectivity
* libMotor
* libROI
* clean.bat
* Frontends
* **Unity/FingerTap**
* **Unity/Medieval**
* Matlab
* Python
* **fsl**
* FRIEND\_License.txt



2. Download example data

Example data of a motor study is provided with the distribution (FES/FES/TestData). Other datasets can be downloaded from the NITRC repository: [http://www.nitrc.org/projects/friend](http://www.nitrc.org/projects/friend) by registered users.


3. Study parameters file

Before trying to execute any study, you should take some time to check/modify study configuration parameters (subject directory, data directory, etc.) in Windows: [FES\Windows\Application\study\_params.txt](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Windows/Application/study_params.txt)

MAC + Linux: [FES\FES\Application\study\_params.txt](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Friend_Engine_Sources/Application/study_params.txt)

4. MAC+Linux: Engine and plug-ins compilation

Execute the corresponding shell scripts in FES/FES:

* OSX 10.9.5 or later (Mavericks): [buildAllmavericks.sh](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Friend_Engine_Sources/buildAllMavericks.sh)
* OSX 10.6 to 10.8 (Mountain Lion): [buildAllMac.sh](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Friend_Engine_Sources/buildAllMac.sh)
* Linux systems installed with fslinstaller.py script provided by FSL  site: [buildAllLinux.sh](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Friend_Engine_Sources/buildAllLinux.sh)
* Debian/Ubuntu systems with NeuroDebian installation: [buildAllUbuntu.sh](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Friend_Engine_Sources/buildAllUbuntu.sh)

After compilation, the engine executable and plug-ins will be created in [FES/FES/Application](https://github.com/InstitutoDOr/FriendENGINE/tree/master/Friend_Engine_Sources/Friend_Engine_Sources/Application) folder.

5. Windows: Engine and plug-ins compilation

To compile on Windows, we suggest the Visual Studio Express 2013 for Windows suite, downloaded from: [http://www.visualstudio.com/en-us/products/visual-studio-express-vs.aspx](http://www.visualstudio.com/en-us/products/visual-studio-express-vs.aspx). After installation, users can open the solution file at [FES\Windows\buildAll\buildAll.sln](https://github.com/InstitutoDOr/FriendENGINE/blob/master/Friend_Engine_Sources/Windows/buildAll/buildAll.sln). After compilation, the engine.exe and plug-ins file will be created in the [FES\Windows\Application](https://github.com/InstitutoDOr/FriendENGINE/tree/master/Friend_Engine_Sources/Windows/Application) directory.

6. Engine execution

After executing, the engine will start listening communication in the default port 5678. To change the port, users can execute the engine passing the number of the desired port as a parameter, e.g.  engine 6000.
