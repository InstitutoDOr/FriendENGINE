//
//  PluginHandler.h
//  
//
//  Created by IDOR on 31/07/13.
//
//

#ifndef ____PluginHandler__
#define ____PluginHandler__
#include "vardb.h"

#ifndef WINDOWS
#include <dlfcn.h>
#define HMODULE void *
#endif

// plug in function pointer definitions
typedef int (*InitializationFunction)(studyParams &vardb, void * &userData);
typedef int (*TrainFunction)(studyParams &vardb, void *&userData);
typedef int (*TestFunction)(studyParams &vardb, int index, float &classnum, float &projection, void * &userData);
typedef int (*FinalizationFunction)(studyParams &vardb, void * &userData);
typedef int (*VolumeFunction)(studyParams &vardb, int index, char *fileName, void * &userData);
typedef int (*AfterPreprocessingFunction)(studyParams &vardb, int index, char *fileName, void * &userData);

// class responsible for handling all aspects of importing functions from an external library
class PluginHandler
{
   // function variables
   InitializationFunction initFunction;
   TrainFunction trainFunction;
   TestFunction testFunction;
   FinalizationFunction finalFunction;
   VolumeFunction volumeFunction;
   AfterPreprocessingFunction afterPreprocFunction;

   // pointer to the extra memory used by the plug-in
   void *userData;
   void* handler;
   // names of the library and functions
   char libraryPath[200], libraryFile[200], initFname[200], trainFname[200], testFname[200], finalFname[200], volumeFname[200], afterPreprocFname[200];
public:
	// variable that holds the logging facility
	LogObject *logObject;

   // loads an external library and import its functions
   void loadFunctions(char *library, char *trainFunc, char *testFunc, char *initFunc, char *finalFunc, char *volumeFunc, char *afterPreprocFunc);
   
   // set the library path
   void setLibraryPath(char *path);
   
   // call loadFunctions with default arguments
   void loadLibrary();
   
   // unloads the library and unreference the function pointers
   void unLoadLibrary();
   
   // functions just to call related plug in functions
   int callInitFunction(studyParams &vdb);
   int callTrainFunction(studyParams &vdb);
   int callTestFunction(studyParams &vdb, int index, float &classnum, float &projection);
   int callFinalFunction(studyParams &vdb);
   int callVolumeFunction(studyParams &vdb, int index, char *fileName);
   int callAfterPreprocessingFunction(studyParams &vdb, int index, char *fileName);

#ifndef __GNUC__
   PluginHandler ()
   {
      // setting default arguments
      libraryPath[0] = 0;
      strcpy(libraryFile, "libsvmp.dylib");

      strcpy(trainFname, "trainSVM");
      strcpy(testFname,  "testSVM");
      strcpy(initFname,  "initSVM");
      strcpy(finalFname,  "finalSVM");
      strcpy(volumeFname,  "no");
      strcpy(afterPreprocFname,  "no");
      
      trainFunction = NULL;
      testFunction = NULL;
      initFunction = NULL;
      finalFunction = NULL;
      volumeFunction = NULL;
      afterPreprocFunction = NULL;
      handler = NULL;
      userData = NULL;
   }
#else //The GCC way
   PluginHandler ()
   {
      // setting default arguments
      libraryPath[0] = 0;
      strcpy(libraryFile, "libsvmp.dylib");

      strcpy(trainFname, "trainSVM");
      strcpy(testFname,  "testSVM");
      strcpy(initFname,  "initSVM");
      strcpy(finalFname,  "finalSVM");
      strcpy(volumeFname,  "no");
      strcpy(afterPreprocFname,  "no");
      
      
      trainFunction = NULL;
      testFunction = NULL;
      initFunction = NULL;
      finalFunction = NULL;
      volumeFunction = NULL;
      afterPreprocFunction = NULL;
      handler = NULL;
      userData = NULL;
   }
   
#endif
};
#endif /* defined(____PluginHandler__) */
