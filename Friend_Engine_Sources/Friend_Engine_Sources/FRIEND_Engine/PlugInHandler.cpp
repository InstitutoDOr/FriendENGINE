//
//  PluginHandler.cpp
//  
//
//  Created by IDOR on 31/07/13.
//
//


#ifdef WINDOWS
#include <windows.h>
#define dlsym GetProcAddress
#endif
#include "PlugInHandler.h"

// call loadFunctions with default arguments
void PluginHandler::loadLibrary()
{
   char library[BUFF_SIZE];
   if (libraryPath)
      sprintf(library, "%s%c%s", libraryPath, PATHSEPCHAR, libraryFile);
   else
      strcpy(library, libraryFile);
   loadFunctions(library, trainFname, testFname, initFname, finalFname, volumeFname, afterPreprocFname);
}

// set the library path
void PluginHandler::setLibraryPath(char *path)
{
   strcpy(libraryPath, path);
}

// unloads the library and unreference the function pointers
void PluginHandler::unLoadLibrary()
{
   if (handler != NULL)
   {
#ifndef WINDOWS
      if (dlclose(handler))
         fprintf(stderr, "dlclose error : %s\n", dlerror());
#else
	   FreeLibrary((HMODULE) handler);
#endif
   }
   trainFunction = NULL;
   testFunction = NULL;
   initFunction = NULL;
   finalFunction = NULL;
   volumeFunction = NULL;
   afterPreprocFunction = NULL;
   userData = NULL;
}

// loads an external library and import its functions
void PluginHandler::loadFunctions(char *library, char *trainFunc, char *testFunc, char *initFunc, char *finalFunc, char *volumeFunc, char *afterPreprocFunc)
{
   char tempLibrary[1024];
   unLoadLibrary();
#ifdef WINDOWS
   sprintf(tempLibrary, "%s.dll", library);
   handler = (HMODULE) LoadLibraryA(tempLibrary);
#else
#ifdef LINUX
   sprintf(tempLibrary, "%s/%s.so", libraryPath, library); 
   handler = dlopen(tempLibrary, RTLD_LAZY);
#else
   sprintf(tempLibrary, "%s.dylib", library);
   handler = dlopen(tempLibrary, RTLD_LAZY);
#endif
#endif
   if (handler != NULL) // library loaded. trying to import the functions
   {
      fprintf(stderr, "Found library : %s\n", library);
      
      trainFunction = (TrainFunction) dlsym((HMODULE)handler, trainFunc);
      
      if (trainFunction == NULL) fprintf(stderr, "train function not loaded.\n");
      
      testFunction = (TestFunction) dlsym((HMODULE)handler, testFunc);
      
      if (testFunction == NULL) fprintf(stderr, "test function not loaded.\n");
      
      finalFunction = (FinalizationFunction) dlsym((HMODULE)handler, finalFunc);
      
      if (finalFunction == NULL) fprintf(stderr, "final function not loaded.\n");
      
      initFunction = (InitializationFunction) dlsym((HMODULE)handler, initFunc);
      
      if (initFunction == NULL) fprintf(stderr, "init function not loaded.\n");
      
      volumeFunction = (VolumeFunction) dlsym((HMODULE)handler, volumeFunc);
      
      if (volumeFunction == NULL) fprintf(stderr, "volume function not loaded.\n");
      
      afterPreprocFunction = (AfterPreprocessingFunction) dlsym((HMODULE)handler, afterPreprocFunc);
      
      if (afterPreprocFunction == NULL) fprintf(stderr, "after Preprocessing function not loaded.\n");
   }
   else // fail in importing library
   {
      fprintf(stderr, "Error loading library : %s\n", library);
#ifndef WINDOWS
      fprintf(stderr, "Error : %s\n", dlerror());
#endif      

      trainFunction = NULL;
      testFunction = NULL;
      initFunction = NULL;
      finalFunction = NULL;
      volumeFunction = NULL;
      afterPreprocFunction = NULL;
	  userData = NULL;
   }
}

// functions just to call related plug in functions
int PluginHandler::callInitFunction(studyParams &vdb)
{
   if (initFunction != NULL)
   {
	   if (userData == NULL)
         return initFunction(vdb, userData);
   }
   else return 0;
}

int PluginHandler::callVolumeFunction(studyParams &vdb, int index, char *fileName)
{
   if (volumeFunction != NULL)
      return volumeFunction(vdb, index, fileName, userData);
   else return 0;
}

int PluginHandler::callAfterPreprocessingFunction(studyParams &vdb, int index, char *fileName)
{
   if (afterPreprocFunction != NULL)
      return afterPreprocFunction(vdb, index, fileName, userData);
   else return 0;
}

int PluginHandler::callTrainFunction(studyParams &vdb)
{
   if (trainFunction != NULL)
      return trainFunction(vdb, userData);
   else return 0;
}

int PluginHandler::callTestFunction(studyParams &vdb, int index, float &classnum, float &projection)
{
   if (testFunction != NULL)
      return testFunction(vdb, index, classnum, projection, userData);
   else return 0;
}

int PluginHandler::callFinalFunction(studyParams &vdb)
{
	int r;
   if (finalFunction != NULL)
   {
	   r = finalFunction(vdb, userData);
	   userData = NULL;
	   return r;
   }
   else return 0;
}
