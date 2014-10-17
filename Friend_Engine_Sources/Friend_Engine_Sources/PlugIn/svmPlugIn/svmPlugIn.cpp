//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//
#ifdef WINDOWS
#include <windows.h>
#include <direct.h>
#define snprintf sprintf_s
#endif

#include "vardb.h"
#include "svmPlugIn.h"
#include <sys/stat.h>
#include <fstream>

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// remove spaces from string
std::string trim(const std::string &s);

// initializes variables and creates the svm directory
void SVMProcessing::initializeVars(studyParams &vdb)
{
   size_t buffSize = BUFF_SIZE-1;
   char svmTestingFile[BUFF_SIZE];
   
   // initializing specific variables
   svmFeatureSelection = vdb.readedIni.GetLongValue("FRIEND", "SVMFeatureSelection", 0);
   
   snprintf(svmDir, buffSize, "%s%s%c", vdb.outputDir, "svm", PATHSEPCHAR);
   
   sprintf(svmTrainingFile, "%s%s%s%s", svmDir, "training", vdb.trainFeatureSuffix, ".txt");
   sprintf(svmTestingFile, "%s%s%s%s", svmDir, "training", vdb.testFeatureSuffix, ".tst");

   fprintf(stderr, "%s\n", svmTrainingFile);
   fprintf(stderr, "%s\n", svmTestingFile);

   extrapolationFactor = 1.2;
   minDistance = 0;
   maxDistance = 0;
   if (fileExists(svmTestingFile))
   {
	   ifstream testFile;
	   testFile.open(svmTestingFile);
	   string row;
	   while (getline(testFile, row))
	   {
		   double sampleClass, hyperplaneDistance, predictedClass;
		   row = trim(row);
		   if (row.empty()) continue;
		   stringstream str(row);
		   str >> sampleClass >> hyperplaneDistance >> predictedClass;
		   minDistance = min(minDistance, hyperplaneDistance);
		   maxDistance = max(maxDistance, hyperplaneDistance);
	   }
	   testFile.close();
	   fprintf(stderr, "Distance limits %f %f\n", minDistance, maxDistance);
   }

   
   sprintf(svmModelFile, "%s%s%s%s", svmDir, vdb.subject, vdb.trainFeatureSuffix, ".model");
   
   sprintf(svmModelPredictFile, "%s%s%s%s", svmDir, vdb.subject, vdb.testFeatureSuffix, ".model");

   sprintf(svmWeightNormFile, "%s%s%s%s",  svmDir, "weights_norm", vdb.trainFeatureSuffix, ".nii");
           
   sprintf(svmWeightFile, "%s%s%s%s",  svmDir, "weights", vdb.trainFeatureSuffix, ".nii");
           
    // creating the svm directory
#ifdef WIN32
           _mkdir(svmDir);
#else
           mkdir(svmDir, 0777); // notice that 777 is different than 0777
#endif
   model = NULL;
   vdbPtr = &vdb;
}

// deallocates the memory used
void SVMProcessing::cleanUp()
{
   if (model!=NULL)
      unloadModel(model);
}

// handles the training
void SVMProcessing::train()
{
   char prefix[BUFF_SIZE];
   vector <int> classes, indices;
   
   std::stringstream CmdLn;
   
   vdbPtr->getFinalVolumeFormat(prefix);
   
   // calculating the sliding window and saving a 4D volume
   estimateActivation(1, vdbPtr->interval.maxIndex(), vdbPtr->slidingWindowSize, prefix, vdbPtr->train4DFile);
   
   // getting the volume indexes, excluding the baseline and the first ones discarted by haemodynamic stabilization
   vdbPtr->interval.getVolumeIndices(vdbPtr->offset, 1, vdbPtr->interval.maxIndex(), indices);
   
   // getting a vector containing the classes for each volume
   vdbPtr->interval.getClassArray(classes);
   
   char svmMask[BUFF_SIZE], svmTestingFile[BUFF_SIZE];
   sprintf(svmMask, "%s%s", vdbPtr->featuresSuffix, vdbPtr->trainFeatureSuffix);
   
   // transforms the 4D volume in a svm like input file
   saveSVMFile(vdbPtr->train4DFile, svmMask, svmTrainingFile, 0, indices, classes);
   
   CmdLn.str("");
   SVMObj svmObject;
   
   // training
   CmdLn << "svmtrain -t 0 " << svmTrainingFile << " " << svmModelFile;
   svmObject.train(CmdLn.str().c_str());

   CmdLn.str("");

   sprintf(svmTestingFile, "%s%s%s%s", svmDir, "training", vdbPtr->trainFeatureSuffix, ".tst");

   // testing the training data
   CmdLn << "svmpredict " << svmTrainingFile << " " << svmModelFile << " " << svmTestingFile;
   svmObject.predict(CmdLn.str().c_str());

   
   // Generating weight map volume
   fprintf(stderr, "Generating weight map volumes\n");
   
   sprintf(svmMask, "%s.nii.gz", vdbPtr->featuresTrainSuffix);

   fprintf(stderr, "using %s \n", svmMask);
   model=svm_load_model(svmModelFile);
   if (model != NULL)
   {
      generateWeightVolume(model, svmMask, 1, svmWeightNormFile);
      generateWeightVolume(model, svmMask, 0, svmWeightFile);
      unloadModel(model);
   }
}

// handles the testing
void SVMProcessing::test(char *volumeFile, float &classnum, float &projection)
{
   // getting the training mask name
   char featuresTestMask[BUFF_SIZE];
   sprintf(featuresTestMask, "%s%s.nii.gz",  vdbPtr->featuresSuffix, vdbPtr->testFeatureSuffix);

   if (model == NULL) model=svm_load_model(svmModelPredictFile);
   if (model == NULL) fprintf(stderr, "model file %s not loaded.\n", svmModelPredictFile);
   else
   {
      fprintf(stderr, "volumeFile tested : %s with mask : %s\n", volumeFile, featuresTestMask);
      predict(model, volumeFile, featuresTestMask, classnum, projection);

	  if (projection < 0) projection = projection / minDistance;
	  else projection = projection / maxDistance;
	  projection *= extrapolationFactor;
   }
}

// plug in initialization function
DLLExport initSVM(studyParams &vdb, void *&userData)
{
   SVMProcessing *svmProcessingVar = (SVMProcessing *) userData;
   if (svmProcessingVar) delete svmProcessingVar;
   svmProcessingVar = new SVMProcessing;
   svmProcessingVar->initializeVars(vdb);
   userData = svmProcessingVar;
   return 0;
}

// plug in finalization function
DLLExport finalSVM(studyParams &vdb, void *&userData)
{
   if (userData != NULL)
   {
      SVMProcessing *svmProcessingVar = (SVMProcessing *) userData;
      svmProcessingVar->cleanUp();
      delete svmProcessingVar;
      userData = NULL;
   }
   return 0;
}

// plug in training function
DLLExport trainSVM(studyParams &vdb, void *&userData)
{
   SVMProcessing *svmProcessingVar = (SVMProcessing *) userData;
   svmProcessingVar->train();
   return 0;
}

// plug in test function
DLLExport testSVM(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
   char tempVolume[BUFF_SIZE], prefix[BUFF_SIZE];
   SVMProcessing *svmProcessingVar = (SVMProcessing *) userData;

   sprintf(tempVolume, "%s%s",  vdb.outputDir, "temp.nii.gz");
   
   vdb.getFinalVolumeFormat(prefix);
   
   remove(tempVolume);
   
   fprintf(stderr, "Generating tempfile.\n");
   estimateActivation(index, index, vdb.slidingWindowSize, prefix, tempVolume);
   
   fprintf(stderr, "Classifying.\n");
   svmProcessingVar->test(tempVolume, classnum, projection);
   
   remove(tempVolume);
   return 1;
}