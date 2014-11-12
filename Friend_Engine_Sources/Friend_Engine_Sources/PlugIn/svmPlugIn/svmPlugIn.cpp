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
   cummulativeTraining = vdb.readedIni.GetLongValue("FRIEND", "SVMCummulativeTraining", 0);

   snprintf(svmDir, buffSize, "%s%s%c", vdb.outputDir, "svm", PATHSEPCHAR);
   
   sprintf(svmTrainingFile, "%s%s%s%s", svmDir, "training", vdb.trainFeatureSuffix, ".txt");

   if (cummulativeTraining)
	   sprintf(svmTestingFile, "%s%s%s%s", svmDir, "cummulative_training", vdb.testFeatureSuffix, ".tst");
   else sprintf(svmTestingFile, "%s%s%s%s", svmDir, "training", vdb.testFeatureSuffix, ".tst");

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
	   double negativeCount = 0, positiveCount = 0;
	   while (getline(testFile, row))
	   {
		   double sampleClass, projection, predictedClass;
		   row = trim(row);
		   if (row.empty()) continue;
		   stringstream str(row);
		   str >> sampleClass >> projection >> predictedClass;

		   // taking the mean of the projections of the two classes as the maximum distance of each class
		   if (projection < 0)
		   {
			   negativeCount++;
			   minDistance += projection;
		   }
		   else if (projection > 0)
		   {
			   positiveCount++;
			   maxDistance += projection;
		   }
	   }
	   if (negativeCount != 0) minDistance /= negativeCount;
	   if (positiveCount != 0) maxDistance /= positiveCount;

	   testFile.close();
	   fprintf(stderr, "Distance limits %f %f\n", minDistance, maxDistance);
   }

   
   sprintf(svmModelFile, "%s%s%s%s", svmDir, vdb.subject, vdb.trainFeatureSuffix, ".model");
   
   if (cummulativeTraining)
	   sprintf(svmModelPredictFile, "%s%s%s%s%s", svmDir, "cummulative_", vdb.subject, vdb.testFeatureSuffix, ".model");
   else sprintf(svmModelPredictFile, "%s%s%s%s", svmDir, vdb.subject, vdb.testFeatureSuffix, ".model");

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
	char prefix[BUFF_SIZE], lastCummulativeTrainingFile[BUFF_SIZE], actualCummulativeTrainingFile[BUFF_SIZE], 
		cummulativeModelFile[BUFF_SIZE], cummulativeTestingFile[BUFF_SIZE];
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

   if (cummulativeTraining)
   {
	   sprintf(lastCummulativeTrainingFile, "%s%s%s%s", svmDir, "cummulative_training", vdbPtr->testFeatureSuffix, ".txt");
	   sprintf(actualCummulativeTrainingFile, "%s%s%s%s", svmDir, "cummulative_training", vdbPtr->trainFeatureSuffix, ".txt");

	   if (fileExists(lastCummulativeTrainingFile))
	   {
		   mergeFiles(svmTrainingFile, lastCummulativeTrainingFile, actualCummulativeTrainingFile);
	   }
	   else copyFile(svmTrainingFile, actualCummulativeTrainingFile);

	   sprintf(cummulativeModelFile, "%s%s%s%s%s", svmDir, "cummulative_", vdbPtr->subject, vdbPtr->trainFeatureSuffix, ".model");

	   // training
	   CmdLn.str("");
	   CmdLn << "svmtrain -t 0 " << actualCummulativeTrainingFile << " " << cummulativeModelFile;
	   svmObject.train(CmdLn.str().c_str());

	   CmdLn.str("");

	   sprintf(cummulativeTestingFile, "%s%s%s%s", svmDir, "cummulative_training", vdbPtr->trainFeatureSuffix, ".tst");

	   // testing the cummulative training data
	   CmdLn << "svmpredict " << actualCummulativeTrainingFile << " " << cummulativeModelFile << " " << cummulativeTestingFile;
	   svmObject.predict(CmdLn.str().c_str());
   }

   
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

   if (cummulativeTraining)
   {
	   char cummulativeSvmWeightNormFile[BUFF_SIZE], cummulativeSvmWeightFile[BUFF_SIZE];
	   sprintf(cummulativeSvmWeightNormFile, "%s%s%s%s", svmDir, "cummulative_weights_norm", vdbPtr->trainFeatureSuffix, ".nii");
	   sprintf(cummulativeSvmWeightFile, "%s%s%s%s", svmDir, "cummulative_weights", vdbPtr->trainFeatureSuffix, ".nii");

	   model = svm_load_model(cummulativeModelFile);
	   if (model != NULL)
	   {
		   generateWeightVolume(model, svmMask, 1, cummulativeSvmWeightNormFile);
		   generateWeightVolume(model, svmMask, 0, cummulativeSvmWeightFile);
		   unloadModel(model);
	   }
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

   sprintf(tempVolume, "%s%s",  vdb.outputDir, "temp.nii");
   
   vdb.getFinalVolumeFormat(prefix);
   
   remove(tempVolume);
   
   fprintf(stderr, "Generating tempfile.\n");
   estimateActivation(index, index, vdb.slidingWindowSize, prefix, tempVolume);
   
   fprintf(stderr, "Classifying.\n");
   svmProcessingVar->test(tempVolume, classnum, projection);
   
   if (fileExists(tempVolume))
      remove(tempVolume);
   return 1;
}