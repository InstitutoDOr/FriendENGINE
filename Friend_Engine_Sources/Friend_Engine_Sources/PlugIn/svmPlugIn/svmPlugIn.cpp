//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//
#include "svmPlugIn.h"

#ifdef WINDOWS
#include <direct.h>
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
   hasPredicted = 0;
   svmFeatureSelection = vdb.readedIni.GetLongValue("FRIEND", "SVMFeatureSelection", 0);
   cummulativeTraining = vdb.readedIni.GetLongValue("FRIEND", "SVMCummulativeTraining", 0);
   adaptTraining = vdb.readedIni.GetLongValue("FRIEND", "AdaptingSVM", 0);

   sprintf(svmDir, "%s%s%c", vdb.outputDir, "svm", PATHSEPCHAR);
   
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
   accuracyResults.setNeutralClassCount(1); // we expect one neutral class
   accuracyResults.setRank(vdb.interval.conditionNames.size());
   for (int t = 0; t < vdb.interval.conditionNames.size(); t++)
	   accuracyResults.setClassName(t + 1, vdb.interval.conditionNames[t]);

   if (adaptTraining)
	   adaptativeTraining.initialize(svmModelPredictFile);

   sprintf(logFileName, "%slogfile%s_performance.txt", vdb.logDir, vdb.trainFeatureSuffix);
   sprintf(projectionsFilename, "%sprojections%s.txt", svmDir, vdb.trainFeatureSuffix);
   projectionsFile = 0;
   svmLog = 0;
}

// deallocates the memory used
void SVMProcessing::cleanUp()
{
	fprintf(stderr, "Unloading model, if any loaded.\n");
   if (model!=NULL) unloadModel(model);

   fprintf(stderr, "Closing log file, if open.\n");
   if (svmLog)
   {
	   fclose(svmLog);
	   svmLog = 0;
   }

   if (projectionsFile)
   {
	   fclose(projectionsFile);
	   svmLog = 0;
   }
   fprintf(stderr, "Clean up finished.\n");
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
   generateProjetionsGraph(svmTestingFile);

   // testing the training data with the prediction model
   if (fileExists(svmModelPredictFile))
   {
	   CmdLn.str("");
	   sprintf(svmTestingFile, "%s%s%s%s", svmDir, "testing", vdbPtr->trainFeatureSuffix, ".tst");
	   CmdLn << "svmpredict " << svmTrainingFile << " " << svmModelPredictFile << " " << svmTestingFile;
	   svmObject.predict(CmdLn.str().c_str());
	   generateProjetionsGraph(svmTestingFile);
   }

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
	   generateProjetionsGraph(cummulativeTestingFile);

	   if (fileExists(svmModelPredictFile))
	   {
		   CmdLn.str("");
		   sprintf(cummulativeTestingFile, "%s%s%s%s", svmDir, "cummulative_testing", vdbPtr->trainFeatureSuffix, ".tst");
		   CmdLn << "svmpredict " << actualCummulativeTrainingFile << " " << svmModelPredictFile << " " << cummulativeTestingFile;
		   svmObject.predict(CmdLn.str().c_str());
		   generateProjetionsGraph(cummulativeTestingFile);
	   }
   }

   
   // Generating weight map volume
   fprintf(stderr, "Generating weight map volumes\n");
   
   sprintf(svmMask, "%s.nii", vdbPtr->featuresTrainSuffix);

   fileExists(svmMask);

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

// initalize the logging
void SVMProcessing::initializeLog()
{
	svmLog = fopen(logFileName, "wt+");
	if (svmLog)
	{
		fprintf(svmLog, "SUBJECT     : %s\n", vdbPtr->subject);
		fprintf(svmLog, "Design File : %s\n", vdbPtr->designFile);
		fprintf(svmLog, "Model File  : %s\n", svmModelPredictFile);
	}
	projectionsFile = fopen(projectionsFilename, "wt+");
}

// write in the log
void SVMProcessing::writeLog(char *msg)
{
	if (svmLog)
	{
		fprintf(svmLog, "%s", msg);
		//fflush(svmLog);
	}
}

// handles the testing
void SVMProcessing::test(int index, char *volumeFile, float &classnum, float &projection)
{
	char msg[500], result[10];
	if (svmLog == 0) initializeLog();

	int idxInterval = vdbPtr->interval.returnInterval(index);
	if (index == vdbPtr->interval.intervals[idxInterval].start)
	{
		sprintf(msg, "\nBeginning of interval %d\n", (idxInterval+1));
		writeLog(msg);
	}
	// getting the training mask name
   char featuresTestMask[BUFF_SIZE];
   sprintf(featuresTestMask, "%s%s.nii",  vdbPtr->featuresSuffix, vdbPtr->testFeatureSuffix);

   // correcting the file name
   fileExists(featuresTestMask);

   if (model == NULL) model=svm_load_model(svmModelPredictFile);
   if (model == NULL) fprintf(stderr, "model file %s not loaded.\n", svmModelPredictFile);
   else
   {
	  //fprintf(stderr, "volumeFile tested : %s with mask : %s\n", volumeFile, featuresTestMask);
	  predict(model, volumeFile, featuresTestMask, classnum, projection);
	  if (adaptTraining)
	  {
		  fprintf(stderr, "Adapting projection.\n");
		  adaptativeTraining.adaptResult(classnum, projection);
	  }

	  if (projectionsFile)
		  fprintf(projectionsFile, "%f\n", projection);

	  hasPredicted = 1;

	  int actualClass = vdbPtr->getClass(index);
	  int predicted = (int)classnum;
	  int idxInterval = vdbPtr->interval.returnInterval(index);
	  if (!vdbPtr->interval.isBaselineCondition(vdbPtr->interval.intervals[idxInterval].condition))
	  {
		  if (index >= vdbPtr->interval.intervals[idxInterval].start + vdbPtr->offset)
			  accuracyResults.reportResult(actualClass, predicted);

		  if (actualClass == predicted) sprintf(result, "OK");
		  else sprintf(result, "NOK");

		  sprintf(msg, "Classification : %s %d %s : %f\n", vdbPtr->interval.conditionNames[predicted - 1].c_str(), index, result, projection);
		  writeLog(msg);
	  }

	  if (index == vdbPtr->interval.intervals[idxInterval].end)
	  {
		  if (accuracyResults.examples()) sprintf(result, "%.2f", accuracyResults.hits());
		  else sprintf(result, "0.00");
		  strcat(result, "%");

		  sprintf(msg, "Ending of interval %d. Accuracy so far %d in %d (%s)\n", (idxInterval + 1), accuracyResults.correctExamples(), accuracyResults.examples(), result);
		  writeLog(msg);
	  }

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
	  if (svmProcessingVar->hasPredicted)
	  {
		  char confusionMatrixFile[500];
		  sprintf(confusionMatrixFile, "%s%s%s.txt", svmProcessingVar->svmDir, "confusionMatrix", vdb.trainFeatureSuffix);
		  svmProcessingVar->accuracyResults.saveMatrixReport(confusionMatrixFile);

		  if (fileExists(svmProcessingVar->projectionsFilename))
		  {
			  stringstream CmdLn;
			  char pngFile[BUFF_SIZE];

			  // generating the projections graph png of the actual run
			  changeFileExt(svmProcessingVar->projectionsFilename, ".png", pngFile);
			  fprintf(stderr, "Generating svm projection graphics of the current run performance\n");
			  CmdLn << "fsl_tsplot -i " << svmProcessingVar->projectionsFilename << " -t \"SVM projections\" -u 1 --start=1 --finish=1 -a Projections -w 640 -h 144 -o " << pngFile;

			  fsl_tsplot((char *)CmdLn.str().c_str());

		  }
	  }
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

/*
// plug in test function
DLLExport testSVM(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
	char tempVolume[BUFF_SIZE], prefix[BUFF_SIZE];
	SVMProcessing *svmProcessingVar = (SVMProcessing *)userData;

	sprintf(tempVolume, "%s%s", vdb.outputDir, "temp.nii");

	vdb.getFinalVolumeFormat(prefix);

	remove(tempVolume);

	vdb.setActivationFile(index);

	fprintf(stderr, "Generating tempfile.\n");
	estimateActivation(index, index, vdb.slidingWindowSize, prefix, tempVolume);

	fprintf(stderr, "Classifying.\n");
	svmProcessingVar->test(index, tempVolume, classnum, projection);

	if (fileExists(tempVolume))
		remove(tempVolume);
	return 1;
}
*/

// plug in test function
DLLExport testSVM(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
   char prefix[BUFF_SIZE];
   SVMProcessing *svmProcessingVar = (SVMProcessing *) userData;

   vdb.getFinalVolumeFormat(prefix);

   fprintf(stderr, "Generating activation file.\n");
   vdb.setActivationFile(index);

   estimateActivation(index, index, vdb.slidingWindowSize, prefix, vdb.maskFile, vdb.activationFile);
   
   fprintf(stderr, "Classifying.\n");
   svmProcessingVar->test(index, vdb.activationFile, classnum, projection);
   
   return 1;
}