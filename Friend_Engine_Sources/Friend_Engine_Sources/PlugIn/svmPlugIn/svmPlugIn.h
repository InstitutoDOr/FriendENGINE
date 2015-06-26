//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____svmPlugIn__
#define ____svmPlugIn__

#include "fslfuncs.h"
#include "svmfuncs.h"
#include "svmobj.h"
#include "confusionmatrix.h"
#include <iostream>
#include <sstream>
#include "vardb.h"
#include <sys/stat.h>
#include <fstream>

// object responsible for handling all the SVM steps needed for
class SVMProcessing
{
   studyParams *vdbPtr;
   public:
   char svmTrainingFile[BUFF_SIZE], svmModelFile[BUFF_SIZE], svmModelPredictFile[BUFF_SIZE], svmWeightNormFile[BUFF_SIZE], svmWeightFile[BUFF_SIZE];
   char svmDir[BUFF_SIZE];
   int adaptTraining;
   int svmFeatureSelection;
   int cummulativeTraining;
   double minDistance, maxDistance, extrapolationFactor;
   int hasPredicted;
   svm_model *model;
   ConfusionMatrix accuracyResults;
   AdaptingSVM adaptativeTraining;
   FILE *svmLog, *projectionsFile;
   char logFileName[BUFF_SIZE], projectionsFilename[BUFF_SIZE];
   
   // initializes variables and creates the svm directory
   void initializeVars(studyParams &vdb);
   
   // handles the train step
   void train();
   
   // handles the train step
   void test(int index, char *volumeFile, float &classnum, float &projection);
   
   // deallocates the memory used
   void cleanUp();

   // initalize the logging
   void initializeLog();

   // write in the log
   void writeLog(char *msg);
};

#endif /* defined(____svmPlugIn__) */
