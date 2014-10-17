//
//  vardb.h
//  
//
//  Created by Radiologia Rededor on 5/21/13.
//
//

#ifndef _vardb_h
#define _vardb_h
#include "defs.h"
#include "intervals.h"
#include "filefuncs.h"
#include "SimpleIni.h"
#include "socket.hpp"
#include "socket2.h"
#include "fslio/fslio.h"
#include "session.h"

// handles all information needed for FRIEND engine and plugins to work properly, like file names, config file parameters, communication socket.
class studyParams
{
public:
   // directory variables
   char studyDir[BUFF_SIZE], inputDir[BUFF_SIZE], logDir[BUFF_SIZE], outputDir[BUFF_SIZE], baseDir[BUFF_SIZE], glmDir[BUFF_SIZE],  preprocDir[BUFF_SIZE];
   
   // prefix and suffix variables
   char preprocVolumePrefix[BUFF_SIZE], mcSuffix[BUFF_SIZE], featuresSuffix[BUFF_SIZE], featuresAllSuffix[BUFF_SIZE], featuresAllTrainSuffix[BUFF_SIZE], featuresTrainSuffix[BUFF_SIZE], trainFeatureSuffix[BUFF_SIZE], testFeatureSuffix[BUFF_SIZE], rawVolumePrefix[BUFF_SIZE];
   
   // filename variables
   char subjectSpaceMask[BUFF_SIZE], parFile[BUFF_SIZE], rmsFile[BUFF_SIZE], fsfFile[BUFF_SIZE], activationFile[BUFF_SIZE], predictionFile[BUFF_SIZE], predictionFileAux[BUFF_SIZE], baseImage[BUFF_SIZE], matrixFile[BUFF_SIZE], maskFile[BUFF_SIZE], baseFunctional[BUFF_SIZE], contrastFile[BUFF_SIZE], train4DFile[BUFF_SIZE], trainGLM4DFile[BUFF_SIZE], glmMatrixFile[BUFF_SIZE], motionRefVolume[BUFF_SIZE], designFile[BUFF_SIZE], glmTOutput[BUFF_SIZE], actualBaseline[BUFF_SIZE];
   
   // parameter variables
   char baselineCondition[BUFF_SIZE];
   
   int actualImg, actualInterval;

   int offset, clusterSize, includeMotionParameters, referenceFirstVolSequence, byCutOff, conditionContrasts, averageMeanOffset, performSUSAN, slidingWindowSize, runSize;
   
   int performNeurofeedback, useWholeSubjectSpaceMask, referenceWholeVolume, storePredictions, feedBackType, numberWidth;
   int invX, invY, invZ;
   
   char mniTemplate[BUFF_SIZE], mniMask[BUFF_SIZE];
   
   double tTestCutOff, FWHM, trainingPercentage, percentileHigherVoxels;
   
   char mcflirtParams[BUFF_SIZE], betParameters[30];

   DesignObject interval;
   CSimpleIniA readedIni;
   
   // object used for sending information to frontend
   Socket2 socks;

   // pointer to the reference volume
   FSLIO *runReferencePtr;
   
   // control variables
   bool rIniRead, rPrepVars, rPreProc, rPipeline, rGLM, rFeatureSel, rTrain;

   char subject[BUFF_SIZE], raiFile[BUFF_SIZE], rfiFile[BUFF_SIZE];

   float TR;
   bool randomRun;
   
   // session variable
   Session *sessionPointer;

   // returns the final volume file name format. To retrieve the volume file name, you have to first call a sprintf function e.g. sprintf(fileName, format, 1);
   void getFinalVolumeFormat(char *format);
   
   // returns the prefix for the raw volumes (internal function)
   void _getPreprocVolumePrefix(char *prefix);
   
   // returns the prefix for the raw volumes
   void getPreprocVolumePrefix(char *prefix);
   
   // returns the sprintf format for the number part
   void getFormat(char *format);
   
   // return the final volume filename
   void getFinalVolumeName(char *outfile, char *number);
   void getFinalVolumeName(char *outfile, int number);

   // returns the motion corrected and gaussian filtered volume file names
   void getMCGVolumeName(char *outfile, char *number);
   void getMCGVolumeName(char *outfile, int number);

   // returns the motion corrected
   void getMCVolumeName(char *outfile, char *number);
   void getMCVolumeName(char *outfile, int number);
   
   // returns the volume file name given the number index in a string. Returns three types :
   // 1 - motion corrected volume,
   // 2 - motion corrected and gaussian filtered
   // 3 - motion corrected, gaussian filtered, baseline mean subtracted
   void getVolume(char *outfile, char *number, int type);
   void getVolume(char *outfile, int number, int type);

   // reads a config file
   void readConfigFile(char *configFile);
   
   // reads a config file from a buffer
   void readConfigBuffer(char *buffer, int size);
   
   // reads a config file in a SimpleIni variable
   void readConfig(CSimpleIniA &ini);

   // reads the information known by the studyParams object into its variables
   void readConfigVars();
   
   // saves a config information in a buffer to disk
   void saveConfigBuffer(char *buffer, int size, char *configFile);
   
   // initiates the internal variables after the config file load
   void prepRealtimeVars();
   
   // sets the value of a variable. Not used right now.
   void setVar(char *var, char *value);
   
   // returns the class of a volume index. It calls the respective DesignObject function
   int getClass(int idx);
   
   // determines, given a volume index, if its in a baseline block. Calls the respective DesignObject function
   bool isBaselineCondition(int idx);
   
   // get the condition name, given the volume index
   const char *getCondition(int idx);
   
   // sets the FRONT END response socket
   void setSocketfd(int Sock);
   
   // initializating control variables
   void initializeStates();
   
   studyParams () { mcflirtParams[0]='\0'; strcpy(mcflirtParams, "-plots -mats -rmsabs"); betParameters[0]='\0'; strcpy(betParameters, "-f 0.3 -o"); sessionPointer = NULL; };

    virtual ~studyParams() { };
    
    
};

#endif
