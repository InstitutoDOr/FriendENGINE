#ifndef ____motivPlugIn__
#define ____motivPlugIn__

#include <iostream>

#include "newimage/newimageall.h"
#include "fslfuncs.h"
#include "vardb.h"
#include "masks.h"
#include <sstream>

using namespace NEWIMAGE;

// class responsible for roi feedback processing, using percent signal change and a target value to calculate the termometer value
class emotionRoiProcessing
{
   RoiMeanCalculation meanCalculation;
   WeightedMean meanActivationLevel;
   double targetValue;
   volume<float> meanbaseline;
   double baseline;
   double pscMean;
   double levelMultiplyer;
   int blocksAboveTarget;
   
public:
   int feedbackType, windowSize;
   // load a volume for computation
   void loadVolume(studyParams &vdb, volume<float> &v, int index);
   // initializes the object variables
   int initialization(studyParams &vdb);
   // calculates the feedback value
   int processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue);
   // just calculate the percent signal change value
   float PSC(float value, float base);
   // creates the subject space ROI mask filtering for the best voxels  
   void createROIVolume(studyParams &vdb);
};

#endif /* defined(____motivPlugIn__) */

