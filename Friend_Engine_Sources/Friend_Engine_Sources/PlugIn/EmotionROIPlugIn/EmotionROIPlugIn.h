//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____emotionPlugIn__
#define ____emotionPlugIn__

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
   WeightedMean positiveActivationLevel, negativeActivationLevel;
   double targetValue;
   volume<float> meanbaseline;
   int positiveIndex, negativeIndex;
   double positiveBaseline, negativeBaseline;
   double positivePSCMean, negativePSCMean;
   double positiveTargetValue, negativeTargetValue, levelMultiplyer;
   int negativeBlocksAboveTarget, positiveBlocksAboveTarget;
   
public:
   // initializes the object variables
   int initialization(studyParams &vdb);
   // calculates the feedback value
   int processVolume(studyParams &vdb, int index, float &classnum, float &projection);
   // just calculate the percent signal change value
   float PSC(float value, float base);
   // creates the subject space ROI mask filtering for the best voxels  
   void createROIVolume(studyParams &vdb);
};

#endif /* defined(____emotionPlugIn__) */

