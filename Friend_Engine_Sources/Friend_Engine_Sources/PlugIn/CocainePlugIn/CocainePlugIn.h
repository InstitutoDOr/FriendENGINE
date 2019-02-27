//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____cocainePlugIn__
#define ____cocainePlugIn__

#include <iostream>

#include "newimage/newimageall.h"
#include "fslfuncs.h"
#include "vardb.h"
#include "masks.h"
#include <sstream>

using namespace NEWIMAGE;

// class responsible for roi feedback processing, using percent signal change and a target value to calculate the termometer value
class roiProcessing
{
   RoiMeanCalculation meanCalculation;
   double lastBaselineValue, targetValue;
   volume<float> meanbaseline;
   
public:
   // initializes the object variables
   int initialization(studyParams &vdb);

   // calculates the feedback value
   int processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue);

   // caculates the minimum and maximum psc values based on the last volumes
   void recalculateLimits();
   // just calculate the percent signal change value
   float PSC(float value, float base);
   int feedbackType;
   int baselineOffset;
   float minPSC, maxPSC;
   int WindowMaxPSC;
   std::vector<float> suprimePSCValues;
   FILE *outputPSC;
};

#endif /* defined(____cocainePlugIn__) */

