//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____svmPlugIn__
#define ____svmPlugIn__

#include <iostream>

#include "newimage/newimageall.h"
#include "fslfuncs.h"
#include "vardb.h"
#include "masks.h"
#include <sstream>

using namespace NEWIMAGE;

// class responsible for roi feedback processing, using percent signal change and a target value to calculate the termometer value
class MotorRoiProcessing
{
   RoiMeanCalculation meanCalculation;
   double firstBaselineValue, secondBaselineValue, targetValue;
   volume<float> meanbaseline;
   int firstRoi, secondRoi;
   int firstRoiIndex, secondRoiIndex;

   
public:
   // initializes the object variables
   int initialization(studyParams &vdb);
   // calculates the feedback value
   int processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue);
   // just calculate the percent signal change value
   float PSC(float value, float base);
};

#endif /* defined(____svmPlugIn__) */

