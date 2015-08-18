//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____AttentionalPlugIn__
#define ____AttentionalPlugIn__

#include <iostream>

#include "newimage/newimageall.h"
#include "fslfuncs.h"
#include "vardb.h"
#include "masks.h"
#include <sstream>

using namespace NEWIMAGE;

// class responsible for roi feedback processing, using percent signal change and a target value to calculate the termometer value
class AttentionalProcessing
{
   double dmnMeanValue, tpnMeanValue, stValue, targetValue;
   volume<float> meanbaseline;
   int offset;
   float pscValues[20], baselineMeans[20];
   bool invertedFeedback;
   int roiIndexes[20];
   RoiMeanCalculation meanCalculation;

   
public:
   // mni mask in native space
   char mniMaskNativeSpace[500];
   // initializes the object variables
   int initialization(studyParams &vdb);
   // calculates the feedback value
   int processVolume(studyParams &vdb, int index, float &classnum, float &projection);
   // just calculate the percent signal change value
   float PSC(float value, float base);
   // load reference mask
   void loadReferenceMask(char *mask);
   // This functions tries to attenuate the effects of motion correction by using as reference the first
   // volume of a run. The mask used is adjusted to the new reference
   void setNewMCReference(studyParams &vdb, char *mask);

   FILE *f;
   char reportFile[500];
};

#endif /* defined(____AttentionalPlugIn__) */

