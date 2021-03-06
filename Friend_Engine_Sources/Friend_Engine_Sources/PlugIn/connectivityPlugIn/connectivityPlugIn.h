//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____connectivityPlugIn__
#define ____connectivityPlugIn__

#include <iostream>
#include "fslfuncs.h"
#include "masks.h"
#include "vardb.h"
#include <sstream>
#include "defs.h"
#include <string>

// handles all the necessary steps for functional connectivity calculation
class FunctionalConnectivity
{
   RegionCorrelation connectivityCalculator;
   double minCorrelationBaseline, maxCorrelationBaseline;
   int initialized;
   int volumeType;
   
public:
   fstream outputReport;

   // function that calculates the sliding window correlation
   void processVolume(studyParams &vdb, int index, float &classnum, float &projection);

   // function that creates a roi map from glm and a previously defined roi in mni
   void createROIVolume(studyParams &vdb);

   // initializes the object with the config file information
   void initialize(studyParams &vdb);
   
   FunctionalConnectivity()
   {
      volumeType = 2; // Gaussian motion correction
      initialized = 0;
   };
};

#endif /* defined(____connectivityPlugIn__) */
