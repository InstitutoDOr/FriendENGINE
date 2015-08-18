//
//  svmPlugIn.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____MDDPlugIn__
#define ____MDDPlugIn__

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
   WeightedMean correlationWeightedMean;
   vector<IncrementalStats> zNormRegionMeans;
   vector<double> zNormRegionMeanValues;
   IncrementalStats correlationStats;
   double correlationMultiplyer;
   double minCorrelationBaseline, maxCorrelationBaseline;
   int initialized;
   int volumeType;
   bool fixedLimit;
   char stabilizationCondition[100];
   
public:
   fstream outputReport;
   void processVolume(studyParams &vdb, int index, float &classnum, float &projection);
   // function that creates a roi map from glm and a previously defined roi in mni
   void createROIVolume(studyParams &vdb);
   
   // initializes the object with the config file information
   void initialize(int correlationWindow, int blockSize, int calculationWindowSize, char *regionfile, char *referencedir, char *correlationMapFile);
   
   // initializes the object with the config file information
   void initialize(studyParams &vdb);
   
   FunctionalConnectivity()
   {
      volumeType = 2; // Gaussian motion correction
      fixedLimit = 1;
      initialized = 0;
      correlationMultiplyer = 1;
   };
};

#endif /* defined(____MDDPlugIn__) */
