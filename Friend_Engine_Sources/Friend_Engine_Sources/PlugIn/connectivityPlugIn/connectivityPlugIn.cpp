//
//  connectivityPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "connectivityPlugIn.h"
#include "parser.h"
#include "filefuncs.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int _stdcall 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object with the config file information
void FunctionalConnectivity::initialize(studyParams &vdb)
{
   if (!initialized)
   {
      char regionfile[500], correlationMapFile[500], reportFileName[500];
      int correlationWindow, calculationWindowSize;
      
      // getting some variables from config information
      strcpy(regionfile, vdb.readedIni.GetValue("FRIEND", "RegionsFile"));
      
      strcpy(correlationMapFile, vdb.readedIni.GetValue("FRIEND", "CorrelationMapFile"));
   
      correlationWindow = vdb.readedIni.GetLongValue("FRIEND", "CorrelationWindowSize");
      
      calculationWindowSize = vdb.readedIni.GetLongValue("FRIEND", "CalculationWindowSize", correlationWindow);
   
      // this is just an example of how you build a plugin and this is a very specific calculation. This is an hard coded assignment.
      strcpy(stabilizationCondition, "INDIGNATION");
      initialize(correlationWindow, vdb.interval.maxIndex(), calculationWindowSize, regionfile, NULL, correlationMapFile);
      sprintf(reportFileName, "%smeans.txt", vdb.outputDir);
      
      outputReport.open(reportFileName, fstream::in | fstream::out | fstream::trunc);
      
      initialized = 1;
   }
}

// initializes the object with the config file information
void FunctionalConnectivity::initialize(int correlationWindow, int runSize, int calculationWindowSize, char *regionfile, char *referencedir, char *correlationMapFile)
{
    // loads the region definition file
    connectivityCalculator.loadRegions(regionfile, referencedir);
   
    // loads the only roi map used
    if (fileExists(correlationMapFile))
    {
        fprintf(stderr, "loading the map file = %s", correlationMapFile);
        connectivityCalculator.loadRegionMap(correlationMapFile);
    }
    else fprintf(stderr, "correlation map file %s not found.\n", correlationMapFile);

    // updates the window size
    connectivityCalculator.setWindowSize(correlationWindow);
   
    // initializating the object responsible for calculating the weighted correlation mean
    correlationWeightedMean.setMeanWidth(calculationWindowSize);
    correlationWeightedMean.setVectorSize(calculationWindowSize);

    correlationWeightedMean.initialize();
    correlationWeightedMean.initializeCoefs();
   
    // initializing the object responsible for the incremental stats (mean and variance) calculation of the correlations values
    correlationStats.initialize();
   
    // initializing the object responsible for calculating incrementally the mean and variance for z values calculations to display (in FRONT END) the mean values of different ROIs
    zNormRegionMeans.resize(connectivityCalculator.numRegions);
    zNormRegionMeanValues.resize(connectivityCalculator.numRegions);
    for (int t=0; t< connectivityCalculator.numRegions; t++) zNormRegionMeans[t].initialize();
}

void FunctionalConnectivity::processVolume(studyParams &vdb, int index, float &classnum, float &projection)
{
    char volumeFile[BUFF_SIZE];
   
    // calculate correlations
    vdb.getVolume(volumeFile, index, volumeType);
    connectivityCalculator.calculateCorrelations(volumeFile);

    // calculate zNorm ROI means for graphic display. This values are meant for FRONT END display. We do a incremental zcore calculation to avoid problems in graphic display causedby differente bold values levels in different rois.
    for (int t=0; t < connectivityCalculator.numRegions; t++)
    {
       double value = connectivityCalculator.means[t];
       outputReport << connectivityCalculator.means[t] << " ";
       zNormRegionMeans[t].addValue(value);
       zNormRegionMeanValues[t] = zNormRegionMeans[t].zValue(value);
    }
   
    // calculate mean correlation
    double correlationMean=0;
    for (int t=0; t < connectivityCalculator.numPairs; t++)
        correlationMean += connectivityCalculator.correlations[t];
   
    // remember here we have only one pair of rois to form a correlation. Other calculations with more than one pair are possible. You have the power to define it in your way.
    correlationMean /= connectivityCalculator.numPairs;
   
   // saving the results in a file, for debug purposes
    outputReport << correlationMean <<  correlationStats.mean <<correlationStats.deviation << "\n";
    outputReport.flush();
   
    // calculate weigted mean
    correlationWeightedMean.addValue(correlationMean);

    // calculates the incremental correlation stats
    correlationStats.addValue(correlationMean);

    if (true)
    {
       // Calculate a dynamic range to be used in coupling task
       minCorrelationBaseline = correlationWeightedMean.mean - correlationMultiplyer * correlationStats.deviation;
       maxCorrelationBaseline = correlationWeightedMean.mean + correlationMultiplyer * correlationStats.deviation;
    
       // get Classe
       classnum = vdb.getClass(index);
    
       // get projection (percentage)
       fprintf(stderr, "condition = %s\n", vdb.getCondition(index));
       bool randomLike = (strcmp(vdb.getCondition(index), stabilizationCondition) == 0) || vdb.randomRun;
       
       if (vdb.interval.isBaselineCondition(index))
       {
          fprintf(stderr, "Baseline condition.\n");
          projection = 0;
       }
       else if (randomLike) // stabilization mode. the participant is encouraged to maintain state.
       {
           fprintf(stderr, "Stabilization task.\n");
           fprintf(stderr, "Mean = %f mean 2 = %f\n", correlationMean, correlationWeightedMean.mean);
           projection = 1-abs(correlationMean-correlationWeightedMean.mean);
       }
       else // coupling task. the participant is encouraged to coupling anterior temporal lobe and subgenual areas. Here we use a dynamic range for transform the correlation in a percentage.
       {
           fprintf(stderr, "Coupling task.\n");
           fprintf(stderr, "Mean = %f min = %f max = %f \n", correlationMean, minCorrelationBaseline, maxCorrelationBaseline);
           projection = (correlationMean-minCorrelationBaseline)  / (maxCorrelationBaseline-minCorrelationBaseline);

       }
    }
   
    // enforcing value limits 0 and 1
    if (projection > 1) projection = 1;
    else if (projection < 0) projection = 0;
    fprintf(stderr, "Projection value = %f \n", projection);
}

// function that creates a roi map from glm and a previously defined roi in mni
void FunctionalConnectivity::createROIVolume(studyParams &vdb)
{
   // Bring MNI Mask to Subject Space
   char outputFile[500], name[500], regionExtractMapFile[500], roiVolumeFile[500];
   double correlationExtractPercent;
   char prefix[30]="_RFI2";
   
   extractFileName(vdb.mniMask, name);
   for (int t=0;t<strlen(name);t++)
      if (name[t] == '.') name[t] = '_';

   sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);
   
   MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);

   // read Region extract map file and percentage
   correlationExtractPercent = vdb.readedIni.GetDoubleValue("FRIEND", "RegionExtractionMapPerc") / 100.0;
   strcpy(regionExtractMapFile, vdb.readedIni.GetValue("FRIEND", "RegionExtractionMapFile"));
   

   RegionExtraction extractor;
   

   // reading the mappings
   std::map<int, int> mappings;
   extractor.readMappings(regionExtractMapFile, mappings);
   
   sprintf(roiVolumeFile, "%s%s%s%s", vdb.outputDir, "ROIsMap", vdb.trainFeatureSuffix, ".nii");
   
   // Execute segmentation
   extractor.regionsExtraction(outputFile, vdb.glmTOutput, vdb.featuresAllTrainSuffix, roiVolumeFile, mappings, correlationExtractPercent);
}

// plug in initialization function
DLLExport initializeFunctionalConectivity(studyParams &vdb, void *&userData)
{
   FunctionalConnectivity *fcconn = (FunctionalConnectivity *) userData;
   if (fcconn) return 1;
   fcconn = new FunctionalConnectivity();
   fcconn->initialize(vdb);
   userData = fcconn;
   return 1;
}

// plug in training function. Calls the createROIVolume function
DLLExport buildROIs(studyParams &vdb, void *&userData)
{
   FunctionalConnectivity *fcconn = (FunctionalConnectivity *) userData;
   fcconn->createROIVolume(vdb);
   return 1;
}

// plug in test function
DLLExport calculateFeedback(studyParams &vdb, int index, float &classnum, float &projection, void *&userData)
{
   FunctionalConnectivity *fcconn = (FunctionalConnectivity *) userData;
   fcconn->processVolume(vdb, index, classnum, projection);
   return 1;
}

// plug in finalization function
DLLExport finalizeFunctionalConectivity(studyParams &vdb, void *&userData)
{
   FunctionalConnectivity *fcconn = (FunctionalConnectivity *) userData;
   fcconn->outputReport.close();
   delete fcconn;
   return 1;
}

// plug in volume function
DLLExport volumeFunctionalConectivity(studyParams &vdb, int index, char *volume, void *&userData)
{
   //if (index==1)
   //   sprintf(vdb.motionRefVolume, "%s", volume);
   return 1;
}