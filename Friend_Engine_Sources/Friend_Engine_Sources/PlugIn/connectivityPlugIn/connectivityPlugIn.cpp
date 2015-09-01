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
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object with the config file information
void FunctionalConnectivity::initialize(studyParams &vdb)
{
   if (!initialized)
   {
	  // loads the region definition file
	  char roiMask[500];
	  strcpy(roiMask, vdb.readedIni.GetValue("FRIEND", "ROIMask"));
	  connectivityCalculator.initialize();
	  connectivityCalculator.loadRoiMask(roiMask);

	  // updates the window size
	  int correlationWindow = vdb.readedIni.GetLongValue("FRIEND", "CorrelationWindowSize");
	  connectivityCalculator.setWindowSize(correlationWindow);

      initialized = 1;
   }
}

void FunctionalConnectivity::processVolume(studyParams &vdb, int index, float &classnum, float &projection)
{
    char volumeFile[BUFF_SIZE];
   
	if (connectivityCalculator.numPairs > 0)
	{
		// calculate correlations based on the gaussian motion corrected volume
		vdb.getVolume(volumeFile, index, volumeType);
		// sliding window correlation
		connectivityCalculator.calculateCorrelations(volumeFile);
		// this object can actually calculate more than one correlation beween regions, but for now, just the first
		projection = connectivityCalculator.correlations[0];
	}
	else projection = 0;
	classnum = vdb.getClass(index);
}

// function that creates a roi map from glm and a previously defined roi in mni
void FunctionalConnectivity::createROIVolume(studyParams &vdb)
{
	RoiMeanCalculation roiMask;
	
	// Bring MNI Mask to Subject Space
	char outputFile[500], name[500], roiVolumeFile[500];
	char ROIIntensities[50];
    char prefix[30]="_RFI2";
   
    extractFileName(vdb.mniMask, name);
    for (int t=0;t<strlen(name);t++)
        if (name[t] == '.') name[t] = '_';

    sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);
   
    MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);

	// loading the roi Mask in native space
	roiMask.loadReference(outputFile);
    
	// reading the extraction percentage and the rois intensities
    float ROIPercent = vdb.readedIni.GetDoubleValue("FRIEND", "ROIPercentage", 100.0) / 100.0;
	strcpy(ROIIntensities, vdb.readedIni.GetValue("FRIEND", "ROIIntensities"));

	int count;
	char **rois;

	// separate the rois number from the string. The seoarator is a comma
	parser(ROIIntensities, count, rois, ',');
	if (count > 1)
	{
		RegionExtraction extractor;

		volume<float>glmResult, outputVolume;
		int region, regionSize;

		// loads the GLM volume
		load_volume(glmResult, string(vdb.featuresAllTrainSuffix));

		outputVolume.reinitialize(glmResult.xsize(), glmResult.ysize(), glmResult.zsize());
		outputVolume.copyproperties(glmResult);

		// this is important. Zeroing the volume
		outputVolume = 0;

		// just picking the first two roi numbers from the list
		for (int t = 0; t < 2; t++)
		{
			// transforming the char roi number into a number
			deblank(rois[t]);
			region = atoi(rois[t]);

			// counting the size of a roi intensities
			regionSize = roiMask.roiSize(roiMask.roiIndex(region));

			// choosing the best voxels of that region
			extractor.regionBestVoxels(roiMask, glmResult, outputVolume, region, regionSize, ROIPercent);
		}
		
		// building the final name
		sprintf(roiVolumeFile, "%s%s%s%s", vdb.outputDir, "ROIsMap", vdb.trainFeatureSuffix, ".nii");

		save_volume(outputVolume, string(roiVolumeFile));
	}
	freeparser(count, rois);
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

// this function creates the roi volume by selecting the best voxels based on the GLM result of the current run (localizer run) 
// specific ROIs of a MNI mask transformed in native space
DLLExport buildROIs(studyParams &vdb, void *&userData)
{
   FunctionalConnectivity *fcconn = (FunctionalConnectivity *) userData;
   fcconn->createROIVolume(vdb);
   return 1;
}

// this function calculates the sliding window correlation
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
   userData = NULL;
   return 1;
}
