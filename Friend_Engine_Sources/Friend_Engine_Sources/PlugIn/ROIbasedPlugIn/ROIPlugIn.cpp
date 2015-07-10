//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "ROIPlugIn.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object variables. This function brings the mni mask to subject space
int roiProcessing::initialization(studyParams &vdb)
{
	if (vdb.readedIni.IsEmpty())
	{
		fprintf(stderr, "the studyparams.txt file was not read.\n");
		return 1;
	}

	strcpy(vdb.mniMask, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMask"));
	strcpy(vdb.mniTemplate, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMaskReference"));
	targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel");
	int masktype = vdb.readedIni.GetLongValue("FRIEND", "ActivationLevelMaskType");

	fprintf(stderr, "mnimask = %s\n", vdb.mniMask);
	fprintf(stderr, "mnitemp = %s\n", vdb.mniTemplate);
	fprintf(stderr, "masktype = %d\n", masktype);
	if ((fileExists(vdb.mniMask)) && (fileExists(vdb.mniTemplate)) && (masktype == 2))
	{
		char outputFile[500], prefix[500] = "_RFI2", name[500];

		extractFileName(vdb.mniMask, name);
		for (int t = 0; t<strlen(name); t++)
			if (name[t] == '.') name[t] = '_';

		sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);

		fprintf(stderr, "Calculating the native template %s\n", outputFile);
		// bringing mni mask to subject space
		MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);

		// loads the reference mask
		meanCalculation.loadReference(outputFile);
	}
	else if ((fileExists(vdb.mniMask)) && (masktype == 1))
	{
		// loads the reference mask
		fprintf(stderr, "Loading native space mask %s\n", vdb.mniMask);
		meanCalculation.loadReference(vdb.mniMask);
	}    
	lastBaselineValue = 0;
    return 0;
}

// calculates the feedback value
int roiProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue)
{
   char processedFile[200];
   int idxInterval = vdb.interval.returnInterval(index);

   volume<float> v;
   // gets the motion corrected and gaussian file
   vdb.getMCGVolumeName(processedFile, index);
   read_volume(v, string(processedFile));

   // if in baseline condition, calculates the mean volume
   classnum = vdb.getClass(index);
   feedbackValue = 0;

   if (vdb.interval.isBaselineCondition(index))
   {
      if (vdb.interval.intervals[idxInterval].start == index) meanbaseline = v;
      else meanbaseline += v;

	  // returning the feedback value in the baseline condition. Nio problem if lasBaslineValue is 0
	  feedbackValue = PSC(meanCalculation.roiMean(0), lastBaselineValue);
	  feedbackValue = feedbackValue / targetValue;

	  // in the end of condition, divides the sum by the size of the current block
      if (vdb.interval.intervals[idxInterval].end == index) 
      {
         meanbaseline /= (vdb.interval.intervals[idxInterval].end-vdb.interval.intervals[idxInterval].start+1);
         meanCalculation.calculateMeans(meanbaseline);
         
         // calculates the mean roi value of the mean volume. We just need this value
         lastBaselineValue = meanCalculation.roiMean(0);
      }
   }
   else // task condition. Taking the mean of the volume and calculating the PSC
   {
	   meanCalculation.calculateMeans(v);
	   feedbackValue = PSC(meanCalculation.roiMean(0), lastBaselineValue);
	   feedbackValue = feedbackValue / targetValue;
      
      // enforcing 0..1 range
      //if (projection > 1) projection = 1;
      //else if (projection < 0) projection = 0;
	   fprintf(stderr, "Feedback value = %f\n", feedbackValue);
   }
   return 0;
}

// just calculate the percent signal change value
float roiProcessing::PSC(float value, float base)
{
   return (base) ? ((value-base) / base) : 0;
}

// plugin function for initializating the roi processing object
DLLExport initializeROIProcessing(studyParams &vdb, void *&userData)
{
   roiProcessing *roiVar = new roiProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

// plugin function for finalizing the object
DLLExport finalizeProcessing(studyParams &vdb, void *&userData)
{
   roiProcessing *roiVar = (roiProcessing *) userData;
   delete roiVar;
   userData = NULL;
   return 0;
}

// plugin function for calculationg feedback value
DLLExport processROI(studyParams &vdb, int index, float &classnum, float &feedbackValue, void * &userData)
{
   roiProcessing *roiVar = (roiProcessing *) userData;
   roiVar->processVolume(vdb, index, classnum, feedbackValue);
   return 0;
}
