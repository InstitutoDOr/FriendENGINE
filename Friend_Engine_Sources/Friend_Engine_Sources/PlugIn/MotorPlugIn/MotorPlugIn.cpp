//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "MotorPlugIn.h"
#include "session.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object variables. This function brings the mni mask to subject space
int MotorRoiProcessing::initialization(studyParams &vdb)
{
   firstRoi = vdb.readedIni.GetLongValue("FRIEND", "firstROI");
   secondRoi = vdb.readedIni.GetLongValue("FRIEND", "secondROI");
   strcpy(vdb.mniMask, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMask"));
   strcpy(vdb.mniTemplate, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMaskReference"));
   fprintf(stderr, "mnimask = %s\n", vdb.mniMask);
   fprintf(stderr, "mnitemp = %s\n", vdb.mniTemplate);
   targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel");
   if ((fileExists(vdb.mniMask)) && (fileExists(vdb.mniTemplate)))
   {
      char outputFile[500], prefix[500]="_RFI2", name[500];
      
      extractFileName(vdb.mniMask, name);
      for (int t=0;t<strlen(name);t++)
      if (name[t] == '.') name[t] = '_';
      
      sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);
      
      fprintf(stderr, "Calculating the native template %s\n", outputFile);
      // bringing mni mask to subject space
      MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);
      
      // loads the reference mask
      meanCalculation.loadReference(outputFile);
	  firstRoiIndex  = meanCalculation.roiIndex(firstRoi);
	  secondRoiIndex = meanCalculation.roiIndex(secondRoi);
   }
   firstBaselineValue=0;
   secondBaselineValue = 0;
   return 0;
}

// calculates the feedback value
int MotorRoiProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &projection)
{
   char processedFile[200];
   int idxInterval = vdb.interval.returnInterval(index);
   float firstRoiMean, secondRoiMean;
   float firstProjection, secondProjection;

   volume<float> v;
   // gets the motion corrected and gaussian file
   vdb.getMCGVolumeName(processedFile, index);
   read_volume(v, string(processedFile));

   // if in baseline condition, calculates the mean volume, one by one
   classnum = vdb.getClass(index);
   projection = 0;

   if (vdb.interval.isBaselineCondition(index))
   {
      if (vdb.interval.intervals[idxInterval].start == index) meanbaseline = v;
      else meanbaseline += v;
      // in the end of condition, divides the sum by the size of the current block
      if (vdb.interval.intervals[idxInterval].end == index) 
      {
         meanbaseline /= (vdb.interval.intervals[idxInterval].end-vdb.interval.intervals[idxInterval].start+1);
         meanCalculation.calculateMeans(meanbaseline);
         
         // calculates the mean roi value of the mean volume. We just need this value
		 firstBaselineValue = meanCalculation.roiMean(firstRoiIndex);
		 secondBaselineValue = meanCalculation.roiMean(secondRoiIndex);
	  }
   }
   else // task condition. Taking the mean of the volume and calculating the PSC
   {
      meanCalculation.calculateMeans(v);
	  firstRoiMean = meanCalculation.roiMean(firstRoiIndex);
	  secondRoiMean = meanCalculation.roiMean(secondRoiIndex);

	  firstProjection = PSC(firstRoiMean, firstBaselineValue) / targetValue;
	  secondProjection = PSC(secondRoiMean, secondBaselineValue) / targetValue;
      projection = firstProjection;
	  char secondRoiString[100];
	  sprintf(secondRoiString, "%f", secondProjection);
	  fprintf(stderr, "%s\n", secondRoiString);

	  Session *session = (Session *)vdb.sessionPointer;
	  session->processAdditionalFeedBackInfo(index, secondRoiString);

      // enforcing 0..1 range
      //if (projection > 1) projection = 1;
      //else if (projection < 0) projection = 0;
      fprintf(stderr, "Projection value = %f\n", projection);
   }
   return 0;
}

// just calculate the percent signal change value
float MotorRoiProcessing::PSC(float value, float base)
{
   return (base) ? ((value-base) / base) : 0;
}

// plugin function for initializating the roi processing object
DLLExport initializeMotorProcessing(studyParams &vdb, void *&userData)
{
   MotorRoiProcessing *roiVar = new MotorRoiProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

// plugin function for finalizing the object
DLLExport finalizeMotorProcessing(studyParams &vdb, void *&userData)
{
   MotorRoiProcessing *roiVar = (MotorRoiProcessing *) userData;
   delete roiVar;
   return 0;
}

// plugin function for calculationg feedback value
DLLExport processMotorROI(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
   MotorRoiProcessing *roiVar = (MotorRoiProcessing *) userData;
   roiVar->processVolume(vdb, index, classnum, projection);
   return 0;
}
