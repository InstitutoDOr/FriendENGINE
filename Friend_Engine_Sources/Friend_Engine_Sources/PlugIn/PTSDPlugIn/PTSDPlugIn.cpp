//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "PTSDPlugIn.h"
#include "session.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object variables. This function brings the mni mask to subject space
int PTSDProcessing::initialization(studyParams &vdb)
{
   // getting roi intensities
   firstRoi = vdb.readedIni.GetLongValue("FRIEND", "firstROI");
   secondRoi = vdb.readedIni.GetLongValue("FRIEND", "secondROI");


   // loads the region definition file
   char roiMask[500];
   strcpy(roiMask, vdb.readedIni.GetValue("FRIEND", "ROIMask"));

   targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel");
   // if exists the two files, bring  the mnimask to native space 
   if (fileExists(roiMask))
   {
      // loads the reference mask
	   meanCalculation.loadReference(roiMask);

	  // Theses indexes are the indexes relative to the roi intensities in the meanCalculation object mean array
	  firstRoiIndex  = meanCalculation.roiIndex(firstRoi);
	  secondRoiIndex = meanCalculation.roiIndex(secondRoi);
   }
   firstBaselineValue=0;
   secondBaselineValue = 0;
   return 0;
}

// calculates the feedback value
int PTSDProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue)
{
   char processedFile[200];
   int idxInterval = vdb.interval.returnInterval(index);
   float firstRoiMean, secondRoiMean;
   float firstRoiFeedbackValue, secondRoiFeedbackValue;

   classnum = vdb.getClass(index);
   feedbackValue = 0;

   // if in baseline condition (rest), calculates the mean volume
   if (vdb.interval.isBaselineCondition(index))
   {
	   volume<float> v;
	   // gets the motion corrected and gaussian smoothed file
	   vdb.getMCGVolumeName(processedFile, index);
	   read_volume(v, string(processedFile));

	   if (vdb.interval.intervals[idxInterval].start == index) meanbaseline = v;
      else meanbaseline += v;
      // in the end of condition, divides the sum by the size of the current block
      if (vdb.interval.intervals[idxInterval].end == index) 
      {
         meanbaseline /= (vdb.interval.intervals[idxInterval].end-vdb.interval.intervals[idxInterval].start+1);
         meanCalculation.calculateMeans(meanbaseline);
         
         // calculates the mean roi values of the mean volume. We just need this values
		 firstBaselineValue = meanCalculation.roiMean(firstRoiIndex);
		 secondBaselineValue = meanCalculation.roiMean(secondRoiIndex);
	  }
   }
   if ((strcmp(vdb.getCondition(index), "FEEDBACK") == 0) && (vdb.interval.intervals[idxInterval].start == index))
   // it could also be in any other posterior block between the task and feedback
   // we do the calculations only in the beginning of the interval
   {
	   int idxPreviousTaskRun = idxInterval - 2; // considering that the feedback and task run has two intervals apart (task, rate 1, rate 2, FB)
	   int endTaskIndex = vdb.interval.intervals[idxPreviousTaskRun].end; // index of the last volume in the previous task run
	   int taskBlockSize = 6; // you can consider not put here 8 to wait for hemodynamic delay
	   volume<float> taskMean;
	   char prefix[500];

	   // gets the motion corrected and gaussian smoothed file prefix
	   vdb.getMCGVolumeFormat(prefix);
	   vdb.setActivationFile(index); // or you could use endTaskIndex instead of index

	   // doing the sliding window in the mc gaussian files to generate the mean volume for the task block
	   // It will use taskBlockSize, not slidingWindow to calculate the mean.
	   estimateActivation(endTaskIndex, endTaskIndex, taskBlockSize, prefix, vdb.maskFile, vdb.activationFile);
	   read_volume(taskMean, string(vdb.activationFile));

	   // calculate the roi means
	   meanCalculation.calculateMeans(taskMean);

	   firstRoiMean = meanCalculation.roiMean(firstRoiIndex);
	   secondRoiMean = meanCalculation.roiMean(secondRoiIndex);

	   firstRoiFeedbackValue = PSC(firstRoiMean, firstBaselineValue) / targetValue;
	   secondRoiFeedbackValue = PSC(secondRoiMean, secondBaselineValue) / targetValue;
	   feedbackValue = firstRoiFeedbackValue;

	   // the first roi feedback value goes directly through the normal channel, but the second roi need another route
	   char secondRoiString[100];
	   sprintf(secondRoiString, "%f", secondRoiFeedbackValue);
	   fprintf(stderr, "%s\n", secondRoiString);

	   Session *session = (Session *)vdb.sessionPointer;
	   session->processAdditionalFeedBackInfo(index, secondRoiString);

	   fprintf(stderr, "Feedback value = %f\n", feedbackValue);
   }
   return 0;
}

// just calculate the percent signal change value
float PTSDProcessing::PSC(float value, float base)
{
   return (base) ? ((value-base) / base) : 0;
}

void PTSDProcessing::createROIVolume(studyParams &vdb)
{
	RoiMeanCalculation roiMask;

	// Bring MNI Mask to Subject Space
	char outputFile[500], name[500], roiVolumeFile[500];
	int roi1 = 1, roi2 = 2;

	// you should see, via contrast.txt file, which contrast index you wanna use. You can also get these values from a new variable you create in study_params.txt
	// like the ROIPercentage. Remeber that these values are zero based.
	int contrastRoi1 = 1, contrastRoi2 = 2;
	char prefix[30] = "_RFI2";

	extractFileName(vdb.mniMask, name);
	for (int t = 0; t<strlen(name); t++)
		if (name[t] == '.') name[t] = '_';

	sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);

	MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);

	// loading the roi Mask in native space
	roiMask.loadReference(outputFile);

	// reading the extraction percentage 
	float ROIPercent = vdb.readedIni.GetDoubleValue("FRIEND", "ROIPercentage", 100.0) / 100.0;

	int count;
	char **rois;


	RegionExtraction extractor;

	volume4D<float>glmResult;
	volume<float>outputVolume;
	int regionSize;

	// loads the GLM volume
	load_volume4D(glmResult, string(vdb.glmTOutput));

	outputVolume.reinitialize(glmResult.xsize(), glmResult.ysize(), glmResult.zsize());
	outputVolume.copyproperties(glmResult[0]);

	// this is important. Zeroing the volume
	outputVolume = 0;

	// counting the size of a roi 1
	regionSize = roiMask.roiSize(roiMask.roiIndex(roi1));

	// choosing the best voxels of that region
	extractor.regionBestVoxels(roiMask, glmResult[contrastRoi1], outputVolume, roi1, regionSize, ROIPercent);

	// counting the size of a roi 2
	regionSize = roiMask.roiSize(roiMask.roiIndex(roi2));

	// choosing the best voxels of that region		
	extractor.regionBestVoxels(roiMask, glmResult[contrastRoi2], outputVolume, roi2, regionSize, ROIPercent);

	// building the final name. You should set this 
	sprintf(roiVolumeFile, "%s%s%s%s", vdb.outputDir, "ROIsMap", vdb.trainFeatureSuffix, ".nii");

	save_volume(outputVolume, string(roiVolumeFile));
}

// specific ROIs of a MNI mask transformed in native space
DLLExport buildROIs(studyParams &vdb, void *&userData)
{
	PTSDProcessing* motor = (PTSDProcessing*)userData;
	motor->createROIVolume(vdb);
	return 1;
}

// plugin function for initializating the roi processing object
DLLExport initializePTSDProcessing(studyParams &vdb, void *&userData)
{
   PTSDProcessing *roiVar = new PTSDProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

// plugin function for finalizing the object
DLLExport finalizePTSDProcessing(studyParams &vdb, void *&userData)
{
	PTSDProcessing *roiVar = (PTSDProcessing *)userData;
   delete roiVar;
   return 0;
}

// plugin function for calculationg feedback value
DLLExport processPTSDROI(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
	PTSDProcessing *roiVar = (PTSDProcessing *)userData;
   roiVar->processVolume(vdb, index, classnum, projection);
   return 0;
}
