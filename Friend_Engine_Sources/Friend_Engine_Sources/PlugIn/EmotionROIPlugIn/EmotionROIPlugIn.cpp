//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "EmotionROIPlugIn.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object variables. This function brings the mni mask to subject space
int emotionRoiProcessing::initialization(studyParams &vdb)
{
	if (vdb.readedIni.IsEmpty())
	{
		fprintf(stderr, "the studyparams.txt file was not read.\n");
		return 1;
	}

   char roiMask[BUFF_SIZE];
   strcpy(roiMask, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMask"));
   int masktype = vdb.readedIni.GetLongValue("FRIEND", "ActivationLevelMaskType");
   targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel");
   int sigmoidSize = vdb.readedIni.GetDoubleValue("FRIEND", "SigmoidCalculationSize", 10);
   levelMultiplyer = vdb.readedIni.GetDoubleValue("FRIEND", "LevelMultiplyer", 0.25);
   negativeTargetValue = vdb.readedIni.GetDoubleValue("FRIEND", "NegativeActivationLevel", targetValue);
   positiveTargetValue = vdb.readedIni.GetDoubleValue("FRIEND", "PositiveActivationLevel", targetValue);

   fprintf(stderr, "mnimask = %s\n", roiMask);
   fprintf(stderr, "masktype = %d\n", masktype);
   if ((fileExists(roiMask)) && (masktype == 1))
   {
	   // loads the reference mask
	   fprintf(stderr, "Loading native space mask %s\n", roiMask);
	   meanCalculation.loadReference(roiMask);
	   positiveIndex = meanCalculation.mapping[3]; // 3 in the intensity value of positive emotion ROI
	   negativeIndex = meanCalculation.mapping[1]; // 1 in the intensity value of negative emotion ROI
	   fprintf(stderr, "Positive Index %d\n", positiveIndex);
	   fprintf(stderr, "Negative Index %d\n", negativeIndex);
   }
   positiveBaseline=0;
   negativeBaseline=0;
   positiveBlocksAboveTarget = 0;
   negativeBlocksAboveTarget = 0;

   // initializating the object responsible for calculating the weighted correlation mean
   positiveActivationLevel.setMeanWidth(sigmoidSize);
   positiveActivationLevel.setVectorSize(vdb.interval.maxIndex());

   positiveActivationLevel.initialize();
   positiveActivationLevel.initializeCoefs();

   negativeActivationLevel.setMeanWidth(sigmoidSize);
   negativeActivationLevel.setVectorSize(vdb.interval.maxIndex());

   negativeActivationLevel.initialize();
   negativeActivationLevel.initializeCoefs();

   return 0;
}

// function that creates a roi map from glm and a previously defined roi in mni
void emotionRoiProcessing::createROIVolume(studyParams &vdb)
{
	// Bring MNI Mask to Subject Space
	char outputFile[500], name[500], regionExtractMapFile[500], roiVolumeFile[500];
	double correlationExtractPercent;
	char prefix[30] = "_RFI2";

	extractFileName(vdb.mniMask, name);
	for (int t = 0; t<strlen(name); t++)
	if (name[t] == '.') name[t] = '_';

	sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);

	MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);

	// read Region extract map file and percentage
	correlationExtractPercent = vdb.readedIni.GetDoubleValue("FRIEND", "PercentageOfBestVoxelsPerROI") / 100.0;
	strcpy(regionExtractMapFile, vdb.readedIni.GetValue("FRIEND", "ROIExtractionMapFile"));


	RegionExtraction extractor;


	// reading the mappings
	std::map<int, int> mappings;
	extractor.readMappings(regionExtractMapFile, mappings);

	sprintf(roiVolumeFile, "%s%s%s%s", vdb.outputDir, "ROIsMap", vdb.trainFeatureSuffix, ".nii");

	// Execute segmentation
	extractor.regionsExtraction(outputFile, vdb.glmTOutput, vdb.featuresAllTrainSuffix, roiVolumeFile, mappings, correlationExtractPercent);

	// calculate initial target values for positive and negative conditions
	if (fileExists(roiVolumeFile))
	{
		meanCalculation.loadReference(roiVolumeFile);
		positiveIndex = meanCalculation.mapping[3]; // 3 in the intensity value of positive emotion ROI
		negativeIndex = meanCalculation.mapping[1]; // 1 in the intensity value of negative emotion ROI

		for (int i=1; i <= vdb.interval.maxIndex(); i++)
		{
			float classnum, projection;
			processVolume(vdb, i, classnum, projection);
		}

		char negativeCurveFile[BUFF_SIZE], positiveCurveFile[BUFF_SIZE];
		sprintf(negativeCurveFile, "%s%c%s%s.txt", vdb.outputDir, PATHSEPCHAR, "negative_curve", vdb.trainFeatureSuffix);
		sprintf(positiveCurveFile, "%s%c%s%s.txt", vdb.outputDir, PATHSEPCHAR, "positive_curve", vdb.trainFeatureSuffix);

		positiveActivationLevel.saveCurves(positiveCurveFile);
		negativeActivationLevel.saveCurves(negativeCurveFile);

		vdb.readedIni.SetDoubleValue("FRIEND", "NegativeActivationLevel", negativeActivationLevel.mean);
		vdb.readedIni.SetDoubleValue("FRIEND", "PositiveActivationLevel", positiveActivationLevel.mean);
		vdb.readedIni.SaveFile(vdb.configFileNameRead);
	}
}

// calculates the feedback value
int emotionRoiProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &projection)
{
	char processedFile[200], prefix[BUFF_SIZE], tempVolume[BUFF_SIZE];
   int idxInterval = vdb.interval.returnInterval(index);
   double positivePSC, negativePSC;

   volume<float> v;
   // gets the motion corrected and gaussian file

   if (0)
   {
	   vdb.getMCGVolumeName(processedFile, index);
	   read_volume(v, string(processedFile));
   }
   else
   {
	   sprintf(tempVolume, "%s%s", vdb.outputDir, "temp.nii");
	   vdb.getMCGVolumeFormat(prefix);
	   estimateActivation(index, index, vdb.slidingWindowSize, prefix, tempVolume);
	   read_volume(v, string(tempVolume));
	   if (fileExists(tempVolume))
		   remove(tempVolume);
   }

   // if in baseline condition, calculates the mean volume, one by one
   classnum = vdb.getClass(index);
   projection = 0;

   double intervalSize = (double)(vdb.interval.intervals[idxInterval].end - vdb.interval.intervals[idxInterval].start + 1);

   if (vdb.interval.isBaselineCondition(index))
   {
      if (vdb.interval.intervals[idxInterval].start == index) meanbaseline = v;
      else meanbaseline += v;
      // in the end of condition, divides the sum by the size of the current block
      if (vdb.interval.intervals[idxInterval].end == index) 
      {
         meanbaseline /= intervalSize;
         meanCalculation.calculateMeans(meanbaseline);
         
         // calculates the mean roi value of the mean volume. We just need this value
         positiveBaseline = meanCalculation.means[positiveIndex];
		 negativeBaseline = meanCalculation.means[negativeIndex];
	  }
   }
   else // task condition. Taking the mean of the volume and calculating the PSC
   {
      meanCalculation.calculateMeans(v);
      positivePSC = PSC(meanCalculation.means[positiveIndex], positiveBaseline);
	  negativePSC = PSC(meanCalculation.means[negativeIndex], negativeBaseline);

	  // negative emotion Feedback
	  if (classnum == 1)
	  {
		  projection = negativePSC / negativeTargetValue;
		  if (vdb.interval.intervals[idxInterval].start == index) negativePSCMean = negativePSC / intervalSize;
		  else  negativePSCMean += negativePSC / intervalSize;

		  if (vdb.interval.intervals[idxInterval].end == index)
		  {
			  if (negativePSCMean >= negativeTargetValue)
				  negativeBlocksAboveTarget++;
		  }

		  if (negativePSC > 0)
			  negativeActivationLevel.addValue(negativePSC, 1);

		  double valueLevel = negativeActivationLevel.mean * (1 + levelMultiplyer);
		  if (valueLevel > 0) negativeTargetValue = valueLevel;
	  }
	  else if (classnum == 3)
	  {
		  projection = positivePSC / positiveTargetValue;
		  if (vdb.interval.intervals[idxInterval].start == index) positivePSCMean = positivePSC / intervalSize;
		  else  positivePSCMean += positivePSC / intervalSize;

		  if (vdb.interval.intervals[idxInterval].end == index)
		  {
			  if (positivePSCMean >= positiveTargetValue)
				  positiveBlocksAboveTarget++;
		  }

		  if (positivePSC > 0)
			  positiveActivationLevel.addValue(positivePSC, 1);

		  double valueLevel = positiveActivationLevel.mean * (1 + levelMultiplyer);
		  if (valueLevel > 0) positiveTargetValue = valueLevel;
	  }
      fprintf(stderr, "Projection value = %f\n", projection);
   }
   return 0;
}

// just calculate the percent signal change value
float emotionRoiProcessing::PSC(float value, float base)
{
   return (base) ? ((value-base) / base) : 0;
}

// plugin function for initializating the roi processing object
DLLExport initializeEmotionROIProcessing(studyParams &vdb, void *&userData)
{
	emotionRoiProcessing *roiVar = new emotionRoiProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

// plugin function for finalizing the object
DLLExport finalizeEmotionROIProcessing(studyParams &vdb, void *&userData)
{
	emotionRoiProcessing *roiVar = (emotionRoiProcessing *)userData;
   delete roiVar;
   return 0;
}

// plugin function for calculationg feedback value
DLLExport processEmotionROI(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
	emotionRoiProcessing *roiVar = (emotionRoiProcessing *)userData;
   roiVar->processVolume(vdb, index, classnum, projection);
   return 0;
}

// plugin function for building the ROI mask
DLLExport buildEmotionROI(studyParams &vdb, void * &userData)
{
	emotionRoiProcessing *roiVar = (emotionRoiProcessing *)userData;
	roiVar->createROIVolume(vdb);
	return 0;
}
