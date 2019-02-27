//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "MotivPlugIn.h"

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
   // these two following values are configured by the frontend via the SET command
   strcpy(roiMask, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMask"));
   int masktype = vdb.readedIni.GetLongValue("FRIEND", "ActivationLevelMaskType");
   feedbackType = vdb.readedIni.GetLongValue("FRIEND", "EmotionROIFeedbackType", 2);
   windowSize = vdb.readedIni.GetLongValue("FRIEND", "MinMaxWindowSize", 30);

   // reading the target values. If not defined, reading a default value ActivationLevel
   targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel", 0.005);

   // size of the last values that     
   int sigmoidSize = vdb.readedIni.GetDoubleValue("FRIEND", "SigmoidCalculationSize", 10);
   levelMultiplyer = vdb.readedIni.GetDoubleValue("FRIEND", "LevelMultiplyer", 0.25);

   fprintf(stderr, "mnimask = %s\n", roiMask);
   fprintf(stderr, "masktype = %d\n", masktype);
   if ((fileExists(roiMask)) && (masktype == 1))
   {
	   // loads the roi mask 
	   fprintf(stderr, "Loading native space mask %s\n", roiMask);
	   meanCalculation.loadReference(roiMask);
   }
   baseline=0;

   // currently we ar wnot using these two values
   blocksAboveTarget = 0;

   // initializating the object responsible for calculating the weighted correlation mean
   meanActivationLevel.setMeanWidth(sigmoidSize);
   meanActivationLevel.setVectorSize(vdb.interval.maxIndex());

   meanActivationLevel.initialize();
   meanActivationLevel.initializeCoefs();

   return 0;
}

// function that creates a roi map from glm and a previously defined roi in mni
void emotionRoiProcessing::createROIVolume(studyParams &vdb)
{
	// Bring MNI Mask to Subject Space
	char compositeFile[500], outputFile[500], name[500], regionExtractMapFile[500], roiVolumeFile[500], meanFile[500];
	char pngFile[500];
	float correlationExtractPercent;
	stringstream CmdLn;
	char prefix[30] = "_RFI2";
	vector<vector<double> > timeseries;

	if ((!fileExists(vdb.mniMask)) || (!fileExists(vdb.mniTemplate)))
	{
		if (!fileExists(vdb.mniTemplate)) fprintf(stderr, "Mni template file %s not found \n", vdb.mniTemplate);
		if (!fileExists(vdb.mniMask)) fprintf(stderr, "Mni mask file %s not found \n", vdb.mniMask);
		return;
	}

	vector<int> roiValues;
	meanCalculation.getRoiValues(roiValues);


	if (roiValues.size() > 0)
	{
		// generating the means file
		volume<float> v;

		// preparing timeseries variable. Two timeseries one for each roi
		timeseries.resize(meanCalculation.roiCount());
		for (int i = 0; i < timeseries.size(); i++)
			timeseries[i].resize(vdb.interval.maxIndex());

		// calculating the means of each roi, passing through all volumes ofthe current run
		for (int i = 0; i < vdb.interval.maxIndex(); i++)
		{
			loadVolume(vdb, v, i + 1);
			meanCalculation.calculateMeans(v);

			for (int j = 0; j < meanCalculation.roiCount(); j++)
				timeseries[j][i] = meanCalculation.roiMean(j);
		}

		// z normalizing the timeseries
		for (int j = 0; j < timeseries.size(); j++)	znormalise(timeseries[j]);

		// saving the result to file
		sprintf(meanFile, "%s%s_%s%s.txt", vdb.logDir, vdb.subject, "means", vdb.trainFeatureSuffix);
		fstream Output(meanFile, fstream::in | fstream::out | fstream::trunc);
		// actually saving the data
		for (int i = 0; i < vdb.interval.maxIndex(); i++)
		{
			for (int j = 0; j < timeseries.size(); j++)
				Output << timeseries[j][i] << '\t';
			Output << '\n';
		}
		// closing the file
		Output.close();

		// generating the roi means graph png of the calculated timeseries.
		fprintf(stderr, "Generating roi means graphic\n");
		changeFileExt(meanFile, ".png", pngFile);

		// First part of the command
		CmdLn << "fsl_tsplot -i " << meanFile << " -t \"z-normalised roi means plot\" -u 1 --start=1 --finish=" << meanCalculation.roiCount() << " -a ";
		// Building the labels with the intensities of the roi volume file
		int counter = 0;
		CmdLn << '\"';
		for (int t = 0; t < roiValues.size(); t++)
		{
			CmdLn << "Intensity " << roiValues[t];
			counter++;
			if (counter < meanCalculation.roiCount())
				CmdLn << ',';
		}

		// completing the command and executing it
		CmdLn << '\"';
		CmdLn << " -w 640 -h 144 -o " << pngFile;
		fsl_tsplot((char *)CmdLn.str().c_str());
		fprintf(stderr, "Creating mean roi timeseries graph.\n %s \n", CmdLn.str().c_str());
	}

	// bringing the reference roi mask in MNI space to native space
	// preparing the name of the outputfile
	extractFileName(vdb.mniMask, name);
	for (int t = 0; t<strlen(name); t++)
	if (name[t] == '.') name[t] = '_';

	sprintf(outputFile, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);

	// Executing the command
	MniToSubject(vdb.maskFile, vdb.mniMask, vdb.mniTemplate, outputFile, prefix);

	// saving a composite file just to assure the side of the roi volume. No real use in computation
	sprintf(compositeFile, "%s%s%s.nii", vdb.inputDir, "compositeFile", vdb.trainFeatureSuffix);
	uniteVolumes(vdb.rfiFile, outputFile, compositeFile);

	// Now create the best voxels roi mask

	// read Region extract map file and percentage
	RegionExtraction extractor;
	correlationExtractPercent = vdb.readedIni.GetDoubleValue("FRIEND", "PercentageOfBestVoxelsPerROI", 20.0) / 100.0;

	// Defining the output file name. This name wiil be used in the Unity frontend via the SET command
	sprintf(roiVolumeFile, "%s%s%s%s", vdb.outputDir, "ROIsMap", vdb.trainFeatureSuffix, ".nii");

	// Execute segmentation
	extractor.regionExtraction(outputFile, vdb.featuresAllTrainSuffix, roiVolumeFile, correlationExtractPercent);

	// calculate initial target values for positive and negative conditions
	if (fileExists(roiVolumeFile))
	{
		meanCalculation.loadReference(roiVolumeFile);

		for (int i=1; i <= vdb.interval.maxIndex(); i++)
		{
			float classnum, projection;
			processVolume(vdb, i, classnum, projection);
		}

		// saving a file with the roi means of the recently calculated roi mask
		char curveFile[BUFF_SIZE];
		sprintf(curveFile, "%s%c%s%s.txt", vdb.outputDir, PATHSEPCHAR, "curve", vdb.trainFeatureSuffix);

		meanActivationLevel.saveCurves(curveFile);

		// saving the values in the configuration file
		vdb.readedIni.SetDoubleValue("FRIEND", "ActivationLevel", meanActivationLevel.mean);
		vdb.readedIni.SaveFile(vdb.configFileNameRead);
	}
}

// load a volume for computation
void emotionRoiProcessing::loadVolume(studyParams &vdb, volume<float> &v, int index)
{
	char processedFile[200], prefix[BUFF_SIZE], tempVolume[BUFF_SIZE];

	// making the activation volume for viewing by applying the FRIEND pipeline sliding window mean. 
	// This file will not be used in the current processing
	vdb.getFinalVolumeFormat(prefix);
	vdb.setActivationFile(index);
	estimateActivation(index, index, vdb.slidingWindowSize, prefix, vdb.maskFile, vdb.activationFile);

	// gets the motion corrected and gaussian file and calculates the FRIEND pipeline sliding window mean
	sprintf(tempVolume, "%s%s", vdb.outputDir, "temp.nii");
	vdb.getMCGVolumeFormat(prefix);
	estimateActivation(index, index, vdb.slidingWindowSize, prefix, tempVolume);

	// if the file exists, read the volume and delete it.
	if (fileExists(tempVolume))
	{
		read_volume(v, string(tempVolume));
		remove(tempVolume);
	}
}

// calculates the feedback value depending of the current block 
int emotionRoiProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue)
{
   int idxInterval = vdb.interval.returnInterval(index);
   double pscValue;

   volume<float> v;


   loadVolume(vdb, v, index);
   // if in baseline condition, calculates the mean volume, one by one
   classnum = vdb.getClass(index);
   feedbackValue = 0;

   double intervalSize = (double)(vdb.interval.intervals[idxInterval].end - vdb.interval.intervals[idxInterval].start + 1);

   // in the baseline condition, just calculate the mean volume
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
         baseline = meanCalculation.roiMean(0);
		 fprintf(stderr, "baseline value = %f \n", baseline);
	  }
   }
   else // task condition. Taking the mean of the volume and calculating the PSC
   {
      meanCalculation.calculateMeans(v);
	  // Calculates the PSC value
	  pscValue = PSC(meanCalculation.roiMean(0), baseline);

	  if (feedbackType == 1)
	  {
		  // Calculates the percentage relative to the target
		  feedbackValue = pscValue / targetValue;

		  /////////////////////////The calculations are not used right now////////////////////////////////////////
		  //if (vdb.interval.intervals[idxInterval].start == index) negativePSCMean = negativePSC / intervalSize;
		  //else  negativePSCMean += negativePSC / intervalSize;
		  //
		  //if (vdb.interval.intervals[idxInterval].end == index)
		  //{
		  //	  if (negativePSCMean >= negativeTargetValue)
		  //		  negativeBlocksAboveTarget++;
		  //}
		  ////////////////////////////////////////////////////////////////////////////////////////////////////////

		  // only add positive PSC values to the sigmoid calculation
		  if (pscValue > 0)
			  meanActivationLevel.addValue(pscValue, 1);

		  // updating the target value for next volume
		  double valueLevel = meanActivationLevel.mean * (1 + levelMultiplyer);
		  if (valueLevel > 0) targetValue = valueLevel;
	  }
	  else if (feedbackType == 2)
	  {
		  meanActivationLevel.addValue(pscValue, 0);
		  fprintf(stderr, "PSC value = %f\n", pscValue);
		  feedbackValue = meanActivationLevel.scaleValue(pscValue, 0, windowSize);
	  }
	  else if (feedbackType == 3)
	  {
		  feedbackValue = pscValue / targetValue;
	  }
      fprintf(stderr, "Feedback value = %f\n", feedbackValue);
   }
   return 0;
}

// just calculate the percent signal change value
float emotionRoiProcessing::PSC(float value, float base)
{
   return (base) ? ((value-base) / base) : 0;
}

// plug-in function for initializing the roi processing object
DLLExport initializeEmotionROIProcessing(studyParams &vdb, void *&userData)
{
	emotionRoiProcessing *roiVar = new emotionRoiProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

// plug-in function for finalizing the object
DLLExport finalizeEmotionROIProcessing(studyParams &vdb, void *&userData)
{
	emotionRoiProcessing *roiVar = (emotionRoiProcessing *)userData;
   delete roiVar;
   userData = NULL;
   return 0;
}

// plug-in function for calculating feedback value
DLLExport processEmotionROI(studyParams &vdb, int index, float &classnum, float &feedbackValue, void * &userData)
{
	emotionRoiProcessing *roiVar = (emotionRoiProcessing *)userData;
   roiVar->processVolume(vdb, index, classnum, feedbackValue);
   return 0;
}

// plug-in function for building the ROI mask
DLLExport buildEmotionROI(studyParams &vdb, void * &userData)
{
	emotionRoiProcessing *roiVar = (emotionRoiProcessing *)userData;
	roiVar->createROIVolume(vdb);
	return 0;
}
