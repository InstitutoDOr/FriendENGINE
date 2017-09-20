//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "CocainePlugIn.h"
#include "SimpleIni.h"
#include "fslfuncs.h"

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

	stringstream CmdLn;
	char auxString[500];

	baselineOffset = vdb.averageMeanOffset;
	minPSC = vdb.readedIni.GetLongValue("FRIEND", "minPSC", 0);
	maxPSC = vdb.readedIni.GetLongValue("FRIEND", "maxPSC", 0);
	WindowMaxPSC = vdb.readedIni.GetLongValue("FRIEND", "WindowMaxPSC", 60);

	fprintf(stderr, "Baseline offset=%d\n", baselineOffset);
	fprintf(stderr, "Minimum PSC=%d\n", minPSC);
	fprintf(stderr, "Maximum PSC=%d\n", maxPSC);

	char outputPSCFileName[500];
	sprintf(outputPSCFileName, "%soutputPSC%s.txt", vdb.outputDir, vdb.trainFeatureSuffix);
	outputPSC = fopen(outputPSCFileName, "w+");

	// reading from the study info file to see the subject feedback type
	char studyInfoFile[500];
	sprintf(studyInfoFile, "%s%cstudyinfo.txt", vdb.studyDir, PATHSEPCHAR);

	CSimpleIniA ini;
	SI_Error rc = ini.LoadFile(studyInfoFile);
	if (rc < 0) fprintf(stderr, "Error reading study Info File = %s\n", studyInfoFile);

	feedbackType = ini.GetLongValue(vdb.subject, "feedbackType", 1);

	switch (feedbackType)
	{
	case 1: fprintf(stderr, "Feedback type is NFB_C.\n");
		break;
	case 2: fprintf(stderr, "Feedback type is Non NFB_C.\n");
		break;
	case 3: fprintf(stderr, "Feedback type is NFB_P.\n");
		break;
	case 4: fprintf(stderr, "Feedback type is Non NFB_P.\n");
		break;
	}

	targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel");
	if ((feedbackType == 1) || (feedbackType == 3))
	{
		strcpy(vdb.mniMask, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMask"));
		strcpy(vdb.mniTemplate, vdb.readedIni.GetValue("FRIEND", "ActivationLevelMaskReference"));
	}
	else
	{
		strcpy(vdb.mniMask, vdb.readedIni.GetValue("FRIEND", "NonNFBLevelMask"));
		strcpy(vdb.mniTemplate, vdb.readedIni.GetValue("FRIEND", "NonNFBLevelMaskReference"));
	};

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

		// binarizes the roi mask
		changeFileExt(outputFile, "_bin.nii", auxString);
		CmdLn << "fslmaths " << outputFile << " -bin " << auxString << " -odt char";
		fslmaths((char *)CmdLn.str().c_str());

		fprintf(stderr, "Binarizing and loading roi mask : %s\n", auxString);
		// loads the reference mask
		meanCalculation.loadReference(auxString);
	}
	else if ((fileExists(vdb.mniMask)) && (masktype == 1))
	{
		// loads the reference mask
		fprintf(stderr, "Loading native space mask %s\n", vdb.mniMask);

		// binarizes the roi mask
		changeFileExt(vdb.mniMask, "_bin.nii", auxString);
		CmdLn << "fslmaths " << vdb.mniMask << " -bin " << auxString << " -odt char";
		fslmaths((char *)CmdLn.str().c_str());

		fprintf(stderr, "Binarizing and loading roi mask : %s\n", auxString);
		// loads the reference mask
		meanCalculation.loadReference(auxString);
	}
	else
	{
		if (!fileExists(vdb.mniMask))
			fprintf(stderr, "!!!!!!!   File not found %s.\n", vdb.mniMask);

		if (!fileExists(vdb.mniTemplate))
			fprintf(stderr, "!!!!!!!   File not found %s.\n", vdb.mniTemplate);
	};
	lastBaselineValue = 0;
	return 0;
}

void roiProcessing::recalculateLimits()
{
	int first = suprimePSCValues.size() - WindowMaxPSC;
	if (first < 0) first = 0;
	float minValue, maxValue;
	for (int i = first; i < suprimePSCValues.size(); i++)
	{
		if (i == first)
		{
			minValue = suprimePSCValues[i];
			maxValue = suprimePSCValues[i];
		}
		else
		{
			if (minValue > suprimePSCValues[i]) minValue = suprimePSCValues[i];
			if (maxValue < suprimePSCValues[i]) maxValue = suprimePSCValues[i];
		}
	};
	minPSC = minValue * 0.9;
	maxPSC = maxValue * 1.1;
}

// calculates the feedback value
int roiProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue)
{
   char processedFile[500], prefix[500];
   int idxInterval = vdb.interval.returnInterval(index);

   volume<float> v;

   // if in baseline condition, calculates the mean volume
   classnum = vdb.getClass(index);
   feedbackValue = 0;

   if (strcmp(vdb.interval.intervals[idxInterval].condition, "exp") == 0) // we are in an experience block
   {
	   if (index >= vdb.interval.intervals[idxInterval].start + baselineOffset)
	   {
		   // gets the motion corrected and gaussian file
		   vdb.getMCGVolumeFormat(prefix);

		   fprintf(stderr, "Generating activation file.\n");

		   // preparing the activation file name
		   vdb.setActivationFile(index);

		   // writing the sliding window file in activation file
		   estimateActivation(index, index, vdb.slidingWindowSize, prefix, vdb.maskFile, vdb.activationFile);

		   // reading this file
		   read_volume(v, string(vdb.activationFile));

		   if ((vdb.interval.intervals[idxInterval].start + baselineOffset) == index) meanbaseline = v;
		   else meanbaseline += v;

		   // in the end of condition, divides the sum by the size of the current block
		   if (vdb.interval.intervals[idxInterval].end == index)
		   {
			   meanbaseline /= (vdb.interval.intervals[idxInterval].end - vdb.interval.intervals[idxInterval].start + 1);
			   meanCalculation.calculateMeans(meanbaseline);

			   // calculates the mean roi value of the mean volume. We just need this value
			   lastBaselineValue = meanCalculation.roiMean(0);
		   }
	   }
   }

   if (strcmp(vdb.interval.intervals[idxInterval].condition, "sup") == 0) // we are in a suppression block
   {
	   if (index >= vdb.interval.intervals[idxInterval].start + vdb.offset)
	   {
		   // gets the motion corrected and gaussian file
		   vdb.getMCGVolumeFormat(prefix);

		   fprintf(stderr, "Generating activation file.\n");

		   // preparing the activation file name
		   vdb.setActivationFile(index);

		   // writing the sliding window file in activation file
		   estimateActivation(index, index, vdb.slidingWindowSize, prefix, vdb.maskFile, vdb.activationFile);

		   // reading this file
		   read_volume(v, string(vdb.activationFile));

		   // calculate the roi
		   meanCalculation.calculateMeans(v);

		   float psc = PSC(meanCalculation.roiMean(0), lastBaselineValue);
		   suprimePSCValues.push_back(psc);
		   feedbackValue = psc;
		   switch (feedbackType)
		   {
		   case 1: { // im just leaving it here to let you change this specific feedbacktype
			   recalculateLimits();
			   break;
		      }

		   case 2: { // im just leaving it here to let you change this specific feedbacktype
			   recalculateLimits();
			   break;
		      }

		   case 3: { // im just leaving it here to let you change this specific feedbacktype
			   break;
		      }


		   case 4: { // im just leaving it here to let you change this specific feedbacktype
			   break;
		      }
		   }
		   float pscSize = maxPSC - minPSC;
		   if (pscSize == 0) feedbackValue = 0;
		   else feedbackValue = (feedbackValue - minPSC) / (pscSize);
	   }
	   fprintf(stderr, "Feedback value = %f\n", feedbackValue);
   }
   fprintf(outputPSC, "%f\n", feedbackValue);
   fflush(outputPSC);
   return 0;
}

// just calculate the percent signal change value
float roiProcessing::PSC(float value, float base)
{
   return (base) ? ((value-base) / base) : 0;
}

// plug-in function for initializing roi processing object
DLLExport initializeROIProcessing(studyParams &vdb, void *&userData)
{
   roiProcessing *roiVar = new roiProcessing;
   roiVar->initialization(vdb);
   userData = roiVar;
   return 0;
}

// plug-in function for finalizing the object
DLLExport finalizeProcessing(studyParams &vdb, void *&userData)
{
   roiProcessing *roiVar = (roiProcessing *) userData;
   fclose(roiVar->outputPSC);
   delete roiVar;
   userData = NULL;
   return 0;
}

// plug-in function for calculating feedback value
DLLExport processROI(studyParams &vdb, int index, float &classnum, float &feedbackValue, void * &userData)
{
   roiProcessing *roiVar = (roiProcessing *) userData;
   roiVar->processVolume(vdb, index, classnum, feedbackValue);
   return 0;
}
