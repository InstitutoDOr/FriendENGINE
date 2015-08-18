//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "AttentionalPlugIn.h"
#include "session.h"
#include "SimpleIni.h"
#include "fslfuncs.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

// initializes the object variables. This function brings the mni mask to subject space
int AttentionalProcessing::initialization(studyParams &vdb)
{
	fprintf(stderr, "mnimask = %s\n", vdb.mniMask);
	fprintf(stderr, "mnitemp = %s\n", vdb.mniTemplate);
	offset = vdb.readedIni.GetDoubleValue("FRIEND", "AverageMeanOffset", 5);
	targetValue = vdb.readedIni.GetDoubleValue("FRIEND", "ActivationLevel");
	invertedFeedback = 0;

	// reading from the study info file to see the subject feedback type
	char studyInfoFile[500];
	sprintf(studyInfoFile, "%s%cstudyinfo.txt", vdb.studyDir, PATHSEPCHAR);

	CSimpleIniA ini;
	SI_Error rc = ini.LoadFile(studyInfoFile);
	if (rc < 0) fprintf(stderr, "Error reading study Info File = %s\n", studyInfoFile);

	invertedFeedback = ini.GetLongValue(vdb.subject, "invertedFeedback", 0);

	if (invertedFeedback) fprintf(stderr, "Feedback type is inverted.\n");
	else fprintf(stderr, "Feedback type is normal.\n");

	if ((fileExists(vdb.mniMask)) && (fileExists(vdb.mniTemplate)))
	{
		char prefix[500] = "_RFI2", name[500];

		extractFileName(vdb.mniMask, name);
		for (int t = 0; t<strlen(name); t++)
			if (name[t] == '.') name[t] = '_';

		sprintf(mniMaskNativeSpace, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);

		fprintf(stderr, "Calculating the native template %s\n", mniMaskNativeSpace);
		// bringing mni mask to subject space
		MniToSubject(vdb.rfiFile, vdb.mniMask, vdb.mniTemplate, mniMaskNativeSpace, prefix);

		loadReferenceMask(mniMaskNativeSpace);
	}

	// open the output file
	sprintf(reportFile, "%s\STvalues%s.txt", vdb.outputDir, vdb.trainFeatureSuffix);
	f = fopen(reportFile, "w+");
	return 0;
}

// calculates the feedback value
int AttentionalProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &projection)
{
	char processedFile[200], prefix[500];
	int idxInterval = vdb.interval.returnInterval(index);
	float firstRoiMean, secondRoiMean;
	float firstProjection, secondProjection;

	volume<float> v;
	// gets the motion corrected and gaussian file
	vdb.getMCGVolumeFormat(prefix);

	fprintf(stderr, "Generating activation file.\n");
	vdb.setActivationFile(index);
	estimateActivation(index, index, vdb.slidingWindowSize, prefix, vdb.maskFile, vdb.activationFile);

	vdb.getMCGVolumeName(processedFile, index);
	if (true) read_volume(v, string(processedFile));
	else read_volume(v, string(vdb.activationFile));

	// if in baseline condition, calculates the mean volume, one by one
	classnum = vdb.getClass(index);
	projection = 0;

	if (vdb.isBaselineCondition(index))
	{
		// this volume is discarded formean calculation
		if (index < (vdb.interval.intervals[idxInterval].start + offset)) return 0;
		// first volume to consider
		else if (index == (vdb.interval.intervals[idxInterval].start + offset)) meanbaseline = v;
		// the following ones
		else meanbaseline += v;

		// in the end of condition, divides the sum by the size of the current block
		if (vdb.interval.intervals[idxInterval].end == index)
		{
			meanbaseline /= (vdb.interval.intervals[idxInterval].end - (vdb.interval.intervals[idxInterval].start + offset) + 1);
			meanCalculation.calculateMeans(meanbaseline);

			// calculates the mean roi values of the mean volume. 
			for (int i = 0; i < 20; i++) baselineMeans[i] = meanCalculation.roiMean(roiIndexes[i]);
		}
	}
	else // task condition. Taking the mean of the volume and calculating the PSC
	{
		float stValue;
		meanCalculation.calculateMeans(v);
		for (int i = 0; i < 20; i++) pscValues[i] = PSC(meanCalculation.roiMean(roiIndexes[i]), baselineMeans[i]);

		dmnMeanValue = 0;
		for (int i = 0; i < 11; i++) dmnMeanValue += pscValues[i]/11;

		tpnMeanValue = 0;
		for (int i = 11; i < 20; i++) tpnMeanValue += pscValues[i]/9;

		if (invertedFeedback) stValue = (tpnMeanValue - dmnMeanValue);
		else stValue = (dmnMeanValue - tpnMeanValue);

		fprintf(stderr, "dmn = %f, tpn = %f, st = %f\n", dmnMeanValue, tpnMeanValue, stValue);
		projection = stValue / targetValue;

		fprintf(f, "Index = %d\n", index);

		fprintf(f, "Baseline Means :\n");
		for (int i = 0; i < 20; i++) fprintf(f, "%f ", baselineMeans[i]);
		fprintf(f, "\n");
		fprintf(f, "\n");

		fprintf(f, "Task Means :\n");
		for (int i = 0; i < 20; i++) fprintf(f, "%f ", meanCalculation.roiMean(roiIndexes[i]));
		fprintf(f, "\n");
		fprintf(f, "\n");

		fprintf(f, "ST value = %f\n", stValue);
		fflush(f);
		projection = (projection + 1) / 2;
	}
	return 0;
}

// load reference mask
void AttentionalProcessing::loadReferenceMask(char *mask)
{
	fprintf(stderr, "Loading mask volume = %s\n", mask);
	// loads the reference mask
	meanCalculation.loadReference(mask);

	// get the indexes of each roi intensity
	for (int i = 1; i <= 20; i++)
		roiIndexes[i - 1] = meanCalculation.roiIndex(i);
}

// This functions tries to attenuate the effects of motion correction by using as reference the first
// volume of a run. The mask used is adjusted to the new reference
void AttentionalProcessing::setNewMCReference(studyParams &vdb, char *fileName)
{
	strcpy(vdb.motionRefVolume, fileName);

	stringstream CmdLn;

	CmdLn << "bet " << fileName << " " << vdb.inputDir << "RFI_sks" << vdb.trainFeatureSuffix << " " << vdb.betParameters;
	bet((char *)CmdLn.str().c_str());

	char outputFile[500];

	changeFileExt(mniMaskNativeSpace, vdb.trainFeatureSuffix, outputFile);
	functionalNormalization(mniMaskNativeSpace, vdb.rfiFile, fileName, outputFile, true);
	loadReferenceMask(outputFile);
}

// just calculate the percent signal change value
float AttentionalProcessing::PSC(float value, float base)
{
	return (base) ? ((value - base) / base) : 0;
}

// plugin function for initializating the roi processing object
DLLExport initializeAttentionalProcessing(studyParams &vdb, void *&userData)
{
	AttentionalProcessing *roiVar = new AttentionalProcessing;
	roiVar->initialization(vdb);
	userData = roiVar;
	return 0;
}

// plugin function for finalizing the object
DLLExport finalizeAttentionalProcessing(studyParams &vdb, void *&userData)
{
	AttentionalProcessing *roiVar = (AttentionalProcessing *)userData;
	fclose(roiVar->f);
	delete roiVar;
	return 0;
}

// plugin function for calculationg feedback value
DLLExport processAttentional(studyParams &vdb, int index, float &classnum, float &projection, void * &userData)
{
	AttentionalProcessing *roiVar = (AttentionalProcessing *)userData;
	roiVar->processVolume(vdb, index, classnum, projection);
	return 0;
}


DLLExport changeMotionCorrectionReference(studyParams &vdb, int index, char *fileName, void * &userData)
{
	AttentionalProcessing *roiVar = (AttentionalProcessing *)userData;
	if (index == 1)
		roiVar->setNewMCReference(vdb, fileName);
	return 0;
}