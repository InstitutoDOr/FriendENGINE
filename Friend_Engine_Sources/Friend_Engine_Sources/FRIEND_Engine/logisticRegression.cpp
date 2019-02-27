//
//  svmPlugIn.cpp
//  
//
//  Created by IDOR on 06/06/13.
//
//

#include "logisticRegression.h"

#ifdef WINDOWS
#define DLLExport extern "C" __declspec(dllexport) int 
#else
#define DLLExport extern "C" int 
#endif

void NormalizeVector(Matrix &lastValues, Matrix &columnVector, int lastValuesCount=30)
{
	for (int t = 1; t <= lastValues.Nrows(); t++)
	{
		double minv, maxv;

		// including in this window the last point
		lastValuesCount++;
		// calculating the firstIndex to consider to build the min max range
		int firstIndex = max(1, lastValues.Ncols() - lastValuesCount);
		for (int i = firstIndex; i <= lastValues.Ncols(); i++)
		{
			if (i == 1)
			{
				minv = lastValues(t, i);
				maxv = lastValues(t, i);
			}
			else
			{
				minv = min(minv, lastValues(t, i));
				maxv = max(maxv, lastValues(t, i));
			}
		}
		if (minv == maxv) columnVector(t, 1) = 0;
		else columnVector(t, 1) = (columnVector(t, 1) - minv) / (maxv - minv);
	}
}

int MaskSize(NEWIMAGE::volume<float> &mask)
{
	int i = 0;
	for (int z = 0; z < mask.zsize(); z++)
		for (int y = 0; y < mask.ysize(); y++)
			for (int x = 0; x < mask.xsize(); x++)
				if (mask.value(x, y, z) > 0)
				{
					i++;
				}
	return i;
}

void Volume2Matrix(volume<float> &vol, volume<float>&mask, int masksize, Matrix &sample)
{
	sample.ReSize(masksize + 1, 1);

	int i = 2;
	// first element of the column vector is the intercept
	sample(1, 1) = 1;

	// building the vector based on the mask
	for (int z = 0; z < mask.zsize(); z++)
		for (int y = 0; y < mask.ysize(); y++)
			for (int x = 0; x < mask.xsize(); x++)
				if (mask.value(x, y, z) > 0)
					sample(i++, 1) = vol.value(x, y, z);
}

// initializes the object variables. This function brings the mni mask to subject space
int logisticRegressionProcessing::initialization(studyParams &vdb)
{
	if (vdb.readedIni.IsEmpty())
	{
		vdb.logObject->writeLog(1, "the studyparams.txt file was not read.\n");
		return 1;
	}

	char theta1FileName[1024], theta2FileName[1024], mask1Filename[1024], mask2Filename[1024];

	strcpy(theta1FileName, vdb.readedIni.GetValue("FRIEND", "Theta1Matrix"));
	strcpy(theta2FileName, vdb.readedIni.GetValue("FRIEND", "Theta2Matrix"));

	strcpy(mask1Filename, vdb.readedIni.GetValue("FRIEND", "Mask1Filename"));
	strcpy(mask2Filename, vdb.readedIni.GetValue("FRIEND", "Mask2Filename"));

	// reading the masks
	if (fileExists(mask1Filename))
	{
		readMask1(mask1Filename);
	}
	else mask1Size = 0;

	if (fileExists(mask2Filename))
	{
		readMask2(mask2Filename);
	}
	else mask2Size = 0;

	// reading the matrices
	if (fileExists(theta1FileName))
		readTheta1Matrix(theta1FileName);

	if (fileExists(theta2FileName))
		readTheta2Matrix(theta2FileName);
	normalization = 1;
	return 0;
}

void logisticRegressionProcessing::readMask1(char *fileName)
{
	read_volume(mask1, string(fileName));
	mask1Size = MaskSize(mask1);
}

void logisticRegressionProcessing::readMask2(char *fileName)
{
	read_volume(mask2, string(fileName));
	mask2Size = MaskSize(mask2);
}

void logisticRegressionProcessing::readTheta1Matrix(char *fileName)
{
	theta1Matrix = read_ascii_matrix(fileName);
}

void logisticRegressionProcessing::readTheta2Matrix(char *fileName)
{
	theta2Matrix = read_ascii_matrix(fileName);
}

// calculates the feedback value
int logisticRegressionProcessing::processVolume(volume<float> &v, int index, float &feedbackValue)
{
	Matrix sample1, sample2;
	Volume2Matrix(v, mask1, mask1Size, sample1);
	Volume2Matrix(v, mask2, mask2Size, sample2);

	// including the actual column vector to the last value matrix so the values dont go beyond [0, 1] range
	if (lastM1Values.Nrows() == 0)
	{
		lastM1Values = sample1.Column(1);
		lastM2Values = sample2.Column(1);
	}
	else
	{
		lastM1Values |= sample1.Column(1);
		lastM2Values |= sample2.Column(1);
	}

	if (normalization)
	{
		NormalizeVector(lastM1Values, sample1);
		NormalizeVector(lastM2Values, sample2);
	}

	sample1(1, 1) = 1;
	sample2(1, 1) = 1;

	Matrix resp1 = sigmoid(theta1Matrix *sample1);
	Matrix resp2 = sigmoid(theta2Matrix * sample2);

	int quadrant = bestIndex(resp1);
	int eccentricity = bestIndex(resp2);

	feedbackValue = quadrant * 10 + eccentricity;
	return 0;
}

// calculates the feedback value
int logisticRegressionProcessing::processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue)
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

	return 0;
}

// just calculate the sigmoid function for a value
double logisticRegressionProcessing::sigmoid(double value)
{
	return 1.0 / (1.0 + exp(-value));
}

// just calculate the sigmoid function for a Matrix
Matrix logisticRegressionProcessing::sigmoid(Matrix value)
{
	Matrix resp;
	resp.ReSize(value.Nrows(), value.Ncols());
	for (int i = 1; i <= value.Ncols(); i++)
	{
		for (int j = 1; j <= value.Nrows(); j++)
		{
			resp(j, i) = sigmoid(value(j, i));
		}
	}
	return resp;
}

// return the row index of a column vector with the max value
int logisticRegressionProcessing::bestIndex(Matrix &value)
{
	int idx = 0;
	float bestValue = -100.0;

	for (int t = 1; t <= value.Nrows(); t++)
	{
		if (bestValue < value(t, 1))
		{
			idx = t;
			bestValue = value(t, 1);
		}
	}
	return idx;
}

// plugin function for initializating the roi processing object
DLLExport initializeProcessing(studyParams &vdb, void *&userData)
{
	logisticRegressionProcessing *roiVar = new logisticRegressionProcessing;
	roiVar->initialization(vdb);
	userData = roiVar;
	return 0;
}

// plugin function for finalizing the object
DLLExport finalizeProcessing(studyParams &vdb, void *&userData)
{
	logisticRegressionProcessing *roiVar = (logisticRegressionProcessing *)userData;
	delete roiVar;
	userData = NULL;
	return 0;
}

// plugin function for calculationg feedback value
DLLExport processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue, void * &userData)
{
	logisticRegressionProcessing *roiVar = (logisticRegressionProcessing *)userData;
	roiVar->processVolume(vdb, index, classnum, feedbackValue);
	return 0;
}
