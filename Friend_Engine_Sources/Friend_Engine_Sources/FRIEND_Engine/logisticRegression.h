//
//  logisticRegression.h
//  
//
//  Created by IDOR on 06/06/13.
//
//

#ifndef ____logisticRegression__
#define ____logisticRegression__

#include <iostream>

#include "newimage/newimageall.h"
#include "fslfuncs.h"
#include "vardb.h"
#include "masks.h"
#include <sstream>

using namespace NEWIMAGE;

// class responsible for roi feedback processing, using percent signal change and a target value to calculate the termometer value
class logisticRegressionProcessing
{
	Matrix theta1Matrix, theta2Matrix;
	Matrix lastM1Values, lastM2Values;
	volume<float>mask1, mask2;
	int mask1Size, mask2Size;

public:

	int normalization;
	// initializes the object variables
	int initialization(studyParams &vdb);
	// calculates the feedback value
	int processVolume(studyParams &vdb, int index, float &classnum, float &feedbackValue);
	// calculates the feedback value
	int processVolume(volume<float> &v, int index, float &feedbackValue);
	// just calculate the sigmoid function for a value
	double sigmoid(double value);
	// just calculate the sigmoid function for a Matrix
	Matrix sigmoid(Matrix value);
	// return the row index of a column vector with the max value
	int bestIndex(Matrix &value);

	void readMask1(char *fileName);
	void readMask2(char *fileName);
	void readTheta1Matrix(char *fileName);
	void readTheta2Matrix(char *fileName);

};

#endif /* defined(____logisticRegression__) */

