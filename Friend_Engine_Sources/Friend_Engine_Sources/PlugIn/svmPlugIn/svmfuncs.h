#ifndef SVMFUNC_H
#define SVMFUNC_H
#include "svm.h"
#include <vector>

struct weightPoint {
	int index;
	double value;
};

using namespace std;

// calculates the svm weights of a model, storing the values in a vector object
void calculateWeightVector(const svm_model *model, vector <double> &weightVector);

// get a svm score based on a vector sample. this functions returns the prediction class as a return value.
double svmSamplePredict(const svm_model *model, const svm_node *sample, double &score);

// loads a model file in memory
svm_model* loadModel(char *modelFileName);

// unloads a model in memory
void unloadModel(svm_model *model);

// transforms a 4D Volume in a SVM samples file, based on a mask
void saveSVMFile(const char *volume4DFileName, const char *maskFileName, const char *outputFileName, float minValue, vector <int > &indexes, vector <int> &classes);

// function to return the prediction class and svm score of a volume.
float predict(svm_model *model, const char *volumeFileName, const char *maskFileName, float &predictedClass, float &svmScore);

// the following functions generates a weight volume from a svm model
void generateWeightVolume(svm_model *model, const char *maskFileName, int normalize, char *outputFileName);

void generateWeightVolume(svm_model *model, const char *maskFileName, int vectorSize, float minValue, int normalize, char *outputFileName);


#endif
