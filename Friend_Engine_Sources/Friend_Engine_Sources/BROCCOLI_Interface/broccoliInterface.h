#ifndef ____broccoliInterface__
#define ____broccoliInterface__

#include "broccoli_lib.h"
#include <opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include "nifti1_io.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "newimage/newimageall.h"
#include "fslio/fslio.h"
#include "fslfuncs.h"
#include "filefuncs.h"
#include <sys/timeb.h>

using namespace NEWIMAGE;

typedef struct
{
	float *h_Quadrature_Filter_1_Linear_Registration_Real, *h_Quadrature_Filter_1_Linear_Registration_Imag,
		  *h_Quadrature_Filter_2_Linear_Registration_Real, *h_Quadrature_Filter_2_Linear_Registration_Imag,
		  *h_Quadrature_Filter_3_Linear_Registration_Real, *h_Quadrature_Filter_3_Linear_Registration_Imag;
} broccoliLinearRegistrationFilters;


void printBufferErrors(BROCCOLI_LIB &BROCCOLI);
void allocateLinearRegistrationFilter(broccoliLinearRegistrationFilters &filter);
void freeLinearRegistrationFilter(broccoliLinearRegistrationFilters &filter);
void prepareLinearRegistrationFilters(broccoliLinearRegistrationFilters &filter);
void printRegistrationMatrix(float *h_Registration_Parameters);
void prepareBROCCOLIRegistration(BROCCOLI_LIB &BROCCOLI, volume4D<float> &T1, volume<float> &mniVolume, volume4D<float> &output);
void prepareBROCCOLIMotionCorrection(BROCCOLI_LIB &BROCCOLI, volume<float> &ref, float *&h_Reference_Volume, float *&h_Motion_Parameters, float *&fmriVolume, broccoliLinearRegistrationFilters &filter);
void vector2volume(volume4D<float> &volume, float *volumeVector);
int64_t volume2vector(volume4D<float> &volume, float *&volumeVector);
int64_t volume2vector(volume<float> &volume, float *&volumeVector);
void copyVolume2vector(volume<float> &volume, float *&volumeVector);
void BROCCOLILinearRegistration(char *T1file, char *MNIfile, char *outputFile);
void BROCCOLIMotionCorrection(char *refFile, char *prefix, char *motionparms, int numFMRIs);
void prepareBROCCOLISmoothing(BROCCOLI_LIB &BROCCOLI, volume<float> &inputData, float fwhm, float *&fMRI_Volume, float *&h_Certainty);
double GetWallTime();


class BROCCOLIEngineInterface
{
   public:
	   int64_t allocatedHostMemory;

	   char rfiName[500];
	   volume<float> rfiVolume;
	   float *epiVolume, *h_Certainty;
	   BROCCOLI_LIB *BROCCOLI;
	   float fwhmValue;

	   broccoliLinearRegistrationFilters filter;
	   float *h_Reference_Volume, h_Motion_Parameters[12];

	   void setRFI(char *rfiFileName);
	   void createBROCCOLIObject(int platform, int device);
	   void deallocateBROCCOLIObject();

	   void prepareRegistrationFilters();

	   void prepareSmoothing(float fwhm);
	   void runSmoothing(volume<float> &epi);

	   void prepareMotionCorrection();
	   void runMotionCorrection(volume<float> &epi);

	   void preparePipeline(float fwhm);

	   void pipelineEngine(char *epiFileName, char *outputFileName);
	   void pipelineNormalEngine(char *epiFileName, char *outputFileName);
	   void pipelineEngine(volume<float> &epi);

	   BROCCOLIEngineInterface() 
	   {  
		   //createBROCCOLIObject();
	   };

	   ~BROCCOLIEngineInterface()
	   {
		   //deallocateBROCCOLIObject();
	   };
};

#endif