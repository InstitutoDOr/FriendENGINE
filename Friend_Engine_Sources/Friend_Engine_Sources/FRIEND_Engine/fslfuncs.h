#ifndef FSLFUNCS_H
#define FSLFUNCS_H
#include "fslio/fslio.h"

// functions handling FSL commands
#ifdef WINDOWS

extern "C" int _stdcall flirt(char *CmdLn);
extern "C" int _stdcall bet(char *CmdLn);
extern "C" int _stdcall tsplot(char *CmdLn);
extern "C" int _stdcall fslmaths(char *CmdLn);
extern "C" int _stdcall mcflirt(char *CmdLn);
extern "C" int _stdcall fsl_tsplot(char *CmdLn);
extern "C" int _stdcall feat_model(char *CmdLn);
extern "C" int _stdcall fsl_glm(char *CmdLn);
extern "C" int _stdcall convert_xfm(char *CmdLn);

#else
int callCommand(char *command);
int fslmaths(char *cmdLn);
int fslroi(char *cmdLn);
int susan(char *cmdLn);
int bet(char *cmdLn);
int flirt(char *cmdLn);
int mcflirt(char *cmdLn);
int fsl_glm(char *cmdLn);
int fsl_tsplot(char *cmdLn);
int feat_model(char *cmdLn);
int convert_xfm(char *cmdLn);
#endif

// verifies if a file is readable by FSL
bool isFSLReadable(char *fileName);

// reorients a volume file to axial
int axial(char *innam, char *outnam);

// accomodate a volume `mask` to the same coordinate space of a volume `base`
int sameFov(char *base, char *mask, char *output);

// reformats a `mask` volume to axial and harmonizes the voxels geometry with `base` volume.
int standardizeVolume(char *base, char *mask, char *output, int NN);

// calculates the final step of the FRIEND pipeline, the sliding window mean
int estimateActivation(int ini, int end, int slidingWindowSize, char *suffix, char *output);
int estimateActivation(int ini, int end, int slidingWindowSize, char *suffix, char *maskFileName, char *output);

// loads the data struture of a volume in memory
FSLIO * fslioopen(char *file);

// filters a volume based on a min cluster size and a minimum value
void clusterSizeFiltering(char *fileName, char *outputName, int minClusterSize, float minValue, int connectionType=6);

// frees the allocated memory
void fslioclose(FSLIO *src);

// resample `volume` file to the same voxel dimensions that `reference`
void equalVoxelDim(char *volume, char *reference, char *output, float TR, int nn);

// same functionality of the fslswapdim command
int fslSwapDimRT(const char *CmdLn, FSLIO *src);

// changes the voxel dimension of a volume. Depends of resample function
int resampleVolume(char *anatomic, char *sAnatomic, float dx, float dy, float dz, float TR, int nn);

// centralizes a volume in FOV
int centralizeVolume(char *Anatomic, char *sAnatomic);

// transforms `mniTemplate` volume in MNI space to native space of `betRFI`
void MniToSubject(char *BetRFI, char *mniTemplate, char *mniStandard, char* RFI2MNI, char* RFI2MNITransf, char* MNI2RFITransf, char* output, char *prefix = NULL);

// transforms `mniTemplate` volume in MNI space to native space of `betRFI`
void MniToSubject(char *BetRFI, char *mniTemplate, char *mniStandard, char* output, char *prefix = NULL);

// this function adjusts a `mask` volume corregistered with a `functional` with another `reference` volume. This function is used to adjust a mask coregistered with the RFI volume with the first volume of a run, to account for movements.
void functionalNormalization(char *mask, char *functional, char *reference, char *output, bool nearestNeighbour);

// this function engraves a roi volume in a RFI volume, to make sure of the side
void uniteVolumes(char *referenceVolume, char *roiVolume, char *outputFile);

#endif
