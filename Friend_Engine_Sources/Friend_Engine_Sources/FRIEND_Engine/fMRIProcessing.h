#ifndef ___fMRIPROCESSING_H
#define ___fMRIPROCESSING_H

#include "fslio/fslio.h"
#include "fslfuncs.h"
#include "filefuncs.h"
#include "direct.h"
#include "newimage/newimageall.h"
#include "masks.h"

using namespace NEWIMAGE;

void FslFree(FSLIO* OP);
extern "C" int _stdcall melodic(char *CmdLn);
extern "C" int _stdcall susan(char *CmdLn);
extern "C" int _stdcall fast(char *CmdLn);
extern "C" int _stdcall fslroi(char *CmdLn);
extern "C" int _stdcall slicer(char *CmdLn);
extern "C" int _stdcall pngappend(char *CmdLn);
extern "C" int _stdcall slicetimer(char *CmdLn);

void regressOutData(volume4D<float> &data, char *regressorFile);
void regressOutDesign(char *inFile, char *designFile, char *outFile);
void regressOutDesign(char *inFile, char *outFile, int initialColumn = 2);
void regressOutFromDataTrends(char *iname, char *regressorFile, char *oname, int NOrdem = 0);
void regressOut(Matrix &regressors, Matrix&data, Matrix &residual);
int applyTransformation(char *volume, char *output, char *reference, char *matrixFile, int nn = 0);
int mniRegistration(char *volume, char *highres, char *reference, char *pngFile, char *matrixFile);
int icaDenoisedSmart(char *iname, char *oname, char *matrix, char *fixIdxs, Matrix& design, int aggressive = 0);
int motionCorrection4D(char *volume4D, char *reference, char *regressorFile, char *mcsuffix = NULL, char *filesuffix = NULL);
void detrendSG(char *filename, char *saida, int w, int degree);
void runMelodic(char *dirSuj, char *dirRun);
void runMelodic4D(char *inName);
int processSubject(char *projectdir, char *subject, char *standard, vector<string>&runs, int window = 0);
void create4Dvolume(char *templateFile, char *output);
void estimateActivation4D(volume4D<float> &data, int slidingWindow);
void applySusan(char *entrada, char *saida, char *mascara, float FWHM = 5.0);
void filtroGaussiano(char *entrada, char *saida, float FWHM = 5.0);
void fsl_prepare_fieldmap(char *phaseImage, char *magImage, char *outputImg, float deltaTE=2.46);
int maskSize(NEWIMAGE::volume<float> &mask);
void getFieldMapInfo(char *fileName, float &dwellTime, int &direction, float &TE);
void preparaDataMelodic(char *entrada, char *saida, float FWHM = 5.0, float highpass = 200.0, float lowpass = 0.0, float brain_thresh = 10.0, float TR = 2.0, int slicetimerType = 0);
void prepareFieldMap(char *reference, char *phaseImage, char *magImage, char *fieldMap, float deltaTE=2.46);
void applyFieldMap(char *fileName, char *fieldMap, char *output, float dwellTime, int direction);
void applyFieldMap(char *fileName, char *fieldMap, char *output);
#endif