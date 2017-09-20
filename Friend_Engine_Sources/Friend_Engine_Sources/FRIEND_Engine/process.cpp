#include "newimage/newimageall.h"
#include "fslfuncs.h"
#include "process.h"
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include "filefuncs.h"
#include "utils.h"
#include "masks.h"
#include "dcm2niiInterface.h"

using namespace std;
using namespace NEWIMAGE;

#ifdef WINDOWS
#include <direct.h>
#define snprintf sprintf_s
#define chdir _chdir
char extension[]="";
#else
char extension[] = "";

#endif

// set socket data for response
void FriendProcess::setSocketfd(int Sock)
{
	vdb.setSocketfd(Sock);
}

// config file functions
void FriendProcess::readConfigFile(char *configFile)
{
   vdb.readConfigFile(configFile);
}

void FriendProcess::readConfigBuffer(char *buffer, int size)
{
   vdb.readConfigBuffer(buffer, size);
}

void FriendProcess::readConfig(CSimpleIniA &ini)
{
   vdb.readConfig(ini);
}

void FriendProcess::saveConfigBuffer(char *buffer, int size, char *configFile)
{
   vdb.saveConfigBuffer(buffer, size, configFile);
}


int isGoodDesign(char *fileName)
{
	DiagonalMatrix eigenvals(1);
	Matrix real_X;
	int result = 0;

	real_X = read_ascii_matrix(string(fileName));
	if (real_X.Ncols()>1)
	{
		SVD(real_X, eigenvals);

		SortDescending(eigenvals);

		float inv_condition = eigenvals.Minimum() / eigenvals.Maximum();

		if (eigenvals.Minimum()<1e-16)
		{
			fprintf(stderr, "Invalid design matrix.\n");
		}
		else result = 1;
	}
	else if (real_X.Ncols() == 1) result = 1;

	eigenvals.Release();
	return result;
}

// prepare files and run fsl_glm on processed volumes
void FriendProcess::glm()
{
   char prefix[BUFF_SIZE], Pref[BUFF_SIZE];
   char fsfFile[BUFF_SIZE];
   char auxString[BUFF_SIZE];
   stringstream CmdLn;

   // confirming that all variables are set
   if (!vdb.rPrepVars) prepRealtimeVars();

   vdb.logObject->writeLog(1, "Generating 4D File\n");
   vdb.getFinalVolumeFormat(prefix);
   vdb.getPreprocVolumePrefix(Pref);
   
   // generating 4D FRIEND pipeline volume for the run. Here we do not performs the final step, sliding window mean
   // as this calculation shifts the voxel time series
   estimateActivation(1, vdb.interval.maxIndex(), 1, prefix, vdb.trainGLM4DFile);

   // demeaning the 4D FRIEND pipeline volume for the run. 
   CmdLn.str("");
   CmdLn << "fslmaths " << vdb.trainGLM4DFile << " -Tmean -mul -1 -add " << vdb.trainGLM4DFile << " " << vdb.trainGLM4DFile;

   fslmaths((char *)CmdLn.str().c_str());


   // concatenating the .par (movements parameters) files
   generateConfoundFile(Pref, 1, vdb.interval.maxIndex(), vdb.parFile);
   
   // concatenating the .rms files
   generateRmsFile(Pref, 1, vdb.interval.maxIndex(), vdb.rmsFile);

   // generating the rotation graph png
   vdb.logObject->writeLog(1, "Generating rotation movement graphics\n");
   CmdLn.str("");
   CmdLn << "fsl_tsplot -i " << vdb.parFile << " -t \"MCFLIRT estimated rotations (radians)\" -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o " << vdb.logDir << "rot" << vdb.trainFeatureSuffix << ".png";

   fsl_tsplot((char *)CmdLn.str().c_str());

   // generating the translations graph png
   vdb.logObject->writeLog(1, "Generating translations movement graphics\n");
   CmdLn.str("");
   CmdLn << "fsl_tsplot -i " << vdb.parFile << " -t \"MCFLIRT estimated translations (mm)\" -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o " << vdb.logDir << "trans" << vdb.trainFeatureSuffix << ".png";

   fsl_tsplot((char *)CmdLn.str().c_str());

   // generating the rms graph png
   vdb.logObject->writeLog(1, "Generating rms graphics\n");
   CmdLn.str("");
   CmdLn << "fsl_tsplot -i " << vdb.rmsFile << " -t \"MCFLIRT estimated mean displacement (mm)\" -u 1 --start=1 --finish=1 -a absolute -w 640 -h 144 -o " << vdb.logDir << "rms" << vdb.trainFeatureSuffix << ".png";

   fsl_tsplot((char *)CmdLn.str().c_str());

   // copying the concatenated movements parameter file to log dir
   vdb.logObject->writeLog(1, "Copying confound files\n");
   CmdLn.str("");
   CmdLn << vdb.logDir << "confounds" <<  vdb.trainFeatureSuffix << ".txt";
   //CmdLn << "/bin/cp -p " << vdb.parFile << " " << vdb.logDir << "confounds" <<  vdb.trainFeatureSuffix << ".txt";
   //system(CmdLn.str().c_str());
   copyFile(vdb.parFile, CmdLn.str().c_str());

   
   // generating the design matrix file
   if (vdb.interval.conditionNames.size() > 0)
   {
	  vdb.logObject->writeLog(1, "Generating conditions files.\n");
	  vdb.interval.generateConditionsBoxCar(vdb.glmDir);
	  vdb.logObject->writeLog(1, "Generating FSF file.\n");

	  strcpy(vdb.interval.glmDir, vdb.glmDir);
	  vdb.interval.generateFSFFile(vdb.fsfFile, vdb.runSize, vdb.includeMotionParameters, 0);

	  vdb.logObject->writeLog(1, "Running feat_model.\n");
	  CmdLn.str("");
	  CmdLn << "feat_model " << vdb.glmDir << vdb.subject << vdb.trainFeatureSuffix;
	  if (vdb.includeMotionParameters) CmdLn << " " << vdb.parFile;
	  chdir(vdb.glmDir);
	  feat_model((char *)CmdLn.str().c_str());

	  remove(vdb.contrastFile);
	  sprintf(auxString, "%s%s%s%s", vdb.glmDir, vdb.subject, vdb.trainFeatureSuffix, ".con");
	  rename(auxString, vdb.contrastFile);

	  remove(vdb.glmMatrixFile);
	  sprintf(auxString, "%s%s%s%s", vdb.glmDir, vdb.subject, vdb.trainFeatureSuffix, ".mat");
	  rename(auxString, vdb.glmMatrixFile);
   }

   sprintf(auxString, "%s%s", vdb.inputDir, "RFI_binmask.nii");

   vdb.logObject->writeLog(1, "Making binary mask.\n");
   CmdLn.str("");
   CmdLn << "fslmaths " << vdb.maskFile << " -bin " << auxString << " -odt char";
   fslmaths((char *)CmdLn.str().c_str());

   vdb.logObject->writeLog(1, "Running fsl_glm\n");
   if (isGoodDesign(vdb.glmMatrixFile))
   {
	   CmdLn.str("");
	   CmdLn << "fsl_glm -i " << vdb.trainGLM4DFile << " -d " << vdb.glmMatrixFile << " -c " << vdb.contrastFile << " -m " << auxString << " -o " << vdb.glmDir << "betas" << vdb.trainFeatureSuffix << " --out_t=" << vdb.glmTOutput << " --out_z=" << vdb.glmZOutput << " --out_p=" << vdb.glmDir << "pvalues" << vdb.trainFeatureSuffix;
	   fsl_glm((char *)CmdLn.str().c_str());
   }
   else
   {
	   vdb.logObject->writeLog(1, "Invalid design matrix.\n");
	   vdb.logObject->writeLog(1, "Problems in design matrix. Error in GLM.\n");
   }
   
   vdb.rGLM=true;
}

// performs the feature selection step
void FriendProcess::featureSelection()
{
   stringstream CmdLn;
   if (!vdb.rPrepVars) prepRealtimeVars();
   if (fileExists(vdb.glmTOutput))
   {
	   // generate file with the maximum T for each voxel in the contrasts made
	   char noInterestConditions[500];
	   noInterestConditions[0]=0;

	   if (vdb.readedIni.GetValue("FRIEND", "noInterestConditions"))
	      strcpy(noInterestConditions, vdb.readedIni.GetValue("FRIEND", "noInterestConditions"));

	   // if there is no interest conditions, generates the Tmax volume special way
	   if (strlen(noInterestConditions) > 0)
	   {
		   vector<int> idxs;
		   vdb.interval.noInterestContrastsIndexes(noInterestConditions, idxs);
		   generateTMaxVoxels(vdb.glmTOutput, vdb.featuresAllTrainSuffix, idxs);
	   }
	   else
	   {
		   vdb.logObject->writeLog(1, "Generating the all contrasts colapsed mask.\n");
		   CmdLn << "fslmaths " << vdb.glmTOutput << " -Tmax " << vdb.featuresAllTrainSuffix;
		   fslmaths((char *)CmdLn.str().c_str());
	   }
   }
   else
   {
	   // reporting error an creating a fake glm t output
	   vdb.logObject->writeLog(1, "Error occurred in glm calculation. Feature selection derived from it is invalid.\n");

	   CmdLn << "fslmaths " << vdb.maskFile << " -mul 0 " << vdb.featuresAllTrainSuffix;
	   fslmaths((char *)CmdLn.str().c_str());

	   CmdLn.str("");
	   CmdLn << "fslmaths " << vdb.maskFile << " -bin " << vdb.glmDir << "pvalues" << vdb.trainFeatureSuffix;
	   fslmaths((char *)CmdLn.str().c_str());
   }

   // generates the mni mask in subject space
   if (fileExists(vdb.mniMask) && fileExists(vdb.mniTemplate))
   {
	  char name[BUFF_SIZE], prefix[30]="_RFI2";
	  
	  vdb.logObject->writeLog(1, "Bringing MNI mask to native space.\n");
	  extractFileName(vdb.mniMask, name);
	  for (int t=0;t<strlen(name);t++)
	  if (name[t] == '.') name[t] = '_';
	  
	  sprintf(vdb.subjectSpaceMask, "%s%s%s.nii", vdb.inputDir, name, vdb.trainFeatureSuffix);
	  
	  // brings the mni mask to subject space
	  MniToSubject(vdb.maskFile, vdb.mniMask, vdb.mniTemplate, vdb.subjectSpaceMask, prefix);
   }
   else if (fileExists(vdb.subjectSpaceMaskUser))
   {
	   char outputFile[BUFF_SIZE];
	   sprintf(outputFile, "%s.nii", vdb.featuresTrainSuffix);
	   vdb.logObject->writeLog(1, "Using all subject mask, defined by user, by copying the file %s to %s.\n", vdb.subjectSpaceMaskUser, outputFile);

	   CmdLn.str("");
	   CmdLn << "fslmaths " << vdb.subjectSpaceMaskUser << " " << outputFile;
	   fslmaths((char *)CmdLn.str().c_str());
   }
   
   if ((vdb.useWholeSubjectSpaceMask) && (fileExists(vdb.subjectSpaceMask)))
   { // just use all mni mask (in native space)
	  char outputFile[BUFF_SIZE];
	  sprintf(outputFile, "%s.nii", vdb.featuresTrainSuffix);
	  vdb.logObject->writeLog(1, "Using all subject mask by copying the file %s to %s.\n", vdb.subjectSpaceMask, outputFile);

	  CmdLn.str("");
	  CmdLn << "fslmaths " << vdb.subjectSpaceMask << " " << outputFile;
	  fslmaths((char *)CmdLn.str().c_str());
   }
   else if (!fileExists(vdb.subjectSpaceMaskUser))
   {
	   switch (vdb.thresholdType)
	   {
		   case 0: // threshold by percentage
		   {
			   /// select a percentage of higher voxels
			   vdb.logObject->writeLog(1, "Selecting the best GLM voxels by percentage.\n");
			   volume<float> features;
			   string featuresFile = vdb.featuresAllSuffix;
			   featuresFile += vdb.trainFeatureSuffix;
			   read_volume(features, featuresFile);
			   float value = features.percentile((double)1.0 - (double)(vdb.percentileHigherVoxels / 100.0));

			   // do threshold
			   CmdLn.str("");
			   CmdLn << "fslmaths " << vdb.featuresAllTrainSuffix << " -thr " << value << " " << vdb.featuresTrainSuffix;
			   fslmaths((char *)CmdLn.str().c_str());

			   break;
		   }

		   case 1: { // threshold by t-value
			   // do threshold
			   vdb.logObject->writeLog(1, "Selecting the best GLM voxels by value.\n");
			   CmdLn.str("");
			   CmdLn << "fslmaths " << vdb.featuresAllTrainSuffix << " -thr " << vdb.tTestCutOff << " " << vdb.featuresTrainSuffix;
			   fslmaths((char *)CmdLn.str().c_str());
			   break;
		   }

		   case 2: { // threshold by p-value
			   // do threshold
			   vdb.logObject->writeLog(1, "Selecting the best GLM voxels by p-value.\n");
			   CmdLn.str("");
			   CmdLn << "fslmaths " << vdb.glmDir << "pvalues" << vdb.trainFeatureSuffix << " -uthr " << vdb.pvalueCutOff << " -bin -mul " << vdb.featuresAllTrainSuffix << " " << vdb.featuresTrainSuffix;
			   fslmaths((char *)CmdLn.str().c_str());
			   break;
		   }
	   }

	   // Saving a copy of the unmasked features
	   CmdLn.str("");
	   CmdLn << "fslmaths " << vdb.featuresTrainSuffix << " " << vdb.featuresTrainSuffix << "_NoMasked";
	   fslmaths((char *)CmdLn.str().c_str());

	   // creates the intersection between the subjectmask and the thresholded mask
	   if (fileExists(vdb.subjectSpaceMask))
	   {
		   CmdLn.str("");
		   CmdLn << "fslmaths " << vdb.featuresTrainSuffix << " -mas " << vdb.subjectSpaceMask << " " << vdb.featuresTrainSuffix;
		   fslmaths((char *)CmdLn.str().c_str());
	   }

	   // erasing clusters with less than clustersize voxels
	   if (vdb.clusterSize > 0)
		   clusterSizeFiltering(vdb.featuresTrainSuffix, vdb.featuresTrainSuffix, vdb.clusterSize, 0);
   }

   // restrict the voxel mask created to the betted rfi 
   CmdLn.str("");
   CmdLn << "fslmaths " << vdb.featuresTrainSuffix << " -mas " << vdb.maskFile << " " << vdb.featuresTrainSuffix;
   fslmaths((char *)CmdLn.str().c_str());

   // create binary version
   CmdLn.str("");
   CmdLn << "fslmaths " << vdb.featuresTrainSuffix << " -bin " << vdb.featuresTrainSuffix << "_bin -odt char";
   fslmaths((char *)CmdLn.str().c_str());

   vdb.rFeatureSel = true;
}

// call plug-in train function
void FriendProcess::train()
{
   if (!vdb.rPrepVars) prepRealtimeVars();
   if (pHandler.callTrainFunction(vdb))
	  vdb.rTrain=true;
}

// call plug-in test function
void FriendProcess::test(int index, float &classNum, float &projection)
{
   if (!vdb.rPrepVars) prepRealtimeVars();
   pHandler.callTestFunction(vdb, index, classNum, projection);
}

// calculates the mean volume of an interval
// just exec : fslmaths volume1 -add volume2 ... volumeN -div n
void FriendProcess::baselineCalculation(int intervalIndex, char *baseline)
{
   char Pref[BUFF_SIZE];
   char suf[20]= "_mc.nii";
   char format[50], number[50];
   
   // facilitates the filename generation. Generates the number string with the correct zeros
   sprintf(format, "%c0%dd", '%', vdb.numberWidth);
   char *prefix = extractFileName(vdb.rawVolumePrefix);
   snprintf(Pref, BUFF_SIZE-1, "%s%s", vdb.preprocDir, prefix);
   free(prefix);

   stringstream osc;
   osc  << "fslmaths ";

   // iterating to add the filenames to command
   for(int t=vdb.averageMeanOffset; t<= (vdb.interval.intervals[intervalIndex].end-vdb.interval.intervals[intervalIndex].start); t++)
   {
	  sprintf(number, format,(vdb.interval.intervals[intervalIndex].start + t));
	  // if the second file, adds the `-add` token
	  if (t > vdb.averageMeanOffset) osc << " -add ";
	  // adds the file in command line
	  osc << " " << Pref << number << suf;
   };

   // generates the -div N part
   sprintf(number, "%d",(vdb.interval.intervals[intervalIndex].end-vdb.interval.intervals[intervalIndex].start+1-vdb.averageMeanOffset));
   osc << " -div " << number << " " << baseline << '\0';

   // run fslmaths
   fslmaths((char *)osc.str().c_str());
};

// verifies if the condition is a baseline condition
BOOL FriendProcess::isBaselineCondition(char * condition)
{
   return vdb.interval.isBaselineCondition(condition);
}

// verifies if the next file is ready to be read by FRIEND and transforms it in nifti accordingly
BOOL FriendProcess::isReadyNextFile(int index, char *rtPrefix, char *format, char *inFile)
{
	BOOL result;
	for (int i = 0; i <= vdb.volumesToSkip; i++)
	{
		result = isReadyNextFileCore(index + i, index, rtPrefix, format, inFile);
		if (result) break;
	}
	return result;
}

// verifies if the next file is ready to be read by FRIEND and transforms it in nifti accordingly
BOOL FriendProcess::isReadyNextFileCore(int indexIn, int indexOut, char *rtPrefix, char *format, char *inFile)
{
	BOOL response = 0, fileChecked = 0, fileFound = 0;
	char  numberIn[50], numberOut[50], outFile[BUFF_SIZE], volumeName[BUFF_SIZE], dcm2niiExe[BUFF_SIZE];

	// forming the output file
	sprintf(numberOut, format, indexOut);
	sprintf(outFile, "%s%s%s", rtPrefix, numberOut, ".nii");

#ifdef WINDOWS
	sprintf(dcm2niiExe, "%s%cdcm2nii.exe", exePath, PATHSEPCHAR);
#else
	sprintf(dcm2niiExe, "%s%cdcm2nii", exePath, PATHSEPCHAR);
#endif
	// forming the input file
	sprintf(numberIn, format, indexIn);

	if (!fileExists(outFile))
	{
		// in DICOM
		sprintf(inFile, "%s%s%s", vdb.rawVolumePrefix, numberIn, ".dcm");
		//if (!fileFound) vdb.logObject->writeLog(1, "Searching file : %s\n", inFile);
#ifdef dcm2niifunction
		if (fileExists(inFile))
		{
			if (isReadableSize(inFile, rifSize))
			{
				try
				{
					transformDicom(inFile, vdb.preprocDir);
					fileFound = 1;
				}
				catch (...)
				{
					fileFound = 0;
				};
			}
		}
#else
		if (fileExists(inFile) && fileExists(dcm2niiExe))
		{
			if (isReadableSize(inFile, rfiSize))
			{
				// executes the dcm2nii tool
				stringstream osc;
				osc << dcm2niiExe << " -b " << exePath << PATHSEPCHAR << "dcm2nii.ini -o " << vdb.preprocDir << " " << inFile;
				vdb.logObject->writeLog(1, "Executting dcm2nii : %s\n", osc.str().c_str());
				system(osc.str().c_str());
				fileFound = 1;
		    }
		}
#endif
		else
		{
			// in PAR/REC
			sprintf(inFile, "%s%s%s", vdb.rawVolumePrefix, numberIn, ".par");
			if (fileExists(inFile) && fileExists(dcm2niiExe))
			{
				// executes the dcm2nii tool
				stringstream osc;
				osc << dcm2niiExe << " -b " << exePath << PATHSEPCHAR << "dcm2nii.ini -o " << vdb.preprocDir << " " << inFile;
				vdb.logObject->writeLog(1, "Executting dcm2nii : %s\n", osc.str().c_str());
				system(osc.str().c_str());

				sprintf(inFile, "%s%s%s", rtPrefix, numberIn, ".nii");
				osc.str("");
				osc << "fslswapdim " << inFile << " x y z " << outFile << '\0';

				fslSwapDimRT(osc.str().c_str(), vdb.runReferencePtr);
				fileFound = 1;
			}

			// in Analyze. The engine converts it to nifti and inverts the axis, if needed
			sprintf(inFile, "%s%s%s", vdb.rawVolumePrefix, numberIn, ".img");
			//if (!fileFound) vdb.logObject->writeLog(1, "Searching file : %s\n", inFile);
			if (fileExists(inFile))
			{
				if (isReadable(inFile))
				{
					stringstream osc;
					osc << "fslswapdim " << inFile;

					if (!vdb.invX) osc << " x";
					else osc << " -x";

					if (!vdb.invY) osc << " y";
					else osc << " -y";

					if (!vdb.invZ) osc << " z";
					else osc << " -z";

					osc << "  " << outFile << '\0';
					// transforms an analyze file into a nii file, inverting axes accordingly
					fslSwapDimRT(osc.str().c_str(), vdb.runReferencePtr);

					osc.str("");
					osc << "fslmaths " << outFile << " " << outFile << " -odt float";
					fslmaths((char *)osc.str().c_str());
					fileFound = 1;
				}
			}

			// in NIFTI. 
			sprintf(inFile, "%s%s%s", vdb.rawVolumePrefix, numberIn, ".nii");
			//if (!fileFound) vdb.logObject->writeLog(1, "Searching file : %s\n", inFile);
			if (fileExists(inFile))
			{
				if (isReadable(inFile))
				{
					stringstream osc;
					osc.str("");
					osc << "fslmaths " << inFile << " " << outFile << " -odt float";
					fslmaths((char *)osc.str().c_str());
					fileFound = 1;
				}
			}

			if (fileFound == 0)
			{
				// issuing a warning in terminal if passed TR seconds
				if (GetWallTime() - lastTimeCheck >= checkTimeThreshold)
				{
					vdb.logObject->writeLog(1, "file not found : %s\n", inFile);
					lastTimeCheck = GetWallTime();
				}
			}


		}
	}

	// verifying if the final file exists.
	if ((fileExists(outFile)) && (isReadable(outFile)))
		response = 1;
	else response = 0;

	if (response)
	{
		strcpy(volumeName, outFile);
		pHandler.callVolumeFunction(vdb, indexIn, volumeName);
	}

	return response;
}

/*
// verifies if the next file is ready to be read by FRIEND and transforms it in nifti accordingly
BOOL FriendProcess::isReadyNextFileCore(int indexIn, int indexOut, char *rtPrefix, char *format, char *inFile)
{
   BOOL response=0, fileChecked=0;
   char  numberIn[50], numberOut[50], outFile[BUFF_SIZE], volumeName[BUFF_SIZE];

   // forming the output file
   sprintf(numberOut, format, indexOut);
   sprintf(outFile, "%s%s%s", rtPrefix, numberOut, ".nii");

   // forming the input file
   sprintf(numberIn, format, indexIn);

   if (!fileExists(outFile))
   {
	   // in DICOM
	   sprintf(inFile, "%s%s%s", vdb.rawVolumePrefix, numberIn, ".dcm");
	   if (fileExists(inFile))
	   {
		   // executes the dcm2nii tool
		   stringstream osc;
		   osc << exePath << PATHSEPCHAR << "dcm2nii -b " << exePath << PATHSEPCHAR << "dcm2nii.ini -o " << vdb.preprocDir << " " << inFile;
		   vdb.logObject->writeLog(1, "Executting dcm2nii : %s\n", osc.str().c_str());
		   system(osc.str().c_str());
	   }
	   else
	   {
		   // in Analyze. The engine converts it to nifti and inerts the axis, if needed
		   sprintf(inFile, "%s%s%s", vdb.rawVolumePrefix, numberIn, ".img");
		   if (fileExists(inFile))
		   {
			   if (isFSLReadable(inFile))
			   {
				   stringstream osc;
				   osc << "fslswapdim " << inFile;

				   if (!vdb.invX) osc << " x";
				   else osc << " -x";

				   if (!vdb.invY) osc << " y";
				   else osc << " -y";

				   if (!vdb.invZ) osc << " z";
				   else osc << " -z";

				   osc << "  " << outFile << '\0';
				   // transforms an analyze file into a nii file, inverting axes accordingly
				   fslSwapDimRT(osc.str().c_str(), vdb.runReferencePtr);

				   osc.str("");
				   osc << "fslmaths " << outFile << " " << outFile << " -odt float";
				   fslmaths((char *)osc.str().c_str());
			   }
		   }
	   }
   }

   // verifying if the final file exists. It servers the NIFTI case
   if ((fileExists(outFile)) && (isFSLReadable(outFile)))
		response = 1;
   else response = 0;
   
   if (response)
   {
	  strcpy(volumeName, outFile);
	  pHandler.callVolumeFunction(vdb, indexIn, volumeName);
   }
	  
   return response;
}
*/

// generate a file concatenating all movement parameters of the processed volumes
void FriendProcess::generateConfoundFile(char *dPrefix, int ini, int end, char *output)
{
   char format[50], number[50], parFileName[BUFF_SIZE];
   sprintf(format, "%c0%dd", '%', vdb.numberWidth);
   fstream Output(output, fstream::in | fstream::out | fstream::trunc);
   for (int t=ini; t<=end;t++)
   {
	   sprintf(number, format, t);
	   sprintf(parFileName, "%s%s%s", dPrefix, number, "_mc.nii.par");
	   fstream parFile(parFileName, fstream::in);
	   
	   float rx, ry, rz, tx, ty, tz;
	   parFile >> rx >> ry >> rz >> tx >> ty >> tz;
	   Output << rx << " " << ry << " " << rz <<  " " << tx << " " << ty << " " << tz << '\n';
   }
   Output.close();
}

// generate a file concatenating all the processed volumes rms values
void FriendProcess::generateRmsFile(char *dPrefix, int ini, int end, char *output)
{
   char format[50], number[50], rmsFilename[BUFF_SIZE];
   sprintf(format, "%c0%dd", '%', vdb.numberWidth);
   fstream Output(output, fstream::in | fstream::out | fstream::trunc);
   for (int t=ini; t<=end;t++)
   {
	   sprintf(number, format, t);
	   sprintf(rmsFilename, "%s%s%s", dPrefix, number, "_mc.nii_abs.rms");
	   fstream rmsFile(rmsFilename, fstream::in);
	   
	   float rms;
	   rmsFile >> rms;
	   Output << rms << '\n';
   }
   Output.close();
}

// runs the real time pipeline at once. Calls the realtimePipelineStep for each volume
void FriendProcess::runRealtimePipeline()
{
	char preprocVolumePrefix[BUFF_SIZE], auxConfigFile[BUFF_SIZE], dcm2niiServerConfig[BUFF_SIZE];
	char format[30];
	char msg[50];
   
	rfiSize = fileSize(vdb.rfiFile);
	sprintf(auxConfigFile, "%sstudy_params%s.txt", vdb.outputDir, vdb.trainFeatureSuffix);
	vdb.readedIni.SaveFile(auxConfigFile);
	vdb.logObject->writeLog(1, "######################### processing pipeline begin ######################\n");

	char logName[500];
	sprintf(logName, "%s%c%s%s.txt", vdb.logDir, PATHSEPCHAR, "outputLog", vdb.trainFeatureSuffix);
	vdb.logObject->initializeLogFile(logName);

	if (vdb.interval.intervals.size())
	{
		if (!vdb.rPrepVars) prepRealtimeVars();
		if ((!fileExists(vdb.raiFile)) || (!fileExists(vdb.rfiFile))) return;
		vdb.actualBaseline[0] = 0;
		vdb.actualImg = 1;
		vdb.actualInterval = 0;
		vdb.getFormat(format);

		removeDirectory(vdb.preprocDir);
#ifdef WIN32
		_mkdir(vdb.preprocDir);
#else
		mkdir(vdb.preprocDir, 0777); // notice that 777 is different than 0777
#endif
		sprintf(dcm2niiServerConfig, "%s%c%s", exePath, PATHSEPCHAR, "dcm2niiserver.ini");
		vdb.readedIni.SaveFile(dcm2niiServerConfig);

		// getting the preprocessing volume prefix
		vdb.getPreprocVolumePrefix(preprocVolumePrefix);

		// setting the reference volume
		sprintf(vdb.motionRefVolume, "%s", vdb.rfiFile);

		// getting the FSLIO pointer of the reference volume for fslwapdimRT calls
		vdb.runReferencePtr = fslioopen(vdb.motionRefVolume);

		// reading the design file
		vdb.interval.readDesignFile(vdb.designFile);
		vdb.interval.generateContrasts(vdb.conditionContrasts, 0);

		// looping
		passes = 0;
		while (vdb.actualImg <= vdb.runSize)
		{
			realtimePipelineStep(preprocVolumePrefix, format, vdb.actualBaseline);
			// getting the termination status. If termination evoked, get out of here.
			if (vdb.sessionPointer)
			{
				int status = vdb.sessionPointer->getTerminateState();
				if (status) break;
			}
		}
		fslioclose(vdb.runReferencePtr);

		// issuing the end of run response
		if (vdb.sessionPointer == NULL)
		{
			sprintf(msg, "%s", "ENDGRAPH\n");
			vdb.socks.writeString(msg);
		}
		vdb.rPipeline = true;
	}
	else vdb.logObject->writeLog(1, "######################### Invalid design file    ######################\n");
	vdb.logObject->writeLog(1, "######################### processing pipeline end   ######################\n");
}


// send Graph params to FRONTEND process. Maybe latter changing to JSON
void FriendProcess::sendGraphParams(char *mcfile, char *number)
{
	fstream gfile;
   
	char parFile[BUFF_SIZE], rmsFile[BUFF_SIZE];
	sprintf(parFile,  "%s%s", mcfile, ".par");
	sprintf(rmsFile,  "%s%s", mcfile, "_abs.rms");
   
	float value;
	stringstream msg;
   
	gfile.open(parFile, fstream::in | fstream::out);
	msg << "GRAPHPARS;" << number << ";";
	gfile >> value;
	msg << value << ";";
   
	gfile >> value;
	msg << value << ";";
   
	gfile >> value;
	msg << value << ";";
   
	 gfile >> value;
	msg << value << ";";
   
	gfile >> value;
	msg << value << ";";
   
	gfile >> value;
	msg << value << ";";
   
	gfile.close();
	gfile.open(rmsFile, fstream::in | fstream::out);
   
	gfile >> value;
	msg << value << '\n';
	gfile.close();
   
	if (vdb.sessionPointer == NULL)
	{
	   vdb.logObject->writeLog(1, "sending : %s\n", msg.str().c_str());
	   vdb.socks.writeString(msg.str().c_str());
	}
	else
	{
	   vdb.sessionPointer->processGraphMessage(msg.str().c_str());
	}
}

// theoretically loading the default plugin library and functions
void FriendProcess::loadLibrary()
{
}

// unloading the plugin library
void FriendProcess::unLoadLibrary()
{
   pHandler.unLoadLibrary();
}

// set the library file and function names of the plugIn
void FriendProcess::loadFunctions(char *library, char *trainFunc, char *testFunc, char *initFunc, char *finalFunc, char *volumeFunc, char *afterPreProcFunc)
{
	pHandler.logObject = vdb.logObject;
	pHandler.loadFunctions(library, trainFunc, testFunc, initFunc, finalFunc, volumeFunc, afterPreProcFunc);
//   if (vdb.rPrepVars) pHandler.callInitFunction(vdb);
}

// set the path for plug in library file
void FriendProcess::setLibraryPath(char *path)
{
	pHandler.logObject = vdb.logObject;
	pHandler.setLibraryPath(path);
}

// set the path for plug in library file
void FriendProcess::setApplicationPath(char *path)
{
	strcpy(exePath, path);
}

// process one volume in the pipeline
void FriendProcess::realtimePipelineStep(char *rtPrefix, char *format, char *actualBaseline)
{
   char mcfile[BUFF_SIZE], mcgfile[BUFF_SIZE], inFile[BUFF_SIZE], outFile[BUFF_SIZE], CmdLn[BUFF_SIZE], matOldName[BUFF_SIZE], matNewName[BUFF_SIZE], number[50];
   if (isReadyNextFile(vdb.actualImg, rtPrefix, format, inFile))
   {
      vdb.logObject->writeLog(1, "Processing file = %s\n", inFile);
	  sprintf(number, format, vdb.actualImg);
	  sprintf(inFile, "%s%s", rtPrefix, number);
	  vdb.getMCVolumeName(mcfile, number);
	  vdb.getMCGVolumeName(mcgfile, number);
	  vdb.getFinalVolumeName(outFile, number);
	  sprintf(CmdLn, "mcflirt -in %s -reffile %s -out %s %s", inFile, vdb.motionRefVolume, mcfile, vdb.mcflirtParams);


	  // mcFlirt
	  mcflirt(CmdLn);
	  
	  // copying the .mat file
	  sprintf(matOldName, "%s%s%c%s", mcfile, ".mat", PATHSEPCHAR, "MAT_0000");
	  sprintf(matNewName,  "%s%s%s", rtPrefix, number, "_mc.mat");
	  rename(matOldName, matNewName);
	  sprintf(matOldName, "%s%s", mcfile, ".mat");
	  rmdir(matOldName);

	  // subtraction process
	  // ending of the block ?
	  if (vdb.actualImg > vdb.interval.intervals[vdb.actualInterval].end)
	  {
		  if (!vdb.skipMeanSubtraction)
		  {
			  // Baseline Calculation, if this is a baseline block
			  if (isBaselineCondition(vdb.interval.intervals[vdb.actualInterval].condition))
			  {
				  vdb.logObject->writeLog(1, "Baseline mean calculation\n");
				  sprintf(actualBaseline, "%s%s%d%s", vdb.preprocDir, "mc_bsl", (vdb.actualInterval + 1), "mean.nii");
				  baselineCalculation(vdb.actualInterval, actualBaseline);
			  }
		  }
		  // going to the next interval
		  vdb.actualInterval++;
	  };
	  // actual subtraction
	  if (actualBaseline[0] != 0)
	  {
		  sprintf(CmdLn, "fslmaths %s -sub %s %s", mcfile, actualBaseline, outFile);
		  fslmaths(CmdLn);

	  }
	  else // if no baseline mean already calculated, zeros volume. Note zeroing the supposed subtracted volume
	  {
		 sprintf(CmdLn, "fslmaths %s -mul 0 %s", mcfile, outFile);
		 fslmaths(CmdLn);
	  }

	  // gaussian filtering the subtracted volume
	  sprintf(CmdLn, "fslmaths %s -kernel gauss %f -fmean %s", outFile, (float) vdb.FWHM/2.3548, outFile);
	  fslmaths(CmdLn);

	  // gaussian filtering the motion corrected volume
	  sprintf(CmdLn, "fslmaths %s -kernel gauss %f -fmean %s", mcfile, (float) vdb.FWHM/2.3548, mcgfile);
	  fslmaths(CmdLn);

	  // send graph params to FRONT END
	  sendGraphParams(mcfile, number);
	  pHandler.callAfterPreprocessingFunction(vdb, vdb.actualImg, outFile);
	  
	  if (vdb.sessionPointer != NULL)
		 if (vdb.sessionPointer->getFeedbackResponses)
		 {
			float classNum, feedBackResponse;
			test(vdb.actualImg, classNum, feedBackResponse);
			vdb.sessionPointer->processFeedback(vdb.actualImg, classNum, feedBackResponse);
		 }
	  vdb.actualImg++;
   }
}

// initializing control variables
void FriendProcess::initializeStates()
{
   vdb.initializeStates();
}

// initializing all other variables
void FriendProcess::prepRealtimeVars()
{
   if (!vdb.rPrepVars) 
   {
	  vdb.prepRealtimeVars();
	  vdb.rPrepVars=1;
   }
}

// change a config file variable value
void FriendProcess::setVar(char *var, char *value)
{
   vdb.setVar(var, value);
}

// clean up memory the code before exit
void FriendProcess::cleanUp()
{

}

// make the last steps after the processing of the run
void FriendProcess::wrapUpRun()
{
   char newpreprocDir[BUFF_SIZE], tempVolume[BUFF_SIZE];
	  
   // renaming the preproc directory
   sprintf(newpreprocDir, "%s%s%s",  vdb.outputDir, "preproc", vdb.trainFeatureSuffix);
   if (vdb.rPipeline)
   {
	  if (fileExists(newpreprocDir)) removeDirectory(newpreprocDir);
	  rename(vdb.preprocDir, newpreprocDir);

	  // copying source dir into subject directory
	  char destDir[BUFF_SIZE], sourceDir[BUFF_SIZE];
	  extractFilePath(vdb.rawVolumePrefix, sourceDir);
	  sprintf(destDir, "%s%ccopied%s", vdb.outputDir, PATHSEPCHAR, vdb.trainFeatureSuffix);
	  copyDirectory(sourceDir, destDir);
	  pHandler.callFinalFunction(vdb);
   }
   else vdb.logObject->writeLog(1, "Pipeline not executed.\n");
   
}

void FriendProcess::setLogObject(LogObject * logO)
{
	vdb.logObject = logO;
}

void FriendProcess::setSessionPointer(Session *sessionPtr)
{
   vdb.sessionPointer = sessionPtr;
   if (sessionPtr != NULL)
	  sessionPtr->setVDBPointer(&vdb);
}

int FriendProcess::isConfigRead()
{
	return vdb.rIniRead;
}

// automatic feedback calculations
void FriendProcess::setFeedbackCalculation(int automaticCalculations)
{
   if (vdb.sessionPointer != NULL)
	  vdb.sessionPointer->getFeedbackResponses = automaticCalculations;
}

// sets the status phase 0 : begin 1 : end
void FriendProcess::setPhaseStatus(string phase, int status)
{
   if (vdb.sessionPointer != NULL)
	  vdb.sessionPointer->setCommandResponse(phase, status);
}

// gets the status phase
void FriendProcess::getPhaseStatus(string phase, char *response)
{
   if (vdb.sessionPointer != NULL)
	  vdb.sessionPointer->getCommandResponse(phase, response);
}

// steps before realtime processing
void FriendProcess::prepRealTime()
{
   char arqAxial[BUFF_SIZE] = { }, CmdLn[BUFF_SIZE] = { }, betAnat[BUFF_SIZE] = { };
   size_t buffSize = BUFF_SIZE-1;
   
   if (!vdb.rPrepVars)
	   prepRealtimeVars();

   if ((!fileExists(vdb.raiFile)) || (!fileExists(vdb.rfiFile))) return;
   vdb.createDirectories();

   pHandler.callInitFunction(vdb);

   if (!fileExists(vdb.baseImage))
   {
	  // Anatomic Processing
	  snprintf(arqAxial, buffSize, "%s%s", vdb.inputDir, "RAI_ax.nii");
	  resampleVolume(vdb.raiFile, arqAxial, 1, 1, 1, vdb.TR, 0);
	  axial(arqAxial, arqAxial);
	  centralizeVolume(arqAxial, vdb.baseImage);
   }

   // creating a resampled anatomic with the same voxel dim as the functional, for registering the pipelines outcomes with this volume for presenting in FRONT END. Actually the FRIEND engine does not do this corregistering
   if (!fileExists(vdb.baseFunctional))
	  equalVoxelDim(vdb.baseImage, vdb.rfiFile, vdb.baseFunctional, vdb.TR, 0);

   // bet Functional
   if (!fileExists(vdb.maskFile))
   {
	  snprintf(CmdLn, buffSize, "bet %s %s %s", vdb.rfiFile, vdb.maskFile, vdb.betParameters);
	  bet(CmdLn);

	  snprintf(CmdLn, buffSize, "bet %s %s %s", vdb.maskFile, vdb.maskFile, vdb.betParameters);
	  bet(CmdLn);
   }
 
   // Coregistering RFI to RAI
   if (!fileExists(vdb.matrixFile))
   {
	  snprintf(betAnat, buffSize, "%s%s", vdb.inputDir, "RAI_ax_cube_sks.nii");

	  snprintf(CmdLn, buffSize, "bet %s %s", vdb.baseImage, betAnat);
	  bet(CmdLn);

	  snprintf(CmdLn, buffSize, "bet %s %s", betAnat, betAnat);
	  bet(CmdLn);

	  snprintf(CmdLn, buffSize, "flirt -ref %s -in %s -dof 7 -omat %s", betAnat, vdb.maskFile, vdb.matrixFile);
	  flirt(CmdLn);
   }
   vdb.rPreProc=true;
}

