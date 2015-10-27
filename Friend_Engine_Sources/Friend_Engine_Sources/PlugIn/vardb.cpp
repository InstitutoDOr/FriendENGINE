//
//  vardb.cpp
//  
//
//  Created by Radiologia Rededor on 5/21/13.
//
//

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>

#ifdef WINDOWS
#define snprintf _snprintf
#include <direct.h>
#include <windows.h>
#endif
#include "vardb.h"

using namespace std;

void replaceAll(string& str, const string& from, const string& to) {
	if (from.empty())
		return;
	string wsRet;
	wsRet.reserve(str.length());
	size_t start_pos = 0, pos;
	while ((pos = str.find(from, start_pos)) != string::npos) {
		wsRet += str.substr(start_pos, pos - start_pos);
		wsRet += to;
		pos += from.length();
		start_pos = pos;
	}
	wsRet += str.substr(start_pos);
	str.swap(wsRet); // faster than str = wsRet;
}

// create log file
void studyParams::initializeLogFile()
{
	char logName[500];
	sprintf(logName, "%s%c%s%s.txt", logDir, PATHSEPCHAR, "outputLog", trainFeatureSuffix);
	outputLog = fopen(logName, "wt+");
}

// closes the log file
void studyParams::closeLogFile()
{
	if (outputLog)
		fclose(outputLog);
}

// writes the message in the log file
void studyParams::writeLog(int inScreen, const char * format, ...)
{
	va_list args;
	va_start(args, format);
	vfprintf(outputLog, format, args);
	if (inScreen)
		vfprintf(stderr, format, args);
	va_end(args);
}

// returns the class of a volume index. It calls the respective DesignObject function
int studyParams::getClass(int idx)
{
    return interval.getClass(idx);
}

// returns the volume file name given the number index in a string. Returns three types :
// 1 - motion corrected volume
// 2 - motion corrected and gaussian filtered
// 3 - motion corrected, gaussian filtered and baseline mean subtracted
void studyParams::getVolume(char *outfile, char *number, int type)
{
    if (type==1) getMCVolumeName(outfile, number);
    else if (type==2) getMCGVolumeName(outfile, number);
    else if (type==3) getFinalVolumeName(outfile, number);
}

// same as above, but given number in int format
void studyParams::getVolume(char *outfile, int number, int type)
{
    if (type==1) getMCVolumeName(outfile, number);
    else if (type==2) getMCGVolumeName(outfile, number);
    else if (type==3) getFinalVolumeName(outfile, number);
}

// returns the printf format for the number part
void studyParams::getFormat(char *format)
{
   sprintf(format, "%c0%dd", '%', numberWidth);
}

// the following three functions returns the concatenation of a volume filename type
// motion corrected
void studyParams::getMCVolumeName(char *outFile, char *number)
{
   sprintf(outFile, "%s%s%s", preprocVolumePrefix, number, "_mc.nii");
}

// motion corrected and gaussian filtered
void studyParams::getMCGVolumeName(char *outFile, char *number)
{
   sprintf(outFile, "%s%s%s", preprocVolumePrefix, number, "_mc_g.nii");
}

// returns the motion corrected gausian file name format. To retrieve the volume file name, you have to first call a sprintf function e.g. sprintf(fileName, format, 1);
void studyParams::getMCGVolumeFormat(char *format)
{
	sprintf(format, "%s%c0%dd%s", preprocVolumePrefix, '%', numberWidth, ".nii");
}

// motion corrected, gaussian filtered, baseline mean subtracted
void studyParams::getFinalVolumeName(char *outFile, char *number)
{
   sprintf(outFile, "%s%s%s%s", preprocDir, "mc_ms_G_3sw_", number, ".nii");
}

// returns the final volume file name format. To retrieve the volume file name, you have to first call a sprintf function e.g. sprintf(fileName, format, 1);
void studyParams::getFinalVolumeFormat(char *format)
{
	sprintf(format, "%s%s%c0%dd%s", preprocDir, "mc_ms_G_3sw_", '%', numberWidth, ".nii");
}

// the following three functions returns the concatenation of a volume filename type. The diference from above functions is in `number` parameter
void studyParams::setActivationFile(int number)
{
	char format[50], numberString[50];
	getFormat(format);
	sprintf(numberString, format, number);
	sprintf(activationFile, "%s%s%s%s", activationsDir, "activation_", numberString, ".nii");
}

void studyParams::getMCVolumeName(char *outFile, int number)
{
   char format[50], numberString[50];
   getFormat(format);
   sprintf(numberString, format, number);
   getMCVolumeName(outFile, numberString);
}

void studyParams::getMCGVolumeName(char *outFile, int number)
{
    char format[50], numberString[50];
    getFormat(format);
    sprintf(numberString, format, number);
    getMCGVolumeName(outFile, numberString);
}

void studyParams::getFinalVolumeName(char *outFile, int number)
{
    char format[50], numberString[50];
    getFormat(format);
    sprintf(numberString, format, number);
    getFinalVolumeName(outFile, numberString);
}

// returns the prefix for the raw volumes (internal function)
void studyParams::_getPreprocVolumePrefix(char *pref)
{
    char *prefix = extractFileName(rawVolumePrefix);
    sprintf(pref, "%s%s", preprocDir, prefix);
    free(prefix);
}

// returns the prefix for the raw volumes
void studyParams::getPreprocVolumePrefix(char *pref)
{
    _getPreprocVolumePrefix(preprocVolumePrefix);
    sprintf(pref, "%s", preprocVolumePrefix);
}

// reads a config file
void studyParams::readConfigFile(char *configFile)
{
    CSimpleIniA ini;
    SI_Error rc = ini.LoadFile(configFile);
    if (rc < 0)
    {
        fprintf(stderr, "Error reading config file = %s\n", configFile);
        return;
    }
	strcpy(configFileNameRead, configFile);
    readConfig(ini);
}

// reads a config file from a buffer
void studyParams::readConfigBuffer(char *buffer, int size)
{
    CSimpleIniA ini;
    SI_Error rc = ini.LoadData(buffer, size);
    if (rc < 0)
    {
        fprintf(stderr, "Error reading config buffer");
        return;
    }
    readConfig(ini);
}

// reads a config file in a SimpleIni variable
void studyParams::readConfig(CSimpleIniA &ini)
{
	std::string buffer;

	ini.Save(buffer);
	readedIni.LoadData(buffer);
	readConfigVars();
}

// reads the information known by the studyParams object
void studyParams::readConfigVars()
  {
	randomRun = 0;
	volumesToSkip = 1;
	skipMeanSubtraction = 0;

	volumesToSkip = readedIni.GetDoubleValue("FRIEND", "VolumesAllowedToSkip", volumesToSkip);
	skipMeanSubtraction = readedIni.GetDoubleValue("FRIEND", "skipMeanSubtraction", skipMeanSubtraction);

	strcpy(subject, readedIni.GetValue("FRIEND", "SUBJECT"));
   
   strcpy(raiFile, readedIni.GetValue("FRIEND", "RAI"));
   expandFilename(raiFile);
   
   strcpy(rfiFile, readedIni.GetValue("FRIEND", "RFI"));
   expandFilename(rfiFile);
   
   strcpy(rawVolumePrefix, readedIni.GetValue("FRIEND", "Prefix"));
   expandFilename(rawVolumePrefix);
   
   strcpy(designFile, readedIni.GetValue("FRIEND", "Design"));
   expandFilename(designFile);
   
   strcpy(studyDir, readedIni.GetValue("FRIEND", "StudyDir"));
   expandFilename(studyDir);
   
   TR = readedIni.GetDoubleValue("FRIEND", "TR", 2);
   slidingWindowSize = readedIni.GetLongValue("FRIEND", "SlidingWindowSize", 3);
   runSize = readedIni.GetLongValue("FRIEND", "FuncVolumes");
   strcpy(baselineCondition, readedIni.GetValue("FRIEND", "BaselineCondition", "NEUTRAL"));
   offset = readedIni.GetLongValue("FRIEND", "Offset", 4);
   tTestCutOff = readedIni.GetDoubleValue("FRIEND", "tTestCutOff", 5.5);
   pvalueCutOff = readedIni.GetDoubleValue("FRIEND", "pvalueCutOff", 0.005);
   clusterSize = readedIni.GetLongValue("FRIEND", "ClusterSize", clusterSize);
   
   sprintf(trainFeatureSuffix, "_%s", readedIni.GetValue("FRIEND", "CurrentRunSuffix"));
   sprintf(testFeatureSuffix,  "_%s", readedIni.GetValue("FRIEND", "ModelRunSuffix"));
   
   thresholdType = readedIni.GetLongValue("FRIEND", "ByCutOff", 0);
   includeMotionParameters = readedIni.GetLongValue("FRIEND", "IncludeMotionParameters", 1);
   referenceFirstVolSequence = readedIni.GetLongValue("FRIEND", "ReferenceFirstVolSequence", 0);
   conditionContrasts = readedIni.GetLongValue("FRIEND", "ConditionContrasts", 1);
   averageMeanOffset = readedIni.GetLongValue("FRIEND", "AverageMeanOffset", 0);
   strcpy(mniTemplate, readedIni.GetValue("FRIEND", "MNITemplate", mniTemplate));
   strcpy(mniMask, readedIni.GetValue("FRIEND", "MNIMask"));
   trainingPercentage = readedIni.GetDoubleValue("FRIEND", "TrainingPercentage", trainingPercentage);
   FWHM = readedIni.GetDoubleValue("FRIEND", "FWHM", FWHM);
   performSUSAN = readedIni.GetLongValue("FRIEND", "PerformSUSAN", 0);
   percentileHigherVoxels = readedIni.GetDoubleValue("FRIEND", "PercentileHigherVoxels", percentileHigherVoxels);
   performNeurofeedback = readedIni.GetLongValue("FRIEND", "PerformNeurofeedback", 1);
   useWholeSubjectSpaceMask = readedIni.GetLongValue("FRIEND", "UseWholesubjectSpaceMask", 1);
   referenceWholeVolume = readedIni.GetLongValue("FRIEND", "ReferenceWholeVolume", 1);
   storePredictions = readedIni.GetLongValue("FRIEND", "StorePredictions", 0);
   feedBackType = readedIni.GetLongValue("FRIEND", "FeedBackType", feedBackType); // 1 : SVM 2 : Bold Level 3 : Functional Correlation
   numberWidth = readedIni.GetLongValue("FRIEND", "SufixNumberWidth", numberWidth);
   
   invX = readedIni.GetLongValue("FRIEND", "InvX", invX);
   invY = readedIni.GetLongValue("FRIEND", "InvY", invY);
   invZ = readedIni.GetLongValue("FRIEND", "InvZ", invZ);
   rIniRead = 1;
}

// saves a config information in a buffer to disk
void studyParams::saveConfigBuffer(char *buffer, int size, char *configfile)
{
    CSimpleIniA ini;
    SI_Error rc = ini.LoadData(buffer, size);
    if (rc < 0)
    {
        fprintf(stderr, "Error reading Config Buffer");
        return;
    }
    readConfig(ini);
    ini.SaveFile(configfile, false);
}

// creates the necessary directories
void studyParams::createDirectories()
{
	fprintf(stderr, "Creating study dir\n");
#ifdef WIN32
	_mkdir(studyDir);
#else
	mkdir(studyDir, 0777); // notice that 777 is different than 0777
#endif

	fprintf(stderr, "Creating subject dir\n");
#ifdef WIN32
	_mkdir(outputDir);
#else
	mkdir(outputDir, 0777); // notice that 777 is different than 0777
#endif

	fprintf(stderr, "Creating activation dir\n");
#ifdef WIN32
	_mkdir(activationsDir);
#else
	mkdir(activationsDir, 0777); // notice that 777 is different than 0777
#endif

	fprintf(stderr, "Creating input dir\n");
#ifdef WIN32
	_mkdir(inputDir);
#else
	mkdir(inputDir, 0777); // notice that 777 is different than 0777
#endif

	fprintf(stderr, "Creating glm dir\n");
#ifdef WIN32
	_mkdir(glmDir);
#else
	mkdir(glmDir, 0777); // notice that 777 is different than 0777
#endif

	fprintf(stderr, "Creating log dir\n");
#ifdef WIN32
	_mkdir(logDir);
#else
	mkdir(logDir, 0777); // notice that 777 is different than 0777
#endif

	fprintf(stderr, "Creating preproc dir\n");
#ifdef WIN32
	_mkdir(preprocDir);
#else
	mkdir(preprocDir, 0777); // notice that 777 is different than 0777
#endif
}

// initializes the internal variables after the config file load
void studyParams::prepRealtimeVars()
{
    size_t buffSize = BUFF_SIZE-1;
    includeTrailingPathDelimiter(studyDir);
    snprintf(outputDir, buffSize, "%s%s%c", studyDir, subject, PATHSEPCHAR);

    snprintf(activationFile, buffSize, "%s%s", outputDir, "temp.nii");
   
    snprintf(predictionFile, buffSize, "%s%s", outputDir, "prediction.nii");
   
    snprintf(predictionFileAux, buffSize, "%s%s", outputDir, "predictionaux.nii");
    
	snprintf(activationsDir, buffSize, "%s%s%s%c", outputDir, "activations", trainFeatureSuffix, PATHSEPCHAR);
	
	snprintf(inputDir, buffSize, "%s%s%c", outputDir, "input", PATHSEPCHAR);

	snprintf(baseImage, buffSize, "%s%s", inputDir, "RAI_ax_cube.nii");
   
    snprintf(matrixFile, buffSize, "%s%s", inputDir, "RFI2RAI_ax_cube_sks.mat");
   
    snprintf(maskFile, buffSize, "%s%s", inputDir, "RFI_sks.nii");
   
    snprintf(baseFunctional, buffSize, "%s%s", inputDir, "RAI_ax_cube_sks_rspl2RFI.nii");
   
    snprintf(activationFile, buffSize, "%s%s", inputDir, "RAI_ax_cube.nii");
    
    snprintf(glmDir, buffSize, "%s%s%c", outputDir, "glm", PATHSEPCHAR);
    
    snprintf(contrastFile, buffSize, "%s%s%s%s", glmDir, "constrast", trainFeatureSuffix, ".txt");
   
    snprintf(glmTOutput, buffSize, "%s%s%s", glmDir, "tstats", trainFeatureSuffix);
   
	snprintf(glmZOutput, buffSize, "%s%s%s", glmDir, "zstats", trainFeatureSuffix);

	snprintf(featuresSuffix, buffSize, "%s%s", glmDir, "tstats_features");
   
    snprintf(featuresAllSuffix, buffSize, "%s%s", glmDir, "tstats_features_ALL");
   
    snprintf(featuresTrainSuffix, buffSize, "%s%s", featuresSuffix, trainFeatureSuffix);
   
    snprintf(featuresAllTrainSuffix, buffSize, "%s%s", featuresAllSuffix, trainFeatureSuffix);
   
    snprintf(train4DFile, buffSize, "%s%s%s%s", glmDir, "4D_mc_ms_G", trainFeatureSuffix, ".nii");
   
    snprintf(trainGLM4DFile, buffSize, "%s%s%s%s", glmDir, "GLM_4D_mc_ms_G", trainFeatureSuffix, ".nii");
   
    snprintf(glmMatrixFile, buffSize, "%s%s%s%s", glmDir, "design_with2Gamma", trainFeatureSuffix, ".txt");
   
    snprintf(logDir, buffSize, "%s%s%c", outputDir, "log", PATHSEPCHAR);
    
   snprintf(preprocDir, buffSize, "%s%s%c", outputDir, "preproc", PATHSEPCHAR);

   // setting Design object internal variables
    fprintf(stderr, "reading design file %s\n", designFile);
    interval.readDesignFile(designFile);
	fprintf(stderr, "Copying information to Design Object\n");
	strcpy(interval.glmDir, glmDir);
    strcpy(interval.baselineCondition, baselineCondition);
	fprintf(stderr, "generating contrasts \n");
    interval.generateContrasts(conditionContrasts, 0);
    sprintf(fsfFile, "%s%s%s%s", glmDir, subject, trainFeatureSuffix, ".fsf");
    sprintf(parFile, "%s%s%s%s", glmDir, "confounds", trainFeatureSuffix, ".txt");
    sprintf(rmsFile, "%s%s%s%s", logDir, "rms_abs", trainFeatureSuffix, ".rms");
    getPreprocVolumePrefix(preprocVolumePrefix);
	fprintf(stderr, "Number of conditions : %d\n", interval.conditionNames.size() > 0);
}

// sets the value of a variable. Not used right now
void studyParams::setVar(char *var, char *value)
{
	string temp = value, from, to;

	from = "glmdir"; to = glmDir;
	replaceAll(temp, from, to);

	from = "outputdir"; to = outputDir;
	replaceAll(temp, from, to);

	from = "studydir"; to = studyDir;
	replaceAll(temp, from, to);

	from = "inputdir"; to = inputDir;
	replaceAll(temp, from, to);

	from = "logdir"; to = logDir;
	replaceAll(temp, from, to);

	from = "preprocdir"; to = preprocDir;
	replaceAll(temp, from, to);

	readedIni.SetValue("FRIEND", var, temp.c_str(), NULL, true);
	/*
    strToUpper(var);
	if (strcmp(var, "MASKFILE") == 0) strcpy(maskFile, temp.c_str());
	if (strcmp(var, "PREFIX") == 0) strcpy(rawVolumePrefix, temp.c_str());
	if (strcmp(var, "DESIGN") == 0) strcpy(designFile, temp.c_str());
	if (strcmp(var, "CURRENTRUNSUFFIX") == 0) strcpy(trainFeatureSuffix, temp.c_str());
	if (strcmp(var, "MODELRUNSUFFIX") == 0) strcpy(testFeatureSuffix, temp.c_str());
	*/
	readConfigVars();
	prepRealtimeVars();
}

// determines, given a volume index, if its in a baseline block. Calls the respective DesignObject function
bool studyParams::isBaselineCondition(int idx)
{
   return interval.isBaselineCondition(idx);
}

// get the condition name, given the volume index
const char *studyParams::getCondition(int idx)
{
   return interval.getCondition(idx);
}

// sets the FRONT END response socket
void studyParams::setSocketfd(int Sock)
{
   socks.setSocketfd(Sock);
}

// initializating control variables
void studyParams::initializeStates()
{
   rIniRead=0;
   rPrepVars=0;
   rPreProc=0;
   rGLM=0;
   rFeatureSel=0;
   rTrain=0;
}

