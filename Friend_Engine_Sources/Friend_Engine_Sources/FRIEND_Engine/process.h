#ifdef CUDAENGINE
#include "realTimeCorrection.cuh"
#endif

#include "cctypes.h"
#include "defs.h"
#include "intervals.h"
#include "vardb.h"
#include "PlugInHandler.h"


#define poly_order 1
#define num_regressors 0

// Class responsible for control all the core FRIEND engine functionality
class FriendProcess
{
private:

	// holds the engine executable path
	char exePath[BUFF_SIZE];
	int passes;

#ifdef CUDAENGINE
	RealTimeCorrection rtc;
	void *mask;
	float *volumeVector;
	float *confounds;
#endif

	// holds plugin information
	PluginHandler pHandler;

	// holds config file information
	studyParams vdb;

	// verifies if the condition is a baseline condition
	BOOL isBaselineCondition(char * condition);

	// calculates the mean volume of a baseline block
	void baselineCalculation(int intervalindex, char *baseline);

	// verifies if the next file is ready to be processed by FRIEND
	BOOL isReadyNextFile(int index, char *rtPrefix, char *format, char *infile);

	// verifies if the next file is ready to be processed by FRIEND
	BOOL isReadyNextFileCore(int indexIn, int indexOut, char *rtPrefix, char *format, char *inFile);

	// generate a file concatenating all the processed volumes rms values
	void generateRmsFile(char *dPrefix, int ini, int end, char *output);

	// generate a file concatenating all movement parameters of the processed volumes
	void generateConfoundFile(char *dPrefix, int ini, int end, char *output);

	// clean up memory the code before exit
	void cleanUp();

	// send params to frontend
	void sendGraphParams(char *mcfile, char *number);

	// plugIn related functions
	void loadLibrary();
	void unLoadLibrary();

	// process one volume in the pipeline
	void realtimePipelineStep(char *rtPrefix, char *format, char *actualbaseline);
public:

	// config file functions
	void readConfigFile(char *configFile);
	void readConfigBuffer(char *buffer, int size);
	void readConfig(CSimpleIniA &ini);
	void saveConfigBuffer(char *buffer, int size, char *configfile);

	// set socket data for response
	void setSocketfd(int Sock);


	// Preprocessing steps
	// initializating control variables
	void initializeStates();

	// initializating all other variables
	void prepRealtimeVars();

	// steps before realtime processing
	void prepRealTime();
	void runRealtimePipeline();

	// preprare files and run fsl_glm on processed volumes
	void glm();

	// performs the feature selection step
	void featureSelection();

	// the following functions need the correct plugin registration to work
	void train();
	void test(int index, float &classnum, float &projection);

	// change a config file variable value
	void setVar(char *var, char *value);

	// make the last steps after the processing of the  run
	void wrapUpRun();

	// set the engine executable directory
	void setApplicationPath(char *path);

	// functions related to the plugIn
	// set the library file and function names of the plugIn
	void loadFunctions(char *library, char *trainFunc, char *testFunc, char *initFunc, char *finalFunc, char *volumeFunc, char *afterPreprocFunc);

	// set the path for plug in library file
	void setLibraryPath(char *path);

	// sets the session pointer
	void setSessionPointer(Session *sessionPtr);

	// automatic feedback calculations
	void setFeedbackCalculation(int automaticCalculations);

	// sets the status phase 0 : begin 1 : end
	void setPhaseStatus(string phase, int status);

	// gets the status phase
	void getPhaseStatus(string phase, char *response);

	// return in the config file was read
	int isConfigRead();


	virtual ~FriendProcess()
	{
		cleanUp();
		unLoadLibrary();
	};

};