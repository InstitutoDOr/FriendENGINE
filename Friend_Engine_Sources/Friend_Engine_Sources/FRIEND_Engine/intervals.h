#ifndef INTERVALS_H
#define INTERVALS_H
#include "defs.h"
#include <vector>
#include <string>
#include <stdlib.h>

typedef struct
{
	int start, end;
	char condition[255];
} Interval;

typedef struct
{   char name[255];
   float  *contrastVector;
} Contrast;

using namespace std;

string trim(const std::string &s);

class DesignObject
{
   public:

         char glmDir[BUFF_SIZE];
         char baselineCondition[100];
         vector <Interval> intervals;
         vector <int> classes;
         vector <string> conditionNames;
		 vector <string> conditionNamesWithoutBaseline;
		 vector <Contrast> contrasts;
   
		 // is the first baseline block? 
		 int IsFirstBaselineBlock(int index);

         // tests if a condition is the baseline condition
         bool isBaselineCondition(const string condition);
   
         // tests if a condition is the baseline condition
         bool isBaselineCondition(const char *condition);
   
         // verifies if condition is the baseline condition
         bool isBaselineCondition(int idx);
   
         // returns the first baseline interval prior to `interval`
         int getLastBaseline(int interval);
   
         // saves the name of the conditions in a file
         void saveConditionNames(char *file);
   
         // get condition name of a volume index
         const char *getCondition(int idx);
   
         // populates the conditionNames vector with the design file information and sort
         int getConditionsList();
   
         // returns the used volume indexes, taking in account offset and baseline conditions
         int getVolumeIndices(int offset, int start, int end, vector <int> &indexes,int NoFirstInterval=0);
   
         // get a condition volume indexes
         int getConditionVolumeIndexes(vector <int> &classVector, int classIndex, vector <int> &indexes, vector <int> &classIndexes);
   
         // returns the box car regressor of a condition
         int getConditionBoxCar(vector <float> &boxcar, const char*condition, int start=0, int end=0);
   
         // generates box car files for each condition for feat_model
         void generateConditionsBoxCar(char *outputdir, int inicio=0, int end=0);
   
         // applies an shift in the array classes. This function will be used for accomodating haemodynamic delays
         void shiftClassArray(vector <int> &classesVector, int shift);
   
         // returns the condition index in the conditionNames vector of the interval of the given volume index
         int getClass(int index);

         // populates the class array
         int getClassArray(vector <int> &classArray);
   
         // returns the index of name in vector list
         int stringIndex(vector <string> &list, char* name);
   
         // return the interval of a volume index
         int returnInterval(int index);

		 // returns the first activation index 
		 int firstActivationIndex();

		 // is this index the start of an activation block? 
		 int IsActivationBlockStart(int index);

		 // is this index the start of a block? 
		 int IsBlockStart(int index);
   
         // returns the max volume index in the design file
         int maxIndex();
   
         // reads a design file
         int readDesignFile(char *designFile);
   
		 // no Interest Contrast
		 void noInterestContrastsIndexes(char *noInteresetConditions, vector<int> &indexes);

         // populates the contrasts and contrast vectors
         void generateContrasts(bool conditionContrasts = 1, bool temporalDerivatives = 0);
   
         // generates a simple design matrix
         void generateDesignMatrix(char *fileName, int inicio=0, int end=0);
   
         // generates the FSF file necessary to run feat_model
         void generateFSFFile(char *fileName, int numDynamics, bool motionParams = 0, bool temporalDerivatives = 0);
   
         // function that cleans up the memory
         void cleanUp();

#ifndef __GNUC__
		DesignObject () {};
#else //The GCC way
		DesignObject () {};
		virtual ~DesignObject() {cleanUp();};
#endif

};

#endif
