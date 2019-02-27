#include <fstream>
#include "defs.h"
#include "stringUtils.h"

using namespace std;

extern "C" void estimate(double initialAlpha, int ntopics, char *settingsFile, char *dataFile, char *directoryName, int runType);
extern "C" void inference(char *settingsFile, char *modelFile, char *dataFile, char *directoryName);
void saveBestTopics(char *outputDir, strings &vocabulary);
string generateLDALine(strings &vocabulary, strings &list);