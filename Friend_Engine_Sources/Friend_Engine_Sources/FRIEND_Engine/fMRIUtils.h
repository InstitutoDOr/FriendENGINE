#include "fMRIProcessing.h"

short getDataType(char *iname);
void volumeProperties(char *iname, int &runSize, short &datatype);
template <class T> int middlevol(char *iname, int middleIdx, char *oname);
int extractVolume(char *iname, char *oname, int idx, short datatype);
void splitVolume(char *fname, char *oname);
void znormalization(Matrix& mat, const int dim = 1);