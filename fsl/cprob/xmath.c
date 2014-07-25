#include <math.h>
#include <float.h>

int Xisnan(double x) { return _isnan(x); }
int Xisfinite(double x) { return _finite(x) && !_isnan(x); }
