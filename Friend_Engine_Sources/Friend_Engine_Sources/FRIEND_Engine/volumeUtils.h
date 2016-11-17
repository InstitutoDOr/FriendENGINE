#include "newimage/newimageall.h"
#include "fslio/fslio.h"


int maskSize(NEWIMAGE::volume<float> &mask);
void volume2Vector(NEWIMAGE::volume<float> &mask, NEWIMAGE::volume<float> &volSample, float *volumeVector);
void vector2Volume(NEWIMAGE::volume<float> &mask, NEWIMAGE::volume<float> &volSample, float *volumeVector);
