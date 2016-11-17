#include "volumeUtils.h"

#define maxvoxels 10

int maskSize(NEWIMAGE::volume<float> &mask)
{
	int i = 0;
	for (int z = 0; z < mask.zsize(); z++)
		for (int y = 0; y < mask.ysize(); y++)
			for (int x = 0; x < mask.xsize(); x++)
				if (mask.value(x, y, z) > 0)
				{
					i++;
					//if (i >= maxvoxels) return i;
				}
	return i;
}

void volume2Vector(NEWIMAGE::volume<float> &mask, NEWIMAGE::volume<float> &volSample, float *volumeVector)
{
	int i = 0;
	for (int z = 0; z < volSample.zsize(); z++)
		for (int y = 0; y < volSample.ysize(); y++)
			for (int x = 0; x < volSample.xsize(); x++)
				if (mask.value(x, y, z) > 0)
				{
					volumeVector[i] = volSample.value(x, y, z);
					i++;
					//if (i >= maxvoxels) return;
				}
}

void vector2Volume(NEWIMAGE::volume<float> &mask, NEWIMAGE::volume<float> &volSample, float *volumeVector)
{
	int i = 0;
	volSample = 0;
	for (int z = 0; z < volSample.zsize(); z++)
		for (int y = 0; y < volSample.ysize(); y++)
			for (int x = 0; x < volSample.xsize(); x++)
				if (mask.value(x, y, z) > 0)
				{
					volSample.value(x, y, z) = volumeVector[i];
					i++;
					//if (i >= maxvoxels) return;
				}
}
